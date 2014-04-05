#!/usr/bin/perl -w

my $isBeta;

my $VERSION = 
    ' $Id: mapLocReporter.pl,v 1.8 2014/01/02 15:50:27 tilfordc Exp $ ';

use strict;
use BMS::ArgumentParser;
use BMS::Utilities::Escape;

my $args = BMS::ArgumentParser->new
    ( -nocgi      => $ENV{HTTP_HOST} ? 0 : 1,
      );


$args->shell_coloring();

my $esc   = BMS::Utilities::Escape->new();

my ($twPre, $twPro, $tiddlers, %anonMeth, %modules);

my @tagOrder = qw
    (title modifier created modified tags changecount
     taggly.excerpts taggly.numCols taggly.sortby taggly.sortorder);

# My POD is generally inspired by BioPerl, so I will try to follow their
# nomenclature:
my %highLvlPod = map { $_ => 1 } split(/\|/, "NAME|DESCRIPTION|MODULE DESCRIPTION|SYNOPSIS|USAGE|NOTE|SEE ALSO|FEEDBACK|AUTHOR|COPYRIGHT|DISCLAIMER|CONTRIBUTORS|APPENDIX");
my %ignorePkg = map { $_ => 1 } qw(strict vars);

foreach my $ig ($args->each_split_val(qw(ignore))) {
    next unless ($ig);
    $ignorePkg{$ig} = 1;
}

my $twfile = $args->val(qw(tw tiddlywiki));

my $auth   = "PodToTiddlyWiki";
my $anTag  = 'AnonymousMethods';
my $unCat  = "Uncategorized Methods";
my $swTok  = "##SWAPOUT##";
my $now    = `date -u +'%Y%m%d%H%M'`; $now =~ s/[\n\r]+$//;

&parse_tw($twfile);
&parse_code();
&write_tw($twfile);

sub parse_code {
    my @dirs = $args->each_split_val(qw(dir directory));
    if ($#dirs == -1) {
        $args->msg("[!]","Please supply the directory path(s) that hold the Perl code of interest with -dir");
        exit;
    }
    foreach my $dir (@dirs) {
        if (-d $dir) {
            &parse_dir($dir);
        } else {
            &parse_file($dir);
        }
    }
    if (my $req = $args->val(qw(path))) {
        # The -dir argument specifies the 'main' directory or files
        # All modules and scripts found there will be incorporated
        # into the TiddlyWiki

        # The -path argument specifies one or more 'salvage' directories
        # that may contain supporting modules. These areas will be inspected
        # for ISA/use-ed files that were not in the 'main' location(s)

        my @paths = ref($req) ? @{$req} : ($req);
        map { s/\/$// } @paths;
        my $working = 1;
        while ($working) {
            $working = 0;
            while (my ($pkg, $found) = each %modules) {
                next if ($found);
                my $subPath = '/'.join('/', split(/::/, $pkg)).".pm";
                for my $p (0..$#paths) {
                    my $fullPath = $paths[$p].$subPath;
                    if (-s $fullPath) {
                        if ($working = &parse_file( $fullPath )) {
                            # Found the package
                            last;
                        }
                    }
                }
            }
        }
    }
}

sub parse_dir {
    my $dir = shift || "";
    $dir =~ s/\/$//;
    return unless ($dir);
    foreach my $file (sort $args->read_dir( -dir => $dir,
                                            -keep => '\.(pm|pl)$')) {
        &parse_file($file);
    }
}

sub parse_file {
    my $file = shift;
    return unless ($file);
    my $short = $file;
    $short =~ s/.+\///;
    return if ($short =~ /^\./);
    return if ($ignorePkg{$short});

    my $type = "UnknownFile";
    my $parent;
    if ($file =~ /([^\/]+\.pm)$/) {
        $type = "PerlModule";
    } elsif ($file =~ /\.pl$/) {
        $type = "PerlScript";
        $parent = &_script_tiddler_struct( $file );
    }
    $args->msg("[FILE]", $file);
    open(PERL, "<$file") || $args->death
        ("Failed to read $type", $file, $!);
    
    my ($pkg, $isaTxt, $meth, $podData, $podCategory);
    my $newPkgs = 0;
    while (<PERL>) {
        if (/^\=(head\d+|pod)\s*$/ ||
            /^\=(head\d+|pod)\s+(\S.+)$/) {
            # This is the start of a new POD header
            my ($pcls, $tit) = ($1, $2 || "");
            $tit =~ s/[\n\r]+$//;
            my $pd = { type => $pcls, title => $tit };
            if ($tit =~ /^\s*([A-Z ]+)/) {
                my $capTit = $1;
                if ($highLvlPod{$capTit}) {
                    # This is a high-level POD section designed to be
                    # attached to the package or script itself
                    if (my $par = $pkg || $parent) {
                        &_check_lost_pod( $podData );
                        $podData = $par->{POD} ||= $pd;
                        $podData->{type} = 'highlevel';
                    }
                }
            }
            if ($podData && $podData->{type} eq 'highlevel' &&
                !$podData->{done}) {
                # This header is still in the high-level block
                # We should add a section block if there is a title
                if ($tit) {
                    my $num = 1;
                    if ($pcls =~ /head(\d+)/) { $num = $1; }
                    push @{$podData->{body}}, ('!' x ($num+1)). " $tit";
                    next;
                }
            }
            if ($pcls eq 'head1' || $pcls eq 'head3') {
                # head1 : High level method groups
                # head3 : Musings
                if ($tit) {
                    &_check_lost_pod( $podData );
                    $podData = $pd;
                    if ($pcls eq 'head1') {
                        if ($podCategory = &_head1_tiddler_struct($tit)) {
                            $podCategory->{POD} = $podData;
                            &_set_tiddler_parent( $podCategory, $parent );
                        }
                    } else {
                        if (my $tid = &_head3_tiddler_struct($tit)) {
                            $tid->{POD} = $podData;
                        }
                    }
                } else {
                    # $args->msg("[-]", "Untitled head1 POD section");
                }
            } elsif ($pcls eq 'head2') {
                # We will associate these with subroutines
                &_check_lost_pod( $podData );
                $podData = $pd;
            } elsif ($pcls eq 'head4') {
                # This is used for markup within another POD block
                push @{$podData->{body}}, "!!!! $tit" if ($podData);
            } else {
                $args->msg("[-]", "Uncaptured POD header '$pcls'");
            }
        } elsif ($podData && !$podData->{done}) {
            # We are in a POD block (hopefully)
            if (/^=cut/) {
                $podData->{done} = 1;
            } elsif (my $body = $podData->{body}) {
                # After the info table
                push @{$body}, $_;
            } elsif (/^\s*$/) {
                # Ignore whitespace prior to the body
            } elsif (/^\s*([A-Z][^:]{2,10}?)\s*:\s*(.+?)\s*$/) {
                # First row of info section, eg
                #  Title   : get_location_by_startend
                my ($h, $d) = (&_pod_info_header($1),
                               &_pod_info_detail($2));
                push @{$podData->{info}}, [ $h, $d ] if ($h && $d);
            } elsif (/^\s+(\S.*)$/) {
                # Presumably an indented section of an info row, eg
                # the second line of:
                # Args    : [0] Chromosome designation
                #           [1] Start coordinate
                
                # Extend the info table row
                &_pod_info_detail($1, $podData->{info});
            } else {
                # Hopefully this means the start of the text body
                $podData->{body} = [ $_ ];
            }
        } elsif (/^\s*package (\S+)/) {
            $pkg = &_package_tiddler_struct($1);
            &_set_tiddler_parent( $pkg, $parent );
            if ($type eq "PerlPackage" || $type eq "PerlModule") {
                &_check_null_tiddler($parent);
                $parent = $pkg;
            }
            my $pkName = $pkg->{ATTR}{title};
            unless ($modules{$pkName}) {
                $modules{$pkName} = 1;
                $newPkgs++;
            }
        } elsif (/^use (\S+)\s*(?:\s+\'[^\']+\')?\s*;/) {
            my $pkName = $1;
            $pkName =~ s/[\;\s]+$//;
            my $par = $pkg || $parent;
            if ($par) {
                unless ($ignorePkg{$pkName}) {
                    $par->{USE}{$pkName}++;
                    $modules{$pkName} ||= 0;
                }
            }
        } elsif (/^\s*\@ISA\s+=\s+qw\((.*)$/) {
            if ($isaTxt) {
                $args->msg("[!!]", "Multiple \@ISA lines!");
            } else {
                $isaTxt = $1;
            }
        } elsif ($isaTxt) {
            $isaTxt .= $_;
        } elsif (/^sub\s+(\S+)/) {
            my $suName = $1;
            my $mPar = $pkg || $parent;
            $meth = &_method_tiddler_struct($suName, $mPar);
            if ($meth) {
                my $pc = $podCategory || &_head1_tiddler_struct($unCat);
                &_set_tiddler_parent( $meth, $pc );
                my $cat = $pc->{ATTR}{title};
                push @{$mPar->{METHS}{$cat}}, $meth if ($mPar);
                push @{$pc->{METHS}{""}}, $meth;
                &_set_tiddler_parent($meth, $mPar);
                if ($podData && $podData->{type} eq 'head2') {
                    $meth->{POD} = $podData;
                }
                
                # die $args->branch($podData) if ($podData && $podData->{info});
                $podData = undef;
            }
        } elsif (/^\*(\S+)\s+=\s+\\&(\S+)\s*\;/) {
            # method alias by type glob
            my ($ali, $meth) = ($1, $2);
            if ($parent) {
                # warn sprintf("%s = %s [%s]\n", $ali, $meth, $parent->{ATTR}{title});
                $parent->{codealias}{$meth}{$ali} = 1;
            }
            
        }
        if ($isaTxt && $isaTxt =~ /\;\s*$/) {
            &_add_isa($pkg, $isaTxt);
            $isaTxt = "";
        }
    }
    &_check_null_tiddler($parent);
    close PERL;
    return $newPkgs;
}

sub _check_lost_pod {
    my $pd = shift;
    return unless ($pd);
    return unless ($pd->{title} || $pd->{body});
    my $typ = $pd->{type} || "";
    return if ($typ ne 'head2');
    my $ti = "Untitled Method";
    # warn $args->branch($pd);
    $args->msg("[-]", "POD Data for '$typ' discarded : $ti");
}

sub _check_null_tiddler {
    my $tid = shift;
    return unless ($tid);
    my $pd = $tid->{POD};
    return if ($pd && ($pd->{body} || $pd->{info}));
    $args->msg("[-]", "No documentation for ".$tid->{ATTR}{title});
}

sub _pod_info_header {
    my $txt = shift;
    $txt =~ s/^\s+//;
    $txt =~ s/\s+$//;
    return $txt;
}

sub _pod_info_detail {
    my ($txt, $info) = @_;
    $txt =~ s/^\s+//;
    $txt =~ s/\s+$//;
    if ($info) {
        # There is already an info table being built
        my $prior = $info->[-1][1];
        if (ref($prior)) {
            # Building a list
            $txt =~ s/^\[\d+\]\s+//;
            push @{$prior}, $txt;
        } else {
            # Extend a simple string
            $info->[-1][1] .= ($info->[-1][0] eq 'Usage' ? "\n" : " ") . $txt;
        }
        return "";
    } elsif ($txt =~ /^\[\d+\]\s+(\S.+)/) {
        # Start of a list
        return [$1];
    } else {
        return $txt;
    }
}

sub _set_tiddler_parent {
    my ($child, $parent) = @_;
    return unless ($child && $parent);
    if (ref($parent)) {
        $parent = $parent->{ATTR}{title};
    }
    my $ptxt = $parent =~ /\s/ ? "[[$parent]]" : $parent;
    $child->{ATTR}{tags}{$ptxt}++;
}

sub _add_isa {
    my ($pkg, $isaTxt) = @_;
    $isaTxt =~ s/[\s\)\;]*$//;
    foreach my $isPkg (split(/[\s\n\r]+/, $isaTxt)) {
        if ($isPkg) {
            $pkg->{ISA} ||= [];
            push @{$pkg->{ISA}}, $isPkg;
            $modules{$isPkg} ||= 0;
        }
    }
}

sub basic_tiddler {
    my $title = shift;
    return $tiddlers->{$title} ||= {
        ATTR => {
            title    => $title,
            tags     => { },
            modifier => $auth,
            created  => $now,
            modified => $now,
        }
    };
}

sub tag_tiddler {
    my $tid = shift;
    foreach my $tag (@_) {
        $tid->{ATTR}{tags}{$tag} = 1 if ($tag);
    }
}

sub _head1_tiddler_struct {
    my $podName = shift;
    $podName =~ s/[;\s]+$//;
    unless ($tiddlers->{$podName}) {
        my $tid = &basic_tiddler( $podName );
        &tag_tiddler( $tid, "MethodClass" );
    }
    return $tiddlers->{$podName};
}

sub _head3_tiddler_struct {
    my $podName = shift;
    $podName =~ s/[;\s]+$//;
    unless ($tiddlers->{$podName}) {
        my $tid = &basic_tiddler( $podName );
        &tag_tiddler( $tid, "SoftwareThoughts" );
    }
    return $tiddlers->{$podName};
}

sub _script_tiddler_struct {
    my $pkName = shift;
    $pkName =~ s/^.+\///;
    unless ($tiddlers->{$pkName}) {
        my $tid = &basic_tiddler( $pkName );
        &tag_tiddler( $tid, "PerlScript" );
    }
    return $tiddlers->{$pkName};
}

sub _package_tiddler_struct {
    my $pkName = shift;
    $pkName =~ s/[;\s]+$//;
    unless ($tiddlers->{$pkName}) {
        my $tid = &basic_tiddler( $pkName );
        &tag_tiddler( $tid, "PerlPackage" );
    }
    return $tiddlers->{$pkName};
}

sub _method_tiddler_struct {
    my ($suName, $par) = @_;
    return "" unless ($suName);
    # return "" if ($suName =~ /^_/);
    $suName .= "()";
    if ($par) {
        $par = $par->{ATTR}{title} if (ref($par));
        $suName = "${par}::$suName" if ($par);
    }
    unless ($tiddlers->{$suName}) {
        my $tid = &basic_tiddler( $suName );
        &tag_tiddler( $tid, "PerlMethod", "excludeLists" );
        $tid->{PARENT} = $par;
    }
    return $tiddlers->{$suName};
}

sub parse_tw {
    my $file = shift;
    if (!$file) {
        $args->msg("[!]", "Please provide the path to the TiddlyWiki file with -tw");
        exit;
    }
    open(TW, "<$file") || $args->death
        ("Failed to read TiddlyWiki file", $file, $!);
    my $bkup = $file;
    ($twPre, $twPro) = ("","");
    my ($tid);
    while (<TW>) {
        if (!$tiddlers) {
            # Head section of TiddlyWiki
            $twPre .= $_;
            if (/^\s*<div id="storeArea">\s*$/) {
                # Ok, now we have entered the Tiddler area
                $tiddlers = {};
            }
        } elsif ($twPro) {
            # Tail section of TiddlyWiki
            $twPro .= $_;
        } elsif (/^<div(.*\stitle=\".+)>\s*$/) {
            # Start of tiddler
            my $attr = $1;
            if ($tid) {
                $args->msg("[!!]","Tiddler not recorded",
                           &tiddler_div_start($tid));
            }
            $tid = { ATTR => {}, BODY => "" };
            while ($attr =~ /(\s+([^=\"]+)=\"([^\"]+)\")/) {
                my ($rep, $key, $val) = ($1, $2, $3);
                # warn "($rep, $key, $val)";
                if ($key eq 'tags') {
                    $tid->{ATTR}{$key}{$val} = 1;
                } else {
                    $tid->{ATTR}{$key} = $val;
                }
                $attr =~ s/\Q$rep\E//;
            }
        } elsif (/^\s*<\/div>\s*$/) {
            if ($tid) {
                # End of tiddler
                if (my $title = $tid->{ATTR}{title}) {
                    $tiddlers ||= {};
                    if ($tiddlers->{$title}) {
                        $args->msg("[!!]","Duplicate '$title' tiddlers");
                    } else {
                        $tiddlers->{$title} = $tid;
                    }
                } else {
                    $args->msg("[!!]","No title for tiddler",
                               &tiddler_div_start($tid));
                }
                $tid = undef;
            } else {
                # Should be the end of the tiddler store
                $twPro .= $_;
            }
        } elsif ($tid) {
            # Middle of tiddler
            $tid->{BODY} .= $_;
        } else {
            $args->msg("[!!]", "Unexpected content in tiddler store",
                       "'$_'");
        }
    }
    close TW;
    $bkup =~ s/\/[^\/]+$//;
    $bkup .= "/TW_Backup/".`date +'%Y-%m-%d.%H:%M:%S'`;
    $bkup =~ s/[\n\r]+$//;
    $bkup .= ".html";
    $args->assure_dir($bkup,'isFile');
    system("cp \"$file\" \"$bkup\"");
    my @tids = sort keys %{$tiddlers};
    $args->msg("[<]","Loaded: $file","Backup: $bkup",
               "Tiddlers: ".($#tids + 1));
    
}

sub tiddler_attr_text {
    my $tid = shift;
    my %atHash   = %{$tid->{ATTR}};
    my @kv;
    foreach my $key (@tagOrder) {
        if (exists $atHash{$key}) {
            my $v = $atHash{$key};
            if ($key eq 'tags') {
                if (my $r = ref($v)) {
                    if ($r eq 'ARRAY') {
                        $v = join(' ', sort @{$v});
                    } elsif ($r eq 'HASH') {
                        $v = join(' ', sort keys %{$v});
                    }
                }
            }
            push @kv, sprintf("%s=\"%s\"", $key, $esc->esc_xml_attr($v));
            delete $atHash{$key};
        }
    }
    foreach my $key (sort keys %atHash) {
        push @kv, sprintf("%s=\"%s\"", $key, $atHash{$key});
        $args->msg_once("[?]", "Unanticipated tiddler key '$key'");
    }
    return join(' ', @kv);
}

sub tiddler_div_start {
    my $tid = shift;
    return sprintf("<div %s>", &tiddler_attr_text($tid)) || "";
}

sub variant_file {
    my ($file, $mod) = @_;
    my $varF = $file;
    my $sfx = "";
    if ($varF =~ /(.+)(\.[^\.]+)$/) {
        ($varF, $sfx) = ($1, $2);
    }
    $varF .= "-$mod" . $sfx;
    return $varF;
}

sub write_tw {
    my $file = shift;
    my $newF = &variant_file($file, 'plusPod');
    &_cull_anonymous_methods();
    &set_standard_tiddlers();
    open(TW, ">$newF") || $args->death
        ("Failed to write TiddlyWiki file", $newF, $!);
    print TW $twPre;
    my @tids = sort { $a->{ATTR}{title} cmp $b->{ATTR}{title} }
    values %{$tiddlers};
    foreach my $tid (@tids) {
        print TW &tiddler_div_start($tid) ."\n";
        print TW &tiddler_body($tid);
        print TW "</div>\n";
    }
    print TW $twPro;
    close TW;
    my @msg = ("[>]", "Altered TW written", $newF, "Tiddlers: ".($#tids+1));
    my $difF = &variant_file($file, 'diff');
    system("diff \"$file\" \"$newF\" > \"$difF\"");
    push @msg, ("Differences:", $difF) if (-s $difF);
    $args->msg(@msg);
}

sub _cull_anonymous_methods {
    # Find methods that do not have any documentation
    my @tnames = keys %{$tiddlers};
    foreach my $title (@tnames) {
        my $tid = $tiddlers->{$title};
        next unless (exists $tid->{ATTR}{tags} &&
                     exists $tid->{ATTR}{tags}{PerlMethod} &&
                     $tid->{ATTR}{tags}{PerlMethod});
        if (!$tid->{POD} ||
            (!$tid->{POD}{info} && !$tid->{POD}{body})) {
            # There does not seem to be any information on this method
            delete $tiddlers->{$title};
            $anonMeth{$title} = 1;
        }
    }
}

sub tiddler_body {
    my $tid   = shift;
    my $title = $tid->{ATTR}{title} || "";
    my $txt   = $tid->{BODY} || "<pre></pre>";
    my @xtra;
    # die $args->branch($tid) if ($title eq 'BMS::SnpTracker::MapLoc');
    if (my $pd = $tid->{POD}) {
        my $root = $title;
        if ($root =~ /(.+)::([^:]+)\(\)$/) {
            # This is a subroutine, get the package as the root
            $root = $1;
            my $meth = $2;
            if (exists $tiddlers->{$root}) {
                if (my $ali = $tiddlers->{$root}{codealias}{$meth}) {
                    # die $args->branch($ali);
                    # Aliases have been defined for this method
                    my @alis = sort keys %{$ali};
                    push @{$pd->{info}}, 
                    ['Aliases',join(', ', map { "{{{$_( )}}}" } @alis)];
                }
            }
            $tid->{POD}{body} ||= 
                ["// PerlMethod found in PerlPackage $root //"];
        }
        if (my $info = $pd->{info}) {
            my @table;
            foreach my $hd (@{$info}) {
                my ($k, $v) = @{$hd};
                next unless ($v);
                next if ($k eq 'Title'); # Redundant in the wiki
                if (ref($v)) {
                    if ($#{$v} == 0) {
                        $v = &standard_crosslinks( $v->[0], $root.'::' );
                    } else {
                        my @bits;
                        for my $i (0..$#{$v}) {
                            my $p = &standard_crosslinks( $v->[$i], $root.'::' );
                            push @bits, "{{{[$i]}}} $p";
                            
                        }
                        $v = join('<br>', @bits);
                    }
                } else {
                    $v =~ s/^\s+//;
                    $v =~ s/\s+$//;
                    if ($k eq 'Usage') {
                        $v = join('<br>', map { "{{{$_}}}" } split(/\n/, $v));
                    } else {
                        # Find other packages to link to:
                        $v =~ s/\s+/ /g;
                        $v = &standard_crosslinks( $v, $root.'::' );
                    }
                }
                push @table, sprintf("| !%s|%s|", 
                                     $args->esc_xml($k),
                                     $args->esc_xml($v) );
            }
            push @xtra, join("\n", @table) unless ($#table == -1);
        }
        my @body;
        my @rawBody = @{$pd->{body} || []};
        for my $bi (0..$#rawBody) {
            my $btxt = $rawBody[$bi];
            next if (!defined $btxt);
            my $txtType = "null";
            if ($btxt =~ /^ (\*.+)/) {
                $btxt = $1 . "\n";
                $txtType = 'list';
            } elsif ($btxt =~ /^!/) {
                $txtType = 'header';
            } elsif ($btxt =~ /^[\s]*$/) {
                $txtType = 'blank';
            } elsif ($btxt =~ /^\s+\S/) {
                $txtType = 'pre';
            } elsif ($btxt =~ /^\S/) {
                $txtType = 'inline';
            } else {
            }
            # die print CORE::length($btxt)."='$btxt'" if ($txtType eq 'null' && $title eq 'BMS::SnpTracker::MapLoc::location_query()');
            if ($#body == -1) {
                # First example
                push @body, [ $txtType, $btxt ] unless
                    ($txtType eq 'blank' || $txtType eq 'null');
            } else {
                my $prior = $body[-1];
                my $pTyp  = $prior->[0];
                if ($pTyp eq $txtType) {
                    # Extend the text of the prior block
                    $prior->[1] .= "\n" if ($pTyp eq 'header');
                    $prior->[1] .= $btxt;
                } elsif ($txtType eq 'blank') {
                    if ($pTyp eq 'pre') {
                        $prior->[1] .= $btxt;
                    } elsif ($pTyp eq 'inline') {
                        $prior->[1] .= "\\n";
                    }
                } else {
                    # Add a new (different) block
                    push @body, [ $txtType, $btxt ];
                }
            }
        }
        # warn $args->branch({raw => \@rawBody, tw => \@body}) if ($title eq 'BMS::ArgumentParser');
        foreach my $bdat (@body) {
            my ($type, $btxt) = @{$bdat};
            $btxt =~ s/[\s\n\r]+$//; # Terminal whitespace is not useful
            if ($type eq 'pre') {
                # PRE formatted text
                $btxt =~ s/^[\n\r]+//; # Remove leading newlines
                if ($btxt =~ /^\s*(https?:\S+)\s*$/) {
                    # This is a single URL
                    push @xtra, '@@'.$1.'@@';
                } else {
                    push @xtra, sprintf("{{{\n%s\n}}}", $args->esc_xml($btxt));
                }
            } else {
                $btxt =~ s/^[\s\n\r]+//; # Remove leading newlines and spaces
                if ($type eq 'list') {
                    # Inline formatted
                    # Find methods we can cross-link to
                    $btxt = &standard_crosslinks( $btxt, $root . '::' );
                    $btxt =~ s/ +/ /g;
                    #$btxt =~ s/[\n\r]+$//;
                    push @xtra, $args->esc_xml($btxt);
                } elsif ($type eq 'inline') {
                    # Inline formatted
                    # Find methods we can cross-link to
                    $btxt =~ s/(\\n)+$//;
                    $btxt = &standard_crosslinks( $btxt, $root . '::' );
                    # Most newlines become spaces:
                    $btxt =~ s/[\s\n\r]+/ /g;
                    # Explicit newlines are kept:
                    $btxt =~ s/(\\n)+/\n/g;
                    push @xtra, $args->esc_xml($btxt);
                } elsif ($type eq 'header') {
                    push @xtra, $args->esc_xml($btxt); #."\n";
                } else {
                    $args->msg("[-]", "Ignoring '$type' POD body block");
                }
            }
        }
    }
    if (my $ar = $tid->{ISA}) {
        my @list;
        foreach my $pkg (sort @{$ar}) {
            push @list, "* [[$pkg]]";
            delete $tid->{USE}{$pkg};
        }
        push @xtra, "!!! ISA Inherited Modules\n".
            $args->esc_xml(join("\n", @list))
            unless ($#list == -1);
    }
    if (my $h = $tid->{USE}) {
        my @list = map { "* [[$_]]" } sort keys %{$h};
        push @xtra, "!!! Used Modules\n".
            $args->esc_xml(join("\n", @list))
            unless ($#list == -1);
    }
    if (my $mh = $tid->{METHS}) {
        my @list;
        my $sep = " &diams; ";
        my @cats = sort keys %{$mh};
        foreach my $cat (@cats) {
            my $ar = $mh->{$cat};
            my @cl;
            foreach my $meth (@{$ar}) {
                my $mtid = $meth->{ATTR}{title};
                my $show = $mtid; $show =~ s/^.*:://;
                if (exists $anonMeth{$mtid}) {
                    push @cl, [ lc($show), &_anon_meth_link($show)]
                        unless ($show =~ /^_/);
                } else {
                    push @cl, [ lc($show), "[[$show|$mtid]]"];
                }
            }
            if ($#cl != -1) {
                @cl = map { $_->[1] } sort { $a->[0] cmp $b->[0] } @cl;
                if ($cat) {
                    my $wiki = ($cat eq $unCat && $#cats == 0) ? "" :
                        "* @@[[$cat]]@@<br>";
                    $wiki .= join($sep, @cl);
                    push @list, $wiki;
                } else {
                    push @list, join("\n", map { "* $_" } @cl);
                }
            }
        }
        push @xtra, "!!! Methods\n".$args->esc_xml(join("\n", @list))
            unless ($#list == -1);
    }
    return $txt if ($#xtra == -1);
    # There are auto-generated data to be added
    if ($txt =~ /<pre>\n(.+---)\n/) {
        # This is a manually entered block
        unshift @xtra, $1;
    }
    return "<pre>".join("\n", @xtra)."</pre>\n";
}

sub _anon_meth_link {
    my $txt = shift;
    return "{{undocumented{[[$txt|$anTag]]}}}";
}

sub standard_crosslinks {
    my $txt  = shift;
    my $prfx = shift || "";

    my @swap;
    # First find methods
    while ($txt =~ /(\S+\(\))/) {
        my $found = $1;
        my $rep   = $swTok . ($#swap+1);
        $txt      =~ s/\Q$found\E/$rep/g;
        # If the method is not already fully specified with a package,
        # add the package prefix to it:
        my $chk   = ($found =~ /\:\:/) ? $found : $prfx . $found;
        # For display, take out the package:
        my $short = $found; $short =~ s/.+\://;
        if (exists $tiddlers->{$chk} && $tiddlers->{$chk}) {
            # This method exists and has documentation
            $found = "[[$short|$chk]]";
            # $found = $found eq $chk ? "[[$found]]" : "[[$found|$chk]]";
        } elsif ($anonMeth{$chk}) {
            # This is an anonymous method (one with no documentation)
            $found = &_anon_meth_link($found);
        }
        push @swap, $found;
    }
   
    # Now find basic packages
    while ($txt =~ /(\S+\:\:\S+)/) {
        my $found = $1;
        my $rep   = $swTok . ($#swap+1);
        $txt      =~ s/\Q$found\E/$rep/g;
        if (exists $tiddlers->{$found} && $tiddlers->{$found}) {
            $found = "[[$found]]";
        }
        push @swap, $found;
    }

    # Reassemble the linked wiki text:
    while ($txt =~ /(\Q$swTok\E(\d+))/) {
        my ($rep, $ind) = ($1, $2);
        my $sw = $swap[$ind];
        $txt =~ s/\Q$rep\E/$sw/g;
    }
    return $txt;
    
}

sub _tw_crosslink {
    my ($txt, $regexp, $prfx) = @_;
    $prfx ||= "";
    my @swap;
    while ($txt =~ /($regexp)/) {
        my $oMeth = $1;
        my $rep = $swTok . ($#swap+1);
        $txt =~ s/\Q$oMeth\E/$rep/g;
        my $chk = $prfx . $oMeth;
        if (exists $tiddlers->{$chk} && $tiddlers->{$chk}) {
            
            $oMeth = $oMeth eq $chk ? "[[$oMeth]]" : "[[$oMeth|$chk]]";
        } elsif ($anonMeth{$chk}) {
            # This is an anonymous method (one with no documentation)
            $oMeth = &_anon_meth_link($oMeth);
        }
        push @swap, $oMeth;
    }
    # warn $args->branch(\@swap) if ($title eq 'BMS::SnpTracker::MapLoc::text_id_to_objects()');
    while ($txt =~ /(\Q$swTok\E(\d+))/) {
        my ($rep, $ind) = ($1, $2);
        my $sw = $swap[$ind];
        $txt =~ s/\Q$rep\E/$sw/g;
    }
    return $txt;
}

sub set_standard_tiddlers {
    my $clobber = $args->val(qw(clobber));
    my $stnd = <<STNDEOF;

$auth Entries tagged with '~$auth' have been generated by an automatic parsing script that analyzes Perl files and extracts out documentation stored in them. If you are a programmer attempting to write software using these packages, then such documentation will hopefully be of great use to you. If you are a researcher that wishes to just use already existing tools to analyze your data, then these topics will likely be completely irrelevant. 

PerlPackage This tag indicates that the topic describe the techinical aspects of a "Perl Package", a piece of software in the Perl programming language. You may safely ignore such topics if you are not interested in programming Perl. { taggly.numCols => '3' }

PerlScript This tag is applied to "Perl scripts", which are programs that can be run from the command line, and sometimes from a web server. Entries in this category may contain information describing how to use the tool. They will also have some technical information for programmers looking to extend or modify the code. { taggly.numCols => '3' }

PerlMethod Entries with this tag describe specific "methods" used by software in the Perl programming language. These will generally only be of interest to programmers. The topics will generally include technical details that explain how the method is used. { taggly.numCols => '3' }

MethodClass These entries are collections of methods that are somehow related, generally because they perform similar tasks, consume similar input, or produce similar output. They are organizational categories to group related code. { taggly.numCols => '3' }

excludeLists This is simply a TiddlyWiki flag that is used to keep some entries from appearing in lists, such as on the right-hand side of the Wiki. { taggly.numCols => '3' }

STNDEOF

    foreach my $line (split(/[\n\r]+/, $stnd)) {
        if ($line =~ /^(\S+) (.+)/) {
            my ($ti, $desc) = ($1, $2);
            next if ($tiddlers->{$ti} && !$clobber);
            my $tid = &basic_tiddler( $ti );
            if ($desc =~ /(.+)\s+\{([^\}]+)\}\s*$/) {
                $desc = $1;
                my $attr = $2;
                while ($attr =~ /((\S+)\s+=>\s+\'([^\']+)\')/) {
                    my ($rep, $k, $v) = ($1, $2, $3);
                    $attr =~ s/\Q$rep\E//;
                    $tid->{ATTR}{$k} = $v;
                }
            }
            $tid->{POD} = { body => [$desc] };
        }
    }

    my @anMeth = sort keys %anonMeth;
    delete $tiddlers->{$anTag};
    unless ($#anMeth == -1) {
        my $anTid = &basic_tiddler( $anTag );
        my @body = ("These methods do not have any documentation assigned to them; method calls will appear {{undocumented{like_this()}}}. This is fundamentally due to a lack of time, but in many cases the methods are either self-explanatory, mundane, or not intended for programatic access");
        my %struct;
        my $sep = " &diams; ";
        foreach my $meth (@anMeth) {
            my $src = "Unknown";
            if ($meth =~ /(.+)::([^:]+)$/) {
                ($src, $meth) = ($1, $2);
            }
            push @{$struct{$src}}, $meth;
        }
        foreach my $cat (sort keys %struct) {
            my @meths = map { "{{undocumented{$_}}}" } sort @{$struct{$cat}};
            push @body, sprintf("* @@[[%s]]@@<br>%s", $cat, join
                                ($sep, @meths));
        }
        $anTid->{BODY} = "<pre>".join("\n", @body)."</pre>";
    }
}
