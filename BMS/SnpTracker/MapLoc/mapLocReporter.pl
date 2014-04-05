#!/usr/bin/perl -w

=head1 DESCRIPTION

=head1 SYNOPSIS

=head1 AUTHOR

Charles Tilford <podmail@biocode.fastmail.fm>

//Subject __must__ include 'Perl' to escape mail filters//

=head1 LICENSE

Copyright 2014 Charles Tilford

 http://mit-license.org/

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
  BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
  ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
  CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.

=cut

my $isBeta;

my $VERSION = 
    ' $Id: mapLocReporter.pl,v 1.8 2014/01/02 15:50:27 tilfordc Exp $ ';

use strict;
use BMS::SnpTracker::MapLoc;
use BMS::FriendlyGraph;
use BMS::ArgumentParser;
use BMS::SnpTracker::MapLoc::Reporter;
use BMS::SnpTracker::MapLoc::PopulationFilter;
use Math::Trig ':pi';

srand( time() ^ ($$ + ($$<<15)) );

my $pfile = "";
foreach my $try ('/tmp', '/scratch', $ENV{HOME}) {
    if ($try && -d $try) {
        $pfile = $try;
        last;
    }
}
if ($pfile) {
    my $user = $ENV{'REMOTE_USER'} || $ENV{'LDAP_USER'} || $ENV{'USER'} ||
        $ENV{'HTTP_CN'} || $ENV{'LOGNAME'} || "DefaultUser";
    $pfile .= sprintf("/MapLocPrefs-%s.param", $user);
}

my $args = BMS::ArgumentParser->new
    ( -nocgi      => $ENV{HTTP_HOST} ? 0 : 1,
      -mode       => 'auto',
      -minimize   => 1,
      -splitfeat  => 0,
      -noanonfeat => 1,
      -maxfeatlen => 20,
      -fastSort   => 10000,
      # -noensembl  => 85,
      -tallywin   => 36,
      -tallymin   => 1.5,
      -tallydec   => 0,
      -assignmode => 'Query if possible',
      -build      => 'GRCh37',
      -errormail  => 'charles.tilford@bms.com',
      -slop       => 5000,
      -nodbsnpundef => 1,
      -instance   => 'maploc2',
      -minfreq    => undef,
      -range      => 500,
      -casesensitive => 0,
      -paramfile  => [ $pfile,  'BMS::SnpTracker::MapLoc'],
      -paramalias => {
          asnmode  => [qw(assgnmode assignmode)],
          clean    => [qw(pure)],
          exononly => [qw(nointron)],
          explain  => [qw(explainsql)],
          freq     => [qw(minfreq impactfreq)],
          fullhtml => [qw(forcehtml)],
          htmltmp  => [qw(htmltemp)],
          maxfeat  => [qw(maxfeatlength maxfeatlen)],
          output   => [qw(outfile)],
          vb       => [qw(verbose)],
          forcebld => [qw(forcebuild)],
          showbench => [qw(bench dobench benchmark)],
          splitfeat => [qw(splitfeature)],
          pretty    => [qw(ispretty)],
          query    => [qw(id ids input
                          snp snps poly polys loc location var variant
                          rna rnas transcript transcripts
                          gene genes locus loci 
                          pop population populations 
                          cell line lines cells cellline celllines
                          footprint 
                          showrna showgene)],
                      },
      );

$args->xss_protect('all');

# These are arguments which can not be changed via GET or POST:
my @protectedArgs = qw(cxurl toolurl errormail htmltmp favicon
                       imagedir urlmap pgport pghost tiddlywiki);

$args->default_only(@protectedArgs);

$args->debug->skip_key( [qw(MAPLOC sequences PARENT)], 'global' );

my $nocgi     = $args->val(qw(nocgi));
my $vb        = $args->val('vb'); 
$vb           = 1 if (!defined $vb && $nocgi);
my $mode      = &_stnd_mode($args->val(qw(mode)));
my $buildReq  = $args->val(qw(build)) || "";
my $forceBld  = $args->val('forcebld') ? $buildReq : undef;
my $limit     = $args->val(qw(limit)) || 0;
my $compact   = $args->val(qw(compact));
my $explSQL   = $args->val('explain') || 0;
my $dumpSQL   = $args->val(qw(dumpsql)) || $explSQL || 0;
my $slop      = $args->val(qw(slop)) || 0;
my $range     = $args->val(qw(range));
my $query     = $args->rawval(qw(query)) || "";
my $splitAcc  = $args->val(qw(splitacc)) || 0;
my $exonOnly  = $args->val(qw(exononly)) || 0;
my $toolUrl   = $args->val(qw(toolurl)) || "mapLocReporter.pl";
my $splitFeat = $args->val('splitfeat') || 0;
my $noAnonFeat = $args->val(qw(noanonfeat)) ? 1 : 0;
my $maxFeatLen = $args->val(qw(maxfeat)) || "";
my $clobber   = $args->val(qw(clobber));
my $isPretty  = $args->val('pretty');
my $mlInt     = $args->val(qw(instance));
my $minFreq   = $args->val('freq');
my $onlyHigh  = $args->val(qw(highonly onlyhigh));
my $outFile   = $args->val(qw(output));
my $howbad    = $args->val(qw(howbad)) || 0;
my $drawTab   = $args->val(qw(drawtable showtable)) || 0;
my $showBench = $args->val(qw(showbench));
my $assgnMode = $args->val('asnmode');
my $keepNull  = 1;
my $showRNA   = $args->val(qw(showrna)) || "";
my $tallyAll  = $args->val(qw(tallyall alltally)) || 0;
my @ccParam   = qw(customcolor customcolors custcol);
my $custColor = $args->val(@ccParam) || "";
my $debug     = $args->val(qw(debug));
my $cxDebug   = $args->val(qw(cxdebug));
my $noEns     = $args->val(qw(noens noensembl)) || 0;
my $tallWin   = $args->val(qw(tallywindow tallywin)) || 0;
my $tallMin   = !$tallWin ? 0 : $args->val(qw(tallyminimum tallymin)) || 0;
my $tallDec   = $args->val(qw(tallydecrease tallydec)) || 0;
my $isCS      = $args->val(qw(casesensitive));
my $noDbAlleles = $args->val(qw(nodballeles));
my $fh        = *STDOUT;
my $domName   = $args->val(qw(showdomain));
my $htmlTmp   = $args->val(qw(htmltmp));
my $fullHTML  = $args->val('fullhtml') ? 1 : 0;
my $clean     = ($args->val('pure') || $mode eq 'Object Report' || $nocgi) 
    ? 1 : 0;
my $noHTML    = 0;
my $allowEns;
my $allowEnsText = "";
my $toolUrlDir   = $toolUrl; 
$toolUrlDir      =~ s/\/?[^\/]+$//;

# Category key aliases:
my $catKeyMap = {
    w        => 'w',
    weight   => 'w',
    priority => 'w',
    imp      => 'imp',
    impact   => 'imp',
    freq     => 'freq',
    minfreq  => 'freq',
    nullok   => 'nullok',
    oknull   => 'nullok',
    nofilter => 'nofilter',
    showall  => 'nofilter',
    track    => 'track',
    ignore   => 'ignore',
    mask     => 'ignore',
};

if (my $aeReq = $args->val(qw(allowensembl allowens))) {
    my @reqs = ref($aeReq) ? @{$aeReq} : ($aeReq);
    my @found;
    foreach my $req (@reqs) {
        foreach my $line (split(/[\n\r]+/, lc($req || ""))) {
            $line =~ s/_/ /g;
            $line =~ s/\s+/ /; $line =~ s/^ //; $line =~ s/ $//;
            push @found, $line if ($line);
        }
    }
    if ($#found != -1) {
        $allowEns = { map { $_ => 1 } @found };
        $allowEnsText = join("\n", @found);
    }
}

my $hidePop;

# my $mltest = &maploc(); die join(" + ", $mltest->impact_list());

my $urlMap;


$args->url_callback( sub {
    my $path = shift;
    if ($path =~ /^\/stf\/(.+)/) {
        return "/$1";
    } elsif ($path =~ /^\/home\/(.+)\/public_html\/(.+)/) {
        return "/~$1/$2";
    } elsif ($path =~ /^\/home\/(tilfordc\/people\/.+)/) {
        return "http://bioinformatics.bms.com/~$1";
    }
    if ($path =~ /^\/stf(.+)/) { return $1; }
    return undef;
});

my $format = &_stnd_format( $args->val(qw(mode)) );
if (!$format && $outFile && $outFile ne 'STDERR' &&
    $outFile =~ /\.([^\.]+)$/) {
    # Determine format from output suffix
    $format = &_stnd_format( $1 );
}
$format ||= $nocgi ? 'Text' : 'HTML';

if ($format =~ /^(CanvasXpress)$/) {
    $noHTML = 1;    
}

if ($outFile) {
    if ($outFile eq 'STDERR') {
        $fh = *STDERR;
        $outFile = "";
    } elsif ($format eq 'Excel') {
        # Handle this with ExcelHelper
    } elsif (open(OUTPUT, ">$outFile")) {
        $fh = *OUTPUT;
    } else {
        $args->death("Failed to write output", $outFile, $!);
    }
}

my $rejFile;
if ($rejFile = $args->val(qw(rejectfile))) {
    if ($rejFile eq '1') {
        if ($outFile) {
            $rejFile = $outFile . "-Reject.tsv";
        } else {
            $rejFile = "Rejected_Locations.tsv";
        }
    }
    if (open( REJ, ">$rejFile")) {
        print REJ "# Rejection details\n";
        close REJ;
    } else {
        $args->err("Failed to clear reject file", $rejFile, $!);
    }
}


my ($categorySets, %excel, $twHash, %mergedCategories);
my @catParents  = ('Polymorphism', 'Mutation');
my $uncatCat    = "Uncategorized";
my $usedMinFreq = $minFreq;

my $outFH     = *STDOUT;
my $rv        = {
    errors => {},
};

my ($maploc, $sortFields, $tmpPrfx, %tmpFiles, $ad, %popcache, %stuff,
    $srpt, %tally, %userQueries, $impAliases,
    %filteredObjects, %filterCache );

my ($bulkPopCat);
my $survivingMafData        = {};
my $popFilterObj            = {};
my $bulkFrequencies         = {};
my $cachedPopulationFilter  = {};
my $bulkImpFilter           = {};
my $popInfo                 = {};

# Relates category filter native keys to human readable information:
# [Human label, Tiddlywiki key, UI column order, reporting column order]
my $cfColDat = {
    cat      => ['Variant Source',   'VariantSource',   1, 2,
                 "A category of variant"],
    freq     => ['MAF',              'MinimumMAF',      5, 4,
                 "Minimum allowed minor allele frequency"],
    imp      => ['Impact Filter',    'ImpactFilter',    7, 0,
                 "Required or prohibited ('!') impact codes"],
    nullok   => ['NOK',              'NullOk',          6, 5,
                 "If variants with no frequency data are allowed"],
    track    => ['Track',            'VariantTrack',    3, 3,
                 "Graphical track to place variants in"],
    w        => ['Priority',         'VariantPriority', 4, 1,
                 "How filters in this category 'compete' with others"],
    keep     => ['Required Impact',  'RequiredImpact',  0, 6,
                 "Only variants with these codes will be shown"],
    toss     => ['Forbidden Impact', 'ForbiddenImpact', 0, 7,
                 "Variants with these codes will be excluded"],
    nofilter => ['No Filters',       'NoFilter',        2, 0,
                 "Force all variants in this category to be displayed"],
};

my $htmlEnv = $nocgi ? 0 : $clean ? 0 : 1;
my $mime = 'plain';
if ($mode eq 'Object Report') {
    $mime = 'plain';
} elsif ($format =~ /^(HTML|CanvasXpress|Excel|Network)$/) {
    $mime = 'html';
}
if ($nocgi) {
    $args->shell_colors();
    $mime = "";
} else {
    if ($mime) { 
        $args->set_mime( -mime => $mime );
        if ($mime eq 'html' || $fullHTML || $outFile) {
            $htmlEnv = 1;
        } else {
            $htmlEnv = 0;
        }
    }
    if ($format eq 'JSON') {
        $SIG{__WARN__} = sub { print STDERR join("\n", @_) };
    }
    if ($mime eq 'plain') {
        $clean ||= 1;
        $args->manage_callback('FormatCallback', sub {
            my $self = shift;
            return shift;
        }, 'isGlobal');
    }
}
# print "$mode = $format [$mime]\n";

&HTML_START(*STDOUT);
# my $m = &maploc(); print "<p>".join(' ', map { $m->impact_html_token($_) } qw(DEL FRM STP NON SYN INT SP3 SP5 UT3 UT5  NCR UNK))."</p>";
if ($args->val(qw(setdefault))) {
    if (open(PFILE, ">$pfile")) {
        &_category_settings();
        map { $args->blockquote($_, 1) } qw(defaultfilter customcolors);
        print PFILE $args->to_text
            ( -ignore => [qw(pause nocgi setdefault cxurl),
                          qw(paramfile valuefile PARAMFILE_ERR),
                          qw(gene genes rna rnas footprint),
                          qw(format),
                          qw(snp var snps poly polys) ]);
        close PFILE;
        chmod(0666, $pfile);
        $args->msg("Your defaults have been set and will be used to set up future analyses.","You may now close this window, or <a href='$toolUrl'>start a new analysis</a>.");
        &HTML_END(*STDOUT);
        exit;
    } else {
        $args->err("Failed to create parameter file", $pfile, $!);
    }
}

my $doPerma = &process();

&extra();
&HTML_PERMALINK() if ($doPerma);
&HTML_FORM();
&HTML_END(*STDOUT);
if ($outFile) {
    close OUTPUT unless ($format eq 'Excel');
    $args->msg("[OUT]", "Output written to file", $outFile);
}

sub _stnd_mode {
    my $req = lc(shift || "");
    if (!$req || $req =~ /auto/) {
        return "Automatic";
    } elsif ($req =~ /obj.*report/) {
        return "Object Report";
    } elsif ($req =~ /(gene|loc)/) {
        if ($req =~ /(samp|pat)/) {
            return "Locus-Sample";
        } else {
            return "Locus";
        }
    } elsif ($req =~ /(rna|trans)/) {
        return "RNA";
    } elsif ($req =~ /(geno)/) {
        return "Genome";
    } elsif ($req =~ /(report)/) {
        return "Report";
    } elsif ($req =~ /(coord)/) {
        return "Coordinates";
    } elsif ($req =~ /(foobar)/) {
    }
    $args->msg("[?]","Unknown report mode request '$req'");
    return "";
}

sub _stnd_format {
    my $req = lc(shift || "");
    return "" unless ($req && $req !~ /auto/);
    if ($req =~ /(htm|browse)/) {
        return 'HTML';
    } elsif ($req =~ /te?xt/) {
        return 'Text';
    } elsif ($req =~ /(net|graph)/) {
        return 'Network';
    } elsif ($req =~ /(can|xpr|cx)/) {
        return 'CanvasXpress';
    } elsif ($req =~ /(xls|excel|sheet)/) {
        return 'Excel';
    } elsif ($req =~ /(debug|dump)/) {
        return 'Debug';
    } elsif ($req =~ /(json)/) {
        return 'JSON';
    } elsif ($req =~ /(tsv|tab)/) {
        return 'TSV';
    }
    $args->msg("[?]","Unknown format request '$req'");
    return "";
}

sub extra {
    my @efiles;
    foreach my $eType (sort keys %excel) {
        my $eh = $excel{$eType};
        next unless ($eh);
        $eh->close();
        my $url = $eh->url();
        push @efiles, $url;
        if ($#efiles == 0) {
            # print $fh "<script>document.location = '".$args->esc_url_keep_slash($url)."'</script>\n";
            my $eurl = $args->esc_url_keep_slash($url);
            if ($nocgi || !$eurl) {
                $args->msg("[FILE]","Spreadsheet report created",
                           $eh->file_path());
            } else {
                print $fh "<p>Spreadsheet can be found at <a href='$eurl'>$eurl</a></p>\n";
                print $fh "<script>window.onload = function() { document.location = '$eurl' }</script>\n";
            }
            
        }
    }
    if ($maploc && $showBench) {
        my $bm = $maploc->showbench
            ( -minfrac => $showBench == 1 ? 0.01 : $showBench,
              -class   => 'tab',
              -shell   => $mime ? 0 : 1,
              -html    => ($nocgi || $mime eq 'plain') ? 0 : 1);
        if ($nocgi) {
            $args->msg("[BENCH]", $bm);
        } else {
            print $bm; #"<pre>$bm</pre>";
        }
    }
}

sub rna_panel {
    my ($rna, $locs) = @_;
    my $ml = &maploc();
    my $su = $ml->sequtils();
    $ml->bench_start();
    $rna->read();
    print $fh &usage();
    print $fh &locus_information( $rna );
    my (@alns, %alnByAnchor);
    foreach my $aln ($rna->alignments_from_parameters
                     ({ BUILD  => $forceBld || 'best',
                        HOWBAD => $howbad})) {
        next if (&_rna_filter($rna, $aln));
        push @alns, $aln;
        push @{$alnByAnchor{ $aln->chr_acc_ids() || ""}}, $aln;
    }
    if ($#alns == -1) {
        print $fh "<p class='note'>No locations found on $buildReq</p>\n";
        &panel_details();
        $ml->bench_end();
        return;
    }
    my %tracks;
    my @ex = $rna->cx_exon_part( -align => \@alns );
    map { $_->{hideName} = 1 } @ex;
    push @{$tracks{Exons}}, @ex;
    my @feats = $rna->cx_local_feature_parts
        ( -maxfrac => $maxFeatLen,
          -featfilter => \&_feature_filter);
    &add_features( \@feats, \%tracks );
    &collapse_features( \%tracks, $ml );
    my $rnaAccId = $rna->accession_id;
    my $rnaAcc   = $rna->accession;
    my $rid      = $rna->pkey();
    map { $_->anchor_to_genome() } @alns;

    while (my ($lid, $dat) = each %{$locs || {}}) {
        my ($loc, $imp, $pop, $oim, $cx, $track) = @{$dat};
        my $locid   = $loc->pkey();
        my $locAnc  = $loc->anchor_id();
        if (my $isQry = $userQueries{$loc->handle()}) {
            if ($isQry < 0) {
                $filteredObjects{"variant"}{""}{"User Excluded"}++;
                next;
            }
        }

        # Locations are being collected in exon-centric coordinates. A
        # single variant might have multiple exon-centric
        # coordinates. This can happen when an RNA is placed several
        # times in a tandem repeat cluster. Generally such a case will
        # only happen if -howbad is not zero (ie, suboptimal locations
        # are allowed):

        #  --- RnaCopy1 ------ RnaCopy2 ---------- RnaCopy3 ---
        #          ^ VAR

        # In the above case, so long as HowBad is liberal enough there
        # will be three alignment objects for the RNA on the same
        # chromosome. The variation VAR is exonic or intronc to
        # RnaCopy1, but genomic to 2 and 3


        my %dHash;
        my @locAlns = @{$alnByAnchor{ $locAnc } || []};
        foreach my $aln (@locAlns) {
            my ($int, $isoA, $isoB) = $aln->intersect
                ( -hsp => $loc->hsp(), -isob => 1 );
            if ($int) {
                my $iStr  = $int->strand();
                my $gind2 = $int->seqindex($rnaAccId) * 2;
                my @data;
                foreach my $coord ($int->coordinates()) {
                    push @data, join('..',  $coord->[$gind2], 
                                     $coord->[$gind2+1] );
                }
                if ($#data != -1) {
                    my $dkey = join("\t", @data);
                    push @{$dHash{$iStr}{$dkey}}, $aln;
                }
            } elsif ($isoB && $#{$isoB} == 0) {
                if (my $fdat = $isoB->[0][3]) {
                    if (my $hsp = $fdat->{$rnaAccId}) {
                        # We were able to slot the SNP in as non-exonic
                        my $dkey = "$hsp->[1]..$hsp->[0]";
                        push @{$dHash{0}{$dkey}}, $aln;
                    }
                }
            }
        }
        my @istrs = keys %dHash;
        if ($#istrs > 0) {
            # There are more than one strand represented.
            # Delete the zero strand (if there). Rational is to
            # present the variant just inside the query
            delete $dHash{0};
        }
        my @dataLocs;
        while (my ($iStr, $keyH) = each %dHash) {
            while (my ($hspTxt, $alnArr) = each %{$keyH}) {
                my $hsps = [map {[ split(/\.\./, $_) ]} split(/\t/, $hspTxt)];
                push @dataLocs, [ $iStr, $hsps, $alnArr ];
            }
        }
        if ($#dataLocs == -1) {
            # &preprint($loc->to_one_line());
            next;
        } elsif ($#dataLocs != 0) {
            warn "\n\nVariant assigned to multiple exon-centric locations:\n".
                join(" + ", map { join('/', @{$_->[1]}) } @dataLocs);
        }
        $cx ||= $loc->cx_genome_part
            ( -impdata   => $imp,
              -popdata   => $pop,
              -ovidata   => $oim,
              -freq      => $usedMinFreq,
              -tosspop   => $hidePop,
              -rnafilter => \&_rna_filter,
              -keepnull  => $keepNull, );
        my $isAnchored;
        my @cxs;
        $track ||= &_pick_variant_track( $cx);
        next unless ($track);
        foreach my $dbit (@dataLocs) {
            my $lcx  = $#dataLocs == 0 ? $cx : { %{$cx} };
            my ($iStr, $hsps) = @{$dbit};
            if ($#{$hsps} == 0 && $hsps->[0][0] == $hsps->[0][1] + 1) {
                $hsps = [[$hsps->[0][1],$hsps->[0][0]]];
                $lcx->{insertion} = 1;
            }
            $isAnchored ||= $iStr;
            $lcx->{data}  = $hsps;
            if ($iStr < 0) {
                &add_cx_revcom($lcx, $su);
            } else {
                &remove_cx_revcom($lcx);
            }
            &_variant_extras( $lcx );
            push @cxs, $lcx;
        }
        if ($isAnchored) {
            map { &_tally_observations( $_ ) } @cxs;
        }

        # We should re-scale the variance metrics depending on the track
        # that we are displaying the data in.
        &add_cx_to_tracks( $cx, \%tracks, $track);
    }
    &_add_tally_track(\%tracks);
    
    my ($maxC) = sort { $b <=> $a } map { $_->{data}[-1][1] } @{$tracks{'Genome'} || []};
    my ($cxHtml, $legend) = $maploc->cx_genome_panel_html
        ( -canvasid => 0,
          -width    => 900,
          -pretty   => $isPretty,
          -params   => {
              setMin => -10,
              setMax => $maxC ? $maxC + 10 : undef,
          },
          &_standard_confs(),
          -tracks   => \%tracks );
    print $fh $cxHtml;
    print $fh &panel_details( $legend );
}


sub locus_information {
    my ($obj, $via) = @_;
    return "" if ($noHTML);
    my $acc  = $obj->accession();
    # &preprint($obj->to_text());
    my $uqt = "<span class='small qrynear'>(near query)</span>";
    if (my $uq = &_is_user_query($acc)) {
        if ($uq >= 1) {
            $uqt = "<span class='qry'>(query)</span>";
        } elsif ($uq < 0) {
            $uqt = "<span class='small qryexc'>(excluded query)</span>";
        } else {
            $uqt = "<span class='small qryrel'>(query related)</span>";
        }
    }
    my @info;
    if (my $t = $obj->list_of_values('Symbol')) {
        push @info, "<span class='sym'>$t</span>";

    }
    if (my $t = $obj->simple_value('Description')) {
        push @info, "<span class='desc'>$t</span>";

    }
    if (my $ts = &taxa_span($obj->list_of_values('Taxa'))) {
        push @info, $ts;
    }
    push @info, "<i>found by $via</i>" if ($via && $acc !~ /^\Q$via\E/);
    push @info, "<a class='small' href='$toolUrl?query=$acc&mode=Gene'>Gene View</a>" if ($obj->can('obj_type') && $obj->obj_type() ne 'Gene');
    push @info, "<a class='small' href='$toolUrl?query=$acc&mode=Genome'>Genomic View</a>";
    
    my $ibits    = "";
    my $afterbit = "";
    if ($clean) {
        $ibits = join(" | ", @info) || "";
    } else {
        $ibits = join(" |\n", @info) || "";
    }
    my $html  = "<h3>$acc $uqt $ibits</h3>\n$afterbit";
    return $html;
}

sub _is_user_query {
    my $req = uc(shift || "");
    return 0 unless ($req);
    my @try = ($req);
    if ($req =~ /^(\S{5,})\.\d{1,2}$/) {
        # Presumptive versioned ID
        push @try, $1;
    }
    my $rv;
    foreach my $t (@try) {
        $rv = $userQueries{$t} if (!$rv || $rv < $userQueries{$t});
    }
    return $rv || 0;
}

sub taxa_span {
    my $tax = shift;
    return "" unless ($tax);
    $maploc ||= &maploc();
    my $col = $maploc->pastel_text_color($tax);
    return "<span class='small taxa' style='background-color: #444; color:$col'>[$tax]</span>";
}

sub panel_details {
    my ($legend) = @_;
    return "" if ($noHTML);
    my $toss = &report_toss();
    my @bits;
    push @bits, 'legend' if ($legend);
    push @bits, 'excluded data' if ($toss);
    my $details =  "<br style='clear:left;' /><div class='toggle' onclick='toggle_hide(this)'>Click to see ".join(' and ', @bits)."</div><div class='hide' style=''>";
    $details .= $legend || "";
    $details .= $toss;
    $details .= "</div>\n";
    return $details;
}

sub process {
    return 0 if (!$query || $query =~ /^[\s\n\r]*$/);
    &_category_settings();
    my $ml   = &maploc();
    my $input = $ml->classify_requests( -request   => $query,
                                        -build     => $buildReq,
                                        -howbad    => $howbad,
                                        -sensitive => $isCS, );
    my $pickMode = $mode eq 'Automatic' ? "" : $mode;
    my @foundRows;
    my $kfTags = {
        '1'  => ['&#10003;', 'text-align: center; color:green'],
        '-1' => ['&times;', 'text-align: center; color:red'],
        '0'  => ['?', 'text-align: center; font-weight: bold; background-color:silver'],
    };
    my $ifTxt = "<b>All alleles</b> will contribute to impact calculations, regardless of how rare they are.";
    if ($minFreq) {
        $ifTxt = "Alleles will only contribute to impact calculations if they represent <b>at least ".sprintf("%.1f%%", $minFreq * 100)."</b> of their assigned population";
    } elsif ($usedMinFreq) {
        my @freqs = sort {$a <=> $b} values %{$usedMinFreq};
        my ($minF, $maxF) = map { int(0.5 + 1000 * $_) / 10 } ($freqs[0], $freqs[-1]);
        my $rng = ($minF == $maxF) ? "is $minF%" : "ranges from $minF% to $maxF%";
        $ifTxt =  "Your <b>population-specific category filters</b> will be used to determine if an allele can contribute to an impact calculation. For your specified categories this $rng";
    }

    my $doObj   = $pickMode eq 'Object Report' ? 1 : 0;
    my $doCoord = $pickMode eq 'Coordinates' ? [] : 0;

    my $sparam = { };
    unless ($doCoord || $doObj) {
        $sparam->{"Location Limit"} = ['QueryLimit', $limit ? "At most <b>$limit</b> locations will be recovered" : "<b>No limit</b>, all matching locations will be reported"];
        $sparam->{"Impact Frequency"} = ['ImpactFrequency', $ifTxt];
    }
    if ($howbad) {
        $sparam->{"Suboptimal RNAs"} = ['HowBad', "Suboptimal RNA-to-genome alignments will be considered within $howbad percentage points of the best genome alignment"];
    }
    if ($noDbAlleles) {
        $sparam->{"No Database Alleles"} = ['UserOnly', "Request to only consider explicit alleles provided in query"];
    }
    if ($allowEnsText) {
        $sparam ||= {};
        $sparam->{"Ensembl RNAs"} = ['AllowedEnsembl', "Ensembl RNAs will be ignored unless their Status and BioType are one of:<br />".join("<br />", map { "<span class='smalllabel'>$_</span>" } split(/[\n\r]+/, $allowEnsText))];
    }
    
    my @typeOrder = ('Gene', 'RNA', 'Protein', 'Range',
                     'Location', 'Population', 'Unknown');
    my %typeCount = map { $_ => 0 } @typeOrder;
    my (%fetchHSP, @explicitLocation, %rnaConnections, %popReq);
    my @keyMethods = qw(acc handle pkey unversioned_accession
                        each_gene_name each_rna_name);
    # $ml->prebranch($input);
    my $foundObj = 0;
    foreach my $type (@typeOrder) {
        my @handles = sort keys %{$input->{$type} || {}};
        my $showType = $type;
        foreach my $hand (@handles) {
            my $dat = $input->{$type}{$hand};
            unless ($dat) {
                warn "Hmmm... nothing in $type '$hand'";
                next;
            }
            my ($via, $kf, $obj, $tagVal) = @{$dat};
            if ($tagVal) {
                while (my ($tag, $vH) = each %{$tagVal}) {
                    map { $obj->set_transient_tag($tag, $_) } keys %{$vH};
                    $obj->set_transient_tag("Useful Tag", $tag);
                }
            }
            
            # warn "$type : $hand = $via";
            my $viaHtml = join("<br />", map { $args->esc_xml($_) } @{$via});
            if ($type eq 'Unknown') {
                push @foundRows, [ $viaHtml, $kfTags->{'0'}, 
                                   $showType, '', "<i>Unrecognized</i>" ];
                next;
            }
            unless ($obj) {
                if (my $meth = $ml->can('get_'.lc($type))) {
                    $obj = &{$meth}($ml, $hand);
                } else {
                    $args->msg_once("Can not recover '$type' objects");
                    next;
                }
            }
            $foundObj++;
            if ($doObj) {
                print $obj->to_text();
                next;
            } elsif ($doCoord) {
                push @{$doCoord}, $obj;
                next;
            }

            $obj->read();
            # &preprint($obj->to_text());
            # print "<pre>".$args->esc_xml($obj->to_xml())."</pre>" if ($obj->can('to_xml'));
            my ($xreq, $xfnd) = ("", "");
            if ($obj->can('extra_alleles')) {
                if (my $xh = $obj->extra_alleles()) {
                    $xreq .= "<br /><span class='xtra'>Alleles += ".
                        join('/', sort keys %{$xh})."</span>";
                }
            }
            if ($obj->is_transient()) {
                $xfnd .= "<br /><span class='transient'>User-defined transient</span>";
            }
            my $desc =  $obj->description();
            if ($obj->can('taxa')) {
                if (my $ts = &taxa_span($obj->taxa())) {
                    $desc .= " ". $ts;
                }
            }
            if ($type eq 'Location') {
                if ($obj->is_transient() && $obj->width() > 3
                    && !$obj->extra_alleles()) {
                    # This is a custom 'wide' location, treat as range request
                    $showType    = 'Range';
                }
            }
            if ($showType eq 'Range') {
                my $bld  = $obj->build();
                my $chr  = $obj->chr();
                my $targ = $fetchHSP{$kf || 1}{$bld}{$chr} ||= [];
                # Pass explicit zeroes in [2] and [3] to suppress adding
                # flanks
                push @{$targ}, [$obj->lft, $obj->rgt, 0, 0, [@{$via}]];
                my $nm = sprintf("Chr%s [%s] %d bp", $chr, $bld, $obj->width);
                #&preprint($obj->to_text());
                push @foundRows,
                [$viaHtml.$xreq, $nm, $showType,
                 $kfTags->{$kf || '0'} || $kfTags->{'0'}, 
                 "$desc = User-supplied range request"];
                $typeCount{$showType}++ if ($kf > 0);
                next;
            }
            push @foundRows,
            [$viaHtml.$xreq, $obj->name(), $obj->obj_type(),
             $kfTags->{$kf || '0'} || $kfTags->{'0'}, $desc];
            if ($type eq 'Gene' || $type eq 'RNA') {
                # Register the related RNA names to this gene
                foreach my $name ($obj->each_rna_name()) {
                    $rnaConnections{$name} = 1;
                    if ($name =~ /(.+)\.\d{1,3}$/) {
                        $rnaConnections{$1} = 1;
                    }
                }
            }
            
            $kf ||= 1;
            $typeCount{$showType}++ if ($kf > 0);
            foreach my $meth ( @keyMethods ) {
                if (my $cb = $obj->can($meth)) {
                    foreach my $key (&{$cb}($obj)) {
                        next unless ($key);
                        # warn $obj->name()." $meth = $key\n";
                        $key = uc($key);
                        my @keys = ($key);
                        if ($key =~ /(.+)\.\d{1,3}$/) {
                            push @keys, $1;
                        }
                        foreach my $uk (@keys) {
                            $userQueries{$uk} = $kf if
                                (!$userQueries{$uk} ||
                                 $userQueries{$uk} < $kf);
                        }
                    }
                }
            }
            # &preprint($obj->to_text());
            if ($showType eq 'Location') {
                push @explicitLocation, $obj if ($kf > 0);
            } elsif ($showType eq 'Population') {
                my $v = ($kf < 0 ? '!' : '') . $obj->pkey(); 
                push @{$popReq{ '-pops' }}, $v;
            } elsif ($obj->can('chr_hsps')) {
                $sparam ||= {};
                $sparam->{"RNA Queries"} ||= ['ExonOnly', "Genome searches via RNAs will use <b>".($exonOnly ? "only exon" : "both exon and intron")."</b> regions"];
                my $hdat = $obj->chr_hsps( -howbad   => $howbad,
                                           -exononly => $exonOnly,
                                           -build    => $forceBld || 'best');
                # &prebranch($hdat);
                while (my ($build, $cH) = each %{$hdat}) {
                    while (my ($chr, $hsps) = each %{$cH}) {
                        my $targ = $fetchHSP{$kf}{$build}{$chr} ||= [];
                        foreach my $hsp (@{$hsps}) {
                            push @{$targ}, [$hsp->[0], $hsp->[1], undef, undef, [@{$via}]];
                        }
                        # push @{$fetchHSP{$kf}{$build}{$chr}}, @{$hsps};
                    }
                }
            } elsif ($showType eq 'Gene') {

            }
        }
    }
    if ($doObj) {
        exit unless ($foundObj);
        return 0;
    }
    if ($doCoord) {
        my $num = $#{$doCoord} + 1;
        if ($num) {
            push @foundRows, 
            [ sprintf("<i>%d entr%s", $num, $num == 1 ? 'y' : 'ies'), "", 
              "User Request", '', "Collected for query coordinate reporting" ];
        }
    } elsif ($#foundRows == -1) {
        # No queries passed
        return 0;
    }
    if (!$pickMode) {
        # Make an automatic assignment based on the kind of queries submitted
        $typeCount{"BioAccessions"} = 
            $typeCount{Gene} + $typeCount{RNA} + $typeCount{Protein};
        my $maxBio = 5;
        if ($typeCount{"BioAccessions"} > $maxBio) {
            $format = "Excel";
            $sparam->{"Visualization Mode"} = ['AutoMode', "<b>Excel Report</b> - more than $maxBio biological accessions requested" ];
        } elsif ($typeCount{Gene}) {
            $pickMode = 'Locus';
            $sparam->{"Visualization Mode"} = ['AutoMode', "<b>$pickMode</b> - User request for locus identifiers" ];
        } elsif ($typeCount{RNA} || $typeCount{RNA}) {
            $pickMode = 'RNA';
            $sparam->{"Visualization Mode"} = ['AutoMode', "<b>$pickMode</b> - User request for transcript identifiers" ];
        } elsif ($typeCount{Location}) {
            $pickMode = "Genome";
            $sparam->{"Visualization Mode"} = ['AutoMode', "<b>$pickMode</b> - User request for polymorphism/mutation identifiers/positions" ];
        } elsif ($typeCount{Range}) {
            $pickMode = "Genome";
            $sparam->{"Visualization Mode"} = ['AutoMode', "<b>$pickMode</b> - User request for one or more genomic ranges" ];
        } else {
            $pickMode = "Genome";
            $format = "Excel";
            $sparam->{"Visualization Mode"} = ['AutoMode', "<b>Excel Report</b> - User intent uncertain from query, and unsure how many discrete genomic locations will result" ];
        }
    }


    my ($segments, %queryRanges);
    if (my $hspDat = $fetchHSP{'1'}) {
        # There are locations defined to keep
        my $flank = $maploc->default_loc_to_rna_distance();
        my $expCount = 0;
        while (my ($build, $chrH) = each %{$hspDat}) {
            while (my ($chr, $hsps) = each %{$chrH}) {
                if ($flank) {
                    foreach my $hsp (@{$hsps}) {
                        # Store left/right flank for each HSP. We do this
                        # so we can discard flank if an edge gets
                        # subtracted
                        unless (defined $hsp->[2]) {
                            $hsp->[2] = $hsp->[3] = $flank;
                            $expCount++;
                        }
                    }
                }
                my $keepHsp = $maploc->add_hsps( $hsps );
                # $args->msg("[KEEP]", "$chr [$build]"); &pre_branch({in => $hsps, out => $keepHsp});
                if (my $delHsp = $fetchHSP{'-1'}{$build}{$chr}) {

                    push @{$queryRanges{'-1'}{$build}{$chr}}, @{$delHsp};
                    my $disc;
                    ($keepHsp, $disc) =
                        $maploc->subtract_hsps( $keepHsp, $delHsp );
                    if ($#{$disc} != -1) {
                        push @{$queryRanges{'-2'}{$build}{$chr}}, @{$disc};
                    }
                    # $args->msg("[TOSS]", "$chr [$build]"); &pre_branch({ toss => $delHsp, survive => $keepHsp});                    next if ($#{$keepHsp} == -1);
                }
                $segments ||= {};
                $segments->{$build}{$chr} = $keepHsp;
                push @{$queryRanges{'1'}{$build}{$chr}}, @{$keepHsp};
                #print "$chr [$build]<br />"; &prebranch($keepHsp);
                
            }
        }
        $sparam->{"Expanded Search"} ||= ['RnaRange', "Genome searches are <b>".($expCount ? "expanded by ".$ml->comma_number($flank)." bp on each side for $expCount of your queries" : "exact" )."</b>"];

    }

    if (!$pickMode || $pickMode eq 'Locus' || $pickMode eq 'RNA') {
        my @allRnas = keys %rnaConnections;
        my $noRnas  = $#allRnas == -1 ? 1 : 0;
        my $amHtml = "";
        if ($assgnMode =~ /all/i) {
            $assgnMode = "All";
            $amHtml    = "Locations will be assigned to all impacted RNAs";
        } elsif ($assgnMode =~ /(query|qry)/i && $assgnMode =~ /(only)/i) {
            $assgnMode = "Query only";
            $amHtml    = "Locations will only be attached to <span class='qry'>query</span> RNAs, otherwise they will be ignored";
            $amHtml .= ". <span class='warn'>CAUTION - this means no results will be shown, as no query RNAs/genes could be recognized</span>" if ($noRnas);
        } elsif (!$noRnas ||
                 ($assgnMode =~ /(query|qry)/i &&
                  $assgnMode =~ /(possible)/i)) {
            $assgnMode = "Query if possible";
            $amHtml    = "Locations will be attached <i>exclusively</i> to <span class='qry'>query</span> RNAs if available, otherwise will be attached to all impacted RNAs";
        } else {
            $assgnMode = "All";
            $amHtml    = "Locations will be assigned to all impacted RNAs";
        }
       
        $sparam->{"Location Assignment"} = ['LocationAssignment', "<b>$assgnMode</b> - $amHtml" ];
    }

    $sparam->{"Genome Build"} ||= ['BuildToken', $forceBld ? "All results must be from genome build <b>$forceBld</b>" : $buildReq ? "Ambiguous queries (eg symbols) will be from genome build <b>$buildReq</b>, others will use the most recent build" : "Your searches will report hits on all encountered genome builds" ];
    my @ps = sort keys %{$sparam || {}};

    my $tw = &tw_hash();
    my $infofh = $fullHTML ? $fh : *STDOUT;
    if ($clean && !$fullHTML) {
        # No details
    } elsif ($nocgi && !$fullHTML) {
        foreach my $param (@ps) {
            my ($twn, $desc) = @{$sparam->{$param}};
            $desc =~ s/<[^>]+>/ /g;
            $desc =~ s/\s+/ /g;
            $desc =~ s/^ //; $desc =~ s/ $//;
            $args->msg("[-]", $param, $desc);
        }
    } elsif ($htmlEnv) {
        print $infofh "<table class='tab'>\n";
        print $infofh "<tbody>\n";
        print $infofh " <tr><th colspan='5' style='color:#390; text-align:left'>$tw->{UserQuery}User search terms</th></tr>\n";
        print $infofh " <tr>".join("", map { "<th>$_</th>" }("Request", "Found in DB", "Type", "Keep?","Description"))."</tr>\n";
        foreach my $fr (@foundRows) {
            print $infofh "<tr>";
            foreach my $c (@{$fr}) {
                my $sty = "";
                if (ref($c)) {
                    ($c, $sty) = @{$c};
                    $sty = $sty ? " style='$sty'" : "";
                }
                print $infofh "<td$sty>$c</td>";
            }
            print $infofh "</tr>\n";
        }
        if (my $errs = $input->{Error}) {
            print $infofh "<tr><th>Error</th><td colspan='4'>".join("<br />", map { $_->[0]} @{$errs})."</td></tr>\n";
        }
        
        print $infofh " <tr><th colspan='5' style='color:#390; text-align:left'>Search Parameters</th></tr>\n";
        print $infofh " <tr><td colspan='5'>";
        unless ($#ps == -1) {
            print $infofh "<table><tbody>\n";
            foreach my $param (@ps) {
                my ($twn, $desc) = @{$sparam->{$param}};
                print $infofh "  <tr><th style='text-align:right'>$param$tw->{$twn}</th><td>$desc</td></tr>\n";
            }
            print $infofh "  </tbody></table>\n";
        }
        
        if (my $cs = &cat_settings_html('hide')) {
            print $infofh "<div class='toggle' onclick='toggle_hide(this)'>Click to show Category Filters</div>$cs" unless ($doCoord);        
        }
        print $infofh "</td></tr>\n";
        
        print $infofh "</tbody>\n";
        print $infofh "</table>\n";
    }

    if ($doCoord) {
        return &format_coordinates( $doCoord );
    }
    # &prebranch($segments);

    my $lq = $maploc->location_query
        ( -segments => $segments,
          %popReq,
          -dumpsql  => $dumpSQL,
          -limit    => $limit, );
    my $locids = $lq->{loc_id} || [];
    foreach my $expLoc (@explicitLocation) {
        if (my $pkey = $expLoc->pkey()) {
            push @{$locids}, $pkey;
        }
    }
    my $locNum = $#{$locids} + 1; #$#explicitLocation + 2;
    unless ($locNum || ($#explicitLocation + 1)) {
        $args->msg("No variant locations were found with your query");
        return 1;
    }
    # $args->msg("[DEBUG]","Primary query: ".($#{$locids} + 1)." loci");
    $locids = &filter_locations( $locids );
    # $args->msg("[DEBUG]","Fast Filter: ".($#{$locids} + 1)." loci");

    my $unkPid = $ml->_unknown_pop_id();
    my (@locs, @reject, %done);
    # Now calculate the precise impact of each location. Previously,
    # locations inside exons were just treated as impact 'RNA'
    # Now we invest the time to decvonvolute to NON, SYN, UT3 etc
    # &prebranch($survivingMafData);
    foreach my $req (@explicitLocation, @{$locids}) {
        my $loc = $req;
        unless (ref($loc)) {
            # loc_id
            $loc = $ml->get_location($req);
        }
        next if ($done{$loc->handle()}++);

        # &preprint( "$req\n" . $loc->to_text() );
        $loc->lock_alleles() if ($noDbAlleles);
        $loc->read();
        my $hand = $loc->handle();
        my $isQry   = &_is_user_query($hand);
        if ($isQry < 0) {
            $filteredObjects{"variant"}{""}{"User Excluded"}++;
            next;
        }
        my $imp = $loc->impact( -freq      => $usedMinFreq,
                                -keepnull  => $keepNull,
                                -tosspop   => $hidePop,
                                -rnafilter => \&_rna_filter, );
        my $pop = $loc->each_population_id( -tosspop => $hidePop,);
        my $oim = $loc->overall_impact( -impdata => $imp,
                                        -tosspop => $hidePop,
                                        -popdata => $pop );

        # Detailed impact filtering is a bit of a problem.  Consider a
        # SNP location with two populations, each having two alleles:

        # Pop1: AAT/AAC = Asn/Asn
        # Pop2: AAA/AAG = Lys/Lys

        # Considered from only one of either populations, the location
        # is synonymous. However, considering all alleles from both
        # populations, the location is non-synonymous. Another case:

        # Pop3: TGT/TGC = Cys/Cys
        # Pop4: TGA/TGG =  * /Trp

        # Location::impact() will calculate an 'overall' worst-case
        # impact, as well as overall impacts for each overlapping
        # RNA. However, it will not calculate per-population impacts,
        # which could be different from the overall impact, or could
        # be different between different populations.

        # I feel that more nuanced population-specific impact
        # filtration should be put in place. However, it is not clear
        # how to do so. If the user is excluding SYN variants, then
        # the first example would be filtered out by a per-population
        # filter, yet the location as a whole is non-synonymous and
        # probably of interest to them. Using a global calculation
        # that considers all available alleles should result in fewer
        # populations being excluded, presuming that the user is
        # interested in "impactful" variants.

        # So for the moment the global impact will be used. HOWEVER,
        # per-population filters will be applied, so if Pop1 allowed
        # SYN variants, but Pop2 does not, then Pop2 would be
        # excluded.

        my $tok = $oim->{impToken} || "UNK";
        my ($smds, $refilter) = @{$survivingMafData->{$hand} || [[]]};
        my %impResults;
        # my $mafs = $loc->population_maf(); &prebranch($mafs);
        if ($#{$smds} != -1) {
            # We have data that have passed the first filter
            # Collate it and refilter, if neccesary
            foreach my $smd (@{$smds}) {
                my ($pid, $maf, $pf) = @{$smd};
                my $r = ($refilter && $pf) ? $pf->test_impact( $tok ) : "";
                $pf ||= &_user_request_filter();
                my $w = $pf->weight();
                push @{$impResults{$w}}, [$r, $pf, $pid, $maf];
                # args->msg_once("$tok != ".$pf->name()) if ($r);
            }
        } else {
            # If we are here it should mean a manually provided location
            # Keep everything
            my $pf   = &_user_request_filter();
            my $mafs = $loc->population_maf();
            my @ir;
            while (my ($pid, $maf) = each %{$mafs}) {
                push @ir, ["", $pf, $pid, $maf];
            }
            if ($#ir == -1) {
                push @ir, ["", $pf, $ml->_temporary_pop_id(), -1];
            }
            $impResults{0} = \@ir;
        }

        # Find the highest-weighted category(ies)
        my ($topW) = sort { $b <=> $a } keys %impResults;
        my %topResults;
        foreach my $smd (@{$impResults{$topW || 0} || []}) {
            my $pf = $smd->[1];
            push @{$topResults{$smd->[0]}{$pf->name()}}, $smd;
        }
        if ($isQry) {
            # User request, always keep
        } elsif (!exists $topResults{""}) {
            # None of the top-weighted populations passed the more precise
            # impact test. Reject the location
            &_reject_location( \%topResults, \@reject, $hand );
            next;
        }

        # The location is "fully surviving"
        # We now need to find what tracks it should be assigned to, and
        # what MAF to use for each track
        my ($okPids, $track);
        # &prebranch(\%topResults);
        while (my ($w, $res) = each %impResults) {
            # Use 1 for the top-weighted categories, -1 for others
            my $tok = ($w == $topW) ? 1 : -1;
            foreach my $rd (@{$res}) {
                next if ($rd->[0]);
                # Passed the filter, make note of it:
                $okPids ||= {};
                my $pid = $rd->[2];
                $okPids->{$pid} = $tok;
                if ($tok > 0) {
                    # This is a top-weighted category, note which
                    # track(s) to record it into, and keep the best MAF
                    my $maf  = $rd->[3];
                    my @cats = $rd->[1]->categories();
                    @cats    = ('Uncategorized') if ($#cats == -1);
                    $track ||= {};
                    foreach my $cat (@cats) {
                        $track->{$cat} = $maf if (!defined $track->{$cat} ||
                                                 $track->{$cat} < $maf);
                    }
                }
            }
        }
        $track ||= { "Unknown" => -1 };
        unless ($isQry) {
            my %okTrk = %{$track};
            if (exists $okTrk{Discard}) {
                # A discard request will cause the location to always be excluded
                &_reject_location( { "User request to filter out category" =>
                                     {"Discard Filter" => []}}, \@reject, $hand );
                next;
            }
            if (exists $okTrk{Ignore}) {
                # An ignored location can be kept so long as at least one other
                # track is available
                delete $okTrk{Ignore};
                my @remain = keys %okTrk;
                if ($#remain == -1) {
                    &_reject_location({ "Passed, but only in ignored categories" =>
                                        {"Ignore Filter" => []}}, \@reject, $hand);
                    next;
                }
            }
        }
        my @trks = keys %{$track};


        # Calculate the CX object here for filtering
        my $cx = $loc->cx_genome_part
            ( -impdata   => $imp,
              -popdata   => $pop,
              -ovidata   => $oim,
              -freq      => $usedMinFreq,
              -tosspop   => $hidePop,
              #-revcom    => $iStr < 0 ? 1 : 0,
              -rnafilter => \&_rna_filter,
              -keepnull  => $keepNull, );
        $cx->{okPids} = $okPids if ($okPids);

        # $track ||= &_pick_variant_track( $cx);
        #unless ($track) {
        #    if ($userQueries{$hand}) {
        #        $track = ["User Variant"];
        #    } else {
        #        next;
        #    }
        #}
        push @locs, [$loc, $imp, $pop, $oim, $cx, $track];
    }
    # $args->msg("[DEBUG]","Careful Filter: ".($#locs + 1)." locations");
    &_note_rejections( \@reject, "Precise impact filters" );

    if ($format eq 'Network') {
        return &format_network( \@locs );
    } elsif ($format eq 'Text') {
        return &format_text( \@locs );
    } elsif ($format eq 'TSV') {
        return &format_tsv( \@locs );
    } elsif ($format eq 'Excel') {
        return &format_excel( \@locs );
    }
    # We are going to be generating a graphical display
    print $fh &usage();
    # Genome is easy
    unless ($pickMode =~ /(Locus|RNA)/i) {
        return &show_on_genome( \@locs, \%queryRanges );   
    }
    # Expand names / IDs in the RNA connections hash
    #foreach my $k (keys %rnaConnections) {
    #    my $nk;
    #    if ($k =~ /^(\d+)$/) {
    #        $nk = $maploc->pkey_to_text($k);
    #    } else {
    #        $nk = $maploc->text_to_pkey($k);
    #    }
    #    $rnaConnections{$nk} ||= $rnaConnections{$k} if ($nk);
    #}
        
    # Otherwise, we need to figure out the objects that are going to
    # hold the locations
    my (%useRna, %orphan);
    foreach my $dat (@locs) {
        my ($loc, $imp) = @{$dat};
        my $pkey = $loc->pkey();
        my @found;
        # See first if impact can lead us to the right RNAs
        my @checkIds   = keys %{$imp};
        # If no impacts are reported, consider all the RNAs under the location:
        @checkIds = $loc->each_rna_id() if ($#checkIds == -1);
        if ($#checkIds == -1) {
            # No nearby RNAs
            my $targ = $orphan{$pkey} ||= { loc => $loc };
            $targ->{reason}{"No nearby RNAs"}++;
            next;
        }
        my @useIds;
        if ($assgnMode eq 'All') {
            # Locations are assigned to every RNA
            @useIds = @checkIds;
        } else {
            # If any of the IDs are associated with a query, use them
            foreach my $id (@checkIds) {
                my ($acc) = $maploc->rna_ids_to_cached_accessions( $id );
                # $args->msg_once("FOO acc_id=$id : $acc");
                next unless ($acc);
                if ($rnaConnections{$acc}) {
                    push @useIds, $id;
                } elsif ($acc =~ /(.+)\.\d{1,3}$/) {
                    push @useIds, $id if ($rnaConnections{$1});
                }
            }
            if ($#useIds == -1) {
                # Could not relate to a queried gene / RNA
                if ($assgnMode eq 'Query if possible' ||
                    &_is_user_query($pkey)) {
                    # Either it is ok to use every RNA, or the variant
                    # is itself a query
                    @useIds = @checkIds;
                } else {
                    my $targ = $orphan{$pkey} ||= { loc => $loc };
                    $targ->{reason}{"Not associated with a query"}++;
                }
            }
        }
        map { $useRna{$_}{$pkey} ||= $dat } @useIds;
    }
    if ($pickMode eq 'Locus') {
        # We need to convert the RNAs to genes
        my (%useGenes, %missing, %mapped);
        while (my ($rid, $pH) = each %useRna) {
            my ($racc) = $maploc->rna_ids_to_cached_accessions( $rid );
            my $uq = $userQueries{$racc} || 0;
            my @gaccs = $maploc->rna_accession_to_gene_accessions( $racc );
            while (my ($pkey, $dat) = each %{$pH}) {
                if ($#gaccs == -1) {
                    $missing{$pkey} ||= $dat;
                } else {
                    $mapped{$pkey}++;
                    foreach my $gacc (@gaccs) {
                        $useGenes{$gacc}{$pkey} ||= $dat;
                        $userQueries{$gacc} = $uq if 
                            (! defined $userQueries{$gacc} ||
                             $userQueries{$gacc} < $uq);
                    }
                }
            }
        }
        while (my ($pkey, $dat) = each %missing) {
            next if ($mapped{$pkey});
            my $loc = $dat->[0];
            my $targ = $orphan{$loc->pkey()} ||= { loc => $loc };
            $targ->{reason}{"Unable to map over to Gene"}++;
        }
        my @gaccs = &sort_accessions( keys %useGenes );
        foreach my $gacc (@gaccs) {
            my $gene = $maploc->get_gene($gacc);
            &gene_panel( $gene, $useGenes{$gacc} );
        }
    } else {
        my @rids = &sort_accessions( keys %useRna );
        foreach my $rid (@rids) {
            my $rna = $maploc->get_rna_by_id($rid);
            &rna_panel( $rna, $useRna{$rid} );
        }
        # &show_on_genome( \@locs, \%queryRanges );
    }
    # &prebranch(\%userQueries);
    return 1 if ($noHTML);
    my @orphs = values %orphan;
    my $oNum  = $#orphs + 1;
    if ($oNum) {
        print $fh "<i>$oNum location".($oNum == 1 ? '' : 's')." could not be placed on your requested framework. </i><span class='toggle' onclick='toggle_hide(this)'>Show Orphan Locations</span><pre class='hide'>";
        foreach my $targ (@orphs) {
            print $fh $targ->{loc}->to_one_line()." (".join(', ', sort keys %{$targ->{reason}}).")\n";
        }
        print "</pre>\n";
    }
    return 1;
}

sub sort_accessions {
    my @sorter;
    foreach my $acc (@_) {
        next unless ($acc);
        push @sorter, [$acc, $userQueries{$acc} || 0,
                       $acc =~ /^ENS/ ? 0 : 1, uc($acc) ];
    }
    return map { $_->[0] } sort {
        $b->[1] <=> $a->[1] || $b->[2] <=> $a->[2] || $a cmp $b } @sorter;
}

sub gene_panel {
    my ($gene, $locs) = @_;
    my $ml = &maploc();
    my $su = $ml->sequtils();
    $ml->bench_start();
    $gene->read();
    print $fh &usage();
    print $fh &locus_information( $gene );
    my @ehsps = $gene->exon_footprints( -build => $forceBld );
    if ($#ehsps == -1) {
        if ($noHTML) {
        } else {
            print $fh "<p class='note'>No locations found".
                ($buildReq ? " on $buildReq" : "")."</p>\n";
            $ml->bench_end();
        }
        return;
    }
    if ($buildReq && $#ehsps != 0) {
        my $num = $#ehsps + 1;
        print $fh "<p class='note'>$num locations found on $buildReq</p>\n"
            unless ($noHTML);
    }
    my $gid = $gene->pkey();
    foreach my $aln (@ehsps) {
        my $chr   = $aln->chr_name();
        my $build = $aln->chr_build();
        my $str   = $aln->str();
        my ($cs, $ce) = $aln->chr_start_end();
        my $nstr = $str > 0 ? "+$str" : $str || '?';

        my @lidList = keys %{$locs || {}};
        unless ($noHTML) {
            print $fh "<h4>$chr.$build";
            print $fh " <span class='coord'>$cs..$ce</span>";
            print $fh " <span class='str'>[$nstr]</span>";
            my $lnum = $#lidList + 1;
            my $ltxt = sprintf("%d variant%s", $lnum, $lnum == 1 ? '' : 's');
            print $fh " <span class='count'>$ltxt</span>";
            print $fh "</h4>\n";
        }
        my $tracks = $gene->cx_exon_parts( -aln => $aln );
        foreach my $lid (@lidList) {
            if (my $isQry = $userQueries{$lid}) {
                if ($isQry < 0) {
                    $filteredObjects{"variant"}{""}{"User Excluded"}++;
                    next;
                }
            }
            my $dat = $locs->{$lid};
            my ($loc, $imp, $pop, $oim, $cx, $track) = @{$dat};
            my ($int, $isoA, $isoB) = $aln->intersect
                ( -hsp => $loc->hsp(), -isob => 1 );
            my $iStr = 1;
            my @data;
            my $isAnchored;
            if ($int) {
                $iStr     = $int->strand();
                my $gind2 = $int->seqindex($gid) * 2;
                foreach my $coord ($int->coordinates()) {
                    push @data, [ $coord->[$gind2], $coord->[$gind2+1] ];
                }
                $isAnchored ||= $iStr;
            } elsif ($isoB && $#{$isoB} == 0) {
                if (my $fdat = $isoB->[0][3]) {
                    if (my $hsp = $fdat->{$gid}) {
                        # We were able to slot the SNP in as non-exonic
                        push @data, $hsp;
                    }
                }
            }
            if ($#data == -1) {
                next;
            }
            $cx ||= $loc->cx_genome_part
                ( -impdata   => $imp,
                  -popdata   => $pop,
                  -ovidata   => $oim,
                  -freq      => $usedMinFreq,
                  -tosspop   => $hidePop,
                  # -revcom    => $iStr < 0 ? 1 : 0,
                  -rnafilter => \&_rna_filter,
                  -keepnull  => $keepNull, );
            $track ||= &_pick_variant_track( $cx);
            next unless ($track);
            if ($#data == 0 && $data[0][0] == $data[0][1] + 1) {
                @data = ([$data[0][1],$data[0][0]]);
                $cx->{insertion} = 1;
            }

            $cx->{data} = \@data;
            if ($iStr < 0) {
                &add_cx_revcom($cx, $su);
            } else {
                &remove_cx_revcom($cx);
            }
            &_variant_extras( $cx );
            &_tally_observations( $cx ) if ($isAnchored);
            &add_cx_to_tracks( $cx, $tracks, $track);
        }
        &_add_tally_track($tracks);

        my ($maxC) = sort { $b <=> $a } map { $_->{data}[-1][1] } @{$tracks->{'Genome'} || []};
        my ($cxHtml, $legend) = $maploc->cx_genome_panel_html
            ( -canvasid => 0,
              -width    => 900,
              -pretty   => $isPretty,
              -params   => {
                  setMin => -10,
                  setMax => $maxC ? $maxC + 10 : undef,
              },
              &_standard_confs(),
              -tracks   => $tracks );
        print $fh $cxHtml;
        print $fh &panel_details( $legend );
    }

}

sub add_cx_to_tracks {
    my ($cx, $tracks, $trackMaf) = @_;
    while (my ($tname, $maf) = each %{$trackMaf}) {
        # We need to make a local copy of each CX object
        # This is because different tracks will be representing
        # different populations, which will have different MAFs
        my %lcx = %{$cx};
        if (defined $maf) {
            $maf = 0 if ($maf < 0);
            my $sh = int(0.5 + 100 * $maf) / 100;
            $sh = 0.1 if ($sh < 0.1);
            $lcx{scaleHeight} = $sh;
            $lcx{MAF} = int(0.5 + 1000 * $maf) / 10;
        }
        push @{$tracks->{$tname}}, \%lcx;
    }
}

sub add_cx_revcom {
    my ($cx, $su) = @_;
    my $aArr = $cx->{alleles};
    return undef unless ($aArr);
    return $cx->{revcom} if ($cx->{revcom});
    my $rcH = $cx->{revcom} = {};
    foreach my $al (@{$aArr->[0] || []}) {
        my $rc = $su->revcom($al);
        $rcH->{$al} = $rc;
    }
    $cx->{alleles} = [[ sort values %{$rcH} ]];
    return $rcH;
}

sub remove_cx_revcom {
    my ($cx) = @_;
    my $rcH = $cx->{revcom};
    return unless ($rcH);
    delete $cx->{revcom};
    $cx->{alleles} = [[ keys %{$rcH} ]];
    return $rcH;
}

sub show_on_genome {
    my ($locData, $queryRanges) = @_;
    my $maploc  = &maploc();
    my $flank   = $maploc->default_loc_to_rna_distance();
    my $allVers = $args->val(qw(allvers allversions));
    my %locByChr;
    foreach my $ldat (@{$locData}) {
        my $loc = $ldat->[0];
        push @{$locByChr{$loc->build}{$loc->chr}}, $ldat;
    }
    print $fh &ref_diff_key();
    my $wsUrl = $maploc->web_support_url();
    foreach my $build (sort keys %locByChr) {
        &html_xtra("<h5>Genome Build $build</h5>\n");
        my @allC = sort keys %{$locByChr{$build}};
        for my $ac (0..$#allC) {
            my $chr = $allC[$ac];
            $maploc->bench_start('Render CX Panel');
            my @locs = @{$locByChr{$build}{$chr}};
            if ($#locs == -1) {
                $args->msg("[CODE ERROR]", "Weird - no locations for $chr");
                next;
            }
            my $chrAnc = $locs[0][0]->location_accession;
            # Tabulate all the observed populations
            my (%pidH, @cols, %tracks );
            my @locRange = map { [ $_->[0]->lft - $flank, 
                                   $_->[0]->rgt + $flank ] } @locs;
            my $chrID    = $locs[0][0]->loc_acc_id();
            
            my ($keptLoc, $minC, $maxC) = (0, 10 ** 15, -10);
            for my $li (0..$#locs) {
                my ($loc, $imp, $pop, $oim, $cx, $track) = @{$locs[$li]};
                $cx ||= $loc->cx_genome_part
                    ( -impdata   => $imp,
                      -popdata   => $pop,
                      -ovidata   => $oim,
                      -tosspop   => $hidePop,
                      -freq      => $usedMinFreq,
                      -rnafilter => \&_rna_filter,
                      -keepnull  => $keepNull, );

                $track ||= &_pick_variant_track( $cx);
                unless ($track) {
                    $filteredObjects{"variant"}{"Failed to find track"}{""}++;
                }
                &_variant_extras( $cx );
                my $hand   = $loc->handle();
                if (my $isQry = $userQueries{$hand}) {
                    if ($isQry < 0) {
                        $filteredObjects{"variant"}{""}{"User Excluded"}++;
                        next;
                    }
                }
                &_tally_observations( $cx );
                &add_cx_to_tracks( $cx, \%tracks, $track);
                $keptLoc++;
                my ($s, $e) = ($loc->start, $loc->end);
                $minC = $s if ($minC > $s);
                $maxC = $e if ($maxC < $e);
            }

            if ($queryRanges) {
                while (my ($kf, $bH) = each %{$queryRanges}) {
                    next unless (exists $bH->{$build} && $bH->{$build});
                    my $cH = $bH->{$build};
                    next unless (exists $cH->{$chr} && $cH->{$chr});
                    if (my $cxq = $maploc->cx_query_part
                        ( -hsps => $cH->{$chr},
                          -flag => $kf )) {
                        push @{$tracks{Query}}, @{$cxq};
                        push @locRange, map { 
                            [$_->[0], $_->[1] ] } @{$cH->{$chr}};
                    }
                }
            }
            
            if (1) {
                # Instead of using isolated ranges, use the whole span
                @locRange = sort { $a->[0] <=> $b->[0] ||
                                       $a->[1] <=> $b->[1] } @locRange;
                @locRange = ([ $locRange[0][0], $locRange[-1][1] ]);
            }
            my $rnaRows  = $maploc->overlapping_rnas_for_flanks
                ( $chrID, \@locRange);
            my %rnas;
            foreach my $rr (@{$rnaRows}) {
                my ($alnId, $seqId, $sc, $ptl, $ptr, $rId) = @{$rr};
                my $rna = $maploc->get_rna_by_id( $rId );
                my $accV = $rna->acc();
                my $accU = $accV;
                my $vers = 0;
                if (!$allVers && $accV =~ /(.+)\.([^\.]+)$/) {
                    ($accU, $vers) = ($1, $2);
                }
                push @{$rnas{$accU}{$vers}}, [ $rna, $alnId ];
            }
            # &prebranch(\%rnas);
            while (my ($accU, $vH) = each %rnas) {
                my ($vers) = sort { $b <=> $a } keys %{$vH};
                my $rns = $accU =~ /ENS/ ? 'Ensembl' : 'RefSeq';
                foreach my $dat (@{$rnas{$accU}{$vers}}) {
                    my ($rna, $alnId) = @{$dat};
                    my $aln = $maploc->get_alignment( $alnId );
                    # next if ($aln->howbad > $howbad);
                    next if (&_rna_filter($rna, $aln));
                    $rna->read();
                    foreach my $cx ($rna->cx_genome_part($alnId)) {
                        # next if ($cx->{anchor} ne $chrAnc);
                        $cx->{hideName} = 1 if 
                            ($rns eq 'Ensembl' && 
                             !$cx->{tags}{"Official Symbol"});
                        # &prebranch($cx);
                        push @{$tracks{$rns}}, $cx;
                    }
                    my @feats = $rna->cx_feature_parts
                        ( -aln => $aln,
                          -maxfrac => $maxFeatLen);
                    &add_features( \@feats, \%tracks );
                }
            }
            
            foreach my $rns ('RefSeq', 'Ensembl') {
                next unless (exists $tracks{$rns});
                my @sorter;
                foreach my  $dat (@{$tracks{$rns}}) {
                    my $sk = sprintf("%06.2f", $dat->{howbad} || 0);
                    if (my $lab = $dat->{label}) {
                        if ($lab =~ /(.+) \([^\d]+(\d+)/) {
                            $sk .= sprintf("%s ZZZZZ %09d", $1, $2);
                        } else {
                            $sk .= $lab;
                        }
                    } else {
                        my $id = $dat->{id};
                        if ($id =~ /^[^\d]+(\d+)/) {
                            $sk .= sprintf("ZZZZZ %09d", $1);
                        } else {
                            $sk .= "ZZZZZ $id";
                        }
                    }
                    push @sorter, [ uc($sk), $dat ];
                }
                
                $tracks{$rns} = [ map { $_->[1] } sort {
                    $a->[0] cmp $b->[0] } @sorter ];
            }
            &collapse_features( \%tracks, $maploc );
            &_add_tally_track(\%tracks);

            my $hBar = "<h3>Chromosome $chr $chrAnc";
            $hBar  .= sprintf( " <span class='mini'>%d locations</span>", 
                               $keptLoc);
            $hBar .= "</h3>\n";
            &html_xtra($hBar);

            my ($cxHtml, $legend) = $maploc->cx_genome_panel_html
                ( -canvasid => $noHTML ? undef : 0,
                  -width    => 900,
                  -pretty   => $isPretty,
                  -params   => {
                      setMin => $minC - 10,
                      setMax => $maxC + 10,
                  },
                  &_standard_confs(),
                  -tracks   => \%tracks );
            print $fh $cxHtml;
            print $fh &panel_details( $legend );
            my $brk = $clean ? "" : $ac == $#allC ? "<a target='_blank' title='Powered by CanvasXpress' href='http://canvasxpress.org'><img src='$wsUrl/images/CXbanner10pxTrans.png' /></a>" : "&hearts;";
            print $fh "<h4 style='text-align:center; clear:both'>$brk</h4>\n";
            $maploc->bench_end('Render CX Panel');
        }
    }
}

sub html_xtra {
    return unless ($htmlEnv);
    print $fh join("\n", @_);
}

sub _find_rnas {
    my $vv = shift;
    my $acc = $vv;
    my $what = "unversioned";
    my @rnas;
    my $ml = &maploc();
    if ($acc =~ /(.+)\.(\d+)$/) {
        # Versioned accesion
        $what = "versioned";
        @rnas = $ml->get_rnas( $acc );
    } else {
        # Find best version currently known
        my @accIds = $ml->text_search("$acc.%" , 'wc');
        
        my %vers;
        foreach my $va (@accIds) {
            if ($va =~ /^\Q$acc\E\.(\d+)$/) {
                $vers{$1} = $va;
            }
        }
        # Take the highest version with information:
        foreach my $v (sort { $b <=> $a } keys %vers) {
            my $accV = $vers{$v};
            @rnas = $ml->get_rnas( $accV );
            last unless ($#rnas == -1);
        }
    }
    return (\@rnas, $what);
}
sub format_coordinates {
    my $objs = shift;
    my @cols = ("Name", "Type", "DB ID", "Build", "Chr", 
                "Start", "End", "Width", "Symbols",
                "Description" );
    my $widths = {
        "Name"          => 16,
        "Type"          => 10,
        "DB ID"         => 8,
        "Build"         => 8,
        "Chr"           => 8, 
        "Start"         => 12,
        "End"           => 12,
        "Width"         => 6,
        "Symbols"       => 10,
        "Description"   => 20,
    };
    my %addCol;
    my $cMap = {
        "Name"        => 'name',
        "Type"        => 'obj_type',
        "DB ID"       => 'pkey',
        "Build"       => 'chr_build',
        "Chr"         => 'chr_names',
        "Start"       => 'chr_start',
        "End"         => 'chr_end',
        "Width"       => 'chr_width',
        "Description" => 'description',
        "Symbols"     => ['all_symbols'],
    };
    my @nullOut = ("DB ID");
    my @rows;
    my $joiner = ' / ';
    foreach my $obj (@{$objs}) {
        my %data;
        while (my ($col, $cbName) = each %{$cMap}) {
            my $isArray = 0;
            if (ref($cbName)) {
                $isArray = 1;
                $cbName = $cbName->[0];
            }
            if (my $cb = $obj->can($cbName)) {
                if ($isArray) {
                    # Get as array
                    my @vals = &{$cb}($obj);
                    $data{$col} = join($joiner, @vals);
                } else {
                    $data{$col} = &{$cb}($obj);
                }
            }
        }
        foreach my $tag ($obj->tag_values("Useful Tag")) {
            my @vals = $obj->tag_values($tag);
            unless ($#vals == -1) {
                unless ($widths->{$tag}) {
                    push @cols, $tag;
                    $widths->{$tag} = 15;
                }
                $data{$tag} = join($joiner, @vals);
            }
            
        }
        map { $data{$_} ||= "" } @nullOut;
        my @row = map { defined $_ ? $_ : "" } map { $data{$_} } @cols;
        push @rows, \@row;
        # $obj->maploc->prebranch($obj); die;
        # print "<pre>".$obj->to_text()."</pre>";
    }
    if ($format eq 'HTML') {
        print $fh "<table class='tab'>\n";
        print $fh " <caption>Coordinate report for your queries</caption>\n";
        print $fh " <tbody>\n";
        print $fh " <tr>".join("", map { "<th>$_</th>" } @cols)."</tr>\n";
        foreach my $row (@rows) {
            print $fh " <tr>".join("", map {"<td>$_</td>"} @{$row})."</tr>\n";
        }
        print $fh "</tbody></table>\n";
    } elsif ($format eq 'Excel') {
        my $sname = "Coordinate Report";
        my $eh = $excel{$sname};
        unless ($eh) {
            $outFile  ||= sprintf
                ("%s/Coordinates-%d-%d.xls", $htmlTmp, time, $$);
            $eh = $excel{$sname} = BMS::ExcelHelper->new( $outFile );
            my $url  = $args->path2url($outFile);
            $eh->url($url);
            my $sheet = $eh->sheet( -name    => $sname,
                                    -freeze  => 1,
                                    -width   => [ map { $widths->{$_} || 10 } @cols],
                                    # -formats => \@fmts,
                                    -columns => \@cols, );
        }
        map { $eh->add_row($sname, $_) } @rows;
    } else {
        print "<pre>" if ($htmlEnv);
        map { print $fh join("\t", @{$_})."\n" } ( \@cols, @rows );
        print "</pre>\n" if ($htmlEnv);
    }
    return 1;
}

sub format_text {
    my $locdata = shift;
    my $lnum   = $#{$locdata} + 1;
    $args->msg("[>]", "Exporting $lnum location".($lnum == 1 ? '':'s')." to Text");
    print "<pre style='background-color:#cc9'>" unless ($nocgi);
    foreach my $ldat (@{$locdata}) {
        my ($loc, $imp, $pop, $oim) = @{$ldat};
        print $fh $loc->to_text();
    }
    print "<pre>" unless ($nocgi);
    return 1;
}

sub format_network {
    my $locdata = shift;
    $args->msg("[!!]","Network not implemented",
               "I need to come up with an objective way to identify the samples / populations / individuals that are 'unusual'");
    return;
    my $ml   = &maploc();
    my %geneIdLookup;
    foreach my $ldat (@{$locdata}) {
        my ($loc, $imp, $pop, $oim, $cx, $track) = @{$ldat};
        my %gids;
        while (my ($rid, $idat) = each %{$cx->{impact}}) {
            my $gid = $geneIdLookup{$rid};
            unless (defined $gid) {
                $gid = 0;
                if (my $rna = $ml->get_rna_by_id( $rid )) {
                    my @genes = $rna->each_gene;
                    $gid = $genes[0]->pkey() if ($#genes == 0);
                }
                $geneIdLookup{$rid} = $gid;
            }
            $gids{$gid}{$idat->{imp}}++ if ($gid);
        }
    }
    return 1;
}

sub format_excel {
    my $locdata = shift;
    my $lnum   = $#{$locdata} + 1;
    $args->msg("[>]", "Exporting $lnum location".($lnum == 1 ? '':'s').
               " to Excel");
    my $ml   = &maploc();
    my $srpt = &reporter();
    my (%rnaLocs, %geneLocs, @allLocs);
    foreach my $ldat (@{$locdata}) {
        my ($loc, $imp, $pop, $oim, $cx) = @{$ldat};
        $cx ||= $loc->cx_genome_part
            ( -impdata   => $imp,
              -popdata   => $pop,
              -tosspop   => $hidePop,
              -ovidata   => $oim,
              -freq      => $usedMinFreq,
              -rnafilter => \&_rna_filter,
              -keepnull  => $keepNull, );
        my $lc = [$loc, $cx];
        push @allLocs, $lc;
        while (my ($rId, $rImp) = each %{$cx->{impact}}) {
            my $aln = $rImp->{align};
            push @{$rnaLocs{$rId}{$aln}}, $lc;
        }
    }
    my @varRows = $srpt->tally_locations
        ( -locs => \@allLocs, -chooserna => 1);

    my @rnaRows;
    while (my ($rid, $alnH) = each %rnaLocs) {
        my $rna = $srpt->rna_obj($rid, 'isID');
        my $acc = $rna->accession();
        next unless ($acc =~ /^[NX][MR]_/);
        my $loc = $rna->unique_value('LL');
        my $ldat = $geneLocs{$loc} ||= {};
        push @{$ldat->{rnas}}, $acc;

        while (my ($alnId, $rlocs) = each %{$alnH}) {
            my $aln = $ml->get_alignment( $alnId );
            $ldat->{alns}{$acc}{$alnId} = $aln;
            push @{$ldat->{locs}}, @{$rlocs};
            my $rrow = $srpt->rna_detail_table
                ( -rna => $rna, -locs => $rlocs, -alns => [$aln]);
            push @rnaRows, @{$rrow};
        }
    }
    my @locRows;
    while (my ($loc, $ldat) = each %geneLocs) {
        next unless ($loc);
        my %u = map { $_->[0]->pkey() => $_ } @{$ldat->{locs}};
        my $lrow = $srpt->locus_detail_table
            ( -locus => $loc,
              -rnas  => $ldat->{rnas},
              -alns  => $ldat->{alns},
              -locs  => [ values %u ] );
        push @locRows, @{$lrow};
    }
    my $eh = &eh($outFile, "Location Report");
    $srpt->add_variant_excel(\@varRows, $eh);
    $srpt->add_rna_excel(\@rnaRows, $eh);
    $srpt->add_gene_excel(\@locRows, $eh);
    return 1;
}

sub format_tsv {
    my $locdata = shift;
    my $lnum   = $#{$locdata} + 1;
    $args->msg("[>]", "Exporting $lnum location".($lnum == 1 ? '':'s')." to TSV");
    my $srpt   = &reporter();
    my $ml     = &maploc();
    my @rows;
    foreach my $ldat (@{$locdata}) {
        my ($loc, $imp, $pop, $oim) = @{$ldat};
        my $cx = $loc->cx_genome_part
            ( -impdata   => $imp,
              -popdata   => $pop,
              -ovidata   => $oim,
              -tosspop   => $hidePop,
              -freq      => $usedMinFreq,
              -rnafilter => \&_rna_filter,
              -keepnull  => $keepNull, );
        push @rows, $srpt->tally_locations
            ( -locs => [[$loc, $cx]], -chooserna => 1);
    }
    @rows = sort { $a->{Chr} cmp $b->{Chr} ||
                       $a->{ChrLft} <=> $b->{ChrLft} ||
                       $a->{ChrRgt} <=> $b->{ChrRgt} } @rows;
    print $fh $srpt->table_to_tsv
        ( -rows => \@rows, -type => 'variant' );
    if (my $toss = &report_toss()) {
        $toss =~ s/<[^>]+>//g;
    }
    return 1;
}

sub add_features {
    my ($feats, $tracks) = @_;
    foreach my $cx (@{$feats}) {
        delete $cx->{hideName} if ($domName);
        if ($splitFeat && $#{$cx->{data}} > 0) {
            $args->msg_once("Features spanning multiple locations will be broken into isolated parts");
            foreach my $loc (@{$cx->{data}}) {
                my %cpy = %{$cx};
                $cpy{data} = [ $loc ];
                push @{$tracks->{feature}}, \%cpy;
            }
        } else {
            push @{$tracks->{feature}}, $cx;
        }
    }
}

sub collapse_features {
    my ($tracks, $ml) = @_;
    return unless (exists $tracks->{feature});
    my @group = $ml->collapse_cx_set( $tracks->{feature} );
    $tracks->{feature} = \@group;
}

sub _feature_filter {
    my ($rng) = @_;
    return 1 unless ($rng);
    my $name = $rng->name();
    return 1 unless ($name);
    if ($name =~ /^[A-Z]{2,6}\:\d+$/) {
        # HPRD:04264, CDD:1234 etc
        my $ml = &maploc();
        my $obj = $ml->get_text( $name );
        $obj->read_tags();
        my $tags = $obj->all_tag_values();
        if (my $desc = $tags->{Description}) {
            # A description has been set, we will assume it's 
            return 0;
        }
        $filteredObjects{"feature"}{"No meaningful description"}{"$name"} = -1;
        return 1;
    }
    return 0;
}

sub _rna_filter {
    my ($rna, $aln) = @_;
    return 1 unless ($rna);
    my $accV = $rna->acc();
    if (my $isQry = $userQueries{uc($accV)}) {
        if ($isQry < 0) {
            $filteredObjects{"RNA"}{"User Excluded"}{$accV}++;
            return $filterCache{$accV} = 1;
        } else {
            return $filterCache{$accV} = 0;
        }
    }

    if ($aln) {
        # Alignment specific filtration
        if ($forceBld) {
            my $bld = $aln->chr_build();
            return 1 if (!$bld || uc($bld) ne uc($forceBld));
        }
        if ($aln->howbad > $howbad) {
            my $tag = $howbad ? "Genome match > ${howbad}% worse than best" :
                "Sub-optimal genome match";
            $filteredObjects{"RNA"}{$tag}{$accV} = -1;
            return 1;
        }
    }
    return $filterCache{$accV} if (defined $filterCache{$accV});
    # Keep if it is not Ensembl:
    return $filterCache{$accV} = 0 unless ($accV =~ /^ENS/);

    if ($allowEns) {
        # We are only keeping certain status and biotype combinations
        $rna->read_tags();
        my @bts = map { lc($_) } $rna->tag_values('Ensembl Biotype');
        my @ess = map { lc($_) } $rna->tag_values('Ensembl Status');
        map { s/_/ /g } @bts;
        my @flags;
        # Calculate the combined flags
        foreach my $es (@ess) {
            foreach my $bt (@bts) {
                push @flags, "$es $bt";
            }
        }
        # Also allow the flags by themselves
        push @flags, @bts, @ess;
        my $ok = 0;
        map { $ok++ if ($allowEns->{$_} ) } @flags;
        unless ($ok) {
            my $fail = $flags[0] ? "'$flags[0]'" :  "-Not Defined-";
            $filteredObjects{"RNA"}{"Status + BioType is $fail"}{$accV} = -1;
            return $filterCache{$accV} = 1;
        }
    }

    # Keep if no Ensembl filtering:
    return $filterCache{$accV} = 0 unless ($noEns);
    my $filtTag = "Less than ${noEns}% match to RefSeq";
    
    # Keep if it is an explicit query:
    my $accU = $accV; $accU =~ s/\.\d+$//;
    return $filterCache{$accV} = 0 if (&_is_user_query($accV, $accU));
    my $near = $rna->nearby_rnas( -align => $aln, -howbad => $howbad );
    return $filterCache{$accV} = 1 unless ($near);
    while (my ($gtag, $data) = each %{$near}) {
        foreach my $nb (@{$data->{nearby}}) {
            last if ($nb->{match} < $noEns);
            return $filterCache{$accV} = 0 if ($nb->{acc} =~ /^[NX][MR]/);
        }
    }
    $filteredObjects{"RNA"}{$filtTag}{$accV} = -1;
    return $filterCache{$accV} = 1;
}

sub cat_settings_html {
    my ($cls) = @_;
    my $maploc  = &maploc();
    &_category_settings();
    my $html = "";
    my @filtRows;
    my %filtCols;
    my @show = qw(w cat freq nullok track);
    foreach my $cs (values %{$categorySets}) {
        my %row = ( w => -1, cat => '' );
        foreach my $key (@show) {
            my $v = $cs->{$key};
            next if (!defined $v || $v eq '');
            if ($key eq 'nullok') {
                $v = $v ? ['&#10003;', 'text-align: center; color:green'] :
                    ['&times;', 'text-align: center; color:red'];
            }
            $row{$key} = $v;
            $filtCols{$key}++;
        }
        if (my $ih = $cs->{_impHash}) {
            while (my ($flag, $hash) = each %{$ih}) {
                my $key = $flag ? 'keep' : 'toss';
                my @bits = map {
                    $maploc->impact_html_token($_) } sort keys %{$hash};
                unless ($#bits == -1) {
                    $row{$key} = join(' ', @bits);
                    $filtCols{$key}++;
                }
            }
        }
        push @filtRows, \%row;
    }
    return $html if ($#filtRows == -1);
    my $tw = &tw_hash();
    my @colOrder = sort { $cfColDat->{$a}[3] <=> 
                              $cfColDat->{$b}[3] } keys %filtCols;
    $cls ||= "";
    $html .= "<table class='tab $cls'><tbody>\n";
    $html   .= " <tr>";
    foreach my $tok (@colOrder) {
        my $cfd  = $cfColDat->{$tok};
        $html .= sprintf("<th style='align:center' title='%s'>%s<br />%s</th>",
                         $args->esc_xml_attr($cfd->[4]), $cfd->[0], 
                         $tw->{$cfd->[1]});
    }
    $html .= "</tr>\n";

    my $naTok = ['N/A', 'color:gray'];
    my %notAp = map { $_ => $naTok } qw(freq imp nullok keep toss);
    foreach my $row (sort { $b->{w} <=> $a->{w} || 
                                $a->{cat} cmp $b->{cat} }
                     @filtRows) {
        my %using = %{$row};
        $html .= " <tr>\n";
        foreach my $tok (@colOrder) {
            $html .= "  <td";
            my $v = $using{$tok};
            if (!defined $v) {
                $v = "";
            } elsif (ref($v)) {
                $html .= " style='$v->[1]'";
                $v = $v->[0];
            } elsif ($tok eq 'track') {
                my $sty;
                if ($v eq 'Ignore') {
                    $sty = "color:white; background-color: black;"
                } elsif ($v eq 'Discard') {
                    $sty = "color:yellow; background-color: red;"
                }
                if ($sty) {
                    $html .= " style='$sty'";
                    while (my ($kk, $vv) = each %notAp) {
                        $using{$kk} = $vv;
                    }
                }
            }
            $html .= ">$v</td>\n";
        }
        $html .= " </tr>\n";
    }
    $html .= "</tbody></table>\n";
    return $html;
}

sub _category_settings {
    my $cat = shift || "";
    unless ($categorySets) {
        my $ml = &maploc();
        $ml->bench_start("Parse Category Filters");
        # Find out if explicit values have been passed, or if we should
        # use defaults
        
        $categorySets = {};
        # First cycle through "bulk text" parameters
        # These can be large chunks of multiline text via textarea
        # or selected options
        # TCGA | w:5 | freq: 0.1
        foreach my $param (qw( defaultfilter
                               catfilt catfilter categoryfilter
                               advancedfilter)) {
            my @catFilt = $args->each_split_val( $param );
            $args->clear_param( $param );
            next if ($#catFilt == -1);
            # $args->msg("[DEBUG]","Category Param: $param");
            foreach my $catDat (@catFilt) {
                next unless ($catDat);
                my @bits = split(/\s*\|\s*/, $catDat);
                my $cat = &_stnd_cat_name( shift @bits );
                next unless ($cat);
                foreach my $bit (@bits) {
                    if ($bit =~ /([^\:]+)\s*\:\s*(.+)/) {
                        &_set_category_keyval($ml, $cat, $1, $2);
                    }
                }
            }
        }
        # Now cycle through named parameters. 
        foreach my $param ($args->each_param()) {
            if ($param =~ /^cat(egory)?filt(er)?_(.+)_(.+)$/) {
                my ($cat, $key) = (&_stnd_cat_name($3),
                                   $catKeyMap->{lc($4)} || lc($4));
                if ($cat && $key) {
                    my $val = $args->val($param);
                    &_set_category_keyval($ml, $cat, $key, $val);
                }
                $args->clear_param( $param );
            }
        }

        my %catPar;
        foreach my $cName (@catParents, "Region") {
            my $cPar = "$cName Categories";
            my $cat  = $ml->get_text( $cPar );
            $cat->read();
            foreach my $pop ($cat->tag_values("Member")) {
                $catPar{$pop} = $cName;
            }
        }
        
        my $userSet = 0;
        foreach my $td (values %{$categorySets}) {
            #if ($td->{track} =~ /^(discard|ignore)$/i) {
            #    # $td->{w} = 0;
            #    # push @ignoreCat, $td->{cat};
            #    next;
            #}
            $userSet++;
        }
        # Set default values for the parent categories:
        foreach my $cat (@catParents) {
            my %defs;
            if ($cat eq 'Polymorphism') {
                %defs = ( freq => 0.05,
                          imp => '!Ugly !NonCoding' );
            } else {
                %defs = ( freq => 0.05,
                          imp => '!Ugly !GEN !INT !UTR' );
            }
            $defs{track} = $cat;
            while (my ($k, $v) = each %defs) {
                &_set_category_keyval($ml, $cat, $k, $v, 1);
            }
        }
        &_set_category_keyval($ml, $uncatCat, 'track', 'Polymorphism', 1);
        while (my ($pop, $cName) = each %catPar) {
            &_set_category_keyval($ml, $pop, 'track', $cName, 1);
        }

        my ($freqHash);
        # Let parent categories pass on values to children:
        foreach my $td (values %{$categorySets}) {
            my $cat = $td->{cat};
            if (my $cName = $catPar{$cat}) {
                # The parent of this category is defined
                if (my $pd = $categorySets->{ $cName } || 
                    $categorySets->{ "$cName Categories" }) {
                    $td->{par} = $cName;
                    # Transfer any parent values to the children
                    while (my ($key, $val) = each %{$pd}) {
                        $td->{$key} = $val unless (defined $td->{$key});
                    }
                }
            } else {
                # No population? Hm.
                $td->{par}   = "Unknown";
                $td->{imp} ||= '!Ugly !NonCoding';
                $td->{freq}  = 0.05;
            }
            if ($td->{nofilter}) {
                # Clear all filters for this category
                $td->{freq}   = 0;
                $td->{imp}    = 'ALL';
                $td->{nullok} = 1;
                $td->{w}      = 10;
                $td->{track}  = $td->{cat};
            }
            unless (defined $td->{w} && $td->{w} ne '') {
                # No weight set
                if ($td->{cat} eq $uncatCat) {
                    $td->{w} = 1;
                } elsif (!$td->{par}) {
                    $td->{par} = 'Unknown';
                    $td->{w} = 2;
                } elsif ($td->{par} eq 'Unknown') {
                    $td->{w} = 2;
                } elsif ($td->{par} eq 'Polymorphism') {
                    $td->{w} = 4;
                } else {
                    $td->{w} = 5;
                }
            }
            my @pkeys = ($td->{tid});
            if (my $pop = $ml->get_population($cat)) {
                # This category is (or is also) a population
                my $pid = $td->{pid} = $pop->pkey();
                push @pkeys, $pid;
            }
            push @pkeys, 0 if ($cat eq $uncatCat);
            if (defined $td->{freq}) {
                # Build frequency hash keyed to object ID
                # This will be used by each_allele() for allele filtering
                $freqHash ||= {};
                foreach my $pk (@pkeys) {
                    $freqHash->{ $pk } = $td->{freq} if
                        (!defined $freqHash->{ $pk } || 
                         $freqHash->{ $pk } < $td->{freq});
                }
            }
        }
        $usedMinFreq = $freqHash if (!defined $usedMinFreq);
        my $catSer = "";
        my $aliTxt = $args->val(qw(impactalias)) || 
            "Ugly UNK COD\nNonCoding GEN INT NCR UTR UT3 UT5";
        if ($aliTxt) {
            my @reqs = ref($aliTxt) ? @{$aliTxt} : ($aliTxt);
            foreach my $req (@reqs) {
                next unless ($req);
                foreach my $line (split(/[\n\r]+/, $req)) {
                    next unless ($line);
                    my @bits = split(/\s+/, uc($line));
                    my $ali  = shift @bits;
                    push @{$impAliases->{$ali}}, @bits;
                }
            }
        }
        my %noSerial = map { $_ => 1 } qw(cat cname par tid);
        foreach my $cat (sort { uc($a) cmp uc($b) } keys %{$categorySets}) {
            $catSer .= "$cat";
            my $td   = $categorySets->{$cat};
            if ($td->{ignore}) {
                delete $categorySets->{$cat};
                next;
            }
            foreach my $k (sort keys %{$td}) {
                next if ($noSerial{$k});
                my $v = $td->{$k};
                next unless (defined $v && $v ne "");
                $v =~ s/\://g;
                $catSer .= " | $k: $v";
            }
            $catSer .= "\n";
            my $imp = $td->{imp};
            if ($imp && $imp !~ /^\s*all\s*$/i) {
                my $ih = $td->{_impHash} = {};
                my $txt = uc($imp);
                foreach my $tt (split(/\s+/, uc($imp))) {
                    my $k = '1';
                    if ($tt =~ /^\!(\S+)/) {
                        ($k, $tt) = (0, $1);
                    }
                    my @toks = ($tt);
                    if (my $ali = $impAliases->{$tt}) {
                        @toks = @{$ali};
                    }
                    foreach my $tok (@toks) {
                        my $id = $ml->impact_details($tok);
                        if (my $tk = $id->{token}) {
                            $ih->{$k}{$tk} = $id->{name};
                        } else {
                            $args->msg_once("[?]", "Unrecognized impact token '$tok'");
                        }
                    }
                }
                &_add_coarse_filter($td);
            }
        }
        $args->set_param('DefaultFilter', $catSer);
        $ml->bench_end("Parse Category Filters");
        #&prebranch($categorySets->{'MGP Polymorphisms'})
    }
    my $rv = exists $categorySets->{$cat} ? 
        $categorySets->{$cat} : { w => 0, track => "" };
    return $rv;
}

sub _stnd_cat_name {
    my $name = shift || "";
    my $lcNm = lc($name);
    if ($lcNm =~ /^ignored?$/) {
        $name = "Ignore";
    } elsif ($lcNm =~ /^discard(ed|ing)?$/) {
        $name = "Discard";
    }
    return $name;
}

sub _set_category_keyval {
    my ($ml, $cat, $kreq, $val, $noReWrite) = @_;
    return unless ($cat);
    my @tids = $ml->text_to_known_pkey( $cat );
    if ($#tids < 1) {
        # Not typo, should be < 1, not < 0
        $args->msg_once("Category filter for unrecognized request '$cat'");
        return;
    }
    unless ($tids[0]) {
        # Normalize the case, since the query did not match
        $cat = $ml->pkey_to_text( $tids[0] = $tids[1] );
    }
    my $key = $catKeyMap->{lc($kreq)};
    unless ($key) {
        $args->msg_once("Unrecognized category filter parameter '$kreq'");
        return;
    }
    my $cd = $categorySets->{ $cat };
    unless ($cd) {
        $cd = $categorySets->{ $cat } ||= {
            cat   => $cat,
            cname => $cat,
            tid   => $tids[0],
            track => "",
            nofilter => "",
        };
    }
    if (defined $cd->{$key} && $cd->{$key} ne "") {
        return if ($noReWrite);
        #$args->msg("[?]","Category filter for '$cat' has redefined:",
        #           "$key:'$cd->{$key}' to $key:'$val'",
        #           "Verify that this is what you desired")
        #    if ($cd->{$key} ne $val);
    }
    # $args->msg("[DEBUG]", "Category $cat : $key = '$val'");
    $cd->{$key} = $val;
}

sub filter_locations {
    my ($lids) = @_;
    my %uLidH  = map { $_ => undef } @{$lids || []};
    my @uLids  = keys %uLidH;
    my $ml     = &maploc();
    $ml->bench_start();
    my $lidTab  = $ml->dbh->list_to_temp_table($lids, 'integer');
    my $newPids = $ml->bulk_MAF_for_location_ids( \@uLids, $bulkFrequencies );
    $ml->bulk_tag_values( \@uLids, 'textToo' ); #$lidTab, 'textToo' );
    $ml->bulk_population_criteria( $newPids );
    # Not using $rRows or $rTab
    my ($imps, $rRows, $rTab, $rids) =
        $ml->bulk_rnas_for_location_ids( \@uLids, undef, $howbad );
    $ml->bulk_features_for_object_ids( $rids, 'textToo');

    my $unkPid = $ml->_unknown_pop_id();
    my (@keep, @rejF, %popFilter, %done);
    foreach my $lid (@uLids) {
        my $isQry = &_is_user_query( $lid );
        if ($isQry && $isQry < 0) {
            $filteredObjects{"variant"}{""}{"User Excluded"}++;
            next;
        }
        next if ($done{$lid}++);
        # Frequency filters get applied on a per-population level, so
        # we test each population independently
        my $ldat = $bulkFrequencies->{$lid};
        my @allPids = keys %{$ldat};
        if ($#allPids == -1) {
            # Location with no population
            push @allPids, $unkPid;
        }

        # We are not considering alleles here, so the loosely-calculated impact
        # will be the same for all populations
        my $predImp = "GEN";
        if (my $impDat = $imps->{$lid}) {
            # TO DO: This should be parameterized for user control...
            # Preferentially use the more conservative RefSeq impact:
            $predImp = $impDat->{imp}{RSR} || $impDat->{imp}{ALL} || $predImp;
        }
        my %freqResults;
        foreach my $pid (@allPids) {
            my $pf = $popFilter{$pid} ||= &_population_filter( $pid );

            # To decide whether we will display this location at all
            # we will consider only the top-weighted population(s).
            # However, if the location is ultimately kept, it may be
            # displayed in multiple tracks. In those cases, we may
            # need to provide different formatting in some tracks
            # (based on MAF) or exclude it all together in some.

            my $w = $pf->weight();
            my $f = $ldat->{$pid};
            my $r = $pf->test_loose($f, $predImp);
            # $args->msg_once("$predImp --> $r") if ($r);
            push @{$freqResults{$w}}, [ $r, $pf, $pid, $f ];
        }
        # &dump_freqResults(\%freqResults) if ($lid == 9027787);
        
        my ($topW) = sort { $b <=> $a } keys %freqResults;
        my %topResults;
        map { $topResults{$_->[0]}{$_->[1]->name()} = 1 } @{$freqResults{$topW}};
        if ($isQry) {
            # User request, always keep
        } elsif (!exists $topResults{""}) {
            # None of the top-weighted populations passed the frequency test
            # We will reject this location
            &_reject_location( \%topResults, \@rejF, $lid );
            next;
        }
        # The location is good!
        push @keep, $lid;

        # We need to note the categories that survived filtration.
        # They will be used once more for one more round of strict
        # impact filtering (if the impact here was 'RNA'), to
        # determine which tracks will display the location, and to
        # adjust MAF-based visualizations on those tracks

        my $smd = $survivingMafData->{$lid} ||= 
            [ [], $predImp eq 'RNA' ? 1 : 0];
        foreach my $pfdats (values %freqResults) {
            foreach my $pfdat (@{$pfdats}) {
                if (!$pfdat->[0]) {
                    my ($pid, $maf, $pf) = ($pfdat->[2], $pfdat->[3], $pfdat->[1]);
                    $maf = -1 unless (defined $maf);
                    push @{$smd->[0]}, [ $pid, $maf, $pf ];
                }
            }
        }
    }
    # warn scalar(@keep)." locations kept, filtered:\n";&prebranch(\%filteredObjects);
    &_note_rejections( \@rejF, "Frequency and loose impact filters" );
    $ml->bulk_accessions_for_object_ids( \@keep );
    $ml->bench_end();
    return \@keep;
}

sub dump_freqResults {
    my $fr = shift;
    my $txt = "";
    foreach my $w (sort { $b <=> $a } keys %{$fr}) {
        $txt .= "WEIGHT $w:\n";
        foreach my $sd (sort { $a->[2] <=> $b->[2] } @{$fr->{$w}}) {
            $txt .= sprintf("   %10d %8s [%s] %s\n", $sd->[2],
                            defined $sd->[3] ? int(0.5 + 1000 * $sd->[3])/10 : '--', $sd->[1]->name(), $sd->[0]);
        }
    }
    &preprint($txt);
}

sub _reject_location {
    my ($results, $rejectArray, $lid) = @_;
    my %rcats;
    while (my ($r, $fnH) = each %{$results}) {
        foreach my $fname (keys %{$fnH}) {
            $rcats{$fname}{$r} = 1;
        }
    }
    while (my ($fname, $rH) = each %rcats) {
        my $r = join(', ', sort keys %{$rH});
        $filteredObjects{"variant"}{$fname}{$r}++;
        push @{$rejectArray}, [$lid, $r, $fname];
    }
}

sub _population_filter {
    my $pid = shift;
    unless ($popFilterObj->{$pid}) {
        my $ml = &maploc();
        $ml->bench_start("Create Population Filter");
        my $pdat = $bulkPopCat->{$pid};
        unless ($pdat) {
            $ml->bulk_population_criteria( [$pid] );
            $pdat = $bulkPopCat->{$pid};
        }
        
        my @cats = @{$pdat->[0]{cat} || []};
        push @cats, @{$pdat->[0]{name} || []};
 
        my %popCatDat;
        foreach my $cat (@cats) {
            if (my $td = &_category_settings( $cat ) ) {
                # Using a hash to avoid duplication
                # Probably not needed
                if ($td->{cat}) {
                    # Ignore the default filter
                    $popCatDat{$td} ||= $td;
                }
            }
        }
        # Sort the relevant categories by user-defined weight
        my @useCat = values %popCatDat;
        @useCat = (&_category_settings($uncatCat)) if ($#useCat == -1);
        $popFilterObj->{$pid} = &_make_filter_object( @useCat );
        $ml->bench_end("Create Population Filter");
    }
    return $popFilterObj->{$pid};
}

sub _filter_for_pid {
    my $pid = shift;
    unless ($cachedPopulationFilter->{$pid}) {
        my $pdat = $bulkPopCat->{$pid};
        unless ($pdat) {
            my $ml = &maploc();
            $ml->bulk_population_criteria( [$pid] );
            $pdat = $bulkPopCat->{$pid};
        }
        my @cats = @{$pdat->[0]{cat} || []};
        push @cats, @{$pdat->[0]{name} || []};
        my %popCatDat;
        foreach my $cat (@cats) {
            if (my $td = &_category_settings( $cat ) ) {
                # Using a hash to avoid duplication
                # Probably not needed
                if ($td->{cat}) {
                    # Ignore the default filter
                    $popCatDat{$td} ||= $td;
                }
            }
        }
        # Sort the relevant categories by user-defined weight
        my @useCat = values %popCatDat;
        @useCat = (&_category_settings($uncatCat)) if ($#useCat == -1);
        $cachedPopulationFilter->{$pid} = &_merge_filters( @useCat );
    }
    return $cachedPopulationFilter->{$pid};
}

sub _user_request_filter {
    unless ($popFilterObj->{-1}) {
        my $ml = &maploc();
        my $pf = $popFilterObj->{-1} = 
            BMS::SnpTracker::MapLoc::PopulationFilter->new
            ( -maploc => $ml );
        $pf->name("Category Filter for User Requests");
        $pf->null_ok( 1 );
        $pf->weight( 99 );
        $pf->categories( [ "User Variant" ] );
    }
    return $popFilterObj->{-1};
}

sub _make_filter_object {

    # This function is called for a specific population (pop_id),
    # called by _filter_for_pid(). It is supplied with one or more
    # hash structures in @useCat. Each structure defines a set of
    # filter criteria that determine if allelic data for the
    # population at a specific location should be included in a query
    # result.

    # The structures are either defined by the defaults or set by the
    # user, and are keyed either to specific populations or to broad
    # categories - of which a population may belong to zero or
    # more. So it is possible for one population to have two or more
    # filter structures assigned to it. The code below is designed to
    # federate the potentially conflicting filter dictates to a single
    # set of values.

    my @useCat    = sort { $b->{w} <=> $a->{w} || $a->{cat} cmp $b->{cat} } @_;
    my @consider;
    foreach my $cat (@useCat) {
        # We will only use the top ranked filter set(s).
        last if ($cat->{w} < $useCat[0]{w});
        push @consider, $cat;
    }
    my @names     = map { $_->{cat} || "Unspecified" } @consider;
    my $mergeCats = "PF::".(join("+", @names) || "CodeError!");
    my $pf;
    unless ($pf = $mergedCategories{$mergeCats}) {
        # Need to make a new filter object
        my $ml = &maploc();
        $pf = $mergedCategories{$mergeCats} = 
            BMS::SnpTracker::MapLoc::PopulationFilter->new
            ( -maploc => $ml );
        # The weight will be the largest one from the categories
        $pf->weight( $consider[0]{w} || 0 );
        $pf->name("Category Filter for ".join(' + ', @names));
        my %tracks = map { ($_->{track} || "Unspecified") => 1 } @consider;
        $pf->categories( [ sort keys %tracks ] );
        my ($impNum, %imps, @freqs, $nullOk) = (0);
        foreach my $td (@consider) {
            push @freqs, $td->{freq} if (defined $td->{freq});
            if (my $ih = $td->{_impHash}) {
                $impNum++;
                while (my ($flag, $hash) = each %{$ih}) {
                    while (my ($tok, $name) = each %{$hash}) {
                        push @{$imps{$flag}{$tok}}, $name;
                    }
                }
            }
            $nullOk++ if ($td->{nullok});
        }
        # Take the most generous frequency
        my ($f) = sort { $a <=> $b } @freqs;
        $pf->min_maf( $f );
        # If any of the categories allowed null frequencies, we will allow them
        $pf->null_ok( $nullOk );
        if ($impNum) {
            # Take the most generous intersection of impacts
            $impNum--;
            while (my ($flag, $hash) = each %imps) {
                my %common;
                while (my ($tok, $arr) = each %{$hash}) {
                    if ($flag == 1 || $#{$arr} == $impNum) {

                        # Either this is a keep flag, or ALL filter
                        # structures agreed to exclude this impact
                        $common{$tok} = 1;
                    }
                }
                my @ctok = keys %common;
                if ($flag == 1) {
                    $pf->require_impact( \@ctok );
                } else {
                    $pf->forbid_impact( \@ctok );
                }
            }
        }
        # $ml->preprint($pf->to_text);
    }
    return $pf;
}

sub _merge_filters {
    return $_[0] if ($#_ == 0);
    my @useCat    = sort { $b->{w} <=> $a->{w} || $a->{cat} cmp $b->{cat} } @_;
    my $mergeCats = join("+", map { $_->{cat} || "UNK" } @useCat) || "";
    my $rv;
    unless ($rv = $mergedCategories{$mergeCats}) {
        $maploc->bench_start();
        $rv = $mergedCategories{$mergeCats} = { w => -1, via => $mergeCats };
        my $topWeight = $rv->{w} = $useCat[0]{w};
        my ($impNum, %imps, @freqs) = (0);
        foreach my $td (@useCat) {
            last if ($td->{w} < $topWeight);
            push @freqs, $td->{freq} if (defined $td->{freq});
            if (my $ih = $td->{_impHash}) {
                $impNum++;
                while (my ($flag, $hash) = each %{$ih}) {
                    while (my ($tok, $name) = each %{$hash}) {
                        push @{$imps{$flag}{$tok}}, $name;
                    }
                }
            }
            if (my $cat = $td->{cat}) {
                my @cats = ref($cat) ? keys %{$cat} : ($cat);
                map { $rv->{cat}{$_} = 1 } @cats;
            }
        }
        # Take the most generous frequency
        my ($f) = sort { $a <=> $b } @freqs;
        $rv->{freq} = $f;
        $rv->{cname} = join(',', sort keys %{$rv->{cat}})
            || "Unknown Category";
        if ($impNum) {
            # Take the most generous intersection of impacts
            $impNum--;
            while (my ($flag, $hash) = each %imps) {
                while (my ($tok, $arr) = each %{$hash}) {
                    if ($flag == 1 || $#{$arr} == $impNum) {
                        # Either this is a keep flag, or
                        # All agree to filter on this item
                        $rv->{_impHash}{$flag}{$tok} = $arr->[0];
                    }
                }
            }
        }
        &_add_coarse_filter($rv);
        # &prebranch($rv); &prebranch(\@useCat);
        $maploc->bench_end();
    }
    return $rv;
}

sub _add_coarse_filter {
    my $mf = shift;
    my $ih = $mf->{_impHash};
    return unless ($ih);
    $maploc->bench_start();
    while (my ($flag, $hash) = each %{$ih}) {
        my $targ = $mf->{_impCoarse}{$flag} ||= {};
        while (my ($tok, $name) = each %{$hash}) {
            $targ->{$tok} = $name;
        }
        if ($flag > 0) {
            # Keep rules
            $targ->{SPL} = $maploc->impact_name('SPL')
                if ($targ->{SP3} || $targ->{SP5});

            $targ->{UTR} = $maploc->impact_name('UTR')
                if ($targ->{UT3} || $targ->{UT5});

            my $rnaCount = 0;
            map { $rnaCount++ if ($targ->{$_}) } 
            qw(CPX FRM STP DEL NON SYN COD UT5 UT3 UTR NCR LOC);
            $targ->{RNA} = $maploc->impact_name('RNA') if ($rnaCount);
        } else {
            # Discard rules
            $targ->{SPL} = $maploc->impact_name('SPL')
                if ($targ->{SP3} && $targ->{SP5});
            $targ->{UTR} = $maploc->impact_name('UTR')
                if ($targ->{UT3} && $targ->{UT5});
        }
    }
    $maploc->bench_end();
}

sub _filter_lid_vs_impact {
    my ($lid, $predImp, $td, $noRequire) = @_;
    return () unless ($predImp && $td);

    # It is advantageous to do a quick-and-dirty initial impact
    # calcualtion That breaks the variants into simply GEN, INT, RNA
    # (genomic, intronic, somewhere in an RNA) In such cases, we do
    # not want to exclude variants that lack a particular impact (eg
    # NON, SYN) becuase they have not been computed yet. The coarse
    # hash takes generic impacts into account

    my $impH = $noRequire ? $td->{_impCoarse} : $td->{_impHash};
    if (my $t = $impH->{1}) {
        unless ($t->{$predImp}) {
            # We failed to match a required impact
            my $req  = join(',', sort values %{$t});
            my $r  = sprintf("Impact $predImp does not match %s", $req);
            # warn "$lid : $r" if ($predImp eq 'GEN' && !$noRequire);
            # $args->msg_once( "$predImp = $r [$td->{cname}]");
            return ($td->{cname}, $r);
        }
        # $args->msg_once( "$predImp = $r [$rcat]");
    }
    # Do not apply exclude filters if full impact not yet calculated:
    return () if ($noRequire);
    if (my $t = $impH->{0}) {
        # &prebranch({pred => $predImp, td => $td }) if ($lid == 44030822);
        if (my $req = $t->{$predImp}) {
            # We matched an excluded impact
            my $r  = sprintf("Impact matches %s", $req);
            return ($td->{cname}, $r);
        }
    }
    return ();
}

sub _note_rejections {
    my ($arr, $note) = @_;
    return unless ($rejFile);
    my $num = $#{$arr} + 1;
    return unless ($num);
    if (open( REJ, ">>$rejFile")) {
        print REJ "\n# $note :\n" if ($note);
        foreach my $rd (@{$arr}) {
            $rd->[0] = "loc_id=".$rd->[0];
            print REJ join("\t", @{$rd})."\n";
        }
        close REJ;
        my $msg = "$num rejected locations ";
        $msg .= "($note) " if ($note);
        $msg .= "written to file";
        $args->msg("[#]", $msg );
        $args->msg_once("[OUT]", $rejFile);
    } else {
        $args->err("Failed to write reject file", $rejFile, $!);
    }
}

sub _pick_variant_track {
    my ($cx) = @_;
    my $ml = &maploc();
    $ml->bench_start();
    my @found = { w => -1, track => "Unknown" };
    my @cats  = @{$cx->{cats} || []};
    push @cats, $uncatCat if ($#cats == -1);
    foreach my $cat (@cats) {
        if (my $td = &_category_settings($cat)) {
            push @found, $td;
        }
    }
    @found     = sort { $b->{w} <=> $a->{w} } @found;
    # die $args->branch(\@found) unless ($#found == -1);
    my $slop   = 0;
    my $weight = $found[0]{w};
    my (%tracks, $reason);
    foreach my $td (@found) {
        last if ($td->{w} + $slop < $weight);
        my $track = $td->{track};
        next unless ($track);
        if (my $r = &_filter_track( $cx, $td )) {
            $reason ||= [$track, $r];
            next;
        }
        $tracks{$track} = 1;
    }
    delete $tracks{Ignore};
    my @t = keys %tracks;
    if ($#t == -1) {
        my $key = uc($cx->{pkey} || $cx->{handle} || "");
        if ($userQueries{$key}) {
            push @t, "User Variant";
        } else {
            my ($track, $reas) = @{$reason || []};
            # warn "($track, $reas)";
            # &prebranch({ found => \@found, cats => \@cats }) unless ($reas);
            $reas ||= "No track for ".join(', ', @cats);
            $filteredObjects{"variant"}{$track || ""}{$reas}++;
            $ml->bench_end();
            return undef;
        }
    }
    $ml->bench_end();
    return \@t;
}

sub _filter_track {
    my ($cx, $td) = @_;
    return 0 if &_is_user_query($cx->{loc_id});
    if (my $ih = $td->{_impHash}) {
        my $tok = $cx->{impToken} || "UNK";
        if (my $name = $ih->{0}{$tok}) {
            return "Impact matches $name";
        } elsif (my $req = $ih->{1}) {
            unless ($req->{$tok}) {
                return "Impact $tok not required";
            }
        }
    }
    
    return 0;
}

sub report_toss {
    my $rv = "";
    return $rv if ($rejFile);
    my $ml = &maploc();
    foreach my $cat (sort keys %filteredObjects) {
        my $cH = $filteredObjects{$cat};
        my @types = sort keys %{$cH};
        next if ($#types == -1);
        my $tot = 0;
        my $txt = "";
        foreach my $t (@types) {
            my $rH = $cH->{$t};
            my @bits;
            foreach my $reason (sort { $rH->{$b} <=> $rH->{$a} || $a cmp $b }
                                keys %{$rH}) {
                my $n = $rH->{$reason};
                next unless ($n);
                $reason = $ml->esc_xml($reason);
                if ($n == -1) {
                    $tot++;
                    if ($cat eq 'feature') {
                        push @bits, sprintf("<a href='http://www.google.com/search?q=%s' target='_blank'>%s</a>", $reason, $reason);
                    } elsif ($cat eq 'RNA') {
                        push @bits, sprintf("<a href='$toolUrl?gene=%s&mode=genome' target='_blank'>%s</a>", $reason, $reason);
                    } else {
                        push @bits, $reason;
                    }
                } else {
                    push @bits, "<br  />&rarr; $reason <span class='mono blue'>($n)</span>";
                    $tot += $n;
                }
            }
            my $bnum = $#bits + 1;
            next unless ($bnum);
            $txt .= "<span class='red'>&times; ";
            $txt .= "<b>".$ml->esc_xml($t)."</b>: " if ($t);
            $txt .= "</span>".($bnum <= 10 ? join(', ', @bits) : "$bnum transcripts"). "<br />\n";
        }
        $rv .= sprintf("<span class='desc'>A total of <span class='mono blue'>%d</span> %s%s were excluded:</span><br />\n", $tot, $cat, $tot == 1 ? '' : 's');
        $rv .= $txt;
    }
    %filteredObjects = ();
    %filterCache     = ();
    return $rv;
}

sub _freqs_for_category {
    my ($cx, $cat, $min) = @_;
    my %alleles;
    while (my ($pid, $fd) = each %{$cx->{freqs} || {}}) {
        my $pdat = &popCache( $pid );
        next unless ($pdat->{cat}{$cat});
        while (my ($al, $v) = each %{$fd}) {
            push @{$alleles{$al}}, defined $v->[0] ? $v->[0] : 1;
        }
    }
    my $pass = 0;
    while (my ($al, $vs) = each %alleles) {
        my ($max) = sort {$b <=> $a} @{$vs};
        next if ($max < $min);
        $pass++;
    }
    return $pass;
}

sub _standard_confs {
    my @rv =  
        (
         -refseqconf => {
             name          => 'RefSeq RNAs',
             subtracksMax => 25,
             trackType    => 'rna',
         },
         -rnaconf => {
             name         => 'Transcripts',
             subtracksMax => 50,
             trackType    => 'rna',
         },
         -genomeconf => {
             # hideName => 1,
             name         => 'Genome',
             outline      => '#ff99ff',
         },
         -ensemblconf => {
             name         => 'Ensembl RNAs',
             height       => 7,
             trackType    => 'rna',
             subtracksMax => 25,
         },
         -featureconf => {
             height       => 3,
             name         => 'Protein Domains',
             trackType    => 'feature',
         },
         -queryconf => {
             height       => 3,
             name         => 'Query',
             trackType    => 'information',
         },
         -populationsconf => {
             height => 3,
             name => 'Cell Lines',
             trackType => 'feature',
             noLegend => 1,
         },
         -observationsconf => {
         },
         -tissuesconf => {
             height => 3,
             noLegend => 1,
         },
         );
    my @sort = ('tally', 'feat', 'var', 'exons', 'refseq','tiss','pop','ensembl');
    &_category_settings();
    my $tCats = $tally{CATS};
    my $ml    = &maploc();
    my @tCols = map { $ml->pastel_text_color($_) } @{$tCats};
    my @oCols = map { $ml->rgba_color( $_, 0 ) } @tCols;
    my @lCols = map { $ml->rgba_color( $_, 0.6 ) } @tCols;
    my $tallyConf = {
        height  => 15,
        displayedFeatures => 1,
        totalFeatures => 1,
        trackType => 'information',
        honorType => 1,
        briefLen  => 1,
        names     => $tCats,
    };
    push @rv, ( -tallyconf => {
        type    => "bar",
        width   => 2,
        fill    => \@tCols,
        outline => \@oCols,
        %{$tallyConf},
    }, -windowed_tallyconf => {
        type    => "line",
        setMaxY => 10,
        outline => \@lCols,
        fill    => \@tCols,
        autowidth => 1,
        window    => 1 + 2 * $tallWin,
        %{$tallyConf},
    });
    %tally = ();
    foreach my $td (sort { $a->{w} <=> $b->{w} } values %{$categorySets}) {
        my $track = $td->{track};
        next unless ($track);
        my $key = lc("-${track}conf");
        $key =~ s/\s+/_/g;
        push @rv, ( $key, {
            name      => $track,
            trackType => 'polymorphism',
            bumpSpace => 1,
            honorType => 1,
        });
        unshift @sort, $track;
    }
    push @rv, (-sort => \@sort);
    return @rv;
}

sub _tally_observations {
    my ($cx) = @_;
    my $ml      = &maploc();
    my ($s, $e) = ($cx->{data}[0][0], $cx->{data}[-1][1]);
    my $show    = $ml->comma_number($s, $e, $cx->{w} == 0 ? '^' : '-');
    my $pos     = $cx->{data}[0][0];
    my $okPids  = $cx->{okPids};
    # &prebranch($cx->{freqs} || {});
    while (my ($pid, $fd) = each %{$cx->{freqs} || {}}) {
        # Do not consider this population if it has failed a filter:
        next if ($okPids && !$okPids->{$pid});
        my $pdat = &popCache( $pid );
        foreach my $t (@{$pdat->{tally}}) {
            my $td = $tally{$pos} ||= {
                pos  => $pos,
                imp  => $cx->{impToken},
                cls  => {},
                show => {},
            };
            $td->{cls}{$t}++;
            $td->{show}{$show} = 1;
        }
    }
}

sub _add_tally_track {
    my ($tracks) = @_;
    my @pos   = sort { $a->{pos} <=> $b->{pos} } values %tally;
    return if ($#pos == -1);
    my %clsH  = map { $_ => 1 } map { keys %{$_->{cls}} } @pos;
    my @cls   = sort keys %clsH;
    my @clInd = (0..$#cls);
    $tally{CATS} = \@cls;
    my $talTrk = "Tally";
    my %byPos;
    foreach my $td (@pos) {
        my @data;
        my $num = 0;
        my $pos = $td->{pos};
        # print "<p>$pos</p>" unless (int($pos) == $pos);
        for my $ci (0..$#cls) {
            my $cnt = $td->{cls}{$cls[$ci]} || 0;
            push @data, $cnt || undef;
            $num += $cnt;
            $byPos{$pos}[$ci] += $cnt || 0;
            
        }
        my $show = join('<br />', sort keys %{$td->{show}});
        push @{$tracks->{$talTrk}}, {
            hideName => 1,
            info     => sprintf("%dbp : %d Observation%s",
                                $pos, $num, $num == 1 ? '' : 's'),
            #fill     => [$rgb],
            #outline  => [$rgb],
            id       => $pos,
            offset   => $pos,
            data     => \@data,
            show     => $show,
            impToken => $td->{imp},
        };
    }
    return unless ($tallMin);

    # See if windowed tallies can be made
    # 100 points along an arbitrary binomial distribution:
    my @binom = (1, 0.999999, 0.998571, 0.995717, 0.99145, 0.985784, 0.978742, 0.970352, 0.960647, 0.949668, 0.937457, 0.924064, 0.909541, 0.893948, 0.877345, 0.859797, 0.841371, 0.822139, 0.802171, 0.781542, 0.760327, 0.738602, 0.716442, 0.693924, 0.671122, 0.648111, 0.624962, 0.601748, 0.578536, 0.555393, 0.532382, 0.509564, 0.486996, 0.464732, 0.442821, 0.421311, 0.400244, 0.379659, 0.35959, 0.340068, 0.32112, 0.302769, 0.285034, 0.267931, 0.251471, 0.235663, 0.220512, 0.206021, 0.192187, 0.179007, 0.166476, 0.154584, 0.143321, 0.132673, 0.122627, 0.113167, 0.104275, 0.095932, 0.08812, 0.080818, 0.074006, 0.067662, 0.061765, 0.056294, 0.051227, 0.046543, 0.042221, 0.03824, 0.03458, 0.03122, 0.028143, 0.025328, 0.022759, 0.020418, 0.018288, 0.016355, 0.014602, 0.013017, 0.011585, 0.010294, 0.009132, 0.008088, 0.007152, 0.006314, 0.005565, 0.004897, 0.004303, 0.003774, 0.003305, 0.002889, 0.002522, 0.002198, 0.001912, 0.00166, 0.00144, 0.001246, 0.001077, 0.000929, 0.0008, 0.000688);
    # There is likely a more elegant way of doing this

    my %win;
    my @ap     = keys %byPos;
    my $winWid = 1 + 2 * $tallWin;
    my @dists  = ((0-$tallWin)..$tallWin);
    my $denom  = 0;
    my @mults;
    foreach my $dist (@dists) {
        #my $m = (1 - ($tallDec || 0)) ** abs($dist);
        # $m = 0.1 if ($m < 0.1);
        # Silly semicircle:
        # my $m = sin(pi * ( 1 - abs($dist) / $tallWin) / 2);
        my $bnInd = int($#binom * abs($dist) / $tallWin);
        my $m     = $binom[$bnInd];
        push @mults, $m;
        $denom += $m;
    }
    # &pre_branch({denom => $denom, dists => \@dists, mults => \@mults });
    foreach my $pos (@ap) {
        my $bp = $byPos{$pos};
        for my $d (0..$#dists) {
            my $p = $pos + $dists[$d];
            my $m = $mults[$d];
            map { $win{$p}[$_] += $m * $bp->[$_] } @clInd;
        }
    }

    my $winTrk  = "Windowed Tally";
    my $minNorm = $tallMin; # * $denom;
    @ap = sort { $a <=> $b } keys %win;
    my %sums;
    foreach my $pos (@ap) {
        my $dat = $win{$pos};
        my $num = 0;
        foreach my $ci (@clInd) {
            my $c = $dat->[$ci];
            $num += $c;
            $dat->[$ci] = int(0.5 + 100 * $c ) / 100;
        }
        next if ($num < $minNorm);
        my $targ = $tracks->{$winTrk} ||= [];
        my $data = $targ->[-1];
        if (!$data || ($data->{offset} + $#{$data->{data}[0]} + 1) < $pos) {
            $data = {
                hideName => 1,
                offset   => $pos,
                data     => [],
                sum      => 0,
            };
            push @{$targ}, $data;
        }
        # push @{$data->{data}}, $num; next;
        my $ind = $pos - $data->{offset};
        for my $d (0..$#cls) {
            $data->{data}[$d][$ind] = $dat->[$d];
        }
    }

    my $trk =$tracks->{$winTrk};
    return unless ($trk);
    my $ml = &maploc();
    foreach my $data (@{$trk}) {
        # $data->{window} = $denom;
        my $dat  = $data->{data};
        my $wid  = $data->{width} = $#{$dat->[0]} + 1;
        my $s    = $data->{offset};
        my $e    = $s + $wid - 1;
        my $show = $ml->comma_number($s, $e, '-');
        my $cnts = $data->{counts} = [];
        $data->{show} = $data->{id} = $show;
        $data->{caption} = "${winWid}bp-windowed tally over ${wid}bp";
        # map { $cnts->[0] += $_ } @{$dat};
        for my $i (0..$#cls) {
            my $cdat  = $dat->[$i];
            my $count = 0;
            for my $ind (0..$#{$cdat}) {
                $cdat->[$ind] = $win{ $ind + $s}[$i]
                    unless (defined $cdat->[$ind]);
                $count += $cdat->[$ind];
            }
            $dat->[$i] = [] unless ($count);
            $cnts->[$i] = int(0.5 + 100 * $count /($denom))/100;
        }
    }
}

sub _variant_extras {
    my ($cx) = @_;
    if ($noAnonFeat) {
        if (my $feats = $cx->{features}) {
            my @ids = keys %{$feats};
            foreach my $id (@ids) {
                my $feat = $feats->{$id};
                delete $feats->{$id} if ($id =~ /^[A-Z]{2,6}\:\d+$/ && 
                                         ! $feat->{name});
            }
        }
    }
    
}

sub usage {
    my $rv = "";
    my @mbits;
    my $filt = &report_toss();
    if ((!$nocgi || $fullHTML) && $filt) {
        $rv .= $filt;
        push @mbits, "filtered locations";
    }
    unless ($nocgi || $stuff{UsageDone}++) {
        $rv .= "<div style='border: blue solid 1px; font-size:0.8em; width:30em;'>".join
            ("<br />", "<span style='font-style:italic; font-weight: bold;'>Usage Instructions</span>",
             "Think 'Google Maps': <b>Click-and-drag</b> to pan left/right, <b>mouse wheel</b> to zoom in/out. <b>Shift-click-drag</b> to draw a box, which will expand when mouse released. <b>Hover</b> over feature for details, <b>click</b> to 'pin' details. <b>Escape key</b> resets view.",
             &help('SoftwareOverview','[FullUsage Instructions]')). "</div>";
        push @mbits, "usage instructions";
    }
    $rv = "<div class='toggle' onclick='toggle_hide(this)'>Click for ".join
        (' and ', @mbits)."</div><div class='hide'>$rv</div>" if ($rv);
    return $rv;
}

sub popCache {
    my $pid = shift;
    unless ($popcache{$pid}) {
        my $ml = &maploc();
        $ml->bench_start('Get Population');
        my $pop = $ml->get_population($pid);
        $pop->read();

        my $name = $pop->name();
        my $pdat = $popcache{$pid} = {
            name => $name,
            pop  => $pop,
            cat  => {},
        };
        foreach my $cat ($pop->tag_values("Category")) {
            $pdat->{cat}{$cat} = 1;
        }
        my @tallied = $pop->tag_values("Tally");
        if ($#tallied == -1 && $tallyAll) {
            my ($ctag) = $pop->colored_tag();
            push @tallied, $ctag;
        }
        $pdat->{tally} = \@tallied;
        $ml->bench_end('Get Population');
    }
    return $popcache{$pid};
}

sub HTML_PERMALINK {
    return unless ($htmlEnv);
    my $runfile = &temp_file("Parameters", 'param');
    if (open(RUNFILE, ">$runfile")) {
        map { $args->blockquote($_, 1) } qw(catfilt);
        print RUNFILE $args->to_text
            ( -ignore => [@protectedArgs, 'pause',  'paramfile', 'valuefile'],
              -blockquote => ['DefaultFilter']);
        close RUNFILE;
        chmod(0666, $runfile);
        print $fh "<span style='color: gray ! important; font-size:0.6em;'>Links to this analysis - may be copied for email: <a href='$toolUrl?valuefile=$runfile'>Run</a> or <a href='$toolUrl?pause=1&valuefile=$runfile'>Allow Changes</a>. You can also <a target='_blank' href='$toolUrl?setDefault=1&valuefile=$runfile'>make these settings your default</a>.</span><br />";
    } else {
        $args->err("Failed to create parameter file", $runfile, $!);
    }
}

sub temp_file {
    my ($tag, $sfx) = @_;
    $tag ||= "Data";
    $sfx ||= "unk";
    unless ($tmpPrfx) {
        $tmpPrfx = "/tmp/MapLoc/MapLoc-$$-".time;
        $args->assure_dir($tmpPrfx, 1);
    }
    my $counter = -1;
    my $path;
    do {
        $counter++;
        $path = sprintf("%s-%s%s.%s", $tmpPrfx, $tag, $counter || "", $sfx);
    } while ($tmpFiles{$path}++);
    return $path;
}

sub ref_diff_key {
    my $html = "";
    return $html if (!$drawTab || $noHTML);
    $html .= "<table class='tab freqtab'>\n";
    $html .= " <caption>Color Key: Difference relative to normal control</caption><tbody>\n";
    $html .= "<tr>\n";
    for my $rd (0..10) {
        my $cl = sprintf(" norm%d", $rd);
        $html .= sprintf(" <td class='norm%s' title='%d%% relative to control'>%s/100</td>",
                         $cl, 10 * $rd, 10 * $rd);
    }
    $html .= "</tr>\n";
    $html .= "</tbody></table>\n";
    return $html;
}

sub prebranch {
    print $args->branch( @_ );
}

sub preprint {
    my $txt = join("\n", map { defined $_ ? $_ : "" } @_);
    return "" unless ($txt);
    if ($nocgi) {
        warn "$txt\n";
    } else {
        $txt = join("\n", map { $args->esc_xml($_) } map { defined $_ ? $_ : "" } @_);
        print $fh "<pre style='border: solid blue 1px; margin:3px; background-color:#eee;'>$txt</pre>";
    }
    return "";
}

sub collapse_ranges {
    my $rangeReq = shift;
    # Nucleate with the first HSP:
    my @sorted = sort { $a->[0] <=> $b->[0] } @{$rangeReq};
    my @rv = ( [ @{shift @sorted} ] );
    foreach my $r (@sorted) {
        if ($r->[0] <= $rv[-1][1]) {
            # Overlap - extend the HSP if needed
            $rv[-1][1] = $r->[1] if ( $rv[-1][1] < $r->[1]);
        } else {
            # Non overlap, 
            push @rv, $r;
        }
    }
    return wantarray ? @rv : \@rv;
}

sub split_req {
    my $text = shift || "";
    my @reqs;
    foreach my $line (split(/[\n\r\,]+/, $text)) {
        $line =~ s/^\s+//;
        $line =~ s/\s+$//;
        next if ($line =~ /^\s*$/);
        push @reqs, [split(/\s+/, $line)];
    }
    return @reqs;
}

sub HTML_START {
    return unless ($htmlEnv);
    my $bxtra = "";
    my $url = $args->val(qw(cxurl)) || "";
    my $hxtra = "";
    if ($debug) {
        $url =~ s/\.min\./\.debug\./;
        $hxtra .= "  <script type='text/javascript' src='http://xpress.pri.bms.com/JAVASCRIPT/canvas/js/canvasXpress.public.min.js'></script>\n";
    }
    if ($format =~ /(xpress|canvas|cx|html)/i) {
        if ($url) {
            my $ieKludge = $url;
            $ieKludge =~ s/\/[^\/]+$//;
            $hxtra .= <<EOF;
  <!--[if IE]>
    <script type='text/javascript' src='$ieKludge/flashcanvas.js'></script>
  <![endif]-->
  <script type='text/javascript' src='$url'></script>
EOF

        } else {
            $bxtra .= "<p class='err'>Could not find URL for canvasXpress</p>\n";
        }
    }
    my $ml   = &maploc();
    my ($mlH, $mlB) = $ml->cx_support_links();
    $bxtra .= $mlB;
    $hxtra .= $mlH;
    if (my $ico = $args->val(qw(favicon))) {
        $hxtra .= sprintf('<link rel="shortcut icon" href="%s">', $ico);
    }

    print $fh <<EOF;
<!DOCTYPE HTML>
<html>
 <head>
  <title>MapLoc Data Explorer</title>
  <meta http-equiv="X-UA-Compatible" content="chrome=1">
$hxtra </head>
 <body>$bxtra
EOF

}

sub help { return $args->tiddly_link( @_); }

sub tw_hash {
    unless ($twHash) {

        my @twTopics = qw
            (FindGenes OnlyExons FindCellLines FindSNPs OnlyHigh
             RnaMatchQuality OutputFormat BuildToken VariantTable RnaView
             FilterEnsembl ValidationStringency CategoryFilter PopulationFilter
             VisualMode TrackConfiguration CustomColors ExcludeEnsembl
             TallyAll QueryLimit RnaRange UserQuery LocationAssignment
             DebugLevel SplitFeatures BenchMarks MaxFeatLen NoAnonFeat
             VariantSource VariantTrack VariantPriority MinimumMAF
             NullOk ImpactFilter ExonOnly AllowedEnsembl AutoMode UserOnly
             HowBad ImpactFrequency NoFilter ForceBuild ShowPopulation
             AdvancedSettings
             );
        push @twTopics , map { $_->[1] } values %{$cfColDat};
        $twHash = {};
        map { $twHash->{$_} ||= &help($_) } @twTopics;
    }
    return $twHash;
}

sub HTML_FORM {
    return unless ($htmlEnv);
    my $ml = &maploc();
    $ml->bench_start();
    my $impExamp = "";

    my %chk = ( exononly    => $exonOnly   ? "CHECKED " : "",
                drawtab     => $drawTab    ? "CHECKED " : "",
                splitfeat   => $splitFeat  ? "CHECKED " : "",
                showbench   => $showBench  ? "CHECKED " : "",
                onlyhigh    => $onlyHigh   ? "CHECKED " : "",
                noanonfeat  => $noAnonFeat ? "CHECKED " : "",
                tallyall    => $tallyAll   ? "CHECKED " : "",
                forcebld    => $forceBld   ? "CHECKED " : "",
                );

    my $tw = &tw_hash();

    &_category_settings();
    my %cpri = ( 'Polymorphism' => 10, 'Mutation' => 20, 'Region' => 30 );
    my @cats = sort { ($cpri{$a->{par} || ""} || 100) <=>
                          ($cpri{$b->{par} || ""} || 100) || 
                          ($b->{w} || 0) <=> ($a->{w} || 0) ||
                      lc($a->{cat} || 'zzz') cmp lc($b->{cat} || 'zzz')
                  } values %{$categorySets};
    my @catTab;
    my %isPar = map { $_ => 1 } @catParents;
    foreach my $td (@cats) {
        my $pop   = $td->{cat}; # Full name of category
        next if ($pop =~ /Categories/);
        my $cName = $td->{par} || ""; # High-level category parent
        my $lab   = $pop;       # Pretty name of category
        $lab      =~ s/\s*${cName}s?\s*/ /gi if ($cName);
        $lab      =~ s/\s+/ /g; $lab =~ s/\s+$//; 
        next if ($pop eq $cName || $lab eq $cName);
        my $slab  = $lab;
        $slab .= " <i>(others)</i>" if ($isPar{$pop} && $cName eq 'Unknown');
        my ($defW, $defT) = ($td->{w} || 1, $td->{track} || "Ignore");
        my %row = ( cat => $slab);
        my @options = ( $cName, $lab, "Ignore", "Discard");
        # Track Selector:
        my $trSel      = "<select name='catfilter'>\n";
        foreach my $track (@options) {
            $trSel .= " <option value='$pop | track:$track'";
            if (lc($track) eq lc($defT)) {
                $trSel .= " selected='SELECTED'";
                $defT = "";
            }
            $trSel .= ">$track</option>\n";
        }
        $trSel .= " <option value='$pop | track:$defT' selected='SELECTED'>$defT</option>\n" if ($defT);
        $trSel .= "</select>";
        $row{track} =  $trSel;

        # Weight Selector:
        my $wSel = "<select name='catfilter'>\n";
        foreach my $w (1..9) {
            $wSel .= " <option value='$pop | w:$w'";
            if (lc($w) eq lc($defW)) {
                $wSel .= " selected='SELECTED'";
                $defW = "";
            }
            $wSel .= ">$w</option>\n";
        }
        $wSel .= " <option value='$pop | w:$defW' selected='SELECTED'>$defW</option>\n" if ($defW);
        $wSel .= "</select>";
        $row{w} = $wSel;

        # Frequency filter
        my $ff = defined $td->{freq} ? $td->{freq} : "";
        $row{freq} = "<input type='text' style='width:3em;' name='catfilt_${pop}_freq' value='$ff' />";

        # Null Ok flag
        $row{nullok} = "<input type='checkbox' name='catfilt_${pop}_nullok' value='1' ".($td->{nullok} ? 'CHECKED' : '')."/>";

        $row{nofilter} = "<input type='checkbox' name='catfilt_${pop}_nofilter' value='1' ".($td->{nofilter} ? 'CHECKED' : '')."/>";

        # Impact filter
        my $imf = defined $td->{imp} ? $td->{imp} : "";
        $row{imp} = "<input type='text' style='width:15em;' name='catfilt_${pop}_imp' value='$imf' />";

        push @catTab, \%row;
    }
    my @colOrder;
    while (my ($tok, $cfd) = each %{$cfColDat}) {
        if (my $cn = $cfd->[2]) {
            $colOrder[$cn-1] = $tok;
        }
    }
    
    my $catHTML = "<table class='tab'><tbody>\n";
    $catHTML   .= " <tr>";
    foreach my $tok (@colOrder) {
        my $cfd = $cfColDat->{$tok};
        my $nm  = $cfd->[0];
        $nm = '&infin;' if ($tok eq 'nofilter');
        $catHTML .= sprintf("<th style='align:center' title='%s'>%s<br />%s</th>",
                            $args->esc_xml_attr($cfd->[4]), $nm, 
                            $tw->{$cfd->[1]});
    }
    $catHTML .= "</tr>\n";
    foreach my $cr (@catTab) {
        my @row;
        for my $c (0..$#colOrder) {
            my $tag = $c ? 'td' : 'th';
            # warn $colOrder[$c]."\n";
            push @row, "<$tag>".$cr->{$colOrder[$c]}."</$tag>";
        }
        $catHTML .= " <tr>".join('', @row)."</tr>\n";
    }
    $catHTML .= "</tbody></table>\n";
    foreach my $ali (sort keys %{$impAliases || {}}) {
        $catHTML .= "<b>$ali</b> = ".join(" ", map {
            $ml->impact_html_token($_)
        } @{$impAliases->{$ali}}). "<br />\n";
    }

    my $defFilter = $args->val(qw(defaultfilter)) || "";
    my $selHTML = "$tw->{OutputFormat}<span class='param'>Format:</span> ".
        "<select name='format'>\n";
    foreach my $f ('Text', 'HTML', 'CanvasXpress', 'Excel', 'Debug') {
        $selHTML .= sprintf("  <option value='%s'%s>%s</option>\n", $f, 
                            (lc($f) eq lc($format)) ? " selected='SELECTED'" : "", $f);
    }
    $selHTML .= "</select><br />\n";
    # Overriding - putting this all into mode
    $selHTML = "";


    my $modeHTML = "$tw->{VisualMode}<span class='param'>Report Format:</span> ".
        "<select name='mode'>\n";
    my $chosenMode = "";
    foreach my $md ('Auto', 'Genome Browser', 'Locus Browser', 'RNA Browser',
                    'Excel Report', 'Coordinate Excel',
                    'Gene-Sample Network', 'Text Report') {
        my $sel = "";
        my $chkMd = &_stnd_mode($md);
        if ($chkMd eq $mode) {
            my $chkFmt = &_stnd_format($md);
            if ($chkFmt eq $format) {
                $sel = " SELECTED";
                $chosenMode = $md;
            }
        }
        $modeHTML .= sprintf
            ("  <option value='%s'%s>%s</option>\n", $md, $sel, $md);
    }
    unless ($chosenMode) {
        # The user choice was not one of the standard ones
        $modeHTML .= sprintf
            ("  <option value='%s'%s>%s</option>\n", $mode," SELECTED", $mode);
    }
    $modeHTML .= "</select><br />\n";

    my $buildHTML = "$tw->{BuildToken}<span class='param'>For non-accessioned queries use build:</span> ".
        "<select name='build'>\n";
    foreach my $bld ('', $maploc->all_builds()) {
        $buildHTML .= sprintf
            ("  <option value='%s'%s>%s</option>\n", $bld, 
             ($buildReq eq $bld) ? " SELECTED" : "", $bld || "All");
    }
    $buildHTML .= "</select><br />\n";
    $buildHTML .= "$tw->{ForceBuild}<input type='checkbox' name='forcebld' value='1' $chk{forcebld}/>Also force the above build for all queries<br />";

    my $dbHTML =  "<span class='param'>Debug:</span> ".
        "<select name='dumpsql'>\n";
    foreach my $ds ([0,'Quiet'], [1,'Show SQL'], [2, 'Explain SQL']) {
        my ($v, $n) = @{$ds};
        $dbHTML .= sprintf
            ("  <option value='%s'%s>%s</option>\n", $v, 
             ($v eq $dumpSQL) ? " SELECTED" : "", $n);
    }
    $dbHTML .= "</select><br />\n";

    my $qtxt = $query;
    if (ref($qtxt) && ref($qtxt) eq 'ARRAY') {
        $qtxt = join("\n", @{$qtxt});
    }
    my $showPopText = join("\n", $args->each_split_val( 'showpop' )) || "";
    
    print $fh <<EOF;
  <form method='post'>
      <table><tbody><tr><td>
      <b>$tw->{UserQuery}Your Query</b> <span class='ifaceNote'>One per line; SNPs, Genes, RNAs, etc</span><br />
       <textarea spellcheck='false' rows='10' cols='40' name='query' style='background-color:#afa'>$qtxt</textarea><br />
</td><td style='vertical-align:top;'>
     $selHTML
     $modeHTML
     <input type='submit' style='background-color:lime; font-weight: bold; font-size: 1.5em;' value='Search' />
</td></tr></tbody></table>

   <div class='toggle' onclick='toggle_hide(this)'>$tw->{AdvancedSettings} Click for advanced settings</div><span class='hide spanDiv' style='background-color:#eef;'>
      $tw->{CategoryFilter}<i>Choose how to organize variants, and relative priority to give to different categories:</i>
$catHTML

   $tw->{CustomColors}<span class='param'>Custom Colors</span> <span class='note'>eg: "phosphorylation site DarkOliveGreen"</span> <a class='toggle' target='_blank' href='http://en.wikipedia.org/wiki/Web_colors#X11_color_names'>Valid Colors</a><br />
    <textarea spellcheck='false' rows='5' cols='70' name='customcolors' wrap='off' style='background-color:#ff9'>$custColor</textarea><br />

   $buildHTML

   <input type='hidden' name='splitfeat' value='0' />
   $tw->{SplitFeatures}<input type='checkbox' name='splitfeat' value='1' $chk{splitfeat}/> Break large features into pieces<br />

   <input type='hidden' name='noanonfeat' value='0' /> 
   $tw->{NoAnonFeat}<input type='checkbox' name='noanonfeat' value='1' $chk{noanonfeat}/> Ignore features without descriptions<br />
    $tw->{MaxFeatLen}Disregard features longer than  
   <input style='width:2em' type='text' name='maxfeatlen' value='$maxFeatLen' />%
    of their hosting RNA.<br />
    $tw->{QueryLimit}Recover at most  
   <input style='width:4em' type='text' name='limit' value='$limit' />
    results for any given query.<br />

   $tw->{RnaMatchQuality}Consider RNA matches up to 
   <input type='text' name='howbad' value='$howbad' style='width:2em' />%
    worse than best genome match.<br />
   $tw->{RnaRange}Extend search by 
   <input type='text' name='range' value='$range' style='width:5em' /> bp
    around RNA/Gene queries.<br />

   $tw->{ExonOnly}<input type='checkbox' name='exononly' value='1' $chk{exononly}/>
   RNA/Gene queries should search the genome only around exons (will still be extended by above value).<br />

   $tw->{AllowedEnsembl}Only consider the following Ensembl status / biotypes:
   <br /><textarea spellcheck='false' rows='3' cols='40' name='allowens' wrap='off' style='background-color:#0ff'>$allowEnsText</textarea><br />
   <input type='hidden' name='allowens' value='0' /> 
   
   $tw->{ShowPopulation}In Excel Report, add MAF columns for these populations:
   <br /><textarea spellcheck='false' rows='3' cols='40' name='showpop' wrap='off' style='background-color:#ff0'>$showPopText</textarea><br />

   $tw->{TallyAll}<input type='checkbox' name='tallyall' value='1' $chk{tallyall}/> Force tallies for all SNP categories<br />

   <i>Nerd Settings</i><br />
   $tw->{DebugLevel}$dbHTML
   $tw->{BenchMarks}<input type='checkbox' name='showbench' value='1' $chk{showbench}/>Show program timing benchmarks<br />

     <input type='submit' style='background-color:lime; font-weight: bold; font-size: 1.5em;' value='Search' />
   </span>

  </form>
EOF

    $ml->bench_end();
}

sub HTML_END {
    return unless ($htmlEnv);
    print $fh <<EOF;
 </body>
</html>
EOF
}

sub maploc {
    return $maploc if ($maploc);
    $maploc = BMS::SnpTracker::MapLoc->new
        ( -build    => $args->val(qw(build)),
          -noenv    => $args->val(qw(noenvironment noenv)),
          -instance => $mlInt, );
    if ($nocgi) {
        # For statically created files, we probably need to make URLs
        # absolute.
        
        my $wsu = $maploc->web_support_url();
        $maploc->web_support_url("$toolUrlDir/$wsu");
    }

    foreach my $assure ($uncatCat, @catParents) { 
        $maploc->text_to_pkey($assure);
    }
    my %cx_global = ( debug => $cxDebug );
    if (my $idir = $args->val(qw(imagedir))) {
        $cx_global{imageDir} = $idir;
    }
    $maploc->cx_panel_settings( \%cx_global );

    $maploc->verbosity( $vb || 0);
    $maploc->default_howbad( $howbad );
    $maploc->default_loc_to_rna_distance( $range );
    # $maploc->population_js_extender( \&_extend_pop_json );
    $custColor = "";
    my @errs;
    foreach my $cc ($args->each_split_val(@ccParam)) {
        next unless $cc;
        $custColor .= "$cc\n";
        if ($cc =~ /^\s*(.+?)\s+(\S+)\s*$/) {
            my ($text, $col) = ($1, $2);
            my $hex = $maploc->hex_color( $col );
            unless ($hex) {
                push @errs, "'$text' = $col";
                next;
            }
            $maploc->user_text_color($text, $hex);
        }
    }
    map { $args->clear_param( $_ ) } @ccParam;
    $args->set_param('CustomColors', $custColor);
    $args->msg_once("[!]","Some color assignments could not be recognized",
                    sort @errs) unless ($#errs == -1);
    $bulkPopCat = $maploc->bulk_cache_source( 'popcat' );
    if ($args->val(qw(nodbsnpundef noundefdbsnp))) {
        if (my $ppk = $maploc->population_name_to_pkey('Unvalidated dbSNP')) {
            $hidePop ||= [];
            push @{$hidePop}, $ppk;
        }
    }
    return $maploc;
}

sub reporter {
    unless ($srpt) {
        $srpt = BMS::SnpTracker::MapLoc::Reporter->new( &maploc );
        $srpt->bulk_frequencies( $bulkFrequencies );
        $srpt->base_url($toolUrl);
        foreach my $typ (qw(gene rna variant)) {
            if (my $st = $args->val($typ.'showtag')) {
                foreach my $tag (split(/\s*[\n\r]+\s*/, $st)) {
                    $srpt->enumerate_tag( $typ, $tag );
                }
            }
        }
        foreach my $req ($args->each_split_val( 'showpop' )) {
            $srpt->add_population_column( $req );
        }
    }
    return $srpt;
}

sub eh {
    my $path  = shift;
    my $type  = shift || "Excel Report";
    my $isNew = 0;
    my $eh;
    unless ($eh = $excel{$type}) {
        my $srpt = &reporter();
        $path  ||= sprintf("%s/%d-%d.xlsx", $htmlTmp, time, $$);
        $args->assure_dir($path, 'isFile');
        $isNew   = 1;
        my $url  = $args->path2url($path);
        $eh      = $excel{$type} ||= $srpt->
            excel_helper( -path => $path, -addgo => 1 );
        $eh->url($url);
    }
    return wantarray ? ($eh, $isNew) : $eh;
}

sub pre_branch {
    print "<pre>".$args->branch(@_)."</pre>";
}
