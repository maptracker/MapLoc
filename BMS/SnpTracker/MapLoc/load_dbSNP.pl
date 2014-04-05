#!/usr/bin/perl -w


=head1 DESCRIPTION

=head1 SYNOPSIS

=head1 Usage

 Some typical command line options:

  Test parse chromosome 1:
  load_dbSNP_XML.pl -chr 1 -limit 10

  Load all chromosomes into the database:
  load_dbSNP_XML.pl -doall

=head1 Options

       -dir Path to the directory containing the NCBI files. Will be
            extracted from the wget script if it is not explicitly given.

  -progress Default 180. Frequency, in seconds, that progress
            information is reported.

      -fork Default 1. Then number of children to assign to
            processing.

     -error Default 'dbSNP_Parse_Errors'. Similar to -prog, but
            contains information on every error encountered during the
            parse. Very useful to see what problems are encountered
            with the data set. Almost all of the errors will deal with
            inconsistencies in the XML file.

       -snp Default undef. Use to specify a specific SNP to load
            (useful for testing purposes).

       -chr Use this option to specify a particular chromosome that
            you want parsed. Examples include 1, X, NotOn.

   -skipchr Default undef. If defined, will be interpreted as a list
            of chromosomes that should NOT be analyzed. Really only
            useful when using -doall.

     -doall If true, then ALL chromosomes found in -dir will be
            parsed. The program will automatically spawn a weightless
            (nice -n 19) -chr process for each file it finds.

            If neither -doall nor -chr are specified, this help will
            be shown.

     -limit Default 0. If zero or undef, all records in the XML file
            will be parsed. If set, then only that number will be
            analyzed.

     -clear Default undef. If true, then all previously generated
            progress and error files will be deleted. If set to 'full'
            then SnpTracker load flags will be cleared. The load flags
            are used by SnpTracker to remember entries it previously
            parsed - if you want to reparse all entries, you should
            -clear full.

   -verbose Default 1. If true, then the program will provide detailed
            information while running.

     -cache Default 5000. Number of MapTracker rows to gather before
            pushing them into the load queue.

  -fulltest Default 0. If true, then show test information even for
            SNPs that are recorded as being processed.

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

use strict;

use BMS::SnpTracker::MapLoc;

use BMS::Utilities::SequenceUtilities;
use BMS::Utilities::FileUtilities;
use BMS::Utilities::Escape;

use BMS::ArgumentParser;
use BMS::ForkCritter;
use BMS::FriendlySAX;
use BMS::TableReader;

my $popFmt  = "Pop:%d";
my $stAuth  = "dbSNP";
my $defPrfx = "dbSNP_";
my $args = BMS::ArgumentParser->new
    ( -nocgi          => $ENV{'HTTP_HOST'} ? 0 : 1,
      -error          => "${defPrfx}Parse_Errors.txt",
      -note           => "${defPrfx}Notes.txt",
      -verb           => "${defPrfx}Verbose_Output.txt",
#      -dir            => "/work5/tilfordc/dbsnp",
      -pops           => 'poplist.xml.gz',
      -subs           => 'database/shared_data/Submitter.bcp.gz',
      -chr            => "",
      -doall          => 0,
      -cache          => 5000,
      -simpledir      => '',
      -clear          => '',
      -fork           => 1,
      -verbose        => 1,
      -progress       => 300,
      -species        => 'human',
      -paramfile      => 'BMS::SnpTracker::MapLoc',
      -wget           => '/work5/tilfordc/WGET/dbsnp_cmd.sh',
      );

$args->shell_colors();
$args->debug->skip_key([qw(PARENT)]);

my ($fc, $mapLoc, $locCategory, $locCatName, $rootPop,
    $nowSnp, %simpPops, %noted, $mlPops, $newPops, $anonPopId,
    $specName, $globalBuild, $textMaps, $cache,
    $cfrm, %statFiles, %errWarnCount);
my $ewcLvl = 50;
my $seqBlock = 80;

my $joiner    = "|";
my $isTrial   = $args->val(qw(trial istrial));
my $vb        = $isTrial || $args->val(qw(vb verbose)) || 0;
my $forkNum   = $args->val(qw(forknum fork)) || 1;
my $dbInst    = $args->val(qw(instance));
my $cname     = "";
my $starttime = time;
my $counter   = 0;
my $dir       = $args->{DIR};
my $simpDir   = $args->val(qw(simpdir simpledir));
my $su        = BMS::Utilities::SequenceUtilities->new();
my $wcmd      = $args->val(qw(wgetcmd wget)) || "";
my $limit     = $args->val(qw(limit)) || 0;
my $progress  = $args->val(qw(progress)) || 0;
my $esc       = BMS::Utilities::Escape->new();
my $mode      = lc($args->val(qw(mode)) || "");
my $doDump    = $args->val(qw(dodump dump)) || 0;

my $noMaintenance = $args->val(qw(nomaint nomaintain nomaintenance)) || 0;
my $minSampleSize = 10;

unless ($dir) {
    # Recover the directory from the WGET update script
    my @dirBits   = split(/[\n\r]+/, `grep '^TARGDIR=' $wcmd`);
    my $db = $dirBits[-1] || "";
    if ($db =~ /TARGDIR=\'([^\']+)\'/) {
        $dir = $1;
        $dir =~ s/\/+$//;
        $dir .= "/ftp.ncbi.nih.gov/snp";
    } else {
        $args->death("Failed to find the data directory from wget command",
                     $wcmd);
    }
}
$dir =~ s/\/+$//;
unless ($simpDir) {
    $simpDir = "$dir/dbSnpSimple";
}
$args->msg("[-]","dbSNP data directory:", $dir, "Simplified XML:", $simpDir);
my $orgdir    = "$dir/organisms";

if (my $upd = $args->val(qw(doupdate))) {
    $args->msg("Mirroring via shell script:", $wcmd, `date`) if ($vb);
    system($wcmd);
    $args->msg("Done ",`date`) if ($vb);
    if ($upd =~ /only/i) {
        $args->msg("Only running update");
        exit;
    }
}

foreach my $key (qw(ERROR VERB NOTE)) {
    if (my $file = $args->val($key)) {
        $statFiles{$key} = $file;
    }
}

&text_maps();

my $dsfrm = '%s/organisms/%s/XML/ds_ch%s.xml.gz';

foreach my $modeArg (qw(loadpop loadgt)) {
    last if ($mode);
    $mode = $modeArg if ($args->val($modeArg));
}

if ($mode =~ /purge/) {
    &purge_alleles();
} elsif ($mode =~ /(simp|load)/) {
    my ($meth, $what);
    if ($mode =~ /load/) {
        if ($mode =~ /(pop)/) {
            ($meth, $what) = (\&load_populations, "Load Populations");
        } elsif ($mode =~ /(geno|gt)/) {
            ($meth, $what) = (\&load_simple_genotypes, "Load Genotypes");
        } elsif ($mode =~ /(asn|pos)/) {
            ($meth, $what) = (\&load_simple_asn_files, "Load ASN Aliases");
        } else {
            $args->death("NEED $mode");
        }
    } elsif ($mode =~ /(pos|asn)/) {
        ($meth, $what) = (\&SIMPLIFY_ASN, "Simplify ASN1");
    } elsif ($mode =~ /(geno|gt)/) {
        ($meth, $what) = (\&SIMPLIFY_GT, "Simplify Genotype XML");
    } else {
        ($meth, $what) = (\&SIMPLIFY, "Simplify Basic XML");
    }
    $args->death("Unknown mode '$mode'") unless ($meth);
    $args->msg("Mode : $what");
    &{$meth};
} else {
    $args->death("Add usage help here");
}

sub load_simple_asn_files {
    my $dDir = $simpDir;
    my $files = &find_files($dDir, '^ds_flat', '(.+)\.tsv\.gz$');
    map { &alter_statfile( $_, 'LoadASN') } ('ERROR','NOTE');
    foreach my $file (@{$files}) {
        $args->msg("[LOAD]", $file);
        $fc = BMS::ForkCritter->new
            ( -inputtype   => 'tsv header hash',
              -init_meth   => \&init_maploc,
              -method      => \&load_simple_asn,
              -finish_meth => \&finalize,
              # -last_item_meth => \&get_simple_item,
              -limit       => $limit,
              -progress    => $progress,
              -verbose     => $vb );
        $fc->reset();
        $fc->input($file);
        $args->msg("[+]","Parsing $file", `date`);
        &run_fork();
    }
    &report_files();
}

sub load_simple_genotypes {
    my $dDir = $simpDir;
    my $files = &find_files($dDir, '^gt_');
    map { &alter_statfile( $_, 'LoadGeno') } ('ERROR','NOTE');
    foreach my $file (@{$files}) {
        my $bld = &extract_build( $file );
        &load_gt_populations($file);
        $args->msg("[LOAD]", $file);
        $fc = BMS::ForkCritter->new
            ( -inputtype   => 'sax',
              -inputargs   => [ -tag => ['rs'] ],
              -init_meth   => \&init_maploc,
              -method      => \&load_simple_gt,
              -finish_meth => \&finalize,
              -last_item_meth => \&get_simple_item,
              -limit       => $limit,
              -progress    => $progress,
              -verbose     => $vb );
        $fc->reset();
        $fc->input($file);
        $args->msg("[+]","Parsing $file", `date`);
        &run_fork();
    }
    &report_files();
}

sub load_gt_populations {
    my $file = shift;
    # Population data reside at the start of the XML file in <popDat> tags,
    # prior to the main mass of <rs> entries. The load_popdat() method runs
    # through the popDat entries then aborts when it hits the first rs.
    &init_maploc();
    $newPops = 0;
    eval {
        my $fs = BMS::FriendlySAX->new
            ( -file    => $file,
              -tag     => ['popDat', 'rs'],
              -method  => \&load_popdat,  );
    };
    
}

sub load_popdat {
    my $node = shift;
    if ($node->{NAME} eq 'rs') {
        # We have finished processing the popDat records
        if ($mlPops) {
            # Write all the populations after we move out of the block
            my @pops = values %{$mlPops};
            foreach my $pop (@pops) {
                $pop->update();
            }
            $mlPops = undef;
            $args->msg("Updated ".scalar(@pops)." population entries", `date`);
        } elsif ($newPops) {
            $args->msg("$newPops new populations noted");
        }
        # Exit from ForkCritter:
        die 0;
    }
    my ($pid, $textName, $handle, $pars)
        = map { $node->{ATTR}{$_} } qw(id name handle class);
    return if ($simpPops{$pid});
    $newPops++;
    # We are just going to make the basic population here:
    my $name  = sprintf($popFmt, $pid);
    my $pop   = $mapLoc->get_population($name, $stAuth);
    my $pkey  = $pop->pkey;
    my $attr  = $simpPops{$pid} = {
        id     => $pid,
        name   => $textName,
        handle => $handle,
        class  => $pars,
        pkey   => $pkey,
    };
    # $mlPops ||= {};
    # $mlPops->{$pkey} ||= $pop;
}

sub purge_alleles {
    unless ($args->val(qw(confirm))) {
        $args->msg("[!!]", "You are about to delete all allele information associated with dbSNP populations",
                   "If you really wish to do this, please include the parameter -confirm");
        exit;
    }
    &init_maploc();
    $args->msg("[<]","Recovering all dbSNP populations");
    $rootPop = $mapLoc->get_population("dbSNP Populations", $stAuth);
    my @allPops = $rootPop->deep_children();
    $args->msg("[+]", "Found ".($#allPops+1)." populations, deleting alleles from each...");
    my $sth = $mapLoc->dbh->prepare
        ( -name => "Delete all alleles for a population",
          -sql  => "DELETE FROM allele WHERE pop_id = ?", );
    foreach my $pop (@allPops) {
        $args->msg("[-]", $pop->to_one_line());
        my $pid = $pop->pkey();
        $sth->execute($pid);
    }
    $args->msg("Finished");
    exit;
}

sub load_populations {
    &init_maploc();
    &set_meta_stuff('setXtra');
    my $dDir = "$dir";
    my $oFiles = &find_files($dDir, '^Pop.+.bcp.gz', "");
    my $kindCols = {
        PopMandLine => [qw(pop_id line_num mandline create_time last_updated_time)],
        PopLine     => [qw(pop_id line_num line create_time last_updated_time)],
        Population  => [qw(pop_id handle loc_pop_id loc_pop_id_upp create_time last_updated_time src_id)],
        PopClass => [qw(pop_id pop_class_id snp_count)],
        PopClassCode => [qw(pop_class_id pop_class pop_class_text)],
        Submitter => [qw(handle name fax phone email lab institution address create_time last_updated_time)],
    };


    my %taxGroup;
    my (%popClass, %pcCode);
    foreach my $path (@{$oFiles}) {
        my @pbits = split(/\//, $path);
        my $kind;
        if ($pbits[-1] =~ /^(.+)\.bcp\.gz$/) {
            $kind = $1;
        } else {
            $args->msg("[!]", "Unknown population file", $path);
            next;
        }
        my $cols = $kindCols->{$kind};
        unless ($cols) {
            $args->msg("[!]", "Failed to find table columns for $kind", $path);
            next;
        }
        unless ($path =~ /shared_data/) {
            my $tax = $pbits[-4];
            if (!$tax || $tax !~ /^[a-z]+_\d+$/) {
                $args->msg("[!]", "Failed to find taxa in file path", $path);
            } else {
                push @{$taxGroup{$tax}}, [ $path, $kind, $cols ];
            }
            next;
        }
        $args->msg("Loading global $kind data");
        my $tr     = BMS::TableReader->new
            ( -format    => 'tsv' );
        $tr->input($path);
        $tr->set_header( $cols );
        if ($kind eq 'PopClass') {
            while (my $row = $tr->next_clean_hash()) {
                my ($pid, $pcid, $snps) = map { $row->{$_} } @{$cols};
                $popClass{$pid}{$pcid} = $snps || 1 if ($pid && $pcid);
            }
        } elsif ($kind eq 'PopClassCode') {
            my $tr     = BMS::TableReader->new
                ( -format    => 'tsv' );
            $tr->input($path);
            $tr->set_header( $cols );
            while (my $row = $tr->next_clean_hash()) {
                my ($pcid, $pc, $desc) = map { $row->{$_} } @{$cols};
                next unless ($pcid && $pc);
                my $name = $pc;
                my $pop  = $mapLoc->get_population($name, $stAuth);
                $pop->parent($rootPop);
                $pcCode{$pcid} = $pop;
                $pop->tag_values("pop_class_id", $pcid);
                $desc = &clean_string($desc);
                $pop->tag_values("Description", $desc) if ($desc);
                $pop->update();
                # warn $pop->to_text();
            }
            my $unkPop = $pcCode{-999}[0] = $mapLoc->get_population
                ("Unclassified dbSNP", $stAuth);
            next;
        } elsif ($kind eq 'Submitter') {
            next;
        } else {
            $args->msg("[!]","Not sure what to do with $kind information");
        }
    }
    
    foreach my $tax (sort keys %taxGroup) {
        $args->msg("Parsing population data for $tax");
        my (%data, %ignored);
        foreach my $tg (@{$taxGroup{$tax}}) {
            my ($path, $kind, $cols) = @{$tg};
            $args->msg("[POP]", $kind);
            my $tr     = BMS::TableReader->new
                ( -format    => 'tsv' );
            $tr->input($path);
            $tr->set_header( $cols );
            while (my $row = $tr->next_clean_hash()) {
                my $pid = $row->{pop_id};
                unless ($pid) {
                    $args->msg("[!]","No pop_id", $args->branch($row));
                    next;
                }
                my $targ = $data{$pid} ||= {
                    pop_id => $pid,
                };
                while (my ($key, $val) = each %{$row}) {
                    if ($key =~ /_time$/) {
                        # Just keep the date, strip out the timestamp:
                        $val =~ s/\s+.+$//;
                        push @{$targ->{$key}}, $val;
                    } elsif ($key =~ /line$/) {
                        my $ind = $row->{line_num};
                        if (defined $ind) {
                            $targ->{$key}[$ind] = $val
                                if (defined $val && $val ne "");
                        } else {
                            $args->msg("[!]","No line number",
                                       $args->branch($row));
                        }
                    } elsif ($key =~ /^(handle|loc_pop_id)$/) {
                        $targ->{$key}{$val} = 1;
                    } else {
                        $ignored{$key}++;
                    }
                }
            }
        }
        map { delete $ignored{$_} } qw(line_num loc_pop_id_upp pop_id);
        my @igk = sort keys %ignored;
        $args->msg("Some population columns were ignored", 
                   map { "$_ $ignored{$_}" } @igk) unless ($#igk == -1);
        my @pops = sort { $a->{pop_id} <=> $b->{pop_id} } values %data;
        my $pFile = "$tax-Populations.txt";
        open(PF, ">$pFile") || $args->death
            ("Failed to write population summary", $pFile, $!);
        my %lids;
        # my $dbi = $mapLoc->dbi();
        foreach my $targ (@pops) {
            my $pid  = $targ->{pop_id};
            my @subm = keys %{$targ->{handle} || {}};
            my $name = sprintf($popFmt, $pid);
            my $pop   = $mapLoc->get_population($name, $stAuth);
            $pop->tag_values("pop_id", $pid);
            $pop->tag_values("Category", $locCatName);
            if ($#subm == 0) {
                $pop->tag_values("Submitter", $subm[0]);
            } else {
                $args->msg("[PopErr]", "pop_id = $pid : No Submitter");
            }
            my @lpids = keys %{$targ->{loc_pop_id} || {}};
            my $handle = $#lpids == 0 ? $lpids[0] : "";
            if ($handle) {
                $pop->tag_values("Handle", $handle);
                push @{$lids{$handle}}, $pid;
            } else {
                $args->msg("[PopErr]", "pop_id = $pid : Non-unique local IDs");
            }
            

            foreach my $key ('line','mandline') {
                if (my $arr = $targ->{$key}) {
                    my $desc = &clean_string
                        (join(' ', map { defined $_ ? $_ : "" } @{$arr}));
                    my $tag = $key eq 'line' ? 
                        "Description" : "Mandatory Description";
                    $pop->tag_values($tag, $desc) if ($desc);
                }
            }
            if (my $pcs = $popClass{$pid}) {
                my @pcids = sort { $pcs->{$b} <=> $pcs->{$a} } keys %{$pcs};
                my @pars  = map  { $pcCode{$_} } @pcids;
                $pop->parent( $pars[0] );
                $pop->tag_values("Population Class", map { $_->name } @pars);
            } else {
                $pop->parent( $pcCode{-999}[0] );
            }
            foreach my $key ('create_time','last_updated_time') {
            }
            $pop->update();
            print PF $pop->to_text();
        }
        close PF;
        while (my ($lid, $pidA) = each %lids) {
            $args->msg("[!]", "Duplicated local ID $lid: ".join
                       (",", @{$pidA})) if ($#{$pidA} != 0);
        }
        $args->msg(scalar(@pops)." populations recovered", $pFile);
    }
}

sub load_simple_asn {
    my $dat = shift;
    my ($chr, $l, $r, $bld) = map { $dat->{$_} } qw(Chr Lft Rgt Build);
    my $loc = $mapLoc->flanks_to_location($chr, $l, $r, $bld);
        
    unless ($loc) {
        &simp_err("Failed to recover location", "$chr.$bld:$l^$r");
        return;
    }
    foreach my $acc (split(/\,/, $dat->{Accs})) {
        $loc->add_accession_to_db( $acc, $stAuth );
    }
    if (my $ttxt = $dat->{Tags}) {
        foreach my $tv (split(/\|/, $ttxt)) {
            if ($tv =~ /^([^=]+)=(.+)$/) {
                $loc->tag_values($1, $2);
            }
        }
        $loc->write_tags();
    }
    if (my $all = $dat->{Alleles}) {
        # This indicates that the SNP is unlikely to have defined populations
        foreach my $base (split(/\//, $all)) {
            $loc->add_allele_to_db($base, $anonPopId);
        }
    }
    $loc->add_category( $locCategory );
    $loc->write_categories();
}

sub load_simple_gt {
    my $node = shift;
    my $snp = $nowSnp = $node->{ATTR}{id};
    my @locs;
    foreach my $lx (@{$node->{BYTAG}{loc} || []}) {
        next if ($lx->{ATTR}{ignore});
        my ($chr, $bld) = ( $lx->{ATTR}{chr},  $lx->{ATTR}{bld} );
        ($chr, $bld) = $mapLoc->standardize_chromosome
            ( $lx->{ATTR}{acc} ) unless ($chr && $bld);
        unless ($chr && $bld) {
            &simp_err("Failed to identify Chr and Build",
                      $lx->{ATTR}{acc} || '-UNDEF-');
            next;
        }
        my $loc = $mapLoc->flanks_to_location
            ($chr, $lx->{ATTR}{l}, $lx->{ATTR}{r}, $bld);
        $loc->add_accession_to_db( $snp, $stAuth );
        # Segregate locations with orientation set to -1 :
        my $ori = $lx->{ATTR}{ori} || 0;
        push @{$locs[ $ori < 0 ? 1 : 0]}, $loc;
        $loc->add_category( $locCategory );
        $loc->write_categories();
    }
    foreach my $px (@{$node->{BYTAG}{pop} || []}) {
        my $pid = $simpPops{ $px->{ATTR}{id} }{pkey};
        unless ($pid) {
            &simp_err("Failed to convert dbSNP population ID to MapLoc",
                      $px->{ATTR}{id} || '-UNDEF-');
            next;
        }
        my $sz = $px->{ATTR}{sz};
        my @bases = map { [ $_->{ATTR}{b}, 
                            $_->{ATTR}{f} ] } @{$px->{BYTAG}{f} || []};
        for my $i (0..$#locs) {
            if ($i) {
                # This location is -1 strand relative to the population report
                # We need to revcom the bases
                foreach my $bd (@bases) {
                    $bd->[0] = &clean_allele( $bd->[0], 1 );
                }
            }
            foreach my $loc (@{$locs[$i]}) {
                foreach my $bd (@bases) {
                    my ($b, $f) = @{$bd};
                    $loc->add_allele_to_db($b, $pid, $f, $sz);
                }
            }
        }
    }
    if ($doDump && $#locs != -1) {
        print "\n------------------\n";
        map { map {  print $_->to_text(); } @{$_} } @locs;
    }
    #  die BMS::FriendlySAX::node_to_text( $node );
}

sub report_files {
    my @files = @_;
    push @files, map { $statFiles{$_} } sort keys %statFiles;
    my @out;
    foreach my $file (@files) {
        if (my $sz = -s $file) {
            my $cmd = "less -S %s";
            $cmd = "gunzip -c %s | less -S" if ($file =~ /\.gz/);
            push @out, sprintf("[%7.3fMb] $cmd", $sz/1000000, $file);
        }
    }
    $args->msg("[FILES]", @out) unless ($#out == -1);
}

sub make_file_unique {
    my $file = shift;
    return unless ($file && -s $file);
    my %u;
    if (open(UNIQ, "<$file")) {
        while (<UNIQ>) {
            $u{$_}++;
        }
        close UNIQ;
        if (open(UNIQ, ">$file")) {
            foreach my $line (sort keys %u) {
                print UNIQ $line;
            }
            close UNIQ;
        } else {
            $args->err("[ERROR]", "Failed to write file for unique filtering",
                       $file, $!);
        }
    } else {
        $args->err("[ERROR]", "Failed to read file for unique filtering",
                   $file, $!);
    }
    
}

sub find_files {
    my $dDir = shift || "$dir/organisms";
    my $match = shift;
    my $mainF = shift;
    $mainF    = '(.+)\.xml\.gz$' unless (defined $mainF);
    my $cReq = uc($args->val(qw(chr chromosome)) || "");
    my $sReq = uc($args->val(qw(skipchr nochr
                                skipchromosome nochromosome)) || "");
    $args->msg("[LIMIT]", "Only analyzing Chromosome $cReq") if ($cReq);
    $args->msg("[LIMIT]", "Skipping Chromosome $sReq") if ($sReq);
    $args->msg("[LIMIT]", "Matching files like $match") if ($match);
    my $tid;
    if (my $sreq = $args->val(qw(taxa species))) {
        my $tti = $textMaps->{taxToId};
        my $sname;
        if ($sreq =~ /^\d+$/) {
            $tid = $sreq;
            $sname = $textMaps->{taxFromId}{$tid} || "/unknown name/";
        } elsif ($tid = $tti->{lc($sreq)}) {
            $sname = $textMaps->{taxFromId}{$tid} || $sreq;
        } else {
            my @known; 
            map { push @known, $_ if ($tti->{$_}) } sort keys %{$tti};
            $args->death("Failed to resolve '$sreq' to an NCBI TaxID",
                         "You can update the \$textMaps->{taxToId} lookup, or pick from the known list below:", map { "$_ [$tti->{$_}]" } @known);
        }
        $args->msg("[LIMIT]", "Request to analyze species $sname [$tid]");
    }
    my @toSort;
    my @toScan = ("");
    while ($#toScan != -1) {
        my $subDir = shift @toScan;
        my $sd     = "$dDir/$subDir";
        $sd =~ s/\/{2,}/\//g;
        # $args->msg("Scanning $sd ...");
        opendir(SCANDIR, $sd) || $args->death
            ("Failed to read data directory", $sd, $!);
        foreach my $file (readdir SCANDIR) {
            next if ($file =~ /^\.+$/);
            my $path = "$sd/$file";
            $path =~ s/\/{2,}/\//g;
            if (-d $path) {
                if ($tid && $file =~ /^([a-z]+)_(\d+)$/) {
                    my ($sn, $tax) = ($1,$2);
                    if ($tax != $tid) {
                        $args->msg_once("[LIMIT]","Skipping $sn [$tax]");
                        next;
                    }
                }
                my $newSD = "$subDir/$file";
                push @toScan, $newSD unless ($newSD eq $simpDir);
            } elsif (!$mainF || $file =~ /$mainF/) {
                my $base = $1 || $file;
                # Do not get testing limit files
                next if ($base =~ /limit/);
                my $sorter = $base;
                if ($base =~ /chr?([a-z0-9]+?)(-simple)?$/i) {
                    $sorter = uc($1);
                }
                # If requested, only get particular chromosomes
                next if ($cReq && $sorter !~ /^($cReq)$/);
                # If requested, ignore particular chromosomes
                next if ($sReq && $sorter =~ /^($sReq)$/i);
                next if ($match && $file !~ /$match/);
                if ($sorter =~ /^\d+$/) { $sorter = sprintf("%03d", $sorter) }
                push @toSort, [ $sorter, $path ];
            }
        }
        closedir(SCANDIR);
    }
    my @files = map { $_->[1] } sort { $a->[0] cmp $b->[0] } @toSort;
    $args->msg("Recovered ".scalar(@files). " files from $dDir");
    return \@files;
    
}

sub SIMPLIFY_ASN {
    my $dDir = "$dir/organisms";
    my $oFiles = &find_files($dDir, '^ds_flat', '(.+)\.flat\.gz$');
    my @files;
    foreach my $path (@{$oFiles}) {
        if ($path =~ /(.*?)([^\/]+)\.flat\.gz$/) {
            my ($prfx, $base) = ($1, $2);
            my $tsd  = $prfx;
            $tsd     =~ s/\/ASN1_flat/\//;
            $tsd     =~ s/.+organisms\///;
            $base .= '-limit' if ($limit);
            my $targ = sprintf("%s/%s/%s-simple.tsv.gz",
                               $simpDir, $tsd, $base);
            $targ =~ s/\/{2,}/\//g;
            push @files, [ $path, $targ ] unless (-s $targ && !$limit);
        } else {
            $args->msg("[ERR]", "File name not recognized", $path);
        }
    }
    $args->msg("Simplifying ".($#files+1)." ASN1 files", 
               "SRC: $dDir", "DST: $simpDir");
    $fc = BMS::ForkCritter->new
        ( -inputtype   => 'array',
          -input       => [ sort { $a->[0] cmp $b->[0] } @files ],
          -method      => \&simplify_asn,
          -progress    => $progress,
          -verbose     => $vb );

    map { &alter_statfile( $_, 'SimplifyASN') } ('ERROR','NOTE');
    &run_fork();
    &make_file_unique( $statFiles{NOTE} );
    &report_files();
}

sub SIMPLIFY_GT {
    my $dDir = "$dir/organisms";
    my $oFiles = &find_files($dDir, '^gt_');
    my @files;
    foreach my $path (@{$oFiles}) {
        if ($path =~ /(.*?)([^\/]+)\.xml\.gz$/) {
            my ($prfx, $base) = ($1, $2);
            my $tsd  = $prfx;
            $tsd     =~ s/\/XML/\//;
            $tsd     =~ s/.+organisms\///;
            $base .= '-limit' if ($limit);
            my $targ = sprintf("%s/%s/%s-simple.xml.gz",
                               $simpDir, $tsd, $base);
            $targ =~ s/\/{2,}/\//g;
            push @files, [ $path, $targ ] unless (-s $targ && !$limit);
        } else {
            $args->msg("[ERR]", "File name not recognized", $path);
        }
    }
    $args->msg("Simplifying ".($#files+1)." Genotype XML files", 
               "SRC: $dDir", "DST: $simpDir");

    $fc = BMS::ForkCritter->new
        ( -inputtype   => 'array',
          -input       => [ sort { $a->[0] cmp $b->[0] } @files ],
          -method      => \&simplify_genotype,
          -progress    => $progress,
          -verbose     => $vb );

    map { &alter_statfile( $_, 'SimplifyGeno') } ('ERROR','NOTE');
    &run_fork();
    &make_file_unique( $statFiles{NOTE} );
    &report_files();
}

sub alter_statfile {
    my ($key, $tag) = @_;
    if (my $f = $statFiles{$key}) {
        my ($base, $sfx) = ($f, "");
        if ($base =~ /(.+)\.([^\.]{2,6})$/) {
            ($base, $sfx) = ($1, $2);
        }
        $f = $base."-$tag";
        $f .= ".$sfx" if ($sfx);
        $statFiles{$key} = $f;
    }
}

sub run_fork {
    $args->death("Can not fork ForkCritter object!")
        unless ($fc);
    while (my ($key, $file) = each %statFiles) {
        $fc->output_file( $key, ">>$file" );
    }
    if (my $failed = $fc->execute( $forkNum || 1)) {
        $args->death("$failed jobs failed to properly convert");
    }
}

sub SIMPLIFY {
    my $dDir = "$dir/organisms";
    $args->msg("Simplifying Original XML", $dDir);
    my $oFiles = &find_files($dDir, '^gs_');
    my @files;
    foreach my $path (@{$oFiles}) {
        if ($path =~ /(.*?)([^\/]+)\.xml\.gz$/) {
            my ($prfx, $base) = ($1, $2);
            my $tsd  = $prfx;
            $tsd     =~ s/\/XML/\//;
            $tsd     =~ s/.+organisms\///;
            $base .= '-limit' if ($limit);
            my $targ = sprintf("%s/%s/%s-simple.xml.gz",
                               $simpDir, $tsd, $base);
            $targ =~ s/\/{2,}/\//g;
            push @files, [ $path, $targ ] unless (-s $targ && !$limit);
        } else {
            $args->msg("[ERR]", "File name not recognized", $path);
        }
    }

    $fc = BMS::ForkCritter->new
        ( -inputtype   => 'array',
          -input       => \@files,
          -method      => \&simplify_xml,
          # -limit       => $limit,
          -progress    => $progress,
          -verbose     => $vb );
    map { &alter_statfile( $_, 'Simplify') } ('ERROR');
    &run_fork();
    &make_file_unique( $statFiles{NOTE} );
    &make_file_unique( $statFiles{ERROR} );
    &report_files();
}

sub simplify_asn {
    my ($fd) = @_;
    my ($src, $trg) = @{$fd};
    my $unz  = $trg;
    $unz     =~ s/\.gz$//;
    $args->msg("[SIMPLIFY]", $src, (map { "$_.gz" } ($unz)));
    my $tr = BMS::TableReader->new();
    my $fh = $tr->_filehandle_for_file( $src );
    my $dat;
    my $done = 0;
    my $txm  = &text_maps();
    my $asmOk = '^('.join('|', values %{$txm->{buildmap}}).')';
    &init_maploc();
    open(SIMPOUT, ">$unz") || $args->death
        ("Failed to write to simplified file", $unz, $!);
    print SIMPOUT join("\t", qw(Chr Build Lft Rgt Pos Accs Type Tags Alleles))."\n";
    while (<$fh>) {
        s/[\n\r]+$//;
        if (/^\s*$/) {
            $done += &_write_simple_loc( $dat );
            last if ($limit && $done >= $limit);
            $dat = {};
        }
        if (/^rs(\d+) \|/) {
            push @{$dat->{acc}}, $1;
        } elsif (/^CTG \| (.+)/) {
            my $ctg = $1;
            my %loc;
            foreach my $bit (split(/\s+\|\s+/, $ctg)) {
                if ($bit =~ /([^=]+)=(.+)/) {
                    my ($k, $v) = ($1, $2);
                    $v =~ s/\.\d+$// if ($k eq 'assembly');
                    $loc{$k} = $v;
                }
            }
            if (my $asm = $loc{assembly}) {
                if ($asm =~ /$asmOk/) {
                    $loc{assembly} = $1;
                    push @{$dat->{ctg}}, \%loc;
                } else {
                    $args->msg_once("Ignoring assembly '$asm'")
                        if ($fc->child() == 1);
                }
            }
        } elsif (/^CLINSIG \| ([^\|]+)/) {
            # Hey, might as well grab it...
            my $cs = $1;
            if ($cs =~ /.+=(.+)/) { $cs = $1; }
            $cs =~ s/\s+$//;
            $dat->{tags}{'Clinical Significance'}{$cs} = 1;
        } elsif (/^VAL .+byFrequency/) {
            $dat->{freqVal} = 1;
        } elsif (/^SNP.+alleles=\'([^\']+)\'/) {
            $dat->{alleles}{$1} = 1;
        }
    }
    close SIMPOUT;
    foreach my $f ($unz) {
        system("gzip -f $f");
    }
}

sub _write_simple_loc {
    my $dat = shift;
    return 0 unless ($dat);
    my @ids = map { "rs$_" } @{$dat->{acc} || []};
    if ($#ids == -1) {
        &simp_err("No accession found for entry");
        return 0;
    }
    my $idtxt = join(',', @ids);
    my @ctgs = @{$dat->{ctg} || []};
    if ($#ctgs == -1) {
        &simp_err("No locations found", $idtxt);
        return 0;
    }
    my $done = 0;
    for my $c (0..$#ctgs) {
        my $ctg = $ctgs[$c];
        my ($chr, $build, $pos, $lt) = map { $ctg->{$_} }
        ('chr', 'assembly', 'chr-pos', 'loctype');
        my ($l, $r);
        if (!defined $pos) {
            &simp_err("No Position defined for CTG $c, loctype $lt", $idtxt);
            next;
        } elsif ($pos !~ /^\d+$/) {
            &simp_err("Ignoring pos '$pos' for CTG $c, loctype $lt", $idtxt)
                unless ($pos eq '?');
            next;
        } elsif (!defined $chr) {
            &simp_err("No CHR for CTG $c, loctype $lt", $idtxt);
            next;
        } elsif ($lt == 2) {
            # This is a classic SNP
            ($l, $r) = ($pos - 1, $pos + 1);
        } elsif ($lt == 1) {
            # This is an insertion == MNP
            my ($s, $e) = ($ctg->{'ctg-start'}, $ctg->{'ctg-end'});
            if ($s && $e && $s <= $e) {
                ($l, $r) = ($pos - 1, $pos + $e - $s + 1);
            } else {
                &simp_err("Failed to resolve loctype $lt", $idtxt);
                next;
            }
        } elsif ($lt == 3) {
            # This is a deletion
            ($l, $r) = ($pos, $pos + 1);
        } else {
            &simp_err("Ignoring loctype $lt", $idtxt);
            next;
        }
        $done++;
        my $ni = $mapLoc->pretty_flanks($l, $r);
        my @tt;
        foreach my $tag (sort keys %{$dat->{tags}}) {
            foreach my $val (sort keys %{$dat->{tags}{$tag}}) {
                push @tt, "$tag=$val";
            }
        }
        my $tags = join('|', @tt) || "";
        my $all = "";
        if (!$dat->{freqVal} && $dat->{alleles}) {
            # This location does not have frequency validation
            # That makes it unlikely that any alleles have been assigned
            my @atxt = keys %{$dat->{alleles}};
            if ($#atxt == 0) {
                $all = $atxt[0];
                if ($ctg->{orient} eq '-') {
                    my @bases = map { &clean_allele($_, 1) } 
                    split(/\//, $atxt[0]);
                    $all = join('/', @bases);
                }
            } else {
                &simp_err("Multiple alleles defined", $idtxt, @atxt);
            }
        }
        print SIMPOUT join
            ("\t", $chr, $build, $l, $r, $ni, $idtxt, $lt, $tags, $all)."\n";
        # warn sprintf("%10s %s.%s:%s\n", $idtxt, $chr, $build, $ni) unless ($ni eq $pos);
    }
    #unless ($done) {
    #    $args->msg("[!]", "Location not inferred for ".$idtxt);
    #
    #}
    return $done;
}

sub simplify_genotype {
    my ($fd) = @_;
    my ($src, $trg) = @{$fd};
    my $unz  = $trg;
    $unz     =~ s/\.gz$//;
    my $popF = $unz;
    $popF    =~ s/\.xml$//;
    $popF   .= "-Populations.tsv";

    &init_maploc();
    &set_taxa_from_file_name( $src );

    my $bld = &extract_build( $src );
    $args->msg("[SIMPLIFY]", $src, (map { "$_.gz" } ($unz)), $popF);
    $args->assure_directory($unz, 'isFile');
    open(SIMPOUT, ">$unz") || $args->death
        ("Failed to write to simplified file", $unz, $!);
    print SIMPOUT '<?xml version="1.0" encoding="UTF-8"?>'."\n";
    print SIMPOUT "<simpleRefSnp dbSnpBuild='$bld'>\n";
    eval {
        my $fs = BMS::FriendlySAX->new
            ( -file    => $src,
              -tag     => ['SnpInfo', 'Population'],
              -limit   => $limit,
              -nolimit => ['Population'],
              -skip    => ['GTypeByInd'],
              -method  => \&simplify_gt_rec,  );
    };
    if ($@ && $@ !~ /^\s+$/) {
        my $expected =  'user limit';
        unless ($@ =~ /\Q$expected\E|\[IGNORE\]/i) {
            $args->err("FriendlySAX error", $@);
        }
    }
    print SIMPOUT "</simpleRefSnp>\n";
    close SIMPOUT;
    if (open(SIMPOP, ">$popF")) {
        my @pops = sort { $a->{id} <=> $b->{id} } values %simpPops;
        my @head = qw(id name handle class size);
        foreach my $sd (@pops) {
            my %locsd = %{$sd || {}};
            # Sizes are stored in a hash and collected during processing
            $locsd{size} = join($joiner, map { "$_:$sd->{size}{$_}" }
                               sort { $a <=> $b } keys %{$sd->{size}});
            my @row = map { $locsd{$_} } @head;
            print SIMPOP join("\t", @row)."\n";
        }
        close SIMPOP;
    } else {
        $args->err
        ("Failed to write to population file", $popF, $!);
    }
    foreach my $f ($unz) {
        system("gzip -f $f");
    }
}

sub simplify_xml {
    my ($fd) = @_;
    my ($src, $trg) = @{$fd};
    my $unz  = $trg;
    $unz     =~ s/\.gz$//;
    my $ssfa = $unz;
    $ssfa    =~ s/\.xml$//;
    my $ssX  = $ssfa . "-SSids.xml";
    $ssfa   .= "-SS.fa";

    eval {
        # Quick pass of file to get species
        my $fs = BMS::FriendlySAX->new
            ( -file       => $src,
              -tag        => 'SourceDatabase',
              -limit      => 1,
              -quietlimit => 1,
              -verbose    => 0,
              -method     => \&set_taxa,  );
    };
    my $bld = &extract_build( $src );



    $args->msg("[SIMPLIFY]", $src, map { "$_.gz" } ($unz, $ssX, $ssfa));
    $args->assure_directory($unz, 'isFile');
    open(SIMPOUT, ">$unz") || $args->death
        ("Failed to write to simplified file", $unz, $!);
    open(SIMPSSX, ">$ssX") || $args->death
        ("Failed to write to SS ID file", $ssX, $!);
    open(SIMPSSFA, ">$ssfa") || $args->death
        ("Failed to write to simplified SS Fasta File", $ssfa, $!);
    print SIMPOUT '<?xml version="1.0" encoding="UTF-8"?>'."\n";
    print SIMPSSX '<?xml version="1.0" encoding="UTF-8"?>'."\n";
    print SIMPOUT "<simpleRefSnp dbSnpBuild='$bld'>\n";
    print SIMPSSX "<simpleSSIds dbSnpBuild='$bld'>\n";
    eval {
        my $fs = BMS::FriendlySAX->new
            ( -file    => $src,
              -tag     => ['Rs'],
              -limit   => $limit,
              # -skip    => $skip,
              -method  => \&simplify_record,  );
    };
    if ($@ && $@ !~ /^\s+$/) {
        my $expected =  'user limit';
        unless ($@ =~ /\Q$expected\E|\[IGNORE\]/i) {
            $args->err("FriendlySAX error", $@);
        }
    }
    print SIMPOUT "</simpleRefSnp>\n";
    print SIMPSSX "</simpleSSIds>\n";
    close SIMPOUT;
    close SIMPSSX;
    close SIMPSSFA;
    foreach my $f ($unz, $ssfa, $ssX) {
        system("gzip -f $f");
    }
}

sub simplify_gt_rec {
    my $node = shift;
    my $fh    = *SIMPOUT;
    my $pad   = " ";
    if ($node->{NAME} eq 'Population') {
        my ($pid, $name, $hand)
            = map { $node->{ATTR}{$_} } qw(popId locPopId handle);
        my %pcs;
        foreach my $pcn (@{$node->{BYTAG}{popClass} || []}) {
            if (my $pc = $pcn->{ATTR}{self}) { $pcs{$pc} = 1}
        }
        my $attr = $simpPops{$pid} = {
            id => $pid,
            name => $name,
            handle => $hand,
            class  => join($joiner, sort keys %pcs),
        };
        print $fh &singleChild('popDat', $attr, $pad);
        $attr->{size} = {};
        return;
    }
    my $snp_acc = $nowSnp = $node->{ATTR}{rsId};
    unless ($snp_acc && $snp_acc =~ /^\d+$/) {
        $nowSnp ||= '-UNDEF-';
        &simp_err("Poorly formed SNP accession");
        return 0;
    }
    my $pad2  = " $pad";
    my $attr  = {id => 'rs'.$snp_acc};

    my %strands;
    my @locs;
    foreach my $sl (@{$node->{BYTAG}{SnpLoc} || []}) {
        my %param = %{$sl->{ATTR} || {}};
        my $assm = $param{genomicAssembly};
        my ($build);
        if ($assm =~ /^[^\:]+\:(.+?)(\.p\d+)?$/) {
            $build = $textMaps->{buildmap}{$1};
        }
        unless ($build) {
            &note_once("Ignoring genome build", $assm);
            next;
        }
        my $chr = $param{chrom};
        my $acc = &chr_to_acc( $chr, $build );
        my $lattr = { chr => $chr, bld => $build };
        if ($acc) {
            $lattr->{acc} = $acc;
        } else {
            &note_once("Ignoring unmappable chromosome", $chr) if ($chr);
            $lattr->{ignore}++;
        }
        # As far as I can tell, the reference allele is in scaffold +1 frame
        # That is, it is the allele as always read on the 'top strand' of
        # the reference contig:
        my $ref = &clean_allele($param{contigAllele});

        # The strand relates the /relative/ orientation of the rs entry
        # to the reference
        my $str = 
            $textMaps->{orient}{$sl->{ATTR}{rsOrientToChrom} || ""} || 0;
        my ($l, $r);


        # Start and End are ZERO INDEXED!
        my ($s, $e, $lt) = ($param{start}, $param{end}, $param{locType});

        # http://www.ncbi.nlm.nih.gov/SNP/specs/alignment_types.htm
        if (!defined $s) {
            # Not sure why entries exist without start coordinates. Example:
            # <SnpLoc genomicAssembly="37:GRCh37.p5" chrom="6" locType="2" contigAllele="G" contig="NT_167244:1"/>
            next;
            &attr_err($lattr, "No start coordinate provided");
            $lattr->{locType} = $lt || '-UNDEF-';
            $lattr->{ignore}++;
        } elsif ($lt == 2) {
            # This is a classic SNP
            ($l, $r) = ($s, $s + 2);
        } elsif ($lt == 1) {
            # This is an insertion == MNP
            if (defined $e) {
                ($l, $r) = ($s, $e + 2);
                $lattr->{typ} = 'MNP';
            } else {
                &attr_err($lattr, "MNP position does not define end");
                $lattr->{ignore}++;
            }
        } elsif ($lt == 3) {
            # This is a deletion
            ($l, $r) = ($s + 1, $e + 1);
            $lattr->{typ} = 'DEL';
        } elsif ($lt == 4) {
            # This is a "Range insertion"
            # AFAICT, it is very hard to clearly define this position
            $lattr->{typ} = 'RNG';
            $lattr->{ignore}++;
        } elsif ($lt == 5) {
            # This is a "Range substitution"
            # AFAICT, it is very hard to clearly define this position
            $lattr->{typ} = 'RNGSUB';
            $lattr->{ignore}++;
        } elsif ($lt == 6) {
            # This is a "Range deletion"
            # AFAICT, it is very hard to clearly define this position
            $lattr->{typ} = 'RNGDEL';
            $lattr->{ignore}++;
        } else {
            my $msg = "Unrecognized location type '$lt'";
            map { $_ ||= '-UNDEF-' } ($s, $e);
            $args->msg_once($msg, "$s..$e");
            $lattr->{typ} = 'UNK';
            $lattr->{ignore}++;
            &attr_err($lattr, $msg); 
        }
        $lattr->{l} = $l;
        $lattr->{r} = $r;
        $lattr->{ref}  = $ref;
        $lattr->{ori}  = $str if ($acc);
        push @locs, $lattr;
        # $slTxt .= &singleChild('loc', $lattr, $pad2);
        # warn BMS::FriendlySAX::node_to_text( $sl );
    }
    # If there is a common strand in the genomic locations, use it and
    # correct the alleles to be normalized on the +1 genome strand
    my %locStrands; map { $locStrands{$_->{ori}}++ } @locs;
    my @allLocStrs = keys %locStrands;
    my $rsMapStr = 0;
    if ($#allLocStrs == 0) {
        # There is a consistent strand for all mapped locations
        # We can use this to normalize information to a common strand
        # and save us time later
        $rsMapStr = $attr->{rsStr} = $allLocStrs[0];
        # We can also delete the orientation information from the data
        map { delete $_->{ori} } @locs;
    }
    my %slTxts;
    map { $slTxts{ &singleChild('loc', $_, $pad2) } = 1 } @locs;
    my $slTxt = join('', sort keys %slTxts) || "";
    $attr->{locN} = $#locs + 1;
    

    my %freqs;
    foreach my $ss (@{$node->{BYTAG}{SsInfo} || []}) {
        my $ssStr  = $textMaps->{orient}{$ss->{ATTR}{ssOrientToRs} || ""} || 0;
        my $mapStr = $ssStr * $rsMapStr;
        my $do_rev = ($mapStr < 0) ? 1 : 0;
        foreach my $pop (@{$ss->{BYTAG}{ByPop} || []}) {
            my $pid = $pop->{ATTR}{popId};
            my $sz  = $pop->{ATTR}{sampleSize} || 0;
            $simpPops{$pid}{size}{$sz}++;
            next if ($sz < $minSampleSize);
            my %sf;
            my $do_rev_here = $do_rev;
            my $needsOri;
            if (!$ssStr) {
                # We do not know the strand of this position
                $sf{unknownStrand} = 1;
                $do_rev_here = 0;
            } elsif (!$mapStr && $ssStr < 0) {
                # If we were unable to unambiguously determine a consistent
                # genomic orientation, at least flip over alleles to make
                # all the SS entries consistent with the RS orientation
                $do_rev_here = 1;
                # $needsOri = 1;
            }
            foreach my $af (@{$pop->{BYTAG}{AlleleFreq} || []}) {
                my $base = &clean_allele($af->{ATTR}{allele}, $do_rev_here);
                my $freq = $af->{ATTR}{freq};
                if (defined $sf{$base}) {
                    my $prior = $sf{$base};
                    &attr_err($attr, "Multiple frequencies for allele", 
                              "$pid:$base - $prior vs $freq [$sz]")
                        unless ($prior eq $freq);
                } else {
                    $sf{$base} = $freq;
                }
            }
            push @{$freqs{$pid}}, [ $sz, \%sf, $needsOri ];
        }
    }
    print $fh &tagOpen('rs', $attr, $pad) . "\n";
    print $fh $slTxt;
    foreach my $pid (sort { $a <=> $b } keys %freqs) {
        my ($largest) = sort { $b->[0] <=> $a->[0] } @{$freqs{$pid}};
        my ($sz, $sf, $needsOri) = @{$largest};
        print $fh &tagOpen('pop', {
            id => $pid, sz => $sz, 
            # needsOri => $needsOri,
        }, $pad2);
        foreach my $base (sort {$sf->{$a} 
                                <=> $sf->{$b}} keys %{$sf}) {
            my $fattr = { b => $base, f => $sf->{$base}, NoCR => 1 };
            print $fh &singleChild('f', $fattr, "");
        }
        print $fh &tagClose('pop');
    }
    print $fh &tagClose('rs', $pad);
    # warn BMS::FriendlySAX::node_to_text( $node );
    # warn $args->branch($node);
}

sub chr_to_acc {
    my $chr = shift;
    return undef unless ($chr && $cfrm);
    if ($chr eq 'Un' || $chr eq 'ChrUn') {
        return undef;
    } elsif ($chr =~ /^(\d+|[WXYZ]|MT)$/) {
        return sprintf($cfrm, $chr, shift || $globalBuild);
    }
}

sub simplify_record {
    my $node = shift;
    my $fh = *SIMPOUT;
    my $snp_acc = $nowSnp = $node->{ATTR}{rsId};
    unless ($snp_acc && $snp_acc =~ /^\d+$/) {
        $nowSnp ||= '-UNDEF-';
        &simp_err("Poorly formed SNP accession");
        return 0;
    }
    my $class = $node->{ATTR}{snpClass} || 'unknown';
    my $pad   = " ";
    my $attr  = {
        id  => $node->{acc} = 'rs'.$snp_acc,
        tax => $node->{ATTR}{taxId},
    };
    if ($class eq 'named') {
        &simp_err("Not a useful variant", "Class $class");
        return;
    } elsif (my $mtclass = $textMaps->{class}{$class}) {
        $attr->{mtClass} = $mtclass;
    }
    &simple_dates( $node, $attr );
    my $ro = &tagOpen('rs', $attr, $pad) . "\n";
    print $fh $ro;
    print SIMPSSX $ro;
    my $type = $node->{ATTR}{snpType};
    if (my $com = $textMaps->{type}{$type}) {
        # The SNP is deemed unworthy
        print $fh "$pad <deprecated>$com</deprecated>\n";
    }
    print $fh &simp_validation( $node, $pad );
    print $fh &simp_sequence( $node, $pad );
    print $fh &simp_ss( $node, $pad );
    print $fh &simp_map( $node, $pad );

    # print BMS::FriendlySAX::node_to_text( $node );
    my $rc = &tagClose('rs', $pad);
    print $fh $rc;
    print SIMPSSX $rc;
}

sub simp_map {
    my ($root, $pad) = @_;
    my $list   = $root->{BYTAG}{Assembly} || [];
    return "" if ($#{$list} == -1);
    $pad .= " ";
    my $pad2 = "$pad ";
    # OCTOBER 2005 - ALL COORDINATES ARE ZERO-BASED
    my $xml = "";
    my (@fxndata, @locdata, %coms, @counts);
    foreach my $asm (@{$list}) {
        # Only really care about reference assemblies
        my $ref = $asm->{ATTR}{reference} || "";
        my $attr = {};
        next unless ($ref eq 'true');
        # <Assembly dbSnpBuild="132" genomeBuild="37_1" groupLabel="GRCh37" current="true" reference="true">
        my $comps = $asm->{BYTAG}{Component};
        unless ($comps) {
            &attr_err($attr, "Assembly with no components");
            next;
        }
        my ($build, $bnum);
        if (my $gl = $asm->{ATTR}{groupLabel}) {
            if ($gl =~ /^([A-Z]+(\d+))(\.p\d+)?$/i) {
                ($gl, $bnum)  = ($1, $2);
                unless ($build = $textMaps->{buildmap}{$gl}) {
                    &attr_err($attr, "Failed to recognize genome build", $gl);
                    $args->msg_once("Please update \$textMaps->{buildmap}");
                }
            } else {
                &attr_err($attr, "Failed to parse genome groupLabel", $gl);
                next;
            }
        } else {
            &attr_err($attr, "No groupLabel provided");
            next;
        }
        $attr->{build} = $build;
        foreach my $stat (@{$asm->{BYTAG}{SnpStat} || []}) {
            my $count = $stat->{ATTR}{seqlocCount} || 0;
            $attr->{copies} = $count if ($count > 1);
        }
        $xml .= &tagOpen( 'map', $attr, $pad)."\n";
        foreach my $comp (@{$comps}) {
            # Chromosome name, Contig accession
            # Contig start, Contig end, Contig/Chr relative orientation
            my ($chr, $acc, $ccs, $cce, $ori) = map { $comp->{ATTR}{$_} }
            qw(chromosome accession start end orientation);
            my $locs = $comp->{BYTAG}{MapLoc};
            unless ($locs) {
                &attr_err($attr, "Component with no MapLocs");
                next;
            }

            $ori = $textMaps->{orient}{$ori || ''} || 0;
            my $cname;
            if (!$chr) {
                &attr_err($attr, "No chromosome defined");
            } elsif ($cfrm) {
                if ($chr eq 'Un' || $chr eq 'ChrUn') {
                    # Do not map to unknown chromosomes
                } elsif ($chr =~ /^(\d+|[WXYZ]|MT)$/) {
                    $cname = sprintf($cfrm, $chr, $build);
                } else {
                    &attr_err($attr, "Unknown chromosome", $chr);
                }
            }
            foreach my $loc (@{$locs}) {
                my ($lf, $rf, $score, $str) = map { $loc->{ATTR}{$_} }
                qw(leftContigNeighborPos rightContigNeighborPos
                   alnQuality orient);
                $str = $textMaps->{orient}{$str} || 0;
                $score *= 100 if ($score);
                # $lf++; $rf--; # Move from flank to absolute positions
                my $cstr;
                if ($cname && $ccs) {
                    # Record the position on the chromosome
                    my ($cs, $ce);
                    if ($ori < 0) {
                        # The contig is oriented reverse to the chromosome
                        ($cs, $ce, $cstr) = ($cce-$rf+1, $cce-$lf+1, 0-$str);
                    } else {
                        ($cs, $ce, $cstr) = ($lf+$ccs+1, $rf+$ccs+1, $str);
                    }
                    $xml .= &singleChild('acc', {
                        id   => $cname,
                        l    => $cs,
                        r    => $ce,
                        str  => $cstr,
                        sc   => $score,
                    }, $pad2);
                    push @locdata, [ $cname, $cs, $ce, $cstr, $score, 'gdna']
                }
                # Record the position on the accession
                push @locdata, [ $acc, $lf, $rf, $str, $score, 'gdna'];
                $xml .= &singleChild('acc', {
                    id   => $acc,
                    l    => $lf,
                    r    => $rf,
                    str  => $str,
                    sc   => $score,
                }, $pad2);
                foreach my $fxnset (@{$loc->{BYTAG}{FxnSet} || []}) {
                    my %fsAttr = %{$fxnset->{ATTR}};
                    if ($cstr  && $fsAttr{allele}) {
                        my $gst = &mrna_orientation
                            ($cname,$fsAttr{mrnaAcc}) || 0;
                        $fsAttr{snpori} = $gst * $cstr;
                    }
                    push @fxndata, \%fsAttr;
                }
            }
        }
        $xml .= &tagClose( 'map', $pad);
    }

    my %seqinfo;
    foreach my $data (@fxndata) {
        my ($racc, $rv, $tc, $sym, $gid, $pacc, $pv) = map { $data->{$_} }
        qw(mrnaAcc mrnaVer fxnClass symbol geneId protAcc protVer);
        my $class = $textMaps->{impact}{$tc || ''};
        unless (defined $class) {
            $args->msg_once("Unknown fxnClass = $tc");
            next;
        }
        if ($racc && $rv) {
            my $raccV = $racc.'.'.$rv;
            $seqinfo{$raccV}{acc}   = $racc;
            $seqinfo{$raccV}{class} ||= $class;
            $seqinfo{$raccV}{type}  = 'mrna';
            $seqinfo{$raccV}{sym}   = $sym;
            if ($pacc && $pv) {
                my $paccV = $pacc.'.'.$pv;
                my ($allele, $aa, $pos) = map { $data->{$_} }
                qw(allele residue aaPosition);
                $seqinfo{$paccV}{acc}   = $pacc;
                $seqinfo{$paccV}{class} ||= $class;
                $seqinfo{$paccV}{type}  = 'prot';
                $seqinfo{$paccV}{sym}   = $sym;
                if ($pos) {
                    $seqinfo{$paccV}{pos}   = $pos;
                    $seqinfo{$paccV}{str} ||= $data->{snpori};
                }
                $seqinfo{$paccV}{rna}   = $raccV;
                if ($allele && $aa) {
                    $seqinfo{$paccV}{impact}{$allele} = $aa;
                    if ($aa eq '*') {
                        $seqinfo{$raccV}{class} = 'STP';
                        $seqinfo{$paccV}{class} = 'STP';
                    }
                }
            }
        }
    }

    my @accs = keys %seqinfo;
    foreach my $accV (@accs) {
        my $data = $seqinfo{$accV};
        my $type = $data->{type};
        my $str  = $data->{str} || undef;
        my $attr = { id => $accV, str => $str, type => $type,
                     sym => $data->{sym}, rna => $data->{rna} };
        if (my $class = $data->{class}) {
            if (my $mti = $textMaps->{mtImpact}{$class}) {
                $attr->{imp} = $mti;
            } else {
            }
        }
        my @aas;
        if ($str) {
            while ( my ($base, $aa) = each %{$data->{impact}}) {
                push @aas, "$base:$aa";
            }
        }
        $attr->{prot} = join($joiner, @aas);
        $xml .= &singleChild('impact', $attr, $pad);
    }
    return $xml;
}

sub simp_ss {
    my ($root, $pad) = @_;
    my $sslist = $root->{BYTAG}{Ss} || [];
    return "" if ($#{$sslist} == -1);
    
    $pad .= " ";
    my $xml = "";
    my $aliAttr = {};
    my %nums;
    my %methods;
    my @ssids;
    foreach my $ssnode (@{$sslist}) {
        my %attrs  = %{$ssnode->{ATTR}};
        my $ssid   = $attrs{ssId};
        my $ssattr = { NoCR => 1 };
        if (!$ssid || $ssid !~ /^\d+$/) {
            $ssid ||= '-UNDEF-';
            &attr_err($ssattr, "Malformed ss ID", $ssid);
            next;
        }
        push @ssids, $ssid;
        $ssid = $ssattr->{id} = 'ss'.$ssid;
        if (my $alias = $attrs{locSnpId}) {
            if ($alias =~ /^\d+$/) {
                $nums{$alias} = 1;
            } else {
                $aliAttr->{a}{$alias} = 1;
            }
        }

        my $orient = $ssattr->{ori} = $textMaps->{orient}{$attrs{orient}};
        unless ($orient) {
            &attr_err($ssattr, "Unknown relative ss orientation", $ssid);
        }
        my @bases;

        # Set up the SS Protocol
        my $metCl   = $attrs{methodClass} || '';
        my $mapmeth = $textMaps->{method}{ lc($metCl) };
        unless ($mapmeth) {
            &attr_err($ssattr, "Unknown ss method", $ssid, $metCl) if ($metCl);
            $mapmeth = 'Unknown';
        }
        if (my $using = $attrs{molType}) {
            $mapmeth .= " using $using";
        }
        my $auth    = $attrs{handle} || "Unknown";
        my $seqlist = $ssnode->{BYTAG}{Sequence} || [];
        push @{$methods{$mapmeth}{$auth}}, $ssattr;
        if ($#{$seqlist} < 0) {
            &attr_err($ssattr, "No sequence for SS", $ssid);
        } elsif ($#{$seqlist} > 0) {
            &attr_err($ssattr, "Multiple sequences for SS", $ssid);
        } else {
            my $seqnode = $seqlist->[0];
            @bases = &simple_observed_to_bases
                ($seqnode, 'Observed', $ssattr);
            my ($l, $r) = ("","");
            for my $side (5, 3) {
                my $nodes = $seqnode->{BYTAG}{"Seq$side"} || [];
                if ($#{$nodes} < 0) {
                    next;
                } elsif ($#{$nodes} > 0) {
                    &attr_err($ssattr, "Multiple ${side}' sequences", $ssid);
                    next;
                }
                my $seqdata = $nodes->[0]{TEXT};
                $seqdata =~ s/\s+//g;
                if ($side == 5 ) {
                    $l = $seqdata;
                } else {
                    $r = $seqdata;
                }
            }
            my $ap   = length($l) + 1;
            my $fa = sprintf(">%s allelePos=%d alleles=%s authority=%s",
                             $ssid, $ap, join('/', @bases), $auth );
            if (my $url = $attrs{linkoutUrl}) {
                $url =~ s/^\s+//; $url =~ s/\s+$//;
                $fa .= " url=$url";
            }
            my $seq = join('', $l, 'N', $r);
            my $sl  = length($seq);
            for (my $i = 0; $i < $sl; $i += $seqBlock) {
                $fa .= "\n".substr($seq, $i, $seqBlock);
            }
            print SIMPSSFA "$fa\n";
        }
    }
    $aliAttr->{numeric} = join($joiner, sort {$a <=> $b} keys %nums);
    $aliAttr->{ssids} = join($joiner, map { "ss$_" } sort {$a <=> $b} @ssids);
    $xml .= &singleChild('alias', $aliAttr, $pad);
    my $pad2 = " $pad";
    my $ssxml = "";
    foreach my $mapmeth (sort keys %methods) {
        $ssxml .= &tagOpen( 'ssM', { method => $mapmeth }, $pad)."\n";
        foreach my $auth (sort keys %{$methods{$mapmeth}}) {
            $ssxml .= &tagOpen( 'ssA', { auth => $auth }, $pad2);
            foreach my $ssattr (@{$methods{$mapmeth}{$auth}}) {
                $ssxml .= &singleChild('ss', $ssattr);
            }
            $ssxml .= &tagClose( 'ssA');
        }
        $ssxml .= &tagClose( 'ssM', $pad);
    }
    print SIMPSSX $ssxml;
    return $xml;
}

sub simple_dates {
    my ($root, $attr) = @_;
    my ($cd, $cb) = &extract_date( $root, 'Create');
    my ($ud, $ub) = &extract_date( $root, 'Update');
    if ($cd) {
        $attr->{add}  = $cd;
        $attr->{addB} = $cb;
    } else {
        &attr_err($attr, "No Create date");
    }
    if ($ud) {
        $attr->{upd}  = $ud;
        $attr->{updB} = $ub;
    } else {
        &attr_err($attr, "No Update information", "Created $cb / $cd");
    }
}

sub simp_sequence {
    my ($node, $pad) = @_;
    my $attr = {};
    my $stags = $node->{BYTAG}{Sequence} || [];
    if ($#{$stags} < 0) {
        &simp_err("No flanking sequence");
    } elsif ($#{$stags} > 0) {
        &simp_err("Multiple flanking sequences");        
    } else {
        my $stag = $stags->[0];
        foreach my $side ( 5, 3 ) {
            my $flanks = $stag->{BYTAG}{"Seq$side"} || [];
            if ($#{$flanks} > 0) {
                &simp_err("Multiple $side' flanking sequence");
            }  elsif ($#{$flanks} == 0) {
                my $seqdata = $flanks->[0]{TEXT};
                $seqdata =~ s/\s+//g;
                if ($side == 5 ) {
                    $attr->{left} = $seqdata;
                } else {
                    $attr->{right} = $seqdata;
                }
            }
        }
        &simple_observed_to_bases($stag, 'Observed', $attr);
        if (my $ex = $stag->{ATTR}{exemplarSs}) {
            $attr->{exemplar} = "ss$ex";
        }
    }
    return &singleChild('sequence', $attr, " $pad");
}

sub simp_validation {
    my ($node, $pad) = @_;
    my $attr = {};
    my @validated;
    foreach my $vald (@{$node->{BYTAG}{Validation} || []}) {
        while (my ($kname, $bool) = each %{$vald->{ATTR}}) {
            if ($bool eq 'true' && $kname =~ /^by(\S+)/) {
                push @validated, $1;
            }
        }
    }
    $attr->{TEXT} = join($joiner, @validated);
    my ($min,$max) = ($node->{ATTR}{validProbMin},
                      $node->{ATTR}{validProbMax});
    if ($min && $max) {
        $attr->{success} = "$min-$max\%";
    } elsif ($min) {
	$attr->{success} = "at least $min\%";
    } elsif ($max) {
	$attr->{success} = "at best $max\%";
    }
    return &singleChild('valid', $attr, " $pad");
}

sub attr_err {
    my $attr = shift;
    $attr->{Error}{$_[0]}++;
    &simp_err(@_);
}

sub note {
    return unless ($fc && $statFiles{NOTE});
    my @list = @_;
    map { $_ = '-UNDEF-' unless (defined $_) } @list;
    my $txt = join("\t", @list);
    $fc->write_output('NOTE',"$txt\n");
}
sub note_once {
    return unless ($fc && $statFiles{NOTE});
    my @list = @_;
    map { $_ = '-UNDEF-' unless (defined $_) } @list;
    my $txt = join("\t", @list);
    return if ($noted{$txt}++);
    $fc->write_output('NOTE',"$txt\n");
}


sub simp_err {
    my @list = @_;
    unshift @list, $nowSnp if (defined $nowSnp);
    map { $_ = '-UNDEF-' unless (defined $_) } @list;
    if ($statFiles{ERROR}) {
        $fc->write_output('ERROR',join("\t", @list). "\n");
    }
}

sub attr_and_kids {
    my $attr = shift || {};
    my $text = $attr->{TEXT} || "";
    delete $attr->{TEXT};
    my (@atts, %kids);
    foreach my $k (sort keys %{$attr}) {
        my $v = $attr->{$k};
        next if (!defined $v || $v eq '');
        if (my $r = ref($v)) {
            if ($r eq 'HASH') {
                $kids{$k} = [ sort keys %{$v} ];
            }
        } else {
            push @atts, sprintf("%s='%s'", $k, $esc->esc_xml_attr($v));
        }
    }
    $text = $esc->esc_xml($text) if ($text);
    foreach my $kt (sort keys %kids) {
        map { $text .= "<$kt>".$esc->esc_xml($_)."</$kt>" } @{$kids{$kt}};
    }
    return (join(' ', @atts), $text);
}

sub singleChild {
    my ($tag, $attr, $pad) = @_;
    my $xml = $pad || "";
    $xml .= "<$tag";
    my $term = "\n";
    if ($attr->{NoCR}) { $term = ""; delete $attr->{NoCR}; }
    my ($atTxt, $text) = &attr_and_kids( $attr );
    $xml .= " $atTxt" if ($atTxt);
    if ($text) {
        # TEXT or inner children
        $xml .= ">$text</$tag>$term";
    } elsif ($atTxt) {
        # Attributes only
        $xml .= " />$term";
    } else {
        # No text and no attributes
        return "";
    }
    return $xml;
}


sub tagOpen {
    my ($tag, $attr, $pad) = @_;
    my $xml = $pad || "";
    $xml .= "<$tag";
    my ($atTxt, $text) = &attr_and_kids( $attr );
    $xml .= " $atTxt" if ($atTxt);
    $xml .= ">$text";
    return $xml;
}

sub tagClose {
    my ($tag, $pad) = @_;
    return ($pad || "") . "</$tag>\n";
}

# snptracker/testST.pl -folder /work/Mirrors/dbsnp/XML -error snptracker/runerrors.txt

sub text_maps {
    $textMaps = {
        buildtag => {
            'NO LONGER USED'    => '',
            'Homo sapiens'      => 'NCBI%d',
            'Mus musculus'      => 'NCBIM%d',
            'Rattus norvegicus' => 'RGSC%d',
        },
        taxFromId => {
            9606  => 'Homo sapiens',
            9615  => 'Canis lupus familiaris',
            10090 => 'Mus musculus',
            10116 => 'Rattus norvegicus',
            0     => '',
        },
        taxToId   => {
            human => 9606,
            mouse => 10090,
            rat   => 10116,
            dog   => 9615,
        },
        buildmap => {
            MGSCv37 => 'NCBIM37',
            GRCh37  => 'GRCh37',
        },
        orient => {
            '+1'    =>  1,
            '-1'    => -1,
            '1'     =>  1,
            forward =>  1,
            fwd     =>  1,
            rev     => -1,
            reverse => -1,
            unknown =>  0,
            top     =>  1,
            bottom  => -1,
        },
        method => {
            'dhplc'     => 'DHPLC',
            'hybridize' => 'Hybridization',
            'computed'  => 'Computation',
            'sscp'      => 'SSCP',
            'other'     => 'NCBI Other',
            'unknown'   => 'Unknown',
            'rflp'      => 'RFLP',
            'sequence'  => 'Sequencing',
        },
        impact => {
            # http://www.ncbi.nlm.nih.gov/books/NBK21088/#ch5.ch5_4_11_1
            'cds-indel'                         => 'DEL',
            'coding'                            => 'COD',
            'coding-exception'                  => 'EXC',
            'coding-nonsynon'                   => 'NON',
            'coding-nonsynonymous'              => 'NON',
            'coding-nonsynonymous-frameshift'   => 'DEL',
            'coding-nonsynonymous-missense'     => 'NON',
            'coding-nonsynonymous-nonsense'     => 'STP',
            'coding-sequence-variant'           => 'COD',
            'coding-synon'                      => 'SYN',
            'coding-synonymous'                 => 'SYN',
            'coding-unknown'                    => 'COD',
            'downstream-variant-500B'           => 'GEN',
            'downstream-variant-5KB'            => 'GEN',
            'exception'                         => 'EXC',
            'frameshift-variant'                => 'DEL',
            'intron'                            => 'INT',
            'intron-variant'                    => 'INT',
            'locus-region'                      => 'LOC',
            'missense'                          => 'NON',
            'mrna-utr'                          => 'UTR',
            'nc-transcript-variant'             => 'EXN',
            'neargene-3'                        => 'LOC',
            'neargene-5'                        => 'LOC',
            'non-synonymous-codon'              => 'NON',
            'splice-3'                          => 'SPL',
            'splice-5'                          => 'SPL',
            'splice-acceptor-variant'           => 'SPL',
            'splice-donor-variant'              => 'SPL',
            'splice-region-variant'             => 'SPL',
            'splice-site'                       => 'SPL',
            'stop-gained'                       => 'STP',
            'stop-lost'                         => 'STP',
            'synonymous-codon'                  => 'SYN',
            'upstream-variant-2KB'              => 'GEN',
            'upstream-variant-5KB'              => 'GEN',
            'utr-3'                             => 'UTR',
            'utr-5'                             => 'UTR',
            'utr-variant-3-prime'               => 'UTR',
            'utr-variant-5-prime'               => 'UTR',
            'near-gene-5'      => 'GEN',
            'near-gene-3'      => 'GEN',
            'frameshift'      => 'DEL',
            'stop-gain'      => 'STP',
            ''      => '',
            ''      => '',
            ''      => '',
            ''      => '',
            'reference'                         => '',
            'contig-reference'                  => '',
            'nmd-transcript-variant'            => '',
            'mature-miRNA-variant'              => '',
            'incomplete-terminal-codon-variant' => '',
            'complex-change-in-transcript'      => '',
        },
        mtImpact => {
            COD => 'Coding',
            DEL => 'Frameshift',
            EXC => 'Coding',
            GEN => 'Genomic',
            INT => 'Intronic',
            LOC => 'Genomic',
            NON => 'Nonsynonymous',
            SPL => 'Splice Site',
            STP => 'Stop',
            UTR => 'UTR',
            EXN => 'Exonic',
        },
        builds => {
            'Homo sapiens' => {
                'reference' => 'Human_Chr_%s.NCBI_%d',
            },
            'Mus musculus' => {
                'C57BL/6J' => 'Mus_musculus_Chr_%s.NCBIM_%d',
            },
            'Rattus norvegicus' => {
                'RGSC_v3.1' => 'Rattus_norvegicus_Chr_%s.NCBI_%d',
            },            
        },
        side => {
            5 => 'Left',
            3 => 'Right',
        },
        type => {
            'artifact'      => 'Determined to be experimental artifact',
            'gene-dup'      => 'Artifact of duplicated gene region',
            'duplicatesub'  => 'Duplicate submission',
            'notspecified'  => 'Withdrawn without stated reason',
            'ambiguousloc'  => 'Excessive number of genomic locations',
            'lowmapquality' => 'Insufficient evidence supporting variation',
        },
        copynum => {
            "two-hits-in-contig" => "REPETITIVE",
            "less-10-hits"       => "REPETITIVE",
            "multiple-hits"      => "VRYREPETITIVE",
        },
        class => {
            'snp' => 'SNP',
            'in-del' => 'InDel',
            # 'het' => 'unknown seq composition, but observed heterozygous',
            'microsat' => 'STR',
            'named' => 'named',
            #'no-variation' => 'submission reports invariant region',
            'mixed' => 'Mixed variant',
        },
        alleles => {
            "NOVARIATION" => "This 'variant' has been labeled with the pseudoallele 'NOVARIATION', indicating no actual variation at this locus",
            'LARGE DELETION' => "Large deletion polymorphism lacking specific allele information",
            'LARGEDELETION' => "Large deletion polymorphism lacking specific allele information",
            'LARGEINSERTION' => "Large insertion polymorphism lacking specific allele information",
            'ALU' => "ALU repeat polymorphism",
        },
    };
    while (my ($tid, $tname) = each %{$textMaps->{taxFromId}}) {
        $textMaps->{taxToId}{lc($tname)} = $tid;
    }
    return $textMaps;
}

sub err {
    my $snp;
    my @list = map { defined $_ ? $_ : '-UNDEF-' } ( $cname, @_ );
    if (ref($list[-1])) {
        my $snp = pop @list;
        unshift @list, $snp->accession;
    } else {
        unshift @list, 'UNK';
    }
    map { $_ = '-UNDEF-' unless (defined $_) } @list;
    if ($statFiles{ERROR}) {
        $fc->write_output('ERROR',join("\t", @list). "\n");
    }
    # Show user the message, unless it is one of the very common problems
    my $err = $list[2];
    unless ($err =~ /(missing allele|Mis\-sized)/) {
        my $seen = ++$errWarnCount{$err};
        if ($seen >= $ewcLvl) {
            $args->msg("[ERR]", "$err : $seen events");
            $errWarnCount{$err} = 0;
        }
    }
}

sub csv {
    # Child Simple Value
    # Given a parent node and the tag name of a child, will return
    # the text of the child if there is one and only one such child
    my ($parent, $childname) = @_;
    my $list = $parent->{BYTAG}{$childname};
    return undef if (!$list || $#{$list} != 0);
    return $list->[0]{TEXT};
}

sub get_simple_item {
    my ($root) = @_;
    return $root->{ATTR}{id} || "Unknown node <".$root->{NAME}.">";
}

sub simple_observed_to_bases {
    my ($parent, $tag, $attr) = @_;
    my $allstr = &csv($parent, $tag);
    return () unless ($allstr);
    $allstr =~ s/^\s+//;
    $allstr =~ s/\s+$//;
    my @parts = split(/\s*\/\s*/, uc($allstr));
    my $rptUnit;
    my $unusual = 0;
    my %errs;
    for my $p (0..$#parts) {
        my $part = $parts[$p];
        if ($part =~ /^[ACTG]+$/) {
            # Nothing needed for 'normal' alleles:
        } elsif ($part eq 'DEL' || $part eq '-') {
            $parts[$p] = '-';
        } elsif ($part =~ /^\d+$/) {
            # A pure integer, presumably specifying a multiple of a prior unit
            if ($rptUnit) {
                # The repeat unit has already been defined
                # (CA)13/14/15/16/18
                $parts[$p] = $rptUnit x $part;
            } else {
                # Hm, an integer without a preceding repeat unit
                $errs{"Repeat count without repeat unit"}{$part} = 1;
                $parts[$p] = '';
                $unusual++;
            }
        } elsif ($part =~ /^\(([^\)]+)\)(\d+)$/) {
            # Simple repeat
            # Make note of the repeat unit in case it is needed above
            $rptUnit = $1;
            # Expand allele to literal:
            $parts[$p] = $rptUnit x $2;
        } else {
            if ($part =~ /\(/) {
                # Not a simple repeat (would have been caught above)
                if ($part =~ /^\(([^\)]+)\)$/) {
                    # Atypical allele, ie LARGEDELETION
                    my $odd = $1;
                    if ($odd =~ /^(\d+)BP$/i) {
                        # Specific length allele with no specific bases
                        $parts[$p] = 'N' x $1;
                    } else {
                        my $base = $parts[$p];
                        my $com = exists $textMaps->{alleles}{$base} ?
                            $textMaps->{alleles}{$base} : undef;
                        $com = "Large deletion polymorphism that only specifies length"
                            if ($base =~ /^\d+ BP (DEL|DELETED)$/);
                        if ($com) {
                            $attr->{Comment}{$com}++;
                        } else {
                            $errs{"Invalid allele"}{$odd} = 1;
                            $unusual++;
                        }
                        $parts[$p] = '';
                    }
                } else {
                    # Complex allele such as (CAA)6(A)7
                    $parts[$p] = &EXPAND_RPT( $part );
                    if ($parts[$p] =~ /\)/) {
                        $errs{"Failed to expand allele"}{$part} = 1;
                        $parts[$p] = '';
                        $unusual++;
                    }
                }
            }
            $rptUnit = undef;
        }
    }
    map { $attr->{Error}{$_} = 1 } keys %errs;

    # A very few entries over-specify the alleles, such that some alleles
    # are repeated in the list. This generally occurs with a numbered repeat,
    # Such as:
    # (AT)2/3/4/5/ATATAT/-
    # I have seen this when the first researcher reports 2-5 units of an AT
    # repeat, and then a second researcher reports 3 units and a deletion.
    # The result is that (AT)3 === ATATAT is duplicated. Code below is
    # designed to eliminate duplicates
    my %nonredun = map { $_ => 1 } @parts;
    foreach my $base (keys %nonredun) {
        if (length($base) > 1000) {
            $attr->{Comment}{"At least one allele over 1kb in length was ignored. SnpTracker lacks the capacity to handle such long alleles."}++;
            delete $nonredun{$base};
        }
    }
    my @bases = sort keys %nonredun;
    $attr->{Alleles} = join($joiner, @bases);
    return @bases;
}


sub observed_to_bases {
    my ($parent, $tag, $snp) = @_;
    my $allstr = &csv($parent, $tag);
    return () unless ($allstr);
    $allstr =~ s/^\s+//;
    $allstr =~ s/\s+$//;
    my @parts = split(/\s*\/\s*/, uc($allstr));
    my $rptUnit;
    my $unusual = 0;
    my %errs;
    for my $p (0..$#parts) {
        my $part = $parts[$p];
        if ($part =~ /^[ACTG]+$/) {
            # Nothing needed for 'normal' alleles:
        } elsif ($part eq 'DEL' || $part eq '-') {
            $parts[$p] = '-';
        } elsif ($part =~ /^\d+$/) {
            # A pure integer, presumably specifying a multiple of a prior unit
            if ($rptUnit) {
                # The repeat unit has already been defined
                # (CA)13/14/15/16/18
                $parts[$p] = $rptUnit x $part;
            } else {
                # Hm, an integer without a preceding repeat unit
                $errs{"Repeat count without repeat unit"}{$part} = 1;
                $parts[$p] = '';
                $unusual++;
            }
        } elsif ($part =~ /^\(([^\)]+)\)(\d+)$/) {
            # Simple repeat
            # Make note of the repeat unit in case it is needed above
            $rptUnit = $1;
            # Expand allele to literal:
            $parts[$p] = $rptUnit x $2;
        } else {
            if ($part =~ /\(/) {
                # Not a simple repeat (would have been caught above)
                if ($part =~ /^\(([^\)]+)\)$/) {
                    # Atypical allele, ie LARGEDELETION
                    my $odd = $1;
                    if ($odd =~ /^(\d+)BP$/i) {
                        # Specific length allele with no specific bases
                        $parts[$p] = 'N' x $1;
                    } else {
                        my $base = $parts[$p];
                        my $com = exists $textMaps->{alleles}{$base} ?
                            $textMaps->{alleles}{$base} : undef;
                        $com = "Large deletion polymorphism that only specifies length"
                            if ($base =~ /^\d+ BP (DEL|DELETED)$/);
                        if ($com) {
                            $snp->add_comment($com);
                        } else {
                            $errs{"Invalid allele"}{$odd} = 1;
                            $unusual++;
                        }
                        $parts[$p] = '';
                    }
                } else {
                    # Complex allele such as (CAA)6(A)7
                    $parts[$p] = &EXPAND_RPT( $part );
                    if ($parts[$p] =~ /\)/) {
                        $errs{"Failed to expand allele"}{$part} = 1;
                        $parts[$p] = '';
                        $unusual++;
                    }
                }
            }
            $rptUnit = undef;
        }
    }
    $snp->{_ALLELE_ERROR_} = \%errs;
    my @eTypes = sort keys %errs;
    unless ($#eTypes == -1) {
        my @bits = ("Errors parsing observed alleles", substr($allstr,0,60));
        foreach my $et (@eTypes) {
            push @bits, "$et : ".join(", ", sort keys %{$errs{$et}});
        }
        &err(@bits, $snp);
    }
    # A very few entries over-specify the alleles, such that some alleles
    # are repeated in the list. This generally occurs with a numbered repeat,
    # Such as:
    # (AT)2/3/4/5/ATATAT/-
    # I have seen this when the first researcher reports 2-5 units of an AT
    # repeat, and then a second researcher reports 3 units and a deletion.
    # The result is that (AT)3 === ATATAT is duplicated. Code below is
    # designed to eliminate duplicates
    my %nonredun = map { $_ => 1 } @parts;
    foreach my $base (keys %nonredun) {
        if (length($base) > 1000) {
            $snp->add_comment
                ("At least one allele over 1kb in length was ignored. ".
                 "SnpTracker lacks the capacity to handle such long alleles.");
            delete $nonredun{$base};
        }
    }
    delete $nonredun{''};
    my @alleles = sort keys %nonredun;
    &err("$unusual unusual alleles",  substr(join('/', @alleles),0,60), $snp)
        if ($unusual);
    return @alleles;
}

sub EXPAND_RPT {
    # Designed to turn:
    # (CA)13CCCCATCTA(TATC)3(TCTG)4
    # into:
    # CACACACACACACACACACACACACACCCCATCTATATCTATCTATCTCTGTCTGTCTGTCTG
    my ($string) = @_;
    while ($string =~/\(([^\)]+)\)(\d+)/) {
        my ($rpt, $num) = ($1, $2);
        my $found   = "($rpt)$num";
        my $replace = $rpt x $num;
        $string =~ s/\Q$found\E/$replace/g;
    }
    return $string;
}

sub clean_allele {
    my ($val, $is_rev) = @_;
    return "" unless ($val);
    $val =~ s/^\s+//; $val =~ s/\s+$//;
    if ($val =~ /\(/) {
        $val = &EXPAND_RPT($val);
    }
    return $val if ($val eq '(INDETERMINATE)' || 
                    $val eq '(HETEROZYGOUS)');
    $val = '-' if ($val eq 'DEL');
    $val = $su->revcom($val) if ($is_rev);
    return $val;
}

sub extract_date {
    my ($root, $tag) = @_;
    my $list = $root->{BYTAG}{$tag};
    return () unless ($list && $#{$list} == 0);
    my $dtag = $list->[0];
    my $dt   = $dtag->{ATTR}{date};
    return ($dt, $dtag->{ATTR}{build});
}

sub note_verbose {
    my ($lvl, $msg) = @_;
    return if ($vb < $lvl || !$msg);
    if ($fc) {
        $fc->write_output('VERB', $msg);
    } else {
        warn $msg;
    }
}

sub mrna_orientation {
    my ($gname, $accU) = @_;
    return undef unless ($gname && $accU);
    unless (defined $cache->{RNA_Orient}{$accU}{$gname}) {

        $args->death("Need to replace mechanism for getting mRNA orientation",
                     "$accU on $gname");

        my %strands = map {$_ => 1} (die 'array of strands goes here');
        my @str = keys %strands;
        $cache->{RNA_Orient}{$accU}{$gname} = ($#str == 0) ? $str[0] : 0;
        # warn "$accU $gname => " . $cache->{RNA_Orient}{$accU}{$gname};
    }
    return $cache->{RNA_Orient}{$accU}{$gname};
}

sub set_taxa_from_file_name {
    my $file = shift;
    my @bits = split(/\//, $file);
    pop @bits;
    while (my $sdir = pop @bits) {
        if ($sdir =~ /^[a-z]+_(\d+)$/) {
            return &set_taxa( { ATTR => { taxId => $1 } } );
        }
    }
    $args->death("Failed to identify TaxID component in file path",
                 $file);
}

sub set_taxa {
    my ($hash) = @_;
    my $tid = $hash->{ATTR}{taxId};
    if ($tid) {
        unless ( $specName = $textMaps->{taxFromId}{$tid} ) {
            $args->death("Failed to resolve NCBI TaxID $tid to binomial name",
                         "You can update \$textMaps->{taxFromId}");
        }
        my $lspec = lc($specName); $lspec =~ s/\s/_/g;
        $cfrm     = "$lspec.chromosome.%s.%s";
    } else {
        &err("Failed to find TaxaID");
    }
}

sub extract_build {
    my $gzf = shift;
    my $line = `gunzip -c $gzf | head -n 10 |  grep -i dbSnpBuild` || "";
    if ($line =~ /dbSnpBuild(No)?=[\"\'](\d+)[\"\']/i) {
        return $2;
    }
    return 0;
}

sub get_build {
    my ($node) = @_;
    unless ($globalBuild) {
        my $par = $node->{PARENT};
        $globalBuild  = $par->{ATTR}{dbSnpBuild};
        if ($globalBuild) {
            $args->msg("$specName dbSNP $globalBuild") if ($fc->child == 1);
        } else {
            $args->death("Failed to recover build information from ". 
                         $par->{NAME});
        }
    }
    return $globalBuild;
}

sub clean_string {
    my $desc = shift;
    $desc =~ s/<[^>]+>/ /g;
    $desc =~ s/\,/\, /g;
    $desc =~ s/\s+/ /g;
    $desc =~ s/\s$//;
    $desc =~ s/^\s//g;
    return $desc =~ /^\s*$/ ? "" : $desc;
}

sub unescape {
    my ($string) = @_;
    while ($string =~ /\%([0-9A-F]{2})/i) {
        my $hex = $1;
        my $char = chr( hex($hex) );
        $string =~ s/\%$hex/$char/g;
    }
    return $string;
}

sub init_maploc {
    my $bld = $globalBuild || $args->val(qw(build));
    # $args->death("Can not initialize MapLoc without build information")
    #    unless ($bld);
    $mapLoc = &get_maploc();
    $mapLoc->build($bld);
    &set_meta_stuff();
}

sub set_meta_stuff {
    my $setXtra = shift || $args->val(qw(addmeta dometa setmeta));

    my $anonPop = $mapLoc->get_population("Unvalidated dbSNP", $stAuth);
    $anonPopId  = $anonPop->pkey();

    $locCatName = "dbSNP Polymorphisms";
    my $cat = $mapLoc->get_text( $locCatName );
    $locCategory = $cat->pkey();

    $rootPop = $mapLoc->get_population("dbSNP Populations", $stAuth);

    return unless ($setXtra);

    $cat->tag("Description", 
              "Sites from the dbSNP reference polymorphism database at the NCBI ");
    $cat->tag("Short", "dbSNP");
    $cat->tag("Useful Tag", 'Clinical Significance');
    $cat->tag("SNP Class Tag", "Population Class");
    $cat->write();

    my $catPar = $mapLoc->get_text( "Polymorphism Categories" );
    $catPar->tag("Description", "Variant locations attributed to germ line polymorphism across populations");
    $catPar->tag("Member", $cat->text());
    $catPar->write();

    $rootPop->tag_values("Description", "This is the root population for all dbSNP populations");
    $rootPop->update();

    my $unkPop = $mapLoc->get_population("Unclassified dbSNP", $stAuth);
    $unkPop->tag_values("Description",
                        "Parent for populations without PopClass");
    $unkPop->parent($rootPop);
    $unkPop->tag_values("Category", $locCatName);
    $unkPop->update();
    $unkPop->clear_tag("Description",
                        "This is a pseudo-population that allows alleles to be recorded for a variant in the absence of actual population assays. Locations flagged ONLY with this population should be considered as potentially non-variant.");

    $anonPop->parent( $unkPop );
    $anonPop->tag_values("Description",
                        "This is a pseudo-population that allows alleles to be recorded for a variant in the absence of actual population assays. Locations flagged ONLY with this population should be considered as potentially non-variant.");
    $anonPop->tag_values("Category", $locCatName);
    $anonPop->update();
}

sub get_maploc {
    my $mapLoc = BMS::SnpTracker::MapLoc->new
        ( # -build  => $bld,
          -noenv    => $args->val(qw(noenvironment noenv)),
          -instance => $dbInst,
          # -makedb => $args->val(qw(makedatabase makedb)),
          );
    return $mapLoc;
}

sub finalize {
    foreach my $err (sort keys %errWarnCount) {
        $args->msg("[ERR]", "$err : $errWarnCount{$err} events");
    }
    %errWarnCount = ();
    # exit;
}

