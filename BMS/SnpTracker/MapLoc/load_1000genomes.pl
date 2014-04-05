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


use strict;
use BMS::SnpTracker::MapLoc;

use BMS::TableReader;
use BMS::Utilities::SequenceUtilities;
use BMS::Utilities::Escape;

use BMS::ArgumentParser;
use BMS::ForkCritter;

=head1 Notes

Chr2  = 3,307,622 rows
Chr20 =   855,196 rows
Chr22 =   494,358 rows

=cut

my ($mapLoc, $fc, $tr, $locCategory, $samp2popid, $samp2pname, $rootPopId,
    $lastPos, $simpFH, $popCols, %textCache, %tkeyCache);

my $build = 'GRCh37';

my $args = BMS::ArgumentParser->new
    ( -nocgi          => $ENV{'HTTP_HOST'} ? 0 : 1,
      -chr            => "",
      -testmode       => 'quiet',
      -fork           => 1,
      -verbose        => 1,
      -progress       => 300,
      -paramfile      => 'BMS::SnpTracker::MapLoc',
      -srcdir         => '/tmp',
      );

$args->shell_coloring();
$args->debug->maxany(50);
$args->debug->skipkey([qw(GENOTYPE_NAMES)]);

my $debug     = $args->val(qw(debug));
my $limit     = $args->val(qw(limit)) || 0;
my $progress  = $args->val(qw(progress)) || 0;
my $isTrial   = $args->val(qw(trial istrial));
my $vb        = $isTrial || $args->val(qw(vb verbose)) || 0;
my $mlInt     = $args->val(qw(instance));
my $forkNum   = $args->val(qw(forknum fork)) || 1;
my $srcDir    = ($args->val(qw(srcdir)) || '.')."/1000genomes";
my $benchDir  = "$srcDir/benchmarks";
my $mirdir    = $args->val(qw(mirrordir)) || "20110521";
my $simpDir   = "$srcDir/Simple";
my $deferCache = $args->val(qw(defer defercache defersize)) || 300;

my $prfx       = "1kG";
my $rootName   = "$prfx.All";
my $locCatName = "1000 Genomes Polymorphisms";
my $stAuth     = "1000 Genomes";

$args->assure_dir($simpDir);
my $dataDir   = &mirror();

if ($debug) {
    my ($min, $max) = (0,99999999999999);
    if ($debug =~ /<(\d+)/) {
        $max = $1;
    } elsif ($debug =~ />(\d+)/) {
        $min = $1;
    } elsif ($debug =~ /(\d+)[\-\.]+(\d+)/) {
        ($min, $max) = ($1, $2);
    } elsif ($debug =~ /(\d+)/) {
        ($min, $max) = ($1, $1);
    }
    $debug = [ $min, $max ];
}

$args->msg("[INSTANCE]", $mlInt);

if (my $file = $args->val(qw(sample patient population))) {
    &load_sample_info($file);
} elsif ($args->val(qw(simplify))) {
    &simplify();
} elsif ($args->val(qw(load process))) {
    &process();
} else {
    &usage();
    exit;
}

$mapLoc->dbh->disconnect() if ($mapLoc);

$args->msg("[DEBUG]", "Collected Parent+child benchmarks:");
print STDERR &benchmark_table( $fc || $mapLoc );

sub benchmark_table {
    my $obj = shift;
    return "" unless ($obj);
    return $obj->showbench( -shell => 1, -minfrac => 0.0001 );
}

sub usage {
    my $txt =<<USAGE;

#1 - Download and parse Sample Excel file

Before you load data, you need to provide the system with information
about the samples (patients), which will be used to make MapLoc
Populations for each. There are two problems here: First, the data are
stored in an Excel workbook, which increases the likelihood that this
parser will break at some point in the future (I have already
rewritten it once). Second, there does not appear to be a consistent
place to find the workbook. I was once able to find a link here:

    http://www.1000genomes.org/about#ProjectSamples

... which directed me here:

    ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/

... but be aware that the above link will eventually be supplanted by
newer information. Save the Excel file locally, and then run the
program with:

    -sample <path_to_xlsx_file>

That will load the database with metadata for the high-level
populations, as well as Sample-to-Population mappings that will be
required for parsing.

# 2 - Download VCF files

You will next need to provide VCF files to the program. By default, it
will use /tmp/, but you can change that with the -srcdir parameter if
desired. Inside the source directory, it will look in the
'1000genomes' subdirectory (which will be automatically created if not
already present). It is here that you should place your 1000 Genomes
VCF files.

This can be done automatically by running this script with the -mirror
parameter. It will use wget to retrieve files from the NCBI. Note that
the default release is set as '$mirdir', but you can override this
with the -mirrordir parameter. You can of course also download the
files through any other means, just be sure to place them in the
srcdir 1000genomes folder.

The program is expecting a particular nomenclature for the files:

   ALL.chr([XY0-9]+).+.genotypes.vcf.gz

Any file not matching this pattern will be ignored. If 1000 Genomes
changes their naming convention, this will break the loader until the
RegExp pattern is changed.

# 3 - Simplify VCF files

VCF files need to be "simplified" before being processed. Do so by
running the program with the -simplify parameter after you have
downloaded the VCF files you want (presumably all chromosomes in the
genome). This operation will shrink the files in preparation for
loading. If you later download new VCF files you can repeat
simplification without affecting previously generated files.

The simplification process is forked, use -fork to specify the number
of children you want working on the task (default is just 1). Each
child should effectively monopolize a CPU. A 3GHz CPU is capable of
simplifying a little over 6000 entries per minute, so Chr2 (about 3M
entries) will take around 8 hours.

After the "-Simple.tsv" files have been created, you may delete their
corresponding VCF files to save disk space. Version 20110521 occupied
143 Gb of disk space in gzipped VCF format, but only 6 Gb after
simplification.

# 4 - Load variations into MapLoc

You can now begin loading the variants into the database. Run this
script with the -load parameter. Like simplification, the process is
forked, and the default of a single child can be changed by passing
the number of CPUs you wish to use with the -fork parameter. Children
will consume only about 160 Mb each, and little CPU - the process
seems to be I/O limitted. Running -fork 30 on a consumer-grade 3 GHz
quad core system had only a minor impact on CPU and RAM utilitzation,
with most children being near 0% CPU. The load starts slow then speeds
up somewhat as it progresses, with estimates after 10 minutes of about
a week to load chromosome 2.

After a few thousand rows have been loaded, it is a good idea to
analyze the database:

    psql -c 'VACUUM ANALYZE VERBOSE' maploc

This will refresh statistics and should improve load times.

Reloading data is not harmful, so the loader may be halted and
restarted if desired. It does not remember where it left off, but
previously loaded entries will be processed much faster as new rows
will not be inserted into the database. If you wish to load particular
chromosomes, you may specify the -chr parameter:

    -chr 7,12,Y
    -chr 13-20

This may be used to selectively load only some files if you need to
restart the process.

USAGE

    $args->msg("[-]","Usage:", split(/[\n\r]/, $txt),'');
}

sub process {
    my @files = &get_files($simpDir, '([XY0-9]+)-Simple.tsv');
    if ($#files == -1) {
        &usage();
        $args->msg("[!]", "It looks like you need to run -simplify");
        return;
    }
    $args->msg("[^]","Parsing chromosomes ".
               join(',', map { $_->{chr} } @files));
    &samp_to_pop();
    $fc = BMS::ForkCritter->new
        ( -init_meth   => \&init_maploc,
          -method      => \&load_simple_vcf,
          -finish_meth => \&finalize,
          -last_item_meth => \&get_simple_item,
          -limit       => $limit,
          -progress    => $progress,
          -verbose     => $vb );

    foreach my $dat (@files) {
        my $chr  = $dat->{chr};
        my $path = $dat->{path};
        $args->msg("[<]","Parsing $chr", $dat->{short});
        if ($isTrial) {
            $args->msg_once("[-]","Trial mode, no loads occuring");
            next;
        }
        $fc->reset();
        $fc->input_type('tsv header hash');
        $fc->input($path);
        $args->assure_dir($benchDir);
        foreach my $bfile ($args->read_dir( -dir => $benchDir,
                                            -keep => '.tsv$' )) {
            unlink($bfile);
            $args->msg_once("[!!]","Could not remove prior benchmark files",
                            "Benchmark calculations may be inaccurate",
                            $benchDir, ) if (-f $bfile);
        }
        
        $forkNum ||= 1;
        if (my $failed = $fc->execute( $forkNum )) {
            $args->death("$failed jobs failed to properly convert");
        }
        foreach my $bfile ($args->read_dir( -dir => $benchDir,
                                            -keep => '.tsv$' )) {
            $fc->benchmarks_from_file( $bfile, $forkNum );
        }
    }
}

sub simplify {
    my @files = &get_files();
    if ($#files == -1) {
        &usage();
        $args->msg("[!]", "It looks like you need to run -mirror");
        return;
    }
    $args->msg("[^]","Simplifying chromosomes ".
               join(',', map { $_->{chr} } @files));
    &samp_to_pop();
    $fc = BMS::ForkCritter->new
        ( -inputtype   => 'array',
          -init_meth   => \&init_maploc,
          -method      => \&simplify_vcf,
          -finish_meth => \&finalize,
          -last_item_meth => \&get_simple_item,
          -verbose     => $vb );
    $build = 'GRCh37';
    $fc->input(\@files);
    if (my $failed = $fc->execute( $forkNum || 1)) {
        $args->death("$failed jobs failed to properly convert");
    }
    
}

sub samp_to_pop {
    unless ($samp2popid) {
        my $ti  = time;
        my $ml  = &get_maploc();
        my $dat = $ml->tag_query( -tag => "In Population", -explain => 0 );
        $samp2popid = {};
        my %pops;
        foreach my $row (@{$dat}) {
            my ($sid, $tid, $vid) = @{$row};
            my $mem   = $ml->pkey_to_text($sid);
            my $pname = $ml->cached_pkey_to_text( $vid );
            my $pop   = $pops{$vid} ||= $ml->get_population($pname, $stAuth);
            $samp2popid->{$mem}     = $pop->pkey();
            $samp2popid->{$pname} ||= $pop->pkey();
            $samp2pname->{$mem}     = $pname;
        }
        $args->msg("[-]","Loaded population map (Sample -> Population)",
                   ($#{$dat} + 1)." records in ".(time - $ti)." seconds");
        my $rootPop = $ml->get_population($rootName, $stAuth);
        $rootPopId = $rootPop->pkey();
        $samp2popid->{$rootName} ||= $rootPopId;
        
    }
    return $samp2popid;
}

sub err {
    $args->msg("[".$tr->row_count()."]", @_);
}

sub load_simple_vcf {
    my $row = shift;
    my ($chr, $lft, $rgt, $bld) = map { $row->{$_} } qw(Chr Lft Rgt Build);
    my $loc = $mapLoc->flanks_to_location( $chr, $lft, $rgt, $bld );
    $loc->add_category( $locCategory );
    my @bases = split(/\|/, $row->{Alleles});
    while (my ($pname, $numTxt) = each %{$row}) {
        if (my $pid = $samp2popid->{$pname}) {
            my @nums = split(/\|/, $numTxt);
            my $depth = 0; map { $depth += $_ } @nums;
            for my $c (0..$#nums) {
                if (my $num = $nums[$c]) {
                    $loc->add_allele_to_db
                        ($bases[$c], $pid, $num / $depth, $depth, 1 );
                }
            }
        } elsif ($pname =~ /1kG/) {
            $args->msg_once("Failed to recover pop_id for '$pname'");
        }
    }
    if (my $tags = $row->{Tags}) {
        foreach my $tvt (split(',', $tags)) {
            my @tv = map { $tkeyCache{$_} ||= $mapLoc->
                               cached_pkey_to_text($_) } split('=', $tvt);
            $loc->tag_values(@tv);
        }
    }
    if ($debug) {
        print $loc->to_text();
    } else {
        $mapLoc->add_deferred_alleles() 
            if ($mapLoc->deferred_allele_count() > $deferCache);
        $loc->update();
    }
}

sub simplify_vcf {
    my $dat  = shift;
    my $chr  = $dat->{chr};
    my $path = $dat->{path};
    my $file = sprintf("%s/%s-Simple%s.tsv", $simpDir, $dat->{chr},
                       $limit ? "-Limit" : "");
    if (-s $file) {
        $args->msg("[^]","Already simplified: $chr", $dat->{short});
        return;
    }
    $args->msg("[<]","Parsing $chr", $dat->{short});
    $tr = BMS::TableReader->new
        ( -format    => 'vcf',
          -limit     => $limit );
    $tr->input($path);
    $tr->select_sheet(1);
    &init_simp($file, $tr);
    # http://www.1000genomes.org/node/101
    $lastPos = 0;
    my $ti = time;
    my $start = $ti;
    while (my $row = $tr->next_hash()) {
        if (my $halt = &simplify_row($row)) {
            $args->msg("[HALT]", $halt);
            last;
        }
        next unless ($fc->child == 1);
        my $t2 = time;
        if ($t2 - $ti > $progress) {
            my $elaps = $t2 - $start;
            my $num   = $tr->row_count();
            my $rate  = 10 * int(0.5 + $num * 60 / $elaps) / 10;
            $args->msg("[TIME]", sprintf
                       ("%8d rows %d/min %s:%s", $num, $rate, $row->{CHROM}, 
                        $mapLoc->comma_number($row->{POS})));
            $ti = $t2;
        }
    }
    &finish_simp();
}

sub init_simp {
    my ($file, $tr) = @_;
    my $usedMems;
    if ($tr) {
        $usedMems = { map { $_ => 1 } @{$tr->{FILEDAT}{GENOTYPE_NAMES}} };
    }
    my %uPop;
    # Find the 1000 genomes populations:
    while (my ($mem, $pname) = each %{$samp2pname}) {
        next if ($usedMems && !$usedMems->{$mem});
        $uPop{$pname} = 1 if ($pname =~ /^1kG\./);
    }
    # Header will start with basic location data, then populations
    my @head = ('Chr', 'Lft', 'Rgt', 'Build', $rootName);
    push @head, sort keys %uPop;
    # End with alleles and finally tags
    push @head, ('Alleles', 'Tags');

    for my $i (0..$#head) {
        $popCols->{$head[$i]} = $i;
    }
    # Now also map the individual members to the relevant columns
    while (my ($mem, $pname) = each %{$samp2pname}) {
        if (my $ind = $popCols->{$pname}) {
            $popCols->{$mem} = $ind;
        }
    }
    open($simpFH, ">$file") || $args->death
        ("Failed to write simplified file", $file, $!);
    print $simpFH join("\t", @head)."\n";
}

sub finish_simp {
    close $simpFH;
}


sub simplify_row {
    my $row = shift;
    my $filt = $row->{FILTER} || "?";
    unless ($filt eq 'PASS') {
        $args->msg_once("Ignoring FILTER = '$filt'");
        return "";
    }
    my ($chr, $pos) = ($row->{CHROM}, $row->{POS});
    if ($debug) {
        return "" if ($pos < $debug->[0]);
        return sprintf("Analyzed positions %d-%d", @{$debug})
            if ($pos > $debug->[1]);
    }

    unless ($chr && $pos =~ /^\d+$/) {
        &err(sprintf("Unclear position : CHR:%s POS:%s", 
                     $chr || '?', $pos || '?'));
        return "";
    }
    # warn "$chr:$pos" if ($pos == $lastPos);
    $lastPos = $pos;

    my (@alleles, %info, %snpTags);
    if (my $v = $row->{REF}) {
        $alleles[0] = $v;
    } else {
        &err("No reference allele $chr:$pos");
        return "";
    }
    if (my $v = $row->{ALT}) {
        my @bases = split(/\,/,$v);
        for my $j (0..$#bases) {
            $alleles[$j+1] = $bases[$j];
        }
    } else {
        &err("No alternate allele $chr:$pos");
        return "";
    }
    for my $i (1..$#alleles) {
        if ($alleles[$i] =~ /^<(.+)>$/) {
            my $tok = $1;
            $args->msg_once("Ignoring imprecise alterations: <$tok>");
            return "";
            if ($tok eq 'DEL') {
                # Is this right?
                $alleles[$i] = "";
            } else {
                $args->msg_once("Skipping <$tok> entry");
                return "";
            }
        }
    }
    foreach my $ibit (split(/\;/, $row->{INFO} || "")) {
        if ($ibit =~ /^(\S+)=(.+)$/) {
            my ($k, $v) = ($1, $2);
            if ($k =~ /^(([A-Z]+)_)?AF$/) {
                my $pop = $2 || "Global";
                $snpTags{"Allele Frequency - $pop"} = $v;
            }
            $info{$k} = $v;
        }
    }
    if (my $qual = $row->{QUAL}) {
        if ($qual =~ /^\d+$/ && $qual < 100) {
            $snpTags{"Reduced Quality"} = $qual;
        }
    }
    
    my $l   = $pos - 1;
    while (1) {
        # What is the first base of each alternative allele?
        my %lftBases;
        for my $ai (1..$#alleles) {
            my $al = $alleles[$ai];
            # Use the index number for empty alleles
            $lftBases{ $al ? substr($al, 0 , 1) : $ai }++;
        }
        my @seen = keys %lftBases;
        unless ($#seen == 0 && $seen[0] eq substr($alleles[0], 0, 1)) {
            # There are either multiple first bases, or a single one
            # that is not the same as the first reference base
            last;
        }
        # All alleles start with the same base. Chop it off
        map { $alleles[$_] = substr($alleles[$_], 1) } (0..$#alleles);
        # Also increment the left flank
        $l++;
    }
    my $len = length($alleles[0]);
    my $r   = $l + $len + 1;

    my %counts;
    # warn $args->branch(-ref => $row, -maxany => 10);
    my $rootInd = $popCols->{$rootName};
    foreach my $dat (@{$row->{GENOTYPES} || []}) {
        if (my $ais = $dat->{GT}) {
            if (my $popInd = $popCols->{$dat->{name} || ""}) {
                my @alleles = split(/[\|\/]/, $ais);
                foreach my $ind ($rootInd, $popInd) {
                    map { $counts{$ind}[$_]++ } @alleles;
                }
            } else {
                $args->msg_once("[-]", "Failed to recognize population '$dat->{name}'"); 
            }
        } else {
            # $args->msg_once("No genotypes for '$dat->{name}'"); 
        }
    }
    my @sRow = ($chr, $l, $r, $build);

    while (my ($ind, $cdat) = each %counts) {
        $sRow[$ind] = join('|', map { $_ || 0 } @{$cdat});
    }
    # Alleles go next (so less -S is still useful with long alleles)
    $sRow[ $popCols->{Alleles} ] = join("|", map { $_ || '-' } @alleles);
    # and finally tags, which are devolved to normtxt pkeys:
    my @pkt;
    while (my ($tag, $val) = each %snpTags) {
        push @pkt, join('=', map { 
            $textCache{$_} ||= $mapLoc->text_to_pkey( $_ ) } ($tag, $val));
    }
    $sRow[ $popCols->{Tags} ] = join(",", @pkt);
    # Note that using pkeys means the simplified files are not portable to 
    # other MapLoc instances.
    print $simpFH join("\t", map { defined $_ ? $_ : "" } @sRow)."\n";
    return 0;
}

sub get_simple_item {
    my ($line) = @_;
    return sprintf("%s:%s^%s", split(/\t/, $line));
}

sub get_vcf_item {
    my ($dat) = @_;
    return $dat->{chr};
}


sub init_maploc {
    my $bld = $build || $args->val(qw(build));
    $args->death("Can not initialize MapLoc without build information")
        unless ($bld);
    $mapLoc = &get_maploc();
    $mapLoc->build($bld);
}

sub get_maploc {
    my $mapLoc = BMS::SnpTracker::MapLoc->new
        ( # -build  => $bld,
          -instance => $mlInt,
          -noenv  => $args->val(qw(noenvironment noenv)),
          -makedb => $args->val(qw(makedatabase makedb rebuild)),
          );
    my $cat = $mapLoc->get_text( $locCatName );
    $locCategory = $cat->pkey();
    return $mapLoc;
}

sub initialize {
    
}

sub finalize {
    my $cn = $fc->child;
    if ($mapLoc) {
        $mapLoc->add_deferred_alleles();
        if ($cn == 1) {
            $args->msg("[DEBUG]", "Forked child benchmarks:",
                       $mapLoc->_allele_update_status());
            print STDERR &benchmark_table( $mapLoc );
        }
        $mapLoc->benchmarks_to_file( "$benchDir/child-$cn.tsv" );
        $mapLoc->dbh->disconnect();
    }
#    foreach my $err (sort keys %errWarnCount) {
#        $args->msg("[ERR]", "$err : $errWarnCount{$err} events");
#    }
#    %errWarnCount = ();
    # exit;
}


sub get_files {
    my $dir    = shift || $dataDir;
    my $regExp = shift || 'ALL.chr([XY0-9]+).+.genotypes.vcf.gz$';
    my $cReq;
    if (my $r = $args->val(qw(chr chrs))) {
        $cReq = {};
        foreach my $bit ( split(/[^XY0-9\-]+/, uc($r)) ) {
            if ($bit =~ /^(\d+)\-(\d+)$/) {
                map { $cReq->{$_} = 1 } ($1..$2);
            } else {
                $cReq->{$bit} = 1;
            }
        }
        # die $args->branch($cReq);
    }
    opendir(TMPDIR, $dir) || $args->death
        ("Failed to read contents of data directory",$dir, $!);
    $args->msg("[+]","Recovering VCF files", $dir);
    my @found;
    foreach my $file (readdir TMPDIR) {
        if ($file =~ /$regExp/) {
            my ($chr) = (uc($1));
            next if ($cReq && ! $cReq->{$chr});
            my $srt = $chr =~ /^\d+$/ ? 
                sprintf("%03d", $chr) : sprintf("%3s", $chr);
            my $path = "$dir/$file";
            push @found, {
                chr   => $chr,
                path  => $path,
                short => $file,
                sort  => $srt,
                size  => -s $path,
            };
        }
    }
    if ($#found == -1) {
        $args->msg("[?]", "No files found matching $regExp", $dir);
    }
    return sort { $b->{size} <=> $a->{size} } @found;
}

sub mirror {
    my $ftp     = "ftp-trace.ncbi.nih.gov";
    my $subdir  = "1000genomes/ftp/release/$mirdir";
    my $dataDir = "$srcDir/$ftp/$subdir";
    return $dataDir unless ($args->val(qw(update mirror wget)));
    my $log = "$srcDir/wget.log";
    my $cmd = sprintf("wget -t 45 -r -A '*.gz' \\\n  ".
                      "-P \"%s\" -o \"%s\" \\\n  ftp://%s/%s", 
                      $srcDir, $log, $ftp, $subdir);
    $args->msg("[+]", "Mirroring '$mirdir' for 1000 genomes",
               $log);
    $args->msg("[CMD]", $cmd);
    system( $cmd );
    $args->msg("[-]", "Finished", $dataDir);
    return $dataDir;
}

sub OLD_parse_sample_excel {
    my $tr = shift;
    # The format of the sample file has changed
    die "Do not call this method unless you know the sample file has changed back!";

    my $cols;
    my %pops;

    $tr->select_sheet(1);
    while (my $row = $tr->next_clean_row()) {
        if ($cols) {
            my %data = map { $cols->[$_] => $row->[$_] } (0..$#{$row});
            my $id = $data{"Coriell Sample ID"};
            my $pop = $data{"Population"};
            next unless ($id && $pop);
            my @bits;
            if (my $g = $data{"Gender"}) {
                push @bits, $g; }
            if (my $c = $data{"Whole Genome Center For full project"}) { 
                push @bits, "sequenced at $c"; }
            $pops{$pop}{Member}{$id} = join(', ', @bits);
        } else {
            # Have not yet found the data
            if ($row->[0] && $row->[0] eq 'Population') {
                # We seem to have found the patient column header
                $cols = $row;
            } else {
                my ($ignore, $k, $v) = @{$row};
                if ($k && length($k) == 3) {
                    $pops{$k}{Description}{$v} = 1;
                }
            }
        }
    }
    return \%pops;
}

sub _parse_sample_excel {
    my $tr = shift;
    my %pops;


    # Need to replace this with the .ped files
    # ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20110521/supporting/phase1_samples_integrated_20101123.ped
    # The above does not have descriptions for the high-level populations,
    # though. They can be found on an HTML page as well:
    # http://www.1000genomes.org/category/frequently-asked-questions/population
    # ... but that would also require special parsing

    # Thanks to Lester Hui for finding both the above links

    
    my @captureCols = ('Platform', 'Center', 
                       'Family ID', 'Gender', 'Relationship');
    my %splitCol = map { $_ => 1 } ('Platform', 'Center');
    foreach my $sname ('Sample Info', 'Phase1', 'Final Phase Sequence Data') {
        $tr->select_sheet($sname);
        my $cols   = $tr->next_clean_row();
        if ($cols->[0] ne "Sample") {
            # Extra header row
            $cols   = $tr->next_clean_row();
            if ($cols->[0] ne "Sample") {
                $args->msg("[!!]", "Failed to find Sample header in first column", "Worksheet $sname");
                next;
            }
        }
        while (my $row = $tr->next_clean_row()) {
            my %data;
            for my $i (0..$#{$row}) {
                my $col = $cols->[$i] || "Column $i";
                my @v   = $splitCol{$col} ?
                    split(/\s*,\s*/, $row->[$i]) : ($row->[$i]);
                map { $data{$col}{$_} = 1 } @v;
            }
            while (my ($k, $v) = each %data) {
                delete $v->{""};
                $data{$k} = [ sort keys %{$v} ];
            }
            my $id  = $data{"Sample"};
            my $pop = $data{"Population"};
            next unless ($id && $pop);
            if ($#{$id} == 0 && $#{$pop} == 0) {
                $id = $id->[0];
                $pop = $pop->[0];
            } else {
                $args->msg("[!!]","Non-unique Sample and Population",
                           "Row ".$tr->rowcount()." = ".(join("+",@{$id}) || "??").
                           " and ".(join("+",@{$pop}) || "??"));
                next;
            }
            foreach my $cap (@captureCols) {
                if (my $arr = $data{$cap}) {
                    map { $pops{$pop}{MemHash}{$id}{$cap}{$_} = 1 } @{$arr};
                }
            }
            if (my $arr = $data{"Population Description"}) {
                map { $pops{$pop}{Description}{$_} = 1 } @{$arr};
            }
        }
    }

    while (my ($pop, $pH) = each %pops) {
        while (my ($id, $kv) = each %{$pH->{MemHash}}) {
            my @bits;
            if (my $v = join('+', sort keys %{$kv->{Gender}})) {
                if ($v =~ /\+/) {
                    $args->msg_once("[!!]",
                                    "Inconsistent gender $v for $id");
                } else {
                    push @bits, $v;
                }
            }
            if (my $v = join('+', sort keys %{$kv->{Center}})) {
                push @bits, "sequenced at $v";
            }
            if (my $v = join('+', sort keys %{$kv->{Platform}})) {
                push @bits, "using $v";
            }
            $pH->{Member}{$id} = join(', ', @bits);

            if (my $v = join('+', sort keys %{$kv->{'Family ID'}})) {
                if ($v =~ /\+/) {
                    $args->msg_once("[!!]",
                                    "Inconsistent family $v for $id");
                } elsif ($v ne $id) {
                    my @rel = sort keys %{$kv->{'Relationship'}};
                    if ($#rel == 0 && $rel[0] ne 'unrel') {
                        my $r = $rel[0];
                        substr($r, 0, 1) = uc(substr($r, 0, 1));
                        $pH->{Family}{$id} = "$r in $v";
                    } else {
                        $pH->{Family}{$id} = "Member in $v";
                    }
                }
            }
        }
        delete $pH->{MemHash};
    }
    $args->branch(-ref => \%pops, -maxany => 10);
    return \%pops;
}

sub load_sample_info {
    my $file = shift;
    return unless ($file);

    # It is possible that I have overlooked a more tractable file, but
    # this is the only place I was able to come across population description
    # data
    $args->msg("[+]", "Loading sample information", $file);
    $tr = BMS::TableReader->new();
    my $fmt = $tr->format_from_file_name( $file );
    $tr->format($fmt);
    $tr->input($file);
    my $pops = &_parse_sample_excel( $tr );

    # warn $args->branch(\%pops);
    my $ml  = &get_maploc();
    my $cat = $ml->get_text( $locCatName );
    $cat->tag("Description", 
              "Sites observed in the 1000 Genomes population survey");
    $cat->tag("Short", "1kG");
    $cat->write();

    my $catPar = $ml->get_text( "Polymorphism Categories" );
    $catPar->tag("Member", $cat->text());
    $catPar->write();

    my $rootPop = $ml->get_population($rootName, $stAuth);
    $rootPop->tag_values("Description", "This is the root population for all 1000 genomes populations");
    $rootPop->tag_values("Category", $locCatName);
    $rootPop->update();
    
    # Set up the populations
    my @allPop;
    foreach my $tok (sort keys %{$pops}) {
        my $meta  = $pops->{$tok};
        my @mems  = keys %{$meta->{Member}};
        my @descs = keys %{$meta->{Description}};
        my $famH  = $meta->{Family};
        my @errs;
        push @errs, "Strange token" unless 
            ($tok eq uc($tok) && length($tok) == 3);
        push @errs, "No members" if ($#mems == -1);
        push @errs, "Non-unique description" unless ($#descs == 0);
        unless ($#errs == -1) {
            $args->msg("[?]","Deciding '$tok' is not a population", @errs);
            next;
        }
        my $name = "$prfx.$tok";
        push @allPop, $name;
        my $pop  = $ml->get_population($name, $stAuth);
        my $desc = $descs[0];
        $pop->parent($rootPop);
        $pop->tag_values("Reference Population", $rootName);
        $pop->tag_values("Description", $desc);
        $pop->chr( 2 * ($#mems + 1) );
        $pop->tag_values("Category", $locCatName);
        $pop->update();
        $args->msg("[POPULATION]", $pop->to_text());
        for my $m (0..$#mems) {
            my $mem = $mems[$m];
            my $txt = $ml->get_text( $mem );
            $txt->tag_values("In Population", $name);
            my @bits = ($desc);
            if (my $x = $meta->{Member}{$mem}) { push @bits, $x; }
            $txt->tag_values("Description", join(', ', @bits));
            if (my $fam = $famH->{$mem}) {
                $txt->tag_values("Family", $fam);
            }
            $txt->update();
            $args->msg("[SAMPLE]", $txt->to_text()) unless ($m);
        }
    }
    $args->msg("[+]", "Finished loading ".($#allPop+1)." populations",
               join(", ", @allPop));
}

=head1 Long Run Benchmarks (Chr2)


Elapsed time 18.341 hours, 18.054 are benchmarked
  %Elaps w/Kids     N       Avg  Avg+Kids Method                                  
  ------ ------ ----- --------- --------- ----------------------------------------
   59.2%  59.2%  377k  103.75ms  103.75ms BMS::SnpTracker::MapLoc::Location::add_allele_to_db-Create Allele Row
   25.4%  25.4%  779k   21.52ms   21.52ms BMS::SnpTracker::MapLoc::Location::add_allele_to_db-Set allele frequency
    6.2%   6.2% 1027k    3.98ms    3.98ms BMS::SnpTracker::MapLoc::text_to_pkey   
    3.8%  68.8%  779k    3.24ms   58.28ms BMS::SnpTracker::MapLoc::Location::add_allele_to_db-Set allele row
    3.4%   3.4%   41k   53.82ms   53.82ms BMS::SnpTracker::MapLoc::Tagged::write_tags
    0.3%   0.3% 3173k   57.86us   57.86us BMS::MapTracker::Shared::benchend       
    0.1%   0.1%   41k    1.48ms    1.48ms BMS::SnpTracker::MapLoc::Common::_generic_get_pkey
    0.0%   0.5%   41k  472.95us    7.51ms BMS::SnpTracker::MapLoc::Tagged::all_tag_id_pairs
    0.0%   0.1%   41k  126.74us    1.61ms BMS::SnpTracker::MapLoc::Location::pkey 
    0.0%   0.0%   41k   61.31us   61.31us BMS::SnpTracker::MapLoc::Location::new  
    0.0%   0.0%  3074  525.50us  525.50us BMS::SnpTracker::MapLoc::pkey_to_text   
    0.0%   0.0%     2   87.22ms   88.75ms BMS::SnpTracker::MapLoc::connect        
    0.0%   0.0%   259   92.89us    1.01ms BMS::SnpTracker::MapLoc::cached_pkey_to_text
    0.0%   0.0%    23  790.05us    3.23ms BMS::SnpTracker::MapLoc::Population::new
    0.0%   0.0%    22  636.20us    1.35ms BMS::SnpTracker::MapLoc::Tagged::read_tags
    0.0%   0.0%     2    1.54ms    1.54ms BMS::SnpTracker::MapLoc::schema         
    0.0%   0.0%    22  128.55us    2.45ms BMS::SnpTracker::MapLoc::Population::pkey
    0.0%   0.0%    44   58.87us  764.85us BMS::SnpTracker::MapLoc::Common::_generic_get_id
    0.0%   0.0%    23   80.85us    3.31ms BMS::SnpTracker::MapLoc::get_population_fast
    0.0%   0.0%    21   71.04us  142.60us BMS::SnpTracker::MapLoc::Population::parent
    0.0%   0.0%    22   46.10us    1.40ms BMS::SnpTracker::MapLoc::Population::read
    0.0%   0.0%     2  351.55us  351.55us BMS::SnpTracker::MapLoc::new            
    0.0%   0.0%    21   28.24us   28.24us BMS::SnpTracker::MapLoc::Population::ancestors
    0.0%   0.0%     2  243.43us  243.43us BMS::SnpTracker::MapLoc::get_text       

=cut
