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

use strict;
use BMS::ArgumentParser;
use BMS::SnpTracker::MapLoc;
use BMS::ForkCritter;
use BMS::TableReader;

my $args = BMS::ArgumentParser->new
    ( -limit    => 0,
      -nocgi    => 1,
      -fork     => 20,
      -verbose  => 1,
      -build    => 'GRCh37',
      -ageall   => '3 May 2012',
      -instance   => 'maploc2',
      -progress => 120,
      -aspercent => 1,
      -binsize   => 10,
      );

$args->shell_coloring( );


# sed filter useful for browsing files:
# sed 's/[^ ]*\/[^ ]*//g' someFile.txt | sed 's/dbSNP_[0-9]*//' | sed 's/unknown//g' | sed 's/NA//g' | sed 's/PASS//g' | less -S

my ($fc, %stuff, $localML, %populations, $locCategory, $locCatName,
    $currentRow, $capCols);

my $debug    = $args->val(qw(debug));
my $vb       = $args->val(qw(vb verbose)) || 0;
my $tm       = $args->val(qw(tm testmode trial)) || 0;
my $prog     = $args->val(qw(prog progress)) || 300;
my $limit    = $args->val(qw(limit));
my $forkNum  = $args->val(qw(forknum fork)) || 1;
my $build    = $args->val(qw(build));
my $dumpCols = $args->val(qw(dumpcol dumpcols));
$dumpCols    = [split(/\s*[\t\r\n\,]+\s*/, $dumpCols)] 
    if ($dumpCols && !ref($dumpCols));
my $doDump   = $args->val(qw(dump)) || ($dumpCols ? 1 : 0);
my $loadPops = $args->val(qw(loadpop loadpops dopop));
my $deferCache = $args->val(qw(defer defercache defersize)) || 300;
if ($clearPop) {
    $clearPop = {};
    $forkNum = 1;
}

my $auth = "NHLBI";
my $fullAcronym = "NHLBI GO ESP";
my $fullDesc    = "National Heart, Lung, and Blood Institute Grand Opportunity Exome Sequencing Project";

my %builds = ( 36 => 'NCBI36',
               37 => 'GRCh37' );

my %nullVals = map { lc($_) => 1 } qw(none unknown na);

my $standardColumns = {
    "AA-EstimatedAge(kyrs)" => "Estimated African American Age (kyrs)",
    "EA-EstimatedAge(kyrs)" => "Estimated European American Age (kyrs)",
    "AvgSampleReadDepth"    => "Average $auth Read Depth",
    # "ChimpAllele"           => "Chimp Allele",
    "ClinicalInfo"          => "Clinical Info",
    "ConservationScoreGERP" => "GERP Conservation",
    "ConservationScorePhastCons" => "PhastCons Conservation",
    "GranthamScore" => "Grantham Score",
    "Polyphen2(Class:Score)" => "Polyphen2",
    "" => "",
    "" => "",
    "" => "",
};

if (0) {
    my $ml = BMS::SnpTracker::MapLoc->new
        ( -build  => $build, );
    my $dbh = $ml->dbh;
    $args->msg("[-]", "Creating / updating database tables");
    $dbh->make_all();
    $args->msg("[+]", "Database tables created");
    exit;
}

my @ignoring = qw();
my $ignoreCol = { map { $_ => 1 } @ignoring };


if (my $mafDir = $args->val(qw(mafdir dir))) {
    &load_dir($mafDir);
} elsif (my $file = $args->val(qw(file maf))) {
    &load_file($file);
} else {
    &usage();
}

sub usage {
    my $txt = <<USAGETXT;

This script loads a MapLoc variant database with variants from the
National Heart, Lung, and Blood Institute. To use the program, you
will need one or more MAF files. I found the files as a gzipped
tarball in the Downloads tab here:

  http://evs.gs.washington.edu/EVS/

... by recovering the '.txt.tar.gz' file. You will need to gunzip and
untar the file. I do not recall seeing VCF files when I first loaded
the data, but I either overlooked them or they were not initially
available. They are stated to contain the same information as the MAF
file.

You can load a single file with:

  load_NHLBI_maf.pl -file blahblah.chr2.snps_indels.txt

If you want to test the program, you can also use debug mode and
specify a handful of locations to process:

  load_NHLBI_maf.pl -file blahblah.chr2.snps_indels.txt \
                    -limit 10 -debug

To load an entire directory, use the -dir parameter, pointing it to
the folder you untarred:

  load_NHLBI_maf.pl -dir /tmp/downloadedNHLBI/

Currently NHLBI is only reporting data for American individuals of
European or African descent. Further populations should be
forthcoming, in which case the $standardColumns hash should be updated
to indicate which additional metadata columns should be captured (ie
allele age). The code should detect the alleles automatically, but at
the moment populations are being reduced to two letters, so
AfricanAmerican becomes AA. This is likely a bad idea in the long
term, since potential future entries such as AsianAmerican would then
colide.

USAGETXT
    
    $args->msg("[USAGE]", split(/\n\r?/, $txt));
    exit;
}


sub load_dir {
    my $dir = shift;
    foreach my $file ($args->read_dir($dir)) {
        # $args->msg("[DEBUG]", $file); next;
        &load_file($file);
    }
}

sub load_file {
    my $file = shift;
    $args->death("Failed to find requested file", $file) unless (-s $file);
    $args->msg("[+]","Reading MAF file", $file,`date`);

    $fc ||= BMS::ForkCritter->new
        ( -init_meth   => \&init,
          -finish_meth => \&finish_fc,
          -method      => \&process_rows,
          -limit       => $limit,
          -progress    => $prog,
          -verbose     => $vb );
    $fc->reset();
    $fc->group_method( \&group_entries );
    $fc->set_column_separator('\s+');
    # $fc->set_colmap( $standardColumns );
    $fc->input_type('maf header clean');
    $fc->input($file);
    if (my $failed = $fc->execute( $forkNum )) {
        $args->err("$failed children failed to finish execution");
    }

    
}

sub group_entries {
    my $row = shift;
    if ($#rowGroup == -1 || $rowGroup[0][0] eq $row->[0]) {
        push @rowGroup, $row;
        return undef;
    }
    # New ID
    my $rv = {
        num => $#rowGroup + 1,
        id  => $rowGroup[0][0],
    };
    map { shift @{$_} } @rowGroup;
    $rv->{data} = [ sort { $b->[0] <=> $a->[0] } @rowGroup ];
    @rowGroup = ();
    return $rv;
}


sub process_rows {
    my $rows = shift;
    my $struct = {};
    foreach my $row (@{$rows}) {
        &process_row($row, $struct);
    }
    foreach my $dat (values %{$struct}) {
        my $loc = $dat->{loc};
        while (my ($base, $pH) = each %{$dat->{alleles}}) {
            while (my ($pid, $fcnt) = each %{$pH}) {
                my ($num, $total) = (@{$fcnt});
                $loc->add_allele_to_db
                    ($base, $pid, $num / $total, $total, 1 );
            }
        }
        if ($debug) {
            print $loc->to_text();
        } else {
            $loc->update();
        }
    }
    $localML->add_deferred_alleles() 
        if (!$debug && $localML->deferred_allele_count() > $deferCache);
}
sub process_row {
    my ($row, $struct) = @_;
    my @data = map { !defined $_ || $nullVals{lc($_)} ? "" : $_ } @{$row};
    my $gtxt = $data[$capCols->{Chr}];
    return unless ($gtxt);
    my ($chr, $pos);
    if ($gtxt =~ /^([XYZW]|\d{1,2})\:(.+)$/) {
        ($chr, $pos) = ($1, $2);
    } else {
        &err("Unrecognized location: $gtxt");
        return;
    }
    my @alleles;
    my $allAlleleTxt = $data[$capCols->{Alleles}];
    foreach my $atxt (split(/\;/, $allAlleleTxt)) {
        if ($atxt =~ /^([^>]+)>([^>]+)$/) {
            my ($r, $o) = ($1, $2);
            if ($#alleles == -1) {
                # Reference
                $alleles[0] = $r;
            } elsif ($alleles[0] ne $r) {
                &err("Multiple reference alleles defined at $gtxt: $alleles[0]/$r");
                return;
            }
            push @alleles, $o;
        }
    }
    if ($#alleles == -1) {
        &err("No alleles defined at $gtxt");
        return;
    }
    my $l = $pos - 1;
    my $r = $pos + length($alleles[0]);
    my $tidy = 1;
    while ($tidy) {
        $tidy = 0;
        # Can we shorten up the alleles?
        my %firstChar;
        foreach my $base (@alleles) {
            if ($base) {
                $firstChar{substr($base, 0, 1)} = 1;
            } else {
                # Gap allele
                $firstChar{""} = 1;
            }
        }
        my @u = keys %firstChar;
        if ($#u == 0) {
            # All the alleles have the same first base
            # This is a mechanism for dealing with gap alleles
            # MapLoc instead uses Left/Right flanks
            # Shorten each allele:
            map { $alleles[$_] = substr($alleles[$_], 1) } (0..$#alleles);
            $l++; # Increment the left flank
            $tidy = 1; # Recurse
            # warn "$allAlleleTxt [$pos] = ".join('/', @alleles)." [$l,$r]\n";
        }
    }
    my %alleleLU;
    for my $i (0..$#alleles) {
        # Need to map "R" to the reference allele and A1, A2, A3 to others
        # Empty strings are gap alleles
        $alleleLU{$i ? "A$i" : "R"} = $alleles[$i] || '-';
    }

    my $loc = $localML->get_location($chr, $l, $r);
    my $targ = $struct{$loc->pkey()} ||= {
        loc     => $loc,
        alleles => {},
    };
    while (my ($pid, $ind) = each %{$capCols->{DataID}}) {
        my %aHash;
        my $total = 0;
        if (my $cntTxt = $data[$ind]) {
            foreach my $allTxt (split(/\//, $cntTxt)) {
                if ($allTxt =~ /^(.+)=(\d+)$/) {
                    my ($base, $num) = ($alleleLU{$1} || $1, $2);
                    $total += $num;
                    $aHash{$base} = $num;
                } else {
                    &err("Could not interpret allele count '$allTxt'",
                         $loc->to_one_line());
                }
            }
        }
        # warn $args->branch(\%aHash);
        while (my ($base, $num) = each %aHash) {
            if (my $old = $targ->{$base}{$pid}) {
            } else {
                $targ->{$base}{$pid} = [$num, $total];
            }
        }
    }
    while (my ($tag, $ind) = each %{$capCols->{Tags}}) {
        my $val = $data[$ind];
        next if ($val eq "");
        $loc->tag($tag, $val);
    }
    $loc->add_category( $locCategory );
    #print $loc->to_text();
}

sub err {
    $args->msg("[!]", @_);
}

sub init {
    #die $args->branch($fc);
    my @head = @{$fc->{HEAD_ARRAY} || []};
    $capCols = {};
    my $isFirst = $fc->child() == 1 ? 1 : 0;
    for my $i (0..$#head) {
        my $col = $head[$i];
        if ($col =~ /^\#?base\((.+)\)/i) {
            # Header reports build
            # Column contains genome coordinates
            my $gb = $capCols->{Build} = $1;
            if ($gb =~ /(\d+)/) {
                $gb = $1;
            }
            if (my $stnd = $builds{$gb}) {
                $build = $stnd;
            } else {
                $args->msg("Cound not recognize genome build '$gb'")
                    if ($isFirst);
                exit;
            }
            $capCols->{Chr} = $i;
        } elsif ($col eq 'Alleles'){
            $capCols->{Alleles} = $i;
        } elsif ($col =~ /(.+)AlleleCount/){
            my $cc = $1;
            # http://www.perlmonks.org/?node_id=639595
            my $nice = join(' ', $cc =~ /([A-Z][a-z]*)/g);
            if (my $oc = $capCols->{Data}{$nice}) {
                $args->msg("Multiple Allele Counts",
                           "$nice = $oc + $i");
            } else {
                $capCols->{Data}{$nice} = $i;
            }
        } elsif (my $nice = $standardColumns->{$col}){
            $capCols->{Tags}{$nice} = $i;
        } else {
            push @{$capCols->{Ignored}}, "$i : $col";
        }
    }
    unless (defined $capCols->{Alleles}) {
        $args->msg("Failed to find allele column")
            if ($isFirst);
        exit;
    }
    my @pops = sort keys %{$capCols->{Data}};
    if ($isFirst) {
        my $file = $fc->input();
        $args->msg("[<]", "$file","Build: $build", "-- Populations --",
                   @pops,"--     Tags    --",
                   (sort keys %{$capCols->{Tags}}));
    }
    $localML = &maploc();
    foreach my $name (@pops) {
        my $pname = "$auth $name";
        my $pop = $localML->get_population( $pname, $fullAcronym );
        my $desc = "$name individuals in the $fullDesc";
        unless ($name eq 'All') {
            my $par = $localML->get_population( "$auth All", "NHLBI GO ESP" );
            $pop->parent($par);
        }
        $pop->tag_values("Description", $desc);
        $pop->tag_values("Category", $locCatName);
        # $pop->to_text();
        $pop->update() if ($isFirst);
        my $pid = $pop->pkey();
        $capCols->{DataID}{$pid} = $capCols->{Data}{$name};
    }
    # die $args->branch( $capCols);
}

sub finish_fc {
    if ($fc && $localML) {
        my $cn = $fc->child;
        $localML->add_deferred_alleles() unless ($debug);
        if ($cn == 1) {
            $args->msg("[DEBUG]", "Forked child benchmarks:",
                       $localML->_allele_update_status());
            print STDERR &benchmark_table( $localML );
        }
    }
}

sub benchmark_table {
    my $obj = shift;
    return "" unless ($obj);
    return $obj->showbench( -shell => 1, -minfrac => 0.0001 );
}

sub maploc {
    my $ml = BMS::SnpTracker::MapLoc->new
        ( -build  => $build, );
    $locCatName = $fullAcronym;
    my $cat = $ml->get_text( $locCatName );
    $cat->tag("Description", $fullDesc );

    $cat->tag("Useful Tag", values %{$standardColumns});

    $locCategory = $cat->pkey();
    $cat->write();

    my $catPar = $ml->get_text( "Polymorphism Categories" );
    $catPar->tag("Member", $cat->text());
    $catPar->write();
    return $ml;
}


=head1 Sample Data

      {#base(NCBI.37)} => 14:19377672
      {AA-EstimatedAge(kyrs)} => 48.2+/-83.0
      {AfricanAmericanAlleleCount} => C=51/T=4345
      {AfricanAmericanGenotypeCount} => CC=3/CT=45/TT=2150
      {AllAlleleCount} => C=51/T=12937
      {AllGenotypeCount} => CC=3/CT=45/TT=6446
      {Alleles} => T>C
      {AvgSampleReadDepth} => 34
      {ChimpAllele} => T
      {ClinicalInfo} => unknown
      {ConservationScoreGERP} => NA
      {ConservationScorePhastCons} => 0.2
      {EA-EstimatedAge(kyrs)} => NA
      {EuropeanAmericanAlleleCount} => C=0/T=8592
      {EuropeanAmericanGenotypeCount} => CC=0/CT=0/TT=4296
      {FilterStatus} => PASS
      {FunctionGVS} => missense
      {GeneAccession} => NM_001013354.1
      {Genes} => OR11H12
      {GranthamScore} => 22
      {GwasPubMedInfo} => unknown
      {MAFinPercent(EA/AA/All)} => 0.0/1.1601/0.3927
      {OnIlluminaHumanExomeChip} => no
      {Polyphen2(Class:Score)} => benign:0.359
      {RefBaseNCBI37} => T
      {codingDnaSize} => 981
      {dbSNPVersion} => dbSNP_131
      {hgvsCdnaVariant} => c.79T>C
      {hgvsProteinVariant} => p.(F27L)
      {rsID} => rs78043369

=cut
