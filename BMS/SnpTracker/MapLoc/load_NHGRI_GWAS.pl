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
use BMS::ArgumentParser;
use BMS::TableReader;
use LWP::UserAgent;

my ($guessBuild, $locCategory, $locCatName, %rejected, %studies, %pops);

my $args = BMS::ArgumentParser->new
    ( -nocgi          => $ENV{'HTTP_HOST'} ? 0 : 1,
      -file           => "gwascatalog.txt",
      -instance       => "maploc2",
      );

$args->shell_coloring( );
my $limit    = $args->val(qw(limit));
my $build    = $args->val(qw(build));
my $doDump   = $args->val(qw(dodump dump));
my $doReport = $args->val(qw(report));
my $dbInst   = $args->val(qw(instance));
my $log10    = 1 / log(10);

my $standardColumns = {
    "PUBMEDID"                  => "PMID",
    "Chr_id"                    => "Chr",
    "Chr_pos"                   => "Pos",
    "Disease/Trait"             => "Trait",
    "Study"                     => "Study",
    "Pvalue_mlog"               => "LOD",
    "Initial Sample Size"       => "Sample",
    "Snp_id_current"            => "SNP",
    "SNPs"                      => "SNPs",
    "Strongest SNP-Risk Allele" => "RiskAllele",
    "Risk Allele Frequency"     => "RiskFreq",
    "p-Value"                   => "pValue",
};
my %keep = map { $_ => 1 } values %{$standardColumns};

my $ml = &init_maploc();
&parse( $args->val(qw(input catalog file)) );
&finish();

sub finish {
    while (my ($pmid, $pdat) = each %studies) {
        my ($study, $samp) = @{$pdat};
        my $obj = $ml->get_text( $pmid );
        $obj->tag("Study", $study);
        $obj->tag("Sample", $samp);
        $obj->write_tags();
    }
    &report_build();
    my @rej = sort keys %rejected;
    unless ($#rej == -1) {
        $args->msg("[!]", scalar(@rej)." studies were rejected:");
        foreach my $study (@rej) {
            my %reas = %{$rejected{$study}};
            $args->msg("[-]", $study, map { "$_ [$reas{$_}]" } sort keys %reas);
        }
    }
}

sub parse {
    my $file = shift;
    return unless ($file);
    &store_raw_data($file) if (!-s $file || $args->val(qw(update)) );
    $args->msg("[+]", "Parsing catalog", $file, "Build : ".($build || "will be inferred"));
    my $tr     = BMS::TableReader->new
        ( -limit     => $limit,
          -colmap    => $standardColumns,
          -format    => 'tsv',
          -hasheader => 1, );
    $tr->input($file);
    unless ($args->val(qw(nomask))) {
        my @ign = $tr->ignore_unmapped_columns();
        $args->msg("[v]", "Ignoring some input columns", @ign)
            unless ($#ign == -1);
    }
    my $mapped = 0;
    while (my $data = $tr->next_clean_hash()) {
        if (!$data->{SNP}) {
            if (my $s = $data->{SNPs}) {
                my @ss;
                foreach my $rs (split(/\s*\,\s*/, $s)) {
                    if ($rs =~ /^rs(\d+)/) {
                        push @ss, $1;
                    }
                }
                my ($low) = sort { $a <=> $b } @ss;
                $data->{SNP} = $low if ($low);
            }
            if (!$data->{SNP} && $data->{RiskAllele} &&
                $data->{RiskAllele} =~ /^rs(\d+)\-/) {
                # Extract rs number from risk allele
                $data->{SNP} = $1;
            }
        }
        warn &_hash_to_text( $data ) if ($doDump);
        unless ($build) {
            &infer_build($data);
        }
        while (my ($key, $val) = each %{$data}) {
            if ($val =~ /^(NR|NS)$/) {
                $data->{$key} = "";
            }
        }
        my ($chr, $pos, $study, $ra, $rf) = map
        { $data->{$_} } qw(Chr Pos Study RiskAllele RiskFreq);
        unless ($study) {
            $rejected{"Unknown Study"}{"Study name was not defined"}++;
            next;
        }
        my $popName = $study;
        my ($lod, $trait, $pmid) = map { $data->{$_} } qw(LOD Trait PMID);
        unless ($trait) {
            $rejected{$study}{"Trait name was not defined"}++;
            next;
        }
        my @dbits = ($trait);
        my $pop;
        if ($pmid) {
            if ($pmid =~ /^\d+$/) {
                $pmid = "PMID:$pmid";
                push @dbits, $pmid;
                $pop = &get_pop( $pmid, $study, $lod);
            } else {
                $args->msg_once("Malformed PubMed ID '$pmid'");
            }
        } else {
            $rejected{$study}{"No PubMed ID"}++;
            next;            
        }
        if (defined $lod && $lod ne '') {
            $lod = int(0.5 + $lod * 100) / 100;
            push @dbits, "LOD:$lod";
        } elsif (my $pv = $data->{pValue}) {
            if ($pv + 0) {
                $lod = 0 - (log($pv) * $log10);
                $lod = int(0.5 + $lod * 100) / 100;
                push @dbits, "LOD:$lod";
            } else {
                push @dbits, "LOD:9999";
            }
            # warn "$pv = $lod";
        } else {
            $lod = undef;
        }

        my @locs;
        if (!$build) {
            my $what = ($chr && $pos) ? "coordinates" :
                $data->{SNPs} ? "dbSNP IDs" : "stuff";
            $rejected{"General"}{"Can not map $what without build"}++;
            next;
        } elsif ($chr && $pos) {
            push @locs, $ml->request_to_location( $chr, $pos );
        } elsif (my $ss = $data->{SNPs}) {
            foreach my $rs (split(/\s*\,\s*/, $ss)) {
                if ($rs =~ /^rs(\d+)$/) {
                    foreach my $loc ($ml->accession_to_location( $rs )) {
                        if ($build && $loc->build eq $build) {
                            push @locs, $loc;
                        }
                    }
                }
            }
        }
        if ($#locs == -1) {
            $rejected{$study}{"Insufficient genomic location data"}++;
            next;
        }
        $mapped++;
        if ($ra) {
            # It would be nice to store these alleles in the DB
            # However, I am not certain I can reliably get the strand right.
            if ($ra =~ /rs\d+\-(.+)/) {
                my $base = $1;
                $ra = "RiskAllele:$base";
                if ($rf && $rf =~ /^0?\.\d+$/) {
                    $ra .= sprintf("-%d%%", 100 * $rf);
                }
                if ($pop) {
                    foreach my $loc (@locs) {
                        $loc->add_allele_to_db($base, $pop->pkey());
                    }
                }
                push @dbits, $ra;
            } else {
                $args->msg_once("Malformed Risk Allele '$ra'");
            }
        }
        my $desc = join(' ', @dbits);
        $studies{$pmid} ||= [ $data->{Study}, $data->{Sample} ];
        foreach my $loc (@locs) {
            $loc->tag("NHGRI GWAS", $desc);
            $loc->tag("Category", $locCatName);
            $loc->tag("Tally", "NHGRI GWAS");
            $loc->add_category( $locCategory );
            $loc->write();
            if ($doReport) {
                $loc->read();
                $args->msg("[+]", $loc->to_text());
            }
        }
    }
    $args->msg("[DONE]","Mapped $mapped");
}

sub get_pop {
    my ($pmid, $study, $lod) = @_;
    my $name = "NHGRI GWAS $pmid";
    unless ($pops{$name}) {
        my $pop = $pops{$name} = $ml->get_population($name, "NHGRI GWAS");
        $pop->tag_values('PubMed', $pmid);
        $pop->tag_values('Description', $study);
        $pop->tag_values('Category', $locCatName);
        $pop->tag_values('Tally', "NHGRI GWAS");
        $pop->update();
    }
    return $pops{$name};
}

sub _hash_to_text {
    my $data = shift;
    my $txt  = "";
    foreach my $key (sort keys %{$data}) {
        my $val = $data->{$key};
        if (!defined $val) {
            $val = '-UNDEF-';
        } elsif ($val eq '') {
            $val = '-EMPTY STRING-';
        } elsif ($val =~ /^\s+/ || $val =~ /\s+$/) {
            $val = "'$val'";
        }
        $txt .= sprintf(" %20s : %s\n", $key, $val);
    }
    return $txt;
}

sub report_build {
    return unless ($guessBuild);
    my @snps = keys %{$guessBuild};
    my %builds;
    map { $builds{$_}++ } map { @{$guessBuild->{$_}} } @snps;
    $args->msg("Catalog scanned in order to infer genome build being used",
               scalar(@snps)." dbSNP accessions compared to their known locations",
               "Consider using the following genome builds:",
               map { "-build $_ [$builds{$_} SNPs]" } sort { $builds{$b} <=> $builds{$a} } keys %builds);
}

sub infer_build {
    my $data = shift;
    my $rsNum = $data->{SNP};
    unless ($rsNum) {
        return;
    }
    my $rsid      = "rs$rsNum";
    $guessBuild ||= {};
    return if ($guessBuild->{$rsid});
    my ($chr, $pos) = ($data->{Chr}, $data->{Pos});
    return unless ($chr && $pos);

    my $lRows     = $ml->accession_to_location( $rsid );
    my @found;
    foreach my $lid (map { $_->[0] } @{$lRows}) {
        my $loc = $ml->get_location($lid);
        next unless ($loc->chr eq $chr);
        next unless ($loc->pos eq $pos);
        push @found, $loc->build();
    }
    if ($#found == -1) {
        push @found, "UNKNOWN";
    } else {
        my %u = map { $_ => undef } @found;
        @found = sort keys %u;
    }
    $guessBuild->{$rsid} = \@found;
}

sub store_raw_data {
    my $file = shift || "gwascatalog.txt";
    my $url = "https://www.genome.gov/admin/gwascatalog.txt";
    my $ua   = LWP::UserAgent->new;
    my $response = $ua->get($url);
    
    if ($response->is_success) {
        open(RESP, ">$file") || $args->death
            ("Failed to store local catalog", $file, $!);
        
        print RESP $response->decoded_content;
        close RESP;
    } else {
        $args->death("Failed to recover GWAS Catalog",
                     $url, $response->status_line );
    }
    return $file;
}

sub init_maploc {
    my $mapLoc = BMS::SnpTracker::MapLoc->new
        ( -noenv    => $args->val(qw(noenvironment noenv)),
          -instance => $dbInst,
          # -makedb => $args->val(qw(makedatabase makedb)),
          );
    $mapLoc->build($build);

    $locCatName = "NHGRI GWAS Polymorphisms";
    my $cat = $mapLoc->get_text( $locCatName );
    $cat->tag("Description", 
              "Collection of polymorphisms that have been reported as significantly linked to a trait in one or more genome-wide association studies");
    $locCategory = $cat->pkey();
    $cat->tag("Useful Tag", "NHGRI GWAS");
    $cat->write();

    my $catPar = $mapLoc->get_text( "Polymorphism Categories" );
    $catPar->tag("Description", "Variant locations attributed to germ line polymorphism across populations");
    $catPar->tag("Member", $cat->text());
    $catPar->write();
    return $mapLoc;
}
