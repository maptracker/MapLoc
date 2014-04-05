#!/usr/bin/perl -w

my $isBeta;

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

my $VERSION = 
    ' $Id: mapLocReporter.pl,v 1.2 2012/12/12 15:39:20 tilfordc Exp $ ';

use strict;
use BMS::ArgumentParser;
use BMS::SnpTracker::MapLoc;
use BMS::TableReader;

my $args = BMS::ArgumentParser->new
    ( -nocgi      => $ENV{HTTP_HOST} ? 0 : 1,
      -build      => 'GRCh37',
      -instance   => 'maploc2',
      -minfreq    => 0.05,
      -range      => 10000,
      -progress   => 60,
      -taxa       => 9606,
      -srcdir     => "",
      -casesensitive => 0,
      # -paramfile  => [ $pfile,  'BMS::SnpTracker::MapLoc'],
      -tiddlywiki => 'MapLocation' );

$args->shell_coloring();
$args->debug->skip_key( [qw(MAPLOC sequences PARENT)], 'global' );

my $nocgi     = $args->val(qw(nocgi));
my $vb        = $args->val(qw(vb verbose)); 
$vb           = 1 if (!defined $vb && $nocgi);
my $limit     = $args->val(qw(limit)) || 0;
my $explSQL   = $args->val(qw(explain explainsql)) || 0;
my $dumpSQL   = $args->val(qw(dumpsql)) || $explSQL || 0;
my $isPretty  = $args->val(qw(pretty ispretty));
my $mlInt     = $args->val(qw(instance));
my $taxReq    = $args->val(qw(species taxid taxa));
my $isTrial   = $args->val(qw(trial istrial testmode tm)) || 0;
my $beQuiet   = $isTrial =~ /quiet/i ? 1 : 0;
my $progress  = $args->val(qw(progress));
my $srcDir    = $args->val(qw(srcdir));

my ($llCatName, $rsCatName, $ml);
my $statTag = "RefSeq Status";
&set_basic();

if ($srcDir) {
    $srcDir =~ s/\/$//;
    $args->msg("[^]", "Finding data in $srcDir/");
    &load( );
} else {
    $args->msg("[USAGE]","This script loads basic Entrez and RefSeq data into the MapLoc database",
               "You need to indicate the source directory of the Entrez flat files",
               "   -srcdir /path/to/directory   (use '.' for local directory)");
}

sub load {
    $args->msg("[!!]","Trial mode - database will not be altered")
        if ($isTrial);
    &load_info();
    &load_refseq();
    $args->msg("Finished");
}

sub load_refseq {
    my $file = shift || "$srcDir/gene2refseq.gz";
    if ($args->val(qw(skiprefseq))) {
        $args->msg("[x]","Request to skip file", $file);
        return;
    }
    my $rsLim = $limit;
    $limit = 0;
    # We are just going to get basic Locus <-> RNA <-> Protein associations

    unless (-s $file) {
        $args->msg("[x]","File not found, skipping", $file);
        return;
    }
    my $tr = &entrez_reader( $file, 'refseq' );
    my $collect = { LL => "" };
    my $done    = 0;
    my $ti      = time;
    while (my $data = &_get_data($tr)) {
        my $ll = $data->{LL};
        if ($collect->{LL} ne $ll) {
            $done += &_do_refseq( $collect );
            if ($progress && 
                time - $ti > $progress) {
                my $num = $tr->rowcount();
                $args->msg("[$num]", sprintf("%s", $ll));
                $ti = time;
            }
            $collect = {
                LL => $ll,
            };
            last if ($rsLim && $done >= $rsLim);
        }
        my $prot = $data->{Protein} || "";
        my $rna  = $data->{RNA} || "";
        my $targ = $collect->{RNA}{$rna} ||= { };
        $targ->{Protein}{$prot} = 1 if ($prot);
        if (my $stat = $data->{Status}) {
            $collect->{Status}{$stat} = 1;
            $targ->{Status}{$stat} = 1;
        }
    }
    &_do_refseq( $collect );
    $limit = $rsLim;
}

sub _do_refseq {
    my $collect = shift;
    my $lid     = $collect->{LL};
    return 0 unless ($lid);
    
    my $gene = $ml->get_gene( $lid );
    $gene->tag("Category", $llCatName);
    unless ($isTrial) {
        map { $gene->clear_tag($_) } ($statTag, "Has RNA", "Has Protein");
    }
    if (my $stats = $collect->{Status}) {
        # The status is primarily at the RNA level.
        # A locus can have multiple statuses because of this eg LOC354
        # Entrez appears to report just one, though.
        my @u = keys %{$stats};
        if ($#u == -1) {
            $args->msg_once("[?]", "No status for $lid");
        } else {
            map { $gene->tag($statTag, $_ ) } @u;
        }
        while (my ($rnaV, $pHash) = each %{$collect->{RNA} || {}}) {
            my $rnaU = $rnaV;
            $rnaU =~ s/\.\d+$//;
            delete $pHash->{""};
            if ($rnaU) {
                $gene->tag("Has RNA", $rnaU);
            }
            my $stat = [ keys %{$pHash->{Status}} ];
            if ($#{$stat} == 0) {
                $stat = $stat->[0];
            } else {
                $args->msg_once("[?]", "Non-unique status for $rnaV :".
                                (join(',', @{$stat}) || '-NULL-'));
                $stat = "";
            }
            my @prts = keys %{$pHash->{Protein}};
            if ($#prts == 0) {
                my $prtV = $prts[0];
                my $prtU = $prtV;
                $prtU =~ s/\.\d+$//;
                $gene->tag("Has Protein", $prtU);
                if ($rnaU) {
                    # You would think that every protein would have an RNA
                    # but for some reason some species do not
                    # (single-celled critters)
                    # If we do have both RNA and protein, associate them:
                    my $rvo = $ml->get_text($rnaV);
                    $rvo->tag("Has Protein", $prtV);
                    $rvo->tag("Accession", $prtV);
                    my $ruo = $ml->get_text($rnaU);
                    $ruo->tag("Has Protein", $prtU);
                    $ruo->tag("Accession", $prtU);
                    map { $_->tag("Category", $rsCatName) } ($ruo, $rvo);
                    unless ($isTrial) {
                        foreach my $ct ($statTag, "Has Protein") {
                            map { $_->clear_tag( $ct ) } ($ruo, $rvo);
                        }
                    }
                    if ($stat) {
                        $ruo->tag($statTag, $stat );
                        $rvo->tag($statTag, $stat );
                    }
                    if ($isTrial) {
                        warn $ruo->to_text() . $rvo->to_text() unless ($beQuiet);
                    } else {
                        $ruo->update();
                        $rvo->update();
                    }
                }
            } elsif ($#prts != -1) {
                $args->msg_once("[?]", "Multiple proteins for $rnaV: ".
                                join(',', @prts));
            }
        }
    } else {
        return 0;
    }
    if ($isTrial) {
        # warn $args->branch($collect); warn "\n";
        print $gene->to_text() unless ($beQuiet);
    } else {
        $gene->update();
    }
    return 1;
}

sub load_info {
    my $file = shift || "$srcDir/gene_info.gz";
    if ($args->val(qw(skipinfo))) {
        $args->msg("[x]","Request to skip file", $file);
        return;
    }
    unless (-s $file) {
        $args->msg("[x]","File not found, skipping", $file);
        return;
    }
    my $tr = &entrez_reader( $file, 'info' );

    my @direct = ('Description','Gene Type','GeneID','Map Location','Symbol',
                  'Official Symbol', 'Other Symbols');
    my $ti = time;
    while (my $data = &_get_data($tr)) {
        my $lid = $data->{LL};
        next unless ($lid);
        my $gene = $ml->get_gene( $lid );
        $gene->tag("Object Type", "Gene");
        $gene->tag("Accession", $lid);

        foreach my $key (@direct) {
            $gene->clear_tag($key) unless ($isTrial);
            foreach my $val (split(/\s*\|\s*/, $data->{$key} || "")) {
                $gene->tag($key, $val) unless (!defined $val || $val eq '');
            }
        }
        foreach my $val (split(/\s*\|\s*/, $data->{"DB XREF"} || "")) {
            if ($val =~ /^([^\:]+)\:(.+)$/) {
                $gene->tag("$1 XREF", $2);
            }
        }
        if ($isTrial) {
            print "\n".$gene->to_text() unless ($beQuiet);
        } else {
            $gene->update();
        }
        if ($progress && 
            time - $ti > $progress) {
            my $num = $tr->rowcount();
            $args->msg("[$num]", sprintf("%s %s %s", $lid, $data->{Symbol} || '-',
                                         $data->{Description} || ""));
            $ti = time;
        }
    }
}

sub _get_data {
    my $tr = shift;
    my $data = $tr->next_clean_hash();
    if ($data) {
        my $gid = $data->{GeneID} || "";
        if ($gid =~ /^\d+$/) {
            $data->{LL} = "LOC$gid";
        } else {
            $args->msg("Malformed GeneID '$gid'");
        }
    }
    return $data;
}

sub entrez_reader {
    my ($file, $what) = @_;
    $args->death("Could not find Entrez $what file",
                 "-file '$file'") unless (-s $file);
    
    my $tr     = BMS::TableReader->new
        ( -limit     => $limit,
          -file      => $file,
          -format    => 'tsv',
          -hasheader => 0, );

    $args->msg("Parsing Entrez $what file", $file);
    # In case the column headers ever change....
    my $colAli = {
        tax_id       => 'TaxID',
        description  => 'Description',
        type_of_gene => 'Gene Type',
        map_location => 'Map Location',
        Synonyms     => 'Other Symbols',
        dbXrefs      => 'DB XREF',
        status       => 'Status',
        orientation  => 'Strand',
        'protein_accession.version' => 'Protein',
        'RNA_nucleotide_accession.version' => 'RNA',
        'genomic_nucleotide_accession.version' => 'gDNA',
        'RNA_nucleotide_gi' => 'GI RNA',
        'genomic_nucleotide_gi' => 'GI gDNA',
        'protein_gi' => 'GI Protein',
        'assembly' => 'Build',
        'end_position_on_the_genomic_accession' => 'gDNA Stop',
        'start_position_on_the_genomic_accession' => 'gDNA Start',
        Symbol_from_nomenclature_authority => 'Official Symbol',
        
    };

    $tr->remap_value('-', "");
    $tr->remap_value(undef, "");

    # Figure out the header
    my $head;
    while (my $row = $tr->next_clean_row()) {
        $tr->extend_limit( );
        if ($row->[0] =~ /^\#format\s*\:?\s+(.+)/i) {
            $head = [];
            foreach my $col (split(/\s+/, $1)) {
                push @{$head}, $colAli->{$col} || $col;
            }
            last;
        }
    }
    unless ($head) {
        $args->death("Failed to find header row in file");
    }
    $tr->set_header($head);
    my @filters;
    if ($taxReq) {
        $args->msg("[+]", "Only capturing tax_id = '$taxReq'");
        push @filters, sub {
            my $row = shift;
            return $row->[0] eq $taxReq ? 0 : 1;
        };
    }
    if ($file =~ /gene2refseq/) {
        # Exclude all but the reference sequence
        my $ind = $tr->column_name_to_index('Build');
        my $match = "Primary Assembly";
        $args->msg("[+]", "Only capturing rows where the build matches '$match'");
        push @filters, sub {
            my $row = shift;
            #die $args->branch($row);
            #warn $row->[$ind];
            return $row->[$ind] =~ /\Q$match\E/ ? 0 : 1;
        };

    }
    unless ($#filters == -1) {
        $tr->toss_filter( sub {
            my $row = shift;
            for my $f (0..$#filters) {
                return 1 if (&{$filters[$f]}( $row ));
            }
            return 0;
        });
    }

    return $tr;
}


sub wget {
    my $ftp = "ftp://ftp.ncbi.nih.gov/gene/DATA/";
    $args->msg("Recovering Entrez data from the NCBI", $ftp);
    $args->death("WORKING HERE");
}

sub set_basic {
    $ml = BMS::SnpTracker::MapLoc->new
        ( -build    => $args->val(qw(build)),
          -noenv    => $args->val(qw(noenvironment noenv)),
          -instance => $mlInt, );

    $llCatName = "NCBI Entrez";
    $rsCatName = "NCBI RefSeq";

    my $llCat = $ml->get_text( $llCatName );
    $llCat->tag("Description", 
              "Data provided by the National Center for Biotechnology Information Entrez database");
    $llCat->tag("URL Tag", 'tag="'.$statTag.'" url="http://www.ncbi.nlm.nih.gov/projects/RefSeq/key.html#status"');
    $llCat->update();

    my $rsCat = $ml->get_text( $rsCatName );
    $rsCat->tag("Description", 
                "Data provided by the National Center for Biotechnology Information RefSeq database");
    $rsCat->tag("URL Tag", 'tag="'.$statTag.'" url="http://www.ncbi.nlm.nih.gov/projects/RefSeq/key.html#status"');
    $rsCat->update();

}
