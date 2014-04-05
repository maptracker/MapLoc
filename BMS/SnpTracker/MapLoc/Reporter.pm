# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::SnpTracker::MapLoc::Reporter;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

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
use BMS::Utilities::Benchmark;
use Scalar::Util qw(weaken);
use BMS::ExcelHelper;
use BMS::SnpTracker::MapLoc;

use vars qw(@ISA);
@ISA      = qw(BMS::Utilities::Benchmark BMS::Utilities::Escape);

our $snpPopName = '1kG.All';

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self = {
        IMP_HIGH => [qw(CPX FRM STP)],
        IMP_MED  => [qw(DEL NON)],
        IMP_LOW  => [qw(SP5 SP3 SPL)],
        FLANK    => 10,
        BUILD    => 'GRCh37',
        BULKFREQ => {},
        POP_COLS => {},
    };
    bless ($self, $class);
    $self->bench_start();
    $self->{MAPLOC} = shift;
    $self->bench_end();
    $self->add_population_column( $snpPopName );

    return $self;
}

sub default_build {
    my $self = shift;
    return $self->{BUILD};
}

sub maploc_instance {
    return 'maploc2';
}

sub DESTROY {
    my $self = shift;
    return unless ($self);
    $self->disconnect();
}

sub maploc {
    my $self = shift;
    unless ($self->{MAPLOC}) {
        $self->{MAPLOC} = BMS::SnpTracker::MapLoc->new
            ( -build    => $self->default_build(),
              -instance => $self->maploc_instance(),);
    }
    return $self->{MAPLOC};
}

sub bulk_frequencies {
    my $self = shift;
    if (my $nv = shift) {
        $self->{BULKFREQ} = $nv;
    }
    return $self->{BULKFREQ} || {};
}

sub is_connected {
    my $self = shift;
    return $self->{MAPLOC} || 0;
}

sub disconnect {
    my $self = shift;
    if (my $ml = $self->{MAPLOC}) {
        $ml->disconnect();
    }
}
sub base_url {
    my $self = shift;
    if (my $nv = shift) {
        $self->{BASE_URL} = $nv;
    }
    return $self->{BASE_URL};
}


sub flank {
    my $self = shift;
    if (my $nv = shift) {
        die "Unimplemented";
    }
    return $self->{FLANK};
}

sub impact_impact {
    my $self = shift;
    unless ($self->{IMPIMP}) {
        $self->bench_start();
        my $rv = $self->{IMPIMP} = {};
        map { $rv->{$_} = 'High' } $self->high_impacts();
        map { $rv->{$_} = 'Med'  } $self->medium_impacts();
        map { $rv->{$_} = 'Low'  } $self->low_impacts();
        $self->bench_end();
    }
    return $self->{IMPIMP};
}

sub impact_class {
    my $self   = shift;
    my $tok    = shift || "";
    my $impimp = $self->impact_impact();
    return $impimp->{$tok} || ($tok =~ /\,/ ? "Mult" : "Unk");
}

sub high_impacts {
    my $self = shift;
    if (my $nv = shift) {
        delete $self->{IMPIMP};
        die "Unimplemented";
    }
    return @{$self->{IMP_HIGH}};
}

*med_impacts = \&medium_impacts;
sub medium_impacts {
    my $self = shift;
    if (my $nv = shift) {
        delete $self->{IMPIMP};
        die "Unimplemented";
    }
    return @{$self->{IMP_MED}};
}

sub low_impacts {
    my $self = shift;
    if (my $nv = shift) {
        delete $self->{IMPIMP};
        die "Unimplemented";
    }
    return @{$self->{IMP_LOW}};
}

sub enumerate_tag {
    my $self = shift;
    my $type = $self->stnd_type(shift);
    if (defined $_[0]) {
        my $nt = $_[0];
        if ($nt) {
            push @{$self->{ENUM_TAG}{$type}}, $nt;
        } else {
            delete $self->{ENUM_TAG}{$type};
        }
    }
    return @{$self->{ENUM_TAG}{$type} || []};
}

sub _formater_count {
    my ($val) = @_;
    return '';
}


our $colDescription = {
    "Poly5"  => {
        fmt   => "AlnCen",
        cat   => "Variant Counting",
        width => 7,
        desc  => "Number of 1000 Genomes populations with MAF >= 5%",
    },
    "Poly10" => {
        fmt   => "AlnCen",
        cat   => "Variant Counting",
        width => 7,
        desc  => "Number of 1000 Genomes populations with MAF >= 10%",
    },
    "Poly30" => {
        fmt   => "AlnCen",
        cat   => "Variant Counting",
        width => 7,
        desc  => "Number of 1000 Genomes populations with MAF >= 30%",
    },

    "Class" => {
        cat   => "Information",
        width => 20,
        desc  => "High-level categories defined for populations. Could be tissues, disease states, etc",
    },
    "Categories" => {
        cat   => "Information",
        width => 20,
        desc  => "High-level categories assigned to variants. Generally labels the source(s) that defined the variant.",
    },
    "DetailFeature"   => {
        cat   => "Cancer Information",
        width => 20,
        desc  => "List of named features that cover at least part of the position",
    },

    "FullLocation" => {
        fmt   => "",
        cat   => "Genome Metrics",
        width => 15,
        desc  => "A single value representing the full location of a variant, eg 5.NCBI36:25..27 indicates chromosome 5, build NCBI36, 3 bases covering 25 through 27",
    },
    "Chr" => {
        fmt   => "AlnCen",
        cat   => "Genome Metrics",
        width => 7,
        desc  => "The chromosome name (eg '18' or 'X')",
    },
    "ChrPos" => {
        cat   => "Genome Metrics",
        width => 10,
        desc  => "A human-readable chromosome position. Single bases will be a single integer, multiple bases will be START..END, and deletions will be LEFT^RIGHT",
    },
    "Len" => {
        fmt   => "AlnCen",
        cat   => "Genome Metrics",
        width => 7,
        desc  => "The length of the reported object. For variants, the length references the reference allele, a value of zero indicates a deletion, 1 is a single nucleotide, and all others a multiple nucleotide location. Note that other alleles may have different lengths",
    },
    "Impact" => {
        cb    => \&_format_impact,
        fmt   => "AlnCen",
        cat   => "Impact of Variant on Gene",
        width => 9,
        desc  => "Impact code predicted for a variant",
    },
    "ImpRNA" => {
        fmt   => "",
        cat   => "Impact of Variant on Gene",
        width => 14,
        desc  => "Coordinate (transcript) and nucleotide sequence details for how a variant impacts the primary RNA (defined in column ImpactID)",
    },
    "ImpCDS" => {
        fmt   => "",
        cat   => "Impact of Variant on Gene",
        width => 14,
        desc  => "Coordinate (CDS reading frame) and nucleotide sequence details for how a variant impacts the primary RNA (defined in column ImpactID)",
    },
    "ImpProt" => {
        fmt   => "",
        cat   => "Impact of Variant on Gene",
        width => 14,
        desc  => "Coordinate (protein) and residue sequence details for how a variant impacts the primary protein (translated from the RNA in column ImpactID)",
    },
    "Accessions" => {
        fmt   => "",
        cat   => "Gene Information",
        width => 15,
        desc  => "Known accessions (IDs) of the variant. In almost all cases these will be dbSNP IDs, but could also include IDs from other data sources like Affymetrix or Incyte",
        IGNOREcb    => sub {
            # Ignoring this because Excel can not handle more than 65530 links
            # and will lock up the spreadsheet if this is exceeded
            # https://groups.google.com/forum/#!msg/spreadsheet-writeexcel/YXwMJj002J4/SMFf1dv70wwJ
            my $self = shift;
            my $v    = shift;
            my $f    = "";
            if ($v && $v =~ /^rs(\d+)/) {
                $v = [ 'url', "http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=$1", $v ];
                $f = 'url';
            }
            return ($f, $v);
        },
    },
    "ImpactID" => {
        fmt   => "",
        cat   => "Impact of Variant on Gene",
        width => 15,
        desc  => "The RNA accession used to calculate reported impacts. Please note that the nucleotide coordinates are ALMOST NEVER portable to other splice variants in the same gene, the protein coordinates are OFTEN NOT portable, and even the impact class (eg 'SYN') is OFTEN NOT conserved across splice forms",
    },
    "High" => {
        fmt   => "AlnCen",
        cat   => "Variant Counting",
        width => 7,
        desc  => "Number of variants that have a 'High' impact (see help sheet for list of impact tokens)",
    },
    "Med" => {
        fmt   => "AlnCen",
        cat   => "Variant Counting",
        width => 7,
        desc  => "Number of variants that have a 'Medium' impact (see help sheet for list of impact tokens)",
    },
    "Low" => {
        fmt   => "AlnCen",
        cat   => "Variant Counting",
        width => 7,
        desc  => "Number of variants that have a 'Low' impact (see help sheet for list of impact tokens)",
    },
    "Build" => {
        fmt   => "AlnCen",
        cat   => "Genome Metrics",
        width => 7,
        desc  => "The genomic build used for the analysis",
    },
    "ChrLft" => {
        fmt   => "",
        cat   => "Genome Metrics",
        width => 7,
        desc  => "The chromosome base position TO THE LEFT of the location. The locations 10, 15..17 and 45^46 will have widths respectively of 1, 3 and 0, and ChrLft = 9, 14 and 45. The last location is a deletion, and is why this nomenclature is used rather than ChrStart.",
    },
    "ChrRgt" => {
        fmt   => "",
        cat   => "Genome Metrics",
        width => 7,
        desc  => "The chromosome base position TO THE RIGHT of the location. The locations 10, 15..17 and 45^46 will have widths respectively of 1, 3 and 0, and ChrRgt = 11, 18 and 46. The last location is a deletion, and is why this nomenclature is used rather than ChrEnd.",
    },
    "Symbol" => {
        fmt   => "AlnCen",
        cat   => "Gene Information",
        width => 10,
        desc  => "The gene symbol, when available. In some cases multiple genes may overlap the same entry.",
    },
    "Reason" => {
        fmt   => "",
        cat   => "Gene Information",
        width => 25,
        desc  => "This column is used in the 'Ignored' worksheet, and will briefly summarize the reason that the locus was excluded from the analysis.",
    },
    "maploc_key" => {
        fmt   => "",
        cat   => "Geek Information",
        width => 12,
        desc  => "A unique identifier used within the MapLoc database system",
    },
    "AllImp" => {
        cb    => \&_format_impact,
        fmt   => "",
        cat   => "Impact of Variant on Gene",
        width => 15,
        desc  => "All observed impacts (as three letter tokens) at this location, for every RNA considered, both RefSeq and Ensembl",
    },
    "RefImp" => {
        cb    => \&_format_impact,
        fmt   => "",
        cat   => "Impact of Variant on Gene",
        width => 15,
        desc  => "All observed impacts (as three letter tokens) at this location, for every RefSeq transcripts considered",
    },
    "EnsImp" => {
        cb    => \&_format_impact,
        fmt   => "",
        cat   => "Impact of Variant on Gene",
        width => 15,
        desc  => "All observed impacts (as three letter tokens) at this location, for every Ensembl transcript considered",
    },
    "RnaLen" => {
        fmt   => "AlnCen",
        cat   => "Gene Information",
        width => 7,
        desc  => "The reported length of the RNA. This is the number of bases of RNA sequence, which may be larger than the number of bases aligned to the genome and reported in RnaUsed",
    },
    "RnaUsed" => {
        fmt   => "AlnCen",
        cat   => "Gene Information",
        width => 7,
        desc  => "The total number of RNA bases considered. For RNAs, this will be the length of the RNA that could be aligned to the genome (and thus intersect with variants). For genes, it will be the sum of discrete bases in all transcripts (so a base that is in three splice variants counts only once)",
    },
    "ProtLen" => {
        fmt   => "AlnCen",
        cat   => "Gene Information",
        width => 7,
        desc  => "The total number of protein coding residues considered. For RNAs, this will be the number of amino acids in the CDS. For genes, it will be the sum of discrete residues across all coding segments (so a residue that is in three splice variants counts only once).",
    },
    "LL" => {
        fmt   => "",
        cat   => "Gene Information",
        width => 15,
        desc  => "The Entrez (aka LocusLink) accession(s) associated with the entry. Should be unique for Genes and RNAs, can have zero or more for variants",
        IGNOREcb    => sub {
            my $self = shift;
            my $v    = shift;
            my $f    = "";
            if ($v) {
                my $esc = $self->esc_url($v);
                $v = [ 'url', "http://genetracker.pri.bms.com/mev-gt/MEV-GT/#rootId/$esc", $v ];
                $f = 'url';
            }
            return ($f, $v);
        },
    },
    "RSR" => {
        fmt   => "",
        cat   => "Gene Information",
        width => 15,
        desc  => "The RefSeq RNA accession(s) associated with an entry. Loci can have zero (in the case of 'phenotypic' loci) or more, variants can have zero or more.",
    },
    "RSP" => {
        fmt   => "",
        cat   => "Gene Information",
        width => 15,
        desc  => "The RefSeq protein accession(s) associated with an entry. Loci can have zero (in the case of non-coding genes or 'phenotypic' loci) or more, variants can have zero or more.",
    },
    "ENSG" => {
        fmt   => "",
        cat   => "Gene Information",
        width => 17,
        desc  => "The Ensembl gene accession(s) associated with an entry, can have zero or more values. The value is derived from the RefSeq RNA used to initiate data recovery.",
    },
    "ENST" => {
        fmt   => "",
        cat   => "Gene Information",
        width => 17,
        desc  => "The Ensembl transcript accession(s) associated with an entry, can have zero or more values. The value is derived from the RefSeq RNA used to initiate data recovery.",
    },
    "ENSP" => {
        fmt   => "",
        cat   => "Gene Information",
        width => 17,
        desc  => "The Ensembl protein accession(s) associated with an entry, can have zero or more values. The value is derived from the RefSeq RNA used to initiate data recovery.",
    },
    "ID" => {
        fmt   => "",
        cat   => "Geek Information",
        width => 15,
        desc  => "The primary accession for the entry. Most calculated values will be based off of this accession, which will be an NCBI RefSeq transcript identifier for RNAs, and an NCBI Entrez gene identifier for genes.",
        cb    => sub {
            my $self = shift;
            my $v    = shift;
            my $f    = "";
            if ($v && $v =~ /^LOC/) {
                my $esc = $self->esc_url($v);
                my $baseUrl = $self->base_url();
                $v = [ 'url', "$baseUrl/mapLocReporter.pl?showrna=$esc", $v ];
                $f = 'url';
            }
            return ($f, $v);
        },
    },
    "Footprint" => {
        fmt   => "",
        cat   => "Genome Metrics",
        width => 25,
        desc  => "The genomic location(s) occupied by the entity. The format will be CHR:Start..Stop [Score%]. Some RNAs and genes can not be unambiguously mapped to a unique location. These will have multiple positions reported.",
        IGNOREcb    => sub {
            my $self = shift;
            my $v    = shift;
            my $f    = "";
            if ($v) {
                my $esc = $v;
                $esc =~ s/\s*\[[^\]]+\]\s*/ /g;
                $esc = $self->esc_url($esc);
                my $baseUrl = $self->base_url();
                $v = [ 'url', "$baseUrl/mapLocReporter.pl?footprint=$esc", $v ];
                $f = 'url';
            }
            return ($f, $v);
        },
    },
    "Copies" => {
        fmt   => "AlnCen",
        cat   => "Genome Metrics",
        width => 7,
        desc  => "The number of genomic locations occupied by the entity. See the Footprint column for precise details",
    },
    "NumRNA" => {
        fmt   => "AlnCen",
        cat   => "Genome Metrics",
        width => 7,
        desc  => "The number of distinct RefSeq RNAs associated with a gene. See the RSR column for specific accessions",
    },
    "Description" => {
        fmt   => "",
        cat   => "Gene Information",
        width => 25,
        desc  => "Human-readable description of the object.",
    },
    "TargetClass" => {
        fmt   => "",
        cat   => "Gene Information",
        width => 15,
        desc  => "Each target class (listed in the 'Target Classification' section) will have its own column, but to aid filtering all matching target classes will be collected in this column as well.",
    },
    "" => {
        fmt   => "",
        cat   => "",
        width => 7,
        desc  => "",
    },
    "" => {
        fmt   => "",
        cat   => "",
        width => 7,
        desc  => "",
    },
};

our @geneDetails = qw(LL RSR RSP Symbol ENSG ENST ENSP);
sub rna_detail_table {
    my $self = shift;
    my $args = $self->parseparams( @_ );
    my $rnaAcc = $args->{RNA};
    my $ml     = $self->maploc();
    my $rna;
    if (ref($rnaAcc)) {
        $rna    = $rnaAcc;
        $rnaAcc = $rna->accession;
    } else {
        $rna = $self->rna_obj( $rnaAcc );
    }
    my $rtags    = $rna->all_tag_values();
    my $rv = { ID => $rnaAcc, maploc_key => $rnaAcc };
    my $used;
    if (my $alns = $args->{ALNS}) {
        # We can calculate the footprint
        foreach my $aln (@{$alns}) {
            push @{$rv->{Footprint}}, $aln->chr_footprint();
        }
        my @genCDS = $ml->HSPs_to_genomic_footprint( $alns );
        $used = 0;
        map { $used += $_->length() } @genCDS;
    }
    $rv->{RnaUsed} = $used;
    # GO terms are anchored off the unversioned accession
    my $unv = $rnaAcc; $unv =~ s/\.\d+$//;
    my $obj = $ml->get_text($unv);
    $self->_add_go_and_desc($rv, $obj);
    foreach my $col (@geneDetails) {
        $rv->{$col} = $rna->comma_tag_value( $col );
    }

    my $vars = $self->tally_locations
        ( -locs => $args->{LOCS}, -rna => $rnaAcc );
    my ($cs, $ce) = $rna->cds();
    my $aa = $rv->{ProtLen} = $cs ? int(($ce - $cs + 1) / 3) : 0;
    my $bp = $rv->{RnaLen}  = $rna->length();
    $self->collapse_locations( $rv, $vars, $used || $bp, $aa );
    return [ $rv ];
}

sub _add_go_and_desc {
    my $self = shift;
    my ($rv, $obj) = @_;
    $obj->read();
    my $dc = $self->druggable_classes();
    while (my ($go, $gi) = each %{$dc}) {
        if (my $v = $obj->unique_value( $go )) {
            my $tok = $gi->{tok};
            $rv->{$tok} = $v;
            $rv->{TargetClass}{$tok}++;
        }
    }
    foreach my $col (qw(Description)) {
        $rv->{$col} = $obj->comma_tag_value( $col );
    }
}

sub locus_detail_table {
    my $self = shift;
    my $args  = $self->parseparams( @_ );
    my $locid = $args->{LOCUS} || $args->{ID};
    my $rnas  = $args->{RNAS} || $args->{RNA};
    my $locs  = $args->{LOCS};
    my $build = $args->{BUILD} || 'GRCh37';
    my $ralns = $args->{ALNS}  || {};
    my $ml    = $self->maploc();
    $rnas   ||= $ml->get_rnas_for_locus( $locid );
    my $rv    = { ID         => $locid, 
                  maploc_key => $locid,
                  NumRNA     => $#{$rnas} + 1 };
    $self->rna_annotate_row( $rv, $rnas );
    my $obj   = $ml->get_text($locid);
    $self->_add_go_and_desc($rv, $obj);    
    my $rkey  = $locid =~ /^ENSG/ ? 'ENST' : 'RSR';
    my @raccs = $obj->tag_values( $rkey );
    # warn $self->branch(\@raccs);
    # We will need to project the RNA onto the genome to tally up the
    # non-redundant bases and coding positions
    my (@ests, @cdss, %notes);
    foreach my $r (@{$rnas}) {
        my $rna  = $self->rna_obj($r);
        my $racc = $rna->accession();
        my $raid = $rna->accession_id();
        my $alns = $ralns->{$racc};
        unless ($alns) {
            $alns = [];
            foreach my $aln ( sort { $b->score <=> $a->score } $rna->
                              each_alignment_for_build( $build )) {
                last if ($aln->howbad());
                push @{$alns}, $aln;
            }
        }
        $alns = [values %{$alns}] if (ref($alns) eq 'HASH');
        if ($#{$alns} == -1) {
            # $self->msg("$locid : No genomic alignments for $racc");
            push @{$notes{"No genome locations"}}, $racc;
            next;
        }
        my $rng  = $rna->cds_range();
        foreach my $aln (@{$alns}) {
            
            push @ests, $aln;
            my $cid = $aln->other_seq_id( $raid );
            if ($rng) {
                my $int = $aln->intersect
                    ( -hsp => $rng, -anchor => $rng->anchor_id() );
                push @cdss, $int;
            }
        }
    }
    foreach my $nt (sort keys %notes) {
        push @{$rv->{Notes}}, "$nt: ".join(",", @{$notes{$nt}});
    }
    # Tally up the unique genomic bases used for all the RNAs, and all
    # their coding segments:
    my ($aa, $bp) = (0,0);
    my @geno   = $ml->HSPs_to_genomic_footprint( \@ests );
    foreach my $gn (@geno) {
        $bp += $gn->length();
        push @{$rv->{Footprint}}, $gn->chr_footprint();
    }
    unless ($locs) {
        my $flank = $args->{FLANK} || $self->flank();
        $locs = $ml->locations_for_HSPs( \@geno, $flank );
    }
    my @genCDS = $ml->HSPs_to_genomic_footprint( \@cdss );
    map { $aa += $_->length() } @genCDS;
    $aa = int($aa/3);
    $rv->{ProtLen} = $aa || undef;
    $rv->{RnaUsed} = $bp || undef;
    $rv->{Copies}  = $#geno + 1;
    my $vars = $self->tally_locations
        ( -locs => $locs, -rnas => \@raccs, -chooserna => 1 );
    $self->collapse_locations( $rv, $vars, $bp, $aa );
    return [ $rv ];
}

sub collapse_locations {
    my $self = shift;
    my ($row, $vars, $bp, $aa, $useRNA) = @_;
    my $impH = $self->impact_impact();
    my (%aa2Imp, %hotspot, %aaImps);
    my @promote = keys %{$self->{PROMOTE_RNA}};
    my @tally   = keys %{$self->{TALLY_COL}};
    foreach my $var (@{$vars}) {
        if ($var->{Poly5}) {
            $row->{Poly5}++;
        }
        while (my ($t, $pH) = each %{$var->{_TALLY_MEMBER} ||= {}}) {
            map { $row->{_TALLY_MEMBER}{$t}{$_} = undef } keys %{$pH};
        }
        foreach my $tag (@promote) {
            if (my $h = $var->{$tag}) {
                map { $row->{$tag}{$_} = 1 } keys %{$h};
            }
        }
        my $imp    = $var->{Impact} || "UNK";
        my $impImp = $impH->{$imp} || 'Low';
        my $pLoc   = $var->{protLocation};
        if ($impImp ne 'Low') {
            if ($pLoc) {
                my ($plTok, $plID, $plPos) = @{$pLoc};
                foreach my $t (@tally) {
                    my $tn = $var->{$t} || 0;
                    $aaImps{$t} += $tn;
                    $aa2Imp{$t}{ $plTok } += $tn;
                }
            }
            if ($var->{Len} <= 1) {
                my $l   = $var->{ChrLft} - 3;
                my $r   = $var->{ChrRgt} + 3;
                my @obs;
                foreach my $tag (@promote) {
                    if (my $h = $var->{$tag}) {
                        push @obs, keys %{$h};
                    }
                }
                for my $pos ($l..$r) {
                    map { $hotspot{$pos}{$_} = 1 } @obs;
                }
            }
        }
        $row->{$impImp}++;
        $row->{$imp}++;
    }
    while (my ($t, $pH) = each %{$row->{_TALLY_MEMBER} ||= {}}) {
        my @u = keys %{$pH};
        $row->{$t} = ($#u + 1) || undef;
    }
    # warn $self->branch($row->{_TALLY_MEMBER});
    foreach my $obsHash (values %hotspot) {
        my @obs = keys %{$obsHash};
        my $onum = $#obs + 1;
        foreach my $lvl (3, 6) {
            last if ($onum < $lvl);
            $row->{"HotSpot$lvl"}++;
        }
    }
    # Add a pseudo count to polys:
    my $polys = ($row->{Poly5} || 0) + 1;
    my @derivedKeys;
    foreach my $t (@tally) {
        if ($aa) {
            # Count the number of residues with more than 2 observations
            my $res = 0;
            map { $res++ if ($_ > 1) } values %{$aa2Imp{$t}};
            if ($res) {
                my $key = $t."-aaImp";
                push @derivedKeys, $key;
                $row->{$key} = int(0.5 + 10000 * $res / $aa)/100;
            }
            if (my $aai = $aaImps{$t}) {
                my $key = $t."/100aa";
                push @derivedKeys, $key;
                $row->{$key}  = int(0.5 + 10000 * $aai / $aa)/100;
            }
        }
        if (my $tal = $row->{$t}) {
            my $key = $t."/Poly";
                push @derivedKeys, $key;
            $row->{$key} = $polys ? int(0.5 + 100 * $tal / $polys)/100 : 999;
            if ($bp) {
                my $key = $t."/kb";
                push @derivedKeys, $key;
                $row->{$key} = int (0.5 + 100000 * $tal / $bp) / 100;
            }
        }
    }
    foreach my $key (@derivedKeys) {
        unless ($self->{DERIVED_TAG}{$key}) {
            $self->{DERIVED_TAG}{$key} = 1;
        }
    }
}

our $refPrefs = { NM => 1, NR => 2, XM => 3, XR => 4 };
sub tally_locations {
    my $self = shift;
    my $args     = $self->parseparams( @_ );
    my $locs     = $args->{LOCS} || [];
    my $rnaAcc   = $args->{RNA};
    my $chooseRNA = $args->{CHOOSERNA};
    my $bulkFreq  = $self->bulk_frequencies();
    my $ml        = $self->maploc();
    my @rows;
    if ($rnaAcc) {
        if (my $r = ref($rnaAcc)) {
            if ($r eq 'ARRAY') {
                $rnaAcc = $#{$rnaAcc} == -1 ? undef : 
                { map { $_ => 1 } @{$rnaAcc} };
            }
        } else {
            $rnaAcc = { $rnaAcc => 1 };
        }
    }
    my %showPop = map { $_ => 1 } keys %{$self->{POP_COLS} || {}};
    foreach my $loc (@{$locs}) {
        my %row; # = %{$defRow};
        my $r = ref($loc);
        my $cx;
        if ($r eq 'ARRAY') {
            ($loc, $cx) = @{$loc};
        }
        # $cx ||= $loc->cx_genome_part( -howbad => 0 );
        # $self->msg("[VERBOSE]", $loc->full_loc());
        my $chr       = $row{Chr}    = $loc->chr;
        my $chrL      = $row{ChrLft} = $loc->lft;
        my $chrR      = $row{ChrRgt} = $loc->rgt;
        my @accs      = $loc->each_accession();
        $row{sortKey} = sprintf("%3s %9d %9d", $chr, $chrL, $chrR);
        $row{Build}   = $loc->build();
        $row{Len}     = $loc->length();
        $row{FullLocation} = $loc->full_loc();
        $row{ChrPos}       = $ml->pretty_flanks($chrL, $chrR);
        my $lid = $row{maploc_key}   = $loc->pkey();
        $row{Accessions}   = \@accs;
        map { $row{Categories}{$_} = 1 } $loc->each_category();
        my ($impData, $popData, $featData);
        if ($cx) {
            # Get already-calculated data from CX object
            $impData  = $cx->{impact};
            $popData  = $cx->{freqs};
            $featData = $cx->{features};
            $row{Impact} = $cx->{impToken};
        } else {
            # Need to calculate impact now
            $impData = $loc->impact( @_ );
            $popData = $loc->each_population_id();
            my $oi   = $loc->overall_impact
                ( -popdata => $popData, -impdata => $impData, @_ );
            $featData = $oi->{features};
            $row{Impact} = $oi->{impToken};
        }
        foreach my $fid (keys %{$featData}) {
            my $fn = $ml->cached_pkey_to_text( $fid );
            my $obj = $ml->get_text($fn);
            $obj->read();
            my ($name) = $obj->tag('Name');
            $row{Feature}++;
            $row{DetailFeature}{ $name || $fn }++;
        }
        foreach my $tag ($self->enumerate_tag('Variant')) {
            map { $row{$tag}{ $_ }++ } $loc->tag_values($tag);
        }
        my @impIds   = keys %{$impData};
        my $racc     = $rnaAcc;
        my (@allAccs, @chooseAccs);
        foreach my $rid (@impIds) {
            my $rna = $ml->get_rna_by_id($rid);
            my $ver = $rna->acc();
            push @allAccs, $rid;
            if ($racc) {
                my $unv = $ver; $unv =~ s/\.\d+$//;
                next unless ($racc->{$unv} || $racc->{$ver});
            }
            push @chooseAccs, $ver;
        }
        if (!$racc || $chooseRNA) {
            # Find the 'best' RNA to use
            # Prefer RefSeq to Ensembl, using $refPrefs
            # Prefer low number IDs to high numbered
            my @sorter;
            foreach my $acc (@chooseAccs) {
                if ($acc =~ /^ENS[A-Z]+(\d+)/) {
                    push @sorter, [$acc, 10, $1 + 0];
                } elsif ($acc =~ /^([NX][MR])_(\d+)/) {
                    push @sorter, [$acc, $refPrefs->{$1} || 9, $2 + 0];
                } else {
                    push @sorter, [$acc, 20, 1];
                }
            }
            my ($best) = sort { $a->[1] <=> $b->[1] || 
                                    $a->[2] <=> $b->[2] || 
                                    $a->[0] cmp $b->[0] } @sorter;
            $racc = { $best ? $best->[0] : 'N/A', 1 };
        }
        $self->rna_annotate_row( \%row, \@allAccs, 'isID' );
        my @protSorter;
        foreach my $idat (values %{$impData}) {
            my $rid = $idat->{rid};
            my $rna = $ml->get_rna_by_id($rid);
            my $ver = $rna->acc();
            my $unv = $ver; $unv =~ s/\.\d+$//;
            my $imp = $idat->{imp} || "UNK";
            $row{AllImp}{$imp} = 1;
            if ($racc->{$ver} || $racc->{$unv}) {
                $row{Impact}   = $imp;
                $row{ImpRNA}   = $idat->{nucNom};
                $row{ImpCDS}   = $idat->{cdsNom};
                $row{ImpProt}  = $idat->{protNom};
                $row{ImpactID} = $ver;
            }
            if ($ver =~ /^ENS/) {
                $row{EnsImp}{$imp} = 1;
            } elsif ($ver =~ /^([NX][MR])_(\d+)/) {
                my ($tok, $rsNum) = ($1, $2);
                $row{RefImp}{$imp} = 1;
                if (my $p = $idat->{protPos}) {
                    $p = sprintf("%010d", $p) if ($p =~ /^\d+$/);
                    push @protSorter, 
                    [ sprintf("%010d %s", $rsNum, $p),
                      $refPrefs->{$tok} || 9,
                      $rsNum, $p ];
                }
            }
        }
        if ($#protSorter != -1) {
            # Try to find a token that can represent the protein position
            # across multiple transcripts
            my ($using) = sort { $a->[1] <=> $b->[1] || 
                                     $a->[0] cmp $b->[0] } @protSorter;
            $row{protLocation} = [ $using->[0], $using->[2], $using->[3] ];
        }

        
       # $self->maploc->prebranch($bulkFreq->{$lid});
        while (my ($pid, $freqs) = each %{$popData}) {
            my $dat  = $self->population_info( $pid );
            foreach my $tally (@{$dat->{INCS}}) {
                $row{$tally}++;
                $row{_TALLY_MEMBER}{$tally}{$pid} = 1;
            }
            while (my ($k, $val) = each %{$dat->{INC2}}) {
                $row{$k}{ $val }++ if ($val);
            }
            my $lf;
            if (my $bfd = $bulkFreq->{$lid}) {
                $lf = $bfd->{$pid};
            }
            if (defined $lf) {
                if ($showPop{$pid}) {
                    $row{$dat->{name}} = int(0.5 + 10000 * $lf) / 10000;
                }
                if ($dat->{IS_SNP}) {
                    # my $mm = 100 * $self->max_minor_from_freq($freqs);
                    my $mm = int(0.5 + 1000 * $lf) / 10 unless ($lf < 0);
                    # $self->msg("$pid = $mm");
                    
                    # TODO: We should not run the tally loop below
                    # if the population is a high-level parent
                    
                    foreach my $f (5,10,30) {
                        last if ($f > $mm);
                        $row{"Poly$f"}++;
                    }
                }
            }
        }
        push @rows, \%row;
    }
    @rows = sort { $a->{sortKey} cmp $b->{sortKey} } @rows;
    return wantarray ? @rows : \@rows;
}

sub add_population_column {
    my $self = shift;
    if (my $req = shift) {
        my $rm = shift;
        my $ml = $self->maploc;
        if (my $pop = $ml->get_population_fast( $req )) {
            my $pid = $pop->pkey();
            my $dat = $self->population_info( $pid );
            if ($rm) {
                delete $dat->{HAS_COL};
                delete $self->{POP_COLS}{$pid};
            } else {
                $dat->{HAS_COL} = 1;
                $self->{POP_COLS}{$pid} = 1;
            }
        } else {
            $self->msg("[!]","Failed to find population '$req'");
        }
    }
}

sub population_info {
    my $self = shift;
    my $pid  = shift;
    unless ($self->{POP_INFO}{$pid}) {
        $self->bench_start();
        $self->{POP_INFO}{$pid} = {};
        my $ml = $self->maploc;
        if (my $pop = $ml->get_population($pid)) {
            $pop->read();
            my $tags  = $pop->all_tag_values();
            my $par   = $pop->parent();
            my $root  = $pop->root_parent();
            my $name  = $pop->name();
            my $dat   = $self->{POP_INFO}{$pid} = { name => $name };
            my $rname = $dat->{ROOT} = $root ? $root->name() : "";
            my $incs  = $dat->{INCS} = [];
            my $inc2  = $dat->{INC2} = {};
            foreach my $tal (@{$tags->{Tally} || []}) {
                push @{$incs}, $tal;
                $self->{XTRA_COL}{$tal}  = 1;
                $self->{TALLY_COL}{$tal} = 1;
            }
            foreach my $cat (@{$tags->{Category} || []}) {
                my $cobj = $ml->get_text($cat);
                $cobj->read();
                my $ctag = $cobj->all_tag_values();
                foreach my $ct (@{$ctag->{"SNP Class Tag"} || []}) {
                    if (my $ot = $self->_obj_tag($tags, $ct)) {
                        $inc2->{$ct} = $ot;
                        $self->{XTRA_COL}{$ct} = 1;
                        $self->{PROMOTE_RNA}{$ct} = 1;
                    }
                }
                foreach my $ct (@{$ctag->{"Enumerate"} || []}) {
                    if (my $ot = $self->_obj_tag($tags, $ct)) {
                        $inc2->{$ct} = $ot;
                        $self->{XTRA_COL}{$ct} = 1;
                        $self->{PROMOTE_RNA}{$ct} = 1;
                    }
                }
            }
            # $self->msg_once($rname);
            if ($rname eq $snpPopName || $name eq $snpPopName) {
                # push @{$incs}, 'MAF1kg';
                $dat->{IS_SNP} = 1;
            }
            # warn $self->branch( $dat ) unless ($dat->{IS_SNP});
        }
        $self->bench_end();
    }
    return $self->{POP_INFO}{$pid};
}

sub rna_obj {
    my $self = shift;
    my $acc  = shift || "";
    my $isID = shift;
    my $rna;
    if (ref($acc)) {
        $rna = $acc;
        $acc = $rna->accession();
    }
    unless (defined $self->{RNAOBJS}{$acc}) {
        if ($isID) {
            $rna = $self->maploc->get_rna_by_id($acc);
            weaken( $self->{RNAOBJS}{$acc} = $rna );
            $acc = $rna->accession();
        } if (!$rna) {
            ($rna) = $self->maploc->get_rnas( $acc );
        }
        if ($rna) {
            $rna->read();
            foreach my $tag (@geneDetails) {
                my @tv = $rna->tag_values( $tag );
                $rna->{SR_TAGS}{$tag} = \@tv unless ($#tv == -1);
            }
        }
        while ($#{$self->{RNACACHE}} > 10000) {
            shift @{$self->{RNACACHE}};
        }
        weaken( $self->{RNAOBJS}{$acc} = $rna );
        push @{$self->{RNACACHE}}, $rna;
    }
    return $self->{RNAOBJS}{$acc};
}

sub rna_annotate_row {
    my $self = shift;
    my ($row, $accs, $isID) = @_;
    foreach my $acc (@{$accs}) {
        if (my $rna = $self->rna_obj( $acc, $isID )) {
            my $srt = $rna->{SR_TAGS};
            while (my ($tag, $tvs) = each %{$rna->{SR_TAGS}}) {
                map { $row->{$tag}{$_} = 1 } @{$tvs};
            }
            # warn $rna->to_text() unless ($rna->{CT__DONE}++);
        }
    }
}

sub _obj_tag {
    my $self = shift;
    my ($tags, $tag) = @_;
    return $tags->{$tag}[0] if ($tags->{$tag});
    return "";
}

sub max_minor_from_freq {
    my $self = shift;
    my $freqs = shift;
    # $self->maploc->prebranch($freqs);
    my @f = sort { $b <=> $a } map { $_->[0] || 0 } values %{$freqs};
    return $f[1] || 0;
}

sub table_to_tsv {
    my $self = shift;
    my $args = $self->parseparams( @_ );
    my $rows = $args->{ROWS};
    return "" if (!$rows || $#{$rows} == -1);
    my $cols = $self->_requested_cols( @_ );
    my $tsv  = $args->{NOHEAD} ? '' : join("\t", @{$cols})."\n";
    foreach my $row (@{$rows}) {
        my $r = $self->_flatten_row( $row, $cols );
        $tsv .= join("\t", @{$r})."\n";
    }
    return $tsv;
}

our $joiner = ',';
sub _flatten_row {
    my $self = shift;
    my ($row, $cols) = @_;
    my @r;
    foreach my $c (@{$cols}) {
        my $v = $row->{$c};
        if (!defined $v) {
            push @r, '';
        } elsif (my $rf = ref($v)) {
            if ($rf eq 'ARRAY') {
                push @r, join($joiner, @{$v});
            } elsif ($rf eq 'HASH') {
                push @r, join($joiner, sort keys %{$v});          
            } else {
                push @r, "??$rf??";
            }
        } else {
            push @r, $v;
        }
    }
    return \@r;
}

sub table_to_text {
    my $self = shift;
    my $args = $self->parseparams( @_ );
    my $rows = $args->{ROWS};
    return "" if (!$rows || $#{$rows} == -1);
    my $cols = $self->choose_cols( @_ );
    my @tab  = ( [@{$cols}] );
    foreach my $row (@{$rows}) {
        my $r = $self->_flatten_row( $row, $cols );
        push @tab, $r;
    }
    my $fmt = "|";
    my $hd  = "+";
    for my $i (0..$#{$cols}) {
        my ($len) = sort { $b <=> $a } map { length($_->[$i]) } @tab;
        if (my $hpad = int(($len - length($tab[0][$i]))/2)) {
            $tab[0][$i] .= ' ' x $hpad;
        }
        $fmt .= sprintf(" %%%ds |", $len);
        $hd  .= ('-' x ($len+2)) .'+';
    }
    $fmt .= "\n";
    $hd  .= "\n";
    my $txt = $hd;
    for my $r (0..$#tab) {
        $txt .= sprintf($fmt, @{$tab[$r]});
        $txt .= $hd unless ($r % 5);
    }
    $txt .= $hd if ($#tab % 5);
    return $txt;
}

sub table_to_excel {
    my $self = shift;
    my $args = $self->parseparams( @_ );
    my $rows = $args->{ROWS};
    return undef if (!$rows || $#{$rows} == -1);
    my $cols = $self->choose_cols( @_ );
    my $eh   = $self->{EH} || $self->{EXCEL};
    $eh ||= $self->default_excel();
}

sub _requested_cols {
    my $self = shift;
    my $args = $self->parseparams( @_ );
    my $cols = $args->{COLS} || $self->choose_cols( @_ );
    my %got  = map { $_ => 1 } @{$cols};
    foreach my $xtra ( sort keys %{$self->{XTRA_COL}} ) {
        push @{$cols}, $xtra unless ($got{$xtra});
    }
    return $cols;
}

sub stnd_type {
    my $self = shift;
    my $req  = lc(shift || "");
    if ($req =~ /(gene|locus)/) {
        return 'Gene';
    } elsif ($req =~ /(rna|trans)/) {
        return 'RNA';
    } elsif ($req =~ /(var|poly|loc)/) {
        return 'Variant';
    }
    return "";
}

sub choose_cols {
    my $self = shift;
    my $args = $self->parseparams( @_ );
    my $cols = $args->{COLS};
    return $cols if ($cols);
    if (my $type = lc($args->{TYPE} || '')) {
        my $styp = $self->stnd_type( $type );
        my $ml   = $self->maploc();
        my @imps = $ml->impact_list();
        my @xtra = sort keys %{$self->{XTRA_COL}};
        push @xtra, sort $self->enumerate_tag( $styp );
        if ($styp eq 'RNA' || $styp eq 'Gene') {
            my @rv = qw(ID Symbol);
            push @rv, @xtra;
            if (my $dt = $self->{DERIVED_TAG}) {
                push @rv, sort keys %{$dt};
            }
            # push @rv, qw(HotSpot3 HotSpot6);
            push @rv, qw(Description Poly5 High Med Low);
            push @rv, @imps if ($styp eq 'RNA');
            push @rv, 'RnaUsed';
            push @rv, 'RnaLen' unless ($styp eq 'Gene');
            push @rv, qw(ProtLen);
            push @rv, qw(Copies NumRNA) if ($styp eq 'Gene');
            push @rv, "TargetClass";
            if ($args->{ADDGO}) {
                my $dc = $self->druggable_classes();
                my @goToks = sort { uc($a) cmp uc($b) } 
                map { $_->{tok} } values %{$dc};
                push @rv, @goToks;
            }
            push @rv, qw(LL RSR RSP ENSG ENST ENSP Footprint maploc_key Notes);
            return \@rv;
        } elsif ($styp eq 'Variant') {
            my @rv = qw(FullLocation Chr ChrPos Len Symbol Categories);
            push @rv, @xtra;
            push @rv, qw(Accessions LL Impact ImpRNA ImpCDS ImpProt ImpactID);
            my @popNames;
            foreach my $pid (keys %{$self->{POP_COLS} || {}}) {
                my $dat = $self->population_info( $pid );
                my $name = $dat->{name};
                $colDescription->{$name} = {
                    fmt   => "AlnCen",
                    cat   => "Variant Counting",
                    width => 7,
                    desc  => "Minor allele frequency for population $name",
                    cb    => \&_format_freq,
                };
                push @rv, $name;
            }
            push @rv, sort @popNames;
            push @rv, qw(Poly5 Poly30 AllImp RefImp EnsImp ChrLft ChrRgt Build DetailFeature maploc_key Notes);
            return \@rv;
        } else {
            $self->death("Can not interpret table type '$type'");
        }
    }
    if (my $rows = $args->{ROWS}) {
        return [ sort keys %{$rows->[0]} ];
    }
    $self->death("Unable to find column names");
}

sub _format_go {
    my $self = shift;
    my $v    = shift;
    my $f    = "";
    if ($v) {
        my $sc = -1;
        if ($v =~ /^(\S+) (\S+)/) {
            ($sc, $v) = (100 * $1, $2);
        }
        if ($sc < 0) {
            $f .= 'Unknown';
        } elsif ($sc < 90) {
            $f .= 'Poor';
        } elsif ($sc < 100) {
            $f .= 'Good';
        } else {
            $f .= 'Perfect';
        }
    }
    return ($f, $v);
}

sub _format_impact {
    my $self = shift;
    my $v    = shift;
    my $f    = $v ? $self->impact_class( $v )."Impact" : "";
    return ($f, $v);
}

sub _format_freq {
    my $self = shift;
    my $v    = shift;
    return (undef, $v) if (!defined $v || $v eq "");
    my $i    = int(0.5 + 10 * $v);
    my $name = sprintf("GrayScale%d", $i);
    return ($name, $v);
}

sub druggable_classes {
    my $self = shift;
    unless ($self->{DC_GO}) {
        my $ml   = $self->maploc();
        my $goPar = $ml->get_text("Druggable GeneOntology");
        $goPar->read();
        my %gos;
        foreach my $go ($goPar->tag_values('Basic')) {
            my $obj = $ml->get_text($go);
            $obj->read();
            my $tok = $obj->simple_value('Token');
            unless ($tok) {
                $self->msg("[?]","No token for $go");
                next;
            }
            my @bits;
            if (my $short = $obj->simple_value('Short')) {
                push @bits, $short;
            }
            push @bits, $go;
            if (my $dt = $obj->simple_value('Description')) {
                push @bits, $dt;
            }
            my $desc = join(' : ', @bits);
            $gos{$go} = {
                desc => $desc,
                tok  => $tok,
                id   => $go,
            };
            $colDescription->{$tok} = {
                cb    => \&_format_go,
                fmt   => "",
                cat   => "Target Classification",
                width => 7,
                desc  => $desc,
            };
        }
        $self->{DC_GO} = \%gos;
    }
    return $self->{DC_GO};
}

sub excel_helper {
    my $self = shift;
    my $args = $self->parseparams( @_ );
    my $path = $args->{PATH};
    my $eh = BMS::ExcelHelper->new( $path );
    my $ml   = $self->maploc();
    my @imps = $ml->impact_list();
    foreach my $imp (@imps) {
        my $info = $ml->impact_details($imp);
        my $cd = $colDescription->{$imp} = {
            fmt   => "AlnCen",
            cat   => "Impact Codes",
            width => 7,
            desc  => "Number of locations predicted to have $info->{name} impact : $info->{desc}",
        };
    }

    $eh->format( -name  => "AlnCen",
                 -bold  => 0,
                 -align => 'center' );

    $eh->format( -name  => "HighImpact",
                 -bg_color => 'pink',
                 -bold  => 0,
                 -align => 'center' );
    $eh->format( -name  => "MedImpact",
                 -bg_color => 'yellow',
                 -bold  => 0,
                 -align => 'center' );
    $eh->format( -name  => "LowImpact",
                 -bg_color => 'lime',
                 -bold  => 0,
                 -align => 'center' );
    $eh->format( -name  => "MultImpact",
                 -bg_color => 27,
                 -bold  => 0,
                 -align => 'center' );
    $eh->format( -name  => "UnkImpact",
                 -bg_color => 'silver',
                 -bold  => 0,
                 -align => 'center' );

    for my $i (0..10) {
        my $ind = 25 * (10 - $i);
        my $cc  = $eh->set_custom_color(40 + $i, $ind, $ind, $ind);
        my $fg  = $i <= 5 ? 'black' : 'white';
        my $name = sprintf("GrayScale%d", $i);
        $eh->format( -name       => $name,
                     -align      => 'center',
                     -color      => $fg,
                     -bg_color   => $cc );
    }

    my $bgGo = { Unknown => 'gray',
                 Poor => 'pink',
                 Good => 'yellow',
                 Perfect => 'lime'};
    while (my ($name, $color) = each %{$bgGo}) {
        $eh->format( -name     => $name,
                     -bg_color => $color,
                     -align    => 'center' );
    }

    $eh->format( -name  => 'helpBold',
                 -bold  => 1,
                 -color => 'navy',
                 -align => 'right' );

    $eh->format( -name  => 'helpCol',
                 -bold  => 1,
                 -color => 'navy',
                 -align => 'top' );

    $eh->format( -name  => 'note',
                 -italic  => 1,
                 -color => 'brown' );

    $eh->format( -name     => 'section',
                 -bold     => 1,
                 -bg_color => 'yellow',
                 -color    => 'navy',
                 -align    => 'top' );

    $eh->format( -name     => 'helpCat',
                 -bold     => 1,
                 -italic   => 1,
                 -bg_color => 'yellow',
                 -color    => 'blue',
                 -align    => 'top' );

    $eh->format( -name      => 'url',
                 -underline => 1,
                 -color    => 'blue' );

    $eh->format( -name  => 'helpSheet',
                 -align => 'top' );

    $eh->format( -name  => 'helpDescr',
                 -align => 'vjustify' );

    $eh->format( -name  => 'bold',
                 -bold  => 1 );
    my %colCats;
    foreach my $type ("Genes", "RNAs", "Variants") {
        my $cols = $args->{uc($type."COL")} ||
            $self->choose_cols( %{$args}, -type => $type );
        my (@fmts, @widths);
        foreach my $col (@{$cols}) {
            my $w = 10;
            my $f;
            if (my $cd = $colDescription->{$col}) {
                $w = $cd->{width};
                $f = $cd->{fmt};
                push @{$cd->{used}}, $type;
                push @{$colCats{$cd->{cat} || 'Unclassified'}}, $col;
            } else {
                my $cd = $colDescription->{"ZZZ $col"} ||= {
                    desc => "NEED TO ADD DESCRIPTION",
                };
                push @{$cd->{used}}, $type;
            }
            push @widths, $w;
            push @fmts, $f;
        }
        my $sheet = $eh->sheet( -name    => $type,
                                -freeze  => 1,
                                -width   => \@widths,
                                # -formats => \@fmts,
                                -columns => $cols, );
        $sheet->{SR_COLS} = $cols;
    }

    my $secFmt    = ['section','section','section','section'];
    $eh->sheet( -name    => 'Help',
                -width   => [15,15,15,15,15],);

    $eh->add_row('Help', ["Please see the 'Column Help' worksheet for details about specific column headers"], ['note']);

    $eh->blank_row('Help');
    $eh->add_row('Help', ["Impact Classes",'','',''], $secFmt);
    $eh->add_row('Help', ["Variants are classified into the following impact 'classes' based on their impact code"], ['note']);
    $eh->add_row('Help', ['High', join(' ', $self->high_impacts())],
                 ['HighImpact']);
    $eh->add_row('Help', ['Medium', join(' ', $self->med_impacts())],
                 ['MedImpact']);
    $eh->add_row('Help', ['Low', join(' ', $self->low_impacts())],
                 ['LowImpact']);
    $eh->add_row('Help', ['Other', 'All other impact codes not iterated above'],
                 ['UnkImpact']);
    $eh->add_row('Help', ['Multiple', 'If more than one impact code is reported this coloration will be used'],
                 ['MultImpact']);
    $eh->add_row('Help', ["See the 'Column Help' worksheet for detailed information on each code."], ['note']);

    my $metrics = {
        ComPerKB   => [    50,    20,     8 ],
        Common     => [   125,    60,    20 ],
        HotSpot3   => [   270,   120,    60 ],
        HotSpot6   => [    40,    20,    10 ],
        STaaImp    => [     6,     4,     2 ],
        STper100aa => [    40,    25,    18 ],
        STperKB    => [   140,   100,    65 ],
        STperPoly  => [   120,    70,    40 ],
        TCGA       => [   250,   125,    70 ],
    };
    $eh->format( -name     => 'MetHi',
                 -bold     => 1,
                 -bg_color => 'red',
                 -align    => 'center' );
    $eh->format( -name     => 'MetMed',
                 -bold     => 1,
                 -bg_color => 'yellow',
                 -align    => 'center' );
    $eh->format( -name     => 'MetLo',
                 -bold     => 1,
                 -bg_color => 'green',
                 -align    => 'center' );

    $eh->blank_row('Help');
    $eh->add_row('Help', ["Cancer Metrics Coloration",'','',''], $secFmt);
    $eh->add_row('Help', ["The following color schemes are designed to highlight potentially 'interesting' values"], ['note']);
    foreach my $met (sort keys %{$metrics}) {
        my $scale = $metrics->{$met};
        if ($#{$scale} !=2) {
            $args->msg("[?]","Metric $met does not have ranges provided");
            next;
        }
        $colDescription->{$met}{cb} = sub {
            my $self = shift;
            my $v    = shift;
            my $f    = "AlnCen";
            if ($v) {
                if ($v >= $scale->[0]) {
                    $f = 'MetHi';
                } elsif ($v >= $scale->[1]) {
                    $f = 'MetMed';
                } elsif ($v >= $scale->[2]) {
                    $f = 'MetLo';
                }
            }
            return ($f, $v);
        };
        $eh->add_row('Help', [$met, map { ">$_" } @{$scale}], 
                     ['helpBold','MetHi','MetMed','MetLo']);
    }
    $eh->add_row('Help', ["They have been chosen somewhat arbitrarily from scanning actual values reported for Genes"], ['note']);
    $eh->add_row('Help', ["Roughly, the High class corresponds to the top 20 values, Medium the top 100, and Low the top 500"], ['note']);
    $eh->add_row('Help', ["There are roughly 750 genes without reported polymorphism, but will still have an STperPoly metric due to a denominator pseudocount of 1"], ['note']);
    $eh->add_row('Help', ["Suggestions for changing the thresholds (or for new metrics) are welcome"], ['note']);

    $eh->blank_row('Help');

    $eh->sheet( -name    => 'Column Help',
                -freeze  => 1,
                -width   => [15,20,70],
                # -formats => ['helpCol'],
                -columns => ['Column Name', 'Sheets', 'Description'], );
    my @cats = sort {uc($a) cmp uc($b)} keys %colCats;
    foreach my $cat (@cats) {
        $eh->add_row('Column Help', [$cat, '', ''],
                     ['helpCat','helpCat','helpCat']);
        foreach my $col (sort @{$colCats{$cat}}) {
            if (my $cd = $colDescription->{$col}) {
                if (my $u = $cd->{used}) {
                    $eh->add_row('Column Help', [$col, join(',', @{$u}), $cd->{desc}], ['helpCol', 'helpSheet', 'helpDescr'] );
                }
            }
        }
        $eh->blank_row('Column Help');
    }
    map { delete $_->{used} } values %{$colDescription};
    return $eh;
}

sub add_excel_rows {
    my $self = shift;
    my ($rows, $sname, $eh) = @_;
    my $sheet = $eh->sheet($sname);
    my $cols  = $sheet->{SR_COLS};
    foreach my $rh (@{$rows}) {
        my (@row, @fmt);
        for my $i (0..$#{$cols}) {
            my $col = $cols->[$i];
            my $cd  = $colDescription->{$col} || {};
            my $v   = $rh->{$col};
            if (my $r = ref($v)) {
                if ($r eq 'ARRAY') {
                    $v = join($joiner, @{$v});
                } elsif ($r eq 'HASH') {
                    $v = join($joiner, sort keys %{$v});
                }
            }
            my $f   = $cd->{fmt};
            if (my $cb = $cd->{cb}) {
                ($f, $v) = &{$cb}( $self, $v );
            }
            push @fmt, $f;
            push @row, $v;
        }
        $eh->add_row_explicit($sheet, \@row, \@fmt);
    }
}

*add_locus_excel = \&add_gene_excel;
*add_loci_excel = \&add_gene_excel;
sub add_gene_excel {
    my $self = shift;
    my ($rows, $eh) = @_;
    return $self->add_excel_rows($rows, 'Genes', $eh);
}
sub add_rna_excel {
    my $self = shift;
    my ($rows, $eh) = @_;
    return $self->add_excel_rows($rows, 'RNAs', $eh);
}
sub add_variant_excel {
    my $self = shift;
    my ($rows, $eh) = @_;
    my %cbStore;
    foreach my $col (qw(LL)) {
        $cbStore{$col} = $colDescription->{$col}{cb};
        delete $colDescription->{$col}{cb};
    }
    my $rv = $self->add_excel_rows($rows, 'Variants', $eh);
    while (my ($col, $cb) = each %cbStore) {
        $colDescription->{$col}{cb} = $cbStore{$col};
    }
    return $rv;
}


1;
