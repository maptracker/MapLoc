# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::SnpTracker::MapLoc::GeneAlignment;
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
use BMS::SnpTracker::MapLoc::TransientAlignment;
use BMS::SnpTracker::MapLoc::RnaToGeneAlignment;

use vars qw(@ISA);
@ISA      = qw(BMS::SnpTracker::MapLoc::TransientAlignment);

# These objects are created by ::Gene::exon_footprints()
# Each is a genome-anchored alignment representing the union of exonic
# regions for one or more RNAs. These locations could represent complex
# interpretations if a single gene is assigned to multiple locations in the
# genome

sub set_rnas {
    my $self = shift;
    if ($self->{RNA_ALIGN}) {
        $self->err("RNA alignment data already set");
    } elsif (my $hsps = shift) {
        $self->{RNA_ALIGN} = $hsps;
    } else {
        $self->err("No RNA information provided");
    }
}

sub rna_genomic_alignments {
    return map { $_->[0] } @{shift->{RNA_ALIGN} || []};
}

sub each_rna {
    return map { $_->[1] } @{shift->{RNA_ALIGN} || []};
}

sub chrAccId {
    my $self = shift;
    unless (defined $self->{GA_CHR_ACC}) {
        #my @csi = $self->chr_seq_info();
        #my ($gAccId, $cind, $chr, $build) = @{$csi[0] || []};
        #$self->{GA_CHR_ACC} = $gAccId || 0;
        # Just require a pre-set order of [gene,chr]: 
        my @sids = $self->seq_ids();
        $self->{GA_CHR_ACC} = $sids[1];
    }
    return $self->{GA_CHR_ACC};
}

sub geneAccId {
    my $self = shift;
    unless (defined $self->{GA_GENE_ACC}) {
        my @sids = $self->seq_ids();
        $self->{GA_GENE_ACC} = $sids[0];
    }
    return $self->{GA_GENE_ACC};
}

sub rna_alignments {
    my $self = shift;
    unless ($self->{RNA_ALIGN_OBJS}) {
        # Generate coordinate alignment (intersection) between Gene and RNA
        # Because the gene was built quickly by collating genomic coordinates
        # we need to tease this out now
        my $arr = $self->{RNA_ALIGN_OBJS} = [];
        my $chrid = $self->chrAccId();
        return () unless ($chrid);
        my $geneid = $self->geneAccId();
        foreach my $rdat (sort { $a->[1]->sortable_value() cmp 
                                     $b->[1]->sortable_value() }
                          @{$self->{RNA_ALIGN} || []}) {
            my ($aln, $rna) = @{$rdat};
            my ($rcx) = $rna->cx_genome_part( $aln->pkey );
            my $int = $self->intersect( -hsp => $aln, -anchor => $chrid );
            $int->anchor_id( $geneid );
            bless($int, "BMS::SnpTracker::MapLoc::RnaToGeneAlignment");
            $int->{RGA_GENE_IDX} = $int->seqindex($geneid) * 2;
            $int->{RGA_RNA_IDX}  = $int->seqindex($rna->acc_id) * 2;
            $int->{RGA_RNA_ALN}  = $aln;
            $int->{RGA_RNA_OBJ}  = $rna;
            push @{$arr}, $int;
        }
    }
    return @{$self->{RNA_ALIGN_OBJS}};
}
