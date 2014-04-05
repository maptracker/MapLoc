# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::SnpTracker::MapLoc::Gene;
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
use BMS::SnpTracker::MapLoc::Text;
use BMS::SnpTracker::MapLoc::HasAlignments;
use BMS::SnpTracker::MapLoc::TransientAlignment;
use BMS::SnpTracker::MapLoc::GeneAlignment;

use vars qw(@ISA);
@ISA      = qw(BMS::SnpTracker::MapLoc::Text
               BMS::SnpTracker::MapLoc::HasAlignments );

sub obj_type       { return "Gene"; }
sub each_gene      { return (shift); }
sub each_gene_name { return (shift->name()) }
*acc = \&accession;
sub accession      { return shift->name(); }

# ::Gene
sub each_rna_name {
    # This is the most direct connection between genes and their transcripts
    my $self = shift;
    unless ($self->{RNA_NAMES}) {
        my @rnas = $self->tag_values('Has RNA');
        $self->{RNA_NAMES} = \@rnas;
    }
    return @{$self->{RNA_NAMES}};
}

# ::Gene - inherit aggregated alignments from RNA
sub read_alignments {
    my $self = shift;
    unless ($self->{ALIGNS}) {
        my $targ = $self->{ALIGNS} = {};
        foreach my $rna ($self->each_rna()) {
            my $rAln = $rna->read_alignments();
            while (my ($aid, $aln) = each %{$rAln}) {
                $targ->{$aid} ||= $aln;
            }
        }
    }
    return $self->{ALIGNS};
}

sub _best_child_exons {
    my $self = shift;
    my $rv = $self->{BEST_EXONS};
    unless ($rv) {
        $rv = $self->{BEST_EXONS} = {};
        foreach my $rna ($self->each_rna()) {
            # 0 = Howbad - ONLY consider the best alignments for each build
            foreach my $aln ($rna->each_alignment( 0 )) {
                push @{$rv->{$aln->build() || ""}}, [$aln, $rna ];
            }
        }
    }
    if (my $build = shift) {
        return $rv->{$self->normalized_build($build)} || [];
    }
    return $rv;
}

sub _pick_build {
    my $self = shift;
    my $args = $self->parseparams( @_ );
    my $build;
    if (my $bReq = $args->{BUILD}) {
        # Explicitly provided
        $build = $self->normalized_build( $bReq );
        unless ($build) {
            $self->msg_once("Can not recover Gene exons for unrecognized build request '$bReq'");
        }
    } else {
        my $bce   = $self->_best_child_exons( );
        my @known = keys %{$bce};
        if ($#known == -1) {
            $self->msg_once("No genomic alignments available for Gene",
                                $self->to_one_line());
        } else {
            $build = $self->best_build( @known);
            unless ($build) {
                $self->msg_once("Failed to find best genome build for Gene",
                                "Recorded builds: ".join(', ', @known),
                                $self->to_one_line());
            }
        }
    }
    return $build;
}

# :: Gene
sub exon_footprints {
    my $self = shift;
    my $args  = $self->parseparams( @_ );
    my $rv;
    if (my $aReq = $args->{ALN} || $args->{ALNS} || $args->{ALIGN}) {
        $self->death("No support yet for custom Gene locations");
    } else {
        my $build = $self->_pick_build( @_ );
        return () unless ($build);
        $rv = $self->{GENE_EXON_HSPS}{$build};
        unless ($rv) {
            $rv = $self->{GENE_EXON_HSPS}{$build} = [];
            my $bce  = $self->_best_child_exons( );
            my $data = $bce->{$build};
            unless ($data) {
                $self->msg_once("No genomic locations for requested build '$build'",
                                $self->to_one_line());
                return @{$rv};
            }
            my $gName    = $self->name();
            my $gId      = $self->pkey();
            my @clusters = $self->_cluster_genomic_hsps( $data );
            my $ml = $self->maploc();
            foreach my $clust (@clusters) {
                my ($chrId, $chr, $build, $str) = map { $clust->{$_} }
                qw(accid chr build str);
                my $chrName = "$chr.$build";
                my @rDat    = @{$clust->{hsps}};
                my @allHsp;
                foreach my $dat (@rDat) {
                    my ($aln, $rna) = @{$dat};
                    my ($hsps)      = $aln->chr_hsps_for_alignment();
                    push @allHsp, @{$hsps};
                }
                my $foot = $ml->add_hsps( \@allHsp );
                if ($str < 0) {
                    $foot = [sort { $b->[1] <=> $a->[1] || 
                                        $b->[0] <=> $a->[0] } @{$foot}];
                } else {
                    $foot = [ sort { $a->[0] <=> $b->[0] ||
                                         $a->[1] <=> $b->[1] } @{$foot} ];
                }
                # die $self->branch($foot);
                my @gHsps;
                my $lastRight = 1;
                foreach my $ghsp (@{$foot}) {
                    my ($gl, $gr) = @{$ghsp};
                    my $l = $lastRight - 1;
                    my $r = $lastRight = $l + $gr - $gl;
                    push @gHsps, [$l, $r, $gl, $gr];
                }
                my $gAln = BMS::SnpTracker::MapLoc::TransientAlignment->new
                    ( $ml );
                $gAln->strand( $str );
                # Order is important: will expect [gene, chr] data:
                $gAln->hsps( \@gHsps );
                $gAln->seq_ids( [$gId, $chrId] );
                $gAln->_seq_lookup();
                $gAln->set_anchor( 1 );
                # $gAln->{RNA_ALIGN} = $clust->{hsps};
                bless ($gAln, 'BMS::SnpTracker::MapLoc::GeneAlignment');
                $gAln->set_rnas( $clust->{hsps} );
                push @{$rv}, $gAln;
            }
        }
    }
    return @{$rv};
}

sub _cluster_genomic_hsps {
    my $self = shift;
    my $data = shift;
    my %bldChr;
    foreach my $dat (@{$data}) {
        my ($aln, $rna) = @{$dat};
        my @csi = $aln->chr_seq_info();
        next if ($#csi != 0);
        my ($accid) = @{$csi[0]};
        my ($s, $e) = $aln->chr_start_end();
        next unless ($s);
        push @{$bldChr{ $accid }{ $aln->str() }}, [ $aln, $rna, $s, $e];
    }
    my @rv;
    while (my ($accid, $strH) = each %bldChr) {
        while (my ($str, $datA) = each %{$strH}) {
            my @dats = sort { $a->[2] <=> $b->[2] || 
                                  $a->[3] <=> $b->[3] } @{$datA};
            my $seed = shift @dats;
            push @rv, {
                accid => $accid,
                str   => $str,
                hsps  => [ [$seed->[0], $seed->[1]] ],
                s     => $seed->[2],
                e     => $seed->[3],
            };
            foreach my $dat (@dats) {
                my ($aln, $rna, $s, $e) = @{$dat};
                my $prior = $rv[-1];
                if ($s <= $prior->{e}) {
                    # Overlaps with last cluster
                    push @{$prior->{hsps}}, [$aln, $rna];
                    $prior->{e} = $e if ($prior->{e} < $e);
                } else {
                    # Add new cluster
                    push @rv, {
                        accid => $accid,
                        str   => $str,
                        hsps  => [ [$aln, $rna] ],
                        s     => $s,
                        e     => $e,
                    };
                }
            }
        }
    }
    # Add in Build and Chr
    foreach my $dat (@rv) {
        my ($aln) = $dat->{hsps}[0][0];
        $dat->{build} = $aln->build();
        $dat->{chr}   = $aln->chr();
    }
    # die $self->branch(\@rv);
    return @rv;
}

# ::Gene
sub cx_exon_parts {
    my $self = shift;
    $self->bench_start();
    my $args   = $self->parseparams( @_, -howbad => 0 );
    my $ml     = $self->maploc();
    my $tracks = $args->{TRACKS} || {};
    # $gAln is the projected alignment of the Gene onto a single genomic
    # location, as represented by the union of exonic components. So the
    # coordinate space will be genomic in one case, and abstracted 'gene'
    # coordinates covering one or more transcripts:
    my $gAln   = $args->{ALN};
    return $tracks unless ($gAln);
    my $pkey  = $self->pkey();
    my $gind  = $gAln->seqindex($pkey);
    my @csi   = $gAln->chr_seq_info();
    return $tracks unless ($#csi == 0);
    my ($gAccId, $cind, $chr, $build) = @{$csi[0]};
    my $gind2 = $gind * 2;
    my $cind2 = $cind * 2;
    # warn "($gAccId, $cind, $chr, $build) $gind2 $cind2";
    my $str   = $gAln->strand();
    my $fl    = $chr || "Unk";
    $fl      .= ".$build" if ($build);
    my @hdat;
    my $hsps = $gAln->hsps();
    my @coords = sort { $a->[$gind2] <=> $b->[$gind2] } $gAln->coordinates();
    my @cols   = ('#669999','#996699');
    for my $co (0..$#coords) {
        my $coord = $coords[$co];
        # $self->maploc->prebranch($coord);
        my ($gs, $ge) = ($coord->[$cind2], $coord->[$cind2+1]);
        push @{$tracks->{'Genome'}}, {
            hideName => 1,
            name  => "HSP ".($co+1),
            data  => [[$coord->[$gind2], $coord->[$gind2+1]]],
            fill  => $cols[$co % 2],
            note  => "$fl:$gs..$ge",
            len   => $ge - $gs + 1,
        };
    }
    # Assuming that the alignment was created with exon_footprints()
    # then each one should have the RNA_ALIGN set with the original
    # RNA-to-genome alignments
    my %cFeat;
    foreach my $int ($gAln->rna_alignments()) {
        my $rna   = $int->rna();
        my $acc   = $rna->acc();
        my $rAln  = $int->rna_aln();
        my ($rcx) = $rna->cx_genome_part( $rAln->pkey );
        my $gind2 = $int->gene_index();
        my $rind2 = $int->rna_index();
        my @inds  = ($gind2, $gind2 + 1, $rind2, $rind2 + 1);
        my @data;
        foreach my $coord ($int->coordinates()) {
            push @data, [ map { $coord->[$_] } @inds ];
        }
        $rcx->{data}  = [ sort { $a->[0] <=> $b->[0] } @data ];
        $rcx->{exons} = [ sort { $a->[0] <=> $b->[0] } @{$rcx->{exons}} ];
        $rcx->{dir}   = 'right';
        push @{$tracks->{'RNA'}}, $rcx;
        # Add in the CDS information
        foreach my $cgr ($rna->cds_genomic_ranges( -aln => $rAln )) {
            my $cint = $cgr->intersect( -hsp => $int );
            my $type = $cgr->simple_value('Key');
            my $loct = $cint->location_text_for_seq( $pkey );
            # $ml->preprint("$type = $loct\n" . $cint->full_text() );
            push @{$cFeat{$type}{$loct}}, $acc;
        }
        # $self->maploc->prebranch($rcx);
    }

    while (my ($type, $locH) = each %cFeat) {
        while (my ($loct, $accs) = each %{$locH}) {
            my ($str, $hspTxt) = split(':', $loct);
            my @data = map { [split(/\.\./, $_)]} split(',', $hspTxt);
            map { $_->[0]++; $_->[1]--; } @data;
            my $cx = $type eq 'ATG' ? {
                type    => 'InitMET',
                fill    => '#33ff33',
                outline => '#33ff33',
            } : $type eq 'Stop' ? {
                type    => 'Stop',
                fill    => '#ff0000',
                outline => '#ff0000',
            } : {
                type    => 'Unknown',
                fill    => '#dddddd',
                outline => '#dddddd',
            };
            $cx->{dir}      = $str < 0 ? 'left' : 'right';
            $cx->{hideName} = 1;
            $cx->{id} = $type;
            $cx->{data} = [ sort { $a->[0] <=> $b->[0] } @data ];
            $cx->{tags}{"Source"} = $accs;
            push @{$tracks->{'Genome'}}, $cx;
        }
    }
    $self->bench_end();
    return $tracks;
}

# ::Gene -  we dynamically re-cluster alignments based on a 
# user-provided howBad value
sub each_alignment {
    my $self = shift;
    my $hb    = shift;
    my $hbkey = defined $hb ? $hb : "";
    if (my $cga = $self->{CACHED_GENE_ALN}{$hbkey}) {
        return @{$cga};
    }
    $self->bench_start();
    my @alns  = $self->SUPER::each_alignment( $hb );
    my %clustered;
    foreach my $aln (@alns) {
        my @ci  = $aln->chr_seq_info();
        next unless ($#ci == 0);
        my ($chrID, $chrInd, $chr, $build) = @{$ci[0]};
        my $str     = $aln->strand();
        my @seqs    = $aln->seqs();
        my $chrSeq  = $seqs[$chrInd];
        my $targ    = $clustered{$chrSeq}{$str} ||= [];
        my $hsps    = $aln->hsps_for_seq( $chrInd );
        # Sort the individual HSPs by ascending coordinate
        my @shsp    = sort { $a->[0] <=> $b->[0] || 
                                 $a->[1] <=> $b->[1] } @{$hsps};
        push @{$targ}, [ \@shsp, [$aln] ]
    }
    my $ml = $self->maploc();
    my @rv;
    while (my ($cName, $sH) = each %clustered) {
        while (my ($str, $hspSet) = each %{$sH}) {
            # We are now considering the set of HSPs from the same genomic
            # region (chromosome) and on the same strand. Order them
            # by coordinate
            my @alnDat = sort { $a->[0][0][0] <=> $b->[0][0][0] || 
                                    $a->[0][-1][1] <=> $b->[0][-1][1] } @{$hspSet};
            # Seed the cluster array
            my @clusts = ( shift @alnDat );
            foreach my $newHsp (@alnDat) {
                my $priHsp = $clusts[$#clusts];
                if ($newHsp->[0][0][0] + 2 <= $priHsp->[0][-1][1]) {
                    # The new HSP overlaps with the previous one
                    # We are NOT counting adjancency! (+1)
                    # Add the current alignment object to the prior data:
                    push @{$priHsp->[1]}, @{$newHsp->[1]};
                    # Stuff together all the HSPs, sort, and collate:
                    my @both = sort {$a->[0] <=> $b->[0] || 
                                         $a->[1] <=> $b->[1] }
                    ( @{$priHsp->[0]}, @{$newHsp->[0]} );
                    my @union = (shift @both);
                    foreach my $hsp (@both) {
                        if ($hsp->[0] + 1 <= $union[-1][1]) {
                            # This HSP overlaps OR is adjacent
                            # Just extend the prior HSP
                            $union[-1][1] = $hsp->[1] if
                                ($union[-1][1] < $hsp->[1]);
                        } else {
                            # No overlap - add as new
                            push @union, $hsp;
                        }
                    }
                    $priHsp->[0] = \@union;
                } else {
                    # The new HSP is "further down" the genome than the last
                    # add it as a new cluster
                    push @clusts, $newHsp;
                }
            }
            # Now turn the clusters into actual alignment objects
            foreach my $clust (@clusts) {
                my @chrHsps = @{$clust->[0]};
                @chrHsps = sort { $b->[1] <=> $a->[1] || $b->[0] <=> $a->[0] }
                @chrHsps if ($str < 0);
                my $lastRight = 1;
                my @hspPairs;
                foreach my $cHsp (@chrHsps) {
                    # Mind the math! We are dealing with FLANKS, not start/end
                    my $l = $lastRight - 1;
                    $lastRight += $cHsp->[1] - $cHsp->[0] - 1;
                    push @hspPairs, [ $l, $lastRight, @{$cHsp} ];
                }
                my ($sc, $hb, %src);
                foreach my $aln (@{$clust->[1]}) {
                    # For score and howbad, take the worst of all worlds:
                    my $asc = $aln->score();
                    $sc = $asc if (!defined $sc || 
                                   (defined $asc && $asc < $sc));
                    my $ahb = $aln->howbad();
                    $hb = $ahb if (!defined $hb ||
                                   (defined $ahb && $ahb > $hb));
                    if (my $auth = $aln->source()) {
                        $src{$auth} = 1;
                    }
                }
                my $locAln = BMS::SnpTracker::MapLoc::TransientAlignment
                    ->new($ml);
                $locAln->seqs( $self->name(), $cName );
                $locAln->hsps( \@hspPairs );
                $locAln->strand( [1, $str] );
                $locAln->score( $sc );
                $locAln->howbad( $hb );
                $locAln->source( join(', ', sort keys %src) );
                push @rv, $locAln;
            }
        }
    }
    # print "<pre>".$self->branch(\%clustered)."</pre>\n";
    $self->{CACHED_GENE_ALN}{$hbkey} = \@rv;
    $self->bench_end();
    return @rv;
}

sub write_alignments {
    # Do not do try to write alignments from the Gene object
    my $self = shift;
    $self->err("Attempt to write alignments from Gene object");
}
