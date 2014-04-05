# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::SnpTracker::MapLoc::RNA;
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
use BMS::SnpTracker::MapLoc::Featured;
use BMS::SnpTracker::MapLoc::HasAlignments;
use BMS::SnpTracker::MapLoc::Tagged;
use BMS::SnpTracker::MapLoc::Range;
use BMS::SnpTracker::MapLoc::TaggedRange;

use vars qw(@ISA);
@ISA      = qw(BMS::SnpTracker::MapLoc::Featured
               BMS::SnpTracker::MapLoc::HasAlignments
               BMS::SnpTracker::MapLoc::Tagged);

# http://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly
our $refSeqTypeMap = {
    AC => 'Complete genomic molecule',
    NC => 'Complete reference genomic molecule',
    NG => 'Incomplete genomic region',
    NT => 'Contig or scaffold',
    NW => 'WGS contig or scaffold',
    NS => 'Environmental sequence',
    NZ => 'Unfinished WGS',
    NM => 'Reference mRNA',
    NR => 'Reference ncRNA',
    XM => 'Predicted model mRNA',
    XR => 'Predicted model ncRNA',
    AP => 'Protein, alternate assembly',
    NP => 'Reference protein',
    YP => 'Predicted protein',
    XP => 'Predicted model protein',
    ZP => 'Predicted WGS protein',
};

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my ($ml, $acc, $src) = @_;
    my $self = {
        MAPLOC => $ml,
        ALIGNS => {},
    };
    bless ($self, $class);
    $self->bench_start();
    if (!$acc) {
        $self->err("Attempt to create RNA object without ID");
        $self->bench_end();
        return undef;
    } elsif ($acc =~ /^\d+$/) {
        # User is passing an ID directly
        unless ($src) {
            # A single ID indicates an RNA ID
            $self->{PKEY} = $acc;
            my $get = $ml->{STH}{ACC_SRC_FOR_RNA_ID} ||= $ml->dbh->prepare
                ( -name => "Get data with rna_id",
                  -sql  => "SELECT acc_id, src_id, seq_id, cds_start, cds_end, anomaly FROM rna WHERE rna_id = ?" );
            $get->execute($acc);
            my $rows = $get->fetchall_arrayref();
            if ($#{$rows} == 0 && $rows->[0][0]) {
                my @v = @{$rows->[0]};
                $self->{ACC_ID} = $v[0];
                $self->{SRC_ID} = $v[1];
                $self->{SEQ_ID} = $v[2];
                $self->{CDS} = [ $v[3], $v[4] ];
                $self->{ANOMALY} = $v[5];
                $self->bench_end();
                return $self;
            } else {
                $self->err("Failed to recover RNA rna_id = $acc");
                $self->bench_end();
                return undef;
            }
        }
        # Otherwise, this is an accession ID
        $self->{ACC_ID} = $acc;
        ($acc, $src) = ();
    } else {
        $self->{ACC} = $acc;
    }
    if (!$src) {
        $self->err("Attempt to create RNA object without source (acc = '$acc')");
        $self->bench_end();
        return undef;
    } elsif ($src =~ /^\d+$/) {
        # User is passing src_id directly
        $self->{SRC_ID} = $src;
    } else {
        $self->{SRC} = $src;
    }
    # At this point we should be assured of having ACC and SRC defined
    $self->bench_end();
    return $self;
}

sub obj_type { return "RNA"; }

# ::RNA
sub each_rna      { return (shift); }
sub each_rna_name { return (shift->name()); }

sub variant_number {
    my $self = shift;
    unless (defined $self->{VAR_NUM}) {
        $self->read();
        my %nums;
        foreach my $desc ($self->tag_values('Description')) {
            if ($desc =~ /transcript variant (\d+)/) {
                $nums{$1} = 1;
            }
        }
        my @n = sort { $b <=> $a } keys %nums;
        $self->{VAR_NUM} = $n[0] || 0;
    }
    return $self->{VAR_NUM};
}

sub is_pseudo {
    my $self = shift;
    unless (defined $self->{IS_PSEUDO}) {
        $self->read();
        my @descs = $self->tag_values('Description');
        my $dtxt  = join(' ', @descs) || "";
        $self->{IS_PSEUDO} = $dtxt =~ /pseudo ?gene/i ? 1 : 0;
    }
    return $self->{IS_PSEUDO};
}

sub sortable_value {
    my $self = shift;
    unless ($self->{SORT_VAL}) {
        $self->read();
        my @bits = ();
        push @bits, $self->is_pseudo();
        my ($sym) = $self->all_symbols();
        push @bits, sprintf("%6s", $sym || "");
        push @bits, sprintf("%03d", $self->variant_number() || 999);
        my $acc = $self->accession();
        if ($acc =~ /^(ENS[A-Z]+|[NX][MRP])_?(\d+)/) {
            my ($prfx, $num) = ($1,$2);
            push @bits, sprintf("%6s", $prfx);
            push @bits, sprintf("%012d", $num);
        }
        push @bits, $acc;
        $self->{SORT_VAL} = join(' ', @bits);
    }
    return $self->{SORT_VAL};
}

sub nice_name {
    my $self = shift;
    unless ($self->{NICE_NAME}) {
        $self->read();
        my $acc   = $self->accession();
        my @bits  = ($acc);
        push @bits, "[PSEUDOGENE]" if $self->is_pseudo();
        if (my $var =  $self->variant_number()) {
            push @bits, "variant $var";
        }
        my ($sym) = $self->tag_values('Symbol');
        if ($sym) {
            # Use the symbol as the first identifier
            $bits[0] = $sym;
            push @bits, "($acc)";
        }
        $self->{NICE_NAME} = join(' ', @bits);
    }
    return $self->{NICE_NAME};
    
}

# ::RNA
sub to_text {
    my $self = shift;
    $self->bench_start();
    my $pre  = shift || "";
    my $txt  = $pre.$self->to_one_line()."\n";
    $txt    .= $self->feature_text($pre."  ");
    $txt    .= $self->tag_text($pre."  ");
    $self->bench_end();
    return $txt;
}

sub to_xml {
    my $self = shift;
    $self->bench_start();
    my $pad = shift || "";
    my $xml = $pad."<rna>\n";
    $xml .= $pad."</rna>\n";
    $self->bench_end();
    return $xml;
}

*fasta = \&to_fasta;
sub to_fasta {
    my $self = shift;
    $self->bench_start();
    my $args  = $self->parseparams( @_ );
    my $fa    = ">".$self->acc();
    if (my $desc = $self->desc()) {
        $fa .= " $desc";
    }
    my $seq   = $args->{VARIANT} ? $self->variant_seq( @_ ) : $self->seq();
    $fa .= "\n";
    my $block = $args->{BLOCK} || 80;
    my $slen  = CORE::length($seq);
    for (my $i = 0; $i < $slen; $i += $block) {
        $fa .= substr($seq, $i, $block)."\n";
    }
    $self->bench_end();
    return $fa;
}

# ::RNA
sub to_one_line {
    my $self = shift;
    $self->bench_start();
    my $len = $self->length();
    if ($len) {
        my ($s, $e) = $self->cds();
        if ($s) {
            $len .= " [$s,$e]";
            if (my $an = $self->anomaly()) {
                $len .= " !ANOMALY! $an";
            }
        }
    } else {
        $len = "?";
    }
    my $txt  = sprintf("%s [%s] %s Len:%s", $self->accession(),
                       $self->handle(), $self->source(), $len);
    $self->bench_end();
    return $txt;
}

# ::RNA
sub cx_exon_part {
    my $self = shift;
    $self->bench_start();
    my $args   = $self->parseparams( @_ );
    my $ml     = $self->maploc();
    my $aid    = $self->acc_id();
    my $acc    = $self->acc();
    my @alns   = $self->alignments_from_parameters( $args );
    my %exons;
    my @rv;
    foreach my $aln (@alns) {
        $aln->set_anchor( $aid );
        my @oseqs = $aln->other_seqs( $acc );
        next unless ($#oseqs == 0);
        my ($chr, $build, $taxa, $type) = $ml->standardize_chromosome
            ( $oseqs[0] );
        # warn "$chr, $build, $taxa, $type)";
        my ($Acrd, $inds, $strs) = $aln->ordered_coordinates( );
        my $sc  = $aln->score();
        my $src = $aln->source();
        $src    =~ s/ Unmasked//;
        my $hb  = $aln->howbad();
        my $str = $aln->strand();
        $str = !$str ? '?' : $str < 0 ? '-1' : '+1';
        for my $h (0..$#{$Acrd}) {
            my ($rs, $re, $gs, $ge) = @{$Acrd->[$h]};
            my $num = $h + 1;
            my $ex = $exons{"$rs\t$re"} ||= {
                s => $rs,
                e => $re,
                chr => {},
                num => {},
                type => 'Exon',
            };
            $ex->{num}{$num} = 1;
            push @{$ex->{chr}{$chr}}, [ $build, $gs, $ge, $num, $str ];
            $ex->{src}{$src} = 1;
            $ex->{sc}{$sc} = 1;
            $ex->{hb}{$hb} = 1;

            next unless ($h);
            my ($ps, $pe, $pgs, $pge) = @{$Acrd->[$h-1]};
            my ($is, $ie) = ($rs, $rs-1);
            my $int = $exons{"$is\t$ie"} ||= {
                s => $is,
                e => $ie,
                type => 'Intron',
                chr => {},
                num => {},
            };
            $int->{num}{$h} = 1;
            my ($igs, $ige) = ($pge+1, $gs-1);
            ($igs, $ige) = ($ige, $igs) if ($igs > $ige);
            push @{$int->{chr}{$chr}}, [ $build, $igs, $ige, $h, $str ];
            $int->{src}{$src} = 1;
            $int->{sc}{$sc} = 1;
            $int->{hb}{$hb} = 1;
        }
        # push @rv, $aln->cx_genome_part();
    }
    my @genes = $self->each_gene_name();
    my $gacc  = $#genes == 0 ? $genes[0] : "";
    my $seq = $self->seq();
    my ($cs, $ce) = $self->cds();
    foreach my $ex (sort { $a->{s} <=> $b->{s} || $a->{e} <=> $b->{e} }
                     values %exons) {
        my $type = $ex->{type};
        my @nums = sort { $a <=> $b } keys %{$ex->{num}};
        my $id   = sprintf("%s%s %s", $type, $#nums == 0 ? '' : 's', join(',', @nums));
        my ($s, $e) = ($ex->{s} || 0, $ex->{e} || 0);
        my $len = $e - $s + 1;
        my ($hb) = sort { $a <=> $b } keys %{$ex->{hb}};
        my ($sc) = sort { $b <=> $a } keys %{$ex->{sc}};
        my $hbc  = $hb < 1 ? "#006600" : $hb < 3 ? "#ff9900" : "#ff0000";
        my $cx = {
            id       => $id,
            type     => $type,
            gene     => $gacc,
            data     => [ $s > $e ? [$e, $s] : [$s, $e] ],
            dir      => 'right',
            fill     => '#99ccff',
            outline  => $hbc, # '#663300',
            source   => join(' / ', sort keys %{$ex->{src}}),
            score    => $sc,
            howbad   => $hb,
            coords   => $ex->{chr},
        };
        push @rv, $cx;
        if ($type eq 'Intron') {
            delete $cx->{dir};
            $cx->{fill} = '#cc9933';
            $cx->{outline} = '#cc9966';
            $cx->{insertion} = 1;
            next;
        }
        my ($ss, $se) = ($s, $e);
        my $cds;
        if ($cs && $cs <= $e && $ce >= $s) {
            # The CDS overlaps here
            # Get the CDS coordinates that are just contained here:
            my ($lcs, $lce) = ($cs, $ce);
            $lcs = $s if ($lcs < $s);
            # Find the phase at the start of the exon:
            my ($ps, $psp) = $self->rna_to_protein_coordinate( $lcs );
            # Calculate the local (within this exon) CDS coordinates:
            my $locCDSstart = $lcs - $s + 1 - $psp;
            $cx->{data}[0][2] = 1;
            if ($locCDSstart < 1) {
                # After correcting for phase we are translating just prior to
                # start of exon. We need to get some extra sequence
                my $xtra   = 1 - $locCDSstart;
                $ss       -= $xtra;
                $locCDSstart = 1;
                $cx->{data}[0][2] = $xtra + 1;
            }
            $cx->{data}[0][3] = $cx->{data}[0][2] + $len - 1;
            # Now check the 3-prime end of the exon CDS:
            $lce = $e if ($lce > $e);
            my ($pe, $pse) = $self->rna_to_protein_coordinate( $lce );
            my $mod3 = ($lce - $lcs + 1) % 3;
            if ($mod3) {
                # We need a few more bases to finish the last exon
                $lce += $mod3;
                $se = $lce if ($lce > $se);
            }
            my $locCDSend = $lce - $ss + 1;
            $cx->{aacoord} = [ $ps, $pe ];
            $cx->{phase}  = $psp;
            $cx->{cds}    = [ $locCDSstart, $locCDSend ];
        }
        $cx->{sequence} = substr($seq, $ss - 1, $se - $ss + 1);
    }
    @rv = sort { $a->{howbad} <=> $b->{howbad} } @rv;
    # warn join(" + ", map { $_->{howbad} } @rv);
    if ($cs) {
        unshift @rv, ( {
            id   => 'ATG',
            type => 'InitMET',
            data => [ [$cs, $cs+2] ],
            fill => '#33ff33',
            outline => '#33ff33',
        } , {
            id   => 'Stop',
            type => 'Stop',
            data => [ [$ce-2, $ce] ],
            fill => '#ff0000',
            outline => '#ff0000',
        });
    }
    $self->bench_end();
    return @rv;
}


# ::RNA
sub cx_genome_part {
    my $self = shift;
    $self->bench_start();
    my $alnReq = shift;
    my @rv;
    my $tags  = $self->all_tag_values();
    my $seq   = $self->seq();
    my $acc   = $self->acc();
    my @genes = $self->each_gene_name();
    my $gacc  = $#genes == 0 ? $genes[0] : "";
    my $ml    = $self->maploc();
    my ($cs, $ce) = $self->cds();
    my $cds   = ($cs && $ce) ? [$cs, $ce] : undef;
    my $llid  = $tags && exists $tags->{LocusID} && $#{$tags->{LocusID}} == 0 ?
        $tags->{LocusID}[0] : "";
    $llid =~ s/^LOC//;
    foreach my $alnId ($self->each_alignment_id()) {
        next if ($alnReq && $alnId != $alnReq);
        my $aln = $ml->get_alignment($alnId);
        unless (defined $aln->anchor_index) {
            my @oids = $aln->other_seq_ids( $self->acc_id );
            if ($#oids == 0) {
                $aln->set_anchor( $oids[0] );
            } else {
                next;
            }
        }
        my @hashes = $aln->cx_genome_part();
        # $ml->preprint($self->name()." ".$aln->strand()." ".$hashes[0]{dir});
        my $nn = $self->nice_name();
        my $type = $self->comma_tag_value( "Ensembl Biotype" );
        unless ($type) {
            if ($acc =~ /^([ANXYZ][CGTWSMRP])_/) {
                $type = $refSeqTypeMap->{$1} || "";
            }
        }
        my $col  = $ml->user_text_color_or_pastel( $type );
        foreach my $hash (@hashes) {
            $hash->{type}     = $type;
            $hash->{label}    = $nn;
            $hash->{tags}     = $tags;
            $hash->{llid}     = $llid;
            $hash->{fill}     = $col;
            $hash->{gene}     = $gacc;
            $hash->{outline}  = "#FFFFFF";
            $hash->{sequence} = $seq;
            if ($cds) {
                $hash->{cds} = $cds;
                # $hash->{translate} = [1];
            }
            if ($seq) {
                my @mem = @{$hash->{members} || []};
                my %mlu = map { $mem[$_] => $_ } (0..$#mem);
                my $ind = $mlu{$acc};
                
                if (defined $ind) {
                    $ind *= 2;
                    $hash->{exons} = [ map { 
                        [ $_->[$ind], $_->[$ind+1] ] } @{$hash->{data}} ]
                    #my $sd = $hash->{sequences} = [];
                    #foreach my $d (@{$hash->{data}}) {
                    #    my ($s, $e) = ($d->[$ind], $d->[$ind+1]);
                    #    push @{$sd}, substr($seq, $s-1, $e - $s + 1);
                    #}
                }
            }
        }
        push @rv, @hashes;
    }
    $self->bench_end();
    return @rv;
}

# ::RNA
*name = \&accession;
*acc = \&accession;
sub accession {
    return shift->_generic_getSet_name( 'ACC', @_ );
}

sub unversioned_accession {
    my $self = shift;
    unless (defined $self->{UNV_ACC}) {
        $self->{UNV_ACC} = $self->unversion( $self->accession() );
    }
    return $self->{UNV_ACC};
}

sub version {
    my $self = shift;
    unless (defined $self->{VERS_NO}) {
        if (my $acc = $self->accession) {
            if ($acc =~ /\.(\d+)$/) {
                $self->{VERS_NO} = $1;
            } else {
                $self->{VERS_NO} = 0;
            }
        }
    }
    return $self->{VERS_NO};
}

# ::RNA
*acc_id = \&accession_id;
sub accession_id {
    return shift->_generic_get_id( 'ACC', @_ );
}

*src = \&source;
sub source {
    return shift->_generic_getSet_name( 'SRC', @_ );
}

sub src_id {
    return shift->_generic_get_id( 'SRC', @_ );
}

sub seq {
    return shift->_generic_getSet_name( 'SEQ', @_ );
}

# ::RNA
sub json_data {
    my $self = shift;
    $self->bench_start();
    $self->read();
    my ($cs, $ce) = $self->cds();
    my ($gene)    = $self->each_gene_name();
    my $rv = {
        acc     => $self->acc(),
        tags    => $self->all_tag_values(),
        rid     => $self->pkey(),
        gene    => $gene || "",
        name    => $self->nice_name(),
        anomaly => $self->anomaly(),
        source  => $self->source(),
        cds     => [$cs, $ce],
    };
    $self->bench_end();
    return $rv;
}

# ::RNA
sub impact {
    my $self = shift;
    $self->bench_start();
    $self->read();
    my $args = $self->parseparams( @_ );
    my $rAcc = $self->acc();
    my $ml   = $self->maploc();
    my $su   = $ml->sequtils();
    my $rv = {
        imp    => "GEN",
        rid    => $self->pkey(),
        var    => {},
        str    => 0,
    };

    # What alleles should we consider?
    my (@tryAllele, %alleles, %badAlleles);
    if (my $req = $args->{ALLELES} || $args->{BASES}) {
        # Manually provided bases
        push @tryAllele, ref($req) ? @{$req} : ($req);
    }
    foreach my $base (keys %alleles) {
        if (lc($base) eq 'unknown') {
            push @{$rv->{ERROR}}, "Ignored allele specified as '$base'";
            delete $alleles{$base};
            # print "<pre>". $self->branch(\%alleles)."</pre>";
        }
    }

    my $snpLoc = $args->{SNP} || $args->{LOC} || $args->{LOCATION};
    if ($snpLoc && !$args->{NOEXPAND}) {
        my $minFreq = $args->{FREQ};
        # $rv->{minFreq} = $minFreq;
        # A ::Location object is provided, extract all the alleles
        push @tryAllele, $snpLoc->each_allele
            ( -minfreq  => $minFreq,
              -tosspop  => $args->{TOSSPOP},
              -keepnull => $args->{KEEPNULL} );
    }
    foreach my $al (@tryAllele) {
        if ($al eq '-' || $al eq $su->cached_safe_dna( $al )) {
            # The allele seems ok
            $alleles{uc($al)} = uc($al);
        } else {
            # Amazing the stuff that can work its way into the system
            # Reinforces the need to aggresively validate alleles on load time
            # ss12675119 (and many others) : Alleles '-' and '+'
            # Some TCGA data report alleles '--T', '-G-', '-A', etc
            # Also poor input checking on fallback for TCGA normalization
            # captured 'Unknown' as an allele
            $badAlleles{$al}++;
        }
    }

    # What is the location on the RNA that we are considering?
    my ($l, $r, $upS, $dnS, $str, $hspNum) =
        ($args->{LEFT}, $args->{RIGHT});

    my $rAln;
    if (defined $l && defined $r) {
        # We have provided the lft / rgt coordinates explicitly
    } elsif (defined $args->{START}) {
        $l = $args->{START} - 1;
        if (defined $args->{END}) {
            $r = $args->{END} + 1;
        } else {
            $r = $args->{START} + 1;
        }
    } else {
        my $genoAlgn;
        if ($snpLoc) {
            # We want the system to project the SNP onto the RNA
            $genoAlgn  = $snpLoc->hsp();
        } elsif (defined $args->{GENOSTART}) {
            $self->death("Need to define mechanism for explicit coords");
        }
        if ($genoAlgn) {
            my $anc = $genoAlgn->anchor_id();
            # What alignments should we use for the RNA?
            my @alns;
            if (my $aln = $args->{ALIGN} || $args->{HSP} || $args->{RANGE}) {
                # A specific alignment object is being provided
                @alns = ($aln);
            } else {
                # Use all available alignments
                @alns = $self->each_alignment();
            }
            
            my $rSid = $self->acc_id();
            my %results;
            my @cpx;
            my $done = 0;
            my (%strs, %hspCounts);
            foreach my $aln (@alns) {
                my $ancInd = $aln->seqindex( $anc );
                # Some alignments may be irrelevant
                next unless (defined $ancInd);
                my @alHsps = $aln->hsps();
                $hspCounts{$#alHsps+1}++;
                $done++;
                $aln->set_anchor( $anc );
                my ($int, $isoA, $isoB) = $aln->intersect
                    ( -hsp => $genoAlgn, -isob => 1 );
                my @rKey;
                if ($int) {
                    # warn join("\n", "Intersection Report for Impact:",$genoAlgn->hsp_text(), $int->hsp_text) if ($rAcc =~ /(002654|001206797)/);
                    if (my $is = $int->strand_for_seq( $rSid )) {
                        $strs{$is}++;
                    }
                    my @hsps = $int->hsps_for_seq( $rSid );
                    # There is an intersection
                    $rKey[4] = $int->strand_for_seq( $rSid );
                    if ($#hsps != 0) {
                        push @{$rv->{ERROR}},
                        "Intersecting HSP is present in two or more HSPs: (".
                            join(',', map { "[$_->[0],$_->[1]]" } @hsps).")";
                        push @cpx, $aln;
                        next;
                    }
                    my $rLen  = $int->length();
                    my $lDiff = $genoAlgn->length() - $rLen;
                    if ($lDiff && $rLen) {
                        my $err = "The RNA occupies a ".int($lDiff)."bp ";
                        if ($lDiff > 0) {
                            $err .= "smaller area than the variant footprint. This is likely because the variant overlaps both an exon and non-exonic sequence";
                        } else {
                            $err .= "larger area than the variant footprint. Huh. I really don't know how that can happen.";
                        }
                        push @{$rv->{ERROR}}, $err;
                        push @cpx, $aln;
                        next;
                    }
                    $rKey[0] = $hsps[0][0];
                    $rKey[1] = $hsps[0][1];

                } elsif ($#{$isoB} == -1) {
                    push @{$rv->{ERROR}},
                    "Failed to calculate intersection to HSPed object";
                    next;
                } elsif ($#{$isoB} != 0) {
                    push @{$rv->{ERROR}},
                    "Non-intersecting HSP is present in two or more HSPs: (".join(',', map { "[$_->[0],$_->[1]]" } @{$isoB}).")";
                    push @cpx, $aln;
                    next;

                } else {
                    # The query does not overlap, and exists in a unique loc
                    $rKey[2] = join(",", @{$isoB->[0][1] || []});
                    $rKey[3] = join(",", @{$isoB->[0][2] || []});
                    my $gStr = $genoAlgn->strand();
                    $rKey[4] = $aln->strand_for_seq( $rSid ) * $gStr;
                }
                push @{$results{join("\t", map {defined $_ ? $_ : ""} @rKey)}},
                [ \@rKey, $aln ];
            }
            # Ok, now check for consistency
            my @res = values %results;
            my ($best) = $#res == -1 ? undef :
                sort { $a->[1]->howbad() <=>
                           $b->[1]->howbad() } @{$res[0]};
            if ($best) {
                $rAln   = $best->[1];
            }
            if ($#res == 0) {
                # Hooray!
                ($l, $r, $upS, $dnS, $str) = @{$best->[0]};
            } else {
                my @ss = keys %strs;
                $str = $ss[0] if ($#ss == 0);
                $rv->{halt} = "Can not calculate RNA location";
                if ($#res != -1) {
                    # Boo
                    push @{$rv->{ERROR}},
                    "Inconsistent impacts detected. Please specify a single alignment for the RNA";
                } elsif (!$done) {
                    $rv->{imp} = "NONE";
                } elsif ($#cpx != -1) {
                    # Complex intersection
                    unless ($best) {
                        ($rAln) = sort { $a->howbad() <=> $b->howbad() } @cpx;
                    }
                    $rv->{imp} = "CPX";
                }
            }
            my @hnums = keys %hspCounts;
            $hspNum = $hnums[0] if ($#hnums == 0);
        }
    }
    $rv->{align} = $rAln ? $rAln->pkey() : 0;
    $rv->{rnaHB} = !$rAln ? 100 : 
        defined $rAln->howbad() ? $rAln->howbad() : 100;
    $rv->{str} = $str ||= 0;
    my @alls = keys %alleles;
    if ($str < 0) {
        # We need to revcom the alleles
        foreach my $al (@alls) {
            $alleles{$al} = $su->revcom( $al );
        }
    }
    $rv->{var} = { map { $_ => [undef, undef, $alleles{$_} ] } @alls };
    my @badA = sort keys %badAlleles;
    unless ($#badA == -1) {
        my $bT = "'".join("','", @badA)."'";
        $self->msg_once
            ("$rAcc [".($l||'?')."^".($r||'?')."] Illegal alleles: $bT")
            if ($ml->verbosity());
        push @{$rv->{ERROR}}, "Illegal Alleles : $bT";
    }

    if (exists $rv->{halt}) {
        $self->bench_end();
        return $rv;
    }

    my $spliceDist = 2;
    unless (defined $l && defined $r) {
        # The location does not fall within the RNA
        if (defined $upS || defined $dnS) {
            ($upS, $dnS) = ($dnS, $upS) if ($str < 0);
            my ($upInd, $dnInd);
            ($upS, $upInd) = split(',', $upS);
            ($dnS, $dnInd) = split(',', $dnS);
            $rv->{upstream}   = $upS if ($upS);
            $rv->{downstream} = $dnS if ($dnS);
            my ($upEx, $dnEx) = map { 
                defined $_ && $_ ne '' ? $_ + 1 : 0 } ($upInd, $dnInd);
            if ($str < 0) {
                if ($hspNum) {
                    map { $_ = $hspNum - $_ + 1} ($upEx, $dnEx);
                }
            }
            my $note;
            if (defined $upS && defined $dnS) {
                # The location is between HSPs
                my ($min, $side);
                if ($upS < $dnS) {
                    $side = '5';
                    $min  = $upS;
                    $note = sprintf("%dbp&larr;Ex%d", $upS, $upEx);
                } else {
                    $side = '3';
                    $min = $dnS;
                    $note = sprintf("Ex%d&rarr;%dbp", $dnEx, $dnS);
                }
                # my ($min) = sort {$a <=> $b} ($upS, $dnS);
                $rv->{spliceDistance} = $spliceDist;
                if ($min < $spliceDist) {
                    $rv->{imp} = "SP$side";
                } else {
                    $rv->{imp} = "INT";
                }
            } else {
                $rv->{imp} = "GEN";
                if (defined $upS) {
                    $rv->{side} = "5'" ;
                    $note = sprintf("%dbp&larr;Ex%d", $upS, $upEx);
                } elsif (defined $dnS) {
                    $rv->{side} = "3'";
                    $note = sprintf("Ex%d&rarr;%dbp", $dnEx, $dnS);
                } else {
                    push @{$rv->{ERROR}}, "Unclear non-interesecting location";
                }
            }
            $rv->{impNote} = $note;
            if ($note && $str < 0 && !$hspNum) {
                push @{$rv->{ERROR}}, "Unable to calculate total number of exons; -1 strand results will report wrong exon number.";
            }
        } else {
            # The location is fully undefined
            push @{$rv->{ERROR}}, "Failed to find RNA boundaries";
            $rv->{imp} = "UNK";
        }
        $self->bench_end();
        return $rv;
    }

    # The location is falling within an exon
    $rv->{lft} = $l;
    $rv->{rgt} = $r;
    my @feats = $self->overlapping_feature_ids( $l, $r );
    foreach my $fdat (@feats) {
        my ($fid, $srcid) = @{$fdat};
        $rv->{features}{$fid}{$srcid} = 1;
    }
    while (my ($fid, $srcs) = each %{$rv->{features} || {}}) {
        $rv->{features}{$fid} = [ keys %{$srcs} ];
    }

    $rv->{nucPos} = $ml->pretty_flanks( $l, $r );
    my $seq = $self->seq();
    my $len = $r - $l - 1;
    my $wt  = uc($len ? substr($seq, $l, $len) : "-");
    # die $self->to_text()."$l,$r $len = ".length($seq) unless ($wt);
    my $wtAl = $rv->{ref} = $str < 0 ? $su->revcom( $wt ) : $wt;
    unless ($alleles{$wtAl}) {
        push @alls, $wtAl;
        $alleles{$wtAl} = $wt;
        $rv->{var}{$wtAl} = [undef, undef, $wt ];
    }
    my @allAlls = values %alleles;
    $rv->{nucNom} = $ml->variant_nomenclature($wt, \@allAlls, $l, $r, 'r');
    my ($cs, $ce) = $self->cds();
    my $fsTok = "FrameShift";
    if ($cs) {
        # We have a CDS defined. These positions are 1-indexed
        if ($ce <= $l) {
            $rv->{imp} = "UT3";
            $rv->{side} = "3'";
        } elsif ($cs >= $r) {
            $rv->{imp} = "UT5";
            $rv->{side} = "5'";
        } else {
            my $coff = $cs - 1;
            $rv->{cdsNom} = $ml->variant_nomenclature
                ($wt, \@allAlls, $l - $coff, $r - $coff, 'c');
            # The variant overlaps the CDS
            # The number of bases in the 5' UTR:
            my $utr5 = $coff;
            # This is the left position in CDS space:
            # (ie $cl = 0 = the position to the left of CDS coord 1)
            my $cl   = $l - $utr5;
            # We also need the left CODON position
            my $cdMod = $cl % 3;
            $rv->{cdpos} = $cdMod + 1;
            my $cdl   = $cl - $cdMod;
            # Get the codon length:
            my $cdLen = $len + $cdMod;
            if (my $lMod = $cdLen % 3) {
                $cdLen += 3 - $lMod;
            }
            # Left coordinate in protein space:
            my $prtL = int(10 * $cdl / 3) / 10;
            my $prtR = int(10 * ($prtL + ($cdLen / 3) + 1)) / 10;

            my $refSeq = substr($seq, $utr5 + $cdl, $cdLen);
            my $refPrt = $su->cached_translation($refSeq);
            my (%lens, %prts, %b2p);
            my $imp = $#alls == 0 ? 'COD' : "SYN";
            foreach my $al (@alls) {
                my $base = $alleles{$al};
                $base = "" if ($base eq '-');
                my $codon = $refSeq;
                substr($codon, $cdMod, $len) = $base;

                my $l = length($codon);
                $lens{$l}++;
                my $prt = $l % 3 ? 
                    $fsTok : $su->cached_translation($codon);
                $prts{$prt}++;
                $b2p{$alleles{$al}} = $prt;
                $rv->{var}{$al} = [ $codon, $prt, $base ];
            }
            if (int($prtL) == $prtL && int($prtR) == $prtR) {
                $rv->{protPos} = $ml->pretty_flanks( $prtL, $prtR );
                $rv->{protNom} = $ml->variant_nomenclature
                    ($b2p{$wt}, [keys %prts], $prtL, $prtR, 'p');
            } else {
                push @{$rv->{ERROR}}, "Error calculating codon position (protein bounds = $prtL..$prtR)";
            }
            if (exists $prts{$fsTok}) {
                $imp = "FRM";
            } elsif (exists $prts{'*'}) {
                $imp = "STP";
            } else {
                my @lens = keys %lens;
                if ($#lens == 0) {
                    my @prts = keys %prts;
                    $imp = "NON" if ($#prts != 0);
                } else {
                    $imp = "DEL";
                }
            }
            $rv->{imp} = $imp;
        }
    } else {
        # No CDS, but still inside the gene
        $rv->{imp} = "NCR";
    }
    $self->bench_end();
    return $rv;
}

sub length {
    my $self = shift;
    if (my $seq = $self->seq()) {
        return length($seq);
    }
    return undef;
}

sub seq_id {
    return shift->_generic_get_id( 'SEQ', @_ );
}

sub make_transient {
    my $self = shift;
    if (my $pk = $self->{PKEY}) {
        $self->{PKEY_MASKED} = $pk;
    }
    $self->{PKEY} = 0;
}

# ::RNA
*handle  = \&pkey;
sub pkey {
    my $self = shift;
    unless (defined $self->{PKEY}) {
        $self->bench_start();
        $self->{PKEY} = $self->_generic_get_pkey
            ([$self->acc_id(), $self->src_id()], 'rna',
             'rna_id', ['acc_id','src_id'] );
        $self->bench_end();
    }
    return $self->{PKEY};
}

sub set_seq {
    my $self = shift;
    if (my $nv = shift) {
        $self->bench_start();
        my $ml   = $self->maploc();
        my $sid  = $self->{SEQ_ID} = $ml->text_to_pkey( $nv );
        my $pkey = $self->pkey();
        if ($pkey && !$self->{READ_ONLY}) {
            my $set  = $ml->{STH}{SET_RNA_SEQ} ||= $ml->dbh->prepare
                ( -name => "Add sequence to RNA",
                  -sql  => "UPDATE rna SET seq_id = ? WHERE rna_id = ?" );
            $set->execute($sid, $pkey);
        }
        $self->bench_end();
    }
}

*set_cds = \&cds;
sub cds {
    my $self = shift;
    my ($s, $e) = @_;
    if (defined $s) {
        $self->{CDS} = [ $s, $e ];
        my $pkey = $self->pkey();
        if ($pkey && !$self->{READ_ONLY}) {
            $self->bench_start();
            my $ml  = $self->maploc();
            my $set = $ml->{STH}{SET_RNA_CDS} ||= $ml->dbh->prepare
                ( -name => "Update RNA CDS",
                  -sql  => "UPDATE rna SET cds_start = ?, cds_end = ? WHERE rna_id = ?" );
            $set->execute( $s, $e, $pkey );
            $self->bench_end();
        }
    }
    return @{$self->{CDS} || [0,0]};
}

sub cds_range {
    my $self = shift;
    unless (defined $self->{CDS_RANGE}) {
        my ($s, $e) = $self->cds();
        if ($s) {
            my $aid  = $self->accession_id();
            my $hsp  = $self->{CDS_RANGE} =
                BMS::SnpTracker::MapLoc::Range->new( $self->maploc );
            $hsp->{STRAND} = [1];
            $hsp->hsps( [[ $s-1, $e+1 ]] );
            $hsp->seq_ids( $aid );
            $hsp->set_anchor( $aid );
        } else {
            $self->{CDS_RANGE} = 0;
        }
    }
    return $self->{CDS_RANGE} || undef;
}

sub cds_genomic_ranges {
    my $self  = shift;
    my ($s, $e) = $self->cds();
    return () unless ($s);
    my $args = $self->parseparams( @_ );
    my $ml   = $self->maploc();
    #my $acc  = $self->acc();
    my $accid = $self->acc_id();
    my @alns = $self->alignments_from_parameters( $args );
    my %pair = ( ATG => [ $s - 1, $s + 3],
                 Stop => [$e - 3, $e + 1] );
    my @rrng;
    foreach my $key (keys %pair) {
        my $hsp = $pair{$key};
        my $rng = BMS::SnpTracker::MapLoc::TaggedRange->new( $ml );
        $rng->{STRAND} = [1];
        $rng->hsps( [$hsp] );
        $rng->seq_ids( $accid );
        $rng->set_anchor( $accid );
        $rng->name( $key );
        push @rrng, $rng;
    }
    my @rv;
    foreach my $galn (@alns) {
        $galn->set_anchor( $accid );
        foreach my $rng (@rrng) {
            my ($int) =  $galn->intersect( -hsp => $rng );
            if ($int) {
                # my @csi = $int->chr_seq_info();
                push @rv, $int;
                $int->tag('Key', $rng->name());
                # $int->name($rng->name);
                # warn $int->full_text();
            }
        }
    }
    return @rv;
}

sub rna_to_protein_coordinate {
    my $self = shift;
    my $rpos = shift;
    return undef unless ($rpos);
    my ($cs, $ce) = $self->cds();
    return undef unless ($cs && $ce);
    return undef if ($rpos < $cs || $rpos > $ce);
    my $rp = $rpos - $cs;
    my $ppos = int($rp / 3) + 1;
    # In scalar context just return the protein coordinate
    return $ppos unless (wantarray);
    my $phase = $rp % 3;
    return ($ppos, $phase);
}

sub clear_cds {
    shift->set_cds('','');
}

*set_anomaly = \&anomaly;
sub anomaly {
    my $self = shift;
    my $nv   = shift;
    if (defined $nv) {
        $self->{ANOMALY} = $nv;
        my $pkey = $self->pkey();
        if ($pkey && !$self->{READ_ONLY}) {
            $self->bench_start();
            my $ml   = $self->maploc();
            my $set  = $ml->{STH}{SET_RNA_ANOMALY} ||= $ml->dbh->prepare
                ( -name => "Set RNA anomaly",
                  -sql  => "UPDATE rna SET anomaly = ? WHERE rna_id = ?" );
            $set->execute($nv, $pkey);
            $self->bench_end();
        }
    }
    return $self->{ANOMALY} || "";
}

sub variant_seq {
    my $self  = shift;
    my $args  = $self->parseparams( @_ );
    my $ml    = $self->maploc();
    my $seq   = $self->seq();
    my @locs  = $self->each_location( @_ );
    foreach my $loc (@locs) {
        $self->death("WORKING!");
    }
    return $seq;
}

sub read {
    my $self = shift;
    unless ($self->{READ_ALREADY}) {
        $self->bench_start();
        $self->read_tags();
        $self->read_alignments();
        $self->bench_end();
        $self->read_features();
        $self->{READ_ALREADY} = 1;
    }
    return $self->{READ_ALREADY};
}

# ::RNA
sub write {
    my $self = shift;
    $self->write_tags();
    $self->write_alignments();
}

