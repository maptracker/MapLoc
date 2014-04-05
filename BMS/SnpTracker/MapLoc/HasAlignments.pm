# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::SnpTracker::MapLoc::HasAlignments;
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
use BMS::SnpTracker::MapLoc::TaggedRange;

sub add_alignment {
    my $self = shift;
    foreach my $aln (@_) {
        $self->{ALIGNS}{ $aln->handle() } ||= $aln;
    }
}

sub has_alignment {
    my $self = shift;
    my $aln  = shift || return 0;
    my $hand = $aln->handle();
    return exists $self->{ALIGNS}{ $hand } && $self->{ALIGNS}{ $hand } ? 1 : 0;
}

sub delete_alignment {
    my $self = shift;
    # $self->bench_start();
    my $aln = shift;
    return 0 unless ($self->has_alignment($aln));
    my $id = $aln->pkey();
    # warn "WOULD KILL:"; warn $aln->full_text(); return;

    delete $self->{ALIGNS}{ $aln->handle() };
    $aln->delete();
    return $aln;
}

sub each_alignment_id {
    return keys %{shift->{ALIGNS}};
}

sub each_alignment {
    my $self = shift;
    $self->bench_start();
    my $hb   = shift;
    my $alns = $self->read_alignments();
    my @rv;
    foreach my $id (keys %{$alns}) {
        if (my $aln = $alns->{$id} ||= $self->maploc->get_alignment($id)) {
            next if (defined $hb && $aln->howbad() > $hb);
            push @rv, $aln;
        } else {
            $self->err("Failed to recover alignment for aln_id = $id");
        }
    }
    $self->bench_end();
    return @rv;
}

sub _specific_alignment {
    my $self = shift;
    my $aid  = shift;
    return undef unless ($aid);
    my $hash = $self->read_alignments();
    return 0 unless (exists $hash->{$aid});
    return $hash->{$aid} ||= $self->maploc->get_alignment($aid);
}

sub each_alignment_for_build {
    my $self  = shift;
    my $build = uc(shift || "");
    my $hb    = shift;
    my @rv;
    foreach my $aln ($self->each_alignment( $hb )) {
        my $sinf = $aln->chr_seq_info();
        next if ($#{$sinf} != 0);
        push @rv, $aln if ($#{$sinf} == 0 && uc($sinf->[0][3]) eq $build);
    }
    return @rv;
}

sub alignments_from_parameters {
    my $self = shift;
    my $args = shift;
    my @rv;
    if (my $alns = $args->{ALN} || $args->{ALNS} || $args->{ALIGN}) {
        my $r = ref($alns);
        @rv = $r eq 'ARRAY' ? @{$alns} : ($alns);
    } else {
        my $hb = $args->{HOWBAD};
        if (my $bld = $args->{BUILD}) {
            if (lc($bld) eq 'best') {
                my %byBld;
                foreach my $aln ($self->each_alignment( $hb )) {
                    my $bld = $aln->chr_builds() || "UNK";
                    # print "<pre>".$aln->to_text()."</pre>";
                    push @{$byBld{$bld}}, $aln;
                }
                my $best = $self->best_build( keys %byBld );
                @rv = @{$byBld{$best || ""} || []};
            } else {
                @rv = $self->each_alignment_for_build( $bld, $hb );
            }
        } else {
            @rv = $self->each_alignment( $hb );
        }
    }
    return @rv;
}

sub chr_hsps {
    my $self = shift;
    my %rv;
    my $args = $self->parseparams( @_ );
    foreach my $aln ($self->alignments_from_parameters( $args )) {
        my ($hsps, $build, $chr) = $aln->chr_hsps_for_alignment( );
        next if (!$hsps || $#{$hsps} == -1);
        if ($args->{EXONONLY}) {
            # Footprint will preserve exon structure
            push @{$rv{$build}{$chr}}, @{$hsps};
        } else {
            # Devolve into a single HSP
            push @{$rv{$build}{$chr}}, [ $hsps->[0][0], $hsps->[-1][1] ];
        }
    }
    return \%rv;
}

sub hsp_text {
    my $self = shift;
    return join("\n", map { $_->hsp_text() } $self->each_alignment());
}

# ::HasAlignments
sub full_text {
    my $self = shift;
    return $self->to_text() . $self->hsp_text();
}

sub each_location {
    my $self = shift;
    $self->bench_start();
    my $args  = $self->parseparams( @_ );
    my @alns  = $self->alignments_from_parameters( $args );
    my $flank = $args->{FLANK};
    my $ml    = $self->maploc();
    # map { warn $_->to_text() } @alns;
    # warn "foo";
    my @locs  = $ml->location_query( -aln => \@alns, @_ );
    $self->bench_end();
    return @locs;
}

# ::HasAlignments
sub derived_locations {
    my $self = shift;
    my $args = $self->parseparams( @_ );
    my $dtxt = $args->{LOC};
    my $vb   = $args->{VERBOSE};
    return wantarray ? () : undef unless ($dtxt);
    my $ml   = $self->maploc();
    my @alns = $self->alignments_from_parameters( $args );
    my ($cr, $indel, $l, $r, $typ, $alleles);
    my $crRE  = '([crp]\.)?';
    my $diRE  = '(del|ins)?';
    my $sepRE = '([\.\-\^]+)';
    my $allRE = '([^\d]*)';
    if ($dtxt =~ /^\s*$crRE$diRE\s*(\d+)\s*$sepRE\s*(\d+)$allRE$/i) {
        ($cr, $indel, $l, $typ, $r, $alleles) = 
            (lc($1 || ''), lc($2 || ''), $3, $4, $5, $6 || "");
    } elsif ($dtxt =~ /^\s*$crRE$diRE\s*(\d+)$allRE$/) {
        $typ = "";
        ($cr, $indel, $r, $l, $alleles) = 
            (lc($1 || ''), lc($2 || ''), $3, $3, $4 || '');
    } else {
        my $acc = $self->accession();
        $self->msg("[?]", "Could not parse coordinates/alleles for '$dtxt' on $acc") if ($vb);
        return wantarray ? () : [];
    }
    if ($alleles =~ /^\s*(del|ins)\s*(.+)$/i) {
        ($indel, $alleles) = (lc($1), $2);
    }
    if ($indel) {
        if ($indel eq 'del') {
            if ($alleles) {
                $alleles = "$alleles>-" unless ($alleles =~ /[>\-]/);
            } else {
                $alleles = "-";
            }
        } elsif ($indel eq 'ins') {
            if ($alleles) {
                $alleles = "->$alleles" unless ($alleles =~ /[>\-]/);
            } else {
                $alleles = "-";
            }
            if ($r = $l + 1) {
                $typ = '^';
            } else {
                my $acc = $self->accession();
                $self->msg("[?]", "Insertion coordinats for '$dtxt' on $acc do not make sense") if ($vb);
                return wantarray ? () : [];
            }
        } else {
            my $acc = $self->accession();
            $self->msg("[?]", "Programming error : InDel code '$indel' computed for '$dtxt' on $acc");
                return wantarray ? () : [];
        }
    }
    if ($cr) {
        $cr =~ s/\.$//;
        if ($cr eq 'c' || $cr eq 'p') {
            # The coordinates are in cDNA (c) or protein (p) space
            if ($cr eq 'p') {
                # Convert to cDNA space
                $l = (($l - 1) * 3) + 1;
                $r = (($r - 1) * 3) + 3;
            }
            if ($self->can('cds')) {
                my ($s) = $self->cds();
                if ($s) {
                    $s--;
                    $l += $s;
                    $r += $s;
                } else {
                    my $acc = $self->accession();
                    $self->msg("[?]", "The CDS was not set for $acc, can not adjust cDNA-referenced allele '$dtxt'") if ($vb);
                    return wantarray ? () : [];
                }
            } else {
                my $acc = $self->accession();
                $self->msg("[?]", "No CDS method defined for $acc, can not adjust cDNA-referenced allele '$dtxt'") if ($vb);
                return wantarray ? () : [];
            }
        }
    }
    unless ($typ eq '^') {
        $l--;
        $r++;
    }
    my $rng   = BMS::SnpTracker::MapLoc::TaggedRange->new( $ml );
    my $accid = $self->acc_id();
    $rng->{STRAND} = [1];
    $rng->hsps( [[$l, $r]] );
    $rng->seq_ids( $accid );
    $rng->set_anchor( $accid );
    $rng->name( $dtxt );
    my @rv;
    foreach my $galn (@alns) {
        $galn->set_anchor( $accid );
        my ($int) = $galn->intersect( -hsp => $rng );
        unless ($int) {
            if ($vb) {
                # print "<pre>". $galn->full_text()."</pre>";
                my $acc = $self->accession();
                $self->msg("[?]", "Failed to intersect '$dtxt' on $acc");
            }
            next;
        }
        # $ml->preprint($int->full_text());
        my $gId  = $int->other_seq_ids( $accid ) || 0;
        my @hsps = $int->hsps_for_seq( $gId );
        my @lrs  = sort { $a <=> $b } map { @{$_} } @hsps;
        my $str  = $int->strand_for_seq( $gId );
        my $cnm  = $int->chr_name();
        my $bld  = $int->chr_build();
        push @rv, [ [ $cnm, $lrs[0], $lrs[-1], $bld ], $alleles, $str ];
        # $int->tag('Key', $rng->name());
        # $self->maploc->preprint( "$cnm=\n".$int->full_text() );
    }
    # $self->maploc->prebranch( \@rv );
    return wantarray ? @rv : \@rv;
}

sub each_local_location {
    my $self = shift;
    my $args = $self->parseparams( @_ );
    my @alns = $self->alignments_from_parameters( $args );
    my $ml   = $self->maploc();
    my $sid  = $self->acc_id();
    my $rv   = {};
    my (%seen, %seenpop);
    foreach my $aln (@alns) {
        $aln->anchor_to_genome();
        my $locs = $ml->location_query( @_, -aln => $aln );
        while (my ($mode, $popids) = each %{$locs->{pop_id}}) {
            map { push @{$rv->{pop_id}{$mode}}, $_
                      unless ($seenpop{$mode}{$_}++) } @{$popids};
        }
        if (my $errs = $locs->{errors}) {
            map { push @{$rv->{errors}}, $_ unless ($seen{$_}++) } @{$errs};
        }
        foreach my $lid ( @{$locs->{loc_id}} ) {
            my $loc = $ml->get_location($lid);
            my ($int, $isoA, $isoB) = $aln->intersect
                ( -hsp => $loc->hsp(), -isob => 1 );
            if (($int->length() != $loc->length()) && !$args->{KEEPODD}) {
                # warn $int->length()." != ".$loc->length();
                push @{$rv->{reject}}, [$loc, $int];
                next;
            }
            # warn $int->to_text().$self->branch({isoA => $isoA, isoB => $isoB});
            push @{$rv->{loc_hsp}}, [$loc, $int];
        }
    }
    return wantarray ? @{$rv->{loc_hsp}} : $rv;
}

# ::HasAlignments;
sub exon_hsps {
    my $self = shift;
    $self->bench_start();
    $self->read();
    my $args  = $self->parseparams( @_ );
    my @alns  = $self->alignments_from_parameters( $args );
    my $range = $args->{RANGE} || 0;
    my @rv;
    foreach my $aln (@alns) {
        my ($hsps, $build, $chr) = $aln->chr_hsps_for_alignment( );
        push @rv, [$hsps, $build, $chr, $aln];
    }
    $self->bench_end();
    return @rv;
}

sub _merged_hsp_ranges {
    my $self = shift;
    $self->bench_start();
    my ($subsets, $range) = @_;
    $range ||= 0;
    # Group the sets by build and chromosome
    my %builds;
    foreach my $hset (@{$subsets}) {
        my ($hsps, $build, $chr, $aln) = @{$hset};
        my $sc  = $aln->score();
        my $str = $aln->strand();
        foreach my $hsp (@{$hsps}) {
            push @{$builds{$build}{$chr}}, 
            [$hsp->[0] - $range, $hsp->[1] + $range, $str, $sc ];
        }
    }
    # Now collapse overlapping ranges
    while (my ($build, $cHash) = each %builds) {
        while (my ($chr, $hArr) = each %{$cHash}) {
            my @sorted = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$hArr};
            my @merged = (shift @sorted);
            foreach my $hsp (@sorted) {
                # The last loaded HSP:
                my ($pl, $pr, $pStr, $psc) = @{$merged[-1]};
                if ($hsp->[0] > $pr) {
                    # This hsp does not overlap
                    push @merged, $hsp;
                    next;
                }
                # We overlap
                $self->msg("[TODO]", "Divergent strands on $chr $pl..$pr for "
                           . $self->accession()) if ($pStr != $hsp->[2]);
                $self->msg("[TODO]", "Divergent scores on $chr $pl..$pr for "
                           . $self->accession()) if ($psc != $hsp->[3]);
                if ($hsp->[0] <= $pr) {
                    # Complete overlap
                    next;
                }
                # Update the right coordinate
                $merged[-1][1] = $pr;
            }
            $cHash->{$chr} = \@merged;
        }
    }
    $self->bench_end();
    return \%builds;
}

sub read_alignments {
    my $self = shift;
    my $id   = $self->acc_id();
    return $self->{ALIGNS} || {} unless ($id);
    unless ($self->{READ_ALN_ALREADY}) {
        $self->bench_start();
        $self->{ALIGNS} ||= {};
        my $ml   = $self->maploc();
        my @alnids = $ml->alignment_ids_for_seqid( $id );
        map { $self->{ALIGNS}{$_} ||= undef } @alnids;
        $self->bench_end();
        $self->{READ_ALN_ALREADY} = 1;
    }
    return $self->{ALIGNS};
}

sub write_alignments {
    my $self = shift;
    $self->bench_start();
    my $clear = shift;
    # read_alignments() will populate the alignment hash with all alignment IDs
    # available from the database, without generating alignment objects
    $self->read_alignments();
    my (@oldids, %grouped);
    my $acc = $self->acc();
    while (my ($alnid, $aln) = each %{$self->{ALIGNS}}) {
        if ($aln) {
            # This is a defined alignment, which we will keep and store
            # We want to organize them by build when possible
            my $group  = "Generic";
            my $gdna   = $aln->other_seqs( $acc );
            if ($gdna  =~ /^([^\.]+)\.([^\.]+)\.([^\.]+)\.([^\.]+)$/) {
                # Genomic sequence in BMS-standard format
                # use the build as the group
                $group = $4;
            }
            push @{$grouped{$group}}, $aln;

            # If we have a 'duplicate key value violates unique constraint'
            # error inside a BEGIN block, it will abort the block even though
            # we are set up to ignore the error. Call seqs() and source()
            # to premptively trigger those if needed:
            $aln->seqs();
            $aln->source();
        } else {
            # This is an alignment in the database that we do not have now
            push @oldids, $alnid;
        }
    }
    # Ok, now we will write the alignments and clear ond ones if requested
    my $ml = $self->maploc();
    my $dbh = $ml->dbh();

    $dbh->begin_work();
    if ($clear) {
        # Clear locations if requested, inside the transaction block
        $ml->delete_alignment_by_id( @oldids );
    }
    while (my ($group, $arr) = each %grouped) {
        # Find the best score for this group
        my ($best) = sort { $b <=> $a } map { $_->score() || -1 } @{$arr};
        # Use it to set howbad
        map { $_->howbad( $best - $_->score ) } @{$arr};
        map { $_->write('carefully') } @{$arr};
    }
    $dbh->commit();
    $self->bench_end();
}

# ::HasAlignments
sub nearby_rnas {
    my $self = shift;
    $self->bench_start();
    $self->read();
    my $args     = $self->parseparams( @_ );
    my $ml       = $self->maploc();
    my @alns     = $self->alignments_from_parameters( $args );
    my $range    = $args->{RANGE};
    my $hb       = $args->{HOWBAD};
    my $defFlank = $ml->default_loc_to_rna_distance();
    my $defHB    = $ml->default_howbad();
    my $isDef    = (defined $hb && $hb != $defHB) ||
        (defined $range && $range != $defFlank) ? 0 : 1;
    $range       = $defFlank unless (defined $range);
    $hb          = $defHB unless (defined $hb);
    my $acc   = $self->acc();
    my $cache = $self->{NEARBY_RNA_CACHE} ||= {};
    my %rv;

    foreach my $aln (@alns) {
        my $aKey = $aln->chr_location_tag();
        next unless ($aKey);
        if ($rv{$aKey}) {
            push @{$rv{$aKey}{err}}, "Multiple alignments in this build with this range";
            next;
        }
        my ($hsps, $build, $chr) = $aln->chr_hsps_for_alignment( );
        my ($chrId) = $aln->chr_seq_info( 1 );
        my $str  = $aln->strand() || 0;
        my $len  = $aln->length();
        my $data = $rv{$aKey} = {
            chrtag => $aKey,
          #  hsps   => $hsps,
            build  => $build,
            chr    => $chr,
            str    => $str,
            acc    => $acc,
            len    => $len,
        };
        my $apkey = $aln->pkey();
        my $ckey = join('-', $apkey, $hb);
        if (my $hits = $cache->{$ckey}) {
            $data->{nearby} = $hits;
            next;
        }
        my ($min, $max) = ($hsps->[0][0], $hsps->[-1][1]);
        my %exons = map { join('-', @{$hsps->[$_]}) => $_ + 1 } (0..$#{$hsps});
        my @near;
        $aln->set_anchor( $chrId );
        my $rnaRows;
        if ($isDef) {
            my $bulkSrc = $ml->bulk_cache_source( 'nearby_rnas' );
            $rnaRows    = $bulkSrc->{$apkey};
        }
        $rnaRows ||= $ml->overlapping_rnas_for_flanks
            ( $chrId, [$min - $range, $max + $range ]);
        foreach my $row (@{$rnaRows}) {
            my ($alnId, $oseqId, $osc, $l, $r, $ornaId) = @{$row};
            my $orna = $ml->get_rna_by_id( $ornaId );
            $orna->read();
            my $oaln = $orna->_specific_alignment( $alnId );
            next unless ($oaln);
            $oaln->set_anchor( $chrId );
            my ($ohsps) = $oaln->chr_hsps_for_alignment( );
            my ($omin, $omax) = ($ohsps->[0][0], $ohsps->[-1][1]);
            my $olen  = $oaln->length();
            my $dist  = 0;
            if ($omin + 2 > $max) {
                # upstream
                $dist = $omin - $max;
            } elsif ($min + 2 > $omax) {
                # downstream
                $dist = $min - $omax;
            }
            my $ostr  = $oaln->strand() || 0;
            my $oacc  = $orna->acc();
            my $odata = $near[$#near + 1 ] = {
               # hsps   => $ohsps,
                str    => $ostr,
                match  => 0,
                dist   => $dist,
                acc    => $oacc,
                len    => $olen,
            };
            if ($dist) {
                $odata->{text} = "Non-overlapping";
                next;
            } elsif ($ostr != $str) {
                $odata->{text} = "Opposite strands";
                next;
            }
            # Ok, the RNAs overlap and are on the same strand
            # Characterize the intersection of the two
            my $int   = $aln->intersect( -hsp => $oaln );
            unless ($int) {
                $odata->{text} = "No overlap found";
                next;
            }
            my $ilen  = $int->length();
            my $sc = $odata->{match} = int(0.5 + (10000 * $ilen * 2) / 
                                           ($len + $olen)) / 100;
            if (!$sc) {
                $odata->{text} = "Overlap but not in exons";
                next;
            } elsif ($len == $olen && $ilen == $olen) {
                if ($oacc eq $acc) {
                    # Found the query itself, do not bother reporting
                    pop @near;
                } else {
                    $odata->{text} = "Exactly the same";
                }
                next;
            }
            my $qind  = $int->seqindex( $acc );
            
            my @osort = sort { $a->[0] <=> $b->[0] } @{$ohsps};
            my %exSub = map { join('-', @{$osort[$_]}) => $_ + 1} (0..$#osort);
        }
        $data->{nearby} = $cache->{$ckey} =
            [ sort { $b->{match} <=> $a->{match} || 
                         $a->{dist} <=> $b->{dist} } @near ];
    }
    $self->bench_end();
    return \%rv;
}

return 1;
