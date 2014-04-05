# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::SnpTracker::MapLoc::Featured;
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
use BMS::Utilities;
use BMS::SnpTracker::MapLoc::Common;
use BMS::SnpTracker::MapLoc::TaggedRange;

use vars qw(@ISA);
@ISA      = qw(BMS::SnpTracker::MapLoc::Common);


sub obj_type { return "Featured"; }

sub set_feature {
    my $self = shift;
    my @coords;
    if (ref($_[0])) {
        # Passing (we hope) 2D array
        @coords = @{shift @_};
    } else {
        # Passing single Start, End Pair
        @coords = [ shift @_, shift @_ ];
    }
    my @flanks;
    my $ml = $self->maploc();
    foreach my $crd (@coords) {
        my ($l, $r) = $ml ->position_to_flanks(@{$crd});
        push @flanks, [$l,$r];
    }
    return $self->set_feature_by_flank( \@flanks, @_);
}

sub clear_features {
    my $self = shift;
    delete $self->{FEAT_RANGE_OBJS};
    return if ($self->{READ_ONLY});
    $self->{FEATURES} = {};
    my $pkey = $self->pkey();
    return unless ($pkey);
    my $auth = shift;
    $self->bench_start();
    my $ml = $self->maploc();
    my $clear =  $ml->{STH}{DELETE_FEATURES} ||= $auth ?
        $ml->dbh->prepare
        ( -name => "Remove all features for an object",
          -sql  => "DELETE FROM feature WHERE host_id = ? AND src_id = ?" )
        : $ml->dbh->prepare
        ( -name => "Remove all features for an object",
          -sql  => "DELETE FROM feature WHERE host_id = ?" );
    my $clearHSP =  $ml->{STH}{DELETE_FEATURE_HSPS} ||= $auth ?
        $ml->dbh->prepare
        ( -name => "Remove all feature HSPs for an object",
          -sql  => "DELETE FROM range WHERE rng_id IN ".
          "( SELECT rng_id FROM feature WHERE host_id = ? AND src_id = ?)" )
        : $ml->dbh->prepare
        ( -name => "Remove all feature HSPs for an object",
          -sql  => "DELETE FROM range WHERE rng_id IN ".
          "( SELECT rng_id FROM feature WHERE host_id = ?)" );
    my @binds = ($pkey);
    push @binds, $ml->text_to_pkey($auth) if ($auth);
    $clearHSP->execute( @binds );
    $clear->execute( @binds );
    $self->bench_end();
}

sub set_feature_by_flank {
    my $self = shift;
    my ($flanks, $nv, $src, $str, $pri) = @_;
    return unless ($nv);
    $self->bench_start();
    unless ($flanks) {
        $self->death("set_feature() must define both left and right flank");
    }
    # Make sure our coordinates are ordered
    my @pts = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]} @{$flanks};

    if (!$str) {
        $str = 0;
    } elsif ($str ne '1' && $str ne '-1' && $str ne '+1') {
        $self->death
            ("If strand is provided to set_feature() it must be 1 or -1");
    }
    $self->death("set_feature() must set the source") unless ($src);
    my $ml   = $self->maploc();
    my ($fid, $srcid, $pid) = map {$ml->text_to_pkey( $_ )} ($nv, $src, $pri);
    my ($p2l, $p2r) = ($pts[0][0], $pts[-1][1]);
    my $rv;
    $self->{FEATURES}{join("\t", $fid, $p2l, $p2r)} = 
        [ \@pts, $str, $srcid, $pid, $fid ];
    if ($self->{READ_ONLY}) {
        $self->bench_end();
        return;
    }
    if (my $pkey = $self->pkey()) {
        my $sths = $ml->{STH};
        my $set = $sths->{ADD_FEATURE} ||= $ml->dbh->prepare
            ( -name => "Add feature row",
              -sql  => "INSERT INTO feature (host_id, feat_id, src_id, pos_to_left, pos_to_right, strand, rng_id) VALUES (?, ?, ?, ?, ?, ?, NULL)",
              -ignore => "duplicate key value" );
        $rv = $set->execute( $pkey, $fid, $srcid, $p2l, $p2r, $str );
        if ($#pts > 0) {
            my $rngid = $ml->dbh->nextval('global_seq');
            my $addRange = $sths->{ADD_RANGE} ||= $ml->dbh->prepare
                ( -name => "Add range",
                  -sql  => "INSERT INTO range (rng_id, ptl, ptr) VALUES (?, ?, ?)");
            foreach my $pt (@pts) {
                $addRange->execute( $rngid, @{$pt} );
            }
            my $setRange = $sths->{SET_RANGE_TO_OBJ} ||= $ml->dbh->prepare
                ( -name => "Set range",
                  -sql  => "UPDATE feature SET rng_id = ? WHERE host_id = ? AND feat_id = ? AND pos_to_left = ? AND pos_to_right = ? AND rng_id != ?");
            $setRange->execute( $rngid, $pkey, $fid, $p2l, $p2r, $rngid);
        } else {
            my $setRange = $sths->{SET_RANGE_TO_OBJ} ||= $ml->dbh->prepare
                ( -name => "Set range",
                  -sql  => "UPDATE feature SET rng_id = ? WHERE host_id = ? AND feat_id = ? AND pos_to_left = ? AND pos_to_right = ? AND rng_id IS NOT NULL");
            $setRange->execute( undef, $pkey, $fid, $p2l, $p2r);
        }
        delete $self->{FEAT_RANGE_OBJS};
    }
    $self->bench_end();
    return $rv;
}

sub read_features {
    my $self = shift;
    my $rv = $self->{FEATURES} ||= {};
    if (!$self->{FEATURES_READ} && $self->pkey()) {
        $self->bench_start();
        my $pkey    = $self->pkey();
        my $ml      = $self->maploc();
        my $bulkSrc = $ml->bulk_cache_source( 'features' );
        $self->{FEATURES_READ}++;
        if (my $cache = $bulkSrc->{$pkey}) {
            foreach my $dat (@{$cache}) {
                my ($fid, $src, $pid, $l, $r, $str, $flanks) = @{$dat};
                $rv->{join("\t", $fid, $l, $r)} = 
                    [$flanks, $str, $src, $pid, $fid];
            }
        } else {
            my $sths = $ml->{STH};
            my $get = $sths->{READ_FEATURES} ||= $ml->dbh->prepare
                ( -name => "Read Features",
                  -sql  => "SELECT feat_id, src_id, pri_id, pos_to_left, pos_to_right, strand, rng_id FROM feature WHERE host_id = ?");
            $get->execute($pkey);
            my $rows = $get->fetchall_arrayref();
            foreach my $row (@{$rows}) {
                my ($fid, $src, $pid, $l, $r, $str, $rng) = @{$row};
                my $flanks = [[$l, $r]];
                if ($rng) {
                    my $getRng = $sths->{READ_FEAT_RNGS} ||= $ml->dbh->prepare
                        ( -name => "Read Feature ranges",
                          -sql  => "SELECT ptl, ptr FROM range where rng_id = ? ORDER BY ptl");
                    $getRng->execute($rng);
                    $flanks = $getRng->fetchall_arrayref();
                }
                $rv->{join("\t", $fid, $l, $r)} = 
                    [$flanks, $str, $src, $pid, $fid];
            }
        }
        delete $self->{FEAT_RANGE_OBJS};
        $self->bench_end();
    }
    return wantarray ? values %{$rv} : $rv;
}

sub overlapping_feature_ids {
    my $self = shift;
    my ($l, $r) = @_;
    my @rv;
    if (defined $l) {
        $r = $l unless (defined $r);
        if (my $w = $r + $l - 1) {
            # As long as the query is not a gap position, pull back the
            # boundaries to get true overlap (not just adjacency)
            $l++;
            $r--;
        }
        my $feats = $self->read_features();
        foreach my $dat (values %{$feats}) {
            my $hit = 0;
            foreach my $fd (@{$dat->[0]}) {
                if ( $fd->[0] < $r && $fd->[1] > $l) {
                    $hit = 1; last;
                }
            }
            push @rv, [ $dat->[4], $dat->[2], $dat->[3] ]
                if ($hit); # FeatID, SRC, PID
        }
    }
    return @rv;
}

sub overlapping_features {
    my $self = shift;
    my @rv   = $self->overlapping_feature_ids( @_ );
    my $ml   = $self->maploc();
    foreach my $row (@rv) {
        $row->[0] = $ml->cached_pkey_to_text( $row->[0] );
        $row->[1] = $ml->cached_pkey_to_text( $row->[1] );
    }
    return @rv;
}

sub each_feature {
    my $self = shift;
    my $rv = $self->{FEAT_RANGE_OBJS};
    unless ($rv) {
        $self->bench_start();
        $rv = $self->{FEAT_RANGE_OBJS} = [];
        my $ml   = $self->maploc();
        my $acc  = $self->accession_id();
        my $feats = $self->read_features();
        while (my ($fkey, $dat) = each %{$feats}) {
            my ($flanks, $str, $src, $pid, $fid) = @{$dat};
            my $rng = BMS::SnpTracker::MapLoc::TaggedRange->new( $ml );
            $rng->{STRAND} = [$str];
            $rng->hsps( $flanks );
            $rng->seq_ids( $acc );
            $rng->set_anchor( $acc );
            $rng->name( $ml->cached_pkey_to_text( $fid ) );
            $rng->source( $ml->cached_pkey_to_text( $src ) );
            $rng->{FEAT_KEY} = $fkey;
            push @{$rv}, $rng;
        }
        $self->bench_end();
    }
    $rv = $self->_filter_features( $rv, @_ );
    return @{$rv};
}

sub _filter_features {
    my $self = shift;
    my $list = shift;
    return $list if ($#_ == -1);
    # At least some requests have been passed
    my $args = $self->parseparams( @_ );
    if (my $maxF = $args->{MAXFRAC}) {
        # Request to perform a filter based on length of feature
        # relative to host object
        $maxF  /= 100 if ($maxF > 1); # In case a percentage was given
        my $len = $args->{LENGTH} || $args->{LEN};
        $len = $self->length() if (!$len && $self->can('length'));
        if ($len && $maxF < 1) {
            $maxF = int( 0.5 + $len * $maxF );
            my @keep;
            foreach my $rng (@{$list}) {
                push @keep, $rng unless ($rng->length() > $maxF);
            }
            $list = \@keep;
        }
    }
    if (my $filt = $args->{FEATFILTER} || $args->{FILTER}) {
        # Request to use a callback function to filter features
        my @keep;
        foreach my $rng (@{$list}) {
            push @keep, $rng unless (&{$filt}( $rng ) );
        }
        $list = \@keep;
    }
    return $list;
}

sub anchored_features {
    my $self = shift;
    $self->bench_start();
    unshift @_, '-aln' if ($#_ == -1);
    my $args = $self->parseparams( @_ );
    my $aln  = $args->{ALN};
    my $anc  = $aln->anchor_id();
    my $acc  = $self->accession();
    #warn $self->pkey();
    # warn "\nAlignment:\n".$aln->hsp_text();

    my (@rv, @efParam);
    if (my $filt = $args->{FEATFILTER} || $args->{FILTER}) {
        push @efParam, ( -filter => $filt );
    }
    if (my $maxF = $args->{MAXFRAC}) {
        push @efParam, ( -maxfrac => $maxF,
                         -length  => $aln->length() );
    }
    
    foreach my $rng ($self->each_feature(@efParam)) {
        # Make sure we are using the anchor for the range:
        my $int = $aln->intersect
            ( -hsp => $rng, -anchor => $rng->anchor_id() );
        next unless ($int);
        # Re-anchor the intersection
        $int->set_anchor( $anc );
        $int->rename_seq( $acc, $rng->name );
        $int->tag("Via", $acc);
        $int->tag("Source", $rng->source());
        # warn $int->to_text() if ($rng->name eq 'Phosphotyrosine');
        #$self->msg("[DEBUG]", "Intersection:", $int->hsp_text());
        push @rv, $int
    }
    $self->bench_end();
    return @rv;
}

sub cx_for_feature {
    my $self = shift;
    my $ml   = $self->maploc();
    my $hsp  = shift;
    my $cx   = $hsp->cx_genome_part();
    my $id   = $cx->{id};
    if ($id =~ /^CDD:(\d+)$/) {
        push @{$cx->{links}}, [ "http://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid=$1", $id ];
    }
    my $obj = $ml->get_text( $id );
    $obj->read_tags();
    my $tags = $obj->all_tag_values();
    if (my $desc = $tags->{Description}) {
        $cx->{acc} = $id;
        push @{$tags->{Accession}}, $id;
        $id = $cx->{id} = $desc->[0];
    }
    my $col = $ml->pastel_text_color( $id );
    $cx->{fill}     = $col;
    $cx->{outline}  = $col;
    $cx->{hideName} = 1;
    delete $cx->{showDir};
    $ml->merge_cx_tags( $cx, $tags );
    return $cx;
}

sub cx_feature_parts {
    my $self = shift;
    my @cxs;
    foreach my $hsp ($self->anchored_features( @_ )) {
        push @cxs, $self->cx_for_feature( $hsp );
    }
    return @cxs;
}

sub cx_local_feature_parts {
    my $self = shift;
    $self->bench_start();
    unshift @_, '-aln' if ($#_ == -1);
    my $args = $self->parseparams( @_ );
    my (@cxs, @efParam);
    if (my $maxF = $args->{MAXFRAC}) {
        push @efParam, ( -maxfrac => $maxF );
    }
    if (my $filt = $args->{FEATFILTER} || $args->{FILTER}) {
        push @efParam, ( -filter => $filt );
    }
    foreach my $rng ($self->each_feature(@efParam)) {
        my $cx = $self->cx_for_feature( $rng );
        push @cxs, $cx;
    }
    $self->bench_end();
    return @cxs;
}

sub feature_text {
    my $self = shift;
    my $pad   = shift || "";
    my @feats = $self->each_feature();
    my $num   = $#feats + 1;
    return "" unless ($num);
    my $txt = sprintf("%s%d feature%s\n", $pad, $num, $num == 1 ? '' : 's');
    foreach my $rng (sort { $a->name() cmp $b->name() } @feats) {
        my @coords = $rng->coordinates();
        my $ft = join(',', map { $self->maploc->pretty_coordinates
                                     ($_->[0], $_->[1]) } @coords);
        $ft = "join($ft)" if ($#coords > 0);
        $txt .= sprintf("%s  %12s %s\n", $pad, $rng->name(), $ft);
    }
    return $txt;
}

