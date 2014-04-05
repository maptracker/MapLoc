# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::SnpTracker::MapLoc::HSPed;
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

# Object which has HSPs assigned to it

use strict;
use BMS::SnpTracker::MapLoc::Common;
use BMS::SnpTracker::MapLoc::IsolatedHSP;
use BMS::SnpTracker::MapLoc::TaggedHSP;

use vars qw(@ISA);
@ISA      = qw(BMS::SnpTracker::MapLoc::Common);

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub obj_type { return "HSPs"; }

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

*str = \&strand;
*strands = \&strand;
sub strand {
    my $self = shift;
    if (defined $_[0]) {
        # Benchmark ~ 76.79us
        my $str = $_[0];
        if (my $r = ref($str)) {
            if ($r eq 'ARRAY') {
                # Being passed as an array
                my $errs = 0;
                map { $errs++ unless (/^([\+\-]?1|0)$/) } @{$str};
                $self->death
                    ("$errs Unrecognized strand entries in strand array",
                     join(',', map { "'$_'" } @{$str})) if ($errs);
                map { $_ += 0 } @{$str};
                $self->{STRAND} = $str;
            } else {
                $self->death("Unrecognized strand object '$str'");
            }
        } elsif ($str !~ /^([\+\-]?1|0)$/) {
            $self->death("Unrecognized strand '$str'");
        } else {
            # Single value. Assign to the second object
            $self->{STRAND} = [ 1, $str + 0 ];
        }
    }
    my $ns = $self->num_of_seqs();
    return wantarray ? @{$self->{STRAND} || [ map { 0 } (1..$ns)]} :
        $self->{STRAND} ? $self->{STRAND}[$ns == 1 ? 0 : 1] : 0;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub strand_for_seq {
    my $self = shift;
    my $ind  = $self->seqindex( shift );
    return 0 unless (defined $ind);
    my @strs = $self->strand();
    if ($#strs == 1) {
        return $strs[0] * $strs[1];
    } elsif ($#strs == 0) {
        # One sequence! 
        return $strs[0];
    }
    # More than one sequence. Return relative to the anchor
    my $aind = $self->anchor_index();
    return 0 unless (defined $aind);
    return $strs[$aind] * $strs[$ind];
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub hsps {
    my $self = shift;
    # Benchmark ~ 39.35us
    if (my $nv = shift) {
        $self->death("HSPs must be provided as 2D array ref")
            unless (ref($nv));

        # For an alignment where we have:
        # Seq1 :    1..2       9..10      50..80
        # Seq2 : 1001..1002 1039..1040  1050..1080
        # Seq3 :   91..92     87..88      15..45

        # Then the HSP array should be formatted as:
        # [ [ 1,2   1001,1002  91,92 ],
        #   [ 9,10  1039,1040  87,88 ],
        #   [ 50,80 1050,1080  15,45 ] ]

        # That is, each subarray is a column in the alignment
        my $num = $self->_validate_hsp_array( $nv, shift );
        if ($num) {
            $self->{HSPS} = $nv;
        } else {
            $self->{NUM_OF_SEQS} = 0;
            delete $self->{HSPS};
        }
        delete $self->{BOUNDS};
    }
    if (!$self->{HSPS} && $self->{PKEY}) {
        my $ml   = $self->maploc();
        my $read = $ml->{STH}{LOAD_HSPS} ||= $ml->dbh->prepare
            ( -name => "Recover HSPs from database",
              -sql  => "SELECT ptl1, ptr1, ptl2 FROM hsp WHERE aln_id = ? ORDER BY ptl1", );
        $read->execute($self->{PKEY});
        my $hsps = $read->fetchall_arrayref();
        if ($#{$hsps} == -1) {
            # No HSPs in database. We will set the array to be explicitly empty
            # The user can still override by passing a new value
            $self->{NUM_OF_SEQS} = 0;
            $self->{HSPS} = [];
        } else {
            map { $_->[3] = $_->[2] + $_->[1] - $_->[0] } @{$hsps};
            $self->{NUM_OF_SEQS} = 2;
            $self->{HSPS} = $hsps;
        }
    }
    return @{$self->{HSPS} || []};
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub length {
    my $self = shift;
    my $len = 0;
    # We will just measure the length of the first HSP:
    my @hsps = $self->hsps();
    map { $len += $_->[1] - $_->[0] - 1 } @hsps;
    return $len;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub score {
    my $self = shift;
    # This is not stored anywhere for ::HSPed - it is used transiently
    if (defined $_[0]) {
        $self->{SCORE} = $_[0];
    }
    return $self->{SCORE};
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub hsps_for_seq {
    my $self = shift;
    my $ind  = $self->seqindex( shift );
    return undef unless (defined $ind);
    $ind *= 2;
    my @rv;
    foreach my $hsp ($self->hsps()) {
        push @rv, [ $hsp->[$ind], $hsp->[$ind+1] ];
    }
    return wantarray ? @rv : \@rv;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub location_text_for_seq {
    my $self = shift;
    my $hsps = $self->hsps_for_seq( @_ );
    my $rv  = $self->strand_for_seq( @_ ) || 0;
    $rv .= ":".join(',', map { "$_->[0]..$_->[1]" } @{$hsps});
    return $rv;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub coordinates {
    my $self = shift;
    $self->bench_start();
    my @coords;
    foreach my $hsp ($self->hsps) {
        my @coord = @{$hsp};
        for my $i (0..$#coord) {
            if ($i % 2) {
                # 'Right' positions are decremented to be 'end'
                $coord[$i]--;
            } else {
                # 'Left' positions are incremented to be 'start'
                $coord[$i]++;
            }
        }
        # Note that this array represents a COLUMN in the alignment
        # not the data for a single sequence
        push @coords, \@coord;
    }
    $self->bench_end();
    return @coords;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub paired_coordinates_by_seq {
    my $self = shift;
    $self->bench_start();
    my @coords = $self->coordinates();
    my $ios    = $self->ind_of_seqs;
    # Each entry will ultimately be [ [s1,e1], [s2,e2], etc ]
    # for a single SEQUENCE
    my @paired = map { [] } (0..$ios);
    foreach my $hsp (@coords) {
        for my $i (0..$ios) {
            my $ind = $i * 2;
            push @{$paired[$i]}, [ $hsp->[$ind], $hsp->[$ind+1] ];
        }
    }
    $self->bench_end();
    return @paired;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub write_hsps {
    my $self = shift;
    return if ($self->{READ_ONLY});
    my $nos  = $self->num_of_seqs();
    return unless ($nos == 2);
    my $pkey = $self->pkey();
    return unless ($pkey);
    $self->bench_start();

    my $ml    = $self->maploc();
    my $dbh   = $ml->dbh();
    my $clearhsp = $ml->{STH}{DELETE_HSP_BY_ID} ||= $dbh->prepare
        ( -name => "Clear HSPs for a single alignment",
          -sql  => "DELETE FROM hsp WHERE aln_id = ?");
    my $addhsp = $ml->{STH}{ADD_HSP} ||= $dbh->prepare
        ( -name => "Add a single HSP",
          -sql  => "INSERT INTO hsp (aln_id, ptl1, ptr1, ptl2) VALUES (?,?,?,?)" );
    
    $clearhsp->execute($pkey);
    foreach my $hsp ($self->hsps()) {
        $addhsp->execute($pkey, $hsp->[0], $hsp->[1], $hsp->[2]);
    }
    $self->bench_end();
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub _validate_hsp_array {
    my $self = shift;
    # Benchmark ~ 50.35us
    my $nv = shift;
    return 0 if ($#{$nv} == -1);
    foreach my $h (@{$nv}) {
        my $len = $h->[1] - $h->[0];
        if (defined $h->[3]) {
            # ptr2 is defined. Sanity check it
            unless ($h->[3] - $h->[2] == $len) {
                my @seqs = $self->seqs();
                $self->err("CAUTION : HSP has inconsistent boundaries",
                           join("+", @seqs)." : ".join(",", @{$h}));
            }
        } else {
            # Automatically fill out ptr2:
            $h->[3] = $h->[2] + $len;
        }
    }
    my $nos = $self->{NUM_OF_SEQS} = ($#{$nv->[0]} + 1) / 2;
    unless ($nos == 2) {
        my $override = shift;
        $self->death("HSPs must be defined as left/right pairs for 2 sequences",
                     "You have provided $nos sequences") unless ($override);
    }
    return $#{$nv} + 1;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub num_of_seqs { return shift->{NUM_OF_SEQS} || 0; }

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub ind_of_seqs { return shift->num_of_seqs() - 1; }

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub bounds {
    my $self = shift;
    unless ($self->{BOUNDS}) {
        # Benchmark ~ 40.18us
        my @hsps = $self->hsps();
        if ($#hsps != -1) {
            my $nos  = $self->num_of_seqs();
            my @inds = (0..($nos - 1));
            my @bnds;
            foreach my $ind (@inds) {
                my $i1 = 2 * $ind;
                my $i2 = $i1 + 1;
                my @s  = sort { $a <=> $b } map { $_->[$i1], $_->[$i2] } @hsps;
                push @bnds, ($s[0], $s[-1]);
            }
            $self->{BOUNDS} = \@bnds;
        }
    }
    return @{$self->{BOUNDS} || []};
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub bounds_for_seq {
    my $self = shift;
    my $ind  = $self->seqindex( shift );
    return undef unless (defined $ind);
    my @bounds = $self->bounds();
    $ind *= 2;
    my @rv = ($bounds[$ind], $bounds[$ind+1]);
    return wantarray ? @rv : \@rv;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub pretty_bounds {
    my $self = shift;
    $self->bench_start();
    my $ml   = $self->maploc;
    my @bnds = $self->bounds();
    my @rv;
    for (my $ind = 0; $ind < $#bnds; $ind += 2) {
        push @rv, $ml->pretty_flanks($bnds[$ind], $bnds[$ind + 1]);
    }
    $self->bench_end();
    return @rv;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub seq_ids {
    my $self = shift;
    if (defined $_[0]) {
        # Benchmark ~ 20.42us
        if ($_[0]) {
            my @arr = ref($_[0]) ? @{$_[0]} : (@_);
            $self->death("You are defining ".scalar(@arr).
                         " sequences but $self->{SEQNUM} are required")
                if ($self->{SEQNUM} && $self->{SEQNUM} != ($#arr + 1));
            $self->{SEQIDS}   = [ @arr ];
        } else {
            # Blank entry (0 or "") clears sequences
            map { delete $self->{$_} } qw(SEQIDS SEQNAMES SEQ_LOOKUP);
        }
    }
    if (!$self->{SEQIDS} && $self->{SEQNAMES}) {
        # Benchmark ~ 31.71us
        # We can get the IDs from the names
        my $ml = $self->maploc;
        $self->{SEQIDS} = [ map { $ml->text_to_pkey( $_ ) 
                                  } @{$self->{SEQNAMES}} ];
    }
    return @{$self->{SEQIDS} || [map { 0 } (1..$self->num_of_seqs())]};
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub seqs {
    my $self = shift;
    if (defined $_[0]) {
        $self->bench_start('set');
        if ($_[0]) {
            $self->death("You are defining ".scalar(@_).
                         " sequences but $self->{SEQNUM} are required")
                if ($self->{SEQNUM} && $self->{SEQNUM} != ($#_ + 1));
            $self->{SEQNAMES}   = [ @_ ];
        } else {
            # Blank entry (0 or "") clears sequences
            map { delete $self->{$_} } qw(SEQIDS SEQNAMES SEQ_LOOKUP);
        }
        $self->bench_end('set');
    }
    if (!$self->{SEQNAMES} && $self->{SEQIDS}) {
        # Benchmark ~ 51.52us
        # We can get the names from the IDs
        my $ml = $self->maploc;
        my $ns = $self->{SEQNAMES} = 
            [ map { $ml->cached_pkey_to_text( $_ ) } @{$self->{SEQIDS}} ];
        my @errs;
        for my $ni (0..$#{$ns}) {
            push @errs, sprintf("seq_id %d '%s' does not map to text",
                                $ni, $self->{SEQIDS}[$ni]) unless ($ns->[$ni]);
        }
        $self->death(@errs) unless ($#errs == -1);
    }
    return @{$self->{SEQNAMES} || [map { "" } (1..$self->num_of_seqs())]};
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub rename_seq {
    my $self = shift;
    my ($old, $new) = @_;
    return 0 unless ($new);
    my $ind = $self->seqindex( $old );
    unless (defined $ind) {
        $self->err("Could not find the sequence request '$old'");
        return 0;
    }
    if ($self->{SEQNAMES}) {
        $self->{SEQNAMES}[$ind] = $new;
    }
    if ($self->{SEQIDS}) {
        $self->{SEQIDS}[$ind] = $self->maploc->text_to_pkey( $new );
    }
    return $new;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

*other_seq_id = \&other_seq_ids;
sub other_seq_ids {
    my $self = shift;
    my $sid  = shift || 0;
    my @rv;
    foreach my $osid ($self->seq_ids) {
        push @rv, $osid unless ($osid == $sid);
    }
    return wantarray ? @rv : $rv[0];
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

*other_seq = \&other_seqs;
sub other_seqs {
    my $self = shift;
    my $ml   = $self->maploc;
    my @sids = $self->other_seq_ids( $ml->text_to_pkey( shift ) );
    my @rv   = map { $ml->pkey_to_text( $_ ) } @sids;
    return wantarray ? @rv : $rv[0];
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub chr_seq_info {
    my $self = shift;
    unless ($self->{CHR_ACC_INFO}) {
        $self->bench_start();
        my $ml   = $self->maploc();
        my $look = $ml->{STH}{CHECK_IF_ACC_IS_CHR} ||= $ml->dbh->prepare
            ( -name => "Check accession ID against accession table",
              -sql  => "SELECT chr, build FROM loc_to_acc WHERE acc_id = ?" );
        my @csids;
        my @sids = $self->seq_ids();
        my $targ = $self->{CHR_ACC_INFO} = [];
        for my $i (0..$#sids) {
            $look->execute( $sids[$i] );
            my $rows = $look->fetchall_arrayref();
            next if ($#{$rows} != 0);
            # ACC_ID, Index, ChrName, Build
            # warn "$sids[$i], $i, @{$rows->[0]}";
            push @{$targ}, [ $sids[$i], $i, @{$rows->[0]} ];
        }
        $self->bench_end();
    }
    my @rv = @{$self->{CHR_ACC_INFO}};
    if (my $col = $_[0]) {
        # Passing a value requests a specific column
        $col--; # Turn 'column' into an array index
        my %vals = map { $_->[$col] => undef } @rv;
        @rv = sort keys %vals;
    }
    return wantarray ? @rv : \@rv;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub chr_hsps_for_alignment {
    my $self = shift;
    unless ($self->{CHR_HSPS_FOR_ALN}) {
        my $dat = $self->{CHR_HSPS_FOR_ALN} = [];
        my @csi = $self->chr_seq_info();
        if ($#csi == 0) {
            my ($accid, $ind, $chr, $build) = @{$csi[0]};
            my $hsps = [ sort { $a->[0] <=> $b->[0] ||
                                    $a->[1] <=> $b->[1] }
                         $self->hsps_for_seq( $ind ) ];
            push @{$dat}, ($hsps, $build, $chr);
        }
    }
    return @{$self->{CHR_HSPS_FOR_ALN}};
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

# Note that this is NOT left / right
sub chr_start_end {
    my $self = shift;
    unless ($self->{CHR_START_END}) {
        my ($hsps) = $self->chr_hsps_for_alignment();
        $self->{CHR_START_END} = $hsps ? 
            [ $hsps->[0][0] + 1, $hsps->[-1][1] - 1] : [0,0];
    }
    return wantarray ? @{$self->{CHR_START_END}} : $self->{CHR_START_END};
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub chr_start {
    my $cse  = shift->chr_start_end();
    return $cse->[0];
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub chr_end {
    my $cse  = shift->chr_start_end();
    return $cse->[1];
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub chr_width {
    my $cse  = shift->chr_start_end();
    return $cse->[0] || $cse->[1] ? $cse->[1] - $cse->[0] + 1 : 0;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub chr_location_tag {
    my $self = shift;
    my ($hsps, $build, $chr) = $self->chr_hsps_for_alignment();
    return "" unless ($chr && $build);
    my $rv = "$chr.$build";
    return $rv if (!$hsps || $#{$hsps} == -1);
    $rv .= join('..', $hsps->[0][0], $hsps->[-1][1]);
    return $rv;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

# Simple text representation:
sub chr_footprint {
    my $self = shift;
    my @rv;
    my $sc = $self->score();
    foreach my $ci ($self->chr_seq_info()) {
        my ($s, $e) = $self->bounds_for_seq( $ci->[1] );
        my $pos = $self->maploc->pretty_flanks($s, $e);
        push @rv, sprintf("%s:%s [%.1f%%]", $ci->[2], $pos, $sc);
    }
    return wantarray ? @rv : $rv[0];
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub chr_acc {
    my $self = shift;
    my @rv = map { $self->maploc->pkey_to_text( $_ ) || ""} $self->chr_acc_ids(@_);
    return wantarray ? @rv : $rv[0];
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub chr_acc_ids {
    my $self = shift;
    my @rv = @{$self->{CHR_ACCIDS} ||= $self->chr_seq_info( 1 )};
    return wantarray ? @rv : $rv[0];
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

*chr_ind = \&chr_inds;
sub chr_inds {
    my $self = shift;
    my @rv = @{$self->{CHR_INDS} ||= $self->chr_seq_info( 2 )};
    return wantarray ? @rv : $rv[0];
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

*chr = \&chr_names;
*chr_name = \&chr_names;
sub chr_names {
    my $self = shift;
    my @rv = @{$self->{CHR_NAMES} ||= $self->chr_seq_info( 3 )};
    return wantarray ? @rv : $rv[0];
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

*build = \&chr_builds;
*chr_build = \&chr_builds;
sub chr_builds {
    my $self = shift;
    my @rv = @{$self->{CHR_BUILDS} ||= $self->chr_seq_info( 4 )};
    return wantarray ? @rv : $rv[0];
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

*seq_index = \&seqindex;
sub seqindex {
    my $self = shift;
    my $req  = shift;
    return undef unless (defined $req);
    # Assume small numbers are already indices. This means we can not use
    # this method to reliably recover info for a sequence with acc_id < 5
    return $req if ($req =~ /^\d+$/ && $req < 5);
    my $slu = $self->_seq_lookup();
    return $slu->{uc($req)};
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub _seq_lookup {
    my $self = shift;
    unless ($self->{SEQ_LOOKUP}) {
        # Benchmark ~ 74.16us
        my @seqs = $self->seqs();
        if ($seqs[0]) {
            # names have been set
            my $sl = $self->{SEQ_LOOKUP} =
            { map { uc($seqs[$_]) => $_ } (0..$#seqs) };
            # Also map the IDs into the lookup. There is a teeny chance
            # that one of the IDs could also be a numeric name. Highly unlikely
            my @sids = $self->seq_ids();
            map { $sl->{$sids[$_]} = $_ } (0..$#sids);
        }
    }
    return $self->{SEQ_LOOKUP} || {};
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

*anchor_index = \&set_anchor;
sub set_anchor {
    my $self = shift;
    if (my $nv = shift) {
        my $ind = $self->seqindex($nv);
        if (defined $ind) {
            $self->{ANCHOR_IND} = $ind;
        } else {
            $self->err("Can not set anchor to unknown member '$nv'");
        }
    }
    return $self->{ANCHOR_IND};
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub anchor_to_genome {
    my $self = shift;
    my @aids = $self->chr_acc_ids();
    if ($#aids == 0) {
        return $self->set_anchor( $aids[0] );
    }
    return 0;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub requested_anchor {
    my $self = shift;
    my $aInd = $self->anchor_index();
    if (my $ancReq = shift) {
        $aInd = $self->seqindex($ancReq);
    }
    return $aInd;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub anchor_id {
    my $self = shift;
    return ($self->{SEQIDS} && defined $self->{ANCHOR_IND}) ?
        $self->{SEQIDS}[ $self->{ANCHOR_IND} ] : undef;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub left_anchor {
    my $self = shift;
    my $ind = $self->{ANCHOR_IND};
    return undef unless (defined $ind);
    my @bnds = $self->bounds();
    return $bnds[$ind * 2];
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub right_anchor {
    my $self = shift;
    my $ind = $self->{ANCHOR_IND};
    return undef unless (defined $ind);
    my @bnds = $self->bounds();
    return $bnds[1 + $ind * 2];
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub intersect {
    my $self = shift;
    unshift @_, '-hsp' if ($#_ == 0);
    my $args = $self->parseparams( @_ );
    my $aid  = $args->{ANCHOR} || $self->anchor_id();
    my @HSPs = $self->_prepare_overlap( $args->{HSP}, $aid );
    my $extr = pop @HSPs;
    # Array context means the user wants information on the distribution
    # (spacing) of the non-overlapping HSPs:
    return wantarray ? (undef, undef) : undef unless ($extr);
    $self->bench_start();
    my ($aHsps, $bHsps) = @HSPs;
    my @ancInds       = map { $_->[1] } @{$extr}; # SeqIndex of Anchor
    my @ancHspInd     = map { 2 * $_ } @ancInds; # HSP Index of Anchor
    my ($aHin, $bHin) = @ancHspInd;
    my $ml            = $self->maploc();
    # Idx = the starting index value for the HSP sets
    my ($aIdx, $bIdx) = (0,0);
    # Halt = the index value for the last entry
    my ($aHalt, $bHalt) = ($#{$aHsps}, $#{$bHsps});
    # _prepare_overlap() has normalzied the anchor HSPs in ascending order
    # We do not need to change increment or starting position
    my @isolated;
    # Need better parameter names for this :/
    my $isoA = $isolated[0] = [] if ($args->{ISOA} || $args->{ISOINTA});
    my $isoB = $isolated[1] = [] if ($args->{ISOB} || $args->{ISOINTB});
    # Isolated data are structured in the following format:
    # [ [ Anchor Left, Anchor Right],
    #   [ Bases between this HSP and other HSP on left, other HSP index ],
    #   [ Bases between this HSP and other HSP on left, other HSP index ] ]
    # So a "bases between" value of 0 means the HSPs are immediately adjacent

    my (@errs, @intHSPs, @isoHSPs);
    my $iLoop  = 0;
    my $intCnt = 0;
    # warn "Halts = ($aHalt, $bHalt)\n";
    while ($aIdx <= $aHalt && $bIdx <= $bHalt) {
        if (++$iLoop > 1000) {
            $self->err("Probable infinite loop ($iLoop iterations) finding intersections");
            last;
        }
        # Get the HSPs being considered:
        my $aHsp = $aHsps->[$aIdx];
        my $bHsp = $bHsps->[$bIdx];
        # Get the left and right flank positions:
        my ($aL, $aR) = ($aHsp->[$aHin], $aHsp->[$aHin + 1]);
        my ($bL, $bR) = ($bHsp->[$bHin], $bHsp->[$bHin + 1]);

        # Note that RightA == LeftB is NOT an overlap
        # Example:  positionA = 10 and positionB = 12
        # RightA = 11, LeftB = 11, there is a single position separating them

        # Also, RightA == LeftB + 1 is an adjacency, not an overlap
        if ($aR <= $bL + 1) {
            # $self->msg("[DEBUG]","$aIdx <= $aHalt ($aL, $aR)","$bIdx <= $bHalt ($bL, $bR)");
            # The A HSP is non-overlapping, and to the left. Increment A
            if ($isoA) {
                # Record the isolated HSP and distance to B
                $self->_extract_isolated_hsp
                    ( $aL, $aR, $aid, $isoA, $extr, \@HSPs, 0,
                      $aIdx, $bIdx, $aHin, $bHin );
            }
            $aIdx++;
            next;
        }
        if ($bR <= $aL + 1) {
            # The B HSP is non-overlapping, and to the left. Increment B
            if ($isoB) {
                # Record the isolated HSP and distance to A
                $self->_extract_isolated_hsp
                    ( $bL, $bR, $aid, $isoB, $extr, \@HSPs, 1,
                      $bIdx, $aIdx, $bHin, $aHin );
            }
            $bIdx++;
            next;
        }
        # At this point we should have A and B HSPs that are overlapping
        # We need to find by how much, and slice through all the HSPs to
        # capture the common overlap.
        
        # Take the 'minimum' coordinates as the anchor intersection:
        my $ancL  = $aL < $bL ? $bL : $aL;
        my $ancR  = $aR > $bR ? $bR : $aR;
        my $ancW  = $ancR - $ancL;
        my @hinds = ($aIdx, $bIdx);
        $intCnt++;
        my %newHSP = ( $aid => [$ancL, $ancR] );
        for my $h (0..1) {
            my $sim = $extr->[$h][0];
            next unless ($sim); # Not all HSPs will have a non-anchor object
            my $hsp = $HSPs[$h][ $hinds[$h] ];
            while (my ($hi, $sd) = each %{$sim}) {
                # $hi is the HSP index (x2) for this non-anchor sequence
                # $sid is the sequence ID for the non-anchor seq
                # $str is the strand relative to the anchor
                my ($sid, $str) = @{$sd};
                my $l = $hsp->[$hi];
                if ($str < 0) {
                    # -1 strand orientation - offset bounds from right
                    $l += ($hsp->[$ancHspInd[$h]+1] - $ancR);
                } else {
                    # +1 strand orientation - offset bounds from left
                    $l += ($ancL - $hsp->[$ancHspInd[$h]]);
                }
                my $r = $l + $ancW;
                # These are the coordinates in anchor space:
                my @coord = ($l,$r);
                if (my $prior = $newHSP{$sid}) {
                    # A sequence is present twice. This could happen if it
                    # is represented on both HSPs. We want to verify that
                    # the coordinates are consistent
                    unless ($prior->[0] == $l && $prior->[1] == $r) {
                        my $sn =  $ml->pkey_to_text( $sid );
                        push @errs, "New HSP # $intCnt [$aIdx vs $bIdx $ancL..$ancR] has inconsistent coordinates for $sn: [$prior->[0], $prior->[1]] vs [$l,$r]";
                    }
                } else {
                    # Almost all the time we expect to be here
                    $newHSP{$sid} = [$l, $r];
                }
            }
        }
        push @intHSPs, \%newHSP;
        # At least one of the HSPs should be incremented, but it is possible
        # that incrementing both is not appropriate
        if ($aIdx < $aHalt && $aHsps->[$aIdx + 1][$aHin] < $bR) {
            # We can increment A and the next A HSP still overlaps the B one
            #$self->msg("[DEBUG]","Increment A $aIdx < $aHalt");
            $aIdx++;
        } elsif ($bIdx < $bHalt && $bHsps->[$bIdx + 1][$bHin] < $aR) {
            # We can increment B and the next B HSP still overlaps the A one
            #$self->msg("[DEBUG]","Increment B $bIdx < $bHalt");
            $bIdx++;
        } else {
            # If we can not find another overlap with either of the A or B
            # HSPs, increment both
            #$self->msg("[DEBUG]","Increment Both");
            $aIdx++;
            $bIdx++;
        }
    }
    # At this point, we should have either exhausted both A and B, or have one
    # or more HSPs left in only /one/ of them
    if ($aIdx <= $aHalt) {
        # If we have any A left, it means that B is all gone
        if ($isoA) {
            # See if there are any isolated A HSPs left on the right:
            while ($aIdx <= $aHalt) {
                my $aHsp = $aHsps->[$aIdx];
                my ($aL, $aR) = ($aHsp->[$aHin], $aHsp->[$aHin + 1]);
                $self->_extract_isolated_hsp
                    ( $aL, $aR, $aid, $isoA, $extr, \@HSPs, 0,
                      $aIdx, $bHalt+1, $aHin, $bHin );
                $aIdx++;
            }
        }
    } elsif ($bIdx <= $bHalt) {
        # If we have any B left, it means that A is all gone
        if ($isoB) {
            # See if there are any isolated B HSPs left on the right:
            while ($bIdx <= $bHalt) {
                my $bHsp = $bHsps->[$bIdx];
                my ($bL, $bR) = ($bHsp->[$bHin], $bHsp->[$bHin + 1]);
                # We know that the only B neighbor is the last one, on left
                $self->_extract_isolated_hsp
                    ( $bL, $bR, $aid, $isoB, $extr, \@HSPs, 1,
                      $bIdx, $aHalt+1, $bHin, $aHin );
                $bIdx++;
            }
        }
    }

    # Validate strands
    my %memStr = ( $aid => {1 => 1 } );
    foreach my $exd (@{$extr}) {
        while (my ($hi, $sd) = each %{$exd->[0]}) {
            my ($sid, $str) = @{$sd};
            $memStr{$sid}{$str} = 1;
        }
    }
    my @sids = sort { $a <=> $b } keys %memStr;

    # Check that each sequence has a unique strand:
    while (my ($sid, $strH) = each %memStr) {
        # Do not consider undefined strands:
        delete $strH->{0};
        my @strs = keys %{$strH};
        if ($#strs != 0) {
            # This sequence is pointing in both directions!
            my $sn =  $ml->pkey_to_text( $sid );
            push @errs, "Sequence $sn is represented in BOTH strands!";
        }
        $memStr{$sid} = $strs[0] || 0;
    }

    # Now we need to capture the sequence IDs and
    # transform @intHSP into HSP arrays
    my @finalHSPs = $self->_sid_hash_to_hsp_array( \@intHSPs, \@sids );

    my $isIsolated;
    if ($#finalHSPs == -1) {
        # my $debug = join("\n","\n#######################\nNon-Intersection Report (H1, H2, Isolated):\n#######################", $self->hsp_text(), $args->{HSP}->hsp_text()); warn $debug . $self->branch(\@isolated); # if ($#isolated == -1);#  if ($debug =~ /rs113708176/);
        if (($args->{ISOINTA} && $#{$isoA} != -1) || 
            ($args->{ISOINTB} && $#{$isoB} != -1)) {
            # Request to turn isolated (non-overlapping) intersections into
            # an HSP
            my $iso = $isoA && $#{$isoA} != -1 ? $isoA : $isoB;
            my @isoInts = map { $_->[3] } @{$iso};
            @finalHSPs = $self->_sid_hash_to_hsp_array( \@isoInts, \@sids );
            if ($#finalHSPs == -1) {
                $self->bench_end();
                return wantarray ? (undef, @isolated) : undef;
            }
            $isIsolated = 1;
        } else {
            $self->bench_end();
            return wantarray ? (undef, @isolated) : undef;
        }
    }
    my $final = $isIsolated ? 
        BMS::SnpTracker::MapLoc::IsolatedHSP->new( $ml ) :
        BMS::SnpTracker::MapLoc::TaggedHSP->new( $ml );
    $final->strand( [ map { $memStr{$_} } @sids ] );
    $final->hsps( \@finalHSPs, 1 );
    $final->seq_ids( \@sids );
    $final->_seq_lookup();
    $final->set_anchor( $aid );
    $self->bench_end();
    # my $debug = join("\n","Intersection Report (H1, H2, Int):", $self->hsp_text(), $args->{HSP}->hsp_text(), $final->hsp_text()); warn $debug . $self->branch(\@HSPs)  if ($debug =~ /\Q72509823\E/);
    return wantarray ? ($final, @isolated) : $final;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub _sid_hash_to_hsp_array {
    my $self = shift;
    my ($intHashes, $sids) = @_;
    return () if ($#{$intHashes} == -1);
    my @finalHSPs;
    my @errs;
    for my $iH (0..$#{$intHashes}) {
        my $newHSP = $intHashes->[$iH];
        my @col;
        foreach my $sid (@{$sids}) {
            if (my $lr = $newHSP->{$sid}) {
                push @col, @{$lr};
            } else {
                my $sn = $self->maploc->pkey_to_text( $sid );
                push @errs, "Sequence $sn is not represented in HSP #".
                    ($iH+1);
            }
        }
        push @finalHSPs, \@col;
    }
    if ($#errs != -1) {
        $self->err("HSP generation errors:", @errs);
        return ();
    }
    return @finalHSPs;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub _extract_isolated_hsp {
    my $self = shift;
    my ( $ancL, $ancR, $aid, $priIsoArray, $extr, $HSPs, $priInd,
         $priIdx, $secIdx, $priHin, $secHin ) = @_;
    my $secInd  = $priInd ? 0 : 1;
    my $secHsps = $HSPs->[$secInd];
    my @iso     = ([$ancL, $ancR]);
    if ($secIdx) {
        # The secondary alignment has an HSP to the left
        my $neigh = $secIdx-1;
        $iso[1] = [ $ancL - $secHsps->[$neigh][$secHin + 1] + 1, $neigh];
    }
    unless ($secIdx > $#{$secHsps}) {
        # The secondary alignment has an HSP to the right
        $iso[2] = [$secHsps->[$secIdx][$secHin] - $ancR + 1, $secIdx ];
    }
    push @{$priIsoArray}, \@iso;

    # Record details of
    my $isoHSP = $iso[3] = {
        $aid => [$ancL, $ancR],
    };
    for my $h (0..1) {
        my $sim    = $extr->[$h][0];
        next unless ($sim); # Not all HSPs will have a non-anchor object
        my $hsps   = $HSPs->[$h]; # [ $hinds[$h] ];
        my $isMain = $h == $priInd ? 1 : 0;
        while (my ($hi, $sd) = each %{$sim}) {
            # $hi is the HSP index (x2) for this non-anchor sequence
            # $sid is the sequence ID for the non-anchor seq
            # $str is the strand relative to the anchor
            my ($sid, $str) = @{$sd};
            if ($isMain) {
                # This is the isolated sequence. We will just take the
                # flanks as-is
                $isoHSP->{$sid} = [ $hsps->[ $priIdx ][$hi], 
                                    $hsps->[ $priIdx ][$hi + 1] ];
            } else {
                # This is the 'other' sequence. We will represent it as a gap
                # position
                my $oH    = $isoHSP->{$sid} = [];
                my ($lftMod, $rgtMod, $strMod) = ($str < 0) ?
                    (1, 0, -1) : (0, 1, 1);
                if (!$secIdx) {
                    # We are to the "far left of the other sequence
                    # Key off the 0 HSP:
                    my $refPos = $hsps->[0][$hi + $lftMod];
                    # +1 
                    $oH->[ $lftMod ] = $refPos;
                    $oH->[ $rgtMod ] = $refPos + $strMod;
                } elsif ($secIdx > $#{$hsps}) {
                    # We are to the "far right of the other sequence
                    my $refPos = $hsps->[$#{$hsps}][$hi + $rgtMod];
                    # +1 
                    $oH->[ $lftMod ] = $refPos;
                    $oH->[ $rgtMod ] = $refPos + $strMod;
                } else {
                    # We are between two of the other HSPs
                    $oH->[ $lftMod ] = $hsps->[$secIdx-1][$hi + $rgtMod] - $strMod;
                    $oH->[ $rgtMod ] = $hsps->[$secIdx][$hi + $lftMod] + $strMod;
                }
            }
        }
    }
    # push @isoHSPs, $isoHSP;
    
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub _prepare_overlap {
    my $self = shift;
    my $other = shift;
    return undef unless ($other);
    # Benchmark ~ 113.60us
    $self->bench_start();
    $self->death("intersect() must pass another HSP object, not '$other'")
        unless ($other->isa('BMS::SnpTracker::MapLoc::HSPed'));
    my ($sAid, $oAid);
    if (my $explAnch = shift) {
        $sAid = $oAid = $explAnch;
    } else {
        $sAid ||= $self->anchor_id  || 0;
        $oAid ||= $other->anchor_id || 0;
    }
    unless ($sAid && $oAid) {
        $self->err("Can not perform intersection unless both HSPs are anchored");
        $self->bench_end();
        return undef;
    }
    unless($sAid == $oAid) {
        $self->err("Can not perform intersection unless both HSPs are anchored to the same object");
        $self->bench_end();
        return undef;
    }
    # Get the anchor index in each HSP set:
    my ($sInd, $oInd) = map { $_->seqindex($sAid) } ($self, $other);
    unless(defined $sInd && defined $oInd) {
        $self->err("One or both anchor indices are not defined - can not perform intersection");
        $self->bench_end();
        return undef;
    }

    my ($sHi,  $oHi)  = map { $_ * 2 } ($sInd, $oInd);
    # Recover HSPs, sorted ascending by the anchor coordinates
    my @sHSPs = sort { $a->[$sHi] <=> $b->[$sHi] ||
                           $a->[$sHi+1] <=> $b->[$sHi+1] } $self->hsps();
    my @oHSPs = sort { $a->[$oHi] <=> $b->[$oHi] ||
                           $a->[$oHi+1] <=> $b->[$oHi+1]} $other->hsps();
    
    # Temporarily piggyback object info with sorted HSP data:
    push @sHSPs, [$self,  $sInd];
    push @oHSPs, [$other, $oInd];
    my @rv = (\@sHSPs, \@oHSPs);
    
    # I do not think that these need to be ordered:
    #
    # Order the data such that the one with the furthest 'left' HSP is first:
    #my @rv = $sHSPs[0][$sHi] > $oHSPs[0][$oHi] ?
    #    (\@oHSPs, \@sHSPs) :(\@sHSPs, \@oHSPs);

    # Organize extraction information
    my @extractor;
    for my $h (0..1) {
        my $oDat = pop @{$rv[$h]}; # take off object info from HSP
        my ($obj, $aInd) = @{$oDat};
        my @str          = $obj->strand();  # The strands for the objects
        my @seqs         = $obj->seq_ids(); # seq_ids of the objects
        my $anStr        = $str[$aInd];     # The strand of the anchor
        my $seqIdMap;
        for my $n (0..$#str) {
            my $sid = $seqs[$n];
            next if ($sid == $sAid); # Ignore the anchor
            # Find the relative strand between anchor and this sequence:
            my $str = $str[$n] * $anStr;
            $seqIdMap ||= {};
            $seqIdMap->{($n * 2)} = [ $sid, $str ];
        }
        $extractor[$h] = [$seqIdMap, $aInd, $obj ];
    }
    push @rv, \@extractor;
    $self->bench_end();
    return @rv;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

*full_text = \&hsp_text;
*to_text = \&hsp_text;
sub hsp_text {
    my $self = shift;
    $self->bench_start();
    my @cols;
    my $ml   = $self->maploc();
    my @seqs = $self->seqs();
    my $aind = $self->anchor_index();
    if (defined $aind) {
        $seqs[$aind] .= ' /A\\';
    }
    push @cols, ['Sequence', @seqs ];
    my @strs = $self->strands();
    push @cols, [ 'Str', map { !$_ ? '' : $_ < 0 ? '-1' : '+1'} @strs ];
    my @hsps = $self->hsps();
    for my $h (0..$#hsps) {
        my $hsp = $hsps[$h];
        my @hcol = ();
        for my $i (0..$#seqs) {
            my $ind = $i * 2;
            push @hcol, $ml->pretty_flanks( $hsp->[$ind], $hsp->[$ind+1]);
        }
        push @cols, [ "HSP-".($h+1), @hcol];
    }
    my @colsets;
    for my $cn (0..$#cols) {
        if ($#colsets == -1 || $colsets[-1]{t} > 100) {
            push @colsets,{
                w => [],
                t => 0,
                c => [],
            };
        }
        my ($w) = sort { $b <=> $a } map { CORE::length($_) } @{$cols[$cn]};
        my $dat = $colsets[-1];
        push @{$dat->{w}}, $w;
        $dat->{t} += $w + 2;
        push @{$dat->{c}}, $cn;
    }
    my $txt = "";
    my @bits;
    foreach my $meth ('obj_type', 'pkey', 'name', 'source', 'score', 'howbad', 'length') {
        if (my $cb = $self->can($meth)) {
            my $t = &{$cb}($self);
            if (defined $t && $t ne '') {
                push @bits, sprintf("%s(%s)", $meth, $t);
            }
        }
    }
    $txt .= join(' ', @bits)."\n" unless ($#bits == -1);
    if ($self->can('all_tag_values')) {
        my $tvs = $self->all_tag_values();
        foreach my $tag (sort keys %{$tvs}) {
            $txt .= sprintf("/%s=%s\n", $tag, join(' && ', @{$tvs->{$tag}}));
        }
    }
    my $rn = $#seqs + 1;
    for my $cs (0..$#colsets) {
        my $dat     = $colsets[$cs];
        my @blkCols = @{$dat->{c}};
        my @widths  = @{$dat->{w}};
        my $frm = join("|", map { "%${_}s" } @widths);
        $frm = "> $frm" if ($cs);
        $frm .= " >" if ($cs != $#colsets);
        $frm .= "\n";
        my $bar = sprintf($frm, map { '-' x $_ } @widths);
        $txt .= $bar if ($cs);
        for my $r (0..$rn) {
            my @row;
            for my $ci (0..$#blkCols) {
                my $c = $blkCols[$ci];
                my $v = $cols[$c][$r];
                if ($c > 1) {
                    # HSP column - center it
                    my $pad = int(($widths[$ci] - CORE::length($v))/2);
                    $v .= " " x $pad if ($pad);
                }
                push @row, $v;
            }
            $txt .= sprintf($frm, @row);
            $txt .= $bar unless ($r);
        }
    }
    $self->bench_end();
    return $txt;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub ordered_coordinates {
    my $self = shift;
    my $aInd = $self->requested_anchor( @_ );
    return () unless (defined $aInd);
    my @crds = $self->paired_coordinates_by_seq();
    my @inds = ($aInd);
    my @strs = $self->strands();
    map { push @inds, $_ unless ($_ == $aInd) } (0..$#strs);
    
    # Extract and fully de-reference
    my @Acrd; #  = map { [@{$_}] } @{$crds[$aInd]};
    foreach my $ind (@inds) {
        my @crd = @{$crds[$ind]};
        map { push @{$Acrd[$_]}, @{$crd[$_]} } (0..$#crd);
    }
    
    # Also, Isaac wants the HSPs to always be in ascending order:
    @Acrd    = sort { $a->[0] <=> $b->[0] } @Acrd
        if ($Acrd[0][0] > $Acrd[-1][0]);
    return ( \@Acrd, \@inds, \@strs );
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub cx_genome_part {
    my $self = shift;
    my ($Acrd, $inds, $strs) = $self->ordered_coordinates( @_ );
    return wantarray ? () : undef unless ($Acrd);
    $self->bench_start();
    my $aInd = $inds->[0];
    my $aStr = $strs->[$aInd];
    my @seqs = $self->seqs();
    my $len  = $self->length();
    my @common = ( outline => 'rgb(0,0,0)',
                   len     => $len,
                   anchor  => $seqs[$aInd],
                   members => [ map { $seqs[$_] } @{$inds} ],
                   data    => $Acrd,
                   showDir => 1, );
    push @common, ( isGap => 1 ) if (!$len && defined $len);
    if ($self->can('all_tag_values')) {
        my $tags = $self->all_tag_values();
        push @common, ( tags => $tags, len => $self->length(), );
    }

    my @parts;
    for my $i (0..$self->ind_of_seqs) {
        # Get each non-anchor sequence
        next if ($i == $aInd);
        my $str = ($aStr * $strs->[$i]) || 0;
        push @parts, {
            id  => $seqs[$i],
            dir => $str < 0 ? 'left' : 'right',
            @common,
        };
    }
    $self->bench_end();
    return wantarray ? @parts : $parts[0];
}
