# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::SnpTracker::MapLoc::Range;
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
use BMS::SnpTracker::MapLoc::Transient;
use BMS::SnpTracker::MapLoc::HSPed;
use vars qw(@ISA);

@ISA      = qw(BMS::SnpTracker::MapLoc::Transient
               BMS::SnpTracker::MapLoc::HSPed );

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my ($ml) = @_;
    my $self = {
    };
    $self->{MAPLOC} = $ml;
    bless ($self, $class);
    return $self;
}

sub obj_type { return "Range"; }

# ::Range
sub _validate_hsp_array {
    my $self = shift;
    # Benchmark ~ 25.99us
    my $nv = shift;
    my %w; map { $w{$#{$_}}++ } @{$nv};
    my @seen = keys %w;
    if ($#seen == -1) {
        # Empty array
        $self->{NUM_OF_SEQS} = 0;
        $self->bench_end();
        return 0;
    } elsif ($#seen != 0) {
        $self->err("Inconsistent HSPs provided for Range.",
                   "A Range entry should have only a single sequence represented by pairs of left/right boundaries",
                   join(" + ", map { "[".join(',', @{$_})."]"} @{$nv}));
        $self->bench_end();
        return 0;
    }
    my $bndCnt = $seen[0] + 1;
    unless ($bndCnt == 2) {
        $self->err("Incorrect number of HSP bounds provided ($bndCnt).",
                   "A Range entry should have only a single sequence represented by pairs of left/right boundaries",
                   join(" + ", map { "[".join(',', @{$_})."]"} @{$nv}));
        $self->bench_end();
        return 0;
    }
    
    $self->{NUM_OF_SEQS} = $bndCnt / 2;
    return $#{$nv} + 1;
}


# Strand is being handled oddly (ie proabably poorly)
# The default for ::HSPed is to always return +1 when requested in scalar
# context and when a single entity is present. We want some ranges
# like (features) to be able to occupy the negative strand
*str = \&strand;
*strands = \&strand;
sub strand {
    my $self = shift;
    if (defined $_[0]) {
        # Benchmark ~ 33.21us
        my $str = $_[0];
        if (my $r = ref($str)) {
            $str = $str->[0];
        }
        if ($str !~ /^([\+\-]?1|0)$/) {
            $self->death("Unrecognized strand '$str'");
        } else {
            $self->{STRAND} = [ $str + 0 ];
        }
    }
    my ($rv) = @{$self->{STRAND}};
    return wantarray ? ( $rv ) :  $rv;
}

# ::Range
sub cx_genome_part {
    my $self = shift;
    $self->bench_start();
    my @hsps = $self->coordinates();
    my @seqs = $self->seqs();
    my $len  = $self->length();
    my $str  = $self->strand() || 0;
    my $anc  = $seqs[0];
    my $id   = $self->can('name') ? $self->name() : $self->obj_type();
    my $cx = {
        id      => $id,
        len     => $len,
        data    => \@hsps,
        members => [$id],
        anchor  => $anc,
        outline => 'rgb(0,0,0)',
        dir     => $str < 0 ? 'left' : 'right',
    };
    $cx->{isGap} = 1 if (!$len && defined $len);
    if ($self->can('all_tag_values')) {
        my $tags = $self->all_tag_values();
        $cx->{tags} = $tags;
    }
    if ($self->can('source')) {
        if (my $src = $self->source()) {
            $cx->{tags}{Source} = [ $src ];
        }
    }
    $cx->{tags}{Via} = [ $anc ];
    $self->bench_end();
    return wantarray ? ($cx) : $cx;
}
