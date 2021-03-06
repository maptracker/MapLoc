# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::SnpTracker::MapLoc::IsolatedHSP;
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
use BMS::SnpTracker::MapLoc::TaggedHSP;
use vars qw(@ISA);

@ISA      = qw(BMS::SnpTracker::MapLoc::TaggedHSP );

sub obj_type { return "IsolatedHSP"; }

# ::IsolatedHSP
sub _validate_hsp_array {
    my $self = shift;
    $self->bench_start();
    my $nv = shift;
    return 0 if ($#{$nv} == -1);
    foreach my $h (@{$nv}) {
        # No sanity check, we do not expect the HSP lengths to be consistent
        $h->[3] = $h->[2] + ($h->[1] - $h->[0])
            unless (defined $h->[3]);
    }
    my $nos = $self->{NUM_OF_SEQS} = ($#{$nv->[0]} + 1) / 2;
    $self->death("HSPs must be defined as left/right pairs for 2 sequences",
                 "You have provided $nos sequences") unless ($nos == 2);
    $self->bench_end();
    return $#{$nv} + 1;
}
