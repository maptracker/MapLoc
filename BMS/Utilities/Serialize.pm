# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::Utilities::Serialize;
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
use BMS::Utilities::Escape;

BEGIN {
    srand( time() ^ ($$ + ($$<<15)) );
}

use vars qw(@ISA);
@ISA   = qw(BMS::Utilities::Escape);

our $jsUndef = 'null';

our $honorOctal = 0;

sub set_literal_token {
    return shift->{LITERAL_TOKEN} = shift;
}

sub set_random_literal_token {
    return shift->{LITERAL_TOKEN} = sprintf("<LITERAL-%f>", rand(10));
}

sub literal_token {
    return shift->{LITERAL_TOKEN};
}

*to_json = \&obj_to_json();
sub obj_to_json {
    my $self = shift;
    my ($obj, $ind, $noNice) = @_;
    my $txt = "";
    if (!defined $obj) {
        $txt = $jsUndef;
    } else {
        my $r = ref($obj);
        if (!$r) {
            # literal
            if ($obj =~ /^0\d+$/ && $honorOctal) {
                # Unquoted leading zeros will cause the number to
                # be interpreted as an octal
                $txt = $obj
            } elsif ($obj =~ /^[\+\-]?(\.?\d+|\d+\.\d+)$/) {
                # Add zero to strip leading zeros
                $txt = $obj + 0;
            } else {
                my $lt = $self->literal_token();
                if ($lt && $obj =~ /^\Q$lt\E(.+)/) {

                    # Request to have the object be literally written
                    # out - no escapes, no quotes. This is needed if
                    # you wish to define a method reference. It
                    # potentially represents a security hole if an
                    # attacker knows the token being used. For that
                    # reason it is suggested that
                    # set_random_literal_token() be used to generate a
                    # new token each time the server runs. The client
                    # should not have knowledge of the token since it
                    # will be stripped from the JSON

                    $txt = $1;

                } else {
                    # Escape and quote the scalar
                    $txt = $self->esc_text( $obj, 'forceQuotes' );
                }
            }
        } elsif ($r eq 'ARRAY') {
            my $joiner = ",";
            my $nxtInd;
            if (defined $ind) {
                $joiner .= "\n" .("  " x $ind);
                $nxtInd  = $ind + 1;
            }
            if ($noNice && $noNice->{basicArray}) {
                # If the array is simple (no nested arrays or hashes)
                # then force simple output
                my $refs = 0; map { $refs++ if (ref($_)) } @{$obj};
                $joiner = ',' unless ($refs);
            }
            $txt = "[ ".join($joiner, map { 
                $self->obj_to_json( $_, $nxtInd, $noNice ) } @{$obj})." ]";
        } elsif ($r eq 'HASH') {
            my $joiner = ",";
            my $nxtInd;
            my @ks     = sort keys %{$obj};
            my $prfx   = "{ ";
            my $sffx   = " }";
            if (defined $ind) {
                my $pad  = "\n". ("  " x $ind);
                $joiner .= $pad;
                $nxtInd  = $ind + 1;
                unless ($#ks < 1) {
                    $prfx .= $pad;
                    # $sffx  = ("  " x $ind) . "}";
                }
            }
            $txt .= $prfx;
            for my $k (0..$#ks) {
                my $key = $ks[$k];
                my $val = $obj->{$key};
                $txt   .= $self->esc_text($key, 'forceQuotes') .": ";
                my @recParam = ($val, $nxtInd, $noNice);
                if ($noNice) {
                    my $d = ++$noNice->{NN_KEY_DEPTH}{$key};
                    if (my $l = $noNice->{$key}) {
                        @recParam = ($val) unless ($d < $l);
                    }
                }
                $txt .= $self->obj_to_json( @recParam );
                $txt .= $joiner unless ($k == $#ks);
                $noNice->{NN_KEY_DEPTH}{$key}-- if ($noNice);
            }
            $txt .= $sffx;
        } else {
            $txt = $jsUndef;
        }
    }
    $txt .= "\n" if (defined $ind && !$ind);
    return $txt;

}
