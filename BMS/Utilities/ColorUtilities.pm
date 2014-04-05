# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::Utilities::ColorUtilities;
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
use vars qw(@ISA);
@ISA   = qw(BMS::Utilities);

# http://www.w3.org/TR/css3-color/#svg-color
our $hexMap               = {
    ALICEBLUE            => '#F0F8FF',
    ANTIQUEWHITE         => '#FAEBD7',
    AQUA                 => '#00FFFF',
    AQUAMARINE           => '#7FFFD4',
    AZURE                => '#F0FFFF',
    BEIGE                => '#F5F5DC',
    BISQUE               => '#FFE4C4',
    BLACK                => '#000000',
    BLANCHEDALMOND       => '#FFEBCD',
    BLUE                 => '#0000FF',
    BLUEVIOLET           => '#8A2BE2',
    BROWN                => '#A52A2A',
    BURLYWOOD            => '#DEB887',
    CADETBLUE            => '#5F9EA0',
    CHARTREUSE           => '#7FFF00',
    CHOCOLATE            => '#D2691E',
    CORAL                => '#FF7F50',
    CORNFLOWERBLUE       => '#6495ED',
    CORNSILK             => '#FFF8DC',
    CRIMSON              => '#DC143C',
    CYAN                 => '#00FFFF',
    DARKBLUE             => '#00008B',
    DARKCYAN             => '#008B8B',
    DARKGOLDENROD        => '#B8860B',
    DARKGRAY             => '#A9A9A9',
    DARKGREEN            => '#006400',
    DARKGREY             => '#A9A9A9',
    DARKKHAKI            => '#BDB76B',
    DARKMAGENTA          => '#8B008B',
    DARKOLIVEGREEN       => '#556B2F',
    DARKORANGE           => '#FF8C00',
    DARKORCHID           => '#9932CC',
    DARKRED              => '#8B0000',
    DARKSALMON           => '#E9967A',
    DARKSEAGREEN         => '#8FBC8F',
    DARKSLATEBLUE        => '#483D8B',
    DARKSLATEGREY        => '#2F4F4F',
    DARKTURQUOISE        => '#00CED1',
    DARKVIOLET           => '#9400D3',
    DEEPPINK             => '#FF1493',
    DEEPSKYBLUE          => '#00BFFF',
    DIMGRAY              => '#696969',
    DODGERBLUE           => '#1E90FF',
    FIREBRICK            => '#B22222',
    FLORALWHITE          => '#FFFAF0',
    FORESTGREEN          => '#228B22',
    FUCHSIA              => '#FF00FF',
    GAINSBORO            => '#DCDCDC',
    GHOSTWHITE           => '#F8F8FF',
    GOLD                 => '#D4A017',
    GOLDENROD            => '#DAA520',
    GRAY                 => '#808080',
    GREEN                => '#008000',
    GREENYELLOW          => '#ADFF2F',
    GREY                 => '#808080',
    HONEYDEW             => '#F0FFF0',
    HOTPINK              => '#FF69B4',
    INDIANRED            => '#CD5C5C',
    INDIGO               => '#4B0082',
    IVORY                => '#FFFFF0',
    KHAKI                => '#F0E68C',
    LAVENDER             => '#E6E6FA',
    LAVENDERBLUSH        => '#FFF0F5',
    LAWNGREEN            => '#7CFC00',
    LEMONCHIFFON         => '#FFFACD',
    LIGHTBLUE            => '#ADD8E6',
    LIGHTCORAL           => '#F08080',
    LIGHTCYAN            => '#E0FFFF',
    LIGHTGOLDENRODYELLOW => '#FAFAD2',
    LIGHTGRAY            => '#D3D3D3',
    LIGHTGREEN           => '#90EE90',
    LIGHTGREY            => '#D3D3D3',
    LIGHTPINK            => '#FFB6C1',
    LIGHTSALMON          => '#FFA07A',
    LIGHTSEAGREEN        => '#20B2AA',
    LIGHTSKYBLUE         => '#87CEFA',
    LIGHTSLATEGREY       => '#778899',
    LIGHTSTEELBLUE       => '#B0C4DE',
    LIGHTYELLOW          => '#FFFFE0',
    LIME                 => '#00FF00',
    LIMEGREEN            => '#32CD32',
    LINEN                => '#FAF0E6',
    MAGENTA              => '#FF00FF',
    MAROON               => '#800000',
    MEDIUMAQUAMARINE     => '#66CDAA',
    MEDIUMBLUE           => '#0000CD',
    MEDIUMORCHID         => '#BA55D3',
    MEDIUMPURPLE         => '#9370DB',
    MEDIUMSEAGREEN       => '#3CB371',
    MEDIUMSLATEBLUE      => '#7B68EE',
    MEDIUMSPRINGGREEN    => '#00FA9A',
    MEDIUMTURQUOISE      => '#48D1CC',
    MEDIUMVIOLETRED      => '#C71585',
    MIDNIGHTBLUE         => '#191970',
    MINTCREAM            => '#F5FFFA',
    MISTYROSE            => '#FFE4E1',
    MOCCASIN             => '#FFE4B5',
    NAVAJOWHITE          => '#FFDEAD',
    NAVY                 => '#000080',
    OLDLACE              => '#FDF5E6',
    OLIVE                => '#808000',
    OLIVEDRAB            => '#6B8E23',
    ORANGE               => '#FFA500',
    ORANGERED            => '#FF4500',
    ORCHID               => '#DA70D6',
    PALEGOLDENROD        => '#EEE8AA',
    PALEGREEN            => '#98FB98',
    PALETURQUOISE        => '#AFEEEE',
    PALEVIOLETRED        => '#DB7093',
    PAPAYAWHIP           => '#FFEFD5',
    PEACHPUFF            => '#FFDAB9',
    PERU                 => '#CD853F',
    PINK                 => '#FFC0CB',
    PLUM                 => '#DDA0DD',
    POWDERBLUE           => '#B0E0E6',
    PURPLE               => '#800080',
    RED                  => '#FF0000',
    ROSYBROWN            => '#BC8F8F',
    ROYALBLUE            => '#4169E1',
    SADDLEBROWN          => '#8B4513',
    SALMON               => '#FA8072',
    SANDYBROWN           => '#F4A460',
    SEAGREEN             => '#2E8B57',
    SEASHELL             => '#FFF5EE',
    SIENNA               => '#A0522D',
    SILVER               => '#C0C0C0',
    SKYBLUE              => '#87CEEB',
    SLATEBLUE            => '#6A5ACD',
    SLATEGREY            => '#708090',
    SNOW                 => '#FFFAFA',
    SPRINGGREEN          => '#00FF7F',
    STEELBLUE            => '#4682B4',
    TAN                  => '#D2B48C',
    TEAL                 => '#008080',
    THISTLE              => '#D8BFD8',
    TOMATO               => '#FF6347',
    TURQUOISE            => '#40E0D0',
    VIOLET               => '#EE82EE',
    WHEAT                => '#F5DEB3',
    WHITE                => '#FFFFFF',
    WHITESMOKE           => '#F5F5F5',
    YELLOW               => '#FFFF00',
    YELLOWGREEN          => '#9ACD32', 
};


our $hexRE = '[A-F0-9]';
our $digRE = "\\d{1,3}";
our $colorChannels = [ qw(Red Green Blue) ];
sub hex_color {
    my $self = shift;
    my ($col, $alpha) = @_;
    return wantarray ? ($col) : $col unless ($col);
    my $rv;
    if (my $r = ref($col)) {
        if ($r eq 'ARRAY') {
            if ($#{$col} == 1) {
                # Color, Alpha
                ($col, $alpha) = defined $alpha ?
                    ($col->[0], $alpha) : @{$col};
            } else {
                # RGB? And A?
                $alpha = $col->[3] unless (defined $alpha);
                for my $c (0..2) {
                    my $v = $col->[$c];
                    if (!$v) {
                        $col->[$c] = 0;
                    } elsif ($v !~ /^\d+$/) {
                        $self->err("Color specified by array has non-numeric value '$v' for $colorChannels->[$c]");
                        return wantarray ? () : undef;
                    } elsif ($v < 0 || $v > 255) {
                        $self->err("Color specified by array has illegal value '$v' for $colorChannels->[$c]");
                        return wantarray ? () : undef;
                    }
                }
                $rv = sprintf("#%02X%02X%02X", @{$col});
                return wantarray ? ($rv, $alpha) : $rv;
            }
        }
        $self->err("Can not interpret color request '$col'");
        return wantarray ? () : undef;
    }
    $col = uc($col);
    if (my $known = $hexMap->{$col}) {
        # Known color name
        # We will not return in order to allow $hexMap to use other formatting
        $col = $known;
    }
    if ($col =~ /^\#?($hexRE{6})$/) { 
        # Standard 6 digit hex code:
        $rv = "#$1";
    } elsif ($col =~ /^\#?($hexRE)($hexRE)($hexRE)$/) { 
        # 3 digit hex code:
        $rv = "#$1$1$2$2$3$3";
    } elsif ($col =~ /RGBA?\(\s*($digRE)\s*\,\s*($digRE)\s*\,\s*($digRE)\s*(,\s*(1|0?\.\d+)\s*)?\)/) {
        # RGB decimal notation
        $rv = sprintf("#%02X%02X%02X", $1, $2, $3);
        $alpha = $5 unless (defined $alpha);
    }
    if (defined $alpha) {
        if ($alpha eq '') {
            $alpha = undef;
        } elsif ($alpha =~ /^(\d+|\d?\.\d+)$/ && $alpha >= 0 && $alpha <= 1) {
            $alpha += 0;
        } else {
            $self->err("Illegal alpha value '$alpha'");
            $alpha = undef;
        }
    }
    # warn sprintf("hex_color('%s') = [ %s, %s ]\n", map { defined $_ ? $_ : '?' } ($_[0], $rv, $alpha));
    return wantarray ? ($rv, $alpha) : $rv;
}

sub rgb_values {
    my $self = shift;
    my ($hex, $alp) = $self->hex_color(@_);
    my @rv;
    if ($hex && $hex =~ /^\#($hexRE{2})($hexRE{2})($hexRE{2})$/) {
        @rv = map { hex($_) } ($1, $2, $3);
        $rv[3] = $alp if (defined $alp);
    }
    return @rv;
}

*rgb_color = \&rgba_color;
sub rgba_color {
    # RGB plus optional alpha transparency value
    my $self = shift;
    my @rgb  = $self->rgb_values(@_);
    my $rv;
    if ($#rgb == 2) {
        $rv = sprintf("rgb(%d,%d,%d)", @rgb);
    } elsif ($#rgb == 3) {
        $rv = sprintf("rgba(%d,%d,%d,%s)", @rgb);
    }
    # warn join('+', map { defined $_ ? $_ : '?' } @_)."=".(defined $rv ? $rv : '?')."\n";
    return $rv;
}

sub user_text_color {
    my $self = shift;
    my ($text, $col) = @_;
    $text = "" unless (defined $text);
    if (defined $col) {
        my $hex = $self->hex_color( $col );
        if ($hex) {
            $self->{USER_TEXT_COLOR}{$text} = $hex;
        } elsif (defined $hex) {
            # A defined false value means to clear the user setting
            delete $self->{USER_TEXT_COLOR}{$text};
        } else {
            $self->err("Failed to set text color for '$text' to '$col'");
        }
    }
    return $self->{USER_TEXT_COLOR}{$text} || "";
}

sub user_text_color_or_pastel {
    my $self = shift;
    my $txt  = shift;
    return $self->user_text_color( $txt ) ||
        $self->user_text_color( $txt,  $self->pastel_text_color($txt) );
}

sub user_color_hash {
    return shift->{USER_TEXT_COLOR} || {};
}

sub text_to_color {
    my $self = shift;
    my ($text, $sep, $min, $max) = @_;
    if (my $col = $self->user_text_color($text)) {
        return $col;
    }
    $text ||= "";
    if (!$min || $min < 0) {
        $min = 0;
    } elsif ($min > 255) {
        $min = 255;
    }
    if (!$max || $max > 255) {
        $max = 255;
    } elsif ($max < $min) {
        $max = $min;
    }
    my $span = $max - $min;
    if (!$sep || $sep > ($span/3) || $sep < 0) {
        $sep = int($span / 3);
    }
    return "#000000" if ($span < 1);
    # Add up the ordinal values of the characters:
    my $num = 17;
    map { $num += ord($_) } split('', $text);
    # Make it wildly decimal with a weird root:
    $num = ($num ** 0.471) . "";
    # Remove decimal place:
    $num =~ s/\.//g;
    $num ||= "0";
    # Make sure we have at least 3 characters per channel:
    while (length($num) < 9) { $num .= $num; }
    my @vals;
    while ($#vals < 2) {
        my $val = substr($num, -3, 3);
        push @vals, $val % $span;
        #push @vals, $val * $span / 1000;
        $num = substr($num, 0, length($num) - 3);
    }
    my $rotate = $vals[0] % 3;
    if ($sep) {
        my $cen = $vals[0];
        for my $i (1..2) {
            my $diff = $vals[$i] - $cen;
            next if (abs($diff) > $sep);
            if ($diff > 0) {
                $vals[$i] += $sep;
                $vals[$i] -= $span if ($vals[$i] > $span);
            } else {
                $vals[$i] -= $sep;
                $vals[$i] += $span if ($vals[$i] < 0);
            }
        }
    }
    for my $i (1..$rotate) { push @vals, shift @vals; }
    map { $_ += $min } @vals;
    my $hex = sprintf("#%02X%02X%02X", @vals);
    return $hex;
}

sub pastel_text_color {
    my $self = shift;
    my $txt = $_[0] || "";
    return $self->{PASTEL_TEXT_COLOR}{$txt} ||=
        $self->text_to_color( $txt, 80, 30, 250 );
}

sub dark_text_color {
    my $self = shift;
    my $txt = $_[0] || "";
    return $self->{PASTEL_TEXT_COLOR}{$txt} ||=
        $self->text_to_color( $txt, 80, 0, 150 );
}

