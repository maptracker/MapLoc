# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::SnpTracker::MapLoc::Text;
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
use BMS::SnpTracker::MapLoc::Tagged;

use vars qw(@ISA);
@ISA      = qw(BMS::SnpTracker::MapLoc::Tagged);

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my ($ml, $txt) = @_;
    return undef unless ($txt);
    my $self = {
        TEXT   => $txt,
        MAPLOC => $ml,
    };
    bless ($self, $class);
    return $self;
}

sub obj_type { return "Text"; }

# ::Text
sub to_one_line {
    my $self = shift;
    return sprintf("[%d] %s", $self->handle(), $self->text());
}

# ::Text
sub to_text {
    my $self = shift;
    $self->bench_start();
    my $pre  = shift || "";
    my $txt  = sprintf("%s%s [%s]\n", $pre, $self->text(), $self->handle());
    $txt    .= $self->tag_text($pre."  ");
    $self->bench_end();
    return $txt;
}

sub to_html {
    my $self = shift;
    $self->bench_start();
    my $html = sprintf("<div class='%s'><span class='medlabel'>%s</span> <span class='pkey'>%s</span><br />\n", $self->obj_type(), $self->text(), $self->handle());
    $html .= $self->tag_html();
    $html .= "</div>\n";
    $self->bench_end();
    return $html;
}

*name = \&text;
*txt = \&text;
sub text {
    return shift->{TEXT};
}

*text_id = \&pkey;
*txt_id  = \&pkey;
*id      = \&pkey;
*handle  = \&pkey;
sub pkey {
    my $self = shift;
    return $self->{TXT_ID} ||= $self->maploc->text_to_pkey( $self->{TEXT} );
}

*is_in_db = \&is_known;
sub is_known {
    my $self = shift;
    return ($self->{TXT_ID}) if ($self->{TXT_ID});
    $self->bench_start();
    my $ml   = $self->maploc();
    my $txt  = $self->text();
    my $sth  = $ml->{STH}{CHECK_FOR_KNOWNTEXT} ||= $ml->dbh->prepare
        ( -name => "Check for existing text in DB",
          -sql  => "SELECT txt_id FROM normtxt ".
          "WHERE md5sum = md5(?) AND txt = ?");
    my @ids = $sth->get_array_for_field( $txt, $txt);
    $self->bench_end();
    return 0 if ($#ids == -1);
    return $self->{TXT_ID} = $ids[0];
}

sub read {
    my $self = shift;
    $self->read_tags();
}

# ::Text
*update = \&write;
sub write {
    my $self = shift;
    $self->write_tags();
}

# ::Text
sub json_data {
    my $self = shift;
    $self->bench_start();
    $self->read();
    my $rv = {
        text    => $self->text(),
        tags    => $self->all_tag_values(),
        pid     => $self->pkey(),
    };
    $self->bench_end();
    return $rv;
}
