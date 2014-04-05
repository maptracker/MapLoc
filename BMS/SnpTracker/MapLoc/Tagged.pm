# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::SnpTracker::MapLoc::Tagged;
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
use BMS::SnpTracker::MapLoc::Common;

use vars qw(@ISA);
@ISA      = qw(BMS::SnpTracker::MapLoc::Common);


sub obj_type { return "Tagged"; }

sub all_tag_values {
    my $self = shift;
    $self->bench_start();
    my %rv;
    my $ml = $self->maploc();
    while (my ($tid, $tagDat) = each %{$self->{TAGIDS} || {}}) {
        my $tag   = $tagDat->[1] ||= $ml->cached_pkey_to_text($tid);
        while (my ($vid, $val) = each %{$tagDat->[0]}) {
            # I was using un-cached pkey_to_text() here under the
            # presumption that values would be much more heterogeneous than
            # tags. However, for my data at least, there were a significant
            # number of values reused as well, so am using the cached method:
            $val ||= $tagDat->[0]{$vid} ||= $ml->cached_pkey_to_text($vid);
            push @{$rv{$tag}}, $val;
        }
    }
    if (my $tt = $self->{TRANSIENT_TAGS}) {
        while (my ($tag, $valH) = each %{$tt}) {
            push @{$rv{$tag}}, keys %{$valH};
        }
    }
    $self->bench_end();
    return \%rv;
}

sub tag_text {
    my $self = shift;
    $self->bench_start();
    my $pre  = shift || "";
    my $txt  = "";
    my $atvs = $self->all_tag_values();
    foreach my $tag (sort keys %{$atvs}) {
        $txt .= sprintf("%s%s\n", $pre, $tag);
        foreach my $val (@{$atvs->{$tag}}) {
            $txt .= sprintf("  %s%s\n", $pre, $val);
        }
    }
    $self->bench_end();
    return $txt;
}

sub tag_html {
    my $self = shift;
    $self->bench_start();
    my $html = "";
    my $atvs = $self->all_tag_values();
    my @tags = sort { $#{$b->[1]} <=> $#{$a->[1]} || 
                          uc($a->[0]) cmp uc($b->[0]) } map { [$_, $atvs->{$_} ] } keys %{$atvs};
    if ($#tags == -1) {
        $html .= "<span class='note'>No tags</span><br />\n";
    } else {
        $html .= "<div class='tagcloud'>\n";
        foreach my $tv (@tags) {
            my ($t, $vs) = @{$tv};
            $html .= " <div class='tagbit'><div class='tag'>$t</div> ";
            $html .= join(" ", map { "  <div class='val'>$_</div>" } @{$vs});
            $html .= " </div>\n";
        }
        $html .= "</div>\n";
    }
    $self->bench_end();
    return $html;

}

sub set_transient_tag {
    my $self = shift;
    if (my $tag = shift) {
        my $val = shift;
        if (defined $val) {
            if ($val ne "") {
                $self->{TRANSIENT_TAGS}{$tag}{$val} = 1;
            } else {
                delete $self->{TRANSIENT_TAGS}{$tag};
            }
        }
    }
}

*tag        = \&set_tag_val;
*set_tag    = \&set_tag_val;
*tag_values = \&set_tag_val;
sub set_tag_val {
    my $self = shift;
    if (my $tag = shift) {
        my $tid  = $self->maploc->cached_text_to_pkey($tag);
        my $targ = $self->{TAGIDS}{$tid} ||= [ {}, $tag ];
        my $vH   = $targ->[0];
        my $ml   = $self->maploc();
        foreach my $val (@_) {
            if (defined $val && $val ne "") {
                my $vid = $ml->text_to_pkey($val);
                $vH->{$vid} ||= $val;
            }
        }
        my @rv = sort map { $vH->{$_} ||= 
                                $ml->pkey_to_text( $_ ) } keys %{$vH};
        if (my $tt = $self->{TRANSIENT_TAGS}) {
            if (exists $tt->{$tag}) {
                push @rv, keys %{$tt->{$tag}};
            }
        }

        return @rv;
    }
    return ();
}

sub set_tag_val_ids {
    my $self = shift;
    if (my $tid = shift) {
        my $targ = $self->{TAGIDS}{$tid} ||= [ {} ];
        foreach my $vid (@_) {
            $targ->[0]{$vid} ||= "";
        }
        return keys %{$targ->[0]};
    }
    return ();
}

sub all_tag_id_pairs {
    my $self = shift;
    $self->bench_start();
    my @rv;
    while (my ($tid, $tarr) = each %{$self->{TAGIDS} || {}}) {
        push @rv, map { [$tid, $_ ] } keys %{$tarr->[0] || {}};
    }
    $self->bench_end();
    return wantarray ? @rv : \@rv;
}

*single_value = \&unique_value;
sub unique_value {
    my $self = shift;
    my @vals = $self->tag_values( $_[0] );
    return $#vals == 0 ? $vals[0] : "";
}

sub comma_tag_value {
    my $self = shift;
    my @vals = $self->tag_values( $_[0] );
    return join(',', sort @vals) || "";
}

sub simple_value {
    my $self = shift;
    my ($rv) = sort { length($a) <=> length($b) } $self->tag_values( $_[0] );
    return $rv || "";
}

*description = \&desc;
sub desc        { return shift->simple_value( 'Description' ); }
sub taxa        { return shift->simple_value( 'Taxa' ); }
sub symbol      { return shift->simple_value( 'Symbol' ); }
sub all_symbols { return shift->tag_values('Symbol'); }

sub list_of_values {
    my $self = shift;
    my @vals = $self->tag_values( $_[0] );
    return join(',', @vals) || "";
}

*remove_tag = \&clear_tag;
*delete_tag = \&clear_tag;
sub clear_tag {
    my $self = shift;
    if (my $tag = shift) {
        return if ($self->{READ_ONLY});
        $self->bench_start();
        my $ml  = $self->maploc();
        my $sth = $ml->{STH}{CLEAR_TAG} ||= $ml->dbh->prepare
            ( -name => "Remove all tagvals for an object + tag",
              -sql  => "DELETE FROM tagval WHERE obj_id = ? AND tag_id = ?");
        my $tid   = $ml->cached_text_to_pkey($tag);
        my @binds = ($tid);
        

        if (my $val = shift) {
            $sth = $ml->{STH}{CLEAR_TAG_FOR_VAL} ||= $ml->dbh->prepare
                ( -name => "Remove all tagvals for an object + tag + val",
                  -sql  => "DELETE FROM tagval WHERE obj_id = ? AND tag_id = ? AND val_id = ?");
            my $vid = $ml->cached_text_to_pkey($val);
            push @binds, $vid;
            if (my $targ = $self->{TAGIDS}{$tid}) {
                delete $targ->[0]{$vid};
            }
        } else {
            delete $self->{TAGIDS}{$tid};
        }
        my $pkey = $self->pkey();
        $sth->execute( $pkey, @binds);
        $self->bench_end();
    }
}

sub read_tags {
    my $self = shift;
    # Pass a true value to force a read:
    my $force = shift;
    my $rows;
    if (!$self->{TAGS_READ} || $force) {
        if (my $pkey = $self->pkey()) {
            $self->bench_start();
            my $ml   = $self->maploc();
            my $bcs  = $ml->bulk_cache_source( 'tagval' );
            my $rows;
            if (!$force && exists $bcs->{$pkey} && $bcs->{$pkey}) {
                $rows = $bcs->{$pkey};
            } else {
                # warn "NEW TAGS FOR ".$self->to_one_line();
                my $get  = $ml->{STH}{READ_TAGS} ||= $ml->dbh->prepare
                    ( -name => "Add a tag-value pair for an object",
                      -sql  => "SELECT tag_id, val_id FROM tagval WHERE obj_id = ?" );
                $get->execute($pkey);
                $rows = $get->fetchall_arrayref();
            }
            foreach my $row (@{$rows}) {
                my ($tid, $vid) = @{$row};
                my $targ = $self->{TAGIDS}{$tid} ||= [ {} ];
                $targ->[0]{$vid} ||= "";
            }
            $self->{TAGS_READ} = 1;
            $self->bench_end();
        }
    }
    return $rows;
}

sub write_tags {
    my $self = shift;
    my $tvs  = $self->all_tag_id_pairs();
    return 0 if ($#{$tvs} == -1);
    return 0 if ($self->{READ_ONLY});
    $self->bench_start();
    my $pkey = $self->pkey();
    $self->death("No PKEY available for ".$self->obj_type()." object!", $self)
        unless ($pkey);
    my $ml   = $self->maploc();
    #my $set  = $ml->{STH}{ADD_TAG_VAL} ||= $ml->dbh->prepare
    #    ( -name => "Add a tag-value pair for an object",
    #      -sql  => "INSERT INTO tagval (obj_id, tag_id, val_id) VALUES (?, ?, ?)",
    #      -ignore => "duplicate key value" );
    my $set  = $ml->{STH}{ADD_TAG_VALFUNC} ||= $ml->dbh->prepare
        ( -name => "Add a tag-value pair for an object by function",
          -sql  => "SELECT quiet_tagval_write(?,?,?)",
          -ignore => "duplicate key value" );
    foreach my $tv (@{$tvs}) {
        my ($t, $v) = @{$tv};
        if ($v) {
            $set->get_single_value( $pkey, $t, $v );
        } else {
            # Do we want to delete null tags here?
        }
    }
    $self->bench_end();
    return 1;
}
