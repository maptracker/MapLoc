# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::SnpTracker::MapLoc::Population;
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
use Scalar::Util qw(weaken);
use BMS::SnpTracker::MapLoc::Tagged;
use vars qw(@ISA);
@ISA      = qw(BMS::SnpTracker::MapLoc::Tagged);

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my ($ml, $name, $src, $allowIntegerName) = @_;
    my $self = {
        MAPLOC => $ml,
        PARENT => undef,
    };
    bless ($self, $class);
    $self->bench_start();
    my @errs;
    if (!$name) {

    } elsif (my $r = ref($name)) {
        if ($r =~ /BMS::SnpTracker::MapLoc::Population/) {
            # Already an object
            return $r;
        } else {
            $self->death("No mechanism to make a Population object from '$name'");
        }
    } elsif ($name =~ /^\d+$/ && !$allowIntegerName) {
        # Raw ID provided
        $self->{NAME_ID} = $name;
    } else {
        # Text provided
        $self->{NAME} = $name;
    }
    if (!$src) {
        if (my $pkey = $self->{NAME_ID}) {
            # A single integer will be treated as a primary key
            $self->{PKEY} = $pkey;
            delete $self->{NAME_ID};
        }
    } elsif ($src =~ /^\d+$/) {
        # Raw ID provided
        $self->{SRC_ID} = $src;
    } else {
        # Text provided
        $self->{SRC} = $src;
    }
    my $pkey = $self->pkey();
    unless ($pkey) {
        $self->bench_end();
        return undef;
    }
    my $read = $ml->{STH}{LOAD_POPULATION} ||= $ml->dbh->prepare
        ( -name => "Recover Population info from database",
          -sql  => "SELECT type, chr, par_id, name_id, src_id ".
          "FROM population WHERE pop_id = ?", );
    $read->execute( $pkey );
    my $rows = $read->fetchall_arrayref();
    $self->death("Failed to recover population for $name/$src")
        if ($#{$rows} != 0);
    $self->type( $rows->[0][0] );
    $self->chr( $rows->[0][1] );
    $self->parent_id( $rows->[0][2] );
    $self->{NAME_ID} ||= $rows->[0][3];
    $self->{SRC_ID}  ||= $rows->[0][4];
    $self->bench_end();
    return $self;
}

sub obj_type { return "Population"; }

*pop_id  = \&pkey;
*handle  = \&pkey;
sub pkey {
    my $self = shift;
    unless (defined $self->{PKEY}) {
        $self->bench_start();
        my ($sid, $nid) = ($self->src_id, $self->name_id);
        my $ml = $self->maploc;
        if ($nid && !$sid) {
            # Can we find the population just by name?
            my $read = $ml->{STH}{LOAD_POP_BY_NAME} ||= $ml->dbh->prepare
                ( -name => "Recover Population info from database by name",
                  -sql  => "SELECT type, chr, par_id, pop_id, src_id ".
                  "FROM population WHERE name_id = ?", );
            my $name = $self->name();
            $read->execute( $nid );
            my $rows = $read->fetchall_arrayref();
            if ($#{$rows} == -1) {
                # $self->death("No populations were found with the name '$name'. Please also provide a population source if you wish to create a new population.");
                $self->bench_end();
                return undef;
            } elsif ($#{$rows} != 0) {
                $self->death("A query for populations named '$name' returns more than one result. Please also provide the source to select a distinct population.");                
            }
            $self->type( $rows->[0][0] );
            $self->chr( $rows->[0][1] );
            $self->parent_id( $rows->[0][2] );
            $self->{PKEY}   = $rows->[0][3];
            $self->{SRC_ID} = $rows->[0][4];
            $self->bench_end();
            return $self->{PKEY};
        } elsif (!$sid) {
            $self->death("Can not generate population pkey without src and name");
        }
        my $getSet = $ml->{STH}{"GET_POP_PKEY_FUNC"} ||= $ml->dbh->prepare
            ( -name   => "Get or create population pkey",
              -sql    => "SELECT population_pkey(?,?)", );
        unless ($self->{PKEY} = $getSet->get_single_value( $nid, $sid)) {
            $self->death("Failed to generate or recover pkey for population",
                         "($sid, $nid)");
        }
        #$self->{PKEY} = $self->_generic_get_pkey
        #    ([$nid, $sid], 'population', 'pop_id', ['name_id','src_id'] );
        $self->bench_end();
    }
    return $self->{PKEY};
}

sub type {
    my $self = shift;
    if (my $nv = shift) {
        $self->bench_start();
        if (length($nv) > 20) {
            $self->err("Truncating population type('$nv') to 20 characters");
            $nv = substr($nv, 0, 20);
        }
        $self->{TYPE} = $nv;
        $self->bench_end();
    }
    return $self->{TYPE};
}

*src = \&source;
sub source {
    return shift->_generic_getSet_name( 'SRC', @_ );
}

sub src_id {
    return shift->_generic_get_id( 'SRC', @_ );
}

# ::Population
sub name {
    return shift->_generic_getSet_name( 'NAME', @_ );
}

sub name_id {
    return shift->_generic_get_id( 'NAME', @_ );
}

*ref_pop_id           = \&reference_population;
*normal_pop_id        = \&reference_population;
*normal_population_id = \&reference_population;
sub reference_pop_id {
    my $self = shift;
    unless (defined $self->{REF_POP_ID}) {
        $self->{REF_POP_ID} = $self->maploc->reference_pop_id($self->pkey);
    }
    return $self->{REF_POP_ID};
}

*ref_pop           = \&reference_population;
*normal_pop        = \&reference_population;
*normal_population = \&reference_population;
sub reference_population {
    my $self = shift;
    unless (defined $self->{REF_POP_OBJ}) {
        $self->{REF_POP_OBJ} = 0;
        if (my $rpid = $self->reference_pop_id()) {
            $self->{REF_POP_OBJ} = $self->maploc->get_population($rpid);
        }
    }
    return $self->{REF_POP_OBJ};
}

sub json_data {
    my $self = shift;
    $self->bench_start();
    my $name = $self->name();
    my $par  = $self->parent();
    my $root = $self->root_parent();
    my $rv = {
        name   => $name,
        parent => $par ? $par->name : "",
        root   => $root ? $root->name : "",
        tags   => $self->all_tag_values(),
        pid    => $self->pkey(),
        
    };
    if ($name =~ /^Pop\:(\d+)$/) {
        # dbSNP pop_id
        $rv->{dbSNPid} = $1;
        my $n = $self->simple_value('Handle');
        $n =~ s/_/ /g;
        $rv->{name} = $n;
    }
    $self->bench_end();
    return $rv;
}

sub chr {
    my $self = shift;
    my $nv   = shift;
    if (defined $nv) {
        $self->{CHR} = $nv;
    }
    return $self->{CHR};
}

*par_id = \&parent_id;
sub parent_id {
    my $self = shift;
    if (my $par = $self->{PARENT}) {
        # Check hash to avoid parent_id > parent > ancestors > pkey loop
        return $par->pkey();
    } elsif (my $nv = shift) {
        if ($nv !~ /^\d+$/) {
            $self->death("Attempt to set parent_id( $nv )");
        } elsif (my $par = $self->parent($nv)) {
            return $par->pkey();
        }
    }
    return undef;
}

sub parent {
    my $self = shift;
    if (defined $_[0]) {
        $self->bench_start();
        if ($_[0]) {
            # Setting a value
            my $ml = $self->maploc();
            my $par = ref($_[0]) ? $_[0] : $ml->get_population_fast( @_ );
            unless ($par) {
                $self->err
                    ("Could not understand '".join(',', @_).
                     "' as a parent Population");
                $self->bench_end();
                return undef;
            }
            $self->{PARENT} = $par;
            # Need to make sure it does not cause a cycle
            $self->ancestors();
        } else {
            # Clearing the parent
            $self->{PARENT} = undef;
        }
        $self->bench_end();
    }
    return $self->{PARENT};
}


sub root_parent {
    my $self = shift;
    my $rv   = $self;
    while (my $par = $rv->parent()) {
        $rv = $par;
    }
    return $rv->pkey == $self->pkey ? undef : $rv;
}

sub ancestors {
    my $self = shift;
    $self->bench_start();
    my %seen = ( $self->pkey => 1);
    my @rv;
    # Include a true argument to indicate the query should also be returned
    push @rv, $self if (shift);
    my $node = $self;
    while (my $par = $node->parent()) {
        push @rv, $par;
        if ($seen{$par->pkey}++) {
            $self->death("Loop structure detected in parents",
                         map { $_->to_one_line() } ($self, @rv));
        }
        $node = $par;
    }
    $self->bench_end();
    return @rv;
}

sub colored_tag {
    my $self = shift;
    $self->read();
    my @ctags = $self->name();
    push @ctags, $self->tag('Tally');
    push @ctags, $self->tag('Category');
    push @ctags, map { $_->name } $self->ancestors();
    my $ml = $self->maploc();
    my @rv;
    foreach my $tag (@ctags) {
        if (my $col = $ml->user_text_color($tag)) {
            @rv = ($tag, $col);
            last;
        }
    }
    $rv[0] ||= $ctags[-1];
    return wantarray ? @rv : $rv[0];
}

sub child_ids {
    my $self = shift;
    my $rv   = $self->{KID_IDS};
    # Pass a true value to force an update:
    if (!$rv || shift) {
        $self->bench_start();
        my $get = $self->dbh->prepare
            ( -name => "Recover child populations",
              -sql  => "SELECT pop_id FROM population WHERE par_id = ?", );
        my @ids = $get->get_array_for_field( $self->pkey );
        $rv = $self->{KID_IDS} ||= {};
        map { $rv->{$_} ||= undef } @ids;
        $self->bench_end();
    }
    return wantarray ? keys %{$rv} : $rv;
}

sub children {
    my $self = shift;
    $self->bench_start();
    my $hash = $self->child_ids( @_ );
    my @rv;
    foreach my $id (keys %{$hash}) {
        my $ml = $self->maploc();
        unless ($hash->{$id}) {
            weaken( $hash->{$id} = $ml->get_population_fast( $id ) );
        }
        push @rv, $hash->{$id};
    }
    $self->bench_end();
    return @rv;
}

sub deep_child_ids {
    my $self = shift;
    my $rv = $self->{ALL_KID_IDS};
    # Pass a true value to force an update:
    if (!$rv || shift) {
        $self->bench_start();
        my $get = $self->dbh->prepare
            ( -name => "Recover child populations",
              -sql  => "SELECT pop_id FROM population WHERE par_id = ?", );
        my @stack = ($self->pkey);
        $rv = $self->{ALL_KID_IDS} ||= {};
        while (my $parid = shift @stack) {
            next if (exists $rv->{$parid});
            $rv->{$parid} = undef;
            my @cc = $get->get_array_for_field( $parid );
            push @stack, @cc;
            # $self->msg("$parid = $#cc") if ($#cc > -1);
        }
        $self->bench_end();
    }
    return wantarray ? keys %{$rv} : $rv;
}

sub deep_children {
    my $self = shift;
    $self->bench_start();
    my $hash = $self->deep_child_ids( @_ );
    my @rv;
    foreach my $id (keys %{$hash}) {
        my $ml = $self->maploc();
        unless ($hash->{$id}) {
            weaken( $hash->{$id} = $ml->get_population_fast( $id ) );
        }
        push @rv, $hash->{$id};
    }
    $self->bench_end();
    return @rv;
}

sub update {
    my $self = shift;
    return if ($self->{READ_ONLY});
    $self->bench_start();
    my $ml   = $self->maploc();
    my $set  = $ml->{STH}{UPDATE_POPULATION} ||= $ml->dbh->prepare
        ( -name => "Update population data",
          -sql  => "UPDATE population SET type = ?, chr = ?, par_id = ? ".
          "WHERE pop_id = ?");
    my $chr = $self->chr();
    $chr = undef if (!$chr || $chr eq '0');
    $set->execute( $self->type(), $chr, $self->parent_id, $self->pkey() );
    $self->write_tags();
    $self->bench_end();
    return $self;
}

sub read {
    my $self = shift;
    $self->bench_start();
    $self->read_tags();
    $self->bench_end();
    return $self;
}

# ::Population
sub to_text {
    my $self = shift;
    $self->bench_start();
    my $pre  = shift || "";
    my $txt  = $pre . $self->to_one_line(). "\n";
    $txt .= $self->tag_text($pre."  ");
    $self->bench_end();
    return $txt;
}

sub to_one_line {
    my $self = shift;
    $self->bench_start();
    my $txt  = sprintf("%s [%d] %s", 
                       $self->name(), $self->handle, $self->src());
    if (my $chr = $self->chr) {
        $txt .= sprintf(" %d chr", $chr);
    }
    if (my $par = $self->parent()) {
        $txt .= sprintf(" (child of %s/%s)", $par->name(), $par->src());
    }
    $self->bench_end();
    return $txt;
}

# ::Population
*to_full_text = \&full_text;
sub full_text {
    my $self = shift;
    $self->bench_start();
    my $pre  = shift || "";
    my $txt  = $self->to_text( $pre );
    my @kids = $self->children();
    my $kn   = $#kids + 1;
    if ($kn) {
        $txt .= sprintf
            ("%s  %d child population%s:\n", $pre, $kn, $kn == 1 ? '' : 's');
        foreach my $kid (sort { $a->name cmp $b->name } @kids) {
            $txt .= "$pre    ". $kid->to_one_line()."\n";
        }
    }
    $self->bench_end();
    return $txt;
}

sub all_locations {
    my $self = shift;
    my $args = $self->parseparams( @_ );
    my @ids  = ($self->pkey);
    push @ids, $self->deep_child_ids( $args->{UPDATE} )
        unless ($args->{NOKIDS} || $args->{NOCHILDREN} );
    my $sql = <<SQL;
SELECT l.loc_id
  FROM location l, allele a
 WHERE 
SQL

    $self->death("NEED TO IMPLEMENT");
    
}
