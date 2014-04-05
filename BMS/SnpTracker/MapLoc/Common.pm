# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::SnpTracker::MapLoc::Common;
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
use BMS::Utilities::Benchmark;

use vars qw(@ISA);
@ISA      = qw(BMS::Utilities 
               BMS::Utilities::Benchmark);

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my ($ml) = @_;
    my $self = {
        MAPLOC => $ml,
    };
    bless ($self, $class);
    return $self;
}

sub dbh     { return shift->{MAPLOC}->{DBH}; }
sub maploc  { return shift->{MAPLOC}; }
sub obj_type { return "Common"; }
sub is_transient { return 0; }
sub text_to_pkey { return shift->maploc->text_to_pkey( @_ ); }
sub cached_text_to_pkey { return shift->maploc->cached_text_to_pkey( @_ ); }

sub read_only {
    my $self = shift;
    if (defined $_[0]) {
        $self->{READ_ONLY} = $_[0] ? 1 : 0;
    }
    return $self->{READ_ONLY};
}

sub is_temp {
    my $self = shift;
    my $nv   = shift;
    if (defined $nv) {
        if ($nv) {
            # Set to zero to prevent this object from writing to database
            $self->{PKEY} = 0;
        } else {
            # delete the PKEY to allow it to be re-calculated
            # Will now allow DB writes
            delete $self->{PKEY};
        }
    }
    return $self->{PKEY} ? 0 : defined $self->{PKEY} ? 1 : 0;
}

# Arguments:
#   Array reference of values unique to the PKEY
#   The name of the table being queried
#   The name of the column holding the PKEY
#   Array reference of column names related to the first argument
sub _generic_get_pkey {
    my $self = shift;
    $self->bench_start();
    my @vals = @{ shift() };
    my $tab  = shift;
    my $ml   = $self->maploc();
    my $rv   = 0;
    my $get  = $ml->{STH}{"GET_PKEY_FOR_$tab"} ||= $ml->dbh->prepare
        ( -name => "Get PKEY for table $tab",
          -sql  => "SELECT $_[0] FROM $tab WHERE ".join
          (" AND ",map { "$_ = ?" } @{$_[1]}));
    # die $get->pretty_print();
    unless ($rv = $get->get_single_value( @vals )) {
        if ($self->{READ_ONLY}) {
            $self->bench_end();
            return 0;
        }
        my $set = $ml->{STH}{"MAKE_PKEY_FOR_$tab"} ||= $ml->dbh->prepare
            ( -name => "Insert PKEY row for $tab",
              -sql  => "INSERT INTO $tab (".join(', ', @{$_[1]}).
              ") VALUES (".join(', ', map { '?' } @{$_[1]}).")",
              -ignore => "duplicate key value" );
        $set->execute(@vals);
        unless ($rv = $get->get_single_value( @vals )) {
            $self->death("Failed to recover PKEY for $tab",
                         map { "$_[1][$_] = $vals[$_]" } (0..$#vals));
        }
    }
    $self->bench_end();
    return $rv;
}

sub _generic_get_id {
    my $self = shift;
    $self->bench_start();
    my ($hkey, $nv) = @_;
    if ($nv) {
        $self->err("Not allowed to set IDs via API!","Request = '$nv'");
    }
    if (!$self->{$hkey."_ID"} && $self->{$hkey}) {
        # The ID is not available but the text is. Just map over to normtxt
        $self->{$hkey."_ID"} = $self->maploc->text_to_pkey( $self->{$hkey} );
    }
    $self->bench_end();
    return $self->{$hkey."_ID"} || 0;
}

sub _generic_getSet_name {
    my $self = shift;
    $self->bench_start();
    my ($hkey, $nv) = @_;
    if ($nv) {
        my $cmp = $self->{$hkey};
        if (!$cmp) {
            if (my $id = $self->_generic_get_id($hkey)) {
                $cmp = $self->maploc->cached_pkey_to_text( $id );
            }
        }
        if ($cmp && $cmp ne $nv) {
            $self->err("Not allowed to reset $hkey",
                       "Request => From '$cmp' to '$nv'");
        } else {
            $self->{$hkey} = $nv;
        }
    }
    if (!$self->{$hkey} && $self->{$hkey."_ID"}) {
        # The name is not available, but the ID is
        $self->{$hkey} = $self->maploc->cached_pkey_to_text
            ( $self->{$hkey."_ID"} )
    }
    $self->bench_end();
    return $self->{$hkey} || "";
}

sub normalized_build {
    my $self = shift;
    my $req  = shift;
    # Ultimately need to put logic here to turn things like '36' into 'NCBI36'
    return $req;
}

our $bestBuildCache = {};
sub best_build {
    my $self = shift;
    my @set = sort map { defined $_ ? $_ : "" } @_;
    my $key = join("\t", @set) || "";
    unless (defined $bestBuildCache->{$key}) {
        my @data;
        foreach my $req (@set) {
            my $bld = $self->normalized_build($req);
            next unless ($bld);
            if ($bld =~ /^[^\d]*(\d+\.\d+)/ || $bld =~ /^[^\d]*(\d+)/) {
                push @data, [$1, $bld];
            } else {
                push @data, [0, $bld];
            }
        }
        my ($best) = sort { $b->[0] <=> $a->[0] || uc($a) cmp uc($b) } @data;
        $bestBuildCache->{$key} = $best->[1] || "";
    }
    return $bestBuildCache->{$key};
}

sub unversion {
    my $self = shift;
    my @rv;
    foreach my $txt (@_) {
        if (defined $txt) {
            $txt =~ s/\.\d+$//;
            push @rv, $txt;
        } else {
            push @rv, "";
        }
    }
    return wantarray ? @rv : $rv[0];
}

# ::Common
sub each_rna_name {
    my $self = shift;
    unless ($self->{RNA_NAMES}) {
        if ($self->{RNA_METHOD_LOOP}++) {
            # Both each_rna() and each_rna_name() are defined generically
            # in ::Common - impossible to resolve without infinite loop
            $self->{RNA_NAMES} = [];
        } else {
            $self->{RNA_NAMES} = [ map { $_->name() } $self->each_rna() ];
        }
    }
    return @{$self->{RNA_NAMES}};
}

# ::Common
sub each_rna_id {
    my $self = shift;
    unless ($self->{RNA_IDS}) {
        $self->{RNA_IDS} = [ map { $_->pkey() } $self->each_rna() ];
    }
    return @{$self->{RNA_IDS}};
}

# ::Common
sub each_unversioned_rna_name {
    my $self = shift;
    $self->{UNV_RNA_NAMES} ||= [ $self->unversion($self->each_rna_name()) ];
    return @{$self->{UNV_RNA_NAMES}};
}

# ::Common
sub each_rna {
    my $self = shift;
    unless ($self->{RNAS}) {
        if ($self->{RNA_METHOD_LOOP}++) {
            # Both each_rna() and each_rna_name() are defined generically
            # in ::Common - impossible to resolve
            $self->{RNAS} = [];
        } else {
            $self->bench_start();
            my $targ = $self->{RNAS} = [];
            my $ml   = $self->maploc();
            foreach my $acc ($self->each_rna_name()) {
                my $vids = $ml->versioned_rna_ids( $acc );
                my ($bestVers) = sort {$b <=> $a } keys %{$vids};
                next unless ($bestVers);
                foreach my $rid (@{$vids->{$bestVers}}) {
                    my $rna = $ml->get_rna_by_id( $rid );
                    push @{$targ}, $rna;
                }
            }
            $self->bench_end();
        }
    }
    return @{$self->{RNAS}};
}

# ::Common
sub each_gene_name {
    my $self = shift;
    unless ($self->{GENE_NAMES}) {
        $self->bench_start();
        my %u;
        my $ml = $self->maploc();
        foreach my $rname ($self->each_rna_name()) {
            map { $u{$_} = 1 } $ml->rna_accession_to_gene_accessions($rname);
        }
        delete $u{""};
        $self->{GENE_NAMES} = [sort keys %u];
        $self->bench_end();
    }
    return @{$self->{GENE_NAMES}};
}

# ::Common
sub each_gene {
    my $self = shift;
    unless ($self->{GENE_OBJS}) {
        my $cache = $_[0] ? $self->maploc()->{GENE_CACHE} ||= {} : {};
        my $targ = $self->{GENE_OBJS} = [ ];
        foreach my $name ($self->each_gene_name()) {
            if (my $gene = $cache->{$name} ||= 
                $self->maploc->get_gene( $name )) {
                push @{$targ}, $gene;
            }
        }
    }
    return @{$self->{GENE_OBJS}};
}
