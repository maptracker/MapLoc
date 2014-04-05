# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::SnpTracker::MapLoc::Accessioned;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

=head1 DESCRIPTION

Inheriting this package allows objects to have "accessions" -
machine-friendly identifiers that can uniquely identify them.

Each accession will also be associated with an "authority" - this is
the source that assigned the accession, for example "NM_001234.1"
would be expected to have an authority along the lines of "RefSeq" or
"Entrez"

=head1 SYNOPSIS

 my $obj; # Recovered by some other means
 my @accStrings = $obj->each_accession();
 my $accDetails = $obj->each_accession();
 my @authStrings = $obj->each_authority();

 my $acc = $obj->best_accession();
 
 $obj->add_accession_to_db("KRZ-00139.31.v1", "Some Authority Name");

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
use vars qw(@ISA);
@ISA      = qw(BMS::SnpTracker::MapLoc::Common);

sub obj_type { return "Accessioned"; }

=head2 

 Title   : each_accession
 Usage   : my @accs = $obj->each_accession()
           my $allData = $obj->each_accession()
 Function: Get a simple list of all accessions assigned to the object,
           or a full list of all accession information
 Returns : In array context, an array of strings.
           In scalar context, a 2D array reference of accessions and authors.

This is the fundamental method for an accessioned object. The first
time it is called, it will query the database with the PKEY for the
object and recover all information, which will be pivoted in several
ways and stored in the object.

In array context, a simple array of strings representing all unique
accessions will be returned. In scalar context, a 2D array reference
containing [ 'accession', 'authority' ] rows will be returned.

=cut

sub each_accession {
    my $self = shift;
    unless ($self->{FULL_ACCS}) {
        $self->bench_start();
        my $ml      = $self->maploc();
        my $bulkSrc = $ml->bulk_cache_source( 'accessions' );
        my $pkey    = $self->pkey();
        my (%byAuth, %byAcc);
        if (my $cache = $bulkSrc->{$pkey}) {
            foreach my $row (@{$cache}) {
                my ($acc, $src) = map { $ml->cached_pkey_to_text( $_ ) } @{$row};
                push @{$byAcc{$acc}}, $src;
                push @{$byAuth{$src}}, $acc;
            }
        } else {
            my $get = $ml->{STH}{GETACCESSION} ||= $ml->dbh->prepare
                ( -name => "Get all accessions for an object",
                  -sql  =>
                  "SELECT nt1.txt, nt2.txt ".
                  "FROM normtxt nt1, normtxt nt2, accession acc ".
                  "WHERE acc.obj_id = ? AND nt1.txt_id = acc.acc_id AND nt2.txt_id = acc.auth_id");
            $get->execute( $pkey );
            my $rows = $self->{FULL_ACCS} = $get->fetchall_arrayref();
            map { push @{$byAcc{$_->[0]}}, $_->[1];
                  push @{$byAuth{$_->[1]}}, $_->[0]; } @{$rows};
        }
        $self->{ACCESSIONS} = [ sort keys %byAcc ];
        $self->{AUTHORITIES} = [ sort keys %byAuth ];
        $self->{AUTH_PLUS_ACC} = \%byAuth;
        $self->{ACC_PLUS_AUTH} = \%byAcc;
        $self->bench_end();
    }
    return wantarray ? @{$self->{ACCESSIONS}} : $self->{FULL_ACCS};
}

=head2 

 Title   : each_authority
 Usage   : my @auths = $obj->each_authority()
           my $allData = $obj->each_authority()
 Function: Get the authorities associated with an object
 Returns : In array context, an array of authority strings.
           In scalar context, a 2D array reference of accessions and authors.

Array context returns a simple non-redundant array of authority
strings. Scalar context returns the same 2D array reference as
each_accession().

=cut

sub each_authority {
    my $self = shift;
    $self->each_accession() unless ($self->{AUTHORITIES});
    return wantarray ? @{$self->{AUTHORITIES}} : $self->{FULL_ACCS};
}

=head2 

 Title   : accessions_with_authorities
 Usage   : my $accHash = $obj->accessions_with_authorities()
 Function: Get accessions and the authorities defining them
 Returns : A hash reference

The method will return a hash reference. The keys will be accessions
(strings) and the values will each be an array reference of
authorities (strings) that have set the corresponding accession.

=cut

sub accessions_with_authorities {
    my $self = shift;
    $self->each_accession() unless ($self->{ACC_PLUS_AUTH});
    return $self->{ACC_PLUS_AUTH};
}

=head2 

 Title   : authors_with_accessions
 Usage   : my $authHash = $self->authors_with_accessions()
 Function: Get authorities and the accessions they have set
 Returns : A hash reference

The returned hash reference will be keyed to authorities (strings) and
each value will be an array reference of accessions (strings) that the
corresponding authority has set.

=cut

sub authors_with_accessions {
    my $self = shift;
    $self->each_accession() unless ($self->{AUTH_PLUS_ACC});
    return $self->{AUTH_PLUS_ACC};
}

=head2 

 Title   : best_accession
 Usage   : my $acc = $obj->best_accession()
 Function: Gets a "good" accession for the object
 Returns : A non-null string

If the objection is a Location, and no accessions are available, then
the full_loc() genomic position will be returned (will look something
like "12.GCRh37:198324")

If the object is a Location with one or more accessions, and at least
one is from dbSNP, then the lowest-numbered dbSNP ID (eg "rs12345")
will be returned.

Otherwise, if at least one accession is known, the lowest
alphabetically-sorted accession is returned.

In all other cases (no accessions, not a Location) the method will
return "MapLoc:PKEY", where PKEY is the internal DB ID (unique within
the database). Note that the value is not portable outside of MapLoc,
nor is it portable to other MapLoc installs.

=cut

sub best_accession {
    my $self = shift;
    $self->bench_start();
    my $type = $self->obj_type();
    my @accs = $self->each_accession();
    my $rv = "";
    if ($type eq 'Location') {
        if ($#accs == -1) {
            # No accessions, use the full location specifier
            $rv = $self->full_loc();
        } else {
            my @sorter;
            foreach my $acc (@accs) {
                my $sc = 0;
                if ($acc =~ /^rs(\d+)$/) {
                    $sc = 200000000 - $1;
                }
                push @sorter, [$acc, $sc];
            }
            my ($best) = sort { $b->[1] <=> $a->[1] } @sorter;
            $rv = $best->[0];
        }
    } elsif ($#accs == -1) {
        $rv = "MapLoc:".$self->pkey();
    } else {
        my ($arbitrary) = sort @accs;
        $rv = $arbitrary;
    }
    $self->bench_end();
    return $rv;
}

=head2 

 Title   : add_accession_to_db
 Usage   : my $err = $obj->add_accession_to_db($acc, $auth)
 Function: Store an accession and authority in the database for the object
 Returns : 0 for success, otherwise error string
 Args    : [0] The accession (string)
           [1] The authority (string)

Used to add an accession to the database.

=cut

sub add_accession_to_db {
    my $self = shift;
    my $ml    = $self->maploc();
    my $accid = $ml->text_to_pkey( shift );
    return "No accession ID recovered" unless ($accid);
    my $lid = $self->pkey();
    return "No PKEY for object" unless ($lid);
    my $authid = $ml->text_to_pkey( shift );
    return "No authority ID recovered" unless ($authid);
    delete $self->{ACCS};
    my $set = $ml->{STH}{SETACCESSION} ||= $ml->dbh->prepare
        ( -name => "Associate accession with an object",
          -sql  =>
          "INSERT INTO accession (obj_id, acc_id, auth_id) VALUES (?,?,?)",
          -ignore => "duplicate key value" );
    $set->execute($lid, $accid, $authid);
    return 0;
}
