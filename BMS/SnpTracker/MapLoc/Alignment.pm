# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::SnpTracker::MapLoc::Alignment;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

=head1 DESCRIPTION

Object representing an alignment between two sequences.

=head1 SYNOPSIS

 my $newAln = BMS::SnpTracker::MapLoc::Alignment->new( $ml );
 my $AlnInDB = BMS::SnpTracker::MapLoc::Alignment->new( $ml, $dbid );

 
 
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
@ISA      = qw(BMS::SnpTracker::MapLoc::TaggedHSP);

=head1 Basic Object Methods

=head2 

 Title   : new
 Usage   : my $aln = BMS::SnpTracker::MapLoc::Alignment->new( $ml );
           my $aln = BMS::SnpTracker::MapLoc::Alignment->new( $ml, $dbid );
 Function: Create a new Alignment object
 Returns : A blessed Perl object
 Args    : [0] A BMS::SnpTracker::MapLoc object
           [1] Optional DB PKEY, if recovering an existing alignment

A MapLoc object must be provided as the first argument (provides DB
connectivity as well as other functions).

In the absence of other arguments, a new, unstructured alignment will
be created.

If a database PKEY (integer) is provided as the second argument, it
will be used to recover information stored in the database under that
ID.

=cut

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my ($ml, $id) = @_;
    my ($s1, $s2, $src, $str, $sc, $hb);
    if ($id) {
        my $read = $ml->{STH}{READ_ALIGNMENT_BY_PKEY} ||= $ml->dbh->prepare
            ( -name => "Read basic alignment information",
              -sql  => "SELECT seq1, seq2, src_id, strand, score, howbad ".
              "FROM alignment WHERE aln_id = ?");
        $read->execute($id);
        my $rows = $read->fetchall_arrayref();
        return undef if ($#{$rows} == -1);
        ($s1, $s2, $src, $str, $sc, $hb) = @{$rows->[0]};
    }
    my $self = {
        PKEY   => $id,
        MAPLOC => $ml,
        SRC_ID => $src,
        SEQIDS => ($s1 ? [$s1,$s2] : undef),
        STRAND => $str ? [1, $str] : undef,
        SCORE  => $sc,
        HOWBAD => $hb,
        SEQNUM => 2,
    };
    bless ($self, $class);
    return $self;
}

sub obj_type { return "Alignment"; }

=head2 

 Title   : pkey
 Usage   : my $dbid = $aln->pkey( $newVal )
 Function: Gets / sets the database ID for the alignment
 Returns : The database primary key (integer), or undef if not set
 Args    : [0] Optional new value

=cut

*handle  = \&pkey;
sub pkey {
    my $self = shift;
    if (my $nv = shift) {
        $self->{PKEY} = $nv;
    }
    return $self->{PKEY};
}

=head2 

 Title   : source
 Usage   : my $src = $aln->source( $newVal )
 Function: Gets/Sets the 'source' of the alignment
 Returns : The source (string)
 Args    : [0] Optional new value

Returns the source as a string.

If an argument is passed, will set the source. Will throw an error if
the source is set and the new value is different.

=cut

*src = \&source;
sub source {
    return shift->_generic_getSet_name( 'SRC', @_ );
}

=head2 

 Title   : src_id
 Usage   : my $txtid = $aln->src_id()
 Function: Gets the txt_id associated with the source of the object
 Returns : A normtxt txt_id (integer)

Gets the database txt_id associated with the string recovered in
source(). Will throw an error if an attempt is made to set the value.

=cut

sub src_id {
    return shift->_generic_get_id( 'SRC', @_ );
}

=head2 

 Title   : score
 Usage   : my $sc = $aln->score( $newVal )
 Function: Gets/Sets the score for the alignment
 Returns : The score
 Args    : [0] Optional new value

The score should be a floating point value. In general, this will be a
percentage between 0 and 100.

=cut

sub score {
    my $self = shift;
    if (defined $_[0]) {
        $self->{SCORE} = $_[0];
    }
    return $self->{SCORE};
}

=head2 

 Title   : howbad
 Usage   : my $hb = $aln->howbad()
 Function: Gets a measure of "how bad" the alignment is
 Returns : undef, or a floating point number
 Args    : [0] Optional new value

This concept can be a bit confusing. In some cases, a "thing" can have
more than one reasonable alignments. As a particular example, an RNA
might have two or more genomic locations with nearly identical (or
even identical) alignment metrics (eg % matched). It can be
advantageous to record multiple alignment locations for the RNA:

 * You are not sure which of the locations is the "real" one
 * You believe each alignment is a distinct locus in a repeat family
 * Each alignment is highlighting a potentially interesting region

So in some cases you may want to roll up your sleaves and dig into
each location in turn, savoring the rich and subtle complexity of the
genome.

In other cases (ie 99% of your users) such multiplicity will be cause
for confusion, distress, or even loud words exchanged in the hallway
with vigorous gesticulation. In these cases, you want to report only
the "best" alignment. However, you can not rely on picking the one
with 100% match to the genome, because frequently even the best
alignment has some mismatches. This is where a "howbad" value comes in
handy.

HowBad is calculated as:

  Best Score - This Score

So the best alignment(s) will have a howbad of zero, regardless of
their actual score. If the best score was 97.3, and an alignment has a
score of 95.1, then it will have howbad=2.2.

These values can not be computed automatically - it rests on alignment
loading algorithms to compute them.

=cut

*how_bad = \&howbad;
sub howbad {
    my $self = shift;
    if (defined $_[0]) {
        $self->{HOWBAD} = $_[0];
    }
    return $self->{HOWBAD};
}

=head1 Database Methods

Methods that read or write the Postgres database

=head2 

 Title   : forced_pkey
 Usage   : my $dbid = $aln->forced_pkey()
 Function: Will return the database PKEY, or create one if it does not exist
 Returns : The database primary key (integer)

=cut

sub forced_pkey {
    my $self = shift;
    return $self->{PKEY} ||= $self->maploc->dbh->nextval('global_seq');
}

=head2 

 Title   : read_or_force_pkey
 Usage   : my $dbid = $aln->read_or_force_pkey( @args )
 Function: Tries to match the alignment to on in DB, or makes a new one
 Returns : The database primary key (integer)
 Args    : See below

This method should be used when the alignment boundaries have been
defined. The method will try to match the coordinates of this
alignment to one in the database. If a match is found, it will use
that database entry - any changes made to the object will then
overwrite the DB entry if a write() call is made.

If a suitable match is not found, then forced_pkey() is used to
generate a new alignment entry in the database.

See read_id_by_seq_and_bounds() for details on the recognized
parameters.

=cut

sub read_or_force_pkey {
    my $self = shift;
    $self->read_id_by_seq_and_bounds(@_);
    return $self->forced_pkey();
}

=head2 

 Title   : read_id_by_seq_and_bounds
 Usage   : my $dbid = $aln->read_id_by_seq_and_bounds( @args )
 Function: Find the aln_id of entries in the database that match the
           properties of the Alignment object
 Returns : In array context, an array of all matching aln_id values
           In scalar context a single aln_id is returned if there are no
           other IDs with the same or better 'slop'
 Args    : [0] A sequence identifier, generally a chromosome or contig
           [1] A slop value
           [2] Optional second sequence identifier
           [3] Optional authority 

For this function to work, the Alignment object must have the
sequences and HSP information defined. These values will be searched
against the database to see if there are already stored alignments
that match the one in the Perl object.

The sequence identifier will be parsed by
BMS::SnpTracker::MapLoc::HSPed::seqindex() - it will accept names,
seq_ids, or index values relative to the used object.

Coordinate matching will be done only on the extremes of the HSPs -
this method does not scrutinize individual HSPs, just the minimum and
maximum coordinates.

Slop is an optional value that allows imperfect boundary matching for
the HSPs. The provided value will be used to expand the coordinate
space considered a match. So a value of 14,302 will "match" 14,321 if
the slop is 19 or greater.

The second sequence identifer can be left out (in which case only the
'main' sequence will be matched), or provided as either an explicit
value, or as a wildcard. Wildcards are intended to match versioned
IDs, such as "NM_001234.%"

An optional authority can be used, such as "Sim4".

=cut

sub read_id_by_seq_and_bounds {
    # Use one of the sequence IDs and its boundaries to recover existing
    # alignments in the database. This will generally be useful when
    # the query is genomic DNA.
    my $self = shift;
    $self->bench_start();
    my ($req, $slop, $useBoth, $auth) = @_;
    my $ind  = $self->seqindex($req);
    unless ($ind) {
        $self->death("Can not set read_id with unrecognized sequence '$req'");
    }
    my @bnds = $self->bounds();
    if ($#bnds == -1) {
        $self->death("Can not set read_id by sequence ($req) until bounds() have been set");
    }
    my $ml   = $self->maploc();
    my @sids = $self->seq_ids();
    my $seq1 = $sids[$ind];
    my ($l1, $r1) = ($bnds[2 * $ind], $bnds[1 + 2 * $ind]);
    my $sname = "FIND_ALIGN_";
    $sname   .= "_SLOPPY" if ($slop);
    $sname   .= !$useBoth ? "BY SEQ" : $useBoth =~ /\%/ ?
        "BY_TWO_WC_SEQS" : "BY_TWO_SEQS";
    my $sth   = $ml->{STH}{$sname};
    unless ($sth) {
        # Build the SQL
        my $desc  = "Find alignment IDs for";
        my $repP  = "REPLACE_PRIMARY";
        my $repO  = "REPLACE_OTHER";
        my $base  = "SELECT a.aln_id, ";
        my $ptTst = "= ?";
        if ($slop) {
            # Need to capture the distance from the current boundaries
            $base .= "abs(a.ptl$repP - ?) + abs(a.ptr$repP - ?)";
            $ptTst = "BETWEEN ? AND ?";
            $desc .= " sloppy";
        } else {
            $base .= "0";
            $desc .= " exact";
        }
        $desc .= " boundary matches against";
        $base .= " FROM alignment a";
        my $ubSql = "";
        if ($useBoth) {
            # We also want to match the second seq
            $desc .= " a primary sequence and a secondary";
            if ($useBoth =~ /\%/) {
                # with an explicit wildcard string
                $base .= ", normtxt nt";
                $ubSql = " AND (upper(substr(nt.txt, 1, 20))) LIKE upper(?) ".
                    "AND a.seq$repO = nt.txt_id";
                $desc .= " wild-carded";
            } else {
                $ubSql = " AND a.seq$repO = ?";
            }
            $desc .= " sequence";
        } else {
            $desc .= " a single sequence";
        }
        if ($auth) {
            my $aid  = $ml->text_to_pkey($auth);
            $ubSql .= " AND a.src_id = $aid";
        }

        $base .= " WHERE ".
            "a.seq$repP = ? AND a.ptl$repP $ptTst AND a.ptr$repP $ptTst".
            $ubSql;
        my $sql = $base;
        # Add the parts where primary is seq1 and secondary is seq2
        $sql    =~ s/$repP/1/g;
        $sql    =~ s/$repO/2/g;
        # Then union with 1/2 switched:
        $sql   .= " UNION ";
        $base   =~ s/$repP/2/g;
        $base   =~ s/$repO/1/g;
        $sql   .= $base;

        $sth    = $ml->{STH}{$sname} = $ml->dbh->prepare
                ( -name => $desc,
                  -sql  => $sql);
    }

    my @binds;
    my @locBind = ($l1, $r1);
    if ($slop) {
        # We need two bind values in the SELECT clause, which are used
        # to see how far off from 'ideal' we are:
        push @binds, @locBind;
        # The = clause is now between:
        @locBind = map { $_ - $slop, $_ + $slop } @locBind;
    }
    push @binds, $seq1, @locBind;
    if ($useBoth) {
        my $seq2 = $sids[$ind ? 0 : 1];
        push @binds, $seq2;
    }
    $sth->execute( @binds, @binds);
    my $found = $sth->fetchall_arrayref();
    my @rv = sort { $a->[1] <=> $b->[1] } @{$found};
    my $sv;
    if ($#rv == -1){
        # Nothing found
    } elsif ($#rv == 0 || $rv[0][1] < $rv[1][1]) {
        # One hit found, or one hit has better slop than the others
        $sv = $self->pkey( $rv[0][0] );
    }
    $self->bench_end();
    return wantarray ? @rv : $sv;
}

=head2 

 Title   : write
 Usage   : $aln->write( $noOverwrite )
 Function: Writes the alignment to the database
 Args    : [0] Flag, if true then will not overwrite existing data

Will commite the object to the database. The flag will prevent changes
to table alignment if an entry matching the PKEY already exists.

The method will also call BMS::SnpTracker::MapLoc::HSPed::write_hsps().

=cut

sub write {
    my $self = shift;
    return if ($self->{READ_ONLY});
    $self->bench_start();
    my ($carefully) = @_;
    my $ml    = $self->maploc();
    my $dbh   = $ml->dbh;
    my $pkey  = $self->forced_pkey();
    my $srcid = $self->src_id();
    my @sids  = $self->seq_ids();
    my @bnds  = $self->bounds();
    my $sc    = $self->score();
    my $hb    = $self->howbad();
    my $str   = $self->strand();
    my @errs;
    push @errs, "Must provide source() or src_id()" unless ($srcid);
    push @errs, "Must provide at least one HSP via hsps()" if ($#bnds == -1);
    push @errs, "Must provide two seqs()" if ($#sids != 1);
    $self->death("Failed to write map", @errs) unless ($#errs == -1);

    my $setUp = $ml->{STH}{INIT_ALIGNMENT} ||= $dbh->prepare
        ( -name => "Make sure alignment row is present",
          -sql  => "INSERT INTO alignment (aln_id) VALUES (?)",
          -ignore => "duplicate key value" );
    my $check = $ml->{STH}{CHECK_ALIGNMENT_EXISTANCE} ||= $dbh->prepare
        ( -name => "See if an alignment is present in the DB",
          -sql  => "SELECT aln_id FROM alignment WHERE aln_id = ?");
    my $update = $ml->{STH}{UPDATE_ALIGNMENT} ||= $dbh->prepare
        ( -name => "Update alignment values",
          -sql  => "UPDATE alignment SET seq1 = ?, seq2 = ?, src_id = ?, ".
          "strand = ?, score = ?, howbad = ?, ptl1 = ?, ".
          "ptr1 = ?, ptl2 = ?, ptr2 = ? WHERE aln_id = ?" );

    my $clearhsp = $ml->{STH}{DELETE_HSP_BY_ID} ||= $dbh->prepare
        ( -name => "Clear HSPs for a single alignment",
          -sql  => "DELETE FROM hsp WHERE aln_id = ?");
    my $addhsp = $ml->{STH}{ADD_HSP} ||= $dbh->prepare
        ( -name => "Add a single HSP",
          -sql  => "INSERT INTO hsp (aln_id, ptl1, ptr1, ptl2) VALUES (?,?,?,?)" );
    $dbh->begin_work();
    $setUp->execute($pkey) unless 
        ($carefully && $check->get_single_value($pkey));
    $clearhsp->execute($pkey);
    $update->execute(@sids, $srcid, $str, $sc, $hb, @bnds, $pkey);
    $self->write_hsps();
    $dbh->commit();
    $self->bench_end();
}

=head2 

 Title   : delete
 Usage   : $aln->delete()
 Function: Clears the alignment from the database

This is a convienence wrapper for
BMS::SnpTracker::MapLoc::delete_alignment_by_id()

=cut

sub delete {
    my $self = shift;
    $self->bench_start();
    my $id = $self->pkey();
    return unless ($id);
    my $rv = $self->maploc->delete_alignment_by_id( $id );
    $self->bench_end();
    return $rv;
}

=head1 Reporting Methods

=head2 

 Title   : to_one_line
 Usage   : my $txt = $aln->to_one_line()
 Function: Get human-readable representation of the alignment
 Returns : A string describing the alignment
 Args    : [0] Optional "prefix" to go before the text

Will report something like:

   <pretext><Seq1> vs. <Seq2> <source> Alignment [<strand>] len=<length> Score <score> [HB <howbad>]

Pretext is intended to be optional padding (one or more white spaces).

=cut

sub to_one_line {
    my $self = shift;
    $self->bench_start();
    my $pre  = shift || "";
    my $str  = $self->str();
    my $sc   = $self->score();
    my @seqs = $self->seqs();
    my $txt  = sprintf("%s%s %s Alignment [%s] len=%s", $pre, 
                       join(' vs. ', @seqs), $self->source(),
                       !defined $str ? "" : $str > 0 ? '+1': $str < 0 ? '-1' :
                       'both strands', $self->length());
    if (defined $sc) {
        $txt .= " Score $sc";
        my $hb = $self->howbad();
        $txt .= $hb ? " [HB $hb]" : " [Best]" if (defined $hb)
    }
    $self->bench_end();
    return $txt;
}

=head2 

 Title   : to_text
 Usage   : my $txt = $aln->to_text()
 Function: Generate moderately detailed information on the alignment
 Returns : A block of text describing the alignment
 Args    : [0] Optional "prefix" to go before the text

This will include the output from to_one_line(), and will also include
extra lines describing the cooridinate bounds of the two sequences.

Additionally, the output from
BMS::SnpTracker::MapLoc::Tagged::tag_text() will be added, to indicate
any TagVal pairs assigned to the alignment.

=cut

sub to_text {
    my $self = shift;
    $self->bench_start();
    my $txt  = $self->to_one_line( @_ )."\n";
    my $pre  = shift || "";
    my @seqs = $self->seqs();
    my @pflk = $self->pretty_bounds();
    $txt .= sprintf("%s  [%-20s] %s\n", $pre, $pflk[0], $seqs[0]);
    $txt .= sprintf("%s  [%-20s] %s\n", $pre, $pflk[1], $seqs[1]);
    $txt .= $self->tag_text($pre."  ");
    $self->bench_end();
    return $txt;
}

=head2 

 Title   : full_text
 Usage   : my $txt = $aln->full_text()
 Function: Generate a full, human-readable report on the alignment
 Returns : A detailed text description of the alignment
 Args    : [0] Optional "prefix" to go before the text

This will generate the full block of text produced by to_text(), and
in addition will include an alignment generated by
BMS::SnpTracker::MapLoc::HSPed::HSPed::hsp_text().

=cut

sub full_text {
    my $self = shift;
    return $self->to_text() . $self->hsp_text();
}

=head1 Export Methods

=head2 

 Title   : json_data
 Usage   : my $json = $aln->json_data()
 Function: Get a JSON data structure suitable for web applications
 Returns : A Perl hash

Will return a hash structure with the following keys and values:

 * alnid : pkey()
 * source : source()
 * score : score()
 * howbad : howbad()
 * seqs : [ seqs() ]

=cut

sub json_data {
    my $self = shift;
    $self->bench_start();
    # $self->read();
    my @seqs = $self->seqs();
    my $rv = {
#        tags    => $self->all_tag_values(),
        alnid   => $self->pkey(),
        source  => $self->source(),
        score   => $self->score(),
        howbad  => $self->howbad(),
        seqs    => \@seqs,
    };
    $self->bench_end();
    return $rv;
}

=head2 

 Title   : cx_genome_part
 Usage   : my $cxHash = $aln->cx_genome_part( @args )
 Function: Generate a CanvasXpress data structure for the alignment
 Returns : In array context, an array of Perl hashes
           In scalar context, just the first hash
 Args    : See below

This method generates a data structure that can be interpreted by the
GenomePanel in
CanvasXpress. BMS::SnpTracker::MapLoc::HSPed::HSPed::ordered_coordinates()
is used to parse the arguments to determine which member of the alignment pair is to be used as the "anchor" sequence.

Because Alignment objects have only two sequences, the array should
always just have one object in it.

The data structure will be built from the HSP CX object generated by
BMS::SnpTracker::MapLoc::HSPed::HSPed::cx_genome_part(), and will in
addition include the following keys:

 * score : score()
 * howbad : howbad()
 * source : source()
 * aln_id : pkey()
 * fill : //fill color, currently fixed at gray//
 * outline : outline color, based on HowBad

=cut

*to_canvasXpress = \&cx_genome_part;
sub cx_genome_part {
    my $self  = shift;
    $self->bench_start();
    # Get generic ::Alignment structures:
    my @parts = $self->SUPER::cx_genome_part( @_ );
    my $sc    = $self->score();
    my $hb    = $self->howbad();
    my $src   = $self->source();
    my $pkey  = $self->pkey();
    my $ml    = $self->maploc();
    my $pSc   = $ml->poor_score();
    my $gSc   = $ml->good_score();
    if (0) {
        my $cval  = 255;
        if (defined $sc) {
            $cval  = int(0.5 + 255 * ($sc - $pSc) / ($gSc - $pSc));
            if ($cval < 0) {
                $cval = 0;
            } elsif ($cval > 255) {
                $cval = 255;
            }
        }
        my $col = sprintf("rgb(%d,%d,255)",255 - $cval, $cval);
    }
    my $col  = "#cccccc";
    foreach my $cx (@parts) {
        $cx->{score}   = $sc;
        $cx->{howbad}  = $hb;
        $cx->{source}  = $src;
        $cx->{aln_id}  = $pkey;
        $cx->{fill}    = $col;
        $cx->{outline} = $hb < 1 ? "#006600" : $hb < 3 ? "#ff9900" : "#ff0000";
    }
    $self->bench_end();
    return wantarray ? @parts : $parts[0];
}

