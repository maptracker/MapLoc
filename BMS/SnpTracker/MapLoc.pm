# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::SnpTracker::MapLoc;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

=head1 DESCRIPTION

This is the primary module used to access the MapLoc data
warehouse. It includes methods for retrieving and adding information
from the database.

=head1 SYNOPSIS

 use BMS::SnpTracker::MapLoc;
 my $ml = BMS::SnpTracker::MapLoc->new();

Recovering information:

 # Get a location object
 my $loc = $ml->get_location("16", 227331, 227333, "GRCh37");

 # Get a population object
 my $pop = $ml->get_population("1kG.JPT")

 # Get RNA objects
 my @rnas = $ml->get_rnas( "NM_001234" );

 # Supporting information about an impact code
 my $infoHash = $ml->impact_details("SYN");

The above calls are relatively low-level. For handling human-readable
queries, see classify_requests().

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
use BMS::Utilities::Serialize;
use BMS::FriendlyDBI;
use BMS::Utilities::Benchmark;
use BMS::Utilities::SequenceUtilities;
use BMS::Utilities::ColorUtilities;
use BMS::ErrorInterceptor;
use BMS::ArgumentParser;

use BMS::SnpTracker::MapLoc::Location;
use BMS::SnpTracker::MapLoc::RNA;
use BMS::SnpTracker::MapLoc::Population;
use BMS::SnpTracker::MapLoc::Text;
use BMS::SnpTracker::MapLoc::Gene;
use BMS::SnpTracker::MapLoc::Range;
use BMS::SnpTracker::MapLoc::Alignment;
use BMS::SnpTracker::MapLoc::TransientAlignment;

$BMS::SnpTracker::MapLoc::VERSION = 
    ' $Id: MapLoc.pm,v 1.18 2014/03/21 16:06:10 tilfordc Exp $ ';

use vars qw(@ISA);
@ISA      = qw(BMS::Utilities::Escape
               BMS::ErrorInterceptor
               BMS::Utilities::FileUtilities
               BMS::Utilities::ColorUtilities
               BMS::Utilities::SequenceUtilities
               BMS::Utilities::Benchmark);


our $iLoop = {};
our $instanceName;
our $fsTok = "FrameShift";

our $impactRanks = {
    CPX => {
        rank => 1,
        name => "Complex",
        color => "#ffff00",
        fgcol => "#ff0000",
        desc  => "The variant is 'complicated', hindering attempts to electronically predict impact. Generally this involves multi-base regions covering an exon plus intronic or genomic sequence",
    },
    FRM => {
        rank => 1.1,
        name => "Frameshift",
        color => "#ffcc33",
        fgcol => "#ff0000",
        desc  => "At least part of the variant is within the CDS, and at least one allele does not have a 'modulus three' number of bases relative to at least one other allele",
    },
    STP => {
        rank => 2,
        name => "Stop",
        color => "#ff0000",
        fgcol => "#ffff00",
        desc  => "At least one allele has a stop codon, and at least one other allele does not",
        alias => "STOP"
    },
    DEL => {
        rank => 3,
        name => "Deletion",
        color => "#000000",
        fgcol => "#ffffff",
        desc  => "The alleles code for differing protein lengths WITHOUT disrupting (frameshifting or introduction of a new stop) the reading frame",
    },
    NON => {
        rank => 4,
        name => "Non-synonymous",
        color => "#ff6600",
        fgcol => "#000000",
        desc  => "The alleles code for different proteins of the SAME length",
    },
    SP5 => {
        rank => 4.9,
        name => "Splice Donor",
        color => "#ff33ff",
        fgcol => "#ffff00",
        desc  => "The location is within the two 'GT' 5-prime splice donor bases",
        alias => "SP-5 5-SP 5SP 5'SP SPL5 SPL-5 5-SPL 5SPL 5'SPL"
    },
    SP3 => {
        rank => 4.9,
        name => "Splice Acceptor",
        color => "#ff33ff",
        fgcol => "#ffff00",
        desc  => "The location is within the two 'AG' 3-prime splice acceptor bases",
        alias => "SP-3 3-SP 3SP 3'SP SPL3 SPL-3 3-SPL 3SPL 3'SPL"
    },
    SPL => {
        rank => 5,
        name => "Splice Site",
        color => "#ff33ff",
        fgcol => "#ffff00",
        desc  => "The location is within a splice site. Information has not been provided if the site is an acceptor or donor.",
        kids  => [ qw(SP3 SP5) ],
    },
    SYN => {
        rank => 6,
        name => "Synonymous",
        color => "#3333ff",
        fgcol => "#ffff00",
        desc  => "The location is within the CDS of the transcript, but all alleles code for the same amino acid content",
    },
    COD => {
        rank => 6.5,
        name => "Coding",
        color => "#999999",
        fgcol => "#000000",
        desc  => "The location falls within the CDS of the transcript, but there are not enough reported alleles to calculate synonymy. This could be because of filters that exclude some alleles, or because the site is bogus (ie not variant)",
        kids  => [ qw(FRM STP DEL NON SYN) ],
    },
    UT5 => {
        rank => 6.9,
        name => "UTR-5",
        color => "#ff99ff",
        fgcol => "#cc3300",
        desc  => "The location is within an exon of the transcript in the 5-prime UnTraslated Region",
        alias => "UTR5 UTR-5 5-UTR 5UTR 5'UTR"
    },
    UT3 => {
        rank => 6.9,
        name => "UTR-3",
        color => "#ff99ff",
        fgcol => "#cc3300",
        desc  => "The location is within an exon of the transcript in the 3-prime UnTraslated Region",
        alias => "UTR3 UTR-3 3-UTR 3UTR 3'UTR"
    },
    UTR => {
        rank => 7,
        name => "UTR",
        color => "#ff99ff",
        fgcol => "#cc3300",
        desc  => "The location is within either the 3- or 5-prime UnTraslated Region",
        kids  => [ qw(UT3 UT5) ],
    },
    NCR => {
        rank => 7.1,
        name => "Non-coding RNA",
        color => "#990099",
        fgcol => "#ff9933",
        desc  => "Falls within an exon of a transcript that does not have a CDS defined",
        kids  => [ qw(SYN COD UT3 UT5 UTR) ],
    },
    RNA => {
        rank => 7.3,
        name => "mRNA",
        color => "#990099",
        fgcol => "#ff9933",
        desc  => "The location is within an exon but no other information was available",
        kids  => [ qw(FRM STP DEL NON SYN COD UT3 UT5 UTR) ],
    },
    INT => {
        rank => 8,
        name => "Intronic",
        color => "#cccc33",
        fgcol => "#000000",
        desc  => "The location is within an intron of the transcript, but not within a splice site",
    },
    LOC => {
        rank => 8.5,
        name => "Locus",
        color => "#cccc33",
        fgcol => "#000000",
        desc  => "The location falls within the pre-mRNA of a locus, but no information has been provided if it is intronic or exonic",
        kids  => [ qw(FRM STP DEL NON SYN COD UT3 UT5 UTR SP3 SP5 SPL UT3 UT5 UTR) ],
    },
    GEN => {
        rank => 9,
        name => "Genomic",
        color => "#cccccc",
        fgcol => "#660066",
        desc  => "The location is within a region with no reported transcripts",
    },
    UNK => {
        rank => 10,
        name => "Unknown",
        color => "#9999cc",
        fgcol => "#ff9933",
        desc  => "The impact of the location could not be calculated",
    },
};


my $impactAlias = {};
my $genericImp  = {};
while (my ($tk, $dat) = each %{$impactRanks}) {
    $dat->{token} = $tk;
    if (my $ali = $dat->{alias}) {
        map { $impactAlias->{$_} = $tk } split(/\s+/, $ali);
    }
    if (my $kids = $dat->{kids}) {
        $genericImp->{$tk} = $kids;
    }
}
# map { $impactRanks->{$_}{token} = $_; } keys %{$impactRanks};

our @allImpacts = map { $_->{token} } sort { 
    $a->{rank} <=> $b->{rank} || $a->{token} cmp $b->{token}
}  values %{$impactRanks};

=head1 Primary Object Instantiation Methods

These basic methods allow specific MapLoc objects to be recovered from
the database.

=head2 new

 Title   : new
 Usage   : my $ml = BMS::SnpTracker::MapLoc->new( @params )
 Function: Create a new MapLoc object
 Returns : A blessed MapLoc object
 Args    : Parameter array, see below

This is the primary API entry point to the MapLoc system. To use
MapLoc within your code, create a single object and use it for
generating all queries or additional objects you may need. It will be
shared by reference with the other objects generated by the
system. When the variable "$ml" is used in the documentation, it is
referencing an object generated by this method.

Currently two parameters are recognized:

  -instance The database instance name. The value can be defined in
            the MapLoc.param parameter file, provided here, or it will
            revert to 'maploc' by default.

     -build Optional. The default genome build, will be set with
            the build() method.

=cut

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self = {
        CACHE_SIZE   => 1000000,
        CACHE_PURGES => 0,
        DEF_HOWBAD   => 0,
        LOC2RNA_DIST => 10000,
        VERBOSITY    => 1,
    };
    bless ($self, $class);
    $self->bench_start();

    
    my $args = $self->parseparams( @_ );
    if (my $inm = $args->{INSTANCE}) {
        $instanceName = $inm;
    }
    if (my $build = $args->{BUILD}) {
        $self->build( $build );
    } else {
        # $self->death("You must provide -build when creating a new object");
    }
    # $self->dbh->no_environment( -noenv => $args->{NOENV} );
    $self->bench_end();
    return $self;
}

=head2 DESTROY

 Title   : DESTROY
 Usage   : Called automatically on object destruction
 Function: Cleanly destroy the object when it goes out of scope

This is the special Perl function that handles object release during
garbage collection. In this case it assures that the database handle
is cleanly disconnected.

=cut

sub DESTROY {
    my $self = shift;
    return unless ($self);
    $self->disconnect();
    $self->SUPER::DESTROY( @_ );
}

=head2 get_location

 Title   : get_location
 Usage   : my $locationObject = $ml->get_location( @args )
 Function: Get a single Location object
 Returns : A BMS::SnpTracker::MapLoc::Location object, or undef
 Args    : [0] Chromosome / contig identifier ... or loc_id
           [1] LeftFlank position
           [2] RightFlank position
           [3] GenomeBuild

Primary mechanism to recover Location objects.

If a single parameter is passed it will be presumed to be an internal
loc_id database identifier.

Otherwise, four arguments will be parsed as the chromosome, left
flank, right flank and genome build. he chromosome name and build will
be standardized using standardize_chromosome(). If the build is not
provided or extracted, the default_build() will be used.

=cut

sub get_location {
    my $self = shift;
    # Treat a single parameter as a loc_id:
    return BMS::SnpTracker::MapLoc::Location->new($self, @_ ) if ($#_ == 0);
    # Otherwise, treat the parameters as flanks. We will still tidy up
    # the chromosome, however:
    my ($chr, $lft, $rgt, $build) = @_;
    my ($niceChr, $foundBuild) = $self->standardize_chromosome( $chr );
    return BMS::SnpTracker::MapLoc::Location->new
        ($self, $niceChr, $lft, $rgt, $build || $foundBuild || $self->{BUILD});
}

=head2 get_location_by_startend

 Title   : get_location_by_startend
 Usage   : my $locObj = $ml->get_location_by_startend($chr, $start, $end)
 Function: Gets a Location object for a chromosome and start/end coordinates
 Returns : A BMS::SnpTracker::MapLoc::Location object
 Args    : [0] Chromosome designation
           [1] Start coordinate
           [2] End coordinate

Essentially the same as get_location(), but will expect start/end coordinates

=cut

*request_to_location = \&get_location_by_startend;
sub get_location_by_startend {
    # This method should be used for human-readable positions
    # That is, if coordinates are passed, they should NOT be flanks
    my $self = shift;
    my @args = $self->request_to_args( @_ );
    return BMS::SnpTracker::MapLoc::Location->new( $self, @args );
}

=head2 get_population_fast

 Title   : get_population_fast
 Usage   : my $popObj = get_population_fast( @opts )
 Function: Gets a population object based on ID or name and source
 Returns : A BMS::SnpTracker::MapLoc::Population object
 Args    : [0] A pop_id (integer), or a population name (string)
           [1] A source (string) when a name is provided
           [2] A flag to indicate that an integer first parameter is
               a name, not a database ID

This method will recover a population, but will not call read() on it.

If the first parameter is an ID and the third value is false, it will
be treated as a database ID. Otherwise it will be treated as a
population name. If a unique population if found with a name, it will
be returned. Otherwise, the population source must match the value
provided in the second option.

=cut

sub get_population_fast {
    my $self = shift;
    my $obj = $self->get_cache("Population", @_);
    if (!defined $obj) {
        $self->bench_start();
        # Set the cache to zero initially - prevents recursion in some cases
        $self->set_cache(0, "Population", @_);
        $obj = BMS::SnpTracker::MapLoc::Population->new($self, @_) || 0;
        $self->set_cache($obj ||= 0, "Population", @_);
        $self->bench_end();
    }
    return $obj;
}

=head2 get_population

 Title   : get_population
 Usage   : my $popObj = get_population( @opts )
 Function: Gets a population object based on ID or name and source,
           and reads associated metadata
 Returns : A BMS::SnpTracker::MapLoc::Population object
 Args    : [0] A pop_id (integer), or a population name (string)
           [1] A source (string) when a name is provided
           [2] A flag to indicate that an integer first parameter is
               a name, not a database ID

Identical to get_population_fast(), but additionally calls read() on
the returned object.

=cut

sub get_population {
    my $obj  = shift->get_population_fast( @_ );
    $obj->read() if ($obj);
    return $obj;
}

=head2 get_text

 Title   : get_text
 Usage   : my $txtObj = $ml->get_text( $text )
 Function: Gets a Text object for a request
 Returns : A BMS::SnpTracker::MapLoc::Text object
 Args    : [0] A string representing the text involved

Used to recover text objects. Because a Text object can have arbitrary
key/value tag pairs, and the values of the tags can themselves be text
objects, they are often used to represent complex structures.

=cut

sub get_text {
    my $self = shift;
    my $obj = $self->get_cache("Text", @_);
    if (!defined $obj) {
        $self->bench_start();
        $obj = BMS::SnpTracker::MapLoc::Text->new($self, @_);
        $self->set_cache($obj ||= 0, "Text", @_);
        $self->bench_end();
    }
    return $obj;    
}

=head2 get_gene

 Title   : get_gene
 Usage   : my $geneObj = $ml->get_gene( $accession )
 Function: Recover a MapLoc Gene object by its accession
 Returns : A BMS::SnpTracker::MapLoc::Gene object
 Args    : [0] The accession for the gene

Genes are just Text objects with additional functions assigned to
them. They are represented solely by their accession.

=cut

sub get_gene {
    my $self = shift;
    my $obj = $self->get_cache("Gene", @_);
    if (!defined $obj) {
        $self->bench_start();
        $obj = BMS::SnpTracker::MapLoc::Gene->new($self, @_);
        $self->set_cache($obj ||= 0, "Gene", @_);
        $self->bench_end();
    }
    return $obj;    
}

=head2 get_gene_by_id

 Title   : get_gene_by_id
 Usage   : my $geneObj = $ml->get_gene_by_id( $pkey )
 Function: Get a Gene object by its internal database ID
 Returns : A BMS::SnpTracker::MapLoc::Gene object
 Args    : [0] The txt_id PKEY value for the gene accession

Same as get_gene(), but uses the accession PKEY as the query

=cut

sub get_gene_by_id {
    my $self = shift;
    my $obj = $self->get_cache("Gene", @_);
    if (!defined $obj) {
        $self->bench_start();
        my $txt = $self->pkey_to_text( shift @_ );
        $obj = BMS::SnpTracker::MapLoc::Gene->new($self, $txt, @_);
        $self->set_cache($obj ||= 0, "Gene", @_);
        $self->bench_end();
    }
    return $obj;    
}

=head2 get_alignment

 Title   : get_alignment
 Usage   : my $alnObj = $ml->get_alignment( $pkey )
 Function: Get an Alignment object for a database PKEY id
 Returns : A BMS::SnpTracker::MapLoc::Alignment object
 Args    : [0] The aln_id PKEY value for the alignment

Recovers an Alignment object from its database aln_id PKEY

=cut

sub get_alignment {
    my $self = shift;
    if ($#_ == -1) {
        # A request for a 'blank' alignment
        return BMS::SnpTracker::MapLoc::Alignment->new($self);
    }
    # At least one query parameter has been provided
    my $obj = $self->get_cache("Alignment", @_);
    if (!defined $obj) {
        $self->bench_start();
        $obj = BMS::SnpTracker::MapLoc::Alignment->new($self, @_) || 0;
        $self->set_cache($obj, "Alignment", @_);
        $self->bench_end();
    }
    return $obj;
}

=head2 get_transient_alignment

 Title   : get_transient_alignment
 Usage   : my $transAlnObj = $ml->get_transient_alignment( )
 Function: Creates an empty TransientAlignment
 Returns : A BMS::SnpTracker::MapLoc::TransientAlignment object

Convienence method that simply generates a new object

=cut

sub get_transient_alignment {
    my $self = shift;
    return BMS::SnpTracker::MapLoc::TransientAlignment->new($self);
}

=head2 get_rna_by_id

 Title   : get_rna_by_id
 Usage   : my $rnaObj = $ml->get_rna_by_id( $pkey )
 Function: Get an RNA object by its database PKEY id
 Returns : An BMS::SnpTracker::MapLoc::RNA object
 Args    : [0] the rna_id database PKEY id

Uses the internal rna_id PKEY value to recover an RNA object

=cut

sub get_rna_by_id {
    my $self = shift;
    my $id = shift;
    my $obj = $self->get_cache("RNAID", $id);
    if (!defined $obj) {
        $self->bench_start();
        $obj = BMS::SnpTracker::MapLoc::RNA->new($self, $id) || 0;
        $self->set_cache($obj, "RNAID", $id);
        $self->bench_end();
    }
    return $obj;
}

=head2 get_rnas

 Title   : get_rnas
 Usage   : my @rnaObjs = $ml->get_rnas( $acc, $src )
 Function: Gets one or more RNA objects based on accession and source
 Returns : An array of BMS::SnpTracker::MapLoc::RNA objects
 Args    : [0] A string representing the accession
           [1] A string representing the source for the RNA

The accession can include wildcards. This allows queries such as
"NM_001234.%", which will recover all versions of a particular
accession.

If called in scalar context, then only one of the recovered RNAs will
be returned, effectively chosen at random. This is probably stupid
behaivor. Calling in array context is advised.

=cut

sub get_rnas {
    my $self = shift;
    my $hits = $self->get_cache("RNA", @_);
    unless ($hits) {
        $self->bench_start();
        $hits    = $self->set_cache( [], "RNA", @_);
        # First parameter is accession, second is source
        my @accIds = $self->text_search( $_[0], 'wc');
        if ($#accIds == -1) {
            $self->bench_end();
            return wantarray ? () : undef;
        }
        my @aids = $self->find_text_ids( $_[1] );
        my ($sth, @binds, $rtab);
        if ($#accIds > 10 || $#aids > 0) {
            # do query in bulk off of temporary table
            my $dbh  = $self->dbh;
            $rtab    = $dbh->list_to_temp_table(\@accIds, 'integer');
            my $sql  = "SELECT r.rna_id FROM rna r WHERE r.acc_id ".
                $self->_subquery_clause( 'col1', $rtab);
            unless ($#aids == -1) {
                my $atab = $dbh->list_to_temp_table(\@aids, 'integer');
                $sql    .= " AND r.src_id ".
                    $self->_subquery_clause( 'col1', $atab);
            }
            $sth = $dbh->prepare
                ( -name => "Get RNAs for accession in temp table",
                  -sql  => $sql );
            @binds = ([]);
        } else {
            # Small number of accids
            if ($#aids == -1) {
                # no source defined
                $sth = $self->{STH}{RNA_FOR_1ACC} ||= $self->dbh->prepare
                    ( -name => "Get RNAs for single accession",
                      -sql  => "SELECT r.rna_id FROM rna r WHERE ".
                      "r.acc_id = ?");
                @binds = map { [$_] } @accIds;
            } else {
                # Source defined
                $sth = $self->{STH}{RNA_FOR_1ACC_1SRC} ||= $self->dbh->prepare
                    ( -name => "Get RNAs for single accession and source",
                      -sql  => "SELECT r.rna_id FROM rna r WHERE ".
                      "r.acc_id = ? AND r.src_id = ?");
                foreach my $aid (@aids) {
                    push @binds, map { [$_, $aid] } @accIds;
                }
            }
        }
        my %rids;
        foreach my $bdat (@binds) {
            map { $rids{$_} = 1 } $sth->get_array_for_field( @{$bdat} );
        }
        $self->dbh->clear_temp_table( $rtab ) if ($rtab);
        foreach my $rid (keys %rids) {
            if (my $rna = $self->get_rna_by_id( $rid )) {
                push @{$hits}, $rna;
            }
        }
        $self->bench_end();
    }
    return wantarray ? @{$hits} : $hits->[0];
}

=head2 get_rnas_for_locus

 Title   : get_rnas_for_locus
 Usage   : my @rnaObjs = $ml->get_rnas_for_locus( $acc )
 Function: Get a list of RNA objects for a gene accession
 Returns : An array of BMS::SnpTracker::MapLoc::RNA objects, or an
           array reference if called in scalar context
 Args    : [0] a string representing the locus accession

RNAs are associated witht their loci by tag/value pairs. This method
uses get_tagged_object_ids() to recover RNA object ids, then calls
get_rna_by_id() to get the RNA objects.

=cut

our $notCurrentTag = 'NotCurrent';
sub get_rnas_for_locus {
    my $self = shift;
    my @rv;
    if (my $locid = shift) {
        my ($tag, $re) = $locid =~ /ENST/ ? 
            ('ENSG', '^ENST') : ('LL', '^[NX][MR]_');
        foreach my $rid ($self->get_tagged_object_ids($locid, $tag)) {
            if (my $rna = $self->get_rna_by_id( $rid )) {
                my $acc = $rna->accession();
                next unless ($acc =~ /$re/);
                $rna->read();
                next if ($rna->simple_value($notCurrentTag));
                push @rv, $rna;
            }
        }
    }
    return wantarray ? @rv : \@rv;
}

=head2 string_to_object

 Title   : string_to_object
 Usage   : my @perlObjects = $ml->string_to_object( $text, $caseSens )
 Function: Turn free text (presumably user queries) into blessed objects
 Returns : An array of MapLoc objects, or array reference if a scalar
           was requested
 Args    : [0] A string
           [1] A flag, indicating that case sensitive matching should
               be used

Used primarily by classify_requests(). Will use text_search() to see
if the string matches any strings stored in the database, then will
use text_id_to_objects() to recover MapLoc Perl objects associated
with the string.

The string can be prefixed with "#NAMESPACE#", where NAMESPACE is text
restricting the search to a particular type of object: see
text_id_to_objects().

=cut

sub string_to_object {
    my $self = shift;
    my ($txt, $isCS) = @_;
    my @txtObjs;
    if ($txt) {
        my $nsPrfx;
        if ($txt =~ /^\#([a-zA-Z]+)\#(.+)/) {
            ($nsPrfx, $txt) = ($1, $2);
        }

        # Has the string been seen?
        my $texts = $self->text_search($txt, undef, $isCS);
        foreach my $row (@{$texts}) {
            my ($tid, $txt) = @{$row};
            my @objs = $self->text_id_to_objects( $tid, $nsPrfx );
            push @txtObjs, @objs;
        }
    }
    return wantarray ? @txtObjs : \@txtObjs;
}

=head2 text_id_to_objects

 Title   : text_id_to_objects
 Usage   : my @perlObjects = $ml->text_id_to_objects( $txtID, $namespace )
 Function: Turn a normtxt txt_id into one or more Perl objects
 Returns : An array of MapLoc objects in array context. In scalar context,
           a hash reference with keys being object PKEYs and values the
           objects themselves
 Args    : [0] An integer, should be a txt_id in table normtxt
           [1] Optional namespace

The default is to find ALL objects associated with the text ID. If a
namespace is provided, then it will restrict the search to just the
requested kind of objects. The methods used [and their recognized
namespace flags] are:

 * Location   : acc_id_to_loc_ids()         [location, var, snp]
 * RNA        : txt_id_to_rna_ids()         [rna, trans]
 * Gene       : txt_id_to_gene_ids()        [gene, locus]
 * Population : txt_id_to_population_ids()  [pop, patient]

Multiple namespaces can be combined, so using "genepop" would allow
matches to both Gene and Population objects.

=cut

sub text_id_to_objects {
    my $self = shift;
    my $id   = shift;
    my $ns   = lc(shift || "");
    my %rvids;
    return wantarray ? () : \%rvids if (!$id || $id !~ /^\d+$/);
    if (!$ns || $ns =~ /(location|var|snp)/) {
        foreach my $lid ($self->acc_id_to_loc_ids( $id )) {
            if (my $loc = $self->get_location( $lid )) {
                $rvids{$lid} ||= $loc;
            }
        }
    }
    if (!$ns || $ns =~ /(rna|trans)/) {
        foreach my $rid ($self->txt_id_to_rna_ids( $id )) {
            if (my $rna = $self->get_rna_by_id( $rid )) {
                $rvids{$rid} ||= $rna;
            }
        }
    }
    if (!$ns || $ns =~ /(gene|locus)/) {
        foreach my $gid ($self->txt_id_to_gene_ids( $id )) {
            if (my $gene = $self->get_gene_by_id( $gid )) {
                $rvids{$gid} ||= $gene;
            }
        }
    }
    if (!$ns || $ns =~ /(pop|patient)/) {
        foreach my $pid ($self->txt_id_to_population_ids( $id )) {
            if (my $pop = $self->get_population_fast( $pid )) {
                $rvids{$pid} ||= $pop;
            }
        }
        # $self->prebranch(\%rvids);
    }
    return wantarray ? values %rvids : \%rvids;
    # Population
}

=head1 Primary Search Methods

These get methods listed above allow objects to be created once their name or ID are known. The search methods below allow more complex searches to be performed using related inforamtion, or an intersection of several criteria.

=head2  classify_requests

 Title   : classify_requests
 Usage   : my $parsed = $maploc->classify_requests( @param );
 Function: Parses user input into recognized MapLoc objects
 Returns : A hash reference.
 Args    : Associative array of parameters

This method will look at a human-provied query and try to "make sense"
of it. It will break the request into one or more recognized objects,
returned as a structured hash. The hash is a nested structure:

  Lvl 1: Keyed to object type (simple string) and has hash reference
         as values
     Lvl 2: Keyed to internal database ID (pkey) and has array
            references as values
        Lvl 3: Arrays are structured as:
               [0] an array reference of user query terms
               [1] a Yes/No flag, 1 = yes, -1 = no
               [2] Perl object based on pkey (not always present)
               [3] optional tag/value pairs provided by user

The recognized query parameters are:

    -request Required. The primary request. Can be a single string, or
             an array reference of strings. Also recognized as -req

  -sensitive Default 0. If true, then the query should be treated as
             case-sensitive. This may be desired for gene
             symbols. Also recognized as -casesensitive

     -howbad Default defined by default_howbad(). Optional HowBad value
             for RNA results.

      -build Default defined by default_build(). Optional GenomeBuild

Each request can be more complex than a simple string. Formal Perl
objects can be provided:

   * Simple arrays will be expanded and each element further analyzed
   * MapLoc objects will be returned as is, under their obj_type()

For string requests, if the request begins with an exclamation point
(!), then it is presumed that the request is to be excluded, and the
Yes/No flag will be set to -1.

If the request contains anything that looks like:

   * TAG=VALTEXT
   * TAG="VAL TEXT"
   * TAG='VAL TEXT'

... then that part will be removed, and the Tag/Val pair will be
associated with the returned data. "TAG" must not contain spaces, VAL
can contain spaces if inside quotes.

Additional alleles can be associated with location by using a dollar
sign. For example, "3.GCRh37:1999432$A/T" will associate the alleles A
and T with the location. Using the format "$A>T" will set A to be the
reference allele. Using two dollar signs (eg "rs12345$$GG/C/T") will
force the location to be a TransientLocation, even if it was a known
location within the database.

Blocks of text will first be broken on new lines and returns. If a
whole line can be matched to an object, that match will be
kept. Otherwise, an attempt will be made to break the line on white
space, and the parts will be matched independently.

Locations are identified using parse_genomic_location_text(), or if
the request is "loc_id=INTEGER" then it will be treated as the PKEY
for a Location.

If the request is not deemed to be a location, it will then be
analyzed with string_to_object().

=cut

sub classify_requests {
    my $self = shift;
    unshift @_, "-req" if ($#_ == 0);
    my $args = $self->parseparams( @_ );
    my $rv = { };
    my $req = $args->{REQ} || $args->{REQUEST};
    my $isCS = $args->{CASESENSITIVE} || $args->{SENSITIVE};
    my @vals;
    if (!$req) {
        my $rt = !defined $req ? 'undefined' : $req eq '' ?
            'empty string' : $req;
        $rv->{'Error'} = ["No request ($rt) provided"];
    } elsif (my $r = ref($req)) {
        @vals = $r eq 'ARRAY' ? @{$req} : ($req);
    } else {
        @vals = ($req);
    }
    return $rv if ($#vals == -1);
    for (my $v = 0; $v <= $#vals; $v++) {
        my $txt = $vals[$v];
        # Still need to check since we may be expanding nested arrays...
        next if (!defined $txt);
        if (my $r = ref($txt)) {
            # An object is being passed
            if ($r eq 'ARRAY') {
                # Maybe I will regret this, but we will put nested
                # arrays back onto the stack
                push @vals, @{$r};
            } elsif ($r =~ /::MapLoc::/ 
                     && $r->can('obj_type') && $r->can('handle')) {
                # This is already a MapLoc object
                my $ot = $txt->obj_type();
                my $targ = $rv->{$ot}{$txt->handle()} ||= [ [], 0 ];
                push @{$targ->[0]}, 'Perl Object';
                $targ->[2] ||= $txt;
            } else {
                # Huh?
                push @{$rv->{'Error'}}, ["Unrecognized request '$txt' [Perl Object $r]"];
            }
            next;
        }
        # Ok, it is a scalar - can we parse it? 
        foreach my $val (split(/[\n\r]+/, $txt)) {
            $val =~ s/^\s+//; $val =~ s/\s+$//;
            my $yn = 1; # YesNo : 1 = include, -1 = exclude
            if ($val =~ /^!\s*(.+)$/) {
                $val = $1;
                $yn = -1;
            }
            my $tagVal;
            while ($val =~ /(\s*\/(\S+)=\"([^\"]+)\"\s*)/ ||
                   $val =~ /(\s*\/(\S+)=\'([^\']+)\'\s*)/ ||
                   $val =~ /(\s*\/(\S+)=(\S+)\s*)/) {
                # /tag="some value" /tag='some value' /tag=someValue
                my ($rep, $t, $v) = ($1, $2, $3);
                $v =~ s/^\s+//;
                $v =~ s/\s+$//;
                $tagVal ||= {};
                $tagVal->{$t}{$v}++;
                $val =~ s/\Q$rep\E/ /;
            }
            if ($tagVal) {
                $val =~ s/^\s+//; $val =~ s/\s+$//;
            }
            next unless ($val);
            my ($xtraVar, $forceTrans);
            if ($val =~ /(.+?)(\$+)(.+)$/) {
                ($val, $forceTrans, $xtraVar) = ($1, $2, $3);
                # warn "($val, $forceTrans, $xtraVar)";
                # More than one $ means force this as a transient location
                $forceTrans = 0 if ($forceTrans eq '$');
            }
            if ($val =~ /^loc_id\s*=\s*(\d+)$/) {
                # Raw locus ID
                $rv ||= {};
                my $pkey = $1;
                my $key  = $forceTrans ? 'Range' : 'Location';
                my $targ = $rv->{$key}{$pkey} ||= [ [], -1 ];
                if ($xtraVar) {
                    # Need to add on extra alleles
                    if (my $loc = $targ->[2] ||= $self->get_location($pkey)) {
                        $loc->extra_alleles( $xtraVar );
                        $loc->lock() if ($forceTrans);
                    } else {
                        push @{$rv->{'Error'}}, "Unknown loc_id = $pkey";
                    }
                }
                &_extend_request_object
                    ( $targ, [$val], $yn, undef, $tagVal);
                next;
            }
            my ($chr, $lft, $rgt, $bld) = 
                $self->parse_genomic_location_text( $val );
            if ($chr) {
                # Coordinate specification
                $rv ||= {};
                push @{$rv->{'genome'}}, 
                [ [$val], $yn, $tagVal, [$chr, $lft, $rgt, $bld],
                  $xtraVar, undef, $forceTrans  ];
                next;
            }
            # Consider the line as a whole, and the line broken into parts
            my ($considered, $hasParts, @unk, @parts) = (0);
            if ($xtraVar) {
                # When extra alleles are specified, do not split the line
                @parts = ([$val, $xtraVar, $forceTrans]);
            } else {
                @parts = map { [$_] } split(/[\s\,]+/, $val);
                if ($#parts > 0) {
                    # Prepend the whole line
                    unshift @parts, [$val];
                    $hasParts = 1;
                }
            }
            for my $vi (0..$#parts) {
                my ($vv, $xv, $ft) = @{$parts[$vi]};
                next unless ($vv);
                my $qVal = $xv ? $vv.($ft ? '$$' : '$').$xv : $vv;
                $considered++;
                my $txtObjs = $self->string_to_object( $vv, $isCS );
                if ($#{$txtObjs} == -1) {
                    # Huh. Could not figure it out
                    push @unk, $qVal;
                    next;
                }
                # We were able to map the query over to a formal object
                my @stack  = @{$txtObjs};
                while ($#stack != -1) {
                    # At some point I was dynamically shifting new objects
                    # back onto the stack, but not anymore.
                    # Leaving the while/shift alone for the moment.
                    my $obj = shift @stack;
                    my %tvCopy = %{$tagVal || {}};
                    $tvCopy{"User Query"}{$qVal} = 1;
                    if ($xv && $obj->can('derived_locations')) {
                        # We actually want to map these to locations
                        my @dl = $obj->derived_locations
                            ( -loc     => $xv,
                              -howbad  => $args->{HOWBAD},
                              -build   => $args->{BUILD},
                              -verbose => 1);
                        next if ($#dl == -1);
                        if ($obj->can('acc')) {
                            $tvCopy{"Derived From"}{$obj->acc()} = 1;
                        }
                        $tvCopy{"User Location"}{$xv} = 1;
                        foreach my $dat (@dl) {
                            my ($ldat, $all, $str) = @{$dat};
                            push @{$rv->{'genome'}}, 
                            [[$qVal], $yn, \%tvCopy, $ldat, $all,$str,$ft];
                        }
                        next;
                    }
                    $rv ||= {};
                    my $targ = $rv->{$obj->obj_type()}{$obj->handle()} ||=
                        [ [], -1 ];
                    &_extend_request_object
                        ( $targ, [$vv], $yn, $obj, \%tvCopy);
                    if ($xv && $obj->can('extra_alleles')) {
                        $obj->extra_alleles($xv);
                        # print $obj->to_html();
                    }
                    # $obj->read(); print $obj->to_html();
                }
                if (!$vi && $hasParts) {
                    # If we could match the whole line by itself then
                    # do not consider the parts of the line broken up
                    last;
                }
            }
            next if ($#unk == -1);
            if ($#unk + 1 == $considered) {
                # Whole line was unintelligable
                my $qVal = $xtraVar ? 
                    $val.($forceTrans ? '$$' : '$').$xtraVar : $val;
                @unk = ($qVal);
            }
            foreach my $u (@unk) {
                $rv ||= {};
                my $targ = $rv->{'Unknown'}{$u} ||= [ [$u], -1 ];
                &_extend_request_object
                    ( $targ, undef, $yn, undef, $tagVal);
            }
        }
    }
    return $rv unless ($rv);
    if (exists $rv->{'genome'}) {
        # Map genomic coordinates over to formal Location objects
        my @keep;
        my %transLoc;
        foreach my $info (@{$rv->{'genome'}}) {
            my ($qrys, $yn, $tagVal, $locDat, $xv, $str, $ft) = @{$info};
            my $loc = $self->get_location( @{$locDat} );
            my ($type, $lkey) = ('Location');
            if ($loc->is_known()) {
                # This is a known location
                if ($ft) {
                    # This is being forced as a transient location / range
                    $lkey = $loc->full_loc();
                    $loc->lock();
                    $type = 'Range';
                } else {
                    $lkey = $loc->pkey();
                }
            } else {
                $lkey = $loc->full_loc();
                $loc->make_transient();
                $type = 'Range' if ($ft);
            }
            if (my $oloc = $transLoc{$lkey}) {
                $loc = $oloc;
            } else {
                $transLoc{$lkey} = $loc;
            }
            if ($xv) {
                #$tagVal ||= {};
                #$tagVal->{"User Alleles"}{$xv} = 1;
                $loc->extra_alleles( $xv, $str );
                splice(@{$locDat}, 3);
            }
            my $targ = $rv->{$type}{$loc->handle()} ||= [ [], -1 ];
            &_extend_request_object( $targ, $qrys, $yn, $loc, $tagVal);
        }
        delete $rv->{'genome'};
    }
    #$self->maploc->prebranch($rv);
    return $rv;    
}

sub _extend_request_object {
    # This method is used to help build the structure generated by
    # classify_requests()
    my ($targ, $qry, $yn, $obj, $tagVal) = @_;
    push @{$targ->[0]}, @{$qry || []};
    if (!defined $targ->[1] || $targ->[1] < $yn) {
        $targ->[1] = $yn;
    }
    $targ->[2] ||= $obj;
    if ($tagVal) {
        my $kv = $targ->[3] ||= {};
        while (my ($k, $vH) = each %{$tagVal}) {
            map { $kv->{$k}{$_}++ } keys %{$vH};
        }
    }
}

=head2 location_query

 Title   : location_query
 Usage   : my @locIds = $ml->location_query( @param )
 Function: Returns a list of Location PKEYs based on user queries
 Returns : An array of loc_id PKEY values
           In scalar context, a hash reference with additional information
 Args    : Associative array of parameters

This method is designed to be the workhorse for recovering
locations. In the list below, all parameters are optional, but at
least one must be provided in order to build a query.

       -hsps Can be a single HSP object, or an array reference of
             them. Will use hsp_footprint() to convert the list into a
             structure of the overall genomic footprint. Alias
             parameters are -hsp, -aln or -alns

   -segments If you already have a genomic segments hash structure, it
             can be provided here. It will be merged with any segments
             found via -hsps.

             Providing either HSPs or segments will define genomic
             regions that the recovered Locations must reside in.

     -popids Array reference of Population pop_id PKEYs
       -pops List of population names. Aliases -pop and -population
     -poptab Database table name holding Population pop_id PKEYs.

             Any of the above will specify one or more populations
             that must be present in any recovered location. That is,
             the population must have assigned at least one allele to
             a variant.

             Population names can be preceded with an exclamation
             point (!), in which case the variant must NOT have an
             allele assigned by the specified population.

Additional options modify the behavior of the primary query:

      -flank Optional, used in conjunction with -hsps to expand the
             genomic footprint of each HSP.

      -limit If set will only recover the requested number of Locations.

    -dumpsql If true then the SQL query will be displayed before
             running. Primarily of debugging interest.

 -nocollapse If provided then any provided HSPs will NOT be collapsed
             using add_hsps(). I can not recall why this is here -
             probably SQL testing.

=cut

sub location_query {
    my $self = shift;
    $self->bench_start();
    my $dbh  = $self->dbh();
    my $args = $self->parseparams( @_ );
    my (%tables, @where, @binds, @tempTabs, %seen);
    my $rv = {};
    my @hsps;
    foreach my $key ('HSP', 'HSPS', 'ALN', 'ALNS') {
        if (my $req = $args->{$key}) {
            if (ref($req) eq 'ARRAY') {
                push @hsps, @{$req};
            } else {
                push @hsps, $req;
            }
        }
    }
    my @segments;

    if (my $segs = $args->{SEGMENTS}) {
        # Explicitly defined segments
        push @segments, $segs;
    }
    if ($#hsps != -1) {
        my $fp = $self->hsp_footprint( \@hsps, $args->{FLANK} );
        push @segments, $fp;
    }
    my ($cn, $sn, @segWhere, @posWhere) = (0,0);
    foreach my $segHash (@segments) {
        while (my ($bld, $chrH) = each %{$segHash}) {
            my $bldQt = $dbh->quote($bld);
            while (my ($chr, $hsps) = each %{$chrH}) {
                $cn++; $sn += $#{$hsps} + 1;
                # This SQL is built explicitly, not with bind variables
                # Postgres query planner is unhappy if given a range query
                # and it does not know in advance what the actual values are
                my @flanked;
                foreach my $hsp (@{$hsps}) {
                    my ($l, $r, $lf, $rf) = @{$hsp};
                    push @flanked, [$l - ($lf || 0), $r + ($rf || 0)];
                }
                my $sHsp = $args->{NOCOLLAPSE} ?
                    [ sort {$a->[0] <=> $b->[0]} @flanked ] :
                    $self->add_hsps( \@flanked );

                push @{$rv->{detailed}}, sprintf
                    ("%s:%d-%d", $chr, $sHsp->[0][0], $sHsp->[-1][1]);
                
                # Tried to build monolithic SQL statements with all coordinate
                # spans built into OR / AND blocks.
                # This was disapointingly slow:

                #my @rng;
                #foreach my $hsp (@{$sHsp}) {
                #    my ($l, $r) = @{$hsp};
                #    push @rng, sprintf
                #        ("(l.pos_to_right > %d AND l.pos_to_left < %d)",
                #         $l + 1, $r - 1);
                #}
                # 
                #my $wc = sprintf("(l.build = %s AND l.chr = %s AND (%s))",
                #                 $dbh->quote($bld), $dbh->quote($chr),
                #                 join(' OR ', @rng));
                #push @segWhere, $wc;
                
                # Instead, build a bunch of sub-SQL statements that
                # will query a single chunk of location at a time
                my $bcWhere = sprintf("l.build = %s AND l.chr = %s", $bldQt,
                                      $dbh->quote($chr));
                foreach my $hsp (@{$sHsp}) {
                    my ($l, $r) = @{$hsp};
                    push @posWhere, sprintf
                        ("%s AND l.pos_to_right > %d AND l.pos_to_left < %d",
                         $bcWhere, $l + 1, $r - 1);
                }
            }
        }
    }
    if ($#posWhere != -1) {
        push @{$rv->{basic}}, sprintf
            ("Locations within %d segment%s on %d chromosome%s", 
             $sn, $sn == 1 ? '' : 's', $cn, $cn == 1 ? '' : 's');
        # push @where, '('.join(' OR ', @segWhere).')';
    }

    my (@popQry);
    if (my $explicit = $args->{POPTAB}) {
        # Population table being provided directly
        my $tt = "ttUser";
        push @popQry, [ $self->_subquery_clause('col1', $explicit)." ", "" ];
        # push @popQry, [ "= $tt.col1", ", $explicit $tt" ];
        push @{$rv->{detailed}}, "User provided population table";
    } elsif (my $popIds = $args->{POPIDS}) {
        # List of pop_ids being provided directly
        my $num = $#{$popIds} + 1;
        if ($num <= 10) {
            # Small number of population IDs
            push @popQry, [ "IN (".join(',', map { '?' } @{$popIds}).")", "", $popIds ];
        } else {
            my $tempTab = $dbh->list_to_temp_table($popIds, 'integer');
            push @tempTabs, $tempTab;
            push @popQry, [ $self->_subquery_clause('col1', $tempTab)." ", "" ];
            # my $tt = "tt$#tempTabs";
            # push @popQry, [ "= $tt.col1", ", $tempTab $tt" ];
            push @{$rv->{detailed}}, "$num user provided populations";
        }
    } elsif (my $popReq = $args->{POP} || $args->{POPS} || 
             $args->{POPULATION}) {
        # More complex request
        $self->bench_start('Build population filter');
        my @reqs = ref($popReq) ? @{$popReq} : ($popReq);
        my $advPopBlock = $args->{ADVPOP} ? 1 : 0;
        my @popSets;
        foreach my $popBloc (@reqs) {
            foreach my $line (split(/\s*[\n\r]+\s*/, $popBloc || "")) {
                my $targ = $advPopBlock ? 
                    $popSets[$#popSets+1] = {} : # Each block is own SQL
                    $popSets[0] ||= {}; # All populations together
                foreach my $pr (split(/\s*[\,]+\s*/, $line || "")) {
                    next unless ($pr);
                    my $mode = 0;
                    if ($pr =~ /^\!(.+)/) {
                        # Exclude the population
                        $pr   = $1;
                        $mode = 1;
                    }
                    my $pop = $self->get_population_fast($pr);
                    unless ($pop) {
                        my $msg = "Unknown population '$pr'";
                        push @{$rv->{errors}}, $msg unless ($seen{$msg}++);
                        next;
                    }
                    $targ->{$mode}{$pop->pkey()} = 1;
                }
            }
        }
        my $otherQry = $#where == -1 ? 0 : 1;
        my (%hum);
        foreach my $popSet (@popSets) {
            while (my ($mode, $idH) = each %{$popSet}) {
                my @leafIds = keys %{$idH};
                $self->bulk_tag_values( \@leafIds, "textToo" );
                # $self->bulk_population_criteria( \@leafIds, $bulkPopCat );
                my @pops = sort { $a->name() cmp $b->name() } 
                map { $self->get_population($_) } @leafIds;
                next if ($#pops == -1);
                map { $_->read(); } @pops;
                my %pids = map { $_ => 1 } 
                map { $_->deep_child_ids(), $_->pkey() } @pops;
                my @ids  = keys %pids;
                my $sp   = $#ids + 1;
                my $what = $mode ? 'Exclude' : 'Require';
                push @{$hum{$what}}, map { $_->name } @pops;
                map { $rv->{pop_id}{$what}{$_} = 1 } @ids;
                if ($sp <= 10) {
                    # Small number of population IDs
                    push @popQry, [ "IN (".join(',', map { '?' } @ids).")",
                                    "", \@ids, $mode ];
                } else {
                    my $tempTab = $dbh->list_to_temp_table(\@ids, 'integer');
                    push @popQry, [ $self->_subquery_clause
                                    ( 'col1', $tempTab)." ", "", $mode ];
                    push @tempTabs, $tempTab;
                }
            }
        }
        foreach my $mode (sort keys %hum) {
            push @{$rv->{detailed}}, "$mode populations with frequencies in ".
                join(', ', @{$hum{$mode}});
        }
        $self->bench_end('Build population filter');
    }
    my $popOnly = 0;
    if ($#popQry != -1) {
        if ($#where == -1 && $#posWhere == -1) {
            # This is ONLY a population query
            $tables{l} = 'allele';
            $popOnly = 1;
        }
        my $msg = "";
        my $exc = $rv->{pop_id}{Exclude};
        if ($exc) {
            $exc = $rv->{pop_id}{Exclude} = [keys %{$exc}];
            my $num = $#{$exc} + 1;
            $msg .= sprintf("Exclude %d population%s",
                            $num, $num == 1 ? '' : 's');
        }
        if (my $inc = $rv->{pop_id}{Require}) {
            map { delete $inc->{$_} } @{$exc || []};
            $inc = [keys %{$inc}];
            my $num = $#{$inc} + 1;
            if ($num) {
                $rv->{pop_id}{Require} = $inc;
                $msg .= ", " if ($msg);
                $msg .= sprintf("Require %d population%s",
                                $num, $num == 1 ? '' : 's');
            } else {
                # Oops, exclude got rid of them all
                delete $rv->{pop_id}{Require};
                $msg .= ", but then also exclude EVERY one!"
            }
        }
        push @{$rv->{basic}}, $msg;
        my $lims = $dbh->limit_syntax();
        for my $pq (0..$#popQry) {
            my ($popTest, $popTab, $bind, $exclude) = @{$popQry[$pq]};
            my $pl = 'popLim'.($pq+1);
            if ($popOnly) {
                push @where, sprintf("l.pop_id %s", $popTest);
                if ($popTab =~ /(\S+) (\S+)$/) {
                    $tables{$2} = $1;
                }
            } else {
                push @where, sprintf
                    ("%sEXISTS (SELECT %s.pop_id FROM allele %s%s WHERE ".
                     "%s.loc_id = l.loc_id AND %s.pop_id %s %s 1)",
                     $exclude ? "NOT " : "", $pl, $pl, $popTab, 
                     $pl, $pl, $popTest, $lims);
            }
            push @binds, @{$bind} if ($bind);
        }
    }

    if ($#where == -1 && $#posWhere == -1) {
        push @{$rv->{basic}}, "No recognizable query provided";
        $rv->{loc_id} = [];
        return wantarray ? () : $rv;
    }
    push @posWhere, undef if ($#posWhere == -1);
    
    $tables{l} ||= 'location';
    # This used to be DISTINCT l.loc_id
    # The distinct clause seemed to be causing some perfomance issues
    # Handling it with Perl instead.

    my $baseSql = "SELECT l.loc_id";
    $baseSql   .= " FROM ".join(', ', map { "$tables{$_} $_" } keys %tables);
    $baseSql   .= " WHERE ". (join(' AND ', @where) || "");
    my (@lids, %got);
    foreach my $pw (@posWhere) {
        my $sql = $baseSql;
        $sql   .= ($#where == -1 ? '' : ' AND ') . $pw if ($pw);
        my $sth = $dbh->prepare
            ( -name  => 'Custom location query : '.join
              (' + ', @{$rv->{basic}}),
              -limit => $args->{LIMIT},
              -sql   => $sql, );
        if (my $ds = $args->{DUMPSQL}) {
            my $txt = $sth->pretty_print( @binds );
            $txt .= $sth->explain_text( [@binds] ) if ($ds eq '2');
            $txt = $self->esc_xml($txt);
            $self->msg("[SQL]", $txt );
        }
        $self->bench_start("DB Query");
        foreach my $lid ($sth->get_array_for_field( @binds )) {
            push @lids, $lid unless ($got{$lid}++);
        }
        $self->bench_end("DB Query");
    }
    foreach my $tempTab (@tempTabs) {
        $dbh->clear_temp_table( $tempTab );
    }
    $rv->{loc_id} = \@lids;
    $self->bench_end();
    return wantarray ? @lids : $rv;
}

=head1 Technical Search Methods

These methods can recover derived information, but generally are not
used directly. Instead they are utilized by other methods, either
because they are using internal IDs as queries, or ar returning
internal IDs as result sets.

=head2 loc_ids_between_flanks

 Title   : loc_ids_between_flanks
 Usage   : my @loc_ids = $ml->loc_ids_between_flanks( $chr, $lft, $rgt, $build )
 Function: Get a list of loc_ids falling within a chromosomal range
 Returns : An array of Location PKEYs
 Args    : [0] The chromosome name
           [1] The left flank coordinate
           [2] The right flank coordinate
           [3] The genome build name

Queries table location to find locations that match:

 * chr = $chr  AND build = $build
 * pos_to_right >= $lft AND l.pos_to_left <= $rgt

=head3

Used by loc_ids_for_HSPs(), responsible for recovering variants given
a chromosomal range, so it is important that it work efficiently. It
has been problematic in Postgres (8.3.7b). If explicit values are used
for lft/rgt (ie no placeholder) then PG will try using
pk_location_loc_id, which is useless and results in a full table
scan. Using placeholders results (oddly) in index location primary
being used.

=cut

sub loc_ids_between_flanks {
    my $self = shift;
    $self->bench_start();
    my ($chr, $lft, $rgt, $build) = @_;
    $build ||= $self->{BUILD};
    my @bnds = ($chr, $build);
    my $whr;
#    if (1) {
        # This is the most elegant represenation:
        # 42.79ms for NM_001234
        $whr  = "pos_to_right >= ? AND l.pos_to_left <= ?";
        push @bnds, ($lft, $rgt);
#    } elsif (0) {
        # Alterantively, can check to see if points are between bounds
        # 1.35sec for NM_001234

#        $whr = "(l.pos_to_right BETWEEN ? AND ? ".
#            "OR l.pos_to_left BETWEEN ? AND ? ".
#            "OR ? BETWEEN l.pos_to_right AND l.pos_to_left)";
#        push @bnds, ($lft, $rgt, $lft, $rgt, $lft);

        # The first two clauses check to see if either side of the
        # *variant* overlaps the range. The third clause checks to see
        # if the left side of the *range* overlaps the variant. This
        # is to check for the final possibility - a very large variant that
        # completely contains the range.

        # When prepared this query uses location_primary but only on
        # chr and build - after recovering all locations for a
        # chromosome, it then filters on the positions. NOT GOOD.
#    }

    my $qry  = $self->{STH}{MATCH_ACC_TO_LOC} ||= $self->dbh->prepare
        ( -name => "Get all locations between flanks",
          -level => 2,
          -sql  => "SELECT l.loc_id FROM location l ".
          "WHERE l.chr = ?  AND l.build = ? AND $whr");
    # &preprint( $qry->pretty_print( @bnds )."\nEXPLAIN:\n".$qry->explain_text( [@bnds] )."\n");
    my @ids = $qry->get_array_for_field( @bnds );
    $self->bench_end();
    return wantarray ? @ids : \@ids;
}

=head2 locations_for_HSPs

 Title   : locations_for_HSPs
 Usage   : my @locationObjects = $ml->locations_for_HSPs( $hsps, $flank )
 Function: Get locations for one or more HSP boundaries
 Returns : An array of Location objects, or an array reference in scalar context
 Args    : [0] An array reference of BMS::SnpTracker::MapLoc::HSPed objects
           [1] Optional flank distance

Will call loc_ids_for_HSPs() to convert the parameters into Location
loc_id values, then get_location() to convert those to Location
objects.

=cut

sub locations_for_HSPs {
    my $self = shift;
    my @locids = $self->loc_ids_for_HSPs( @_ );
    my @locs = map { $self->get_location( $_ ) } @locids;
    return wantarray ? @locs : \@locs;
}

=head2 loc_ids_for_HSPs

 Title   : loc_ids_for_HSPs
 Usage   : my @locIds = $ml->loc_ids_for_HSPs( $hsps, $flank )
 Function: Finds locations that overlap a set of HSP objects
 Returns : An array of loc_id values, or an array ref if in scalar context
 Args    : [0] an array reference of BMS::SnpTracker::MapLoc::HSPed objects
           [1] Optional flank expansion value

Used by locations_for_HSPs()

=cut

sub loc_ids_for_HSPs {
    my $self = shift;
    my ($hsps, $flank) = @_;
    my $fp   = $self->hsp_footprint( $hsps, $flank );
    my %locids;
    while (my ($bld, $chrH) = each %{$fp}) {
        while (my ($chr, $hsps) = each %{$chrH}) {
            foreach my $hsp (@{$hsps}) {
                map { $locids{$_} = 1 } $self->loc_ids_between_flanks
                    ( $chr, $hsp->[0], $hsp->[1], $bld);
            }
        }
    }
    my @ids = sort { $a <=> $b } keys %locids;
    return wantarray ? @ids : \@ids;
}

=head2 acc_to_locid

 Title   : acc_to_locid
 Usage   : my @locObjs = $ml->acc_to_locid( $string )
           my $locAndAuthList = $ml->acc_to_locid( $string )
 Function: Convert a variant accession to Location objects or a list
           of loc_id and auth_id PKEYs
 Returns : An array of Location objects or (scalar context) a 2D array
           reference of [loc_id, auth_id]
 Args    : [0] A string representing the accession, eg "rs12345"

If the request exists in table accesion, any location that has been
assigned to it will be returned. Called in array context will return
fully blessed Location objects. Scalar context will return a 2D array
reference, with each element being [loc_id, auth_id] taken from the
accession table.

=cut

*accession_to_location = \&acc_to_locid;
sub acc_to_locid {
    my $self = shift;
    $self->bench_start();
    my $acc  = shift || "";
    my $qry  = $self->{STH}{MATCH_ACC_TO_LOC} ||= $self->dbh->prepare
        ( -name => "Get location for accession",
          -sql  => "SELECT l.loc_id, a.auth_id ".
          "FROM location l, accession a, normtxt nt ".
          "WHERE nt.md5upper = md5(?) AND upper(nt.txt) = ? ".
          "AND a.acc_id = nt.txt_id AND l.loc_id = a.obj_id");
    my $uAcc = uc($acc);
    $qry->execute( $uAcc, $uAcc );
    my $rows = $qry->fetchall_arrayref();
    if (wantarray) {
        my %ids = map { $_->[0] => 1 } @{$rows};
        my @rv = map { $self->get_location( $_ ) } keys %ids;
        $self->bench_end();
        return @rv;
    } else {
        $self->bench_end();
        return $rows;
    }
}

=head1 Impact Information

These methods standardize and normalize impact information or
requests, and recover metadata for each impact code.

=head2 impact_details

 Title   : impact_details
 Usage   : my $impData = $ml->impact_details( $impRequest )
 Function: Get detailed information about an impact by name or token
 Returns : A hash reference,or undef for unrecognized requests
 Args    : [0] The name or token related to the impact

Provides full information for the requested impact. The data are
organized in a small hash structure, with the following keys:

  token : Three letter upper-case token for the impact ("NON", "INT")
   rank : relative 'importance' of the impact, smaller values are
          generally more 'impactful'
   name : A short, human-readable name like "Non-synonymous"
  color : A hexidecimal color, applied to the background for text and
          the fill for graphical objects
  fgcol : A hexidecimal color, applied to the foreground for text and
          the border for graphical objects
   desc : A long, human-readable description of the impact
  alias : Optional. Space-separated list of recognized aliases for the impact.
          Example: UT3 = UTR3 UTR-3 3-UTR 3UTR

=cut

sub impact_details {
    my $self = shift;
    my $tok  = uc(shift || "");
    if (my $realTok = $impactAlias->{$tok}) { $tok = $realTok; }
    return exists $impactRanks->{$tok} ? $impactRanks->{$tok} : undef;
}

=head2 impact_list

 Title   : impact_list
 Usage   : my @impTokens = $ml->impact_list( )
 Function: Get a list of all impact tokens 
 Returns : An array of short strings

This function simply provides a list of all the primary impact token
strings (eg "NON", "UT5") recognized by the system. It will not
provide aliases.

=cut

sub impact_list { return @allImpacts; }

=head2 all_impacts

 Title   : all_impacts
 Usage   : my @impData = $ml->all_impacts( )
 Function: Get full details for each impact code
 Returns : An array of hashes

Similar to impact_list(), but will return an array of hash structures
generated by impact_details().

=cut

sub all_impacts {
    my $self = shift;
    return map { $self->impact_details( $_ ) } $self->impact_list();
}

=head2 impact_name

 Title   : impact_name
 Usage   : my $name = $ml->impact_name( $impRequest )
 Function: Get the human-readable name for an impact
 Returns : A short string
 Args    : [0] A name, token or alias for an impact

If the request is unrecognized the empty string will be returned.

=cut

sub impact_name {
    my $self = shift;
    my $id   = $self->impact_details( @_ );
    return $id ? $id->{name} : "";
}

=head2 impact_html_token

 Title   : impact_html_token
 Usage   : my $htmlSpan = $ml->impact_html_token( $impRequest )
 Function: Get a small block of HTML representing a colorized impact token
 Returns : A modest HTML string
 Args    : [0] A name, token or alias for an impact

Will return an HTML string encoding a span element. The element will
have background and foreground colors appropriately set for the
impact, and the title will be the full description.

Null requests will return empty string. Unrecognized requests will
return a span with the request plus "?".

=cut

sub impact_html_token {
    my $self = shift;
    my $tok  = shift || "";
    return "" unless ($tok);
    my $i    = $self->impact_details($tok);
    return "<span color='gray'>&iquest;$tok?</span>" unless ($i);
    my @sty;
    if (my $c = $i->{fgcol}) {
        push @sty, "color: $c"; }
    if (my $c = $i->{color}) {
        push @sty, "background-color: $c"; }
    my $html = "<span";
    my @db;
    foreach my $t ('name','desc') {
        if (my $dv = $i->{$t}) { push @db, $dv; }
    }
    if (my $d = join(': ', @db)) {
        $html .= " title='".$self->esc_xml_attr( $d )."'";
        push @sty, "cursor: help";
    }
    $html .= " style='".join('; ',@sty)."'" 
        unless ($#sty == -1);
    $html .= ">$tok</span>";
    return $html;
}

=head2 heaviest_impact

 Title   : heaviest_impact
 Usage   : my $mostImpactful = $ml->heaviest_impact( @listOfImpacts )
 Function: Select the most disruptive impact from the provided list
 Returns : A short string
 Args    : [0..] Array of impact names, tokens or aliases

Used to find the "biggest" impact in a list. This is typically used to
determine what impact should be used for a variant with more than one
(depending on the transcript it falls in). For example, a variant that
could be INT, UT5 or SYN will return "SYN".

Biologically, it is not always clear which impacts are "heavier". For
example, Deletions (rank 3) are considered "heavier" than
Non-synonymous variants (rank 4). There will be cases where this is
NOT the case (a missing alanine might go unnoticed, but would be very
disruptive if changed to a histidine). There are no universal rules or
rankings to reliably predict such effects, but the default rankings
should generally be reasonable in most circumstances.

=cut

sub heaviest_impact {
    my $self = shift;
    my @ranked = 
        sort { $a->{rank} <=> $b->{rank} }
    map { $impactRanks->{$_ || "" } ||  $impactRanks->{UNK} }
    map { split(/\//, $_ || "") } @_;
    $ranked[0] ||= $impactRanks->{UNK};
    return wantarray ? @ranked : $ranked[0]{token};
}

=head2 expand_impacts

 Title   : expand_impacts
 Usage   : my @fullImp = $ml->expand_impacts
              ( $listOfImpacts, $parentMode, $kidMode )
 Function: Expand a list of impacts into more specific or generic categories
 Returns : An array of impact tokens (strings), or an array reference of same
 Args    : [0] An array reference of one or more impact requests
           [1] Parent mode - default 0
           [2] Child mode - default 0

In addition to 'specific' impacts, such as NON, INT, etc, there are
'generic' parents that represent lower specificity classes, such as
UTR (either 3' or 5') or COD (coding, not certain if synonymous or
non-synonymous). These generic impacts could occur either because
insufficient data are available (for example, a location reported
within a coding frame, but without alleles) or because software is not
fully calculating impact (for example, fast filtering that can
determine that a location is within an exon, but is not attempting to
determine if it is UTR, SYN, etc).

This method will take one or more impacts, and expand it into more
generic impacts, more specific impacts or both. All the originally
provided impacts will be returned as well as any additional parents or
children that get added.

By default, parent nodes will be added to the list only if ALL their
children are present. So "UTR" will be included only if the input list
contains both "UT3" and "UT5". However, if parent mode is true, then a
parent will be added if ANY of its children are present.

By default, ALL child nodes will be added if their parent is in the
input. If child mode is true, then NO children will be included.

=cut

sub expand_impacts {
    my $self = shift;
    my $reqs  = shift || [];
    my $pMode = shift || 0; # Parent behavior mode
    my $kMode = shift || 0; # Kid behavior mode
    my $uniq;
    foreach my $req (@{$reqs}) {
        if (my $id = $self->impact_details($req)) {
            $uniq ||= {};
            $uniq->{$id->{token}} = 1;
        } else {
            $self->msg_once("[!]", "Unrecognized impact request '$req'");
        }
    }
    unless ($uniq) {
        return wantarray ? () : [];
    }
    my %expanded = %{$uniq};
    while (my ($par, $kids) = each %{$genericImp}) {
        if ($uniq->{$par}) {
            # The general parent token already exists
            if ($kMode) {
                # Do not do anything extra about children
            } else {
                # By default, add in all the more specific children
                map { $expanded{$_} = 1 } @{$kids};
            }
        } else {
            # This parent does not exist. Should we add it?
            my $kidsFound = -1;
            map { $kidsFound++ if ($uniq->{$_}) } @{$kids};
            # my @tmp; map { push @tmp, $_ if ($uniq->{$_}) } keys %{$uniq}; warn "$par ($kidsFound) has [".join(' ', @{$kids})."] ($#{$kids}) via [".join(' ', @tmp)."]\n";
            if ($pMode) {
                # Add general parent if ANY of the kids are present
                unless ($kidsFound == 0) {
                    $expanded{$par} = 1;
                }
            } else {
                # Add general parent only if ALL the kids are present
                $expanded{$par} = 1 if ($kidsFound == $#{$kids});
            }
        }
    }
    my @rv = sort keys %expanded;
    return wantarray ? @rv : \@rv;
}

=head1 Settings

Various methods that do not fit in elswhere. Generally boring, but may
be useful or critical in some circumstances.

=head2 instance_name

 Title   : instance_name
 Usage   : my $instance = $ml->instance_name( )
 Function: Gets the database instance name
 Returns : A string

Must be set when new() is called, can not be changed later.

=cut

*db_name = \&instance_name;
sub instance_name { return $instanceName;  }

=head2 build

 Title   : build
 Usage   : my $genomeBuild = $ml->build( $newVal )
 Function: Gets / Sets the default genome build
 Returns : A string
 Args    : [0] Optional new value

A DefaultBuild is not neccesary, but will be used when ambiguous
queries are provided. For example "4:199246" will recover all known
Chromosome "4" entries, but if the build has been set to "GRCh37" then
only "4.GRCh37:199246" will be considered.

=cut

sub build   {
    my $self = shift;
    if (my $nv = shift) {
        $self->{BUILD} = $nv;
    }
    return $self->{BUILD};
}

=head2 verbosity

 Title   : verbosity
 Usage   : my $currentVerbosity = $ml->verbosity( $optionalNewValue )
 Function: Gets / Sets how chatty the program is
 Returns : The current setting
 Args    : [0] Optional new value

Not heavily used. If set to a true value, some additional information
will be reported as the program runs.

=cut

sub verbosity {
    my $self = shift;
    if (defined $_[0]) {
        $self->{VERBOSITY} = $_[0];
    }
    return $self->{VERBOSITY};
}

=head2 default_loc_to_rna_distance

 Title   : default_loc_to_rna_distance
 Usage   : my $genomicDistance = $ml->default_loc_to_rna_distance( $newVal )
 Function: Gets / sets the maximum "associated" distance between
           Locations and RNAs
 Returns : An integer
 Args    : [0] Optional new value

Sets the default distance a variant can be from an RNA and still be
"associated" with it. The default is typically 500 bases. Variants
that are further than this distance upstream or downstream of the RNA
will not be associated with them. If a variant is associated with no
RNAs, it will be assigned a GEN (genomic) impact.

=cut

sub default_loc_to_rna_distance {
    my $self = shift;
    if (defined $_[0]) {
        if ($_[0] =~ /^\d+$/) {
            $self->{LOC2RNA_DIST} = $_[0];
        } else {
            $self->err("Can not set loc-to-RNA distance to '$_[0]'",
                       "Please provide an integer");
        }
    }
    return $self->{LOC2RNA_DIST};
}

=head2 default_howbad

 Title   : default_howbad
 Usage   : my $howbad = $ml->default_howbad( $newVal )
 Function: Gets / Sets the default "HowBad" value for transcripts
 Returns : A value between zero and 100
 Args    : [0] Optional new  value

Transcripts loaded into the system can be associated with more than
one genomic location, if the administrator so wishes. Each transcript
has an overal match percentage associated with each genomic location
it is assigned to. The highest scoring location has a "HowBad" value
of zero - it is the best location. All other locations have a HowBad
set to the best match score minus their own. So if the best match was
98%, and the second best was 95.5%, then the HowBad for the second
entry would be 2.5.

The default value for HowBad is zero - only the best position(s) will
be reported for any transcript. Please note that the genome is
repetitive - highly so in some regions! This means that a transcipt
may not have a unique HowBad=0 location - two or more locations might
be "equally best".

=cut

sub default_howbad {
    my $self = shift;
    if (defined $_[0]) {
        $self->{DEF_HOWBAD} = $_[0];
    }
    return $self->{DEF_HOWBAD};
}

=head1 Utility Objects

A handful of objects that provide supporting functionality.

=head2 connect

 Title   : connect
 Usage   : my $dbh = $ml->connect( )
 Function: Connect to the MapLoc database
 Returns : A DBI database handle

The MapLoc parameter file will be read to recover database
settings. This file is located in the same folder as MapLoc.pm, but
will be named MapLoc.param. If PGPORT and PGHOST are defined there
they will be used to set the current environment.

If instance_name() has been set, it will be used for the DB
instance. Otherwise INSTANCENAME or INSTANCE from the parameter file
will be used. Failing that, 'maploc' will be used.

=cut

*dbh = \&connect;
*dbi = \&connect;
sub connect {
    my $self = shift;
    unless ($self->{DBH}) {
        $self->bench_start();
        # Try to read connection details from parameter file
        my $params = BMS::ArgumentParser->new
            ( -paramfile  => 'BMS::SnpTracker::MapLoc' );
        foreach my $key ('PGPORT', 'PGHOST') {
            if (my $val = $params->val($key)) {
                $ENV{$key} = $val;
            }
        }
        $instanceName ||= $params->val(qw(instancename instance)) || 'maploc';
        my $dbType   = "dbi:Pg:dbname=$instanceName";
        my $dbn      = 'tilfordc';
        my $dbh      = BMS::FriendlyDBI->connect
            ($dbType, $dbn,
             undef, { RaiseError  => 0,
                      PrintError  => 0,
                      LongReadLen => 100000,
                      RowCacheSize => 1,
                      AutoCommit  => 1, },
             -errorfile => '/tmp/MapLocErrors.err',
             -adminmail => 'tilfordc@bms.com', );
        $dbh->schema( $self->schema() );

        $self->{private_STH} = {
        };
        $self->{DBH} = $dbh;
        $self->bench_end();
    }
    return $self->{DBH};
}

=head2 maploc

 Title   : maploc
 Usage   : my $pointlessSelfReference = $ml->maploc( )
 Function: Returns itself
 Returns : The same object used to call the method

This method is included to allow some generic calls from shared
packages to work properly when called using the MapLoc object.

=cut

sub maploc  { return shift; }

=head2 serializer

 Title   : serializer
 Usage   : my $ser = $ml->serializer( )
 Function: Get an object for data serialization
 Returns : A BMS::Utilities::Serialize object

The serializer is used internally, mostly for turning Perl data
structures into JSON text.

=cut

sub serializer {
    my $self = shift;
    unless ($self->{SERIALIZER}) {
        $self->bench_start();
        my $ser = $self->{SERIALIZER} = BMS::Utilities::Serialize->new();
        $ser->set_random_literal_token();
        $self->bench_end();
    }
    return $self->{SERIALIZER};
}

=head2 sequtils

 Title   : sequtils
 Usage   : my $su = $ml->sequtils( )
 Function: Get an object for basic biological sequence manipulation
 Returns : A BMS::Utilities::SequenceUtilities object

Used internally for manipulating nucleotide sequences, primarily for
reverse complementation.

=cut

sub sequtils {
    my $self = shift;
    unless ($self->{SEQUTILS}) {
        $self->bench_start();
        my $ser = $self->{SEQUTILS} = BMS::Utilities::SequenceUtilities->new();
        $self->bench_end();
    }
    return $self->{SEQUTILS};
}

=head1 Utility Methods

Method calls that are useful but secondary to the primary role of MapLoc

=head2 disconnect

 Title   : disconnect
 Usage   : $ml->disconnect( )
 Function: Disconnect from the database

=cut

sub disconnect {
    my $self = shift;
    if (my $dbh = $self->{DBH}) {
        $dbh->disconnect();# || $self->err( "Disconnection failed" );
    }
}

=head2 fork_safe

 Title   : fork_safe
 Usage   : $ml->dbh->fork_safe( )
 Function: Set the DBH into "fork safe" mode

This method simply calls $ml->dbh->fork_safe()

Perl forking does horrible things to database handles. The methods
fork_safe() and fork_unsafe() are used to make sure that a handle is
not shared between two child processes.

=cut

sub fork_safe   { shift->dbh->fork_safe(); }

=head2 fork_unsafe

 Title   : fork_unsafe
 Usage   : $ml->dbh->fork_unsafe( )
 Function: Set the DBH into "fork unsafe" mode

This method simply calls $ml->dbh->fork_unsafe(). See above.

=cut

sub fork_unsafe { shift->dbh->fork_unsafe(); }

=head2 revcom

 Title   : revcom
 Usage   : my $rcDNA = $ml->revcom( $nucSequence )
 Function: Reverse complement a sequence
 Returns : A string
 Args    : [0] A string of DNA or RNA characters

This is just a convienence call of $ml->sequtils->revcom()

=cut

sub revcom {
    return shift->sequtils->revcom( @_ );
}

=head1 Input Normalization

Methods that standardize or sanity check external input

=head2 standardize_chromosome

 Title   : standardize_chromosome
 Usage   : my ($chr, $build, $taxa, $type) = $ml->standardize_chromosome( $req )
 Function: Parse a chromosome request into standardized parts
 Returns : Four strings, or a single string when called in scalar context
 Args    : [0] A string describing the request

MapLoc uses several internal and legacy mechanisms for carefully
defining ChromosomeNomenclature. This method will parse the various
nomenclatures and return the four primary components that completely
define a chromosome or contig:

  * Chromosome name (eg "X", "14")
  * Genome BuildToken (eg "NCBIM36", "GRCh37")
  * Species name (eg "homo_sapiens", "canis_lupus_familiaris")
  * Sequence type (eg "chromosome", "contig")

If called in a scalar context, only the chromosome name will be returned

=cut

sub standardize_chromosome {
    my $self = shift;
    $self->bench_start();
    my ($chr) = @_;
    my ($taxa, $type, $build) = ("", "", "");
    if ($chr =~ /^([^\.]+)\.([a-z]+)\.([XYZW]|MT|\d{1,2})\.([^\.]+)$/i) {
        # homo_sapiens.chromosome.10.GRCh37
        $taxa  = $1;
        $type  = lc($2);
        $chr   = uc($3);
        $build = $4;
    } elsif ($chr =~ /^([^\.]+)\.([a-z]+)\.(NT_\d{6,9})\.([^\.]+)$/i) {
        # homo_sapiens.supercontig.NT_113884.NCBI36
        $taxa  = $1;
        $type  = lc($2);
        $chr   = uc($3);
        $build = $4;
    } elsif ($chr =~ /^([^_]+)_([a-z]+)_([XYZW]|MT|\d{1,2})\.([^\.]+)$/i) {
        # Human_Chr_11.NCBI_31 - Older format
        $taxa  = $1;
        $type  = lc($2);
        $chr   = uc($3);
        $build = $4;
    } elsif ($chr =~ /^([XYZW]|MT|\d{1,2})(\.(\S+))?$/i) {
        $chr = uc($1);
        $build = $3 || "";
    } else {
        $self->death("Unrecognized chromosome '$chr' from '$_[0]'");
    }
    $self->bench_end();
    return wantarray ? ($chr, $build, $taxa, $type) : $chr;
}

=head2 position_to_flanks

 Title   : position_to_flanks
 Usage   : my ($left, $right) = $ml->position_to_flanks( @request )
 Function: Converts Start/Stop nomenclature to Left/Right flanks
 Returns : Two integers (left, right)
 Args    : One or two position definitions

Most API calls utilize LeftRightCoordinates, where a position is
defined by the bases to the left and right. This method is used to
convert from more typical Start/Stop coordinates. A few common single
string formats are recognized:

  L^R  : eg "4^5", an insertion position. Left = 4, Right = 5
  S..E : eg "10..12", a range of three bases. Left = 9, Right = 13
  S-E  : As above, but with a dash instead of ".."
  S    : A single integer, eg "104". Left = 103, Right = 105

If two parameters are passed, and the first does not match any of the
non-single formats described above, then the second parameter is
presumed to be the end coordinate, and the first the start.

=cut

sub position_to_flanks {
    # Standardizes "human-readable" positions to
    # "position at left" + "position at right"
    my $self = shift;
    my ($lft, $rgt) = @_;
    if ($lft =~ /^(\d+)^(\d+)$/) {
        # Tokenized gap coordinates
        ($lft, $rgt) = ($1, $2);
    } elsif ($lft =~ /^(\d+)\.\.(\d+)$/) {
        # Coordinate range
        ($lft, $rgt) = ($1 - 1, $2 + 1);
    } elsif (!defined $rgt) {
        # User is passing a single position
        $rgt = $lft + 1;
        $lft--;
    } else {
        # Both defined
        $lft--;
        $rgt++;
    }
    return ($lft, $rgt);
}

=head2 request_to_args

 Title   : request_to_args
 Usage   : my ($chr, $lft, $rgt, $build) = $ml->request_to_args( $c, $s, $e )
 Function: Standardize chromosome, start and end
 Returns : Four strings
 Args    : [0] Chromosome designation
           [1] Start coordinate
           [2] End coordinate

This method combines standardize_chromosome(), which parses the first
argument into a chromosome name and GenomeBuild, with
position_to_flanks(), which parses the second and optionally third
value into Left/Right coordinates.

=cut

sub request_to_args {
    my $self = shift;
    my ($chr, $build) = $self->standardize_chromosome( shift );
    my ($s, $e) = @_;
    $build = $_[2] if ($_[2]); # Allow user to over-ride the build
    my ($lft, $rgt) = $self->position_to_flanks( $s, $e );
    return ($chr, $lft, $rgt, $build || $self->{BUILD});
}

our $chrRE     = "(chr|chromosome)?\\.?(x|y|\\d{1,2})";
our $buildRE   = "([A-Z]{2,5}\\d{1,2})";
our $rangeRE   = "(\\d+)(\\.\\.|-)(\\d+)";
our $flankRE   = "(\\d+)(_|\\^)(\\d+)";

=head2 

 Title   : parse_genomic_location_text
 Usage   : my ($chr, $s, $e, $bld) = $ml->parse_genomic_location_text( $txt )
 Function: Parse a string to chromosome coordinate information
 Returns : Array of (Chromosome, start, end, build)
 Args    : The text defining the location

Recognizes a variety of formats that define the chromosome, start, end
and build for a genomic location. Not all formats will identify all
four components. Recognition is based on regular expressions, some
formats will not be recognized, for example 2L and 2R (fly) or 2A/2B
(chimp).

=cut

sub parse_genomic_location_text {
    my $self = shift;
    my $txt = shift;
    return () unless ($txt);
    if ($txt =~ /^$chrRE\.$buildRE\:(\d+)(\/(\d+))?$/i) {
        # This is a CHR:POS request with a build specified
        my ($chr, $bld, $s, $xtra) = ($2, $3, $4, $6 || 1);
        return ($chr, $s - $xtra, $s + $xtra, $bld);
    } elsif ($txt =~ /^$chrRE\:(g\.)?(\d+)[ACTGU]?(\/(\d+))?$/i) {
        # This is a CHR:POS request, no build
        my ($chr, $s, $xtra) = ($2, $4, $6 || 1);
        return ($chr, $s - $xtra, $s + $xtra);
    } elsif ($txt =~ /^$chrRE\.$buildRE\:(\d+)\^(\d+)$/i) {
        # This is a CHR:POS request, with build
        return ($2, $4, $5, $3);
    } elsif ($txt =~ /^(g\.)?$chrRE\:(\d+)[ACTGU]?(>[ACTGU])?$/i) {
        # This is a CHR:POS request
        return ($3, $4-1, $4+1);
    } elsif ($txt =~ /^$chrRE\.$buildRE\:$rangeRE[ACTGU]?(\/(\d+))?$/i ||
             $txt =~ /^$chrRE\:(g\.)?$rangeRE[ACTGU]?(\/(\d+))?$/i) {
        my ($chr, $bld, $s, $e, $xtra) = ($2, $3, $4, $6, $8 || 1);
        return ($chr, $s - $xtra, $e + $xtra, $bld);
    } elsif ($txt =~ /^$chrRE\.$buildRE\:$flankRE$/i) {
        # Flank-specified, with build
        return ($2, $4, $6, $3);
    } elsif ($txt =~ /^$chrRE\:(g\.)?$flankRE[ACTGU]?$/i) {
        # Flank-specified, no build
        return ($2, $4, $6);
    }
    return ();
}

=head2 pretty_flanks

 Title   : pretty_flanks
 Usage   : my $niceText = $ml->pretty_flanks( $lft, $rgt )
 Function: Convert Lft/Rgt coordinates into something humans appreciate
 Returns : A string describing a location
 Args    : [0] The left flanking coordinate
           [1] The right flanking coordinate

The database and much of the code logic works with left/right
coordinates, but these are not very human friendly. The method will convert:

    * []      -> "?"      # Undefined
    * [4,6]   -> "5"      # SNP
    * [10,11] -> "10^11"  # InDel
    * [10,19] -> "11..18" # MNP / Range

=cut

sub pretty_flanks {
    my $self = shift;
    my ($lft, $rgt) = @_;
    return "?" unless (defined $lft);
    if ($rgt == $lft + 1) {
        # Gap position
        return "$lft^$rgt";
    } elsif ($rgt == $lft + 2) {
        # Single base
        return $lft + 1;
    }
    # Range
    return join('..', $lft + 1, $rgt -1 );
}

=head2 pretty_coordinates

 Title   : pretty_coordinates
 Usage   : my $niceText = $ml->pretty_flanks( $s, $e )
 Function: Convert Start/End coordinates into something humans appreciate
 Returns : A string describing a location
 Args    : [0] The start coordinate
           [1] The end coordinate

Similar to pretty_flanks() but takes start/end coordinates as
input. Will return a single human-readable string defining the
location.

    * []      -> "?"      # Undefined
    * [5]     -> "5"      # SNP
    * [5,5]   -> "5"      # SNP
    * [11,10] -> "10^11"  # InDel - note that end = start - 1
    * [11,18] -> "11..18" # MNP / Range

=cut

sub pretty_coordinates {
    my $self = shift;
    my ($s, $e) = @_;
    return "?" unless (defined $s);
    return $s unless (defined $e);
    if ($e == $s - 1) {
        # Gap position
        return "$s^$e";
    } elsif ($e == $s) {
        # Single base
        return $s;
    }
    # Range
    return "$s..$e";
}

=head2 comma_number

 Title   : comma_number
 Usage   : my $txt = $ml->comma_number( $number )
 Function: Injects commas every three digits in a number
 Returns : A string
 Args    : An integer

Just adds commas. Kind of silly, but I find it easier to read large
genomic coordinates with the commas present. Appologies to European
friends who use decimal points.

    * 140453136 -> "140,453,136"

=cut

sub comma_number {
    my $self = shift;
    my ($s, $e, $tok) = @_;
    return "" unless ($s);
    my $rv  = $s;
    my $nxt = 3;
    while ($nxt < length($rv)) {
        substr($rv, length($rv) - $nxt, 0) = ',';
        $nxt += 4;
    }
    if ($e) {
        my $sv = $self->comma_number($e);
        my $sl = CORE::length($sv);
        if ($sl == CORE::length($rv)) {
            for (my $p = 0; $p < $sl; $p++) {
                if (substr($sv, $p, 1) eq substr($rv, $p, 1)) {
                    substr($sv, $p, 1) = " ";
                } else {
                    last;
                }
            }
            $sv =~ s/^ +//;
        }
        $rv .= ($tok || '-') . $sv if ($sv);
    }
    return $rv;
}

=head1 Coordinate Manipulation

Functions used to do coordinate math

=head2 hsp_footprint

 Title   : hsp_footprint
 Usage   : my $genomeSegments = $ml->hsp_footprint( $hsps, $flank )
 Function: Coallesce one or more HSP-compliant objects into a genomic structure
 Returns : A hash reference of genome locations
 Args    : [0] An array reference of one or more HSP objects
           [1] Optional flanking distance to expand HSP footprints

The method will use HSPs_to_genomic_footprint() to generate a set of
genomic Range objects. These will then be built into a nested hash
structure:

  Lvl 1: Keyed to the genome build token
     Lvl 2: Keyed to chromosome / contig names
        Lvl 3: Array of [ left, right ] coordinates

The structure will represent the full "shadow" that the list of HSPs
casts across the genome. This is primarily used for building SQL
queries to searches.

=cut

sub hsp_footprint {
    my $self = shift;
    my ($hsps, $flank) = @_;
    my @geno = $self->HSPs_to_genomic_footprint( $hsps, $flank );
    my %fp;
    foreach my $rng (@geno) {
        my $cid  = $rng->anchor_index();
        my $chr  = $rng->chr();
        my $bld  = $rng->build();
        my $hsps = $rng->hsps_for_seq( $rng->anchor_index() );
        push @{$fp{$bld}{$chr}}, $rng->hsps_for_seq( $cid );
    }
    return \%fp;
}

=head2 HSPs_to_genomic_footprint

 Title   : HSPs_to_genomic_footprint
 Usage   : my @ranges = $ml->HSPs_to_genomic_footprint( $hsps, $flank )
 Function: Organize a list of HSPs into discrete genomic segments
 Returns : An array of Range objects, or an array ref if called in scalar context
 Args    : [0] Array reference of one or more HSP-compliant objects
           [1] Optional flank to use for expanding HSP boundaries

This method takes a list of BMS::SnpTracker::MapLoc::HSPed compliant
objects and generates discrete, non-overlaping genomic segments from
them. Each segment will be encoded by a BMS::SnpTracker::MapLoc::Range
object.

The second parameter is an optional flank value. It will be used to
expand both the left and right side of the HSP (left decreased by
flank, right increased by flank).

The method is primarily used to reduce a set of HSPs from one or more
transcripts into the set of genomic ranges the minimally represents
them.

=cut

sub HSPs_to_genomic_footprint {
    my $self = shift;
    my $hsps = shift;
    my $flank = shift || 0;
    my %chrs;
    foreach my $hsp (@{$hsps || []}) {
        my @cids = $hsp->chr_acc_ids();
        unless ($#cids == 0) {
            $self->msg("[!]","Unable to get genomic footprint: ".
                       scalar(@cids)." recognized genomic sequences",
                       $hsp->to_one_line());
            next;
        }
        my $sc  = $hsp->score() || 0;
        my $cid = $cids[0];
        my @LR = map { [ @{$_} ] } $hsp->hsps_for_seq( $cid );
        map { $_->[2] = $sc } @LR;
        push @{$chrs{$cid}}, @LR;
    }
    my @rv;
    while (my ($cid, $hsps) = each %chrs) {
        if ($flank) {
            map { $_->[0] -= $flank;
                  $_->[1] += $flank; } @{$hsps};
        }
        my @sorted = sort { $a->[0] <=> $b->[0] || 
                                $a->[1] <=> $b->[1] } @{$hsps};
        my @sets = ([ shift @sorted ]);
        while (my $hsp = shift @sorted) {
            my $prior = $sets[-1][-1];
            my $pRgt  = $prior->[1];
            my $hLft  = $hsp->[0];
            if ($hLft > $pRgt + 1) {
                # Non overlapping
                push @{$sets[-1]}, $hsp;
            } else {
                # They overlap
                if ($hsp->[2] < $prior->[2]) {
                    # Keep the worst score
                    $prior->[2] = $hsp->[2];
                }
                if ($hsp->[1] > $pRgt) {
                    # Overlapping, and need to extend the right coordinate
                    $prior->[1] = $hsp->[1];
                }
            }
        }
        foreach my $set (@sets) {
            my $hsp  = BMS::SnpTracker::MapLoc::Range->new( $self );
            my $sc = 100;
            my @LR;
            foreach my $lr (@{$set}) {
                push @LR, [ $lr->[0], $lr->[1] ];
                # Keep the worst score:
                $sc = $lr->[2] if ($lr->[2] < $sc);
            }
            $hsp->score($sc);
            $hsp->{STRAND} = [1];
            $hsp->hsps( \@LR );
            $hsp->seq_ids( $cid );
            $hsp->set_anchor( $cid );
            push @rv, $hsp;
        }
    }
    return wantarray ? @rv : \@rv;
}

=head1 Nomenclature Methods

These methods are used to normalize various names or controlled vocabularies

=head2 acc_to_build_and_chr

 Title   : acc_to_build_and_chr
 Usage   : my ($build, $chr) = $ml->acc_to_build_and_chr( $accession )
 Function: Convert a genomic accession into chromosome name and build
 Returns : Array of chromosome name and build, both strings
 Args    : [0] A string representing the genomic accession

Checks the request against table loc_to_acc to see if it is assigned
as an accession for a particular chromosome / build pair.

=cut

sub acc_to_build_and_chr {
    my $self = shift;
    my $acc  = shift || "";
    unless ($self->{ACC2LOC_CACHED}{$acc}) {
        $self->bench_start();
        # do not let cache get too big
        my @check = keys %{$self->{ACC2LOC_CACHED} || {}};
        $self->{ACC2LOC_CACHED} = {} if ($#check > 1000);
        my $rows = [];
        if ($acc) {
            if (my $accid = $self->text_to_pkey( $acc, 'noCreate' )) {
                my $get  = $self->{STH}{ACC2BUILDCHR} ||= $self->dbh->prepare
                    ( -name => "Get build and chromosome given an accession",
                      -sql  => "SELECT build, chr FROM loc_to_acc WHERE acc_id = ?", );
                $get->execute( $accid );
                $rows = $get->fetchall_arrayref();
            }
        }
        $self->{ACC2LOC_CACHED}{$acc} = $rows->[0] || ["",""];
        $self->bench_end();
    }
    return @{$self->{ACC2LOC_CACHED}{$acc}}
}

=head2 

 Title   : all_builds
 Usage   : my @buildTokens = $ml->all_builds()
 Function: Get a list of known genome build tokens from the database
 Returns : In array context, a list of strings
           In scalar context, a hash reference

A scalar call will return a hash reference keyed to build token, with
values indicating the number of accessions associated with the build.


This method is a full database search against the loc_to_acc
table. For mature genome projects, there will be relatively few
accessions for any given build, and this query will complete quickly,
even though it is by nature a seq_scan.  However, for unfinished
genomes there are typically tens or hundreds of thousands of contig
IDs, which will greatly increase the number of rows in a table.

PLANNED: A more surgical method that recovers values from a
pre-calculated (or pre-stored) list.

=cut

sub all_builds {
    my $self = shift;
    my $rv;
    unless ($rv = $self->{ALL_BUILDS}) {
        $self->bench_start();
        my $get = $self->{STH}{GET_ALL_BUILDS} ||= $self->dbh->prepare
            ( -name => "Get all builds mentioned in the loc_to_acc table",
              -sql  => "SELECT build, count(build) FROM loc_to_acc GROUP BY build",);              
        $get->execute( );
        my $rows = $get->fetchall_arrayref();
        $rv = $self->{ALL_BUILDS} = {};
        map { $rv->{$_->[0]} = $_->[1] } @{$rows};
        $self->bench_end();
    }
    return wantarray ? sort keys %{$rv} : $rv;
}

=head2 

 Title   : accession_to_object
 Usage   : my @oids = $ml->accession_to_object( $someTex )
 Function: Get a list of database IDs for a text query
 Returns : A list of obj_id values in array context
           In scalar context, a hash reference of obj_ids pointing
           to authorities
 Args    : [0] A string

Will find object IDs in the accession table matching the requested
string. Wildcards (%) can be used. Array context returns a list of the
IDs, while scalar context returns a hash reference keyed to the IDs,
with values being the authorities (strings) associated with each ID.

=cut

*acc_to_obj = \&accession_to_object;
sub accession_to_object {
    my $self = shift;
    my $rv = {};
    my $req = shift;
    if ($req =~ /\%/) {
        $self->bench_start('wildcard');
        my $get = $self->dbh->prepare
            ( -name => "Get accession for wildcard name",
              -limit => shift,
              -sql  => "SELECT a.obj_id, a.auth_id FROM accession a, normtxt nt WHERE a.acc_id = nt.txt_id AND (upper(substr(nt.txt, 1, 20))) LIKE ?",);
        $get->execute( $req );
        my $rows = $get->fetchall_arrayref();
        map { push @{$rv->{$_->[0]}}, 
              $self->cached_pkey_to_text( $_->[1] ) } @{$rows};
        $self->bench_end('wildcard');
    } elsif (my $oid = $self->text_to_pkey( $req, 'noCreate' )) {
        $self->bench_start();
        my $get = $self->{STH}{ACCID_TO_OBJECT} ||= $self->dbh->prepare
            ( -name => "Get accession for object name",
              -sql  => "SELECT obj_id, auth_id FROM accession WHERE acc_id = ?",);
        $get->execute( $oid );
        my $rows = $get->fetchall_arrayref();
        map { push @{$rv->{$_->[0]}}, 
              $self->cached_pkey_to_text( $_->[1] ) } @{$rows};
        $self->bench_end();
    }
    return wantarray ? keys %{$rv} : $rv;
}

=head2 

 Title   : chrbld_to_accid
 Usage   : my $accid = $ml->chrbld_to_accid( $chrName, $bldToken )
 Function: Get the database acc_id for a specific chromosome and build
 Returns : An acc_id (integer)
 Args    : [0] Chromosome or conting name ("4", "NT_123456.3")
           [1] Build ("GRCh37")

Searches table loc_to_acc to find accession IDs that match the
requested chromosome name and build ID.

=cut

*loc_to_accid = \&chrbld_to_accid;
sub chrbld_to_accid {
    my $self = shift;
    my $build = shift || "";
    my $chr   = shift || "";
    unless (defined $self->{LOC2ACCID_CACHED}{$build}{$chr}) {
        $self->bench_start();
        # do not let cache get too big
        my @check = keys %{$self->{LOC2ACCID_CACHED}{$build} || {}};
        $self->{LOC2ACCID_CACHED}{$build} = {} if ($#check > 1000);
        my $accid;
        if ($build && $chr) {
            my $get  = $self->{STH}{BUILDCHR2ACC} ||= $self->dbh->prepare
                ( -name => "Get build and chromosome given an accession",
                  -sql  => "SELECT acc_id FROM loc_to_acc WHERE build = ? AND chr = ?", );
            $accid = $get->get_single_value( $build, $chr);
        }
        $self->{LOC2ACCID_CACHED}{$build}{$chr} = $accid || 0;
        $self->bench_end();
    }
    return $self->{LOC2ACCID_CACHED}{$build}{$chr};
}

=head2 

 Title   : chrbld_to_acc
 Usage   : my $accession = $ml->chrbld_to_acc( $chrName, $bldToken )
 Function: Get an accession for a specific chromosome and build
 Returns : An accession (string)
 Args    : [0] Chromosome or conting name ("4", "NT_123456.3")
           [1] Build ("GRCh37")

This method simply calls chrbld_to_accid() and then converts the
acc_id into the accession (eg "homo_sapiens.chromosome.7.GRCh37")

=cut

*loc_to_acc = \&chrbld_to_acc;
sub chrbld_to_acc {
    my $self = shift;
    my $build = shift || "";
    my $chr   = shift || "";
    unless (defined $self->{LOC2ACC_CACHED}{$build}{$chr}) {
        $self->bench_start();
        # do not let cache get too big
        my @check = keys %{$self->{LOC2ACC_CACHED}{$build} || {}};
        $self->{LOC2ACC_CACHED}{$build} = {} if ($#check > 1000);
        my $accid = $self->loc_to_accid( $build, $chr );
        $self->{LOC2ACC_CACHED}{$build}{$chr} = $accid ?
            $self->cached_pkey_to_text( $accid ) : "";
        $self->bench_end();
    }
    return $self->{LOC2ACC_CACHED}{$build}{$chr};
}

=head1 Database Methods

These methods modify the database by changing existing rows or
inserting new ones.

=head2 

 Title   : define_acc_for_loc
 Usage   : $ml->define_acc_for_loc( $acc, $build, $chr )
 Function: Associates an accession with a chromosome and build within
           the database
 Returns : The acc_id (integer) if successful
 Args    : [0] The accession (string)
           [1] The genome build token
           [2] The chromosome name or contig

Used to assign an accession to a chromosome (or contig) entry. The
database will only associate a single accession with a genomic
element, and the method will die if an attempt is made to add a second
accession.

=cut

sub define_acc_for_loc {
    my $self = shift;
    $self->bench_start();
    my ($acc, $build, $chr) = @_;
    my $set = $self->{STH}{SET_ACC2BUILDCHR} ||= $self->dbh->prepare
        ( -name => "Set build and chromosome for an accession",
          -sql => "INSERT INTO loc_to_acc (build, chr, acc_id) VALUES (?,?,?)",
          -ignore => "duplicate key value" );

    my $accid = $self->text_to_pkey( $acc );
    $set->execute( $build, $chr, $accid );
    # Check to see that it was set
    my $ca = $self->loc_to_acc( $build, $chr);
    if ($ca eq $acc) {
        $self->bench_end();
        return $accid;
    }
    my @errs;
    push @errs, "$build + $acc is already defined as '$ca'" if ($ca);
    my ($cb, $cc) = $self->acc_to_build_and_chr( $acc );
    push @errs, "$acc is already defined as $cb + $cc"
        if ($cb && $cc && ($cb ne $build || $cc ne $chr));
    $self->death("You can not define $acc = $build + $chr", @errs);
}

=head2 

 Title   : schema
 Usage   : my $schemaHash = $ml->schema()
 Function: Get a data structure representing the database schema
 Returns : A rich Perl hash reference

This method does not query the database. Rather, it contains a large,
static hash structure that defines the tables, indices, views and
functions that are used by the DB. This structure is used when
creating a new database, and is parsed by BMS::FriendlyDBI .

=cut

sub schema {
    my $self = shift;
    $self->bench_start();
    my %tables;

    $tables{"location"} =
    { name  => 'location',
      com   => 'Normalization of coordinate location.',
      index => {
          location_primary   => {
              cols => [ 'chr', 'pos_to_left','pos_to_right', 'build' ],
              unique => 1,
          },
          location_coord   => {
              cols => [ 'pos_to_left','pos_to_right', 'chr', 'build' ],
          },
      },
      pkey  => 'loc_id',
      cols  => 
          [['loc_id', 'integer',
            'Internal pkey representing this location', {
                SEQUENCE => 'global_seq'
            } ],

           ['build', 'varchar(10)',
            'The genome build the location is anchored on' ],

           ['chr', 'varchar(20)',
            'The chromosome or contig holding the location' ],

           ['pos_to_left', 'integer',
            'The position _next to_ the low coordinate of the location' ],

           ['pos_to_right', 'integer',
            'The position _next to_ the high coordinate of the location' ],

           ['len', 'integer',
            'The base width occupied by this location. A true SNP will be 1' ],

           ] };

    $tables{"loc_cat"} =
    { name  => 'loc_cat',
      com   => 'Location categories',
      index => {
          loc_cat_primary   => {
              cols => [ 'loc_id', 'cat_id' ],
              unique => 1,
          },
          loc_cat_by_cat   => {
              cols => [ 'cat_id' ],
          },
      },
      fkey  => {
          loc_id  => 'location.loc_id',
          cat_id  => 'normtxt.txt_id',
      },
      cols  => 
          [['loc_id', 'integer',
            'Fkey pointing to location, indicating the location being described', ],

           ['cat_id', 'integer',
            'Fkey pointing to normtxt, indicating a category name that is being assigned to the location' ],

           ] };

    $tables{"population"} =
    { name  => 'population',
      com   => 'Can represent true populations, cell lines, or individuals',
      index => {
          population_primary   => {
              cols => [ 'name_id', 'src_id' ],
              unique => 1,
          },
          poptypes => {
              cols => [ 'type' ],
          }
      },
      pkey  => 'pop_id',
      fkey  => {
          name_id => 'normtxt.txt_id',
          src_id  => 'normtxt.txt_id',
          par_id  => 'population.pop_id',
      },
      cols  => 
          [['pop_id', 'integer',
            'Primary key for the population', {
                SEQUENCE => 'global_seq'
            } ],

           ['name_id', 'integer',
            'ID (pointing to normtxt) for the name of the population' ],

           ['src_id', 'integer',
            'ID (pointing to normtxt) for the source of the population' ],

           ['type', 'varchar(20)',
            'The kind of population this is' ],

           ['chr', 'integer',
            'The number of chromosomes represented in the population' ],

           ['par_id', 'integer',
            'Optional parent population' ],

           ] };

    $tables{ "v_pop" } =
    { name  => 'v_pop',
      com   => 'Denormalized view of population',
      view  =>
"
SELECT p.pop_id, nt1.txt AS name, nt2.txt AS source, type, chr
FROM normtxt nt1, normtxt nt2, population p
WHERE nt1.txt_id = p.name_id
AND nt2.txt_id = p.src_id
"
};

    $tables{ "v_stat" } =
    { name  => 'v_stat',
      idea  => 'http://stackoverflow.com/questions/6903938/how-do-i-know-if-the-statistics-of-a-postgres-table-are-up-to-date',
      com   => 'Table access statistics',
      view  =>
"
SELECT relname,  seq_scan, idx_scan, 
       n_tup_ins AS Inserts, n_tup_upd AS Updates, n_tup_del AS Deletes,
       to_char(last_analyze, 'YYYY Mon DD') AS Analyzed,
       to_char(last_vacuum, 'YYYY Mon DD') AS Vacuumed
FROM pg_stat_all_tables where schemaname = 'public'
ORDER BY relname
"
};

    $tables{"allele"} =
    { name  => 'allele',
      com   => 'Iteration of observed alleles at specific locations',
      index => {
          allele_pkey   => {
              cols => [ 'loc_id', 'base_id', 'pop_id' ],
              unique => 1,
          },
          # This did not help updates at all, since allele_pkey was always
          # chosen by the query planner:
          # allele_full   => {
          #     cols => [ 'pop_id', 'freq', 'loc_id', 'base_id', 'count'],
          # },
      },
      fkey  => {
          loc_id  => 'location.loc_id',
          base_id => 'normtxt.txt_id',
          pop_id  => 'population.pop_id',
      },
      cols  => 
          [['loc_id', 'integer',
            'Foreign key pointing to table location'],

           ['base_id', 'integer',
            'ID (pointing to normtxt) for the base composition of the allele' ],
           ['pop_id', 'integer',
            'FKEY pointing to population for the source of frequency data' ],

           ['freq', 'real',
            'Frequency, between 0 and 1, or undefined' ],

           ['count', 'integer',
            'For discrete frequency counts, will indicate the total number of entities counted. The number of entities for this allele will then be count x freq' ],

           ] };

    $tables{ "v_allele" } =
    { name  => 'v_allele',
      com   => 'Denormalized view of allele',
      view  =>
"
SELECT l.loc_id, l.build, l.chr, CASE 
  WHEN l.pos_to_left + 2 = l.pos_to_right
  THEN (l.pos_to_left + 1)::text
  WHEN l.pos_to_left + 1 = l.pos_to_right
  THEN l.pos_to_left::text || '^' || l.pos_to_right::text
  ELSE (l.pos_to_left + 1)::text || '..' || (l.pos_to_right - 1)::text
   END AS Location,
  substr(nt1.txt, 1, 30) AS BASE30,
  substr(nt2.txt,1,50) AS POP50, a.freq, a.count
FROM normtxt nt1, location l, allele a
LEFT OUTER JOIN population pop ON (pop.pop_id = a.pop_id)
LEFT OUTER JOIN normtxt nt2 ON (nt2.txt_id = pop.name_id)
WHERE nt1.txt_id = a.base_id
AND l.loc_id = a.loc_id
"
};

    $tables{"min_major_allele"} =
    { name  => 'min_major_allele',
      com   => 'Reports the smallest frequency major allele for a location. Null values indicate no frequency information. The maximum is 1 and indicates no apparent variation. The minimum for a two allele SNP is 0.5',
      index => {
          mm_primary   => {
              cols => [ 'loc_id' ],
              unique => 1,
          },
      },
      fkey  => {
          loc_id  => 'location.loc_id',
          pop_id  => 'normtxt.txt_id',
      },
      cols  => 
          [['loc_id', 'integer',
            'Foreign key pointing to table location'],

           ['pop_id', 'integer',
            'ID (pointing to normtxt) for the source of the assignment' ],

           ['freq', 'real',
            'Frequency, greater than 0 and less than or equal 1, or undefined' ],

           ] };


    $tables{ "v_minmaj" } =
    { name  => 'v_minmaj',
      com   => 'Denormalized view of min_major_allele',
      view  =>
"
SELECT DISTINCT mm.loc_id, l.build, l.chr, CASE 
  WHEN l.pos_to_left + 2 = l.pos_to_right
  THEN (l.pos_to_left + 1)::text
  WHEN l.pos_to_left + 1 = l.pos_to_right
  THEN l.pos_to_left::text || '^' || l.pos_to_right::text
  ELSE (l.pos_to_left + 1)::text || '..' || (l.pos_to_right - 1)::text
   END AS Location,
  substr(nt1.txt, 1, 30) AS ACCESSION,
  mm.freq, substr(nt2.txt, 1, 30) AS POPULATION
FROM normtxt nt1, min_major_allele mm
JOIN location  l ON ( l.loc_id   = mm.loc_id )
JOIN accession a ON (a.obj_id    = mm.loc_id )
LEFT OUTER JOIN population pop ON (pop.pop_id = mm.pop_id)
LEFT OUTER JOIN normtxt nt2 ON (nt2.txt_id = pop.name_id)
WHERE nt1.txt_id = a.acc_id
  AND a.obj_id   = mm.loc_id
"
};

    $tables{"loc_to_acc"} =
    { name  => 'loc_to_acc',
      com   => 'Provides a mapping of build + chr <-> accession. A one-to-one mapping is enforced',
      index => {
          loc2acc_buildchr   => {
              cols => [ 'chr', 'build' ],
              unique => 1,
          },
          loc2acc_acc => {
              cols => ['acc_id'],
              unique => 1,
          },
      },
      fkey  => {
          acc_id  => 'normtxt.txt_id',
      },
      cols  => 
          [['build', 'varchar(10)',
            'The genome build the location is anchored on' ],

           ['chr', 'varchar(20)',
            'The chromosome or contig holding the location' ],

           ['acc_id', 'integer',
            'PKEY for the accession text' ],

           ] };

    $tables{ "v_locacc" } =
    { name  => 'v_locacc',
      com   => 'View of location to accession',
      view  =>
"
SELECT la.build, la.chr, nt.txt as accession, la.acc_id
FROM loc_to_acc la, normtxt nt
WHERE nt.txt_id = la.acc_id
"
};

    $tables{"alignment"} =
    { name  => 'alignment',
      com   => 'Set of one or more HSPs between two sequences',
      pkey  => 'aln_id',
      index => {
          align_byseq1   => {
              cols => [ 'seq1', 'ptl1', 'ptr1' ],
          },
          align_byseq2   => {
              cols => [ 'seq2', 'ptl2', 'ptr2' ],
          },
          align_bysrc1   => {
              cols => [ 'src_id', 'seq1' ],
          },
          align_bysrc2   => {
              cols => [ 'src_id', 'seq2' ],
          },
      },
      fkey  => {
          seq1   => 'normtxt.txt_id',
          seq2   => 'normtxt.txt_id',
          src_id => 'normtxt.txt_id',
      },
      cols  => 
          [['aln_id', 'integer',
            'Internal pkey representing the alignment', {
                SEQUENCE => 'global_seq'
            } ],

           ['seq1', 'integer',
            'ID (pointing to normtxt) for the first sequences alignment' ],

           ['seq2', 'integer',
            'ID (pointing to normtxt) for the second sequences alignment' ],

           ['src_id', 'integer',
            'ID (pointing to normtxt) for the source of the alignment' ],

           ['strand', 'smallint',
            'The strand of the feature, should be 1, -1 or zero' ],

           ['score', 'real',
            'A score for the alignment' ],

           ['howbad', 'real',
            'How bad this alignment is compared to the \'best\' alignment. A value of zero means this is the top scoring alignment. For what depends on how it was loaded, generally it means the world of the non-genomic sequence.' ],

           ['ptl1', 'integer',
            'The position _to the left of_ sequence 1' ],

           ['ptr1', 'integer',
            'The position _to the right of_ sequence 1' ],

           ['ptl2', 'integer',
            'The position _to the left of_ sequence 2' ],

           ['ptr2', 'integer',
            'The position _to the right of_ sequence 2' ],

           ] };

    $tables{"hsp"} =
    { name  => 'hsp',
      com   => 'Ungapped aligned region between two sequences',
      index => {
          hsp_primary   => {
              cols => [ 'aln_id' ],
          },
      },
      fkey  => {
          aln_id => 'alignment.aln_id',
      },
      cols  => 
          [['aln_id', 'integer',
            'FKEY to the alignment table'],

           ['ptl1', 'integer',
            'The position _to the left of_ sequence 1' ],

           ['ptr1', 'integer',
            'The position _to the right of_ sequence 1' ],

           ['ptl2', 'integer',
            'The position _to the left of_ sequence 2' ],

           # NOTE THAT WE DO NOT NEED ptr2. Because there are no gaps
           # it can be infered from the other three values as:
           # ptl2 + ptr1 - ptl1

           ] };

    $tables{"range"} =
    { name  => 'range',
      com   => 'Like hsp, but a set of coordinates for only a single sequence',
      index => {
          range_primary   => {
              cols => [ 'rng_id' ],
          },
      },
      cols  => 
          [['rng_id', 'integer',
            'The primary key grouping together segments in the range'],

           ['ptl', 'integer',
            'The position _to the left of_ the range segment' ],

           ['ptr', 'integer',
            'The position _to the right of_ the range segment' ],

           ] };

    $tables{"normtxt"} =
    { name  => 'normtxt',
      com   => 'Normalization of arbitrary text.',
      index => {
          text_primary   => {
              cols => [ 'md5sum', 'mduniq' ],
              unique => 1,
          },
          text_upper   => {
              cols => [ 'md5upper' ],
          },
          text_partial => {
              cols => [ 'upper(substr(txt,1,20))' ],
          }
      },
      pkey  => 'txt_id',
      cols  => 
          [['txt_id', 'integer',
            'Internal pkey representing the text', {
                SEQUENCE => 'global_seq'
            } ],
           
           ['md5sum', 'varchar(32)',
            'MD5 checksum of txt' ],

           ['mduniq', 'smallint',
            'Integer iterator to allow for checksum collisions.', {
                DEFAULT => 1,
            } ],

           ['md5upper', 'varchar(32)',
            'MD5 checksum of upper-cased txt' ],

           ['txt', 'text',
            'The text being represented' ],

           ] };

    $tables{ "v_text" } =
    { name  => 'v_text',
      com   => 'Compact view of normtxt',
      view  =>
"
SELECT txt_id, length(txt) as Len, CASE WHEN md5sum = md5upper THEN 'Yes' ELSE '' END AS Upper, substring(txt,1,50) AS text_50chars
FROM normtxt
"
};

    $tables{"tagval"} =
    { name  => 'tagval',
      com   => 'Tag-value pairs assigned to arbitrary object',
      index => {
          tagval_primary   => {
              cols => [ qw(obj_id) ],
          },
          tagval_bytag   => {
              cols => [ qw(tag_id obj_id val_id) ],
              unique => 1,
          },
          tagval_byval   => {
              cols => [ qw(val_id obj_id) ],
          },
          tagval_byvaltag   => {
              cols => [ qw(val_id tag_id) ],
          },
      },
      fkey  => {
          tag_id  => 'normtxt.txt_id',
          val_id  => 'normtxt.txt_id',
      },
      cols  => 
          [['obj_id', 'integer',
            'The PKEY for the object being tagged', ],

           ['tag_id', 'integer',
            'ID (pointing to normtxt) for the tag name' ],

           ['val_id', 'integer',
            'ID (pointing to normtxt) for the tag value' ],

           ] };


    $tables{ "v_tagval" } =
    { name  => 'v_tagval',
      com   => 'Map out IDs to text for Accession table',
      view  =>
"
SELECT 'Text' AS Type, tv.obj_id, nt1.txt AS Object, nt2.txt AS Tag, substring(nt3.txt,1,50) AS Value
FROM normtxt nt1, normtxt nt2, normtxt nt3, tagval tv
WHERE nt1.txt_id = tv.obj_id
AND nt2.txt_id = tv.tag_id
AND nt3.txt_id = tv.val_id
UNION
SELECT 'Population' AS Type, tv.obj_id, nt1.txt AS Object, nt2.txt AS Tag, substring(nt3.txt,1,50) AS Value
FROM population p, normtxt nt1, normtxt nt2, normtxt nt3, tagval tv
WHERE p.pop_id = tv.obj_id
AND nt1.txt_id = p.name_id
AND nt2.txt_id = tv.tag_id
AND nt3.txt_id = tv.val_id
UNION
SELECT 'Location' AS Type, tv.obj_id, l.chr || '.' || l.build || ':' ||
        l.pos_to_left || '^' || l.pos_to_right AS Object, 
        nt2.txt AS Tag, substring(nt3.txt,1,50) AS Value
FROM location l, normtxt nt2, normtxt nt3, tagval tv
WHERE l.loc_id = tv.obj_id
AND nt2.txt_id = tv.tag_id
AND nt3.txt_id = tv.val_id
"
};

    $tables{"accession"} =
    { name  => 'accession',
      com   => 'Accessions assigned to a location',
      index => {
          accession_primary   => {
              cols => [ qw(obj_id acc_id auth_id) ],
              unique => 1,
          },
          acc_byacc  => {
              cols => [ 'acc_id'],
          },
      },
      fkey  => {
          acc_id  => 'normtxt.txt_id',
          auth_id => 'normtxt.txt_id',
      },
      cols  => 
          [['obj_id', 'integer',
            'The PKEY for the object holding the accession', ],

           ['acc_id', 'integer',
            'PKEY for the accession text' ],

           ['auth_id', 'integer',
            'PKEY for the authority text' ],

           ] };


    $tables{ "v_acc" } =
    { name  => 'v_acc',
      com   => 'Map out IDs to text for Accession table',
      view  =>
"
SELECT a.obj_id, nt1.txt AS Accession, nt2.txt AS Authority
FROM accession a, normtxt nt1, normtxt nt2
WHERE nt1.txt_id = a.acc_id
AND nt2.txt_id = a.auth_id
"
};

    $tables{"rna"} =
    { name  => 'rna',
      com   => 'RNA sequence and CDS data',
      index => {
          rna_primary   => {
              cols => [ qw(acc_id src_id) ],
              unique => 1,
          },
      },
      pkey  => 'rna_id',
      fkey  => {
          acc_id  => 'normtxt.txt_id',
          src_id  => 'normtxt.txt_id',
          seq_id  => 'normtxt.txt_id',
      },
      cols  => 
          [['rna_id', 'integer',
            'Internal pkey representing this RNA', {
                SEQUENCE => 'global_seq'
            } ],
           
           ['acc_id', 'integer',
            'ID (pointing to normtxt) for the RNA accession' ],
           
           ['src_id', 'integer',
            'ID (pointing to normtxt) indicating the source of the RNA' ],

           ['cds_start', 'integer',
            'The RNA coordinate marking the start of the coding sequence' ],

           ['cds_end', 'integer',
            'The RNA coordinate marking the end of the coding sequence' ],

           ['seq_id', 'integer',
            'ID (pointing to normtxt) defining the actual DNA sequence' ],

           ['anomaly', 'text',
            'CDS feature coordinates if they are not simple start..end' ],

           ] };

    $tables{ "v_rna" } =
    { name  => 'v_rna',
      com   => 'Map out IDs to tex for RNA table',
      view  =>
"
SELECT r.rna_id, nt1.txt AS Accession, nt2.txt AS Source, r.cds_start, r.cds_end, length(nt3.txt) as Length, r.anomaly
FROM rna r, normtxt nt1, normtxt nt2, normtxt nt3
WHERE nt1.txt_id = r.acc_id
AND nt2.txt_id = r.src_id
AND nt3.txt_id = r.seq_id
"
};

    $tables{"feature"} =
    { name  => 'feature',
      com   => 'Positional annotation on an object',
      index => {
          feat_primary   => {
              cols => [ qw(host_id pos_to_left pos_to_right feat_id) ],
              unique => 1,
          },
          feat_name_ind => {
              cols => [ qw(feat_id src_id) ],
          }
      },
      fkey  => {
          feat_id => 'normtxt.txt_id',
          src_id  => 'normtxt.txt_id',
          pri_id  => 'normtxt.txt_id',
          # rng_id  => 'range.rng_id',
      },
      cols  => 
          [['host_id', 'integer',
            'ID for the object hosting the feature' ],

           ['feat_id', 'integer',
            'ID (pointing to normtxt) for the name / label of the feature' ],

           ['src_id', 'integer',
            'ID (pointing to normtxt) for the source of the feature assignment' ],

           ['pri_id', 'integer',
            'ID (pointing to normtxt) for the primary source of the feature if this assignment is indirect' ],

           ['pos_to_left', 'integer',
            'The position _next to_ the low coordinate of the location' ],

           ['pos_to_right', 'integer',
            'The position _next to_ the high coordinate of the location' ],

           ['rng_id', 'integer',
            'FKEY pointing to range, if more than one segment are present' ],

           ['strand', 'smallint',
            'The strand of the feature, should be 1, -1 or zero' ],

           ] };

    $tables{ "v_feature" } =
    { name  => 'v_feature',
      com   => 'Map out IDs to text for feature table',
      view  =>
"
SELECT f.host_id, nt1.txt AS Host, nt2.txt AS Feature, nt3.txt AS Source, f.pos_to_left as left, f.pos_to_right as right, f.strand
FROM feature f, rna r, normtxt nt1, normtxt nt2, normtxt nt3
WHERE nt1.txt_id = r.acc_id
AND r.rna_id = f.host_id
AND nt2.txt_id = f.feat_id
AND nt3.txt_id = f.src_id
"
};

    $tables{ "queries" } =
    { name  => 'queries',
      com   => 'Shows Postgres SQL statements currently running for ALL databases',
      db    => 'postgres',
      view  =>
"
 SELECT pg_stat_activity.datname, pg_stat_activity.usename, date_trunc('second'::text, now() - pg_stat_activity.query_start) AS query_age, date_trunc('second'::text, now() - pg_stat_activity.backend_start) AS backend_age, btrim(pg_stat_activity.current_query) AS current_query
   FROM pg_stat_activity
  WHERE pg_stat_activity.current_query <> '<IDLE>'::text
  ORDER BY date_trunc('second'::text, now() - pg_stat_activity.query_start), date_trunc('second'::text, now() - pg_stat_activity.backend_start)
"
};

    $tables{ "v_size" } =
    { name  => 'v_size',
      com   => 'Show size of installed postgres databases',
      db    => 'postgres',
      view  =>
"
SELECT datid, datname, 
       pg_size_pretty(pg_database_size(datname)) AS size_on_disk
  FROM pg_stat_database
 ORDER BY pg_database_size(datname) DESC;
"
};

    $tables{ "v_wait" } =
    { name  => 'v_wait',
      com   => 'Find queries that are not immediately returning',
      db    => 'postgres',
      requires => [ 'queries' ],
      view  =>
"
SELECT count(queries.current_query) AS count,
       floor(100::double precision * (avg(date_part('minutes'::text, queries.query_age) * 60::double precision + date_part('seconds'::text, queries.query_age)) / 60::double precision)) / 100::double precision AS minutes,
        queries.current_query
   FROM queries
  GROUP BY queries.current_query
  ORDER BY floor(100::double precision * (avg(date_part('minutes'::text, queries.query_age) * 60::double precision + date_part('seconds'::text, queries.query_age)) / 60::double precision)) / 100::double precision DESC;
"
};

    $tables{ "v_xid" } =
    { name  => 'v_xid',
      com   => 'Report on the transaction IDs for each table',
      db    => 'postgres',
      requires => [ 'queries' ],
      view  =>
"
 SELECT c.relname, ts.spcname AS tablespace,
  c.relpages::double precision / 1000::double precision AS kilopages,
  floor(c.reltuples / 1000::double precision) AS kilotuples,
  age(c.relfrozenxid)::double precision /1000000::double precision AS mega_xid,
  pg_size_pretty(pg_total_relation_size(c.relname::text)) AS disk
   FROM pg_namespace ns, pg_class c
   LEFT OUTER JOIN pg_tablespace ts
     ON (c.reltablespace = ts.oid)
  WHERE ns.oid = c.relnamespace
    AND ns.nspname NOT IN ('pg_catalog', 'information_schema', 'pg_toast')
    AND c.relkind  = 'r'
  ORDER BY c.reltuples DESC
"
};

    # http://stackoverflow.com/questions/2204058/show-which-columns-an-index-is-on-in-postgresql
    # array_agg() not available prior to 8.4
    $tables{ "v_ind" } =
    { name  => 'v_ind',
      com   => 'Summarizes size and location of indices',
      db    => 'postgres',
      view  =>
"
 SELECT c.relname AS Index, tc.relname AS Table, 
        ts.spcname AS tablespace,
  c.relpages::double precision / 1000::double precision AS kilopages,
  floor(c.reltuples / 1000::double precision) AS kilotuples,
  pg_size_pretty(pg_total_relation_size(c.oid)) AS disk
   FROM pg_class tc, pg_namespace ns, pg_index ix, pg_class c
   LEFT OUTER JOIN pg_tablespace ts
     ON (c.reltablespace = ts.oid)
  WHERE ns.oid = c.relnamespace 
    AND ns.nspname NOT IN ('pg_catalog', 'information_schema', 'pg_toast')
    AND tc.oid     = ix.indrelid
    AND c.oid      = ix.indexrelid
  ORDER BY c.reltuples DESC
"
};

    $tables{ "v_tab" } =
    { name  => 'v_tab',
      com   => 'Summarizes activities on tables',
      db    => 'postgres',
      view  =>
"
SELECT relname,  seq_scan, idx_scan, 
       n_tup_ins AS Inserts, n_tup_upd AS Updates, n_tup_del AS Deletes,
       to_char(last_analyze, 'YYYY Mon DD') AS Analyzed,
       to_char(last_vacuum, 'YYYY Mon DD') AS Vacuumed
  FROM pg_stat_all_tables where schemaname = 'public'
 ORDER BY relname
"
};

# Modified from:
# http://stackoverflow.com/questions/1559008/list-stored-functions-using-a-table-in-postgresql

    $tables{ "v_func" } =
    { name  => 'v_func',
      com   => 'Summarizes functions',
      db    => 'postgres',
      view  =>
"
SELECT proname || '(\n  ' || array_to_string(proargnames, ',\n  ') || ')' AS Function, 
       prosrc AS Source
  FROM pg_catalog.pg_namespace n
  JOIN pg_catalog.pg_proc p
    ON pronamespace = n.oid
 WHERE nspname = 'public';
"
};

# http://stackoverflow.com/questions/1109061/insert-on-duplicate-update-postgresql

    $tables{ "add_allele_to_db" } =
    { name  => 'add_allele_to_db',
      com   => 'Perform an upsert on allele data',
      db    => 'postgres',
      args  => [ lid   => 'INT',
                 bid   => 'INT',
                 pid   => 'INT',
                 f     => 'REAL',
                 c     => 'INT',
                 defer => 'INT' ],
      retval    => 'INT',
      language  => 'plpgsql',
      function  =>
"
 BEGIN
   -- Check if there is already an entry in DB
   PERFORM 0 FROM allele WHERE loc_id = lid AND base_id = bid AND pop_id = pid;
 IF FOUND THEN
   -- A row for this location / population / allele base exists
   -- We may need to update it
   IF f IS NULL THEN
     -- Both freq and count should be null
     UPDATE allele SET freq = NULL, count = NULL
      WHERE loc_id = lid AND base_id = bid AND pop_id = pid
        AND freq IS NOT NULL;
     RETURN 1;
   ELSIF c IS NULL THEN
     -- count is null
     UPDATE allele SET freq = f, count = NULL
      WHERE loc_id = lid AND base_id = bid AND pop_id = pid
        AND (freq IS NULL OR freq != f OR count IS NOT NULL);
     RETURN 2;
   ELSE
     -- Both freq and count are not null
     UPDATE allele SET freq = f, count = c
      WHERE loc_id = lid AND base_id = bid AND pop_id = pid
        AND (freq IS NULL OR freq != f OR count IS NULL OR count != c);
     RETURN 3;
   END IF;
 ELSE
   -- We simply need to insert. It is possible that a race condition
   -- occurs here. That is ok. One process wins, the other generates a
   -- key violation warning in the log. We do not have any way to
   -- determine which was 'right' anyway.
   IF defer = 1 THEN
     -- The user wants to do a bulk insert later with COPY
     RETURN 0;
   ELSE
     INSERT INTO allele (loc_id, base_id, pop_id, freq, count)
     VALUES (lid, bid, pid, f, c);
     RETURN 4;
   END IF;
 END IF;
   -- Should end up here with an error
   RETURN -1;
 END;
"
};

    $tables{ "quiet_category_write" } =
    { name  => 'quiet_category_write',
      com   => 'Write a category to the DB without complaining',
      db    => 'postgres',
      args  => [ lid => 'INT',
                 cid => 'INT' ],
      retval    => 'INT',
      language  => 'plpgsql',
      function  =>
"
BEGIN
  PERFORM 0 FROM loc_cat WHERE loc_id = lid AND cat_id = cid;
IF NOT FOUND THEN
  -- It is still possible for a violation to occur here, but we do not care
  -- The unique constraint prevents duplicate rows in DB, and this function
  -- minimizes complaints in the error log
  INSERT INTO loc_cat (loc_id, cat_id) VALUES (lid, cid);
  RETURN 1;
END IF;
  -- The value already exists in database, leave happy
  RETURN 0;
END;
"
};

    $tables{ "quiet_tagval_write" } =
    { name  => 'quiet_tagval_write',
      com   => 'Write a tag/value pair to the DB without complaining',
      db    => 'postgres',
      args  => [ oid => 'INT',
                 tid => 'INT',
                 vid => 'INT' ],
      retval    => 'INT',
      language  => 'plpgsql',
      function  =>
"
BEGIN
  PERFORM 0 FROM tagval WHERE obj_id = oid AND tag_id = tid AND val_id = vid;
IF NOT FOUND THEN
  -- It is still possible for a violation to occur here, but we do not care
  -- The unique constraint prevents duplicate rows in DB, and this function
  -- minimizes complaints in the error log
  INSERT INTO tagval (obj_id, tag_id, val_id) VALUES (oid, tid, vid);
  RETURN 1;
END IF;
  -- The value already exists in database, leave happy
  RETURN 0;
END;
"
};

    $tables{ "quiet_loctoacc" } =
    { name  => 'quiet_loctoacc',
      com   => 'Write a location accession to the DB without complaining',
      db    => 'postgres',
      args  => [ bld => 'VARCHAR(10)',
                 cn  => 'VARCHAR(20)',
                 aid => 'INT' ],
      retval    => 'INT',
      language  => 'plpgsql',
      function  =>
"
BEGIN
  PERFORM 0 FROM loc_to_acc WHERE build = bld AND chr = cn AND acc_id = aid;
  IF FOUND THEN
    -- The value already exists in database, leave happy
    RETURN 0;
  END IF;

  PERFORM 0 FROM loc_to_acc WHERE build = bld AND chr = cn AND acc_id != aid;
  IF FOUND THEN
    -- Oops. A different accession has been assigned to this chromosome
    RETURN 1;
  END IF;
  
  PERFORM 0 FROM loc_to_acc WHERE (build != bld OR chr != cn) AND acc_id = aid;
  IF FOUND THEN
    -- Oops. The accession has been assigned to a different chromosome
    RETURN 2;
  END IF;

  -- Add the accession to the database
  INSERT INTO loc_to_acc (build, chr, acc_id) VALUES (bld, cn, aid);
  RETURN 0;
END;
"
};

    $tables{ "population_pkey" } =
    { name  => 'population_pkey',
      com   => 'Fetch or create a pkey entry for a population',
      db    => 'postgres',
      args  => [ nid => 'INT',
                 sid => 'INT' ],
      retval    => 'INT',
      language  => 'plpgsql',
      function  =>
"
DECLARE rv integer;
BEGIN
  SELECT pop_id INTO rv FROM population WHERE name_id = nid AND src_id = sid;
  IF rv > 0 THEN
    -- done, nothing further to do
    RETURN rv;
  END IF;

  -- Not there - Add a new entry for this name and source
  BEGIN
    INSERT INTO population (name_id, src_id) VALUES (nid, sid);
  EXCEPTION WHEN unique_violation THEN
    -- Do nothing. This happens if another process does the insert before us.
    -- This block is expensive, but allows failover to recover the pkey
  END;

  -- Get the pop_id that should have been auto-generated
  SELECT pop_id INTO rv FROM population WHERE name_id = nid AND src_id = sid;
  IF rv > 0 THEN
    RETURN rv;
  END IF;
  -- Something is horribly wrong if you ever get here!
  RAISE NOTICE 'Failed to generate pop_id';
  RETURN 0;
END;
"
};

    $tables{ "location_pkey" } =
    { name  => 'location_pkey',
      com   => 'Fetch or create a pkey entry for a location',
      db    => 'postgres',
      args  => [ bld => 'VARCHAR(10)',
                 cn  => 'VARCHAR(20)',
                 ptl => 'INT',
                 ptr => 'INT' ],
      retval    => 'INT',
      language  => 'plpgsql',
      function  =>
"
DECLARE rv integer;
BEGIN
  SELECT loc_id INTO rv FROM location WHERE build = bld AND chr = cn
                         AND pos_to_left = ptl AND pos_to_right = ptr;
  IF rv > 0 THEN
    -- done, nothing further to do
    RETURN rv;
  END IF;

  -- Not there - Add a new entry for this position
  BEGIN
    INSERT INTO location (build, chr, pos_to_left, pos_to_right)
                  VALUES (bld, cn, ptl, ptr);
  EXCEPTION WHEN unique_violation THEN
    -- Do nothing. This happens if another process does the insert before us.
    -- This block is expensive, but allows failover to recover the pkey
  END;

  -- Get the loc_id that should have been auto-generated
  SELECT loc_id INTO rv FROM location WHERE build = bld AND chr = cn
                         AND pos_to_left = ptl AND pos_to_right = ptr;
  IF rv > 0 THEN
    RETURN rv;
  END IF;
  -- Something is horribly wrong if you ever get here!
  RAISE NOTICE 'Failed to generate loc_id';
  RETURN 0;
END;
"
};

    $tables{ "text_to_pkey" } =
    { name  => 'text_to_pkey',
      com   => 'Fetch or create a pkey entry for a string',
      db    => 'postgres',
      args  => [ req => 'TEXT' ],
      retval    => 'INT',
      language  => 'plpgsql',
      function  =>
"

-- The text data type is wonderous, allowing arbitrary length
-- strings. However, that is not entirely true ... you can store strings
-- of any length, but you can not INDEX them. Eventually, a large
-- string will cause an index based on it to fail. Table normtxt is designed
-- to allow arbitary string storage with md5() columns for full text searching
-- plus substr() indices for LIKE searches. There is also an integer

DECLARE rv integer;
BEGIN
  SELECT txt_id INTO rv FROM normtxt WHERE md5sum = md5(req) AND txt = req;
  IF rv > 0 THEN
    -- An entry already exists in the database, return it
    RETURN rv;
  END IF;

  -- We need to make an entry, and return the value
  -- We will allow up to 1000 md5 collisions. In reality there will
  -- probably never be a collision
  -- Start by just inserting into mduniq = 1
  INSERT INTO normtxt (md5sum, mduniq, md5upper, txt)
       VALUES (md5(req), 1, md5(upper(req)), req);
  SELECT txt_id INTO rv FROM normtxt WHERE md5sum = md5(req) AND txt = req;
  IF rv > 0 THEN
    -- This job or another managed to record the entry as the first row with
    -- mduniq = 1. Expected 99.some-large-number-of-nines % of the time.
    RETURN rv;
  END IF;

  -- There is another, distinct text entry that shares the MD5 hash!
  -- Find the next available mduniq slot (up to 1000) to place it in

  FOR tryUniq IN 2..1000 LOOP
    BEGIN
      -- Lock to prevent (?) duplicates under different mduniq values
      -- I think locking may not be neccessary
      -- But it should not hurt, and we should rarely ever get here
      LOCK TABLE normtxt IN EXCLUSIVE MODE;
      INSERT INTO normtxt (md5sum, mduniq, md5upper, txt)
           VALUES (md5(req), tryUniq, md5(upper(req)), req);
      SELECT txt_id INTO rv FROM normtxt WHERE md5sum = md5(req) AND txt = req;
      IF rv > 0 THEN
        -- The insert succeeded, return
        RETURN rv;
      END IF;
    END;
  END LOOP;
  -- I NEVER expect to get here, at least with real-life data
  -- If this point is reached, it means 1000 strings with identical md5
  -- hashes have been entered in the database. Tough, you do not get 1001
  -- If you can find 1000 *real-world* strings with identical MD5 hashes
  -- you probably have a celebrated cryptography paper anyway.
  RETURN 0;
END;
"
};

=head3 SQL Performance Thoughts

The database uses unique constraints in many cases to prevent
duplicate entries from building up when sources are reloaded. The goal
was to allow rows to just be slammed into the database and have the
constraint assure that only one instance would be present. This was a
fairly speedy mechanism for update, BUT the log file was getting
hammered with 'violates unique constraint' errors.

The functions below are experiments in trying to find mechanisms that
would efficiently perform inserts only when appropriate.

=head4  add_allele_to_db()

Used to add an allele (with optional frequency and count) to a location

                   Alter DB   No DB Change  
 TIME null+null  : 40 msec    150 -  900 usec
 TIME freq+null  : 43 msec     73 - 1200 usec
 TIME freq+count : 43 msec     75 - 1000 usec

 BEGIN
   -- Check if there is already an entry in DB
   PERFORM 0 FROM allele WHERE loc_id = lid AND base_id = bid AND pop_id = pid;
 IF FOUND THEN
   -- A row for this location / population / allele base exists
   -- We may need to update it
   IF f IS NULL THEN
     -- Both freq and count should be null
     UPDATE allele SET freq = NULL, count = NULL
      WHERE loc_id = lid AND base_id = bid AND pop_id = pid
        AND freq IS NOT NULL;
     RETURN 1;
   ELSIF c IS NULL THEN
     -- count is null
     UPDATE allele SET freq = f, count = NULL
      WHERE loc_id = lid AND base_id = bid AND pop_id = pid
        AND (freq IS NULL OR freq != f OR count IS NOT NULL);
     RETURN 2;
   ELSE
     -- Both freq and count are not null
     UPDATE allele SET freq = f, count = c
      WHERE loc_id = lid AND base_id = bid AND pop_id = pid
        AND (freq IS NULL OR freq != f OR count IS NULL OR count != c);
     RETURN 3;
   END IF;
 ELSE
   -- We simply need to insert. It is possible that a race condition
   -- occurs here. That is ok. One process wins, the other generates a
   -- key violation warning in the log. We do not have any way to
   -- determine which was 'right' anyway.
   INSERT INTO allele (loc_id, base_id, pop_id, freq, count)
   VALUES (lid, bid, pid, f, c);
   RETURN 0;
 END IF;
   -- Should end up here with an error
   RETURN -1;
 END;



 TIME null+null  : 40 msec
 TIME freq+null  : 40 msec
 TIME freq+count : 40 msec
 
 BEGIN
   INSERT INTO allele (loc_id, base_id, pop_id, freq, count)
   VALUES (lid, bid, pid, f, c);
   RETURN 0;
 EXCEPTION WHEN unique_violation THEN
   IF f IS NULL THEN
     UPDATE allele SET freq = NULL, count = NULL
      WHERE loc_id = lid AND base_id = bid AND pop_id = pid
        AND freq IS NOT NULL;
     RETURN 1;
   ELSIF c IS NULL THEN
     UPDATE allele SET freq = f, count = NULL
      WHERE loc_id = lid AND base_id = bid AND pop_id = pid
        AND (freq IS NULL OR freq != f OR count IS NOT NULL);
     RETURN 2;
   ELSE
     UPDATE allele SET freq = f, count = c
      WHERE loc_id = lid AND base_id = bid AND pop_id = pid
        AND (freq IS NULL OR freq != f OR count IS NULL OR count != c);
     RETURN 3;
   END IF;
 END;



=head4  quiet_category_write(lid INT, cid INT)

Used to add a category tag to a location in table loc_cat

http://www.postgresql.org/docs/9.1/static/plpgsql-statements.html#PLPGSQL-STATEMENTS-SQL-NORESULT

 TIME: 130 - 180 us
 WINNER - about 250-500 times faster than error catching

 BEGIN
  PERFORM 0 FROM loc_cat WHERE loc_id = lid AND cat_id = cid;
 IF NOT FOUND THEN
  INSERT INTO loc_cat (loc_id, cat_id) VALUES (lid, cid);
  RETURN 1;
 END IF;
  RETURN 0;
 END;

I really thought this would be the logical way to go...
Apparently exception checking has a large overhead

TIME: 49 msec

 BEGIN
  INSERT INTO loc_cat (loc_id, cat_id) VALUES (lid, cid);
  RETURN 1;
 EXCEPTION WHEN unique_violation THEN
  -- Already in database, no need to do anything
  return 0;
 END;

TIME: 46 msec

 BEGIN
  EXECUTE 'INSERT INTO loc_cat (loc_id, cat_id) VALUES ($1,$2)' 
     USING lid, cid;
  RETURN 1;
 EXCEPTION WHEN unique_violation THEN
  -- Already in database, no need to do anything
  return 0;
 END;

=cut


    $self->bench_end();
    return \%tables;
    
}

=head2 

 Title   : pkey_to_text
 Usage   : my $string = $ml->pkey_to_text( $id )
 Function: Converts a txt_id (PKEY) to a string
 Returns : A string
 Args    : [0] The txt_id of interest

=cut

sub pkey_to_text {
    my $self = shift;
    my $rv = "";
    if (my $id = shift) {
        # 70 usec
        # $self->count_stack();
        $self->bench_start();
        my $sth = $self->{STH}{PKEY2TEXT} ||= $self->dbh->prepare
            ( -name => "Get text for PKEY",
              -sql  => "SELECT txt FROM normtxt WHERE txt_id = ?");
        # die $sth->pretty_print( )."\nEXPLAIN:\n".$sth->explain_text([$id])."\n  ";
        # $self->count_stack();
        $rv = $sth->get_single_value( $id ) || "";
        $self->bench_end();
    }
    return $rv;
}

=head2 

 Title   : cached_pkey_to_text
 Usage   : my $string = $ml->cached_pkey_to_text( $id )
 Function: Converts a txt_id (PKEY) to a string
 Returns : A string
 Args    : [0] The txt_id of interest

=cut

sub cached_pkey_to_text {
    my $self = shift;
    my $pk   = shift || 0;
    unless (defined $self->{OBJ_CACHE}{TEXT}{$pk}) {
        $self->bench_start();
        if (++$self->{CACHE_COUNT}{TEXT} > $self->{CACHE_SIZE}) {
            # This block purges the cache if it grows beyond CACHE_SIZE
            delete $self->{OBJ_CACHE}{TEXT};
            $self->{CACHE_COUNT}{TEXT} = 0;
            $self->{CACHE_PURGES}++;
        }
        my $v = $self->pkey_to_text( $pk );
        $self->{OBJ_CACHE}{TEXT}{$pk} = defined $v ? $v : "";
        $self->bench_end();
    }
    return $self->{OBJ_CACHE}{TEXT}{$pk};
}


# Verify that no duplicates exist:

# select a.txt_id, b.txt_id, a.md5sum, a.txt from normtxt a, normtxt b where a.md5sum = b.md5sum and a.txt_id < b.txt_id;

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub find_text_ids {
    my $self = shift;
    my @rv;
    if (my $txt = shift) {
        $self->bench_start();
        my $get  = $self->{STH}{"GET_PKEY_FOR_UC_TEXT"} ||= $self->dbh->prepare
            ( -name => "Get PKEY for uppercased normtxt",
              -sql  => "SELECT txt_id FROM normtxt ".
              "WHERE md5upper = md5(upper(?)) AND upper(txt) = upper(?)");
        # $get->execute( $txt, $txt );
        @rv = $get->get_array_for_field( $txt, $txt);
        $self->bench_end();
    }
    return @rv;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub text_to_pkey {
    my $self = shift;
    my $txt  = shift;
    my $rv   = 0;
    if (defined $txt && $txt ne '') {
        $self->bench_start();

        # Perl plus several queries = 90 - 120 us (removed code)

        # Single function = 100 - 100 us
        my $get  = $self->{STH}{"TEXT2PKEY_BY_FUNC"} ||= $self->dbh->prepare
            ( -name => "Get PKEY for normtxt",
              -sql  => "SELECT text_to_pkey(?)",
              -ignore => "duplicate key value");
        # The second call is when a duplicate key occurs
        # It is faster this way then putting exception checking in function
        $rv = $get->get_single_value( $txt ) ||
            $get->get_single_value( $txt );
        $self->bench_end();
    }
    return $rv;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub text_to_known_pkey {
    # Will only get already existing entries
    my $self = shift;
    my $txt  = shift;
    # When requested in array context, do a case-insensitive search
    my $noCase = wantarray ? 1 : 0;
    my @rv = ();
    if (defined $txt && $txt ne '') {
        $self->bench_start();
        if ($noCase) {
            my $get  = $self->{STH}{"GET_PKEY_FOR_UCTEXT"} ||=
                $self->dbh->prepare
                ( -name => "Get PKEY for upper normtxt",
                  -sql  => "SELECT txt_id, txt FROM normtxt ".
                  "WHERE upper(substr(txt, 1, 20)) = upper(substr(?, 1, 20)) ".
                  "AND upper(txt) = upper(?)");
            my $rows = $get->selectall_arrayref($txt, $txt);
            # warn "$txt = ".(join(' + ', map { join(',', @{$_}) } @{$rows}) || '-UNDEF-')."\n";
            foreach my $row (@{$rows}) {
                my ($id, $itxt) = @{$row};
                if ($txt eq $itxt) {
                    # Case match
                    $rv[0] = $id;
                }
                # All ids, including any exact matches, are added following
                # This allows iterating from (1..$#rv) to get all IDs
                # or accessing just [0] to see if there was a case match
                push @rv, $id;
            }
        } else {
            my $get  = $self->{STH}{"GET_PKEY_FOR_TEXT"} ||= $self->dbh->prepare
                ( -name => "Get PKEY for normtxt",
                  -sql  => "SELECT txt_id FROM normtxt ".
                  "WHERE md5sum = md5(?) AND txt = ?");
            $rv[0] = $get->get_single_value( $txt );
        }
        $self->bench_end();
    }
    return wantarray ? @rv : $rv[0];
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub cached_text_to_pkey {
    my $self = shift;
    my $txt  = shift || "";
    unless (defined $self->{OBJ_CACHE}{TEXTID}{$txt}) {
        $self->{OBJ_CACHE}{TEXTID}{$txt} = $self->text_to_pkey($txt) || 0;
    }
    return $self->{OBJ_CACHE}{TEXTID}{$txt};
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub tag_query {
    my $self = shift;
    my $args = $self->parseparams( @_ );
    my (@sbits, @binds);
    foreach my $key (qw(TAG VAL)) {
        my $r = $args->{$key};
        next unless ($r);
        my @reqs = ref($r) ? @{$r} : ($r);
        my %ids;
        foreach my $rq (@reqs) {
            if (my $id = $self->text_to_pkey($rq)) {
                $ids{$id} = 1;
            }
        }
        my @idA = keys %ids;
        unless ($#idA == -1) {
            my $col = $key eq 'TAG' ? 'tag_id' : $key eq 'VAL' ? 'val_id' :
                $self->death("Unknown tagval column '$key'");
            push @sbits, "tv.$col IN (".join(', ', @idA).")";
        }
    }
    if ($#sbits == -1) {
        $self->err("No query criteria provided");
        return wantarray ? () : [];
    }
    my $sql = "SELECT tv.obj_id, tv.tag_id, tv.val_id FROM tagval tv WHERE ".
        join(' AND ', @sbits);
    my $sth = $self->dbh->prepare
        ( -name => "tagval query",
          -sql  => $sql,
          -limit => $args->{LIMIT} || 0);
    warn $sth->explain_text() if ($args->{EXPLAIN});
    $sth->execute();
    my $rows = $sth->fetchall_arrayref();
    return wantarray ? map { $_->[0] } @{$rows} : $rows;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub alignment_ids_for_seqid {
    my $self = shift;
    my $id   = shift;
    return () unless ($id);
    $self->bench_start();
    my $get = $self->{STH}{ALIGNMENT_IDS_FOR_SEQID} ||= $self->dbh->prepare
        ( -name => "Get aln_ids for a seq_id",
          -sql  => "SELECT aln_id FROM alignment WHERE seq1 = ? OR seq2 = ?");
    my @rv = $get->get_array_for_field( $id, $id );
    $self->bench_end();
    return wantarray ? @rv : \@rv;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub alignment_ids_for_accession {
    my $self = shift;
    my $acc  = shift;
    return wantarray ? () : [] unless ($acc);
    my $id = $self->text_to_pkey($acc);
    return $self->alignment_ids_for_seqid( $id );
}

=head2 

 Title   : delete_alignment_by_id
 Usage   : $ml->delete_alignment_by_id( $pkey )
 Function: Deletes information from tables alignment and hsp based on aln_id
 Args    : [0] the aln_id PKEY to be deleted

Will clear all entries from alignment and hsp where aln_id = $pkey

=cut

sub delete_alignment_by_id {
    my $self = shift;
    return if ($self->{READ_ONLY});
    $self->bench_start();
    my $dbh  = $self->dbh();
    my $delHsp = $self->{STH}{DELETE_HSP_BY_ID} ||= $dbh->prepare
        ( -name => "Delete HSPs by aln_id",
          -sql  => "DELETE FROM hsp WHERE aln_id = ?");
    my $delAln = $self->{STH}{DELETE_ALIGNMENT_BY_ID} ||= $dbh->prepare
        ( -name => "Delete alignments by aln_id",
          -sql  => "DELETE FROM alignment WHERE aln_id = ?");

    $dbh->begin_work();
    foreach my $id (@_) {
        next unless ($id);
        $delHsp->execute($id);
        $delAln->execute($id);
    }
    $dbh->commit();
    $self->bench_end();
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub delete_tagval_by_id {
    my $self = shift;
    return if ($self->{READ_ONLY});
    $self->bench_start();
    my $dbh = $self->dbh();
    my $del = $self->{STH}{DELETE_TAGVAL_BY_ID} ||= $dbh->prepare
        ( -name => "Delete tag/values by aln_id",
          -sql  => "DELETE FROM tagval WHERE obj_id = ?");
    $dbh->begin_work();
    foreach my $id (@_) {
        next unless ($id);
        $del->execute($id);
    }
    $dbh->commit();
    $self->bench_end();
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub overlapping_locations {
    my $self = shift;
    my ($lft, $rgt, $chr, $build) = @_;
    my @binds;
    my $sTag = "GET_LOCS";
    unless (defined $lft) {
        $self->err("Can not get overlapping locations without at least one coordinate");
        return undef;
    }
    $self->bench_start();
    if (defined $rgt) {
        @binds = ($lft, $rgt);
        if ($chr) {
            $sTag .= "_C";
            push @binds, $chr;
            if ($build) {
                $sTag .= "_B";
                push @binds, $build;
            }
        }
    } else {
        # Single point request
        @binds = ($lft-1, $lft+1);
    }

    my $get = $self->{STH}{$sTag};
    unless ($get) {
        my $sql = "SELECT pos_to_left, pos_to_right, loc_id FROM location ".
            "WHERE pos_to_right > ? AND pos_to_left < ?";
        if (defined $chr) {
            $sql .= " AND chr = ?";
            if (defined $build) {
                $sql .= " AND build = ?";
            }
        }
        $get = $self->{STH}{$sTag} = $self->dbh->prepare
            ( -name => "Get overlapping locations",
              -sql  => $sql );
    }
    $get->execute( @binds );
    my $rows = $get->fetchall_arrayref();
    $self->bench_end();
    return $rows;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub overlapping_alignments {
    my $self = shift;
    my ($l, $r, $a) = @_;
    die "NEEDS WORK";
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub overlapping_alignments_by_id {
    my $self = shift;
    $self->bench_start();
    my ($lft, $rgt, $objid) = @_;
    my $get = $self->{STH}{GET_OVERLAPPING_ALIGN} ||=  $self->dbh->prepare
            ( -name => "Get overlapping alignments",
              -sql  =>
"SELECT a.aln_id, a.seq1, a.score, a.ptl2, a.ptr2
   FROM alignment a
  WHERE a.seq2 = ? AND a.ptr2 > ? AND a.ptl2 < ?
  UNION
 SELECT a.aln_id, a.seq2, a.score, a.ptl1, a.ptr1
   FROM alignment a
  WHERE a.seq1 = ? AND a.ptr1 > ? AND a.ptl1 < ?"

);


    $get->execute($objid, $lft, $rgt, $objid, $lft, $rgt);
    my $rows = $get->fetchall_arrayref();
    $self->bench_end();
    return $rows;
}

# This method just adds an explicit distance to overlapping alignment results
=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub add_distance_to_alignments {
    my $self = shift;
    my ($rows, $l, $r) = @_;
    return if (!$rows || $#{$rows} == -1);
    $self->bench_start();
    $r ||= $l;
    foreach my $row (@{$rows}) {
        # We are expecing the row to be of format
        # [ aln_id, seq_id, ptl, ptr, ... ]
        if ($r <= $row->[2]) {
            # Our range is to the left of the alignment
            push @{$row}, $row->[2] - $r + 1;
        } elsif ($l >= $row->[3]) {
            # Our range is to the left of the alignment
            # We will use a negative distance
            push @{$row}, $row->[3] - $l - 1;
        } else {
            # We overlap
            push @{$row}, 0;
        }
    }
    $self->bench_end();
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub _not_current_text_id {
    my $self = shift;
    unless ($self->{NC_TXT_ID}) {
        my $dtxt = $self->get_text($notCurrentTag);
        $self->{NC_TXT_ID} = $dtxt->pkey();
    }
    return $self->{NC_TXT_ID};
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub overlapping_rnas_by_id {
    my $self = shift;
    $self->bench_start();
    my ($lft, $rgt, $objid, $howbad) = @_;
    # This is the same as overlapping_alignments_by_id(), but will
    # additionally require that the alignment be associated with an RNA
    # and will report the rna ID
    $howbad  = $self->default_howbad() unless (defined $howbad);
    # aln_id, rnaAccID, Score, genoL, genoR, rna_id
    my $get = $self->{STH}{GET_OVERLAPPING_RNAS} ||=  $self->dbh->prepare
        ( -name => "Get overlapping RNAs",
          -sql  =>
"SELECT a.aln_id, a.seq1, a.score, a.ptl2, a.ptr2, r.rna_id
  FROM alignment a, rna r
 WHERE a.seq2 = ? AND r.acc_id = a.seq1 AND a.ptr2 > ? AND a.ptl2 < ?
   AND a.howbad <= ? AND
NOT EXISTS (SELECT nc.val_id FROM tagval nc 
             WHERE nc.obj_id = r.rna_id AND nc.tag_id = ".
          $self->_not_current_text_id().")
 UNION
SELECT a.aln_id, a.seq2, a.score, a.ptl1, a.ptr1, r.rna_id
  FROM alignment a, rna r
 WHERE a.seq1 = ? AND r.acc_id = a.seq2  AND a.ptr1 > ? AND a.ptl1 < ?
   AND a.howbad <= ? AND
NOT EXISTS (SELECT nc.val_id FROM tagval nc 
        WHERE nc.obj_id = r.rna_id AND nc.tag_id = ".
          $self->_not_current_text_id().")"

);

    $get->execute($objid, $lft, $rgt, $howbad,
                  $objid, $lft, $rgt, $howbad);
    my $rows = $get->fetchall_arrayref();
    $self->bench_end();
    return $rows;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub overlapping_rnas_for_flanks {
    my $self = shift;
    $self->bench_start();
    my ($objid, $flanks) = @_;
    # This is the same as overlapping_rnas_by_id(), but takes an array
    # of LFT / RGT pairs as input
    my $dbh = $self->dbh();
    $flanks = [ $flanks ] unless (ref($flanks->[0]));
    my $flankTable = $dbh->list_to_temp_table($flanks, 'integer');
    my $lstxt = $dbh->limit_syntax();
    my $ncid  = $self->_not_current_text_id();
    my $get = $dbh->prepare
        ( -name => "Get overlapping RNAs",
          -sql  =>
"SELECT a.aln_id, a.seq1, a.score, a.ptl2, a.ptr2, r.rna_id
  FROM alignment a, rna r, $flankTable tt
 WHERE a.seq2 = ? AND r.acc_id = a.seq1 AND
EXISTS (SELECT tt.col1 FROM $flankTable tt
        WHERE a.ptr2 > tt.col1 AND a.ptl2 < tt.col2 $lstxt 1) AND
NOT EXISTS (SELECT nc.val_id FROM tagval nc 
        WHERE nc.obj_id = r.rna_id AND nc.tag_id = $ncid)
 UNION
SELECT a.aln_id, a.seq2, a.score, a.ptl1, a.ptr1, r.rna_id
  FROM alignment a, rna r
 WHERE a.seq1 = ? AND r.acc_id = a.seq2 AND
EXISTS (SELECT tt.col1 FROM $flankTable tt
        WHERE a.ptr1 > tt.col1 AND a.ptl1 < tt.col2 $lstxt 1) AND
NOT EXISTS (SELECT nc.val_id FROM tagval nc 
        WHERE nc.obj_id = r.rna_id AND nc.tag_id = $ncid)"

);

    $get->execute($objid, $objid);
    my $rows = $get->fetchall_arrayref();
    $dbh->clear_temp_table( $flankTable );
    $self->bench_end();
    return wantarray ? @{$rows} : $rows;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub closest_rnas {
    my $self = shift;
    $self->bench_start();
    my $args = $self->parseparams( @_ );
    my $objid = $args->{ACCID};
    my $lft   = $args->{LFT};
    my $rgt   = $args->{RGT}      || $lft;
    my $minD  = $args->{MINDIST}  || 0;
    my $maxR  = $args->{MAXRANGE} || 1000000;
    my $block = $args->{BLOCK}    || 50000;
    my @found;
    # First see if there are any RNAs 'on' the location
    my $ontop = $self->overlapping_rnas_by_id($lft, $rgt, $objid);
    push @found, @{$ontop};
    # Now see if there are any nearby but not overlapping
    my ($s, $e) = (1, $block);
    while ($#found == -1 || $minD > $s) {
        # As long as we have no locations found, or if we requested all RNAs
        # within a minimum distance, keep looking
        if ($e > $maxR) {
            # If we have looked beyond the max range, finish searching
            last if ($s > $maxR);
            $e = $maxR;
        }
        # $e = $minD if ($e < $minD);
        my ($ls, $le) = ($lft - $e, $lft - $s);
        my ($rs, $re) = ($rgt + $s, $rgt + $e);
        my @tasks = ([$rs, $re]);
        if ($le >= 0) {
            $ls = 0 if ($ls < 0);
            push @tasks, [$ls, $le];
        }
        my @newFound;
        foreach my $td (@tasks) {
            my ($ts, $te) = @{$td};
            my $alns = $self->overlapping_rnas_by_id($ts, $te, $objid);
            push @newFound, @{$alns};
        }
        unless ($#newFound == -1) {
            $self->add_distance_to_alignments( \@newFound, $lft, $rgt );
            @newFound = sort { abs($a->[6]) <=> abs($b->[6]) } @newFound;
            # No matter what, we will take the closest RNA
            my $d = abs($newFound[0][6]);
            while ($#newFound != -1 && abs($newFound[0][6]) == $d) {
                push @found, shift @newFound;
            }
            if ($minD) {
                # We also want to keep ALL RNAs within a minimum distance
                while ($#newFound != -1 && abs($newFound[0][6]) <= $minD) {
                    push @found, shift @newFound;
                }
            }
        }
        $s += $block;
        $e += $block;            
    }
    $self->bench_end();
    return wantarray ? @found : \@found;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub overlapping_features_by_id {
    my $self = shift;
    $self->bench_start();
    my ($lft, $rgt, $objid) = @_;
    my $get = $self->{STH}{GET_OVERLAPPING_FEATURES} ||=  $self->dbh->prepare
            ( -name => "Get overlapping features",
              -sql  => "SELECT feat_id, src_id, str FROM feature ".
              "WHERE host_id = ? AND pos_to_left < ? AND pos_to_right > ?");

    if (my $w = $rgt + $lft - 1) {
        # As long as the query is not a gap position, pull back the
        # boundaries to get true overlap (not just adjacency)
        $lft++;
        $rgt--;
    }
    $get->execute($objid, $lft, $rgt);
    my $rows = $get->fetchall_arrayref();
    $self->bench_end();
    return $rows;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub cluster_objects {
    my $self  = shift;
    $self->bench_start();
    my $objs  = shift;
    my $range = shift || 0;
    my (@anchored, @unAnchorable);
    foreach my $obj (@{$objs}) {
        if ($obj->can('anchor_id') && $obj->anchor_id) {
            push @anchored, [ $obj->anchor_id,
                              $obj->left_anchor, $obj->right_anchor, $obj ];
        } else {
            push @unAnchorable, $obj;
        }
    }
    @anchored = sort { $a->[0] <=> $b->[0] ||
                           $a->[1] <=> $b->[1] || 
                           $a->[2] <=> $b->[2]  } @anchored;
    my @groups = ( [ shift @anchored ] );
    while (my $adat = shift @anchored) {
        my $prior = $groups[-1][-1];
        if ($prior->[0] == $adat->[0] && $prior->[2] > $adat->[1] - $range) {
            # This features should be in the same group as the last one
            push @{$groups[-1]}, $adat;
        } else {
            # Start a new group
            push @groups, [ $adat ];
        }
    }
    $self->bench_end();
    return wantarray ? @groups : [ \@groups, \@unAnchorable ];
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub cache_size {
    my $self = shift;
    if (my $nv = shift) {
        $self->{CACHE_SIZE} = $nv;
    }
    return $self->{CACHE_SIZE};
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub get_cache {
    my $self = shift;
    my $qkey = join("\t", map { defined $_ ? $_ : "" } @_);
    return $self->{OBJ_CACHE}{$qkey};
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub set_cache {
    my $self = shift;
    my $obj  = shift;
    my $type = shift;
    if (++$self->{CACHE_COUNT}{$type} > $self->{CACHE_SIZE}) {
        # $self->msg("Purging $type cache");
        delete $self->{OBJ_CACHE}{$type};
        $self->{CACHE_COUNT}{$type} = 0;
    }
    my $qkey = join("\t", $type, map { defined $_ ? $_ : "" } @_);
    return $self->{OBJ_CACHE}{$qkey} = $obj || 0;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub population_name_to_pkey {
    my $self = shift;
    my $req  = shift || "";
    unless (defined $self->{POPNAME_PKEY_CACHE}{$req}) {
        # Unvalidated dbSNP = 54669022
        my $arr = $self->{POPNAME_PKEY_CACHE}{$req} = [];
        if (my $popObj = $self->get_population_fast( $req )) {
            push @{$arr}, $popObj->pkey();
        }
    }
    return wantarray ? @{$self->{POPNAME_PKEY_CACHE}{$req}} : 
        $self->{POPNAME_PKEY_CACHE}{$req}[0] || 0;
}

=head2 get_tagged_object_ids

 Title   : get_tagged_object_ids
 Usage   : my @objIds = $ml->get_tagged_object_ids( $value, $tagname )
 Function: Returns database object PKEY ids for tag values / names
 Returns : An array of PKEY integer IDs or an array ref if in scalar context
 Args    : [0] The value of the tag
           [1] Optional name of the tag

The method will return all object IDs matching your query. The name is
optional, with the presumption that tag values will generally be more
specific and therefore of greater utility.

=cut

sub get_tagged_object_ids {
    my $self = shift;
    my ($val, $tag) = @_;
    my @rv;
    return wantarray ? @rv : \@rv if (!defined $val || $val eq '');
    my $valText = $self->get_text( $val );
    my @binds = ($valText->pkey);
    my $sth;
    if ($tag) {
        my $tagText = $self->get_text( $tag );
        push @binds, $tagText->pkey;
        $sth = $self->{STH}{OBJS_FOR_VAL_AND_TAG} ||= $self->dbh->prepare
            ( -name => "Get object ids associated with value and tag",
              -level => 2,
              -sql  => "SELECT obj_id FROM tagval WHERE ".
              "val_id = ? AND tag_id = ?");
    } else {
        $sth = $self->{STH}{OBJS_FOR_VAL} ||= $self->dbh->prepare
            ( -name => "Get object ids associated with value",
              -level => 2,
              -sql  => "SELECT obj_id FROM tagval WHERE val_id = ?");
    }
    @rv = $sth->get_array_for_field( @binds );
    return wantarray ? @rv : \@rv;
}


=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub text_search {
    my $self = shift;
    my ($qry, $isWC, $isCS) = @_;
    return wantarray ? () : [] unless (defined $qry && $qry ne '');
    $self->bench_start();

    # I can not use a prepared statement with placeholders here, the
    # resulting query plan is horribly malformed - basically it
    # defaults to SeqScan. I am instead manually quoting the string
    # myself (via the DBI API).

    # warn "foo";
    my $ucq   = uc($qry);
    my $qlen  = CORE::length($qry);
    my $dbh   = $self->dbh();
    my $b1    = $dbh->quote( ($qlen > 20 ? substr($ucq, 0, 20) : $ucq)."%");

    # Manually escape underscores, so NM_001234 does not get treated
    # as a wildcard

    my $esc   = $dbh->get_info( 14 );
    # $esc should be '\'
    # We need to escape underscores (eg NM_001234)
    $b1 =~ s/([_])/$esc$esc$1/g;

    my $sql;
    if ($isWC) {
        # Looking for a leading match
        my $b2    = $dbh->quote($ucq."%");
        $b2 =~ s/([_])/$esc$1/g;
        $sql = "SELECT txt_id, txt FROM normtxt WHERE ".
            "upper(substr(txt, 1, 20)) LIKE $b1 AND upper(txt) LIKE $b2";
    } else {
        # Looking for an exact match
        my $qq =  $qry;
        $qq    =~ s/\'/_/g;
        my $b2    = $dbh->quote($qq);
        $sql = "SELECT txt_id, txt FROM normtxt WHERE ".
            "upper(substr(txt, 1, 20)) LIKE $b1 AND upper(txt) = ".uc($b2);
        $sql .= " AND txt = $b2" if ($isCS);
    }
    my $sth = $self->dbh->prepare( $sql );
    # $self->preprint( $sth->show_and_explain() );

    $sth->execute();
    my $rv = $sth->fetchall_arrayref();
    $self->bench_end();
    return wantarray ? map { $_->[0] } @{$rv} : $rv;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub poor_score {
    return 80;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub good_score {
    return 100;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub cx_panel_settings {
    my $self = shift;
    my $rv = $self->{CX_PANEL_SETTINGS} ||= {
        featureStaggered     => 'true',
        graphType            => 'Genome',
        featureNameFontColor => 'rgb(29,34,43)',
        trackNameFontColor   => 'rgb(256,100,0)',
        trackNameFontSize    => 10,
        xAxisTickColor       => 'rgb(29,34,43)',
        wireColor            => 'rgba(29,34,43,0.1)',
        autoExtend           => 'true',
        sequenceFill         => '#cccccc',
        infoTimeOut          => 300 * 1000,
        margin               => 2,
        filterSkipNullKeys   => 1,
    };
    if (my $nv = shift) {
        while (my ($key, $val) = each %{$nv}) {
            $rv->{$key} = $val;
        }
    }
    return wantarray ? %{$rv} : $rv;
}

# ::MapLoc - general formatting
=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub cx_genome_panel_data {
    my $self = shift;
    $self->bench_start();
    my $args = $self->parseparams( @_ );
    # Global parameters for the canvas
    my %global = $self->cx_panel_settings();
    my $canvasData = \%global;
    if (my $cdParam = $args->{PARAMS}) {
        while (my ($k, $v) = each %{$cdParam}) {
            $canvasData->{$k} = $v;
        }
    }
    # Panel-specific data:
    my $dataData = { };
    my $tracks = $dataData->{tracks} = [];
    my @tKeys  = map { uc($_) } qw
        (location locations loc locs snp snps
         rna rnas alignment alignments);
    my $src    = $args;
    if (my $tdat = $args->{TRACKS}) {
        $src = $tdat;
        @tKeys = sort keys %{$tdat};
    }
    my %legend;
    foreach my $req (@tKeys) {
        my $dat = $src->{$req};
        next unless $dat;
        my $confTag = uc($req)."CONF";
        $confTag    =~ s/\s/_/g;
        my %conf;
        if (my $confReq = $args->{$confTag}) {
            # User is providing track data
            # Configuration data provided
            my $r = ref($confReq) || "";
            if ($r eq 'ARRAY') {
                %conf = @{$confReq};
            } elsif ($r eq 'HASH') {
                %conf = %{$confReq};
            } else {
                $self->err("Unsure what to do with track configuration for '$req'", $confReq);
            }
        }
        my $name = $args->{uc($req.'NAME')} || $conf{name} || "";
        my $type;
        foreach my $typeCheck ($name, $req, $conf{trackType}) {
            next unless ($typeCheck);
            $type = 
                $typeCheck =~ /(SNP|LOC|VAR|Poly|Mutat)/i ? "Variants" :
                $typeCheck =~ /(RNA|TRANSC)/i ? "RNAs" :
                $typeCheck =~ /(ALIGN)/i ? "Alignments" :
                $typeCheck =~ /(FEAT|MOTIF|DOMAIN)/i ? "Features" :
                $typeCheck =~ /(INFO)/i ? "Information" : "";
            last if $type;
        }
        $type ||= $req;
        $name ||= $req;
        
        $conf{name}    ||= $name;
        $conf{type}    ||= 'box';
        $conf{connect} ||= 'true';
        $conf{data}      = $dat;

        push @{$tracks}, \%conf;

        if ($conf{noLegend}) {
            # Do not record in legend
        } elsif ($type eq 'Features') {
            my $targ = $legend{$type} ||= {};
            foreach my $cx (@{$dat}) {
                my $id = $cx->{id};
                if ($targ->{$id}) {
                    $targ->{$id}{num}++;
                } else {
                    my $sty = $self->style_for_cx($cx);
                    $targ->{$id} = {
                        name => $id,
                        html => sprintf("style='%s'", $sty),
                        num  => 1,
                    };
                }
            }
        } elsif ($type eq 'Information' && $conf{fill}) {
            my @fill = @{$conf{fill}};
            my $targ = $legend{$type} ||= {};
            my $names = $conf{names} || [];
            my @tall;
            for my $c (0..$#fill) {
                my $name = $names->[$c] || "Category ".($c+1);
                my $sty  = "";
                if (my $col = $fill[$c]) {
                    $sty = "style='background-color: $col'";
                }
                my $l = $targ->{$name} ||= {
                    name => $name,
                    html => $sty,
                    num  => 0,
                };
                push @tall, $l;
            }
            foreach my $cx (@{$dat}) {
                if (my $dat = $cx->{data}) {
                    for my $c (0..$#fill) {
                        my $v = $dat->[$c] || 0;
                        # Tally plots use an array of values - ignore them
                        $tall[$c]{num} += $v unless (ref($v));
                    }
                }
            }
            # $self->prebranch($targ);
            
        } elsif ( $type eq 'RNAs') {
            my $targ = $legend{$type} ||= {};
            foreach my $cx (@{$dat}) {
                my $id  = $cx->{type} || 'No Particular Type';
                my $col = $cx->{fill};
                my $l = $targ->{$id} ||= {
                    name => $id,
                    html => "style='background-color: $col'",
                    num  => 0,
                };
                $l->{num}++;
            }
        } elsif ( $type eq 'Variants') {
            my $targ = $legend{$type} ||= {};
            foreach my $cx (@{$dat}) {
                my $id = $cx->{impName} || $cx->{impToken};
                next unless ($id);
                if ($targ->{$id}) {
                    $targ->{$id}{num}++;
                } else {
                    my $sty = $self->style_for_cx($cx);
                    my $desc = "";
                    if (exists $impactRanks->{$cx->{impToken} || ""}) {
                        $desc = $impactRanks->{$cx->{impToken} || ""}{desc};
                    }
                    $targ->{$id} = {
                        name => $id,
                        html => sprintf("style='%s'", $sty),
                        desc => $desc,
                        num  => 1,
                    };
                }
            }
        }
    }
    if (my $srt = $args->{SORT}) {
        my @sArr = ref($srt) ? @{$srt} : split(/\s*[\,\n\r\t]+\s*/, $srt);
        map { $_ = lc($_); s/[^a-z0-9]//gi; } @sArr;
        my @sorter;
        my $def = $#{$tracks} + $#sArr + 10;
        for my $t (0..$#{$tracks}) {
            my $val   = $def + $t;
            my $track = $tracks->[$t];
            my $tn    = lc($track->{name});
            $tn       =~ s/[^a-z0-9]//g;
            for my $s (0..$#sArr) {
                if ($tn =~ /$sArr[$s]/) {
                    $val = $s;
                    last;
                }
            }
            push @sorter, [ $val, $track ];
        }
        $dataData->{tracks} = [ map { $_->[1] } sort 
                                { $a->[0] <=> $b->[0] } @sorter ];
    }
    $self->bench_end();
    return wantarray ? ($dataData, $canvasData, \%legend) : $dataData;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub style_for_cx {
    my $self = shift;
    my $cx   = shift;
    return "" unless ($cx);
    my @bits;
    if (exists $cx->{fill} && $cx->{fill}) {
        push @bits, sprintf("background-color: %s", $cx->{fill});
    }
    if (exists $cx->{outline} && $cx->{outline}) {
        push @bits, sprintf("border: %s solid 1px", $cx->{outline});
    }
    return join('; ', @bits) || "";
}


=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub population_js_extender {
    my $self = shift;
    if (my $cb = shift) {
        $self->{POP_JS_XTND} = $cb;
    }
    return $self->{POP_JS_XTND};
}

our $doneSupportingData;
our $mlSupportName = 'MLsupportData';
=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub serialize_supporting_data {
    my $self = shift;
    my ($cat, $data, $keyReq) = @_;
    my $ser  = $self->serializer();
    my $js   = "";
    unless ($doneSupportingData) {
        # Variable is defined in supporting JS files
        # $js .= "var $mlSupportName = new Array();\n";
        $doneSupportingData = {};
    }
    my $jcat = "${mlSupportName}[".$ser->esc_text( $cat, 'forceQuotes' )."]";
    unless ($keyReq) {
        # $data should be a hash with keys pointing to the data
        my @keyz = sort keys %{$data};
        if ($doneSupportingData->{$cat}) {
            # The hash has already been generated
            # Set each sub key individually
            foreach my $key (@keyz) {
                $js .= $self->serialize_supporting_data
                    ( $cat, $data->{$key}, $key);
            }
        } else {
            # First time seen, we can serialize it all at once
            $js .= "$jcat = " . $ser->obj_to_json
                ($data, 0, { basicArray => 1, tags => 1 } ).";\n";
            map { $doneSupportingData->{$cat}{$_}++ } @keyz;
        }
        return $js;
    }
    unless ($doneSupportingData->{$cat}) {
        # Initialize the object
        $js .= "$jcat = new Object();\n";
    }
    my %keyz   = map { $_ => 1 } (ref($keyReq) ? @{$keyReq} : ($keyReq));
    my $varTxt = "";
    # All the keys will be pointing at the same data
    foreach my $key (sort keys %keyz) {
        unless ($doneSupportingData->{$cat}{$key}++) {
            # This means that we can not update values once they have been
            # set. These values are global to the HTML page, so either way
            # will cause problems. Do not serialize the data until it is
            # ready to be used globally.
            $varTxt .= $jcat."[".$ser->esc_text( $key, 'forceQuotes' )."] = ";
        }
    }
    if ($varTxt) {
        my $jsDat = $ser->obj_to_json($data, 0, {basicArray => 1, tags => 1} );
        $jsDat =~ s/\n$//;
        $js .= $varTxt . $jsDat .";\n";
    }
    return $js;
}
    
our $idsAlreadySupported = {};
=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub population_js_data {
    my $self = shift;
    my $args = $self->parseparams( @_ );
    my $ids  = $args->{IDS} || [];
    my $ser  = $self->serializer();
    my $js   = "";
    my $cats = $args->{CATS};
    my $docat = 0;
    unless ($cats) {
        $cats = {};
        $docat = 1;
    }
    my @stack = @{$ids};
    while (my $id = shift @stack) {
        next if ($idsAlreadySupported->{$id}++);
        my $pop  = $self->get_population( $id );
        unless ($pop) {
            $self->msg_once("Failed to recover population for id $id");
            next;
        }
        my $json = $pop->json_data();
        if (my $cb = $self->population_js_extender()) {
            &{$cb}($self, $pop, $json);
        }
        
        map { $cats->{$_} = 1 } $pop->tag('Category');
        my ($ctag, $col) = $pop->colored_tag();
        $json->{colorTag} = $col || $self->pastel_text_color( $ctag );
        $js .= $self->serialize_supporting_data
            ("Pop", $json, [$id, $pop->name()]);
        if (my $par = $pop->parent()) {
            push @stack, $par->pkey();
        }
    }
    $js .= $self->category_js_data( -cats => $cats ) if ($docat);
    return $js;
}

our $catsAlreadySupported = {};
=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub category_js_data {
    my $self = shift;
    my $args = $self->parseparams( @_ );
    my $cats = $args->{IDS} || $args->{CATS} || [];
    my $ser  = $self->serializer();
    my $js   = "";
    $cats    = [ keys %{$cats} ] if (ref($cats) eq 'HASH');
    foreach my $cat (@{$cats}) {
        next if ($catsAlreadySupported->{$cat}++);
        my $cText = $self->get_text( $cat );
        my $cJson = $cText->json_data();
        $cJson->{colorTag} = $self->user_text_color_or_pastel( $cat );
        $js .= $self->serialize_supporting_data
            ("Cat", $cJson, [$cText->pkey, $cat]);
    }
    return $js;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub feature_js_data {
    my $self = shift;
    my $args = $self->parseparams( @_ );
    my $ids = $args->{IDS} || [];
    my $ser  = $self->serializer();
    my $js   = "";
    foreach my $id (@{$ids}) {
        next if ($idsAlreadySupported->{"F".$id}++);
        my $feat = $self->pkey_to_text( $id );
        my $obj  = $self->get_text($feat);
        unless ($obj) {
            $self->msg_once("Failed to recover Feature for id $id");
            next;
        }
        my $json = $obj->json_data();
        my @name = $obj->tag_values( 'Description' );
        my @coms = $obj->tag_values( 'Comment' );
        my $col  = $self->pastel_text_color( $name[0] || $feat );
        $json->{color} = $col;
        $json->{name}  = $name[0];
        $js .= $self->serialize_supporting_data("Feat", $json, $id);
    }
    return $js;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub rna_js_data {
    my $self = shift;
    my $args = $self->parseparams( @_ );
    my $ids = $args->{IDS} || [];
    my $ser  = $self->serializer();
    my $js   = "";
    foreach my $id (@{$ids}) {
        next if ($idsAlreadySupported->{$id}++);
        my $obj  = $self->get_rna_by_id($id);
        unless ($obj) {
            $self->msg_once("Failed to recover RNA for id $id");
            next;
        }
        foreach my $gene ($obj->each_gene) {
            my $gid  = $gene->pkey();
            next if ($idsAlreadySupported->{$gid}++);
            my $gson = $gene->json_data();
            $js .= $self->serialize_supporting_data
                ("Gene", $gson, [$gid, $gene->acc]);
        }
        my $json = $obj->json_data();
        $js .= $self->serialize_supporting_data
            ("RNA", $json, [$id, $obj->acc]);
    }
    return $js;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub align_js_data {
    my $self = shift;
    my $args = $self->parseparams( @_ );
    my $ids = $args->{IDS} || [];
    my $ser  = $self->serializer();
    my $js   = "";
    foreach my $id (@{$ids}) {
        next if ($idsAlreadySupported->{$id}++);
        my $obj  = $self->get_alignment($id);
        unless ($obj) {
            $self->msg_once("Failed to recover Alignment for id $id");
            next;
        }
        my $json = $obj->json_data();
        $js .= $self->serialize_supporting_data("Align", $json, $id);
    }
    return $js;
}

our $jsVarCounter = 0;
# ::MapLoc - full HTML generation
=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub cx_genome_panel_html {
    my $self = shift;
    $self->bench_start();
    my $args = $self->parseparams( @_ );
    my $cid    = $args->{CANVASID};
    my %vars;
    my ($dd, $cd, $legend) = $self->cx_genome_panel_data( @_ );
    my $pretty    = $args->{PRETTY} ? 0 : undef;
    my $ser = $self->serializer();
    my $lt  = $ser->literal_token();
    my $ed  = {
        mousemove => $lt . "locMouseOver",
        click => $lt . "locMouseClick",
        
    };
    my (%popids, %rnaids, %alnids, %featids, %cats, %genes);
    foreach my $tdat (@{$dd->{tracks} || []}) {
        foreach my $ddat (@{$tdat->{data} || []}) {
            map { $popids{$_} = 1 } keys %{$ddat->{freqs} || {}};
            map { $cats{$_} = 1 } @{$ddat->{cats} || []};
            while (my ($rid, $rdat) = each %{$ddat->{impact} || {}}) {
                $rnaids{$rid} = 1;
                $alnids{ $rdat->{align} || 0} = 1;
                map { $featids{$_} = 1 } keys %{$rdat->{features} || {}};
            }
        }
    }
    my $xtraJs = "";
    $xtraJs .= $self->rna_js_data
        ( -ids => [ sort {$a <=> $b} keys %rnaids ] );
    $xtraJs .= $self->feature_js_data
        ( -ids => [ sort {$a <=> $b} keys %featids ] );
    $xtraJs .= $self->align_js_data
        ( -ids => [ sort {$a <=> $b} keys %alnids ] );
    $xtraJs .= $self->population_js_data
        ( -ids => [ sort {$a <=> $b} keys %popids ], -cats => \%cats );
    $xtraJs .= $self->category_js_data( -cats => [ sort keys %cats ] );
    my $custCol   = $self->user_color_hash();
    my $colStyles = {};
    while (my ($txt, $col) = each %{$custCol}) {
        my $sty = "background-color: $col";
        my $sum = 0; map { $sum += $_ } $self->rgb_values( $col );
        if ($sum < 300) {
            # dark color
            $sty .= "; color: #fff";
        }
        $colStyles->{$txt} = $sty;
    }
    my $supportData = $self->serialize_supporting_data("Color", $colStyles)
        .$xtraJs;
    my $djs = $ser->obj_to_json($dd, 0, {
        basicArray => 1, freqs => 1, impact => 2, tags => 1 } );
    my $cjs = $ser->obj_to_json($cd, $pretty);
    my $ejs = $ser->obj_to_json($ed, $pretty);
    my $w   = $args->{WIDTH} || 600;
    my $h   = $args->{HEIGHT} || 400;
    my $html   = "";
    if (defined $cid) {
        # If the user provides a canvas ID, then generate a full html block
        $cid ||= "MapLocCanvas".++$jsVarCounter;
        $html .= "<canvas id='$cid' width='$w' height='$h'></canvas>\n";
        $html .= "<script>\n$supportData\nnew CanvasXpress('$cid'"
            . ",\n  $djs"
            . ",\n  $cjs"
            . ",\n  $ejs"
            . "\n );</script>\n";
    } else {
        # We will just be building javascript data structures
        my $jv = ++$jsVarCounter;
        my $dVar = $vars{data}   = "cxData$jv";
        my $cVar = $vars{global} = "cxGlobal$jv";
        my $eVar = $vars{events} = "cxEvents$jv";
        $html = join
            ("\n", "var $dVar = $djs;", "var $cVar = $cjs;", "var $eVar = $ejs;")."\n";
    }
    my @lkeys = sort keys %{$legend || {}};
    my $legHTML = "";
    if ($#lkeys != -1 && !$args->{NOLEGEND}) {
        my $spfmt = "<td class='infonum'><span %s>&nbsp;%d&nbsp;&nbsp;</span></td>";
        my $cols  = 4;
        my ($maxC) = sort { $b <=> $a } map { $#{$_} + 1 } map { [values %{$legend->{$_}}] } @lkeys;
        $cols = $maxC if ($cols > $maxC);
        my $tCol  = $cols * 2;
        my $dSty  = sprintf("style='width:%d%%'", 80 / $cols);
        $legHTML .= "<table style='clear:left; font-size:0.8em; width: ${w}px;' class='tab'><tbody>\n";
        foreach my $l (@lkeys) {
            $legHTML .= "<tr><th colspan='$tCol'>$l</th></tr>\n";
            my @dats = sort { uc($a->{name}) cmp uc($b->{name}) }
            values %{$legend->{$l}};
            for my $ld (0..$#dats) {
                my $dat = $dats[$ld];
                $legHTML .= "<tr>\n" unless ($ld % $cols);
                $legHTML .= sprintf($spfmt, $dat->{html}, $dat->{num});
                
                $legHTML .= "  <td $dSty";
                if (my $d = $dat->{desc}) {
                    $legHTML .= sprintf(" title='%s'", $self->esc_xml($d));
                }
                $legHTML .= ">$dat->{name}";
                $legHTML .= "</td>\n";
                $legHTML .= "</tr>\n" if ($ld == $#dats || !(($ld+1) % $cols));
            }
        }
        $legHTML .= "<tbody></table>\n";
        
    }
    $self->bench_end();
    return wantarray ? ($html, $legHTML, \%vars) : $html . $legHTML;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub collapse_cx_set {
    my $self = shift;
    my $cxs  = shift;
    return () unless ($cxs);
    $self->bench_start();
    my $keyz = shift || [qw(id data)];
    my $tagKey = '__tags';
    my %hash;
    foreach my $cx (@{$cxs}) {
        my $key = "";
        foreach my $k (@{$keyz}) {
            my $val = $k eq 'data' ? join
                ('+', map { $_->[0], $_->[1] } @{$cx->{$k}}) : $cx->{$k};
            $key .= "\t" . (defined $val ? $val : '-undef-');
        }
        my $targ = $hash{$key} ||= $cx;
        if ( my $tags = $cx->{tags} ) {
            # Merge the tag hash using a temporary sub-hash
            my $tt = $targ->{$tagKey} ||= {};
            while (my ($k, $arr) = each %{$tags}) {
                map { $tt->{$k}{$_} = 1 } @{$arr};
            }
        }
    }
    my @rv;
    foreach my $cx (values %hash) {
        if ( my $tags = $cx->{$tagKey} ) {
            my $targ =$cx->{tags} = {};
            while (my ($tag, $hash) = each %{$tags}) {
                $targ->{$tag} = [ sort keys %{$hash} ];
            }
            delete $cx->{$tagKey};
        }
        push @rv, $cx;
    }
    $self->bench_end();
    return @rv;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub merge_cx_tags {
    my $self = shift;
    my ($par, $donor) = @_;
    return unless ($donor);
    $self->bench_start();
    my $tags = $par->{tags} ||= {};
    while (my ($k, $arr) = each %{$donor}) {
        if (my $pArr = $par->{$k}) {
            my %u = map { $_ = 1 } (@{$pArr}, @{$arr});
            $tags->{$k} = [ sort keys %u ];
        } else {
            $tags->{$k} = [ @{$arr} ];
        }
    }
    $self->bench_end();
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub variant_nomenclature {
    my $self = shift;
    my ($wt, $otherReq, $l, $r, $type) = @_;
    my $rv = "";
    $type = lc($type || "");
    return $rv unless ($otherReq);
    $wt = uc($wt);
    my $wLen = length($wt);
    if ($wLen > 20) {
        $wt = $wLen . ($type eq 'p' ? 'aa' : 'bp');
    } else {
        $wt =~ s/T/U/g if ($type eq 'r');
    }
    my @other;
    $otherReq = [ $otherReq ] unless (ref($otherReq));
    my %u;
    foreach my $o (@{$otherReq}) {
        my $uco = uc($o);
        $uco =~ s/T/U/g if ($type eq 'r');
        my $oLen = length($uco);
        if ($oLen > 20) {
            $uco = $oLen . ($type eq 'p' ? 'aa' : 'bp');
        }
        $u{$uco} = 1;
    }
    delete $u{$wt};
    @other = sort keys %u;

    return $rv if ($#other == -1 && $type ne 'p');
    my $otxt = join('/', @other);
    my $sz   = 1;
    if (defined $l) {
        # Set up coordinate text
        $sz = $r - $l - 1;
        if (!defined $r || $r == $l + 2) {
            # Single base
            $rv = $l + 1;
        } elsif ($r == $l + 1) {
            # Insertion
            $rv = $l."_".$r."ins";
        } else {
            # Range
            $rv = ($l+1)."_".($r-1);
        }
    }
    if ($otxt eq uc($fsTok)) {
        # Frameshift
        $rv .= 'fs*';
   # } elsif ($wt eq $otxt) {
   #     # The position does not appear variant
   #     $rv .= $wt;
    } elsif ($wt eq '-') {
        if ($sz) {
            # Hm. Should be zero size
            $rv .= '?ins?';
        }
        $rv .= $otxt;
    } elsif (!$sz) {
        # Hm. We have zero size, but a discrete wildtype location
        my $sep = $rv || ">";
        $rv = "?$wt?$sep$otxt";
    } elsif ($otxt eq '-') {
        # We have a deletion
        $rv .= "del$wt";
    } else {
        $otxt ||= ""; # $wt;
        if ($type eq 'p') {
            $rv = "$wt$rv$otxt";
        } else {
            $rv .= "$wt>$otxt";
        }
    }
    $rv = "$type.$rv" if ($type);
    return $rv;
}

# For a query population pop_id recover all category IDs as well as 
# category IDs of parents
=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub all_parent_cat_ids {
    my $self = shift;
    my $pid  = shift || 0;
    my $rv;
    unless ($rv = $self->{CATID_PARENT_ARRAY}{$pid}) {
        $self->benchstart();
        $rv = $self->{CATID_PARENT_ARRAY}{$pid} = [];
        my $get = $self->{STH}{ALLCAT_FOR_POPID} ||= $self->dbh->prepare
            ( -name => "Find category IDs for a pop ID",
              -sql  => "SELECT val_id FROM tagval WHERE obj_id = ? AND tag_id = ".$self->text_to_pkey("Category"), );
        my %seen;
        my @popids = ($pid);
        push @popids, $self->all_parent_pop_ids($pid);
        foreach my $id (@popids ) {
            foreach my $cid ($get->get_array_for_field($id)) {
                if ($cid && !$seen{$cid}++) {
                    # warn "$id = CAT [$cid] ".$self->pkey_to_text($cid)."\n";
                    push @{$rv}, $cid;
                }
            }
        }
        $self->benchend();
    }
    return wantarray ? @{$rv} : $rv;
}

# For a query population pop_id recover all parent pop_ids
=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub all_parent_pop_ids {
    my $self = shift;
    my $pid  = shift || 0;
    my $rv;
    unless ($rv = $self->{POPID_PARENT_ARRAY}{$pid}) {
        $self->benchstart();
        $rv = $self->{POPID_PARENT_ARRAY}{$pid} = [];
        my $get = $self->{STH}{ALLPAR_FOR_POPID} ||= $self->dbh->prepare
            ( -name => "Find parent ID for population ID",
              -sql  => "SELECT par_id FROM population WHERE pop_id = ?", );
        my %seen;
        while ($pid) {
            last if ($seen{$pid}++);
            if ($pid = $get->get_single_value($pid)) {
                push @{$rv}, $pid;
            }
        }
        $self->benchend();
    }
    return wantarray ? @{$rv} : $rv;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub popid_match {
    my $self = shift;
    my $qry  = shift || 0;
    my $lu;
    unless ($lu = $self->{POPID_MATCH_LOOKUP}{$qry}) {
        $self->bench_start();
        $lu = $self->{POPID_MATCH_LOOKUP}{$qry} = {};
        my $get = $self->{STH}{PAR_FOR_POPID} ||= $self->dbh->prepare
            ( -name => "Find parent ID for population ID",
              -sql  => "SELECT par_id FROM population WHERE pop_id = ?", );
        my $pid = $qry;
        while ($pid) {
            last if ($lu->{$pid});
            $lu->{$pid} = 1;
            $pid = $get->get_single_value($pid);
        }
        $self->bench_end;
    }
    foreach my $pid (@_) {
        return 1 if (exists $lu->{$pid});
    }
    return 0;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub css_styles {
    my $styles = "/* Impact Code Styles */";
    # Add styles for the impact codes
    foreach my $imp (sort keys %{$impactRanks}) {
        my $ir = $impactRanks->{$imp};
        my @params = ("");
        if (my $bgc = $ir->{color}) {
            push @params, "background-color: $bgc";
        }
        if (my $fgc = $ir->{fgcol}) {
            push @params, "color: $fgc";
        }
        my $sty = join('; ', @params);
        if ($sty) {
            $sty    .= ";";
            $styles .= sprintf(" .Imp%s { %s }\n", $imp, $sty);
        }
        $styles .= sprintf(" span.Imp%s { %s padding-left: 3px; padding-right: 3px; }\n", $imp, $sty);
    }
    return $styles;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub bulk_frequencies_for_location_ids {
    my $self = shift;


    # NOT BEING USED CURRENTLY - was used by bulk_MAF_for_location_ids()

    $self->bench_start();
    my $lids = shift || [];
    my $rv   = [];
    if (ref($lids)) {
        # Array of loc_ids
        # For ~6000 IDs this works out to ~ 0.7-0.8 msec per ID
        # Does not appear to be cached
        my $sth  = $self->{STH}{FREQ_DATA_FOR_ONE_LID} ||= $self->dbh->prepare
            ( -name => "Get frequency data for a single location ID",
              -sql  => "SELECT a.loc_id, a.pop_id, a.base_id, a.freq, a.count".
              " FROM allele a WHERE a.loc_id = ?", );
        my %done;
        foreach my $lid (@{$lids}) {
            unless ($done{$lid}++) {
                my $rows = $sth->selectall_arrayref( $lid );
                push @{$rv}, @{$rows};
            }
        }
    } else {
        # The request is the name of a temporary table
        # The database appears to heavily cache results from this query
        # (makes benchmarking harder)
        # For ~3000 IDs works out to 11 msec per ID for the first query
        # and 0.1 msec for re-queries.
        my $sth  = $self->dbh->prepare
            ( -name => "Get bulk frequency data for location IDs",
              -sql  => "SELECT a.loc_id, a.pop_id, a.base_id, a.freq, a.count".
              " FROM allele a, $lids tt WHERE a.loc_id = tt.col1", );
        $sth->execute();
        $rv = $sth->fetchall_arrayref();
    }
    $self->bench_end();
    # my $plan = $sth->explain(); warn join("\n", map { $_->[0] } @{$plan}); die $self->show_bench();
    return wantarray ? ($rv) : $rv;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub acceptable_null_frequency {
    my $self = shift;
    if (my $cat = shift) {
        my $cid = $cat =~ /^\d+$/ ? $cat : $self->text_to_pkey($cat);
        my $val = shift;
        if (defined $val && !$val) {
            delete $self->{OK_NULL_FREQ}{$cid};
        } else {
            $self->{OK_NULL_FREQ}{$cid} = $cat;
        }
    }
    return wantarray ? 
        values %{$self->{OK_NULL_FREQ} || {}} : $self->{OK_NULL_FREQ};
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub _maf_from_popBaseHash {
    my $self = shift;
    # $self->bench_start();
    # ~ 150 us, 420 us including child calls
    my ($pbhash, $refPopCache, $rv) = @_;
    $refPopCache ||= {};
    $rv          ||= {};
    my @pids = keys %{$pbhash};
    foreach my $pid (@pids) {
        my $bHash = $pbhash->{$pid};
        my $rbid  = 0;
        # First find out if this population has a reference
        my $rPid = $refPopCache->{$pid};
        # Determine reference ID if this is the first encounter for pid
        $rPid = $refPopCache->{$pid} = $self->reference_pop_id($pid)
            unless (defined $rPid);
        if ($rPid) {
            # Yes, there is a reference population
            if (my $rbHash = $pbhash->{$rPid}) {
                # And we do have allele data. Find the most populous
                # allele. That is what we will compare to
                $rbid = [ sort { $rbHash->{$b} <=> $rbHash->{$a} ||
                                     $a <=> $b } keys %{$rbHash} ];
            }
        }
        # Use -1 for undefined frequencies
        my $maxF = -1;
        if ($rbid) {
            # If a reference allele is defined, take its frequency
            $maxF = $bHash->{$rbid->[0]} || 0;
        } else {
            # We should take the largest allele reported
            ($maxF) = sort { $b <=> $a } values %{$bHash};
        }
        if ($maxF < 0) {
            # Undefined frequency
            if ($rbid) {
                # In some data sets, the reference population has a single
                # allele with undefined frequency, and the sample datasets
                # have one or more alleles with undefined frequency
                # (some TCGA samples, eg STAD TCGA-CG-5721-01)

                # In other data sets, the reference population has TWO
                # frequency-undefined alleles 
                # (7.GRCh37:140449150 GBM TCGA-02-0115-10)

                # I think it is best that if frequencies are not
                # defined, that the reported MAF should be zero if
                # the enumerated alleles are the same as those of
                # the reference.

                my $atxt = join(' ', sort {
                    $bHash->{$b} <=> $bHash->{$a} ||
                        $a <=> $b } keys %{$bHash}) || "";
                my $rtxt = join(' ', @{$rbid}) || "";
                
                if ($atxt eq $rtxt) {
                    # This population has the same alleles assigned
                    # as the reference population. Consider it
                    # 'uninteresting', with a zero MAF
                    # Setting the maximum frequency to 1 results in MAF=0
                    $maxF = 1;
                } else {
                    # A different allele assignment, at least in part
                    # Set the MAF as -1 to indicate undefined
                    # which (in some cases) can be allowed through filters
                    $maxF = 2;
                }
            } else {
                # A maximum frequency of 2 will result in a minor frequency
                # of -1
                $maxF = 2;
            }
        }
        $rv->{$pid} = 1 - $maxF;
    }
    # $self->bench_end();
    return $rv;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub bulk_MAF_for_location_ids {
    my $self = shift;
    my ($lids, $freqs) = @_;
    if ($freqs) {
        # The user is keeping a cache of recovered frequencies
        #my @newlids;
        #foreach my $lid (@{$lids}) {
        #    unless ($freqs->{$lid}) {
        #        push @newlids, $lid;
        #        # $freqs->{$lid} = {};
        #    }
        #}
        #return wantarray ? ([], $freqs) : [] if ($#newlids == -1);
        #$lids = \@newlids;
    } else {
        # Set up new hash
        $freqs = { }; # map { $_ => {} } @{$lids} };
    }
    $self->bench_start("Pre-cached");
    my $bulkSrc  = $self->bulk_cache_source( 'allele' );
    # ?? I believe this is now being handled by PopulationFilter:
    # my %okIfNull = %{$self->acceptable_null_frequency() || {}};
    
    my (%refPops, %needLid);
    foreach my $lid (@{$lids}) {
        if ($freqs->{$lid}) {
            # Already calculated, nothing to do
        } elsif (my $alleleRows = $bulkSrc->{$lid}) {
            # Basic information already recovered
            my %pbHash;
            foreach my $row (@{$alleleRows}) {
                my ($bid, $pid, $f, $c) = @{$row};
                $pbHash{$pid}{$bid} = defined $f ? $f : -1;
            }
            $self->_maf_from_popBaseHash
                (\%pbHash, \%refPops, $freqs->{$lid} = {} );
        } else {
            # Need to get frequency info
            $needLid{$lid} = 1;
        }
    }
    $self->bench_end("Pre-cached");
    my @newLid = keys %needLid;
    unless ($#newLid == -1) {
        # warn scalar(@newLid);
        $self->bench_start("Queried");
        # There are some loc_ids that we need frequency data for
        my $dbh     = $self->dbh();
        my $tempTab = $dbh->list_to_temp_table(\@newLid, 'integer');
        my $get     = $self->dbh->prepare
            ( -name => "Get alleles for location in bulk",
              -sql  => "SELECT loc_id, base_id, pop_id, freq, count ".
              "FROM allele WHERE loc_id ".$self->_subquery_clause
              ( 'col1', $tempTab));
        my $rows = $get->selectall_arrayref( );
        my %pbHashes;
        foreach my $row (@{$rows}) {
            my $lid = shift @{$row};
            push @{$bulkSrc->{$lid}}, $row;
            my ($bid, $pid, $f, $c) = @{$row};
            $pbHashes{$lid}{$pid}{$bid} = defined $f ? $f : -1;
        }
        foreach my $lid (@newLid) {
            $self->_maf_from_popBaseHash
                ($pbHashes{$lid} || {}, \%refPops, $freqs->{$lid} = {} );
        }
        $self->bench_end("Queried");
    }

    # This was originally done in bulk from a temporary table
    # But itterative queries for each loc_id is generally faster
    #my $get = $self->{STH}{ALLELES_FOR_LOC} ||= $self->dbh->prepare
    #    ( -name => "Get alleles for location",
    #      -sql  => "SELECT base_id, pop_id, freq, count ".
    #      "FROM allele WHERE loc_id = ?");
    #my (%refPops, %doneLid);
    #foreach my $lid (@{$lids}) {
    #    next if ($doneLid{$lid}++);
    #    my $alleleRows = $bulkSrc->{$lid} ||= $get->selectall_arrayref( $lid );
    #    my %pbHash;
    #    foreach my $row (@{$alleleRows}) {
    #        my ($bid, $pid, $f, $c) = @{$row};
    #        $pbHash{$pid}{$bid} = defined $f ? $f : -1;
    #    }
    #    $self->_maf_from_popBaseHash(\%pbHash, \%refPops, $freqs->{$lid} ||= {});
    #}
    my $newPids = [keys %refPops];
    return wantarray ? ($newPids, $freqs) : $newPids;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub reference_pop_id {
    my $self = shift;
    $self->bench_start();
    # 415 us
    my $pid  = shift || 0;
    my $sth  = $self->{STH}{FIND_REF_POP_ID} ||= $self->dbh->prepare
        ( -name => "Get refrence population ID for a query pop_id",
          -sql  => "SELECT rp.pop_id FROM population rp, population qp, tagval tv WHERE rp.name_id = tv.val_id AND qp.pop_id = ? AND rp.src_id = qp.src_id AND tv.obj_id = qp.pop_id AND tv.tag_id IN (".join(',', map { $self->text_to_pkey($_) } ('Normal Population','Reference Population') ).")");
    my @rv = $sth->get_array_for_field( $pid );
    #$self->preprint("reference_pop_id($pid) =\n".$sth->pretty_print( $pid )."= ".join(',', @rv));
    if ($#rv == -1) {
        # If the population has no reference, but is itself a reference,
        # then return its own ID
        my $selfCheck  = $self->{STH}{CHECK_SELF_REF_POP} ||=
            $self->dbh->prepare
            ( -name => "Check if a population is a reference for others",
              -limit => 1,
              -sql  => "SELECT qp.pop_id FROM population qp, tagval tv WHERE qp.name_id = tv.val_id AND qp.pop_id = ? AND tv.tag_id IN (".join(',', map { $self->text_to_pkey($_) } ('Normal Population','Reference Population') ).")");
        #$self->preprint("SELF CHECK=\n".$selfCheck->pretty_print( $pid ));
        push @rv, $pid if ($selfCheck->get_single_value( $pid ));
    }
    $self->bench_end();
    return wantarray ? @rv : $rv[0] || 0;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub bulk_rnas_for_location_ids {
    my $self = shift;
    my ($lids, $rnaH, $howbad) = @_;
    if ($rnaH) {
        # The user is keeping a cache of recovered overlaps
        my @newlids;
        foreach my $lid (@{$lids}) {
            unless ($rnaH->{$lid}) {
                push @newlids, $lid;
            }
        }
        return wantarray ? ($rnaH, [], []) : $rnaH if ($#newlids == -1);
        $lids = \@newlids;
    } else {
        # Set up new hash
        $rnaH = {  };
    }
    my $defHB   = $self->default_howbad();
    my $bulkSrc = $self->bulk_cache_source( 'nearby_rnas' );
    if (defined $howbad) {
        $bulkSrc = undef if ($howbad != $defHB);
    } else {
        $howbad = $defHB;
    }
    $self->bench_start();
    my $dbh   = $self->dbh;
    my $getLR = $dbh->prepare
        ( -name => "Get LFT/RGT information for LIDs",
          -sql  => "SELECT build, chr, pos_to_left, pos_to_right".
          " FROM location where loc_id = ?");
    my %clusters;

    foreach my $lid (@{$lids}) {
        my $rows = $getLR->selectall_arrayref( $lid );
        if (my $r = $rows->[0]) {
            push @{$clusters{$r->[0]}{$r->[1]}}, $rnaH->{$lid} ||= {
                l   => $r->[2], 
                r   => $r->[3],
                lid => $lid,
                imp => {},
            };
        }
    }
    # "Not Current" filtration will be done in Perl
    my $ncid = $self->_not_current_text_id();
    my $isNC = $dbh->prepare
        ( -name => "Find not current RNAs",
          -sql  => "SELECT val_id FROM tagval ".
          "WHERE obj_id = ? AND tag_id = $ncid" );

    # Recover HSPs depending on if seq1 or seq2 is the genomic sequence:
    my @ghsp;
    $ghsp[0] = $dbh->prepare
        ( -name => "Get genomic HSPs for 0-index",
          -sql  => "SELECT ptl1, ptr1 FROM hsp WHERE aln_id = ?" );
    $ghsp[1] = $dbh->prepare
        ( -name => "Get genomic HSPs for 1-index",
          -sql  => "SELECT ptl2, ptl2+ptr1-ptl1 FROM hsp WHERE aln_id = ?" );

    my $flank   = $self->default_loc_to_rna_distance();
    my $flPlus1 = $flank + 1; # Used in L/R coordinate space for proximity

    my $spliceSz = 2; # number of base pairs to consider as being "splice site"
    my $splDelta = $spliceSz - 1; # If L/R difference is less than this value
                                  # then we are within the $spliceSz region
    my (%textIds, %rids);
    while (my ($bld, $chrH) = each %clusters) {
        while (my ($chr, $rArr) = each %{$chrH}) {
            my $gaccid = $self->loc_to_accid( $bld, $chr );
            $textIds{$gaccid} = 1;
            # Organize the loci in ascending genomic coordinate space
            my @lrls = sort { $a->{l} <=> $b->{l} ||
                                  $a->{r} <=> $b->{r} } @{$rArr};
            my @hsps = map { [$_->{l}, $_->{r} ] } @lrls;
            # Get the alignments overlapping this region:
            my $rngAlns = $self->_alignments_overlapping_regions
                ( \@hsps, $gaccid, $howbad );

            # Filter old RNA entries
            # Add genomic HSP ranges
            my %rnaAlns;
            while (my ($aid, $adat) = each %{$rngAlns}) {
                my $rid  = $adat->{rid};
                my $rdat = $rnaH->{RNAS}{$rid};
                unless ($rdat) {
                    my $raccid = $adat->{acc};
                    my $acc    = $self->pkey_to_text( $raccid );
                    my $ns     = 
                        $acc =~ /^ENS/       ? 'ENST' : 
                        $acc =~ /^[NX][MR]_/ ? 'RSR'  : 'AR';
                    $rdat = $rnaH->{RNAS}{$rid} = {
                        rid   => $rid,
                        acc   => $acc,
                        accid => $raccid,
                        aln   => {},
                        ns    => $ns,
                        NC    => $isNC->get_single_value( $raccid ) || 0,
                    };
                    # $self->msg_once("FOUND: $acc, rna_id = $rid [$rdat->{NC}]");
                }
                # Skip if RNA is listed as not current
                next if ($rdat->{NC});
                $rids{$rid}  = 1;
                my $shareAln = $rnaH->{ALNS}{$aid};
                unless ($shareAln) {
                    $shareAln = $rnaH->{ALNS}{$aid} = $adat;
                    my $hspSth = $ghsp[ $adat->{gin} ];
                    $hspSth->execute( $aid );
                    my $hspR = $hspSth->fetchall_arrayref();
                    my @hsps = sort { $a->[0] <=> $b->[0] ||
                                          $a->[1] <=> $b->[1] } @{$hspR};
                    # $self->prebranch(\@hsps);
                    $adat->{hsps} = \@hsps;
                    my ($l, $r) = ($hsps[0][0], $hsps[-1][1]);
                    $rdat->{l} = $l if (!defined $rdat->{l} ||
                                        $rdat->{l} > $l);
                    $rdat->{r} = $r if (!defined $rdat->{r} ||
                                        $rdat->{r} < $r);
                }
                $rdat->{aid}{$aid} = 1;
                my $locRdat = $rnaAlns{$rid} ||= {
                    l    => $shareAln->{l},
                    r    => $shareAln->{r},
                    str  => $shareAln->{str},
                    rid  => $rid,
                };
                $locRdat->{l} = $shareAln->{l} if
                    ($locRdat->{l} > $shareAln->{l});
                $locRdat->{r} = $shareAln->{r} if
                    ($locRdat->{r} < $shareAln->{r});
                #warn "aln_id=$aid [$shareAln->{l}..$shareAln->{r}] rna_id=$shareAln->{rid}\n";
            }
            #warn "SELECT * FROM alignment WHERE aln_id IN (".join(',', keys %{$rngAlns}).");\n";
            # Sort RNAs by position on chromosome
            my @rnaInfo = sort { $a->{l} <=> $b->{l} ||
                                     $a->{r} <=> $b->{r} } values %rnaAlns;
            # $self->prebranch(\@rnaInfo);
            $self->bench_start("Estimate Impact");
            foreach my $ldat (@lrls) {
                my ($ll, $lr, $lid) = ($ldat->{l}, $ldat->{r}, $ldat->{lid});

                $ldat->{gid} ||= $gaccid;
                my $limps = $ldat->{imp} ||= {};
                my @consider = @rnaInfo;
                my @remove;
                for my $i (0..$#consider) {
                    my $locRdat = $consider[$i];
                    my $rid  = $locRdat->{rid};
                    # $self->msg("[!]", "$rid [$ll..$lr] vs [$locRdat->{r}..$locRdat->{l}]") if ($lid == 44030822);
                    if ($locRdat->{r} - 2 < $ll) {
                        # The RNA is "to the left" and non-overlapping
                        # This RNA will not overlap with any more locations
                        # Remove the RNA from consideration
                        unshift @remove, $i;
                        next;
                    } elsif ($lr - 2 < $locRdat->{l}) {
                        # The RNA is "to the right" and non-overlapping
                        # The RNA may overlap with subsequent locations
                        # This location will not overlap with any other RNAs
                        last;
                    }

                    # We are likely overlapping with the RNA, though
                    # we may not be if there are two alignments on the
                    # chromosome (we could be between them). Consider
                    # each alignment in turn

                    my $alns = $locRdat->{aln};
                    unless ($alns) {
                        my $rdat = $rnaH->{RNAS}{$rid};
                        my @adats = sort { $a->{l} <=> $b->{l} ||
                                               $a->{r} <=> $b->{r} } 
                        map { $rnaH->{ALNS}{$_} } keys %{$rdat->{aid}};
                        $alns = $locRdat->{aln} = \@adats;
                    }

                    my $isNear = 0;
                    for my $j (0..$#{$alns}) {
                        my $adat = $alns->[$j];
                        my $str  = $adat->{str} || 0;
                        my @hsps = @{$adat->{hsps}};
                        if ($hsps[0][0] - 2 > $lr) {
                            # Location is to the left, not overlapping
                            if ($hsps[0][0] - $flPlus1 < $lr &&
                                (!$isNear || ref($isNear))) {
                                # We are close to the alignment, though
                                $isNear ||= [];
                                push @{$isNear}, $adat->{aid};
                            }
                            # The loc will also not overlap with
                            # following alns. However, it might still
                            # be NEAR to them, so we will next instead
                            # of last.
                            next;
                        } elsif ($ll - 2 > $hsps[-1][1]) {
                            # Location is to the right, not overlapping
                            # It may overlap with a following alignment
                            if ($ll - $flPlus1 > $hsps[-1][1] &&
                                (!$isNear || ref($isNear))) {
                                # We are close to the alignment, though
                                $isNear ||= [];
                                push @{$isNear}, $adat->{aid};
                            }
                            next;
                        }
                        # We are within the RNA boundaries
                        $isNear = -1;
                        my $imp = "INT";
                        for my $k (0..$#hsps) {
                            my $delta1 = $hsps[$k][0] - $lr;
                            if ($delta1 > -2) {
                                # Does not overlap, this exon and all
                                # remaining exons are to the right of
                                # the location
                                if ($delta1 < $splDelta) {
                                    # But we are in the splice site (2bp)
                                    $imp = $str > 0 ?
                                        'SP3' : $str < 0 ? 'SP5' : 'SPL';
                                    # warn "$imp: $hsps[$k][0] [$ll^$lr] $hsps[$k][1]\n";
                                }
                                last;
                            }
                            my $delta2 = $ll - $hsps[$k][1];
                            if ($delta2 > -2) {
                                # Does not overlap, this exon is to the left
                                # of the location. Other exons may overlap
                                if ($delta2 < $splDelta) {
                                    # But we are in the splice site
                                    $imp = $str > 0 ?
                                        'SP5' : $str < 0 ? 'SP3' : 'SPL';
                                    # warn "$imp: $hsps[$k][0] [$ll^$lr] $hsps[$k][1]\n";
                                    last;
                                }
                                next;
                            }
                            # We are inside an exon
                            # TO DO: Resolve this down to UT3 and UT5

                            # Without alleles - and a fair bit more
                            # CPU time - the best we can do here is
                            # label the impact as RNA. For this
                            # reason, when these impacts are passed to
                            # a PopulationFilter they should be tested
                            # with 'loose' methods

                            $imp = 'RNA';
                            last;
                        }
                        $ldat->{aln}{$adat->{aid}} = $imp;
                    }
                    if (ref($isNear)) {
                        # We did not find any alignments that overlapped the
                        # location. However, it was close to one or more
                        map { $ldat->{aln}{$_} = 'GEN' } @{$isNear};
                    }
                }
                # If any RNAs are no longer capable of overlapping locations
                # remove them
                if ($#remove != -1) {
                    # my @rid = map { $rnaInfo[$_]{rid} } @remove; my @acc = map { $rnaH->{RNAS}{$_}{acc} } @rid; $self->msg("[DEBUG]","Removing RNAs ".join(',', @acc). " indices [".join(',',@remove)."] rna_id IN (".join(',',@rid).") after clearing loc_id=$lid [$ll,$lr]. Was = ".join(' + ', map { $rnaH->{RNAS}{$_->{rid}}{acc}  } @rnaInfo));
                    # These should already be sorted, but bad things happen
                    # if they are not. Safety sort desceding:
                    foreach my $ind (sort { $b <=> $a } @remove) {
                        splice(@rnaInfo, $ind, 1);
                    }
                }
                if (my $alnImps = $ldat->{aln}) {
                    # This location is near an RNA
                    my %byNS;
                    while (my ($aid, $tok) = each %{$alnImps}) {
                        my $adat = $rnaH->{ALNS}{$aid};
                        my $rid  = $adat->{rid};
                        my $rdat = $rnaH->{RNAS}{ $rid };
                        foreach my $ns ($rdat->{ns}, 'ALL') {
                            $byNS{$ns}{$tok} = 1;
                        }
                        # Populate the bulk cache for this result in case
                        # it is needed later
                        push @{$bulkSrc->{$lid}},
                        [ $aid, $rdat->{accid}, $adat->{sc}, 
                          $adat->{l}, $adat->{r}, $rid ] if ($bulkSrc);
                    }
                    while (my ($ns, $h) = each %byNS) {
                        $ldat->{imp}{$ns} = $self->heaviest_impact(keys %{$h});
                    }
                } else {
                    # Note that there are no nearby RNAs in the bulk source:
                    $bulkSrc->{$lid} ||= [] if ($bulkSrc);
                }
                $ldat->{imp}{ALL} ||= 'GEN';
            }
            $self->bench_end("Estimate Impact");
        }
    }
    # $self->prebranch($bulkSrc->{44030822});$self->prebranch($bulkSrc->{45863557});

    # $self->prebranch($rnaH);
    my @allRids = keys %rids;
    $self->bulk_tag_values( \@allRids, 'textToo' );
    # $self->bulk_text_values( [keys %textIds] );
    # $self->msg_once("[TODO]","Need to port bulk caching methods from _bulk_rnas_for_flank_table()");
    $self->bench_end();
    return ($rnaH, ['NOT USED'], 'NOT_USED', \@allRids );
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub _alignments_overlapping_regions {
    my $self = shift;
    $self->bench_start();
    my ($hsps, $gaccid, $howbad, $alns) = @_;
    my $dbh = $self->dbh();
    $alns ||= {};
    my @srtHsp  = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$hsps};
    my @ranges  = ( [$srtHsp[0][0], $srtHsp[0][1]] );
    my $flank   = $self->default_loc_to_rna_distance();
    my $twoFl   = $flank * 2; # Twice the flank distance, used for clustering
    for my $l (1..$#srtHsp) {
        if ($srtHsp[$l][0] - $twoFl > $ranges[-1][1]) {
            # New range
            push @ranges, [ $srtHsp[$l][0], $srtHsp[$l][1] ];
        } else {
            # Extend prior range
            $ranges[-1][1] = $srtHsp[$l][1];
        }
    }
    foreach my $rdat (@ranges) {
        my ($l, $r) = @{$rdat};
        $l -= $flank;
        $r += $flank;
        
        # Rather than preparing this statement, we are
        # creating it explicitly. This appears to help the
        # query planner make a proper plan

    my $sql = <<SQL;
SELECT a.aln_id, r.acc_id, a.score, a.ptl1, a.ptr1, r.rna_id, a.strand, 0 AS genind
  FROM alignment a, rna r
 WHERE a.seq1    = $gaccid
   AND a.ptr1   >= $l
   AND a.ptl1   <= $r
   AND a.howbad <= $howbad
   AND a.seq2    = r.acc_id
UNION
SELECT a.aln_id, r.acc_id, a.score, a.ptl2, a.ptr2, r.rna_id, a.strand, 1 AS genind
  FROM alignment a, rna r
 WHERE a.seq2    = $gaccid
   AND a.ptr2   >= $l
   AND a.ptl2   <= $r
   AND a.howbad <= $howbad
   AND a.seq1    = r.acc_id
SQL
                
        my $sth  = $dbh->prepare
            ( -name => "Get RNAs overlapping a genomic range",
              -sql  => $sql );
        $sth->execute( );
        my $rows = $sth->fetchall_arrayref();
        foreach my $row (@{$rows}) {
            my ($alnid, $raccid, $sc, $rl, $rr, $rid, $str, $gind) = @{$row};
            $alns->{$alnid} ||= {
                aid => $alnid,
                gin => $gind,
                sc  => $sc,
                l   => $rl,
                r   => $rr,
                rid => $rid,
                str => $str,
                acc => $raccid,
            };
        }
    }
    $self->bench_end();
    return $alns;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub OLD_bulk_rnas_for_location_ids {
    my $self = shift;
    $self->msg("[DEBUG]","Testing new method"); $self->NEW_bulk_rnas_for_location_ids(@_);
    my ($lids, $rnaH, $howbad) = @_;
    if ($rnaH) {
        # The user is keeping a cache of recovered overlaps
        my @newlids;
        foreach my $lid (@{$lids}) {
            unless ($rnaH->{$lid}) {
                push @newlids, $lid;
            }
        }
        return wantarray ? ($rnaH, [], []) : $rnaH if ($#newlids == -1);
        $lids = \@newlids;
    } else {
        # Set up new hash
        $rnaH = {  };
    }
    $self->bench_start();

    my $dbh   = $self->dbh;
    my $ttab  = $dbh->list_to_temp_table($lids, 'integer');
    my %alns = ();
    my ($rows, $allRids, $alnDat);
    if (0) {
        # Lots of trouble getting the bulk query to run efficiently
        # Postgresql becomes confounded for range queries when it does
        # not have information about the ranges in advance
        $self->bench_start("Location Details");
        my $flank   = $self->default_loc_to_rna_distance();
        my $lSth = $dbh->prepare
            ( -name => "Get location information for LIDs",
              -sql  => "SELECT l.loc_id, la.acc_id, l.pos_to_left, l.pos_to_right FROM $ttab tt, location l, loc_to_acc la WHERE l.loc_id  = tt.col1 AND la.chr = l.chr AND la.build = l.build");
        $lSth->execute();
        my $lrows = $lSth->fetchall_arrayref();
        $self->bench_end("Location Details");
        my %byAcc;
        foreach my $row (@{$lrows}) {
            my ($lid, $accid, $l, $r) = @{$row};
            my $targ = $byAcc{$accid} ||= { };
            push @{$targ->{bait}}, [$l - $flank, $r + $flank];
            push @{$targ->{locid}}, [$l, $r, $lid];
        }
        $self->prebranch(\%byAcc); die;
        while (my ($accid, $data) = each %byAcc) {
            my $baits = $self->add_hsps( $data->{bait} );
            foreach my $bait (@{$baits}) {
                
            }
            warn "$accid = ".$self->branch($baits);
        }
        die "Testing";
    } else {
        # Breaking out the query in two parts to reduce SQL complexity
        # This seems to help Postgres pick the full [ID, LFT, RGT] index
        # in query planning
        $self->bench_start('Normalize locations');
        my @pcols = ("l.loc_id AS obj_id",
                     "la.acc_id AS acc_id",
                     "l.pos_to_left AS ptl",
                     "l.pos_to_right AS ptr", );
        my $prep  = <<PREPSQL;
location l, loc_to_acc la, $ttab tt
 WHERE l.loc_id  = tt.col1
   AND la.chr    = l.chr
   AND la.build  = l.build
PREPSQL

        my $ptab  = $dbh->select_to_temp_table( \@pcols, $prep );
        $dbh->_fast_add_index( $ptab, [ 'obj_id' ] );
        $dbh->_fast_add_index( $ptab, [ 'acc_id' ] );
        $self->bench_end('Normalize locations');
        
        ($rows, $allRids, $alnDat) = $self->_bulk_rnas_for_flank_table
            ( $ptab, $howbad, $rnaH );
        
        unless ($#{$alnDat} == -1) {
            $self->bench_start('Get HSPs');
            my $att  = $dbh->list_to_temp_table
                ($alnDat, ['integer','integer','integer','integer']);
            my $asql = <<SQL;
SELECT h.aln_id, h.ptl1 AS l, h.ptr1 AS r
  FROM hsp h, alignment a, $att tt
 WHERE a.aln_id = tt.col1
   AND a.seq1   = tt.col2
   AND h.aln_id = a.aln_id
UNION
SELECT h.aln_id, h.ptl2 AS l, h.ptl2 + h.ptr1 - h.ptl1 AS r
  FROM hsp h, alignment a, $att tt
 WHERE a.aln_id = tt.col1
   AND a.seq2   = tt.col2
   AND h.aln_id = a.aln_id
SQL

            my $asth = $dbh->prepare
            ( -name => "Get genomic locations for alignment IDs",
              -sql  => $asql);
            # warn $asth->pretty_print().$asth->explain_text();
            $asth->execute();
            my $arows = $asth->fetchall_arrayref();
            foreach my $row (@{$arows}) {
                my ($aid, $l, $r) = @{$row};
                push @{$alns{$aid}}, [$l, $r];
            }
            while (my ($aid, $arr) = each %alns) {
                $alns{$aid} = [ sort { $a->[0] <=> $b->[0] } @{$arr} ];
            }
            $self->bench_end('Get HSPs');
        }
    }
    my $range = 3;
    foreach my $lid (@{$lids}) {
        my $ldat = $rnaH->{$lid};
        unless ($ldat) {
            $rnaH->{$lid} = {
                imp => { ALL => 'GEN' },
            };
            next;
        }
        my ($l, $r) = ($ldat->{l}, $ldat->{r});
        my $limps = $ldat->{imp} = {};
        while (my ($rid, $rdat) = each %{$ldat->{rna}}) {
            my $acc = $rdat->{acc} = 
                $self->cached_pkey_to_text( $rdat->{acc} );
            my $ns  = $acc =~ /^ENS/ ? 'ENST' : 
                $acc =~ /^[NX][MR]_/ ? 'RSR' : 'AR';
            foreach my $aid (keys %{$rdat->{aln}}) {
                my $imp = "UNK";
                if (my $hsp = $alns{$aid}) {
                    $imp  = "";
                    my $h = 0;
                    while ($h <= $#{$hsp}) {
                        my ($s, $e) = @{$hsp->[$h]};
                        # Stop if we are 'to the left of current HSP
                        last if ($r < $s);
                        $h++;
                        # Keep looking if we are still 'to the right'
                        next if ($l > $e);
                        # We overlap
                        if ($r - 2 < $s || $l + 2 > $e) {
                            # .. but only in the splice site
                            $imp = 'SPL';
                        } else {
                            # We should be in an exon
                            $imp = 'RNA';
                        }
                    }
                    unless ($imp) {
                        if ($h == 0 || $h > $#{$hsp}) {
                            # We were not within the alignment
                            $imp = 'GEN';
                        } else {
                            $imp = 'INT';
                        }
                    }
                }
                $rdat->{aln}{$aid} = $imp;
                $limps->{$ns}{$imp} = 1;
                $limps->{ALL}{$imp} = 1;
            }
        }
        while (my ($ns, $h) = each %{$limps}) {
            $limps->{$ns} = $self->heaviest_impact( keys %{$h} );
        }
    }
    $self->prebranch($rnaH);
    $self->bench_end();
    if (wantarray) {
        return ($rnaH, $rows, $ttab, $allRids);
    } else {
        $dbh->clear_temp_table( $ttab );
        return $rnaH;
    }
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub _bulk_rnas_for_flank_table {
    my $self = shift;
    $self->bench_start();
    my ($ptab, $howbad, $rnaH, $lvl2) = @_;
    my $defHB   = $self->default_howbad();
    my $flank   = $self->default_loc_to_rna_distance();
    my $bulkSrc = $self->bulk_cache_source( 'nearby_rnas' );
    my $ncid    = $self->_not_current_text_id();
    my $dbh     = $self->dbh;
    if (defined $howbad) {
        $bulkSrc = undef if ($howbad != $defHB);
    } else {
        $howbad = $defHB;
    }

    my ($cols);
    my ($colSel, $ttW, $ttAcc, $ttL, $ttR) = ("","", "?","?","?");
    if (ref($ptab)) {
        # Explicit values are being provided
    } else {
        # A temp table with the values has been passed.
        $cols    = [qw(obj_id acc_id ptl ptr)];
        $colSel  = ", " .join(', ', map { "tt.$_" } @{$cols});
        $ttW     = ", $ptab tt";
        ($ttAcc, $ttL, $ttR) = 
            ("tt.acc_id","tt.ptl - $flank","tt.ptr + $flank");

    }

    my $sql = <<SQL;
SELECT a.aln_id, r.acc_id, a.score, a.ptl1, a.ptr1, r.rna_id$colSel
  FROM alignment a, rna r$ttW
 WHERE a.seq1    = $ttAcc
   AND a.ptr1   >= $ttL
   AND a.ptl1   <= $ttR
   AND a.howbad <= $howbad
   AND a.seq2    = r.acc_id
   AND NOT EXISTS
      (SELECT nc.val_id FROM tagval nc 
        WHERE nc.obj_id = r.rna_id AND nc.tag_id = $ncid)
UNION
SELECT a.aln_id, r.acc_id, a.score, a.ptl2, a.ptr2, r.rna_id$colSel
  FROM alignment a, rna r$ttW
 WHERE a.seq2    = $ttAcc
   AND a.ptr2   >= $ttL
   AND a.ptl2   <= $ttR
   AND a.howbad <= $howbad
   AND a.seq1    = r.acc_id
   AND NOT EXISTS
      (SELECT nc.val_id FROM tagval nc 
        WHERE nc.obj_id = r.rna_id AND nc.tag_id = $ncid)

SQL

    my $sth  = $dbh->prepare
        ( -name => "Get RNAs overlapping bulk location IDs",
          -sql  => $sql );

    my $rows = [];
    if (ref($ptab)) {
        # Cycle through explicit values. Should be a 2D array of format
        # [ Object ID, GenomicAccessionID, Left, Right ]
        $self->bench_start('SingleEntryQueries');
        for my $p (0..$#{$ptab}) {
            my @xtra = @{$ptab->[$p]};
            my @binds = ($xtra[1], $xtra[2] - $flank, $xtra[3] + $flank);
            # warn $sth->show_and_explain([@binds, @binds]) unless ($p);
            $sth->execute(@binds, @binds);
            my $prows = $sth->fetchall_arrayref();
            push @{$rows}, map { [@{$_}, @xtra] } @{$prows};
        }
        $self->bench_end('SingleEntryQueries');
    } else {
        # Get it all in one go from the temp table
        # This is being painfully slow
        # $self->preprint( $sth->show_and_explain() );
        $self->bench_start('TempTableMerge');
        $sth->execute();
        $rows = $sth->fetchall_arrayref();
        $dbh->clear_temp_table( $ptab );
        $self->bench_end('TempTableMerge');
    }

    my (%alns, %rids, %textIds, %bulk);
    foreach my $row (@{$rows}) {
        my ($aid, $raccid, $sc, $genL, $genR, $rid,
            $oid, $gid, $l, $r) = @{$row};
        # $aid = aln_id in alignment
        # $rid = rna_id in rna
        # $oid = loc_id in location
        # $raccid = the txt_id for the RNA accession
        # $gid = The txt_id for the genomic accession
        # my ($oid, $l, $r, $gid, $rid, $raccid, $aid) = @{$row};
        my $targ = $rnaH->{$oid} ||= {
            l   => $l,
            r   => $r,
            gid => $gid,
            rna => {},
        };
        $rids{$rid}    = 1;
        map { $textIds{$_} = 1 } ($gid, $raccid);
        $alns{$aid}   ||= [ $gid, $genL, $genR ];
        $targ->{rna}{$rid}{aln}{$aid} = "";
        $targ->{rna}{$rid}{acc} = $raccid;
        push @{$bulk{$oid}}, $row;
    }
    my @allRids = keys %rids;
    # warn "RNAs: $#allRids\n";
    $self->bulk_tag_values( \@allRids, 'textToo' );
    $self->bulk_text_values( [keys %textIds] );
    if ($bulkSrc) {
        if ($bulkSrc) {
            while (my ($oid, $arr) = each %bulk) {
                $bulkSrc->{$oid} = $arr;
            }
        }
    }
    if ($lvl2) {
        $self->bench_end();
        return;
    }
    my @alnDat;
    foreach my $aid (sort keys %alns) {
        push @alnDat, [$aid, @{$alns{$aid}}];
    }
    # Now also find all RNAs that are near this RNA, and get boundary
    # information for them as well
    # warn "Alignments: $#alnDat\n";
    my $l2Data = \@alnDat;
    if ($cols) {
        $l2Data = $dbh->list_to_temp_table
            ( \@alnDat, ['integer','integer','integer','integer'], $cols);
        $dbh->_fast_add_index( $l2Data, [ 'obj_id' ] );
        $dbh->_fast_add_index( $l2Data, [ 'acc_id' ] );
    }
    $self->_bulk_rnas_for_flank_table( $l2Data, $howbad, {}, 'Level2' );

    $self->bench_end();
    return ($rows, \@allRids, \@alnDat);
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub bulk_accessions_for_object_ids {
    my $self = shift;
    my ($oids, $tagsToo) = @_;
    my $bulkSrc = $self->bulk_cache_source( 'accessions' );
    my @newids;
    foreach my $id (@{$oids}) {
        unless ($bulkSrc->{$id}) {
            push @newids, $id;
            $bulkSrc->{$id} = [];
        }
        return $bulkSrc if ($#newids == -1);
    }
    $self->bench_start();
    my $dbh   = $self->dbh;
    my $ttab  = $dbh->list_to_temp_table(\@newids, 'integer');
    #my $sth   = $dbh->prepare
    #    ( -name => "Get bulk accessions for object IDs",
    #      -sql  => "SELECT a.obj_id, a.acc_id, a.auth_id FROM accession a, $ttab tt WHERE a.obj_id = tt.col1", );
    my $sth   = $dbh->prepare
        ( -name => "Get bulk accessions for object IDs",
          -sql  => "SELECT a.obj_id, a.acc_id, a.auth_id FROM accession a ".
          "WHERE a.obj_id " . $self->_subquery_clause( 'col1', $ttab));
    # my $plan = $sth->explain(); die join("\n", map { $_->[0] } @{$plan});
    $sth->execute();
    my $rows = $sth->fetchall_arrayref();
    my (%tids);
    foreach my $row (@{$rows}) {
        my ($oid, $aid, $sid) = @{$row};
        push @{$bulkSrc->{$oid}}, [$aid, $sid];
        $tids{$aid} = $tids{$sid} = undef;
    }
    $self->bulk_text_values( [keys %tids] );
    $self->bench_end();
}
=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub bulk_features_for_object_ids {
    my $self = shift;
    my ($oids, $tagsToo) = @_;
    my $bulkSrc = $self->bulk_cache_source( 'features' );
    my @newids;
    foreach my $id (@{$oids}) {
        unless ($bulkSrc->{$id}) {
            push @newids, $id;
            $bulkSrc->{$id} = [];
        }
        return $bulkSrc if ($#newids == -1);
    }
    $self->bench_start();
    my $dbh     = $self->dbh;
    my $ttab  = $dbh->list_to_temp_table(\@newids, 'integer');
    #my $sth  = $dbh->prepare
    #    ( -name => "Get bulk feature data for object IDs",
    #      -sql  => "SELECT f.feat_id, f.src_id, f.pri_id, f.pos_to_left, f.pos_to_right, f.strand, f.rng_id FROM feature f, $ttab tt WHERE f.host_id = tt.col1", );
    my $sth  = $dbh->prepare
        ( -name => "Get bulk feature data for object IDs",
          -sql  => "SELECT f.feat_id, f.src_id, f.pri_id, f.pos_to_left, f.pos_to_right, f.strand, f.rng_id FROM feature f WHERE f.host_id ".
          $self->_subquery_clause( 'col1', $ttab));
    # my $plan = $sth->explain(); die join("\n", map { $_->[0] } @{$plan});
    $sth->execute();
    my $rows = $sth->fetchall_arrayref();
    my (%fids, %ridH);
    foreach my $row (@{$rows}) {
        my $fid = $row->[0];
        $fids{$fid} = undef;
        push @{$bulkSrc->{$fid}}, $row;
        if (my $rid = $row->[6]) {
            $ridH{$rid} = $row;
        } else {
            $row->[6] = [[ $row->[3], $row->[4]]];
        }
    }
    my @rids = keys %ridH;
    unless ($#rids == -1) {
        # We need to get range information, too
        my $rtab  = $dbh->list_to_temp_table(\@rids, 'integer');
        #my $rsth  = $dbh->prepare
        #    ( -name => "Get bulk range data for range IDs",
        #      -sql  => "SELECT r.rng_id, r.ptl, r.ptr FROM range r, $rtab tt WHERE rng_id = tt.col1", );
        my $rsth  = $dbh->prepare
            ( -name => "Get bulk range data for range IDs",
              -sql  => "SELECT r.rng_id, r.ptl, r.ptr FROM range r ".
              "WHERE r.rng_id ". $self->_subquery_clause( 'col1', $rtab));
        $rsth->execute();
        my $rows = $rsth->fetchall_arrayref();
        my %found;
        foreach my $row (@{$rows}) {
            my ($rid, $l, $r) = @{$row};
            push @{$found{$rid}}, [$l, $r];
        }
        while (my ($rid, $row) = each %ridH) {
            if (my $arr = $found{$rid}) {
                $row->[6] = [ sort { $a->[0] <=> $b->[0] } @{$arr} ];
            } else {
                $self->msg("[!!]", "Failed to locate range data for rng_id = $rid");
                $row->[6] = [[ $row->[3], $row->[4]]];
            }
        }
        $dbh->clear_temp_table( $rtab );
    }
    if ($tagsToo) {
        my @fs = keys %fids;
        $self->bulk_text_values( \@fs );
        $self->bulk_tag_values( \@fs, 'textToo' );
    }
    $dbh->clear_temp_table( $ttab );
    $self->bench_end();
    return $bulkSrc;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub bulk_cache_source {
    my $self = shift;
    my $key  = uc(shift || "");
    if (my $hash = shift) {
        $self->{BULK_SOURCE}{$key} = $hash;
    }
    return $self->{BULK_SOURCE}{$key} ||= {};
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub bulk_population_criteria {
    my $self = shift;
    my ($pids) = @_;
    my $pops  = $self->bulk_cache_source( 'popcat' );
    return $pops if (!$pids || $#{$pids} == -1);
    $self->bench_start();
    my @newpids;
    foreach my $pid (@{$pids}) {
        unless ($pops->{$pid}) {
            push @newpids, $pid;
            $pops->{$pid} = [];
        }
    }
    if ($#newpids == -1) {
        $self->bench_end();
        return $pops;
    }
    $pids = \@newpids;

    my (%text, %pars, @newPids, %catIds);
    my @stack = @{$pids};
    my $dbh   = $self->dbh;
    my $catID = $self->text_to_pkey('Category');
    while ($#stack != -1) {
        my %pids;
        foreach my $pid (@stack) {
            $pids{$pid} = 1 unless ($pops->{$pid}[0]);
        }
        my @newIds = keys %pids;
        last if ($#newIds == -1);
        push @newPids, @newIds;
        my $ttab = $dbh->list_to_temp_table(\@newIds, 'integer');
        @stack = ();
        #my $names = $dbh->prepare
        #    ( -name => "Get bulk names for population IDs",
        #      -sql  => "SELECT DISTINCT p.pop_id, p.name_id, p.par_id ".
        #      "FROM population p, $ttab tt WHERE p.pop_id = tt.col1", );
        my $names = $dbh->prepare
            ( -name => "Get bulk names for population IDs",
              -sql  => "SELECT DISTINCT p.pop_id, p.name_id, p.par_id ".
              "FROM population p WHERE p.pop_id ".
              $self->_subquery_clause( 'col1', $ttab));
        my $tags = $self->bulk_tag_values( $ttab, 'textToo' );
        $names->execute();
        my $nRows = $names->fetchall_arrayref();
        my @nids = map { $_->[1] } @{$nRows};
        $self->bulk_text_values( \@nids );
        # Set all the direct values
        foreach my $nrow (@{$nRows}) {
            my ($pid, $nid, $parid) = @{$nrow};
            my $name = $text{ $nid } ||= $self->cached_pkey_to_text( $nid );
            my $dirDat = $pops->{$pid}[0] ||= {};
            $dirDat->{name} = [$name];
            if ($parid) {
                $dirDat->{parid} = [$parid];
                push @stack, $parid;
            }
            if (my $tr = $tags->{$pid}) {
                foreach my $row (@{$tr}) {
                    if ($row->[0] == $catID) {
                        my $cid = $row->[1];
                        $catIds{$cid} = 1;
                        push @{$pops->{$pid}[0]{cat}}, $self->
                            cached_pkey_to_text( $cid );
                    }
                }
            }
        }
        $dbh->clear_temp_table( $ttab );
    }
    $self->bulk_tag_values( [keys %catIds], 'textToo' );
    # Recurse to populate full values into the {all} category
    my $recAll;
    $recAll = sub {
        my ($pid) = @_;
        my $pself = $pops->{$pid};
        my $all   = $pself->[1];
        unless ($all) {
            $all = $pself->[1] = {};
            my @via = ($pself->[0]);
            if (my $par = $pself->[0]{parid}) { 
                push @via, &{$recAll}( $par->[0] );
            }
            my %done;
            foreach my $hash (@via) {
                while (my ($k, $v) = each %{$hash}) {
                    map { push @{$all->{$k}}, $_ unless ($done{$k}{$_}++) } @{$v};
                }
            }
        }
        return $all;
    };
    foreach my $pid (@newPids) {
        &{$recAll}( $pid );
    }
    $self->bench_end();
    return $pops;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub bulk_text_values {
    my $self = shift;
    my $oids = shift;
    return if (!$oids || $#{$oids} == -1);
    $self->bench_start();
    my $cache = $self->{OBJ_CACHE}{TEXT} ||= {};
    my $dbh   = $self->dbh;
    my @keep;
    map { push @keep, $_ unless ($cache->{$_}) } @{$oids};
    my $ttab = $dbh->list_to_temp_table(\@keep, 'integer');
#    my $get = $dbh->prepare
#        ( -name => "Map pkey to text in bulk",
#          -sql  => "SELECT t.txt_id, t.txt FROM normtxt t, $ttab tt WHERE t.txt_id = tt.col1" );
    my $get = $dbh->prepare
        ( -name => "Map pkey to text in bulk",
          -sql  => "SELECT t.txt_id, t.txt FROM normtxt t WHERE t.txt_id ".
          $self->_subquery_clause( 'col1', $ttab));
    $get->execute();
    my $rows = $get->fetchall_arrayref();
    foreach my $row (@{$rows}) {
        $cache->{$row->[0]} ||= $row->[1];
    }
    $dbh->clear_temp_table( $ttab );
    $self->bench_end();
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub bulk_tag_values {
    my $self = shift;
    my $ttab = shift;
    my $textToo = shift;
    my $bcs  = $self->bulk_cache_source( 'tagval' );
    return $bcs unless ($ttab);
    $self->bench_start();
    my $dbh   = $self->dbh;
    my $clear;
    if (ref($ttab)) {
        my @keep;
        foreach my $id (@{$ttab}) {
            unless ($bcs->{$id}) {
                $bcs->{$id} = [];
                push @keep, $id;
            }
        }
        $ttab = $dbh->list_to_temp_table(\@keep, 'integer');
        $clear = $ttab;
    }
    #my $get = $dbh->prepare
    #    ( -name => "Get tag-values in bulk",
    #      -sql  => "SELECT t.obj_id, t.tag_id, t.val_id ".
    #      "FROM tagval t, $ttab tt WHERE t.obj_id = tt.col1" );
    my $get = $dbh->prepare
        ( -name => "Get tag-values in bulk",
          -sql  => "SELECT t.obj_id, t.tag_id, t.val_id ".
          "FROM tagval t WHERE t.obj_id ".
          $self->_subquery_clause( 'col1', $ttab));
    $get->execute();
    my $rows = $get->fetchall_arrayref();
    my %ids;
    foreach my $row (@{$rows}) {
        my ($id, $t, $v) = @{$row};
        push @{$bcs->{$id}}, [$t, $v];
        map { $ids{$_} = undef } ($t, $v);
    }
    $self->bulk_text_values( [keys %ids] ) if ($textToo);
    $dbh->clear_temp_table( $clear ) if ($clear);
    $self->bench_end();
    return $bcs;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub parenthetical_tag_parser {
    my $self = shift;
    my $args = $self->parseparams( @_ );
    my $qry  = $args->{QUERY} || "";
    my $cnt  = $args->{COUNT} || 1;
    $qry     =~ s/[\n\r]+/ /g;
    my (@sqlBits, @parens, $started, $needAndOr, @binds);
    my $rv = {
        request => $qry,
        tables  => {},
        binds   => [],
    };
    while ($qry) {
        if ($qry =~ /^\s+(.*)/) {
            # Whitespace
            $qry = $1;
            next;
        } elsif ($qry =~ /^([\&\|]+)(.*)/) {
            # AND / OR
            $qry = $2;
            my $ao = $1;
            if ($ao =~ /^\&+$/) {
                push @sqlBits, 'AND';
            } elsif ($ao =~ /^\|+$/) {
                push @sqlBits, 'OR';
            } else {
                $rv->{error} = "Could not interpret AND/OR component '$ao'";
                last;
            }
            $needAndOr = 0;
        } elsif ($qry =~ /^([\(\)]+)(.*)/) {
            # Parentheses
            $qry = $2;
            foreach my $par (split('', $1)) {
                if ($needAndOr && $par eq '(') {
                    push @sqlBits, 'AND';
                    $needAndOr = 0;
                }
                push @sqlBits, $par;
                push @parens,  $par;
            }
        } elsif ($qry =~ /^\"([^\"]+)\"(.*)/ ||
                 $qry =~ /^\'([^\']+)\'(.*)/ ||
                 $qry =~ /^([^\(\)\&\|]+)(.*)/) {
            # Quoted value, or non-whitespace value
            $qry = $2;
            my $tv = $1;
            push @sqlBits, 'AND' if ($needAndOr);
            my $tag = $tv;
            my $val;
            my $test = '=';
            if ($tv =~ /^(.*?)\s*!=+\s*(.+?)$/) {
                ($tag, $val) = ($1, $2);
                $test = '!=';
            } elsif ($tv =~ /^(.*?)\s*=+\s*(.+?)$/) {
                ($tag, $val) = ($1, $2);
            }
            my @tvbits;
            my $tvt = "tvp" . $cnt;
            if ($tag) {
                if (my $pk = $self->text_to_pkey($tag, "noCreate")) {
                    push @tvbits, "$tvt.tag_id = $pk";
                }
            }
            if ($val) {
                if (my $pk = $self->text_to_pkey($val, "noCreate")) {
                    push @tvbits, "$tvt.val_id $test $pk";
                }
            }
            if ($#tvbits == -1) {
                if ($test eq '!=') {
                    push @sqlBits, '1 = 1';
                } else {
                    push @sqlBits, '1 = 0';
                }
            } else {
                my $tvs = join(' AND ', @tvbits);
                $tvs = "( $tvs )" if ($#tvbits != 0);
                push @sqlBits, $tvs;
                $rv->{tables}{$tvt} = "tagval";
                $cnt++;
            }
            $needAndOr = 1;
        } else {
            $rv->{error} = "Parse failure at '$qry'";
            last;
        }
    }
    return $rv if ($rv->{error});
    my $chkParen = join('', @parens);
    while ($chkParen =~ /\(\)/) { $chkParen =~ s/\(\)//; }
    if ($chkParen) {
        $rv->{error} = "Parentheses are mismatched: '".join('', @parens)."'";
        return $rv;
    }
    $rv->{sql} = join(' ', @sqlBits);
    return $rv;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub acc_id_to_loc_ids {
    my $self = shift;
    my $sth  = $self->{STH}{ACCID_TO_LOC} ||= $self->dbh->prepare
        ( -name => "Get location for accession ID",
          -sql  => "SELECT l.loc_id, a.auth_id FROM location l, accession a ".
          "WHERE l.loc_id = a.obj_id AND a.acc_id = ?");
    $sth->execute( shift);
    my $rows = $sth->fetchall_arrayref();
    if (wantarray) {
        my %u = map { $_->[0] => 1 } @{$rows};
        return keys %u;
    } else {
        return $rows;
    }
    my @rv;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub txt_id_to_rna_ids {
    my $self = shift;
    my $id   = shift;
    # First, does the ID directly map to an RNA?
    my $get = $self->{STH}{RNA_DIRECT_FROM_ACC} ||= $self->dbh->prepare
        ( -name => "Get RNA with explicit ID",
          -sql  => "SELECT r.rna_id FROM rna r ".
          "WHERE r.acc_id = ?");
    my @rv = $get->get_array_for_field( $id );
    # If we got a hit, take it and run:
    return @rv unless ($#rv == -1);
    # Perchance this is an unversioned accession?
    my $txt  = $self->pkey_to_text( $id );
    return @rv if (CORE::length($txt) <= 6); # Eh, too short
    my $vids = $self->versioned_rna_ids($txt);
    my %uacc = map { $_ => 1 } map { @{$_} } values %{$vids};
    push @rv, keys %uacc;
    return @rv;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub versioned_rna_ids {
    my $self = shift;
    my $uacc = shift; # Unversioned accession
    return wantarray ? () : {} unless ($uacc);
    $self->bench_start();
    my %rv;
    my @accids;
    if ($uacc =~ /\.(\d{1,3})$/) {
        # This is actually a versioned ID, or at least looks like one
        # We will treat it as such. Note that there are some speceies with
        # internal ".###" that will fail here if provided an unversioned ID
        my $vers = $1;
        my $id   = $self->text_to_pkey($uacc);
        push @accids, [$id, $vers];
    } else {
        # We need to use the unversioned ID in a leading wildcard search
        # to see if the DB is aware of any versioned ones

        # I tried this with a joined statement:
        # SELECT r.rna_id, nt.txt FROM rna r, normtxt nt WHERE 
        #       r.acc_id = nt.txt_id AND upper(substr(nt.txt, 1, 20)) LIKE $qry
        # ... but it took 230 ms to run. Performing the accession
        # and RID lookup independently is 100x faster

        # Placeholders are NOT used because query planning is very poor
        # without knowledge of the string being searched.

        my $uca  = uc($uacc);
        my $ulen = CORE::length($uca);
        my $dbh  = $self->dbh();
        my @queries;
        if (1 || $ulen > 16) {
            # We will not be able to do a surgical wildcard query because
            # the finite substr(20) index is too short for the unversioned
            # accession plus the dot plus up to three digits
            my $qry = substr($uca, 0, 20);
            $qry    =~ s/_/\\_/g;
            $qry   .= '.' if ($ulen < 20);
            $qry   .= '%';
            push @queries, $dbh->quote($qry);
            
        } else {
            # The unversioned accession is short enough that we can use
            # single-character wildcards to check for discrete matches
            # Each search takes about 2.5 ms

            # These have the advantage that in certain cases they could
            # constrain returned matches if there were many DB entries
            # such as 'NM_001234.SAMPLE-72913'

            # However, they also introduce three times as many searches
            # so I am overriding this block.
            
            $uca =~ s/_/\\_/g;
            for my $d (1..3) {
                push @queries, $dbh->quote($uca. '.' . ('_' x $d)) ;
            }
        }
    
        foreach my $qry (@queries) {
            # $self->bench_start('Get ACCID');
            # Roughly 2.5 msec
            my $sth = $dbh->prepare
                ( -name => "Get RNA with unversioned accession",
                  -sql  => "SELECT nt.txt_id, nt.txt FROM normtxt nt ".
                  "WHERE upper(substr(nt.txt, 1, 20)) LIKE $qry");
            $sth->execute();
            my $rows = $sth->fetchall_arrayref();
            foreach my $row (@{$rows}) {
                # Need to validate each accession
                my ($id, $vacc) = @{$row};
                if ($vacc =~ /^\Q$uacc\E\.(\d+)$/i) {
                    # The recovered string does, in fact, look like a
                    # versioned accession of the query
                    push @accids, [$id, $1];
                }
            }
            # $self->bench_end('Get ACCID');
        }
    }
    #self->bench_start('Get RID');
    # Roughly 630 usec
    my $acc2rid = $self->{STH}{GET_RID_FOR_VACC} ||= $self->dbh->prepare
        ( -name => "Get RNA ID via versioned accession",
          -sql  => "SELECT r.rna_id FROM rna r WHERE r.acc_id = ?");
    foreach my $adat (@accids) {
        my ($id, $vers) = @{$adat};
        foreach my $rid ($acc2rid->get_array_for_field( $id )) {
            push @{$rv{$vers}}, $rid;
        }
    }
    #$self->bench_end('Get RID');
    $self->bench_end();
    return wantarray ? (sort { $b <=> $a } keys %rv) : \%rv;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub txt_id_to_gene_ids {
    my $self = shift;
    my $id   = shift;
    # First, is the ID already a gene ID? Note: Genes are just
    # text::Text objects that have some methods added to them
    my $get = $self->{STH}{GENE_DIRECT_FROM_ACC} ||= $self->dbh->prepare
        ( -name => "Get Gene with explicit ID",
          -sql  => "SELECT tv1.obj_id FROM tagval tv1 WHERE tv1.obj_id = ? AND tv1.tag_id = ".$self->text_to_pkey("Object Type")." AND tv1.val_id = ". $self->text_to_pkey("Gene") );
    my @rv = $get->get_array_for_field( $id );
    # If we got a hit, take it and run:
    return @rv unless ($#rv == -1);

    my $getAli = $self->{STH}{GENE_FROM_ALIASES} ||= $self->dbh->prepare
        ( -name => "Get Gene from alias",
          -sql  => "SELECT tv1.obj_id FROM tagval tv1, tagval tv2 WHERE tv2.tag_id = ? and tv2.val_id = ? AND tv1.obj_id = tv2.obj_id AND tv1.tag_id = ".$self->text_to_pkey("Object Type")." AND tv1.val_id = ". $self->text_to_pkey("Gene") );
    my $sets = $self->{GENE_LOOKUP_SETS} ||=
        [ map { [ map { $self->text_to_pkey($_) } @{$_} ] } (['Official Symbol'],['Other Symbols']) ];

    # Whatever it was, it was not a simple accession
    for my $s (0..$#{$sets}) {
        # This loop is designed to exit as soon as a match is found
        # So if a hit is made against the official symbol, it will NOT
        # see if hits exist against unofficial symbols
        my @found;
        foreach my $tagid (@{$sets->[$s]}) {
            push @found, $getAli->get_array_for_field( $tagid, $id );
            # warn $getAli->pretty_print( $tagid, $id );
        }
        push @rv, @found;
        last unless ($#found == -1);
    }

    return @rv;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub txt_id_to_population_ids {
    my $self = shift;
    my $id   = shift;
    # First check primary name assigned to populations
    my $get = $self->{STH}{POP_ID_FOR_NAME_TXT_ID} ||= $self->dbh->prepare
        ( -name => "Get Population with name ID",
          -sql  => "SELECT pop.pop_id FROM population pop ".
          "WHERE pop.name_id = ?" );
    my @rv = $get->get_array_for_field( $id );
    # If we got a hit, take it and run:
    return @rv unless ($#rv == -1);

    # dbSNP also has names stored under the Handle tag
    # This is done because not all RefSnp handles are unique
    # Any population also has categories defined
    foreach my $tag ('Handle', 'Category') {
        my $getAli = $self->{STH}{"POP_FROM_TAG_".$tag} ||= $self->dbh->prepare
            ( -name => "Get population from handle",
              -sql  => "SELECT pop.pop_id FROM population pop, tagval tv1 ".
              "WHERE tv1.tag_id = ".$self->text_to_pkey($tag).
              " and tv1.val_id = ? AND tv1.obj_id = pop.pop_id" );
        push @rv, $getAli->get_array_for_field( $id );
    }
    # Anything else I should check for here?
    return @rv;
}

=head2 add_hsps

 Title   : add_hsps
 Usage   : my @grouped = $ml->add_hsps( $hspList1, $hspList2 )
 Function: Simplify list of HSPs, add two HSP lists, and/or expand HSPs
           by extra distance
 Returns : 2D array reference of hsp coordinates
 Args    : [0] A 2D array reference
           [1] Optional second 2D array reference

This method does simple aggregation of HSP coordinates. Each HSP
should be represented by:

     [0] = left coord
     [1] = right coord
     [2] = optional left extra flank
     [3] = optional right extra flank
     [4] = optional array of notes

Each HSP coordinate can have an independent flank value. This was
implemented to allow terminal coordinates (those fully 5' and 3' on a
set of HSPs) to have different flanks than the internal
(intron-facing) coordinates. Additionally, the fourth index can hold
an array of notes (used for tracking distinct query terms on
individual HSPs).

If a second HSP list is provided, it will simply be combined with the
first. The HSPs are sorted, and merged such that the final list covers
all the provided HSPs, and is internally non-overlapping. So:

   [100,120] [150,200] [160,400]

will become:

   [100,120] [150,400]

=cut

sub add_hsps {
    my $self = shift;
    # Consider one or two inputs
    # Flanks are tracked as they are needed for subtractive operations
    my @input = ( @{$_[0] || []}, @{$_[1] || []} );
    map { $_->[2] ||= 0;
          $_->[3] ||= 0;
          $_->[4] ||= []; } @input;

    # If the third argument is defined, then it will be presumed that
    # it is a flank value and all coordinates should be extended by
    # that value AND the values stored in [2] and [3] for each
    # HSP. Pass a third parameter of zero to only use the HSP-specific
    # flank values.
    if (defined $_[2]) {
        my $expandFlank = $_[2];
        foreach my $hsp (@input) {
            $hsp->[0] -= ($expandFlank + $hsp->[2]);
            $hsp->[1] += ($expandFlank + $hsp->[3]);
            # Reset the stored flanks to zero
            $hsp->[2] = 0;
            $hsp->[3] = 0;
        }
    }
    @input = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @input;
    my @rv = shift @input;
    foreach my $hspR (@input) {
        my $hspL = $rv[-1];
        if ($hspL->[1] <= $hspR->[0]) {
            # The right HSP is non-overlapping, add it as a new HSP
            push @rv, $hspR;
            next;
        }
        # The HSPs overlap
        # Extend the right side if needed
        $hspL->[1] = $hspR->[1] if ($hspL->[1] < $hspR->[1]);
        # Extend the extra flanks if needed
        $hspL->[2] = $hspL->[0] - ($hspR->[0] - $hspR->[2]) if
            (($hspL->[0] - $hspL->[2]) > ($hspR->[0] - $hspR->[2]));

        $hspL->[3] = ($hspR->[1] + $hspR->[3]) - $hspL->[1] if
            (($hspL->[1] + $hspL->[3]) < ($hspR->[1] + $hspR->[3]));

        # Combine notes
        push @{$hspL->[4]}, @{$hspR->[4]};
    }
    return wantarray ? (\@rv) : \@rv;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub subtract_hsps {
    my $self = shift;
    # Consider two inputs
    # First normalize the 'keep' and 'toss' arrays:
    my $keep = $self->add_hsps( $_[0], undef, $_[3] );
    my $toss = $self->add_hsps( $_[1] );
    # We need to make sure that the keep array is sorted including the flank:
    $keep = [sort {($a->[0] - $a->[2]) <=> 
                       ($b->[0] - $b->[2]) ||
                       ($a->[1] + $a->[3]) <=> 
                       ($b->[1] + $b->[3])} @{$keep}];
    #$self->prebranch({ keep => \@keep, toss => \@toss });
    my $tind = 0;
    my (@rv, @discard);
    my $iloop = 0;
    while ($#{$keep} != -1) {
        if ($#{$toss} < $tind) {
            # No more subtractive HSPs, keep everything
            push @rv, @{$keep};
            last;
        }
        # Boundaries of the current toss HSP
        my ($tl, $tr) = @{$toss->[$tind]};
        # Boundaries and flanks of current keep HSP:
        my ($kl, $kr, $fl, $fr, $notes) = @{$keep->[0]};
        # Flank-adjusted coordinates:
        my ($klf, $krf) = ($kl - $fl, $kr + $fr);
        # $self->msg('[>]', "($kl, $kr, $fl, $fr [($klf, $krf)]) vs ($tl, $tr)");
        if ($klf >= $tr - 1) {
            # The toss HSP is to the left and completely non-overlapping
            # We can ignore the toss HSP
            # $self->msg('[^]', "TOSS COMPLETE - ($tl, $tr) Dist = ".($klf - $tr + 1));
            $tind++;
            next;
        }
        my $kHsp = shift @{$keep};
        if ($krf <= $tl + 1) {
            # The toss HSP is to the right and completely non-overlapping
            # We can keep the full HSP
            # $self->msg('[+]', "SURVIVE - ".join('+', @{$kHsp}));
            push @rv, $kHsp;
            next;
        }
        # There is some overlap. Find out what we need to trim
        if ($kl >= $tl && $kr <= $tr) {
            push @discard, [$kl, $kr, undef, undef, $notes];
            # HSP is completely lost
            # $self->msg('[-]',"FULL LOSS - "..join('+', map { ref($_) ? '['.join(',',@{$_}).']' : $_ } @{$kHsp}));
            next;
        } elsif ($kl < $tl && $kr > $tr) {
            # The HSP is split in two by the toss HSP
            # We will keep the truncated left HSP, discarding right flank:
            push @rv, [$kl, $tl + 1, $fl, 0];
            push @discard, [$tl, $tr, undef, undef, $notes];
            # Add add a new truncated right HSP back to the stack
            my $newHsp = [$tr - 1, $kr, 0, $fr ];
            # $self->msg('[?]', "OVERLAP - HSP Split");
            if ($#{$keep} == -1 || ($keep->[0][0] - $keep->[0][2] >
                                    $newHsp->[0] - $newHsp->[2])) {
                # We can just put the HSP at the front of the stack
                unshift @{$keep}, $newHsp;
            } else {
                # The HSP is no longer the left-most one!
                unshift @{$keep}, $newHsp;
                $keep = [sort {($a->[0] - $a->[2]) <=> ($b->[0] - $b->[2]) ||
                                   ($a->[1] + $a->[3]) <=> ($b->[1] + $b->[3])} @{$keep}];

                # $self->msg("[!!]", "RE-ORDER SPLIT HSP");
            }
            next;
        }
        
        if ($kl >= $tl) {
            # We need to shave off some of the left side
            my $shaved = $tr - $klf - 1; # == $tr - 1 - ($klf + 1) + 1;
            # $self->msg('[x]', "OVERLAP - HSP Left Trim [$shaved of $fl]");
            # First shave off the flank
            $fl -= $shaved;
            if ($fl < 0) {
                # The entire flank is gone
                # Also need to shave off the actual HSP
                push @discard, [$kHsp->[0], undef, undef, undef, $notes];
                $kHsp->[0] -= $fl;
                $discard[-1][1] = $kHsp->[0] + 1;
                $fl = 0;
            }
            $kHsp->[2] = $fl;
        }
        if ($kr <= $tr) {
            # We need to shave off right side
            my $shaved = $krf - $tl - 1; # == $krf - 1 - ($tl + 1) + 1;
            # $self->msg('[x]', "OVERLAP - HSP Right Trim [$shaved off $fr]");
            $fr -= $shaved;
            if ($fr < 0) {
                # Also need to shave off the actual HSP
                push @discard, [undef, $kHsp->[1], undef, undef, $notes];
                $kHsp->[1] += $fr;
                $discard[-1][0] = $kHsp->[1] - 1;
                $fr = 0;
            }
            $kHsp->[3] = $fr;
        }
        # $self->msg("[#]", "TRIMMED - ".join('+', map { ref($_) ? '['.join(',',@{$_}).']' : $_ } @{$kHsp}));
        if ($kHsp->[1] - $kHsp->[0] > 1) {
            # We have an HSP with some information
            if ($#{$keep} == -1 || ($keep->[0][0] - $keep->[0][2] >
                                    $kHsp->[0] - $kHsp->[2])) {
                # We can just put the HSP at the front of the stack
                # warn "$keep->[0][0] - $keep->[0][2] > $kHsp->[0] - $kHsp->[2]" if ($#{$keep} != -1);
                unshift @{$keep}, $kHsp;
            } else {
                # $self->msg("[!!]", "RE-ORDER TRIMMED HSP");
                # The HSP is no longer the left-most one!
                unshift @{$keep}, $kHsp;
                $keep = [sort {($a->[0] - $a->[2]) <=> ($b->[0] - $b->[2]) ||
                                   ($a->[1] + $a->[3]) <=> ($b->[1] + $b->[3])} @{$keep}];

            }
        }
        # What may get ignored:
        # If a toss HSP overlaps ONLY the padded-flank portion of a keep HSP
        # then the "extra bit" of un-masked flank will not be considered
        # if (++$iloop > 1000);
    }
    return wantarray ? (\@rv, \@discard) : \@rv;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub preprint {
    my $self = shift;
    my $txt = join("\n", map { $self->esc_xml($_) } map { defined $_ ? $_ : "" } @_);
    return "" unless ($txt);
    print "<pre style='margin:3px; background-color:#eee;'>$txt</pre>";
    return "";
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub prebranch {
    my $self = shift;
    print "<pre style='border: solid blue 1px; margin:3px; background-color:#eee;'>".$self->branch(@_)."</pre>";
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub cx_query_part {
    my $self = shift;
    my $args = $self->parseparams( @_ );
    my $hsps = $args->{HSPS};
    my $type = $args->{TYPE};
    my $fill;
    if (my $c = $args->{COLOR}) {
        $fill = $c;
    } elsif (my $f = $args->{FLAG}) {
        $type ||= 
            $f eq  '1' ? 'Requested region' : 
            $f eq '-1' ? 'Excluded region' : 
            $f eq '-2' ? 'Requested, but overlapped with excluded region' :
            "Unknown Region '$f'";
        $fill = 
            $f eq  '1' ? '#33ff33' :  # Requested region
            $f eq '-1' ? '#ff0000' :  # Excluded region
            $f eq '-2' ? '#888888' :  # Requested, but overlapped excluded
            '#cc99ff';                # Default fallback color
    } else {
        $fill = '#ffcc99';
    }
    
    
    
    my $mainFill  = "rgb(".join(",", $self->rgb_values( $fill )).")";
    my $flankFill = "rgba(".join(",", $self->rgb_values( $fill ), 0.3).")";
    my @parts;
    for my $h (0..$#{$hsps}) {
        my $hsp = 
        my ($l, $r, $lf, $rf, $notes) = @{$hsps->[$h]};
        my %un = map { $_ => 1 } @{$notes || []};
        my $note = join(', ', sort keys %un);
        next unless (defined $l);
        my %tags;
        $tags{'Defined By'} = $notes if ($notes);
        $tags{'Region Type'} = [$type] if ($type);
        $tags{'Region Size'} = [($r - $l -1)." bp"];
        push @parts, {
            hideName => 1,
            fill     => $mainFill,
            data     => [[$l + 1, $r -1]],
            tags     => { %tags },
            name     => $note,
        };
        $tags{'Extended Region'} = ["This region represents extra search space requested around your primary query"];
        if ($lf) {
            $tags{'Region Size'} = ["$lf bp"];
            push @parts, {
                hideName => 1,
                fill     => $flankFill,
                data     => [[$l - $lf + 1, $l]],
                name     => $note,
                tags     => { %tags },
                height   => 1,
            };
        }
        if ($rf) {
            my $e = $r + $rf - 1;
            $tags{'Region Size'} = ["$rf bp"];
            if (0 && $h < $#{$hsps}) {
                # Do we overlap with any part of the next HSP?
                my $nxt = $hsps->[$h + 1];
                my $nl  = $nxt->[0] - ($nxt->[2] || 0) + 1;
                my $overlap = $e - $nl + 1;
                if ($overlap > 0) {
                    my $nrf = $rf - $overlap;
                    if ($nrf < 1) {
                        $e = 0;
                    } else {
                        $e = $r + $nrf - 1;
                        $tags{'Shared Region'} = ["Shown region is $overlap bp shorter to simplify shared area with upstream query region"];
                    }
                }
            }
            push @parts, {
                hideName => 1,
                fill     => $flankFill,
                data     => [[$r, $e]],
                name     => $note,
                tags     => {%tags},
                height   => 1,
            } if ($e);
        }
    }
    return wantarray ? @parts : \@parts;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub rna_ids_to_accessions {
    my $self = shift;
    my @rv;
    my $sth = $self->{STH}{"RNAID_TO_ACC"} ||= $self->dbh->prepare
        ( -name => "Convert rna_id to accession",
          -sql  => "SELECT nt.txt FROM normtxt nt, rna r ".
          "WHERE nt.txt_id = r.acc_id AND r.rna_id = ?" );
    foreach my $id (@_) {
        push @rv, $sth->get_single_value( $id ) || "";
    }
    return @rv;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub rna_ids_to_cached_accessions {
    my $self = shift;
    # For many applications we expect the same RNA IDs to occur repeatedly
    my $sth = $self->{STH}{"RNAID_TO_ACC"} ||= $self->dbh->prepare
        ( -name => "Convert rna_id to accession",
          -sql  => "SELECT nt.txt FROM normtxt nt, rna r ".
          "WHERE nt.txt_id = r.acc_id AND r.rna_id = ?" );
    my $cache = $self->{OBJ_CACHE}{RNA_ID_ACC} ||= {};
    my %rv;
    foreach my $id (@_) {
        my $acc = $cache->{$id};
        unless (defined $acc) {
            $acc = $cache->{$id} = $sth->get_single_value( $id ) || "";
            # warn "$id = $acc\n";$self->preprint($sth->show_and_explain([$id]) );
        }
        $rv{$id} = $acc if ($acc);
    }
    return wantarray ? values %rv : \%rv;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub rna_accession_to_gene_accessions {
    my $self = shift;
    my $racc = shift;
    my @rv;
    return @rv unless ($racc);
    $self->bench_start();
    my $sth = $self->{STH}{GENE_NAMES_FROM_RNA_NAME} ||=
        $self->dbh->prepare
        ( -name => "Get gene names for RNA names",
          -sql => "SELECT nt1.txt FROM normtxt nt1, tagval tv WHERE ".
          "tv.tag_id = ".$self->cached_text_to_pkey('Has RNA').
          " AND nt1.txt_id = tv.obj_id AND tv.val_id = ?");
    my $accid = $self->text_to_pkey($racc);
    @rv = $sth->get_array_for_field( $accid );
    if ($#rv == -1 && $racc =~ /(.+)\.\d{1,3}$/) {
        # Try the unversioned accession
        my $accid = $self->text_to_pkey($1);
        @rv = $sth->get_array_for_field( $accid );
    }
    $self->bench_end();
    return @rv;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub web_support_dir {
    my $self = shift;
    if (my $nv = shift) {
        $self->{WEB_SUPPORT_DIR} = $nv;
    } elsif (!$self->{WEB_SUPPORT_DIR}) {
        my $dir = $self->module_path( $self );
        $dir    =~ s/\/[^\/]+$//;
        $dir   .= "/MapLoc/webSupportFiles";
        $self->{WEB_SUPPORT_DIR} = $dir;
    }
    return $self->{WEB_SUPPORT_DIR} || ".";
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub web_support_url {
    my $self = shift;
    if (my $nv = shift) {
        $self->{WEB_SUPPORT_URL} = $nv;
    } elsif (!$self->{WEB_SUPPORT_URL}) {
        $self->{WEB_SUPPORT_URL} = "webSupportFiles";
    }
    return $self->{WEB_SUPPORT_URL} || ".";
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub cx_support_links {
    my $self = shift;
    my $args = $self->parseparams( @_ );
    my $wsUrl = $self->web_support_url();
    my $xsd  = $args->{DIR} || $self->web_support_dir();
    my ($headHtml, $bodyHtml) = ("","");
    return wantarray ? ($headHtml, $bodyHtml) : $headHtml unless (-d $xsd);
    my $dynamicStyles = "$xsd/dynamicStyles.css";
    if ($self->file_is_older_than_module( $dynamicStyles, $self)) {
        # Need to make the CSS file
        if (open(CSS, ">$dynamicStyles")) {
            print CSS "/* This file is automatically generated and may be rewritten without notice! */\n";
            print CSS $self->css_styles();
            close CSS;
            chmod( 0666, $dynamicStyles);
        } else {
            $bodyHtml .= "<span style='color:red'>Failed to write stylesheet '"
                .$self->esc_xml($dynamicStyles)."' : $!</span><br />\n";
        }
    }
    
    foreach my $ojs (sort $self->read_dir( -dir  => $xsd,
                                           -recurse => 1,
                                           -keep => '.(js|css)$')) {
        if ($ojs =~ /^\Q$xsd\E\/(.+)$/) {
            # Found an additional script
            my $xs = $1;
            if ($xs =~ /^[a-z0-9_\/]+\.(js|css)$/i) {
                my $type = $1;
                if ($type eq 'js') {
                    $headHtml .= "  <script src='$wsUrl/$xs' ".
                        "type='text/javascript'></script>\n";
                } elsif ($type eq 'css') {
                    $headHtml .= "  <link href='$wsUrl/$xs' ".
                        "type='text/css' rel='stylesheet' />\n";
                }
            } else {
                $bodyHtml .= "<span style='color:red'>Ignoring extra script '".$self->esc_xml($xs)."' - files can only contain letters, numbers and underscores</span><br />\n";
                next;
            }
        }
    }
    return wantarray ? ($headHtml, $bodyHtml) : $headHtml;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub defer_alleles_to_db {
    my $self = shift;
    my $rows = $self->{DEFERRED_ALLELES} ||= [];
    push @{$rows}, @_;
    return $#{$rows};
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub deferred_allele_count {
    return $#{shift->{DEFERRED_ALLELES} || []} + 1;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub add_deferred_alleles {
    my $self = shift;
    my $rows = $self->{DEFERRED_ALLELES} || [];
    if ($#{$rows} != -1) {
        # For 8202 rows in 27 batches this averaged 52.8ms = 170us per row
        $self->bench_start();
        my $ia = $self->dbh->insert_array( 'allele', $rows);
        if (my $err = $DBI::errstr) {
            if ($err =~ /^ERROR:  duplicate key value/) {
                # Expected Postgres error code for duplicate row
                # We could try to work out which rows are old, but
                # will likely be most robust (and probably as fast)
                # to just fallback to row-by-row addition
                $self->_add_allele_rows( $rows );
            } else {
                $self->death("Unexpected error while COPYing population alleles to database", $err);
            }
        }
        $self->bench_end();
    }
    delete $self->{DEFERRED_ALLELES};
}

*_add_allele_rows = \&_add_allele_rows_function;
=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub _add_allele_rows_function {
    my $self = shift;
    $self->bench_start();
    my $rows = shift || [];
    my $sth = $self->{STH}{"CREATE_ALLELE_BY_FUNCTION"} ||= $self->dbh->prepare
        ( -name => "Create an allele row with function call",
          -sql  => "SELECT add_allele_to_db(?,?,?,?,?,?)",
          -ignore => "duplicate key value" );

    foreach my $binds (@{$rows}) {
        #$self->bench_start('execute');
        my $rv = $sth->execute(@{$binds}, 0);
        #$self->bench_end('execute');
    }
    $self->bench_end();
}


=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub _unknown_pop_id {
    my $self = shift;
    return $self->{UNK_POP_ID} ||= $self->get_population
        ("Undefined Population", 'Undefined')->pkey();
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub _temporary_pop_id {
    my $self = shift;
    return $self->{TEMP_POP_ID} ||= $self->get_population
        ("Temporary Data", 'Temporary')->pkey();
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub _allele_update_status {
    # This is mainly for debugging
    my $self = shift;
    my $nice = {
        0 => "Defer for COPY",
        1 => "Update NULL + NULL",
        2 => "Update Freq + NULL",
        3 => "Update Freq + Count",
        4 => "Insert new row",
        -1 => "ERROR",
    };
    my $txt = "";
    my $vals = $self->{ALLELE_UPDATE_STATUS} || {};
    foreach my $k ( sort {$vals->{$b} <=> $vals->{$a}} keys %{$vals}) {
        $txt .= sprintf(" %5d %s\n", $vals->{$k}, $nice->{$k} || $k);
    }
    return $txt;
}

=head2 

 Title   : 
 Usage   : 
 Function: 
 Returns : 
 Args    : 

=cut

sub _subquery_clause {
    my $self = shift;
    my ($col, $table) = @_;
    return sprintf("= ANY ( (SELECT array(SELECT %s FROM %s))::integer[])",
                   $col, $table);

    # I suffered very bad performance for subqueries. It appears that
    # this is a known problem with postgres (8.1.11). The ANY
    # solution above was found on StackOverflow at:

    # http://stackoverflow.com/questions/14987321/postgresql-in-operator-with-subquery-poor-performance
    
    # Queries that performed poorly:

    # return sprintf("= %s.%s", $table, $col);
    # return sprintf("IN (SELECT %s FROM %s)", $col, $table);

    # EXPLAINing these queries showed that Postgres was running a
    # seqscan on the main query, even if the subplan was returning
    # only a handful of primary key IDs for the main query. The query
    # planner does not seem capable of considering the contents of the
    # subquery.

    # My presumption is that either the ANY or the typecast forces
    # Postgres to collect the values of the subquery in advance, and
    # it is then smart enough to inspect them before planning the rest
    # of the query.

}



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
1;
