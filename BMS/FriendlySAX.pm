
=head1 NAME

 BMS::FriendlySAX.pm - Easy-to-use interface to SAX XML parsing

=head1 DESCRIPTION

A wrapper for XML::Parser::PerlSAX. Provides a very simple module to
extract data from an XML file quickly and with minimal memory
footprint, while maintaining the object structure of the file.

=head1 SYNOPSIS

    use BMS::FriendlySAX;

    # Define the path to the XML file you wish to parse:
    my $file = "/tmp/coolData.xml";
    # Will automatically handle gz files if your system has gunzip:
    $file = "/tmp/lotsOfCoolData.xml.gz";

    # Define the XML unit you wish to capture by tag name. For each
    # such tag, the returned value will include all data for that tag,
    # and all children, recovered recursively. To get <author> nodes
    # and their children:

    my $tags = "author";

    # You can also get more than one type of tag:

    $tags = [ 'author', 'publisher' ];

    # You are free to ignore any tags that do not interest you. If
    # these tags tend to contain a lot of data, you can choose to
    # explicitly skip them. This will slightly speed up parsing, and
    # will minimize the memory footprint:

    my $skip = undef;   # to skip nothing
    $skip = 'fulltext'; # skip one tag
    $skip = [ 'fulltext', 'comments' ]; # 2 or more tags to skip

    # If you do not wish to parse the whole document, you can set a
    # limit; parsing will halt after that many target tags have been
    # returned.

    my $limit = 10;

    # You want to do something with the data, right?
    # Defined a method that will accept each parsed node and analyze it

    my $meth = \&my_cool_method;

    eval {
        my $fs = BMS::FriendlySAX->new
            ( -file    => $file,
              -tag     => $tags,
              -limit   => $limit,
              -skip    => $skip,
              -method  => $meth,  );
    };

    # Note that the eval block is needed if you wish to continue your
    # program flow after reaching a $limit.

    # Every node is a hash reference organized as follows:
    my $example_node = {
        NAME => "XML_NODE_NAME", # eg <table> node name is "table"
        TEXT => "String with text (character) content of node",
        ATTR => {  }, # Hash reference of attributes, attr_name => attr_value
        KIDS => [  ], # Array of all child nodes, in order of appearance
        BYTAG => {
            # Hash reference keyed to the XML tag names of the children
            # Each key points to an array of child nodes
        },
        PARENT => $parentNode, # weakened reference to parent node
    };

    # Child nodes are structured just like above, and are nested under
    # both KIDS *and* BYTAG

    # Null values you may encounter:

    # NAME should always be a string (hopefully a valid XML tag name)
    # TEXT will never be undef, but will be "" (empty string) for no text
    #   Leading and trailing white space will be stripped by default
    # ATTR will *always* return a hash ref
    #   As per any hash, you can request values for any key; you will get
    #   a string if the key was a populated attribute, otherwise undef
    # KIDS will *always* return an array ref, which will be empty for no kids
    # BYTAG will *always* return a hash ref
    #   Again, you can call for any key. If the key corresponds to the tag
    #   name of a child of the node, you will get an array reference,
    #   otherwise you will get undef
    # The root node will have undef for PARENT (it is root, after all)
    
    # Example XML:

    # <cookie type='tasty'>
    #    <source>Bob</source>
    #    <source>MainStreetGrocery</source>
    #    <flavor>lemon</flavor>
    # </cookie>

    sub my_cool_method {
        # This assumes that you set -tags => 'cookie' ...
        my ($node) = @_;
        my $name = $node->{NAME};          # Should be "cookie"
        my $type = $node->{ATTR}{'type'};  # Should be "tasty"
        my $foo  = $node->{ATTR}{'animal'} # Should be undef (no such attr)
        # Do something with *every* child node:
        foreach my $child ( @{$node->{KIDS}} ) {
            # Should cycle through 3 nodes: source, source, flavor
            # You would probably call another method here, and pass $child
        }

        # BYTAG for any given tag will either be undef or an array ref:
        my $sources = $node->{BYTAG}{'source'};
        if ($sources) {
            # Always should test to make sure undef was not returned
            my @list;
            foreach my $source ( @{$sources} ) {
                # Each $source will be another node (<source>)
                push @list, $source->{TEXT};
            }
            # You could also have done the above as:
            @list = map { $_->{TEXT} } @{$sources};
            
            # so @list == ('Bob','MainStreetGrocery');
        }
        # This is a convienence method to allow you to recursively visualize
        # a node as text output:
        print BMS::FriendlySAX::node_to_text( $node );
    }

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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::FriendlySAX;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

my $VERSION = 
    ' $Id: FriendlySAX.pm,v 1.17 2012/04/10 13:50:54 tilfordc Exp $ ';


use strict;
use XML::Parser::PerlSAX;
use Scalar::Util qw(weaken);
use BMS::ErrorInterceptor;
use BMS::Utilities::Escape;

use vars qw(@ISA);
@ISA   = qw(BMS::ErrorInterceptor BMS::Utilities::Escape);

=head1 Methods

=head2 new

 Title   : new
 Usage   : my $obj = BMS::FriendlySAX->new(@arguments)

           Note that the new call should be within an eval {}; block
           if you are planning to use -limit *and* continue on with
           the flow of the calling program. This is because the only
           way I have found to force expat to stop parsing is to
           die().

 Returns : A blessed BMS::FriendlySAX object - which is not terribly
           useful, because instantiating the object also causes all
           parsing to occur.

 Args    : Associative array of arguments. Recognized keys [Default]:

     -file Required. Path to the XML file. Note that the file can be
           gzipped (.gz) - FriendlySAX will transparently gunzip the
           file for parsing (assuming your system has gunzip and the
           shell can find it).

      -tag Required. One or more tag names. The 'tag name' is the name
           of the XML node. For example, the node <dog>poodle</dog>
           has a tag name of 'dog'.

           In most cases you will be getting data for one tag only;
           you can then pass a simple string. You can also pass an
           array reference of strings if you wish to recover more than
           one tag type.

   -method Required. A subroutine reference. Every time a completed
           tag is parsed, the resulting hash structure will be passed
           to this subroutine. See the Synopsis above for more
           information on how you might structure your method to parse
           returned nodes.

     -skip Optional. Similar to -tag, used to specify one or more tag
           names that you wish to *ignore*. Such tags will be absent
           from the returned hash structure. This is useful for
           limiting memory footprint, and will slightly speed up
           parsing.

  -skipall Default 0. If true, the ONLY the tags specified by -tag
           will be captured. Note that if you use this option that you
           MUST specify any potential intermediate tag names with
           -tag, otherwise the deeper nested children will be ignored.

    -limit Optional limit. If undef or zero, then the whole XML
           document will be parsed. If you pass an integer, then only
           that number of tags will be returned, then parsing will be
           halted. See note on eval under Usage above...

  -nolimit Optional array of tag names. Tags listed here will not
           count toward a -limit. Useful if a document has dictionary
           nodes that you want to read in total but still apply the
           limit to the 'main' tags. Note that if the limit is reached
           for other tags, then parsing will still halt (so it is only
           useful if the non-limitted tags occur before the limit is
           reached).

    -final In order to break out of the SAX parser when the
           user-specified limit is reached, I had to die(). There is
           probably a better interupt to send, but rather than spend
           the time researching that, I am including this option; pass
           a method reference, and the method will be called before
           die.

 -textmeth Optional subroutine reference. The subroutine should accept
           an array reference of strings as an argument and return a
           single string scalar. It will be applied to the TEXT
           component of a node when the end of an element is
           reached. The purpose is to allow post-processing of the
           character data; this is particularly useful for escaping
           character content eg &apos;, \n etc.

           Note also that if a node contains both text and child
           nodes, no information on how the two relate will be
           preserved. The child nodes will be accessible in the order
           they occur in, but the text will be mashed into a single
           string. So the following node:

           <exclaim>Boy, <b>that</b> hurt!</exclaim>

           ... will have a single <b> child, and a single text entry
           "Boy,  hurt!".

=cut

our $parser;

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self = { };    
    bless ($self, $class);
    my $args = $self->parseparams( -file  => "",
				   -error => undef,
				   -bench => undef,
                                   -verbose => 0,
                                   -skipall => 0,
				   @_);


    foreach my $key ('FILE', 'VERBOSE', 'TAG', 'METHOD', 'TEXTMETH',
                     'QUIETLIMIT') {
        $self->{$key} = $args->{$key};
    }

    my $file = $args->{FILE} || "";
    my @bits = split(/\//, $file);
    my $shortfile = $bits[-1] || $file;
    unless (-e $file) {
        $file ||= '-UNDEFINED-';
	$self->death("I could not find input -file '$file'");
    }

    my $method = $self->{METHOD} || $self->{CALLBACK};
    if (!$method) {
	$self->death("You must define -method (your XML tag handler)");
    } elsif (ref($method) ne 'CODE') {
	$self->death("-method (your XML tag handler) is not a code reference");
    }

    if ($self->{TEXTMETH} && ref($self->{TEXTMETH}) ne 'CODE') {
	$self->death("-textmeth (your text post-processing handler) ".
                     "is not a code reference");
    }

    $self->err("Initiating parse of '$shortfile'") if ($self->{VERBOSE});

    my $handler = GenericHandler->new( %{$args} );
    my $parser  = XML::Parser::PerlSAX->new( Handler => $handler );

    if ($file =~ /\.gz$/) {
	open(GZIP, "gunzip -c $file|") || 
	    $self->death("Failed to establish gunzip pipe for '$file'");
        $self->ignore_error('gunzip: stdout: Broken pipe');
        $self->{OPENFH} = *GZIP;
	$parser->parse( Source => { ByteStream => *GZIP } );
    } else {
	$parser->parse( Source => { SystemId => $file } );
    }
    $self->err("Finished parse of '$shortfile'") if ($self->{VERBOSE});
    if (my $fcb = $self->{FINAL}) {
        &{$fcb}( );
    }
    return $self;
}

sub DESTROY {
    my $self = shift;
    my $ofh  = $self->{OPENFH};
    if ($ofh) {
        close $ofh;
    }
}

=head2 node_to_text

 Title   : node_to_text
 Usage   : print BMS::FriendlySAX::node_to_text( $node );
 Returns : A string, formatted to represent the hierarchical strucutre of
           the passed node.
 Args    : A single FriendlySAX XML node.

This method is primarily useful for debugging

=cut

sub node_to_text {
    my ($hash, $level, $opts) = @_;
    $level ||= 0;
    $opts    = lc($opts || "");
    my $name = $hash->{NAME};
    my $pad  = "   " x $level;
    my @kids    = @{$hash->{KIDS}};
    my $text    = $hash->{TEXT};
    my $hasText = defined $text && $text ne "" ? 1 : 0;
    my @tags    = sort keys %{$hash->{ATTR}};
    if ($#kids == -1 && $#tags == -1 && !$hasText && $opts !~ /showempty/) {
        # This tag has no children, no attributes and no text
        # And the user has not indicated that they want empty tags
        return "";
    }
    # Open the XML tag
    my $string = "$pad$name";
    # Add attributes, if any, to the tag
    if ($#tags > -1) {
        my @tvs = map {sprintf("%s='%s'", $_, $hash->{ATTR}{$_})} @tags;
        if ($opts =~ /expand/) {
            map { $string .= "\n$pad  ATTR $_" } @tvs;
        } else {
            $string .= " (".join(" ", @tvs).")";
        }
    }
    if ($hasText) {
        $text = substr($text, 0, 70) . '...' if (length($text) > 70);
        $string .= "\n$pad  TEXT" if ($opts =~ /expand/);
        $string .= " : $text";
    }
    $string .= "\n";
    foreach my $kid (@kids) {
        # Recursively print children:
        $string .= &node_to_text( $kid, $level + 1);
    }
    return $string;
}

sub node_to_xml {
    my ($hash, $level, $opts) = @_;
    $level ||= 0;
    my $name = $hash->{NAME};
    my $pad  = "   " x $level;
    my @kids = @{$hash->{KIDS}};
    my $text = $hash->{TEXT};
    # Open the XML tag
    my $xml = "$pad<$name";
    # Add attributes, if any, to the tag
    my $tn = 0;
    foreach my $tag (sort keys %{$hash->{ATTR}}) {
        $xml .= sprintf(" %s='%s'", $tag, &esc_xml_attr
                        ($hash->{ATTR}{$tag}));
        $tn++;
    }
    if ((!defined $text || $text eq '') && $#kids == -1) {
        # This tag only has no text or children
        return "" unless ($tn); # Return nothing if no tags, either
        # Can put in option to allow empty nodes to be preserved...
        return "$xml />\n";
    } elsif ($#kids == -1) {
        # Text-only tag
        return "$xml>".&esc_xml($text)."</$name>\n";
    }
    # At least one child
    $xml .= ">\n";
    $xml .= "  $pad".&esc_xml($text)."\n" if ($text);
    foreach my $kid (@kids) {
        # Recursively print children:
        $xml .= &node_to_xml( $kid, $level + 1);
    }
    $xml .= "$pad</$name>\n";
    return $xml;
}

sub esc_xml { BMS::Utilities::Escape::esc_xml( undef, @_ ) }
sub esc_xml_attr { BMS::Utilities::Escape::esc_xml_attr( undef, @_ ) }

1;

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package GenericHandler;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

use strict;
use BMS::Utilities;
use Scalar::Util qw(weaken);

use vars qw(@ISA);
@ISA      = qw(BMS::Utilities);

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self = {
        ACTIVE => 0,
        COUNT  => 0,
        SKIP   => {},
        TAGS   => {},
	STACK  => [ {
            NAME   => "",
            KIDS   => [],
            PARENT => undef,
            _NOTE_ => 'Null root element',
        } ],
    };
    bless ($self, $class);
    my $args = $self->parseparams
        ( -textmeth => \&_default_text_meth,
          @_);
    $self->{LIMITMSG} = $args->{LIMITMSG};
    if (my $tagkey = $args->{TAG} || $args->{TAGS}) {
        # We are only getting a particular tag
        my @tags   = ref($tagkey) ? @{$tagkey} : split(/\s+/, $tagkey);
        if ($#tags < 0) {
            $self->death("You must define -tag (the XML tagname for the unit ".
                         "you wish to process)");
        }
        $self->{TAGS} = { map { $_ => 1 } @tags };
    } else {
        # We are getting all tags
        $self->{GETALL} = 1;
    }
    my @noLimTags;
    if (my $nolim = $args->{NOLIMIT}) {
        # Some tags should not be counted towards a limit
        @noLimTags   = ref($nolim) ? @{$nolim} : split(/\s+/, $nolim);
    }
    $self->{NOLIMIT} = { map { $_ => 1 } @noLimTags };
    
    if (my $skipkey = $args->{SKIP}) {
        # The user wants to skip information for specific nodes
        my @skip = ref($skipkey) ? @{$skipkey} : split(/\s+/, $skipkey);
        $self->{SKIP} = { map { $_ => 1 } @skip };
    }

    foreach my $key ('METHOD', 'LIMIT', 'QUIETLIMIT',
                     'SKIPALL', 'FINAL', 'TEXTMETH') {
	$self->{$key} = $args->{$key};
    }
    return $self;
}


sub start_element {
    my ($self, $element) = @_;
    my $name   = $element->{Name};
    my $parent = $self->{STACK}[-1];
    my $node   = {
        NAME   => $name,
        ATTR   => $element->{Attributes},
        TEXT   => [],
        KIDS   => [],
        BYTAG  => {},
    };
    weaken($node->{PARENT} = $parent);

    # At a minimum, maintain a thin stack of the full node hierarchy
    # with attributes. Members on the stack that are not a target (or
    # do not have a target parent) can ONLY be counted on to have
    # attribute data associated with them. Child nodes will be
    # stripped from these entries when a target hits end_element().

    push @{$self->{STACK}}, $node;

    return unless ($self->{ACTIVE} || $self->{TAGS}{$name} || $self->{GETALL});
    # Do nothing further unless:
    # 1. We are already processing a requested tag
    # 2. This node is a specifically requested tag
    # 3. We are getting ALL tags
    unless ($self->{ACTIVE}) {
        $self->{ACTIVE} = 1;
        # Specify this node as a requested target
        $node->{IS_TARGET} = 1;
    }
    if ($self->{SKIPALL}) {
        # We are skipping all nodes - continue only if this tag is
        # explicitly set to be kept
        return unless ($self->{TAGS}{$name});
    } elsif ($self->{SKIP}{$name}) {
        # Request to ignore this tag
        return;
    }
    # Add this node to as a child to the previous parent node
    push @{$parent->{KIDS}}, $node;
    # Also add it to the BYTAG referenced structure:
    push @{$parent->{BYTAG}{$name}}, $node;
}

sub end_element {
    my ($self, $element) = @_;
    # Remove the completed node from the stack
    my $node   = pop @{$self->{STACK}};

    # Do nothing further if we are not actively parsing a target
    return unless ($self->{ACTIVE});

    if ($node->{TEXT}) {
        $node->{TEXT} = &{$self->{TEXTMETH}}( $node->{TEXT} );
    } else {
        $node->{TEXT} = '';
    }
    
    if ($node->{IS_TARGET}) {
        # This was a requested target node, we should act on it:
        my $method = $self->{METHOD};
        &{$method}( $node );
        # Check parents to see if any are themselves target nodes:
        my $active = 0; my @parents;
        my $focus  = $node->{PARENT};
        while ($focus) {
            push @parents, $focus;
            $active++ if ($focus->{IS_TARGET});
            $focus = $focus->{PARENT};
        }
        unless ($active) {
            # None of the parent nodes were targets
            # Flag that we are no longer collecting data:
            $self->{ACTIVE} = 0;
            foreach my $parent (@parents) {
                # Clear all parent structures (prevent memory bloat)
                $parent->{KIDS}  = [];
                $parent->{BYTAG} = {};
            }
        }
        # Keep track of records analyzed, in case we are limiting analysis
        my ($name, $lim) = ($node->{NAME}, $self->{LIMIT});
        $self->{COUNT}++ unless ($self->{NOLIMIT}{$name});
        if ($lim && $self->{COUNT} >= $lim) {
            # We have analyzed the user requested parse limit
            $! = 0;
            if (my $fcb = $self->{FINAL}) {
                &{$fcb}( $node );
            }
            if ($self->{QUIETLIMIT}) {
                die "User limit ($lim) halts processing on $name node\n";
            } elsif (my $msg = $self->{LIMITMSG}) {
                die $msg;
            } else {
                die "\n";
            }
        }
    }
}

sub characters {
    my ($self, $chars) = @_;
    # Do nothing if we are not actively parsing a target
    return unless ($self->{ACTIVE});
    my $parent = $self->{STACK}[-1];
    my $txt    = $chars->{Data};
    push @{$parent->{TEXT}}, $txt;
}

sub _default_text_meth {
    my ($arr) = @_;
    # Remove leading and trailing whitespace:
    map { s/\s+$//; s/^\s+//; } @{$arr};
    # Join WITHOUT spaces
    return join('', @{$arr});
}

1;
