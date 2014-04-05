# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::FriendlyGraph;
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

$BMS::FriendlyGraph::VERSION = 
    ' $Id: FriendlyGraph.pm,v 1.14 2013/12/13 14:45:39 tilfordc Exp $ ';

use strict;
use BMS::ArgumentParser;
use Math::Trig;

use vars qw(@ISA);
@ISA   = qw(BMS::FriendlyGraph::Common);

our $idCounter = 0;
our (%attrTypes, %attrCallBacks, %attrCBoverWrite, %xmlCallBacks);

our $cxShapeMap = {
    sphere   => 'sphere',
    circle   => 'sphere',
    square   => 'square',
    triangle => 'triangle',
    star     => 'star',
    rhombus  => 'rhombus',
    diamond  => 'rhombus',
    octagon  => 'octagon',
    oval     => 'oval',
    plus     => 'plus',
    minus    => 'minus',
    pacman   => 'pacman',
    mdavid   => 'mdavid',
};

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self = {
        NODES => {},
        EDGES => {},
        CLASSES => {},
        GNAME => 'Universe',
        LIVE_WRITE => 0,
        LIVE_READ  => 0,
        DB_CURRENT => undef,
    };
    bless ($self, $class);
    return $self;
}

sub type  { return 'Graph' }
sub graph { return shift }

sub name {
    my $self = shift;
    if (my $nv = $_[0]) {
        unless ($self->{GNAME} eq $nv) {
            # We have changed the name of the database
            # We need to clear all 'current' flags
            map { delete $_->{DB_CURRENT} } $self->all_objects();
            delete $self->{DBID};
        }
        $self->{GNAME} = $nv;
    }
    return $self->{GNAME};
}

*dbi = \&database;
sub database {
    my $self = shift;
    unless ($self->{DB}) {
        my $ap   = BMS::ArgumentParser->new( -nocgi => 1 );
        my @defs = $ap->module_parameters( "BMS::FriendlyGraph" );
        my $args = $self->parseparams( @defs, -graph => $self, @_ );
        eval {
            $self->{DB} = BMS::FriendlyGraph::DBH->new( $args );
        };
    }
    return $self->{DB};
}

sub dbh {
    my $self = shift;
    if (my $dbi = $self->dbi) {
        # This used to set $dbi to $self->{DBI}
        # That failed to automatically connect to the database if it had
        # not been done already. I have changed this to force connection
        # by calling the method, but I can't recall why I had not done
        # so in the first place; my change may cause problems. 11 Jan 2011
        return $dbi->dbh;
    } else {
        return undef;
    }
}

sub write {
    my $self = shift;
    $self->write_attributes();
}

sub read {
    my $self = shift;
    $self->read_attributes( );
}

sub read_all {
    my $self = shift;
    map { $_->read() } $self->all_objects;
}

sub write_all {
    my $self = shift;
    map { $_->write() } $self->all_objects;
}

sub live_read {
    my $self = shift;
    if (defined $_[0]) {
        $self->{LIVE_READ} = $_[0] ? 1 : 0;
    }
    return $self->{LIVE_READ};
}

sub live_write {
    my $self = shift;
    if (defined $_[0]) {
        $self->{LIVE_WRITE} = $_[0] ? 1 : 0;
    }
    return $self->{LIVE_WRITE};
}

sub live_database {
    my $self = shift;
    my $lr   = $self->live_read(@_);
    my $lw   = $self->live_write(@_);
    if ($lr && $lw) {
        return 1;
    } elsif (!$lr && !$lw) {
        return 0;
    } else {
        # Partial live
        return -1;
    }
}

sub all_objects {
    my $self = shift;
    return ($self, $self->each_node, $self->each_edge );
}

sub node_sort_attribute {
    my $self = shift;
    if (defined $_[0]) {
        $self->{NODE_SORT_ATTR} = $_[0];
    }
    return $self->{NODE_SORT_ATTR} || "";
}

*all_nodes = \&each_node;
sub each_node {
    my $self = shift;
    my @rv = values %{$self->{NODES}};
    if (my $nsa = $self->node_sort_attribute()) {
        @rv = sort { ($a->attribute($nsa) || "") 
                         cmp ($b->attribute($nsa) || "") } @rv;
    }
    # sort { $a->name cmp $b->name }
    return $self->filter_nodes( \@rv, @_ );
}

sub node {
    my $self = shift;
    my $name = shift || '';
    # Accept node objects as input. We will degrade them to name, in case they
    # are being brought in from a different graph object:
    if (my $r = ref($name)) {
        if ($r =~ /BMS::FriendlyGraph::Node/) {
            $name = $name->name;
        } else {
            $self->death("Not sure how to handle node() request via reference",
                         $name);
        }
    }
    $name =~ s/^\s+//; $name =~ s/\s+$//;
    return unless ($name);
    unless ($self->{NODES}{uc($name)}) {
        my $node = BMS::FriendlyGraph::Node->new($name, $self);
        $self->{NODES}{uc($name)} = $node;
        $self->{NODE_BY_ID}{ $node->numeric_id() } = $node;
    }
    return $self->{NODES}{uc($name)};
}

sub node_by_numeric_id {
    my $self = shift;
    my $id   = shift;
    $id = '?' unless (defined $id);
    return exists $self->{NODE_BY_ID}{$id} ? $self->{NODE_BY_ID}{$id} : undef;
}

sub has_node {
    my $self = shift;
    my $name = shift || '';
    # Accept node objects as input. We will degrade them to name, in case they
    # are being brought in from a different graph object:
    $name = $name->name if (ref($name));
    $name =~ s/^\s+//; $name =~ s/\s+$//;
    return 0 unless ($name);
    return exists $self->{NODES}{uc($name)} ? 1 : 0;
}

*delete_node = \&remove_node;
sub remove_node {
    my $self = shift;
    foreach my $req (@_) {
        my $node = $req;
        my $name = uc(ref($node) ? $node->name : $node);
        next unless (exists $self->{NODES}{$name} && $self->{NODES}{$name});

        # WE PROBABLY SHOULD RESET {graph} TO UNDEF FOR THE NODE
        # It is possible that the node is being held by another structure
        # and may still be utilized later. However, currently a lot of
        # methods will attempt to call ->graph off the node; we should
        # first put in safety checking or fallback behavior.


        if ($node = $self->{NODES}{$name}) {
            foreach my $edge ($node->each_edge) {
                $self->remove_edge( $edge );
            }
        }
        delete $self->{NODES}{$name};
    }
}

sub remove_isolated_nodes {
    my $self = shift;
    my @removed;
    foreach my $node ($self->each_node) {
        next if ($node->edge_count);
        push @removed, $node;
        $self->remove_node($node);
    }
    return @removed;
}

sub node_count {
    my $self = shift;
    my @nodes = $self->each_node();
    return $#nodes + 1;
}

sub edge_count {
    my $self = shift;
    my @edges = $self->each_edge();
    return $#edges + 1;
}

sub set_node_attributes {
    my $self = shift;
    my $node = $self->node( shift );
    return $node->set_attributes( @_ ) if ($node);
}

sub each_edge {
    my $self = shift;
    my @rv;
    foreach my $hash (values %{$self->{EDGES}}) {
        push @rv, values %{$hash};
    }
    return $self->filter_edges( \@rv, @_ );
}

sub edge {
    my $self = shift;
    return undef if ($#_ < 1);
    my $node1 = $self->node(shift);
    my $node2 = $self->node(shift);
    return undef unless ($node1 && $node2);
    my $isDir = shift;
    $isDir    = $self->attribute('directed') unless (defined $isDir);
    if (!$isDir && $node1->uckey() gt $node2->uckey()) {
        ($node1, $node2) = ($node2, $node1);
    }
    my ($k1, $k2) = map { $_->uckey() } ($node1, $node2);
    
    unless ($self->{EDGES}{$k1}{$k2}) {
        $self->{EDGES}{$k1}{$k2} = BMS::FriendlyGraph::Edge->new
            ($node1, $node2, $self, $isDir);
    }
    return $self->{EDGES}{$k1}{$k2};
}

sub remove_edge {
    my $self = shift;
    my ($req1, $req2, $isDir) = @_;
    return undef unless ($req1);
    my ($name1, $name2) = (uc($req1));
    if (my $r = ref($req1)) {
        if ($r =~ /FriendlyGraph::Node/) {
            $name1 = $req1->uckey;
        } elsif ($r =~ /FriendlyGraph::Edge/) {
            ($name1, $name2) = map { $_->uckey } $req1->nodes();
            $isDir = $req1->directed;
        } else {
            $self->err("Unable to remove edge for unknown object", $r);
            return;
        }
    }
    return undef unless 
        (exists $self->{NODES}{$name1} && $self->{NODES}{$name1});
    unless ($name2) {
        return undef unless ($req2);
        if (my $r = ref($req2)) {
            if ($r =~ /FriendlyGraph::Node/) {
                $name1 = $req1->uckey;
            } else {
                $self->err("Unable to remove edge for unknown second object", $r);
                return;
            }
        } else {
            $name2 = uc($req2);
        }
    }
    return undef 
        unless (exists $self->{NODES}{$name2} && $self->{NODES}{$name2});
    my ($node1, $node2) = map { $self->{NODES}{$_} } ($name1, $name2);
    if (!$isDir && $node1->uckey() gt $node2->uckey()) {
        ($node1, $node2) = ($node2, $node1);
    }
    my ($k1, $k2) = map { $_->uckey() } ($node1, $node2);
    my $edge = $self->{EDGES}{$k1}{$k2};
    delete $self->{EDGES}{$k1}{$k2};
    map { $_->unlink($edge) } ($node1, $node2);
    return $edge;
}

sub class {
    my $self = shift;
    my ($name, $nv) = @_;
    my $rv;
    if ($name) {
        $name = uc($name);
        if (defined $nv) {
            if (!$nv) {
                # Clear the class
                delete $self->{CLASSES}{$name};
            } elsif (my $r = ref($nv)) {
                if ($r eq 'HASH') {
                    $self->{CLASSES}{$name} = $nv;
                } else {
                    die "No logic defined to set classes by $r values";
                }
            } else {
                die "No logic defined yet to define classes on scalar values";
            }
        }
        $rv = $self->{CLASSES}{$name};
    }
    return $rv || {};
}

sub rdf_time {
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) =
        localtime(time);
    my $rdfTime = sprintf("%04d-%02d-%02d %2d:%02d:%02d", $year+1900,$mon+1,
                          $mday , $hour,$min,$sec);
    return $rdfTime;
}

sub guess_attribute_types {
    my $self = shift;
    my %types;
    my %possible;
    foreach my $obj ($self->all_objects) {
        foreach my $key ($obj->each_attribute) {
            next if ($types{$key});
            my $val = $obj->attribute($key);
            next unless (defined $val);
            $types{$key} = &BMS::FriendlyGraph::DBH::is_num($val) ?
                undef : 'string';
        }
    }
    map { $types{$_} ||= 'real' } keys %types;
    return \%types;
}

sub to_graphviz {
    my $self = shift;
    my $label   = $self->attribute('label') || "Friendly Graph Output";
    $label     .= " ".&rdf_time() if ($self->attribute('stampLabel'));
    
    my $gv = sprintf("%sgraph %s {\n", 
                     $self->attribute('directed') ? 'di' : '',
                     $self->esc_gv($label));
    $gv .= "  // GraphViz specification - for more information see:\n";
    $gv .= "  // http://www.research.att.com/sw/tools/graphviz/\n";
    $gv .= "  // Built using BMS::FriendlyGraph\n\n";

    # Set up global default values for all types
    my @globals;
    foreach my $type qw(Graph Node Edge) {
        my %attHash;
        foreach my $agdat ($self->each_attribute_generator($type)) {
            my ($key, $val, $ow) = @{$agdat};
            next if (!defined $val || ref($val));
            $attHash{$key} = $self->esc_gv($val);
        }
        $self->remap_attributes( \%attHash, 'graphviz', $type );
        my @attrs = map { $_.'='.$attHash{$_} } sort keys %attHash;
        push @globals, sprintf("  %s [%s]\n", lc($type),join(',', @attrs))
            unless ($#attrs == -1);
    }
    unless ($#globals == -1) {
        $gv .= "  // Global variables\n";
        $gv .= join('', @globals);
        $gv .= "\n";
    }

    my @nodes = $self->each_node();
    if ($#nodes == -1 ) {
        $gv .= "  // Empty Graph, no nodes\n";
    } else {
        $gv .= sprintf("  // Total of %d Node%s\n", $#nodes + 1,
                        $#nodes == 0 ? '' : 's');
        foreach my $node (@nodes) {
            $gv .= $node->to_graphviz();
        }
    }
    $gv .= "\n";

    my @edges = $self->each_edge();
    if ($#edges == -1 ) {
        $gv .= "  // No edges defined\n";
    } elsif ($#nodes == -1) {
        $gv .= "  // ERROR - edges are defined, but nodes are not\n";
    } else {
        $gv .= sprintf("  // Total of %d Edge%s\n", $#edges + 1,
                        $#edges == 0 ? '' : 's');
        foreach my $edge (@edges) {
            $gv .= $edge->to_graphviz();
        }
    }
    $gv .= "\n}\n";
    return $gv;
}

sub to_png {
    my $self = shift;
    my $args = $self->parseparams( @_ );
    my $cmd  = $self->attribute('graphviz');
    unless ($cmd) {
        $self->err("You must provide the path to the graphviz executable in the 'graphviz' attribute");
        return undef;
    }
    my $png = $args->{FILE} || $args->{PNG} ||
        "/tmp/FriendlyGraph-$$-".time.".png";
    my $base = $png;
    $base    =~ s/\.png$//i;
    my $map  = $base.".imap";
    my $gv   = $base.".gv";
    if (open(GVF, ">$gv")) {
        print GVF $self->to_graphviz();
        close GVF;
    } else {
        $self->err("Failed to create graphviz data file", $gv, $!);
        return undef;
    }
    my $cmdLine = sprintf("%s -Timap -o %s -Tpng -o %s %s",
                          $cmd, $map, $png, $gv);
    system($cmdLine);
    unless (-s $png) {
        $self->err("Possible failure generating png file",
                   $png, "Via $gv");
    }
    return wantarray ? ($png, $map, $gv) : $png;
}

sub to_xgmml {
    my $self = shift;

    my $rdfTime = &rdf_time();
    my $id      = $self->id;
    my $label   = $self->attribute('label') || "Friendly Graph Output";
    $label     .= " $rdfTime" if ($self->attribute('stampLabel'));
    my $isDir   = $self->attribute('directed') ? 1 : 0;
    my $desc    = $self->attribute('description') || 'N/A';
    my $type    = $self->attribute('type') || 'N/A';
    my $url     = $self->attribute('url')  || 'http://bioinformatics.bms.com/';

    my $xml =  <<EOF;
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
 <graph label="$label" directed="$isDir" id="$id"
    xmlns:dc="http://purl.org/dc/elements/1.1/" 
 xmlns:xlink="http://www.w3.org/1999/xlink" 
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" 
    xmlns:cy="http://www.cytoscape.org" 
       xmlns="http://www.cs.rpi.edu/XGMML">

 <att name="networkMetadata">
  <rdf:RDF>
   <rdf:Description rdf:about="http://www.cytoscape.org/">
    <dc:date>$rdfTime</dc:date>
    <dc:type>$type</dc:type>
    <dc:description>$desc</dc:description>
    <dc:identifier>N/A</dc:identifier>
    <dc:title>$url</dc:title>
    <dc:source>http://www.cytoscape.org/</dc:source>
    <dc:format>Cytoscape-XGMML</dc:format>
   </rdf:Description>
  </rdf:RDF>
 </att>

EOF

    $xml .= $self->attribute_xgmml(1);
    $xml .= "\n";

    my @nodes = $self->each_node();
    if ($#nodes == -1 ) {
        $xml .= "  <!-- Empty Graph, no nodes -->\n";
    } else {
        $xml .= sprintf("  <!-- Total of %d Node%s -->\n", $#nodes + 1,
                        $#nodes == 0 ? '' : 's');
        foreach my $node (@nodes) {
            $xml .= $node->to_xgmml(1);
        }
    }
    $xml .= "\n";

    my @edges = $self->each_edge();
    if ($#edges == -1 ) {
        $xml .= "  <!-- No edges defined -->\n";
    } elsif ($#nodes == -1) {
        $xml .= "  <!-- ERROR - edges are defined, but nodes are not -->\n";
    } else {
        $xml .= sprintf("  <!-- Total of %d Edge%s -->\n", $#edges + 1,
                        $#edges == 0 ? '' : 's');
        foreach my $edge (@edges) {
            $xml .= $edge->to_xgmml(1);
        }
    }
    $xml .= "\n";

    $xml .= "</graph>\n";
    return $xml;
}

sub to_gml {
    my $self = shift;

    my $rv = "graph [\n";
    $rv   .= $self->attribute_gml(1);

    foreach my $node ($self->each_node()) {
        $rv .= $node->to_gml(1);
    }
    foreach my $edge ($self->each_edge()) {
        $rv .= $edge->to_gml(1);
    }
    $rv .= "]\n";
    return $rv;
}

*to_canvasXpress   = \&to_canvasXpress;
*to_canvasExpress  = \&to_canvasXpress;
*to_canvasexpress  = \&to_canvasXpress;
*to_canvas_express = \&to_canvasXpress;
*to_canvasxpress   = \&to_canvasXpress;
*to_canvas_xpress  = \&to_canvasXpress;
sub to_canvasXpress {
    my $self = shift;
    my $args = $self->parseparams( @_ );
    my $domId = $args->{ID} || $args->{DOMID};
    $domId = $domId ? "'$domId'" : 'null';
    my $attr = $self->all_attribute_values('format=canvasxpress');
    my $expt = $args->{EXPORT};
    my $bigZ = 10 ** 9;
    my @nodes = map { $_->to_canvasXpress( $expt ) } sort {
        ($a->attribute('zIndex') || $bigZ) <=>
            ($b->attribute('zIndex') || $bigZ) ||
            $a->id cmp $b->id
    } $self->each_node;
    my @edges = map { $_->to_canvasXpress( $expt ) } $self->each_edge;
    my $data = {
        nodes => \@nodes,
        edges => \@edges,
        # legend => { nodes => [], edges => [], text => [] },
    };
    if (my $lnks = $attr->{links}) {
        $data->{links} = $lnks;
        delete $attr->{links};
    }
    
    if (my $leg = $self->attribute('legend')) {
        my @errs;
        if (ref($leg) eq 'HASH') {
            my %legend;
            while (my ($ltype, $tdat) = each %{$leg}) {
                my $fgType = 
                    $ltype eq 'nodes' ? 'Node' :
                    $ltype eq 'edges' ? 'Edge' :
                    $ltype eq 'text'  ? 'Graph' : undef;
                unless ($fgType) {
                    push @errs, "Unknown legend part '$ltype'";
                    next;
                }
                if (ref($tdat) eq 'ARRAY') {
                    my @parts;
                    foreach my $part (@{$tdat}) {
                        unless (ref($part) eq 'HASH') {
                            push @errs, "Non-hash component for legend part '$ltype' = '$part'";
                        }
                        my $attr = { %{$part} };
                        $self->remap_attributes
                            ($attr, 'canvasxpress', $fgType);
                        $attr->{color} = $self->rgba_color
                            ( $attr->{color}, $attr->{coloralpha} )
                            if ($attr->{color});
                        
                        $attr->{shape} = $cxShapeMap->{ lc($attr->{shape} || '') } || 'square' if (my $shape = $attr->{shape});
                        push @parts, $attr;
                    }
                    if ($#parts == -1) {
                        push @errs, "No legend components for '$ltype'";
                    } else {
                        $legend{$ltype} = \@parts;
                    }
                } else {
                    push @errs, "Legend part '$ltype' is not an array";
                }
            }
            $data->{legend} = \%legend;
        } else {
            push @errs, "CanvasXpress legend is not a hash";
        }
        $self->err("Errors building CanvasXpress legend:", @errs) unless ($#errs == -1);
    }

    $attr->{graphType}        = 'Network';
    # $attr->{eventKeys}        = 0;
    $attr->{backgroundType} ||= 'solid';
    $attr->{calculateLayout}  = 0 if ($args->{NOLAYOUT});
    

#    my $opts = {
#        showAnimation   => 1,
#        # preScaleNetwork => 0,
#        autoScaleFont   => 1,
#       #  nodeFontSize    => $attr->{nodeFontSize} || 8,
#    };
    require BMS::Utilities::Serialize;
    my $ser = $self->serializer();
    my $ind = $args->{PRETTY} ? 0 : undef;
    
    my $js  = 
        $ser->obj_to_json($data, $ind).",\n".
        $ser->obj_to_json($attr, $ind);

    if (my $cbh = $args->{CALLBACKS}) {
        my $lt  = $ser->literal_token();
        my %loc = %{$cbh};
        while (my ($evt, $method) = each %loc) {
            $loc{$evt} = $lt . $method;
        }
        $js .= ",\n". $ser->obj_to_json(\%loc, $ind);
    }

    if ($args->{BRACKET}) {
        $js = "[\n$js]\n";
    } else {
        $js ="new CanvasXpress( $domId,\n$js);\n";
        if (my $func = $args->{DEFER}) {
            $js = sprintf("function %s() {\n%s}\n", $func, $js);
            # $js = sprintf("function %s() {\n%s}\n", $func, "alert('$func');");
        }
    }

    return $js;
}

sub serializer {
    my $self = shift;
    unless ($self->{SERIALIZER}) {
        my $ser = $self->{SERIALIZER} = BMS::Utilities::Serialize->new();
        $ser->set_random_literal_token();
    }
    return $self->{SERIALIZER};
}

sub to_text {
    my $self = shift;
    return join("\n", sprintf("Graph %s [%d]", $self->name, $self->dbid),
                "# NODES:",$self->to_text_nodes(),
                "# EDGES:", $self->to_text_edges() );
}

sub to_text_nodes {
    my $self = shift;
    return join('', map { $_->to_text } $self->each_node);
}

sub to_text_edges {
    my $self = shift;
    return join('', map { $_->to_text } $self->each_edge);
}

sub layout {
    my $self = shift;
    my $args = $self->parseparams
        ( -scale => 50,
          -x     => 'x',
          -y     => 'y',
          -type  => 'grid',
          -integer => 1,
          @_);
    
    my $type = lc($args->{TYPE} || '');
    my $rv;
    if ($type =~ /rand/) {
        $rv = $self->_random_layout(%{$args});
    } elsif ($type =~ /spring|force/) {
        $rv = $self->_spring_layout(%{$args});
    } elsif ($type =~ /arc|radi/) {
        $rv = $self->_arc_layout(%{$args});
    } else {
        $rv = $self->_grid_layout(%{$args});
    }
    if ($args->{INTEGER}) {
        my @tags = ( $args->{X}, $args->{Y} );
        foreach my $node ($self->each_node) {
            map {
                $node->attribute($_, int(0.5 + $node->attribute($_))) } @tags;
        }
    }
    return $rv;
}

sub _arc_layout {
    my $self     = shift;
    my $args     = $self->parseparams( @_ );
    my $roots    = $args->{ROOT} || $args->{ROOTS} || [];
    my $rootNum  = $#{$roots} + 1;
    return unless ($rootNum);
    $self->benchstart();

    my $wedge    = $args->{WEDGE} || $args->{ARC} || 160;
    my $maxW     = $args->{MAXWEDGE} || $args->{MAXARC} || 15;
    my $radius   = $args->{RADIUS} || 150;
    my $iter     = $args->{ITERATIONS} || 4000;
    my $xtag     = $args->{X} || 'x';
    my $ytag     = $args->{Y} || 'y';
    my $seed     = $args->{SEED} || time() ^ ($$ + ($$<<15));
    my $maxLvl   = 0;
    my %done     = ();
    my @stack    = @{$roots};

    my %byLvl;
    while ($#stack != -1) {
        $maxLvl++;
        my %found;
        foreach my $node (@stack) {
            my $id = $node->numeric_id();
            next if ($done{$id}++);
            my @neigh = $node->neighbors();
            $node->attribute('ArcLevel', $maxLvl);
            push @{$byLvl{$maxLvl}}, $node;
            foreach my $n (@neigh) {
                $found{$n->numeric_id()} ||= $n;
            }
        }
        @stack = values %found;
    }
    my $aAx = $rootNum * $radius;
    my $bAx = $radius;
    my $prior;
    srand( $seed );
    my (%levels, %lookup);
    for my $lvl (1..$maxLvl) {
        my @nodes   = @{$byLvl{$lvl} || []};
        my $nodeNum = $#nodes + 1;
        next unless ($nodeNum);
        my $ldats = $levels{$lvl} = [];
        if ($lvl == 1) {
            # First set of points (the roots)
            # Just slap them in across three locations
            my $step  = 2 * $aAx / $rootNum;
            my $start = 0 - $aAx + $step / 2;
            for my $n (0..$#nodes) {
                my $node = $nodes[$n];
                my $nid  = $node->numeric_id();
                my $x = $start + $step * $n;
                my $y = $#nodes == 0 ? $bAx : 0;
                my $ldat = $lookup{$nid}{now} = {
                    lvl  => $lvl,
                    x    => $x,
                    y    => $y,
                    xy   => [$x, $y],
                    node => $nid,
                };
                push @{$ldats}, $ldat;
                $ldat->{i} = $#{$ldats};
            }
            next;
        }
        # all other levels, partition out discrete locations to use
        # First, calculate some fixed spots that nodes can inhabit
        # Provide more spaces, up to 25% beyond the number of nodes:
        my $xtra  = int($nodeNum * 0.50);
        # Always have at least 5 extra spaces:
        $xtra     = 5 if ($xtra < 5);
        my $slots = $nodeNum + $xtra;
        # Always have at least 8 spaces total:
        $slots    = 8 if ($slots < 8);
        # Determine how many degrees each slot should hold:
        my $step  = 180 / $slots;
        my $lvAng = 2 * atan2($bAx / 2, $bAx * $lvl) * 180 / pi;
        # warn "$lvl = $lvAng\n";
        $step     = $lvAng if ($step > $lvAng);
        $step = int(0.5 + 10 * $step) / 10 unless ($step < 2);
        my $rotate = 90;
        my $cnter  = 0;
        for (my $ang = -90; $ang <= 90; $ang += $step) {
            my $rad = 2 * pi * ($ang + $rotate) / 360;
            my $x = int(0.5 + 10 * $lvl * $aAx * cos($rad)) / 10;
            my $y = int(0.5 + 10 * $lvl * $bAx * sin($rad)) / 10;
            my $ldat = {
                ang  => int(0.5 + 100 * $ang) / 100,
                lvl  => $lvl,
                x    => $x,
                y    => $y,
                xy   => [$x, $y],
                r    => rand(1),
            };
            push @{$ldats}, $ldat;
            my $i = $ldat->{i} = $#{$ldats};
            if (0) {
                # Add visual arcs into the graph
                my $col  = "rgba(150,150,0,0.3)";
                my $node = $ldat->{arc} = $self->node("ArcLayout.$lvl.$i");
                $node->set_attributes( size => 0.01,
                                       color => $col,
                                       x    => $ldat->{x},
                                       y    => $ldat->{y}, );
                if ($i) {
                    my $prior = $ldats->[$i-1]{arc};
                    my $edge  = $self->edge($node, $prior);
                    $edge->set_attributes( color => $col );
                }
            }
        }
        # Now scatter the nodes across those points
        my @scramble = sort { $a->{r} <=> $b->{r} } @{$ldats};
        $self->death("Too few arc positions to hold all the nodes",
                     "$#scramble vs. $#nodes")
            if ($#scramble < $#nodes);
        for my $n (0..$#nodes) {
            my $node = $nodes[$n];
            my $ldat = $scramble[$n];
            my $nid  = $node->numeric_id();
            $ldat->{node}      = $nid;
            $lookup{$nid}{now} = $ldat;
        }
    }
    # Calculate all distances
    my @considering;
    while (my ($nid, $lu) = each %lookup) {
        my $node    = $self->node_by_numeric_id($nid);
        my @neigh   = map { $_->numeric_id() } $node->neighbors();
        $lu->{id}   = $nid;
        $lu->{conn} = \@neigh;
        $lu->{impr} = 999999;
        $lu->{r}    = rand(1);
        my $now     = $lu->{now};
        my $lvl     = $now->{lvl};
        unless ($lvl) {
            next
        }
        my @spots   = @{$levels{$lvl}};
        my $sNum    = $lu->{opts} = $#spots + 1;
        push @considering, $nid if ($sNum > 1);
        my ($x, $y) = @{$now->{xy}};
        $lu->{dist} = &_arc_distance( $nid, \%lookup, $x, $y);
    }

    # Shuffle the nodes around
    my $totNum  = $#considering + 1;
    my @lus     = map { $lookup{$_} } @considering;
    my $halfNum = int($totNum / 2);
    my $ti      = 15;
    my $start   = time;
    my $end     = $start + $ti;
    $iter       = 4000;

    my $getRandom = sub {
        my ($lu)  = @_;
        my $sNum  = $lu->{opts};
        my $i     = $lu->{now}{i};
        my $s     = $i;
        while ($s == $i) {
            $s = int(rand($sNum));
        }
        return $s;
    };

    my $getCloser = sub {
        my ($lu, $lookup)  = @_;
        my $nid       = $lu->{id};
        my $now       = $lu->{now};
        my $lvl       = $now->{lvl};
        my $i         = $now->{i};
        my @spots     = @{$levels{$lvl}};
        my ($x, $y)   = @{$now->{xy}};
        my ($dx, $dy) = &_arc_pull( $nid, $lookup, $x, $y );
        my ($px, $py) = ($x + $dx, $y + $dy);
        my @dists;
        my $nowDist = 0;
        for my $s (0..$#spots) {
            next if ($s == $i);
            my $spot = $spots[$s];
            my $dist = sqrt(($px - $spot->{x}) ** 2 + ($py - $spot->{y}) ** 2);
            $nowDist = $dist if ($s == $i);
            push @dists, [$dist, $s];
        }
        @dists = sort { $a->[0] <=> $b->[0] } @dists;
        @dists = splice(@dists, 0, 5);
        my $pick = int(rand($#dists+1));
        my $best = $dists[$pick];
        # warn join(' -> ', map { int($_) } $nowDist, $best->[0])."\n";
        return $best->[1];
    };

    my $sWindow = 0;
    for my $it (0..$iter) {
        my ($cutOff, $tryNum, $fetchNew) =
            (-1, $#lus, $getRandom);
        my @try;
        if (! ($it % 5)) {
            # Do a full scan of all nodes
            @try = @lus;
#        } elsif ($it < $nearStart) {
#            @try = @lus;
#            $fetchNew = $getCloser;
#        } elsif (! ($it % 13)) {
#            @try = @lus;
#            $fetchNew = $getCloser;
        } elsif (! ($it % 7)) {
            # Scan the longest nodes
            @try = sort { $b->{dist} <=> $a->{dist} || 
                              $a->{r} <=> $b->{r} } @lus;
            $tryNum = int($#try / 5);
        } else {
            # The remainder of the time, focus on the nodes that have
            # shown the most improvement, in the hopes we can make
            # them better still
            @try = sort { $b->{impr} <=> $a->{impr} || 
                              $a->{r} <=> $b->{r} } @lus;
            $cutOff = $try[0]{impr} * 0.7;
        }

        for my $t (0..$tryNum) {
            my $lu    = $try[$t];
            last if ($lu->{impr} < $cutOff);
            my $sNum  = $lu->{opts};
            my $nid   = $lu->{id};
            my $now   = $lu->{now};
            my $lvl   = $now->{lvl};
            my @spots = @{$levels{$lvl}};
            my $dist  = $lu->{dist};
            # Find a new spot
            my $cb      = $getRandom;#((rand(10) + $it / $iter) < 0.7) ? $getCloser : $getRandom;
            my $sCen    = &{$cb}( $lu, \%lookup );
            my $smin    = $sCen - $sWindow;
            my $smax    = $sCen + $sWindow;
            $smin = 0 if ($smin < 0);
            $smax = $#spots if ($smax > $#spots);
            my @attempts;
            for my $s ($smin..$smax) {
                my $targ    = $spots[$s];
                my $oid     = $targ->{node};
                my $newDist = &_arc_distance( $nid, \%lookup, @{$targ->{xy}});
                my $totDist = $newDist;
                my ($olu, $odist);
                if (defined $oid) {
                    # We are swapping with another node
                    $olu     = $lookup{ $oid };
                    $dist    += $olu->{dist};
                    $odist    = &_arc_distance( $oid, \%lookup, @{$now->{xy}});
                    $totDist += $odist;
                }
                my $improve   = $dist - $totDist;
                push @attempts, [ $improve, $targ, $newDist, $olu, $odist ];
            }
            @attempts = sort { $b->[0] <=> $a->[0] } @attempts;
            my ($improve, $targ, $newDist, $olu, $odist) = @{$attempts[0]};

            $lu->{r}      = rand(1);
            if ($improve <= 0) {
                unless (rand(10) < 0.7 - $it / $iter) {
                    $lu->{impr} *= 0.25;
                    next;
                }
            } else {
                $lu->{impr} = $improve;
            }
            # Ok we are moving the node
            $lu->{now}    = $targ;
            $lu->{dist}   = $newDist;
            $targ->{node} = $nid;
            my %doneDist = ( $nid => 1 );
            if ($olu) {
                # We also need to swap out the node that it replaced
                my $oid = $now->{node} = $olu->{id};
                $doneDist{$oid} = 1;
                $olu->{now}  = $now;
                $olu->{dist} = $odist;
                $olu->{impr} = $improve;
                $olu->{r}    = rand(1);
            } else {
                # Swapped into a blank spot
                delete $now->{node};
            }

            # We will also need to update the distance to all the
            # other nodes it is connected to!
            foreach my $look ($lu, $olu) {
                next unless ($look);
                foreach my $oid (@{$look->{conn}}) {
                    next if ($doneDist{$oid}++);
                    my $klu = $lookup{ $oid };
                    $klu->{dist} = &_arc_distance
                        ( $oid, \%lookup, @{$klu->{now}{xy}});
                }
            }
        }
        my $tot = 0; map { $tot += $_->{dist} } @lus;
        print join("\t", $it, int($tot))."\n";
        last if (time > $end);
    }

    while (my ($lvl, $spots) = each %levels) {
        foreach my $spot (@{$spots}) {
            last if ($spot->{node});
            if (my $node = $spot->{arc}) {
                $self->delete_node($node);
            }
        }
        for (my $s = $#{$spots}; $s >= 0; $s--) {
            my $spot = $spots->[$s];
            last if ($spot->{node});
            if (my $node = $spot->{arc}) {
                $self->delete_node($node);
            }
        }
    }
    
    # Finally, assign locations to each node
    while (my ($nid, $lu) = each %lookup) {
        my $node  = $self->node_by_numeric_id($nid);
        my $now   = $lu->{now};
        $node->attribute($xtag, $now->{x});
        $node->attribute($ytag, $now->{y});
    }
    $self->benchend();
}

sub _arc_distance {
    my ($nid, $lookup, $x, $y) = @_;
    my $lu = $lookup->{ $nid };
    my $dist    = 0;
    foreach my $oid (@{$lu->{conn}}) {
        my $xy = $lookup->{$oid}{now}{xy};
        unless ($xy) {
            # warn "Failed to find XY for $nid vs $oid\n";
            next;
        }
        my ($x2, $y2) = @{$xy};
        $dist += sqrt( ($x - $x2) ** 2 + ($y - $y2) ** 2);
    }
    return $dist;
}

sub _arc_pull {
    my ($nid, $lookup, $x, $y) = @_;
    my $lu = $lookup->{ $nid };
    my ($dX, $dY) = (0, 0);
    foreach my $oid (@{$lu->{conn}}) {
        my $xy = $lookup->{$oid}{now}{xy};
        unless ($xy) {
            # warn "Failed to find XY for $nid vs $oid\n";
            next;
        }
        my ($x2, $y2) = @{$xy};
        $dX += $x2 - $x;
        $dY += $y2 - $y;
    }
    return ($dX, $dY);
}

sub _random_layout {
    my $self     = shift;
    my $args     = $self->parseparams( @_ );
    my $seed     = $args->{SEED} || time() ^ ($$ + ($$<<15));
    my @nodes    = $self->each_node();
    my $nn       = $#nodes + 1;
    my $xtag     = $args->{X} || 'x';
    my $ytag     = $args->{Y} || 'y';
    my $maxRad   = $args->{RADIUS} || 
        int(3 * sqrt($nn) * ($args->{SCALE} || 0)) || 50;

    srand( $seed );
    my (@xs, @ys);
    foreach my $node (@nodes) {
        my $rad = rand($maxRad);
        my $deg = rand(2 * pi);
        my $x   = $rad * cos($deg);
        my $y   = $rad * sin($deg);
        $node->attribute($xtag, $x);
        $node->attribute($ytag, $y);
        push @xs, $x;
        push @ys, $y;
    }
    @xs = sort { $a <=> $b } @xs;
    @ys = sort { $a <=> $b } @ys;
    return { alg => 'grid',
             xmin => $xs[0],
             xmax => $xs[-1],
             ymin => $ys[0],
             ymax => $ys[-1],
         };
}

sub _spring_layout_OLD {
    my $self     = shift;
    my $args     = $self->parseparams( @_ );
    my $xtag     = $args->{X};
    my $ytag     = $args->{Y};
    my $restLen  = $args->{SCALE} || 50;
    my $iter     = $args->{ITERATIONS} || 1000;
    my $maxRepul = $args->{MAXREPULSE} || $restLen;
    my $moveSize = 0.8;

    my @nodes    = $self->each_node();
    my %id2ind   = map { $nodes[$_]->id => $_ } (0..$#nodes);
    $self->_random_layout(%{$args});
    my (@data, %edges);
    my $totEdges = 0;
    for my $i (0..$#nodes) {
        my $node = $nodes[$i];
        my @connections;
        foreach my $edge ($node->each_edge) {
            my $onode = $edge->other_node($node);
            push @connections, $id2ind{ $onode->id };
        }
        my $enum = $#connections + 1;
        $totEdges += $enum;
        $data[$i] = {
            x => $node->attribute( $xtag ) || 0,
            y => $node->attribute( $ytag ) || 0,
            i => $i,
            c => $enum + 0, # do not allow zero here
        };
        map { $edges{$i}{$_} = $edges{$_}{$i} = $restLen } @connections;
        # warn sprintf("%s : %s\n", $node->name, join(' | ', map { $nodes[$_]->name} @connections) || '');
    }
    my $test = 0;
    warn sprintf("* %d\n",$totEdges * $restLen) if ($test);

    my @small2large = sort { $data[$a]{c} <=> $data[$b]{c} } (0..$#data);
    for my $rep (1..$iter) {
        my $totLen = 0;
        foreach my $ei (@small2large) {
            my $idat = $data[$ei];
            my ($xi, $yi, $i) = ($idat->{x}, $idat->{y}, $idat->{i});
            my ($moveX, $moveY) = (0,0);
            # Consider this node in relation to all others:
            for my $ej (0..$#small2large) {
                next if ($ei == $ej);
                my $jdat = $data[$ej];
                my ($xj, $yj, $j) = ($jdat->{x}, $jdat->{y}, $jdat->{i});
                my $dx   = $xj - $xi;
                my $dy   = $xj - $yi;
                my $len  = sqrt( ($dx*$dx) + ($dy*$dy) );
                my $force = 0;
                if (my $ideal = $edges{$j}{$i}) {
                    # The j node is attached to the i node
                    # Move in an attempt to set spring to rest length
                    my $diff = $len - $ideal;
                    $force   = $diff * $moveSize;
                    $totLen += $len;
                } elsif ($len < $maxRepul) {
                    # Unconnected nodes that are too close
                    $force   = sqrt($maxRepul ** 2 - $len ** 2);
                }
                if ($force) {
                    # Relative weight - force less connected nodes to move more
                    my $relW = $jdat->{c} / ($jdat->{c} + $idat->{c});
                    $force  *= $relW;
                    my $rad  = $dx && $dy ? atan2($dy, $dx) : rand(2 * pi);
                    $moveX  += $force * cos($rad);
                    $moveY  += $force * sin($rad);
                    warn sprintf("  {%d,%d} = %d : %d = [%d,%d]\n",
                                 $dx, $dy, $rad, $len, $moveX, $moveY)
                        if ($test);
                }
            }
            $idat->{x} += $moveX;
            $idat->{y} += $moveY;
        }
        warn sprintf("* %d\n", $totLen) if ($test);
    }
    foreach my $dat (@data) {
        my $node = $nodes[$dat->{i}];
        $node->attribute($xtag, $dat->{x});
        $node->attribute($ytag, $dat->{y});

    }
    return { alg => 'grid',
             iter => $iter,
         };
}

sub _spring_layout {
    my $self     = shift;
    my $args     = $self->parseparams( @_ );
    my $xtag     = $args->{X};
    my $ytag     = $args->{Y};
    my $restLen  = $args->{SCALE} || 50;
    my $iter     = $args->{ITERATIONS} || 1000;
    my $moveSize = 0.8;
    my @nodes    = $self->each_node();
    my %id2ind   = map { $nodes[$_]->id => $_ } (0..$#nodes);
    my $randinf  = $self->_random_layout(%{$args});

    my ($startWid, $startHgt) =
        ($randinf->{xmax} - $randinf->{xmin}, 
         $randinf->{ymax} - $randinf->{ymin} );
    my @nInfo;
    my %data = ( repulse => $args->{MAXREPULSE} || $restLen,
                 forbid  => $args->{FORBID}     || 0,
                 fbcost  => $args->{FORBIDCOST} || $restLen ** 2);
    my $totEdges = 0;
    for my $i (0..$#nodes) {
        my $node = $nodes[$i];
        my @connections;
        foreach my $edge ($node->each_edge) {
            my $onode = $edge->other_node($node);
            push @connections, $id2ind{ $onode->id };
        }
        my $enum = $#connections + 1;
        $totEdges += $enum;
        $nInfo[$i] = {
            x => $node->attribute( $xtag ) || 0,
            y => $node->attribute( $ytag ) || 0,
            i => $i,
            c => $enum + 0, # do not allow zero here
        };
        map { $data{edges}{$i}{$_} = 
                  $data{edges}{$_}{$i} = $restLen } @connections;
    }
    $data{nodes} = \@nInfo;
    my $test = 0;

    my $circle    = 2 * pi;
    my $trialNum  = $args->{WEDGES} || 10;
    my $wedge     = $circle / $trialNum;
    my $minSpread = $restLen / 10;
    my $spread    = sqrt($startWid ** 2 + $startHgt ** 2) || 
        100 * $restLen || 1000;
    $self->set_attributes("FG Layout"      => 'spring',
                          "FG Iterations"  => $iter,
                          "FG Scale"       => $restLen,
                          "FG Wedge Count" => $trialNum );
    my $baseForce;
    for my $rep (1..$iter) {
        my $totForce = 0;
        my $frac     = (1 + $iter - $rep) / $iter;
        my $sMax     = $spread * $frac;
        $sMax        = $minSpread if ($sMax < $minSpread);
        my $sMin     = $sMax / 5;
        my $sRange   = $sMax - $sMin;
        for my $ei (0..$#nInfo) {
            my $idat = $nInfo[$ei];
            my $x    = $idat->{x};
            my $y    = $idat->{y};
            my $now = &_calc_force(\%data, $ei, $x,  $y);
            my $rad  = rand($circle);
            my @trials;
            $totForce += $now;
            for my $rot (1..$trialNum) {
                my $ang  = $rad + ($rot * $wedge);
                my $len  = $sMin + rand($sRange);
                my $nx   = $x + $len * cos($ang);
                my $ny   = $y + $len * sin($ang);
                my $frc  = &_calc_force(\%data, $ei, $nx, $ny);
                push @trials, [$frc, $nx, $ny];
            }
            my ($best) = sort { $a->[0] <=> $b->[0] } @trials;
            # warn join(" | ", map { int($_) } ($now, @trials));
            if ($now > $best->[0]) {
                $idat->{x} = $best->[1];
                $idat->{y} = $best->[2];
                #warn "$fp > $fn\n";
            }
        }
        $baseForce ||= $totForce;
        warn sprintf("[%3d] %12d (%6.2f%%)\n", $rep, $totForce, 100 * $totForce / $baseForce) if ($test);
    }
    for my $i (0..$#nInfo) {
        my $node = $nodes[$i];
        my $dat  = $nInfo[$i];
        $node->attribute($xtag, $dat->{x});
        $node->attribute($ytag, $dat->{y});
    }
    return { alg => 'grid',
             iter => $iter,
         };
}
sub _calc_force {
    my ($data, $i, $x, $y) = @_;
    my $idat  = $data->{nodes}[$i];
    my $edges = $data->{edges};
    my $repul = $data->{repulse};
    my $fb    = $data->{forbid};
    my $fbc   = $data->{fbcost};
    
    my $force = 0;
    for my $j (0..$#{$data->{nodes}}) {
        next if ($j == $i);
        my $jdat = $data->{nodes}[$j];
        my ($xj, $yj) = ($jdat->{x}, $jdat->{y});
        my $dx   = $xj - $x;
        my $dy   = $yj - $y;
        my $len  = sqrt( ($dx*$dx) + ($dy*$dy) );
        if (my $ideal = $edges->{$j}{$i}) {
            # The j node is attached to the i node
            # Move in an attempt to set spring to rest length
            my $diff = $len - $ideal;
            $force += $diff * $diff;
        } elsif ($len < $repul) {
            # Unconnected nodes that are too close
            if ($len < $fb) {
                $force   += $fbc;
            } else {
                my $diff  = $repul - $len;
                $force   += $diff ** 0.5;
            }
        } else {
            $force += 1 / $len;
        }
    }
    return $force;
}

sub _grid_layout {
    my $self = shift;
    my $args     = $self->parseparams( @_ );
    my $step     = $args->{SCALE};
    my $xtag     = $args->{X};
    my $ytag     = $args->{Y};
    my @nodes    = $self->each_node();
    my $nn       = $#nodes + 1;
    my $numCols  = int(0.99 + sqrt($nn));
    my (@xs, @ys);
    for my $i (0..$#nodes) {
        my $node = $nodes[$i];
        my $x = $step * ($i % $numCols);
        my $y = $step * int($i / $numCols);
        $node->attribute($xtag, $x);
        $node->attribute($ytag, $y);
        push @xs, $x;
        push @ys, $y;
    }
    @xs = sort { $a <=> $b } @xs;
    @ys = sort { $a <=> $b } @ys;
    return { alg => 'grid',
             xmin => $xs[0],
             xmax => $xs[-1],
             ymin => $ys[0],
             ymax => $ys[-1],
         };
}

sub to_text_file {
    my $self = shift;
    my $txt = <<EOF;
# Friendly Graph Text File
#   This file can be parsed by BMS::FriendlyGraph.pm to import rich metadata
#   for a graph or subgraph. It is broken into sections defining different
#   kinds of data describing the graph, its nodes, and its edges.

#   Comment lines (such as this one) are ignored and discarded
#   (so they will not survive a round-trip export)

#   Data lines will be broken into "phrases" which are separated by white space
#   If you wish to use a phrase with white space, it may be enclosed in quotes
#   Quoted phrases should have quotes and slashes escaped, as \" \' \\
#   If a phrase does not contain whitespace, it may be provided unquoted

#   For the Node Section, the first phrase defines a node by its name
#   The Edge Section uses the first two phrases to define the edge
#   All other phrases in a line will be interpreted as attribute key/val pairs
#   Attribues can be of the format KeyPhrase:ValPhrase or KeyPhrase=ValPhrase

EOF

    my @gattr = $self->attributes_for_text_file();
    unless ($#gattr == -1) {
        $txt .= "FGSECTION GRAPH : Contains attributes for the graph as a whole\n\n";
        $txt .= join("\n", @gattr)."\n\n";
    }
    
    $txt .= "FGSECTION NODES : Lists each node in the graph, as well as attributes they may have\n\n";
    foreach my $node ($self->each_node()) {
        $txt .= $self->_esc_text_file($node->name());
        if (my $nattr = $node->attributes_for_text_file()) {
            $txt .= " $nattr";
        }
        $txt .= "\n";
    }
    
    $txt .= "FGSECTION EDGES : Lists each edge in the graph, as well as attributes they may have\n\n";
    foreach my $edge ($self->each_edge()) {
        $txt .= join(' ', map { $self->_esc_text_file($_->name()) } $edge->nodes());
        if (my $nattr = $edge->attributes_for_text_file()) {
            $txt .= " $nattr";
        }
        $txt .= "\n";
    }
    return $txt;
}


our @textEscQuotes  = ('\\\\' => '__SLASH__',
                       '"' => '__DOUBLE_QUOTE__', 
                       "'" => '__SINGLE_QUOTE__' );
our $textQuoteTok   = "QUOTED_REGION<%d>";
# $textQuoteTok =~ s/\\$textEscQuotes[0]/foo/; die;
our $textFindTok    = "($textQuoteTok)"; $textFindTok =~ s/\%d/(\\d+)/;
our $textMeth = {
    NODES => \&_parse_text_nodes,
    EDGES => \&_parse_text_edges,
    GRAPH => \&_parse_text_graph,
    CLASS => \&_parse_text_class,
};

sub parse_text_file {
    my $self = shift;
    my $txt  = shift;
    $self->text_section("");
    my $method = undef;
    my $miss   = 0;
    foreach my $line (split(/[\n\r]+/, $txt)) {
        my $arr = $self->parse_text_line( $line );
        next if (!$arr || $#{$arr} == -1);
        if ($arr->[0] eq 'FGSECTION') {
            if (my $sect = $self->text_section( $arr->[1] )) {
                $method = $textMeth->{$sect};
            }
            next;
        }
        if ($method) {
            &{$method}($self, $arr);
        } else {
            $miss++;
        }
    }
    $self->err("$miss lines were ignored due to missing FGSECTION lines")
        if ($miss);
}

sub _parse_text_nodes {
    my $self = shift;
    my $arr  = shift;
    my $node = $self->node( shift @{$arr});
    if ($node) {
        $node->set_attributes( &data2params( @{$arr} ) );
    }
    return $node;
}

sub _parse_text_edges {
    my $self = shift;
    my $arr  = shift;
    my $edge = $self->edge( shift @{$arr}, shift @{$arr});
    if ($edge) {
        $edge->set_attributes( &data2params( @{$arr} ) );
    }
    return $edge;
}

sub _parse_text_graph {
    my $self = shift;
    my $arr  = shift;
    $self->set_attributes( &data2params( @{$arr} ) );
    return $self;
}

sub _parse_text_class {
    my $self = shift;
    my $arr  = shift;
    my $name = shift @{$arr};
    return undef unless ($name);
    my $params = &data2params($arr);
    my $hash   = {};
    while (my ($ks, $v) = each %{$params}) {
        my @lvls  = split(/\:\:/, $ks);
        my $final = pop @lvls;
        my $trg   = $hash;
        map { $trg = $trg->{$_} ||= {} } @lvls;
        $trg->{$final} = $v;
    }
    return $self->class($name, $hash);
}

sub text_section {
    my $self = shift;
    if (defined $_[0]) {
        my $nv = lc($_[0]);
        if ($nv =~ /node/) {
            $nv = 'NODES';
        } elsif ($nv =~ /edge/) {
            $nv = 'EDGES';
        } elsif ($nv =~ /graph/) {
            $nv = 'GRAPH';
        } elsif ($nv =~ /class/) {
            $nv = 'CLASS';
        } elsif (!$nv) {
            $nv = "";
        } else {
            $self->err("Failed to recognize text section '$nv'");
            return undef;
        }
        $self->{TEXT_SECTION} = $nv;
    }
    return $self->{TEXT_SECTION};
}

sub parse_line {
    my $self = shift;
    my $txt  = shift;
    return undef unless (defined $txt);
    # Ignore comments:
    return undef if ($txt =~ /^\#/);
    # Ignore leading and trailing whitespace:
    $txt =~ s/^\s+//; $txt =~ s/\s+$//;
    return undef if ($txt eq '');

    # Find and tokenize escaped quotes:
    for (my $i = 0; $i < $#textEscQuotes; $i +=2) {
        $txt =~ s/\\$textEscQuotes[$i]/$textEscQuotes[$i+1]/g;
    }
    # Find quoted sections and block them off:
    my @quoted;
    while ($txt =~ /(\"([^\"]+)\")/ ||
           $txt =~ /(\'([^\']+)\')/) {
        my ($orig, $targ) = ($1, $2);
        push @quoted, $targ;
        my $rep = sprintf($textQuoteTok, $#quoted);
        $txt =~ s/\Q$orig\E/$rep/g;
    }

    my @rv;
    # Split on whitespace 
    foreach my $bit (split(/[\s\n\r\t]+/, $txt)) {
        while ($bit =~ /$textFindTok/) {
            my ($orig, $ind) = ($1, $2);
            $bit =~ s/\Q$orig\E/$quoted[$ind]/g;
        }
        # Un-escape quotes:
        for (my $i = 0; $i < $#textEscQuotes; $i +=2) {
            $bit =~ s/\\$textEscQuotes[$i+1]/$textEscQuotes[$i]/g;
        }
        push @rv, $bit;
    }
    return wantarray ? @rv : \@rv;
    
}

sub data2params {
    my $self = shift;
    my ($arr) = @_;
    my %rv;
    foreach my $txt (@{$arr}) {
        next unless ($txt);
        if ($txt =~ /^([^\:]+)\:(.+)$/ ||
            $txt =~ /^([^\=]+)\=(.+)$/) {
            $rv{$1} = $2;
        }
    }
    return wantarray ? %rv : \%rv;
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::FriendlyGraph::Common;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

use strict;
use BMS::Utilities::ColorUtilities;
use BMS::Utilities::Benchmark;
use vars qw(@ISA);
@ISA   = qw(BMS::Utilities::ColorUtilities BMS::Utilities::Benchmark);


sub refText  { my $self = shift; return $self."" }
sub graph    { return shift->{graph} }
sub dbh      { return shift->graph->dbh }
sub database { return shift->graph->database }
sub format   { return lc(shift->graph->attribute('format') || 'generic'); }
sub label    { return shift->attribute('label', @_); }

*dbi = \&database;

*db_id = \&dbid;
sub dbid {
    my $self = shift;
    unless (defined $self->{DBID} || !$self->dbh) {
        $self->{DBID} = $self->database->text2id($self->name) || 0
            if ($self->dbh);
    }
    return $self->{DBID};
}

our $attrMap   = {
    Graph => {
        fontsize => {
            type => 'real',
            desc => 'The size of the text used when drawing the node label',
            map  => {
                canvasxpress => 'labelSize',
            },
        },
        labelcolor => {
            type => 'string',
            desc => 'The color used when drawing the label for the node',
            map  => {
                canvasxpress => 'labelColor',
            },
        },
        color => {
            type => 'string',
            desc => 'The background color of the graph',
            map  => {
                canvasxpress => 'background',
            },
        },
    },
    Node => {
        fill => {
            type => 'string',
            desc => 'Fill / background color of node',
            map  => {
                # cytoscape => 'node.fillColor',
                cytoscape => 'graphics::fill',
                graphviz  => 'fillcolor',
                canvasxpress => 'color',
            },
            aliases => ['color','bgcolor','backgroundcolor'],
        },
        x => {
            type => 'real',
            desc => 'X-coordinate for node placement',
            map  => {
                cytoscape => 'graphics::x',
                canvasxpress => 'x',
            },
        },
        y => {
            type => 'real',
            desc => 'Y-coordinate for node placement',
            map  => {
                cytoscape => 'graphics::y',
                canvasxpress => 'y',
            },
        },
        z => {
            type => 'real',
            desc => 'Z layer for stacking calculation',
            map  => {
                canvasxpress => 'z',
            },
        },
        border => {
            type => 'string',
            desc => 'The color applied to the border of the node',
            map  => {
                cytoscape => 'graphics::outline',
                graphviz  => 'color',
                canvasxpress => 'outline',
            },
            disallow => { 'color' => 1 },
        },
        penwidth => {
            type => 'real',
            desc => 'The thickness of the border around the node',
            map  => {
                graphviz  => 'penwidth',
                canvasxpress => 'outlineWidth',
            },
            disallow => {  },
            aliases => [ 'linewidth', 'borderwidth' ],
        },
        size => {
            type => 'real',
            desc => 'Size / scale of node',
            map  => {
                # cytoscape => 'node.size',
                cytoscape => ['graphics::w','graphics::h'],
                graphviz  => ['width', 'height'],
                canvasxpress => 'size',
            },
        },
        width => {
            type => 'real',
            desc => 'Size of the width of the node',
            map  => {
                cytoscape => 'graphics::w',
                graphviz  => 'width',
                canvasxpress => 'width',
            },
            aliases => [ 'w' ],
        },
        height => {
            type => 'real',
            desc => 'Size of the height of the node',
            map  => {
                cytoscape => 'graphics::h',
                graphviz  => 'height',
                canvasxpress => 'height',
            },
            aliases => [ 'h' ],
        },
        style => {
            type => 'string',
            desc => 'The style rendering of the border eg solid, dashed etc',
            map  => {
                cytoscape => 'graphics::cy:borderLineType',
                graphviz  => 'style',
            },
        },
        tooltip => {
            type => 'string',
            desc => 'Hover message to show when mouse is over the node',
            map  => {
                cytoscape => 'node.toolTip',
                graphviz  => 'tooltip',
                canvasxpress => 'tooltip',
            },
            disallow => { 'label' => 1 },
        },
        shape => {
            type => 'string',
            desc => 'Geometric shape of the node, for example "oval"',
            map  => {
                cytoscape => 'graphics::type',
                graphviz  => 'shape',
                canvasxpress => 'shape',
            },
        },
        label => {
            type => 'string',
            desc => 'The displayed label for the node',
            map  => {
                cytoscape => 'node.label',
                graphviz  => 'label',
                canvasxpress => 'label',
            },
            disallow => { 'name' => 1 },
        },
        labelcolor => {
            type => 'string',
            desc => 'The color used when drawing the label for the node',
            map  => {
                cytoscape => 'node.labelColor',
                graphviz  => 'fontcolor',
                canvasxpress => 'labelColor',
            },
        },
        fontsize => {
            type => 'real',
            desc => 'The size of the text used when drawing the node label',
            map  => {
                cytoscape => 'node.fontSize',
                graphviz  => 'fontsize',
                canvasxpress => 'labelSize',
            },
        },
        fontname => {
            type => 'string',
            desc => 'The font family to use',
            map  => {
                canvasxpress => 'nodeFont',
            },
            aliases => [ 'font','fontfamily' ],
        },
        foo => {
            type => 'string',
            desc => '',
            map  => {
                cytoscape => '',
                graphviz  => '',
            },
        },
        foo => {
            type => 'string',
            desc => '',
            map  => {
                cytoscape => '',
                graphviz  => '',
            },
        },
        foo => {
            type => 'string',
            desc => '',
            map  => {
                cytoscape => '',
                graphviz  => '',
            },
        },
    },
    Edge => {
        color => {
            type => 'string',
            desc => 'Fill / background color of edge',
            map  => {
                cytoscape => 'edge.color',
                graphviz  => 'color',
                canvasxpress => 'color',
            },
            aliases => [ 'border' ],
        },
        labelcolor => {
            type => 'string',
            desc => 'The color used when drawing the label for the edge',
            map  => {
                cytoscape => 'edge.labelColor',
                graphviz  => 'fontcolor',
            },
        },
        label => {
            type => 'string',
            desc => 'The displayed label for the edge',
            map  => {
                cytoscape => 'edge.label',
                graphviz  => 'label',
            },
        },
        fontsize => {
            type => 'real',
            desc => 'The size of the text used when drawing the edge label',
            map  => {
                cytoscape => 'edge.fontSize',
                graphviz  => 'fontsize',
            },
        },
        style => {
            type => 'string',
            desc => 'The style rendering of the line eg solid, dashed etc',
            map  => {
                cytoscape => 'edge.lineType',
                graphviz  => 'style',
                canvasxpress => 'type',
            },
        },
        width => {
            type => 'real',
            desc => 'Thickness / width of the line making up the edge',
            map  => {
                cytoscape => 'edge.lineWidth',
                graphviz  => 'style::setlinewidth',
                canvasxpress => 'width',
            },
            aliases => [ 'linewidth', 'penwidth' ],
        },
        tooltip => {
            type => 'string',
            desc => 'Hover message to show when mouse is over the edge',
            map  => {
                cytoscape => 'edge.toolTip',
                graphviz  => 'tooltip',
                canvasxpress => 'tooltip',
            },
            disallow => { 'label' => 1 },
        },
    },
};

# Expand name lookup for attribute map
while (my ($type, $tdat) = each %{$attrMap}) {
    my @generic = sort keys %{$tdat};
    my $map     = $tdat->{MAP} = { };
    my %special;
    foreach my $gAtt (@generic) {
        # $gAtt is the "generic" name for the attribute that FriendlyGraph
        # has chosen to use.
        # $gdat is the FG data structure for the generic attribute
        my $gdat = $tdat->{$gAtt};
        # $adat holds the mappings from the generic term to specific formats
        my $adat = $gdat->{map};
        # @src will contain all known synonyms for this object
        my @src  = ($gAtt);
        my $dis  = $gdat->{disallow}  || {};
        push @src, @{$gdat->{aliases} || []};

        # Note all the specific attribute values for this generic name:
        foreach my $targ (values %{$adat}) {
            push @src, $targ unless (ref($targ) || $dis->{$targ});
        }
        my %uniq = map {lc($_) => 1 } @src; @src = keys %uniq;

        # Map over to the generic attribute:
        map { $map->{generic}{$_} = $gAtt } @src;
        while (my ($frm, $targ) = each %{$adat}) {
            # $frm is a specific format, $trg is the specific name to use
            # Map the synonyms over to specific values
            if ($dis->{$targ}) {
                # The target attribute for this format is disallowed from
                # being used as a general alias, but we want to be sure
                # that within this format it can still be used;
                $special{$frm}{$targ} = $targ;
            }
            foreach my $src ( @src ) {
                if ($map->{$frm}{$src} && $map->{$frm}{$src} ne $targ) {
                    warn "Multiple $type attribute mappings for $frm:$src\n".
                        "  $map->{$frm}{$src} + $targ\n";
                } else {
                    $map->{$frm}{$src} = $targ;
                }
            }
        }
    }
    while (my ($frm, $srcs) = each %special) {
        while (my ($src, $targ) = each %{$srcs}) {
            $map->{$frm}{$src} = $targ;
        }
    }
}

sub map_attribute {
    my $self = shift;
    my ($attr, $format, $type) = @_;
    return undef unless ($attr);
    $format ||= $self->format();
    $type   ||= $self->type;
    my $rv    = $attrMap->{$type}{MAP}{lc($format)}{lc($attr)};
    $rv       = $attr unless (defined $rv);
    return ref($rv) ? @{$rv} : ($rv);
}

sub generic_attributes {
    my $self = shift;
    my @types = $_[0] ? ($_[0]) : qw(Graph Node Edge);
    my @rv;
    foreach my $type (@types) {
        my $tdat = $attrMap->{$type};
        next unless ($tdat);
        foreach my $gAttr (sort keys %{$tdat}) {
            next if ($gAttr eq 'MAP');
            my $gdat = $tdat->{$gAttr};
            push @rv, [$gAttr, $gdat->{desc}, $type, $gdat->{type}];
        }
    }
    return wantarray ? @rv : \@rv;
}

our $colorValues = {};
map { $colorValues->{$_}{canvasxpress} = 'rgba' }
qw(color background backgroundGradient1Color backgroundGradient2Color outline foreground shadowColor smpLabelColor smpTitleColor varLabelColor nodeHighlightColor smpHighlightColor varHighlightColor overlayFontColor legendColor decorationsColor axisTickColor axisTitleColor xAxisTickColor yAxisTickColor zAxisTickColor blockContrastEvenColor blockContrastOddColor labelColor nodeFontColor smpHairlineColor featureNameFontColor sequenceAColor sequenceCColor sequenceGColor sequenceMColor sequenceTColor trackNameFontColor wireColor dragAreaColor infoAreaColor resizeAreaColor resizeAreaColorCurrent showAnimationFontColor);

sub normalize_colors {
    my $self = shift;
    my ($hash, $format) = @_;
    return unless ($hash);
    $format ||= "";
    while (my ($key, $val) = each %{$hash}) {
        if (my $ch = $colorValues->{$key}) {
            my ($fmt) = $format ? ($ch->{$format}) : sort values %{$ch};
            if ($fmt) {
                if ($fmt eq 'rgb') {
                    $hash->{$key} = $self->rgb_color($key);
                } elsif ($fmt eq 'rgba') {
                    $hash->{$key} = $self->rgba_color
                        ($key, $hash->{$key.'alpha'});
                    delete $hash->{$key.'alpha'};
                }
            }
        }
    }
}

sub id {
    my $self = shift;
    if ($_[0]) {
        $self->{id} = $_[0];
    }
    return $self->{id} ||= $self->numeric_id;
}

sub numeric_id {
    my $self = shift;
    unless ($self->{num_id}) {
        $self->{num_id} = ++$idCounter;
    }
    return $self->{num_id};
}

*set_attribute = \&attribute;
sub attribute {
    my $self = shift;
    my ($key, $val, $setAsEmpty) = @_;
    return undef unless ($key);
    if (defined $val) {
        $self->{DB_CURRENT}{ATTR} = 0
            unless (defined $self->{attr}{$key} && 
                    $self->{attr}{$key} eq $val);
        if ($val eq '' && !$setAsEmpty) {
            delete $self->{attr}{$key};
        } else {
            $self->{attr}{$key} = $val;
        }
    }
    return $self->{attr}{$key};
}

sub mapped_attribute {
    my $self = shift;
    my ($key, $format, $type, $rv) = @_;
    my $val = $self->attribute($key);
    $rv  ||= {};
    if (!defined $val || $val eq '') {
        # Nothing to do
    } elsif ($format) {
        foreach my $mkey ($self->map_attribute( $key, $format, $type )) {
            $rv->{$mkey} = $val;
        }
    } else {
        # No format defined
        $rv->{$key} = $val;
    }
    return $rv;
}

sub write_attributes {
    my $self = shift;
    return -1 if ($self->{DB_CURRENT}{ATTR});
    my $gid  = $self->graph->dbid;
    return 0 unless ($gid);
    my $dbi  = $self->database;
    my $rv   = 1;
    my $id   = $self->dbid;
    while (my ($key, $val) = each %{$self->{attr}}) {
        if (! defined $val || $val eq '') {
            # A null value is a request to delete this attribute on write
            $rv = 0 unless ($dbi->delete_attribute($gid, $id, $key));
        } elsif(ref($val)) {
            $rv = -1 if ($rv);
        } else {
            $rv = 0 unless ($dbi->write_attribute($gid, $id, $key, $val));
        }
    }
    $self->{DB_CURRENT}{ATTR} = 1;
    return $rv;
}

sub read_attributes {
    my $self   = shift;
    my $dbname = shift;
    my $gid    = $dbname ? 
        $self->database->text2id($dbname) : $self->graph->dbid;

    # NEED TO DO
    # Add a mechanism to read attributes from ALL graphs, not just the
    # current one or a specific one. By passing an explicit 0 ?
    my $rv = {};
    if ($gid) {
        $rv = $self->database->read_attributes( $gid, $self->dbid( $dbname ) );
        # warn join("\n", $self->name, map { "$_ => $rv->{$_}" } keys %{$rv});
        $self->set_attributes( %{$rv} );
    }
    return wantarray ? %{$rv} : $rv;
}

sub delete_all_attributes {
    my $self = shift;
    my $deep = shift;
    my $gid  = $self->graph->dbid;
    return 0 unless ($gid);
    my $id   = $self->dbid;
    my $dbi  = $self->database;
    my $rv   = 1;
    if ($deep) {
        # If we are doing a deep delete, make sure we delete everything
        # that exists in the database
        $rv = 0 unless ($dbi->delete_attributes_deep($gid, $id));
    } else {
        # Otherwise, just delete the keys that are currently known to object
        foreach my $key ($self->each_attribute) {
            $rv = 0 unless ($dbi->delete_attribute($gid, $id, $key));
        }
    }
    $self->{attr} = {};
    return $rv;
}

# Many of these are listed as being numeric in the documentation, but
# Cytoscape appears to expect string arguments... hmmm...
our $defaultTypes = {
    'node.size'     => 'string',
    'node.width'    => 'string',
    'node.height'   => 'string',
    'node.fontSize' => 'string',
    'edge.fontSize' => 'string',
    'graphics'      => 'XML',
};

sub attribute_type {
    my $self = shift;
    my $key  = shift || '';
    if (my $val = shift) {
        $attrTypes{$key} = $val;
    }
    return $attrTypes{$key} || $defaultTypes->{$key} || 'string';
}

sub attribute_type_is_set {
    my $self = shift;
    my $key  = shift || '';
    return $attrTypes{$key} ? 1 : 0;
}

sub set_attributes {
    my $self = shift;
    my $atNum = $#_ + 1;
    if ($atNum % 2) {
        $self->err("Odd number of attributes passed to set_attributes()");
        return undef;
    }
    for (my $i = 0; $i <= $#_; $i += 2) {
        $self->attribute($_[$i], $_[$i+1]);
    }
    return $atNum / 2;
}

sub each_attribute {
    my $self = shift;
    return sort keys %{$self->{attr}};
}

sub all_attributes {
    # As above, but includes the keys for the generators
    my $self = shift;
    my @keys = sort keys %{$self->{attr}};
    foreach my $cbdat ($self->each_attribute_generator) {
        my ($key, $cb, $ow) = @{$cbdat};
        push @keys, $key unless (exists $self->{attr}{$key});
    }
    return @keys;
}

sub all_attribute_values {
    my $self = shift;
    my $opts = shift;
    my %kvs;
    foreach my $key ($self->all_attributes) {
        my $val = $self->generated_attribute($key, $opts);
        next unless (defined $val);
        $kvs{$key} = $val;
    }
    for (my $i = 0; $i <= $#_; $i += 2) {
        # User-passed values
        my $key = $_[$i];
        next if (!defined $key || defined $kvs{$key});
        my $val = $_[$i+1];
        $kvs{$key} = $val if (defined $val);
    }
    my $format = $self->format();
    if ($opts && $opts =~ /format=(\S+)/) {
        $format = lc($1);
    }
    $self->remap_attributes( \%kvs, $format ) unless ($format eq 'generic');
    return \%kvs;
}

sub remap_attributes {
    my $self = shift;
    my ($hash, $format, $type) = @_;
    if (my $cname = $hash->{class}) {
        # User is specifying a visual class for the node
        foreach my $cl (split(/\s+/, $cname)) {
            my $dat = $self->graph->class($cl);
            my @targets;
            while (my ($key, $val) = each %{$dat}) {
                push @targets, [ $hash, $key, $val ];
            }
            while (my $dat = shift @targets) {
                my ($targ, $key, $val) = @{$dat};
                if (ref($val)) {
                    # Nested target
                    while (my ($k, $v) = each %{$val}) {
                        push @targets, [ $targ->{$key} ||= {}, $k, $v ];
                    }
                } else {
                    $targ->{$key} = $val unless (defined $targ->{$key});
                }
            }
        }
    }

    return $hash unless ($hash && $format);
    $type ||= $self->type;
    foreach my $atIn (keys %{$hash}) {
        my $val = $hash->{$atIn};
        delete $hash->{$atIn};
        # Is there value in keeping hash keys that are undefined?
        next unless (defined $val);
        # warn "IN: $atIn\n";
        foreach my $atOut ($self->map_attribute($atIn, $format, $type)) {
            # warn "   OUT: $atOut\n";
            my @lvls  = split(/\:\:/, $atOut);
            my $final = pop @lvls;
            my $trg = $hash;
            foreach my $lvl (@lvls) {
                # Drill down a level in the hash
                my $exst = $trg->{$lvl};
                if ($exst && !ref($exst)) {
                    # The target level already exists as a scalar
                    my ($nk, $nv) = ($exst, "");
                    if ($format eq 'graphviz') {
                        if ($exst =~ /^(.+)\((.+)\)$/) {
                            ($nk, $nv) = ($1, $2);
                        }
                    }
                    $trg = $trg->{$lvl} = { $nk => $nv };
                } else {
                    # Does not exist, or is a hash (a reference, at least)
                    $trg = $trg->{$lvl} ||= {};
                }
            }
            # If this is a color attribute, we want to make sure we are
            # using the appropriate color markup for the format
            if (my $ch = $colorValues->{$final}) {
                if (my $fmt = $ch->{$format}) {
                    # For this key in this format, a specific markup is needed
                    if ($fmt eq 'rgb') {
                        $val = $self->rgb_color( $val );
                    } elsif ($fmt eq 'rgba') {
                        $val = $self->rgba_color
                            ( $val, $hash->{$atIn."alpha"} );
                        delete $hash->{$atIn."alpha"};
                    } else {
                        # Oops. We should probably complain somehow...
                    }
                }
            }
            $trg->{$final} = $val;
        }
    }
    if ($format && $format eq 'graphviz' && exists $hash->{fillcolor}) {
        # Graphviz needs odd flags to make some features work
        $hash->{style}{filled} ||= ""
            if (!$hash->{style} || ref($hash->{style}) eq 'HASH');
    }
    # $self->normalize_colors( $hash, $format );
}

sub generator_keys {
    my $self = shift;
    # General class, then lastly any class
    # This method is here in case we want to make a key for individual objects
    # ($self)
    return map { uc($_) } ($self->type, '');
}

sub each_attribute_generator {
    my $self  = shift;
    my @tkeys = ($#_ == -1) ? $self->generator_keys() : map { uc($_) } @_;
    my @rv;
    foreach my $type (@tkeys) {
        while (my ($key, $dat) = each %{$attrCallBacks{$type}}) {
            push @rv, [ $key, @{$dat} ] if ($dat && $dat->[0]);
        }
    }
    return @rv;
}

sub attribute_generator {
    my $self = shift;
    my ($key, $type, $callBack, $overWrite) = @_;
    return undef unless ($key);
    $type = uc($type || '');
    if (defined $callBack) {
        if ($callBack eq '') {
            # An empty string is a request to delete the generator
            delete $attrCallBacks{$type}{$key};
            &clear_AGOW($key, $type);
        } else {
            $attrCallBacks{$type}{$key}   = [$callBack, $overWrite ];
            if ($overWrite) {
                $attrCBoverWrite{$key}{$type} = $overWrite;
            } else {
                &clear_AGOW($key, $type);
            }
        }
    } elsif (defined $overWrite) {
        # Only updating the overwrite flag
        $attrCallBacks{$type}{$key}[1] = $overWrite;
        if ($overWrite) {
            $attrCBoverWrite{$key}{$type} = $overWrite;
        } else {
            &clear_AGOW($key, $type);
        }
    }
    my $rv = $attrCallBacks{$type}{$key} || [];
    return wantarray ? @{$rv} : $rv->[0];
}

sub clear_AGOW {
    my ($key, $type) = @_;
    delete $attrCBoverWrite{$key}{$type};
    my @remain = keys %{$attrCBoverWrite{$key}};
    return unless ($#remain == -1);
    delete $attrCBoverWrite{$key};
}

sub generated_attribute {
    my $self = shift;
    my ($key, $exclude) = @_;
    return undef unless ($key);

    my $direct = $self->attribute($key);
    my $types;
    if (defined $direct) {
        # If an explicit, direct value is defined use it if not overwritten
        return $direct unless (exists $attrCBoverWrite{$key});

        # At least one type has an overwrite generator set for this key
        $types = [ $self->generator_keys ];
        my $ow = 0;
        map {$ow++ if (exists $attrCBoverWrite{$key}{$_} &&
                       $attrCBoverWrite{$key}{$_}) } @{$types};
        return $direct unless ($ow);
        # This type is being overwritten by a generator
    }

    # By default, we will recover all types of generated values:
    my ($doCode, $doRef, $doScalar) = (1,1,1);
    if ($exclude) {
        # User wants to exclude certain types
        $exclude = lc($exclude);
        $doScalar = 0 if ($exclude =~ /scalar/);
        $doCode   = 0 if ($exclude =~ /code/);
        $doRef    = 0 if ($exclude =~ /ref/);
    }
    $types ||= [ $self->generator_keys ];
    foreach my $type (@{$types}) {
        if (my $dat = $attrCallBacks{$type}{$key}) {
            my ($cb, $ow) = @{$dat};
            if (my $ref = ref($cb)) {
                if ($ref eq 'CODE') {
                    return &{$cb}($self) if ($doCode);
                } else {
                    return $cb if ($doRef);
                }
            } else {
                return $cb if ($doScalar);
            }
        }
    }
    return undef;
}

sub esc_gv { 
    my $self = shift;
    my $val  = shift;
    return $val if ($val =~ /^(true|false)$/i ||
                    $val =~ /^[\+\-]?\d+$/  ||
                    $val =~ /^[\+\-]?\d?\.\d+$/);
    $val =~ s/\"/\\\"/g;
    $val =~ s/\n/\\n/g;
    return "\"$val\"";
}

my @xesc = ( amp  => '&',
             quot => '"',
             apos => "'",
             gt   => '>',
             lt   => '<', );

sub esc_xml {
    my $self = shift;
    my $val  = shift;
    return "" unless (defined $val);
    $val =~ s/[\n\r\t]/ /g;
    for (my $i = 0; $i < $#xesc; $i += 2) {
        $val =~ s/$xesc[$i+1]/\&$xesc[$i]\;/g;
    }
    return $val;
}

sub esc_text {
    my $self = shift;
    my $val  = shift;
    return "" unless (defined $val);
    $val =~ s/\"/\\\"/g;
    $val =~ s/\n/\\n/g;
    return $val =~ /[\'\"\s]/ ? "\"$val\"" : $val;
}

sub attribute_xgmml {
    my $self = shift;
    my $indent = shift || 0;
    my $pad    = "  " x $indent;
    my $xml    = "";
    my $attr   = $self->all_attribute_values('format=cytoscape');
    foreach my $key (sort keys %{$attr}) {
        my $val  = $attr->{$key};
        my $type = $self->attribute_type($key);
        if (my $rt = ref($val)) {
            # The attribute is not a string value
            if ($type eq 'XML') {
                # Type is defined as an embedded XML node
                if ($rt eq 'HASH') {
                    my @gattr;
                    foreach my $gkey (sort keys %{$val}) {
                        my $gval = $val->{$gkey};
                        if ($key =~ /graphics/i && $gval !~ /^\d+$/) {
                            if (my $hex = $self->hex_color($gval)) {
                                $gval = $hex;
                            }
                        }

                        push @gattr, sprintf
                            ("%s='%s'", $gkey, $self->esc_xml($gval))
                            if (defined $gval);
                    }
                    $xml .= sprintf("%s<%s %s />\n", $pad, $key,
                                    join(' ', @gattr)) unless ($#gattr == -1);
                }
            }
        } else {
            next if ($type eq 'real' && $val eq '');
            if ($key =~ /color/i) {
                if (my $hex = $self->hex_color($val)) { $val = $hex }
            }
            $xml .= sprintf("%s<att name='%s' value='%s' type='%s' />\n", $pad,
                            $self->esc_xml($key), $self->esc_xml($val),
                            $type);
        }
    }
    return $xml;
}

sub attribute_gml {
    my $self = shift;
    my $indent = shift || 0;
    my $pad    = "  " x $indent;
    my $rv     = "";
    my $attr   = $self->all_attribute_values();
    foreach my $key (sort keys %{$attr}) {
        my $val  = $attr->{$key};
        next if ($val eq '');
        my $type = $self->attribute_type($key);
        if ($type eq 'real') {
            next if ($val eq '');
        } else {
            # Have not been able to find what to do for quote escaping
            $val =~ s/\"//g;
            $val = "\"$val\"";
        }
        $rv .= sprintf("%s%s %s\n", $pad, $key, $val);
    }
    return $rv;
}

sub attribute_gv {
    my $self  = shift;
    my $opts  = shift || "";
    my $attr  = $self->all_attribute_values( "$opts format=graphviz", @_ );
    my @attrs = ();
    foreach my $key (sort keys %{$attr}) {
        my $val  = $attr->{$key};
        if (ref($val)) {
            my @bits;
            while (my ($k,$v) = each %{$val}) {
                next unless ($k);
                if ($v) {
                    push @bits, sprintf("%s(%s)", $k, $v);
                } else {
                    push @bits, $k;
                }
            }
            next if ($#bits == -1);
            $val = join(',', @bits);
        }
        if ($key =~ /color/i) {
            if (my $hex = $self->hex_color($val)) {
                $val = $hex;
            }
        }
        push @attrs, "$key=". $self->esc_gv($val);
    }
    return ($#attrs == -1) ? "" : " [".join(',', @attrs)."]";
}

sub attribute_text {
    my $self  = shift;
    my $attr  = $self->all_attribute_values();
    my @attrs = ();
    foreach my $key (sort keys %{$attr}) {
        my $val  = $attr->{$key};
        next if (ref($val));
        push @attrs, $self->esc_text($key).'='.$self->esc_text($val);
    }
    return ($#attrs == -1) ? "" : " ".join(' ', @attrs);
}

sub attributes_for_text_file {
    my $self = shift;
    my $attr  = $self->all_attribute_values();
    my @attrs;
    foreach my $key (sort keys %{$attr}) {
        my $val  = $attr->{$key};
        next if (!defined $val || ref($val));
        push @attrs, $self->_esc_text_file($key).
            ($val =~ /\:/ || $key =~ /\:/ ? "=" : ':').
            $self->_esc_text_file($val);
    }
    return wantarray ? @attrs : join(" ", @attrs) || "";
}

sub _esc_text_file {
    my $self = shift;
    my $txt = shift;
    return "" unless (defined $txt);
    $txt =~ s/[\n\r]+/ /g;
    return $txt unless ($txt =~ /[\s\"\\]/);
    for (my $i = 0; $i < $#textEscQuotes; $i +=2) {
        $txt =~ s/\\$textEscQuotes[$i]/$textEscQuotes[$i+1]/g;
    }
    return "\"$txt\"";
}
    

sub filter_nodes {
    my $self = shift;
    my $nodes = shift;
    return wantarray ? @{$nodes} : $nodes if ($#_ == -1);

    my $args = $self->parseparams( @_ );
    my $nodeCB = $self->dbi->extract_parenthetical_test
        ( $args->{NODEFILTER} || $args->{NODEATTR} || $args->{FILTER});

    my @keep;
    for my $n (0..$#{$nodes}) {
        my $node = $nodes->[$n];
        next if ($nodeCB && !&{$nodeCB}( $node ));
        push @keep, $node;
    }
    return wantarray ? @keep : \@keep;
}

sub filter_edges {
    my $self = shift;
    my $edges = shift || [];
    return wantarray ? @{$edges} : $edges if ($#_ == -1);

    my $args = $self->parseparams( @_ );
    my $nodeCB = $self->dbi->extract_parenthetical_test
        ( $args->{NODEFILTER} || $args->{NODEATTR} );
    my $edgeCB = $self->dbi->extract_parenthetical_test
        ( $args->{EDGEFILTER} || $args->{EDGEATTR} || $args->{FILTER});

    my @keep;
    for my $e (0..$#{$edges}) {
        my $edge = $edges->[$e];
        next if ($edgeCB && !&{$edgeCB}( $edge ));
        if ($nodeCB) {
            my $node = $edge->other_node($self);
            next unless ($node && &{$nodeCB}( $node ));
        }
        push @keep, $edge;
    }
    return wantarray ? @keep : \@keep;
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::FriendlyGraph::Node;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

use strict;
use vars qw(@ISA);
use Scalar::Util qw(weaken);
@ISA      = qw(BMS::FriendlyGraph::Common);

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my ($name, $graph) = @_;
    return undef unless ($name);
    my $self = {
        name   => $name,
        ucname => uc($name),
        attr   => {},
        edges  => {},
        eByN   => {},
        DB_CURRENT => undef,
    };
    bless ($self, $class);
    weaken($self->{graph} = $graph);
    return $self;
}

sub type  { return 'Node' }
sub name  { return shift->{name} }
sub uckey { return shift->{ucname} }

sub to_xgmml {
    my $self = shift;
    my $indent = shift || 0;
    my $pad    = "  " x $indent;
    my $lab    = $self->name(); # $self->attribute('label') || $self->name;
    my $xml    = sprintf("%s<node id='%s' label='%s'",$pad,
                         $self->esc_xml( $self->id ), $self->esc_xml( $lab ) );
    if (my $attr = $self->attribute_xgmml($indent + 1)) {
        $xml .= ">\n$attr$pad</node>\n";
    } else {
        $xml .= " />\n";
    }
    return $xml;
}

sub to_gml {
    my $self = shift;
    my $indent = shift || 0;
    my $pad    = "  " x $indent;
    my $rv     = sprintf("%snode [\n%s  id \"%s\"\n", $pad, $pad, $self->id());
    $rv       .= $self->attribute_gml( $indent + 1 );
    $rv .= "$pad]\n";
    return $rv;
}

sub edge {
    my $self = shift;
    return $self->graph( $self, @_);
}

sub to_graphviz {
    my $self = shift;
    return sprintf(" %4d%s\n", $self->numeric_id, 
                   $self->attribute_gv('scalar', label => $self->name));
}

sub to_canvasXpress {
    my $self = shift;
    my $expt = shift;
    my $attr = $self->all_attribute_values('format=canvasxpress');
    $attr->{id} = $self->canvasXpressId();
    return $attr;
}

sub canvasXpressId {
    my $self = shift;
    my $attr = $self->all_attribute_values('format=canvasxpress');
    return $attr->{id} || $self->id();
}

sub to_text {
    my $self = shift;
    return $self->esc_text($self->name()).$self->attribute_text()."\n";
}

sub _register_edge {
    my $self = shift;
    my $edge = shift;
    weaken($self->{edges}{$edge->refText} = $edge);
    my ($n1, $n2) = map { $_->uckey } $edge->nodes();
    my $key  = $self->uckey;
    if ($n1 eq $key) {
        weaken($self->{eByN}{$n2} = $edge);
    } else {
        weaken($self->{eByN}{$n1} = $edge);        
    }
}

sub unlink {
    my $self = shift;
    my ($oNode, $otherAlreadyDone) = @_;
    return undef unless ($oNode);
    my $name  = $self->name;
    my $graph = $self->graph;
    if (my $r = ref($oNode)) {
        if ($r =~ /FriendlyGraph::Node/) {
            $oNode = $oNode->name;
        } elsif ($r =~ /FriendlyGraph::Edge/) {
            my ($name1, $name2) = map { $_->name } $oNode->nodes();
            if ($name1 eq $name) {
                $oNode = $name2;
            } elsif ($name2 eq $name) {
                $oNode = $name1;
            } else {
                return undef;
            }
        } else {
            $self->err("Unable to unlink edge from node using $r", $name);
            return undef;
        }
    }
    my $edge;
    if (exists $self->{eByN}{uc($oNode)}) {
        $edge = $self->{eByN}{uc($oNode)};
        delete $self->{eByN}{uc($oNode)};
        delete $self->{edges}{$edge->refText};
        $edge->other_node($self)->unlink($name, 1) unless ($otherAlreadyDone);
    }
    return $edge;
}

sub each_edge {
    my $self = shift;
    my @rv = values %{$self->{edges}};
    return $self->filter_edges( \@rv, @_ );
}

sub neighbors {
    my $self = shift;
    my @rv;
    foreach my $edge ($self->each_edge(@_)) {
        push @rv, $edge->other_node($self);
        
    }
    return @rv;
}

sub all_connected_nodes {
    my $self = shift;
    my $isFirst = 0;
    
    if ($self->{CONNECT_LOOP}) {
        return if $self->{CONNECT_LOOP}{$self->id};
    } else {
        $isFirst = 1;
        $self->{CONNECT_LOOP} = {};
    }
    $self->{CONNECT_LOOP}{$self->id} = $self;
    foreach my $node ($self->neighbors( @_ )) {
        $node->all_connected_nodes( @_ );
    }
    if ($isFirst) {
        my @rv = values %{$self->{CONNECT_LOOP}};
        delete $self->{CONNECT_LOOP};
        return @rv;
    }
}

sub edge_count {
    my $self = shift;
    my @edges = $self->each_edge;
    return $#edges + 1;
}

sub write {
    my $self = shift;
    $self->write_attributes();
    return -1 if ($self->{DB_CURRENT}{WRITE});
    my $gid  = $self->graph->dbid;
    return 0 unless ($gid);
    my $id   = $self->dbid;
    if ($self->database->write_node( $gid, $id )) {
        $self->{DB_CURRENT}{WRITE}  = 1;
        $self->{DB_CURRENT}{DELETE} = 0;
        return $id;
    } else {
        return 0;
    }
}

sub read {
    my $self = shift;
    # warn "[FOO] ".$self->name;
    return $self->read_attributes(  );
}

sub read_edges {
    # Load all edges for a node from the database
    my $self  = shift;
    my $graph = $self->graph;
    my $dat   = $self->database->query_edges
        ( -dumpsql => 0,
          -node => $self->name(),
          -graph => $graph->name() );
    foreach my $set (@{$dat}) {
        my ($n1, $n2) = @{$set};
        $graph->edge($n1, $n2);
    }
    return $#{$dat} + 1;
}

sub delete {
    my $self = shift;
    my $deep = shift;
    return -1 if ($self->{DB_CURRENT}{DELETE});
    my $gid  = $self->graph->dbid;
    return 0 unless ($gid);
    $self->delete_all_attributes( $deep );
    $self->delete_all_edges( $deep );
    my $id   = $self->dbid;
    if ($self->database->delete_node( $gid, $id )) {
        $self->{DB_CURRENT}{WRITE}  = 0;
        $self->{DB_CURRENT}{DELETE} = 1;
        return $id;
    }
    return 0;
}

sub delete_edge {
    # Delete a specific named edge
    my $self = shift;
    my $onode = shift;
    my $edge = $self->graph->edge( $self, $onode );
    return undef unless ($edge);
    my $deep = shift;
    $self->graph->remove_edge($edge);
    $self->unlink( $edge->other_node($self) );
    return $edge->delete( $deep );
}

sub delete_all_edges {
    my $self = shift;
    my $deep = shift;
    $self->read_edges() if ($deep);
    my $rv = 1;
    my $graph = $self->graph();
    foreach my $edge ($self->each_edge()) {
        $rv = 0 unless ($edge->delete( $deep ));
        $graph->remove_edge($edge);
        $self->unlink( $edge->other_node($self) );
    }
    return $rv;
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::FriendlyGraph::Edge;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

use strict;
use Scalar::Util qw(weaken);
use vars qw(@ISA);
@ISA   = qw(BMS::FriendlyGraph::Common);

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my ($node1, $node2, $graph, $isDir) = @_;
    return undef unless ($node1 && $node2);
    $isDir    = $graph->attribute('directed') unless (defined $isDir);
    ($node1, $node2) = ($node2, $node1) if 
        (!$isDir && $node1->name gt $node2->name);
    my $self = {
        node1 => $node1,
        node2 => $node2,
        isDir => $isDir ? 1 : 0,
        DB_CURRENT => undef,
    };
    bless ($self, $class);
    weaken($self->{graph} = $graph);
    map { $_->_register_edge( $self ) } ($node1, $node2);
    return $self;
}


sub type     { return 'Edge' }
sub graph    { return shift->{graph} }
sub dbh      { return shift->graph->dbh }
sub directed { return shift->{isDir} }
sub name     {
    my $self = shift;
    my $tok  = $self->directed ? ' -> ' : ' -- ';
    return '{'.join($tok, map { $_->name } $self->nodes).'}';
}

*source = \&node1;
*target = \&node2;
sub node1 { return shift->{node1} }
sub node2 { return shift->{node2} }
sub nodes { my $self = shift; return ($self->{node1}, $self->{node2}) }

sub to_graphviz {
    my $self = shift;
    my ($n1, $n2) = $self->nodes();
    # Note that edge direction depends on if the whole graph is directed
    my $tok       = $self->graph->attribute('directed') ? '->' : '--';
    return sprintf(" %4d %s %4d%s\n", $n1->numeric_id, $tok, $n2->numeric_id,
                   $self->attribute_gv('scalar'));
}

sub to_xgmml {
    my $self = shift;
    my $indent = shift || 0;
    my $pad    = "  " x $indent;
    my $lab    = $self->attribute('label') || join
        (" : ", $self->node1->name, $self->node2->name);
    my $xml    = sprintf
        ("%s<edge id='%s' label='%s' source='%s' target='%s'",$pad,
         $self->esc_xml( $self->id ), $self->esc_xml( $lab ),
         $self->esc_xml($self->node1->id), $self->esc_xml($self->node2->id));
    if (my $attr = $self->attribute_xgmml($indent + 1)) {
        $xml .= ">\n$attr$pad</edge>\n";
    } else {
        $xml .= " />\n";
    }

    return $xml;
}

sub to_gml {
    my $self = shift;
    my $indent = shift || 0;
    my $pad    = "  " x $indent;
    my $rv     = sprintf("%sedge [\n%s  source \"%s\"\n%s  target \"%s\"\n", $pad,
                         $pad, $self->node1->id,
                         $pad, $self->node2->id,);
    $rv       .= $self->attribute_gml( $indent + 1 );
    $rv .= "$pad]\n";
    return $rv;
}

our $canvasXpressEdgeAttributes = 
    [qw(arrowPointSize dashLength dotLength lineDecoration type)];
sub to_canvasXpress {
    my $self = shift;
    my $expt = shift;
    my $attr = $self->all_attribute_values('format=canvasxpress');
    # die $self->branch($self);
    unless ($attr->{type}) {
        my $isDir = $attr->{directed};
        $isDir = $self->graph->attribute('directed') unless (defined $isDir);
        $attr->{type} = $isDir ? 'arrowHeadLine' : 'line';
    }
    $attr->{width} ||= 1; # Do we want to allow zero width lines?
    $attr->{id1}     = $self->node1->canvasXpressId();
    $attr->{id2}     = $self->node2->canvasXpressId();
    return $attr;
}

sub to_text {
    my $self = shift;
    return join(' ', map { $self->esc_text($_->name()) } $self->nodes).
        $self->attribute_text()."\n";
}

sub other_node {
    my $self = shift;
    my $node = shift;
    return undef unless ($node);
    $node = $self->graph->node($node) unless (ref($node));
    my $nid  = $node->id;
    my ($n1, $n2) = $self->nodes;
    if ($n2->id eq $nid) {
        return $n1;
    } elsif ($n1->id eq $nid) {
        return $n2;
    } else {
        return undef;
    }
}

sub write {
    my $self = shift;
    my $id  = $self->dbid;
    # warn sprintf("Edge %s = %d\n", $self->name, $id);
    return 0 unless ($id);
    $self->write_attributes();
    return $id;
}

sub read {
    my $self = shift;
    # Read the nodes, too:
    map { $_->read() } $self->nodes;
    return $self->read_attributes();
}

sub dbid {
    my $self = shift;
    my $dbname = shift || "";
    unless (defined $self->{DBID}{$dbname}) {
        my $dbi = $self->database();
        my $gid = $dbname ? $dbi->text2id( $dbname ) : $self->graph->dbid;
        if ($gid) {
            my ($n1, $n2) = map { $_->dbid } $self->nodes();
            $self->{DBID}{$dbname} = $self->database->edge2id($gid, $n1, $n2);
        }
    }
    return $self->{DBID}{$dbname};
}

sub delete {
    my $self = shift;
    return -1 if ($self->{DB_CURRENT}{DELETE});
    my $gid  = $self->graph->dbid;
    return 0 unless ($gid);
    my $deep = shift;
    $self->delete_all_attributes( $deep );
    my @nids = map { $_->dbid } $self->nodes;
    if ($self->database->delete_edge( $gid, @nids )) {
        $self->{DB_CURRENT}{WRITE}  = 0;
        $self->{DB_CURRENT}{DELETE} = 1;
        return @nids;
    } else {
        return 0;
    }
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::FriendlyGraph::DBH;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

use strict;
use BMS::FriendlyDBI;
use BMS::Utilities;
use Scalar::Util qw(weaken);

use vars qw(@ISA);
@ISA   = qw(BMS::Utilities BMS::FriendlyDBI);

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self = {
        LOCKED     => 1,
        SINGMSG    => {},
    };
    bless ($self, $class); 
    my $args   = shift;
    $ENV{PGPORT} = $args->{PORT} if ($args->{PORT});
    $ENV{PGHOST} = $args->{HOST} if ($args->{HOST});

    my $dbh    = $self->connect( $args );
    if ($dbh) {
        $dbh->schema( $self->schema );
        $dbh->make_all();
    } else {
        my $err = $fdbLastError;
        if ($err =~ //) {
            # The database does not exist, try to create it.
            $dbh = $self->create( $args );
        }
        warn "FriendlyGraph Failed to connect to database:\n  ".
            "$fdbLastError\n  " unless ($dbh);
    }
    $self->{DBH} = $dbh;
    weaken($self->{graph} = $args->{GRAPH});
    return $self;
}

sub type      { return "Database Handle"; }
sub unlock    { shift->{LOCKED} = 0; }
sub lock      { shift->{LOCKED} = 1; }
sub is_locked { return shift->{LOCKED}; }

sub single_message {
    my $self  = shift;
    my $msg   = shift;
    return unless ($msg);
    unless ($self->{SINGMSG}{$msg}++) {
        $self->msg($msg);
    }
}


sub dbh      { return shift->{DBH} }
sub graph    { return shift->{graph} }

our $idCache    = {};
our $cacheCount = 0;
sub text2id {
    my $self = shift;
    my $text = shift;
    my $id   = 0;
    return $id unless ($text);
    return $idCache->{$text} if (defined $idCache->{$text});
    # This SQL is also defined below in bulk_text2id()
    my $sth = $self->{STHS}{TEXT2ID} ||= $self->dbh->prepare
        ( -name => "Recover the numeric ID for an arbitrary piece of text",
          -sql  => "SELECT id FROM text2id WHERE name = ?",
          -level => 3);
    unless ($id = $sth->get_single_value($text)) {
        if ($self->is_locked) {
            $self->single_message("Can not create text id for '$text' - database is locked");
            return 0;
        }
        my $dbh   = $self->dbh;
        my $newID = $dbh->nextval('t2i_seq');
        my $upd   = $self->{STHS}{ADD_T2I};
        unless ($upd) {
            $upd = $self->{STHS}{ADD_T2I} = $dbh->prepare
                ( -name => "Create a new text / id pair in text2id",
                  -sql  => "INSERT INTO text2id (id, name) VALUES (?,?)",
                  -ignore => 'duplicate key value',
                  -level => 3);
        }
        # There is the small posibility that this will fail the unique
        # constraint if two procs are simultaneously trying to add the
        # same text.
        $upd->execute($newID, $text);
        # If the above failed, we should get the value from the other
        # process that won the race to insert:
        $id = $sth->get_single_value($text);
    }
    return &_set_id_cache( $text, $id );
}

sub _set_id_cache {
    my ($text, $id) = @_;
    if (++$cacheCount > 20000) {
        # Do not let the ID cache grow too large
        $idCache    = {};
        $cacheCount = 0;
    }
    return $idCache->{$text} = $id;
}

sub bulk_text2id {
    # WILL NOT CREATE NEW ENTRIES
    my $self = shift;
    my %rv;
    # This SQL is also defined above in text2id()
    my $sth = $self->{STHS}{TEXT2ID} ||= $self->dbh->prepare
        ( -name => "Recover the numeric ID for an arbitrary piece of text",
          -sql  => "SELECT id FROM text2id WHERE name = ?",
          -level => 3);
    foreach my $text (@_) {
        next unless ($text);
        my $id = $idCache->{$text};
        unless (defined $id) {
            if ($id = $sth->get_single_value($text)) {
                &_set_id_cache( $text, $id );
            }
        }
        $rv{$text} = $id || 0;
    }
    return \%rv;
}

sub id2text {
    my $self = shift;
    my $id   = shift;
    my $text = "";
    return $text unless ($id && $id =~ /^\d+$/);
    my $sth = $self->{STHS}{ID2TEXT};
    unless ($sth) {
        $sth = $self->{STHS}{ID2TEXT} = $self->dbh->prepare
            ( -name => "Recover the text associated with an id",
              -sql  => "SELECT name FROM text2id WHERE id = ?",
              -level => 3);
    }
    return $sth->get_single_value($id) || "";
}

sub edge2id {
    my $self = shift;
    my ($gid, $n1, $n2) = @_;
    my $id   = 0;
    return $id unless ($gid && $n1 && $n2);

    my $sth = $self->{STHS}{EDGE2ID};
    unless ($sth) {
        $sth = $self->{STHS}{EDGE2ID} = $self->dbh->prepare
            ( -name => "Recover the numeric ID for an edge",
              -sql  => "SELECT edge_id FROM edge WHERE graph_id = ? AND node1 = ? AND node2 = ?",
              -level => 3);
    }
    unless ($id = $sth->get_single_value($gid, $n1, $n2)) {
        if ($self->is_locked) {
            $self->single_message("Can not create edge id - database is locked");
            return 0;
        }
        my $dbh   = $self->dbh;
        my $newID = $dbh->nextval('t2i_seq');
        my $upd   = $self->{STHS}{ADD_E2I};
        unless ($upd) {
            $upd = $self->{STHS}{ADD_E2I} = $dbh->prepare
                ( -name => "Create a new edge entry",
                  -sql  => "INSERT INTO edge (graph_id, node1, node2, edge_id) VALUES (?,?,?,?)",
                  -ignore => 'duplicate key value',
                  -level => 3);
        }
        $upd->execute($gid, $n1, $n2, $newID);
        # There is the small posibility that this will fail the unique
        # constraint if two procs are simultaneously trying to add the
        # same edge.
        # If the above failed, we should get the value from the other
        # process that won the race to insert:
        $id = $sth->get_single_value($gid, $n1, $n2);
    }
    return $id;
}

sub bulkLoadText {
    my $self = shift;
    return 0 if ($#_ == -1);
    my %uniq;
    if (ref($_[0])) {
        map { $uniq{$_ || ''} = 1 } @{$_[0]};
    } else {
        map { $uniq{$_ || ''} = 1 } @_;
    }
    delete $uniq{''};
    my @texts = sort keys %uniq;
    if ($self->is_locked) {
        $self->single_message("Can not bulkLoadText - database is locked");
        return $#texts + 1;
    }
    my $dbh   = $self->dbh;
    my $tab   = "temp_list";
    my $col   = "member";
    $dbh->execute("CREATE TEMPORARY TABLE $tab ( $col text );");
    my $clear = $dbh->prepare("TRUNCATE TABLE $tab");
    my $add   = $dbh->prepare("INSERT INTO $tab ($col) values (?)");
    my $need  = $dbh->prepare
        ("SELECT $col FROM $tab WHERE NOT EXISTS ( ".
         "SELECT t.name FROM text2id t WHERE t.name = $col )");
    my $needed = 0;
    while (my @batch = splice(@texts, 0, 300)) {
        $clear->execute();
        map { $add->execute( $_ ) } @batch;
        my @missing = $need->selectcol_array();
        next if ($#missing == -1);
        my $data = [ map { [$dbh->nextval('t2i_seq'), $_] } @missing ];
        $dbh->insert_array('text2id', $data);
        $needed += $#{$data} + 1;
    }
    $dbh->execute("DROP TABLE $tab;");
    return $needed;
}

sub write_node {
    my $self = shift;
    my ($gid, $nid) = @_;
    return 0 unless ($gid && $nid);
    if ($self->is_locked) {
        $self->single_message("Can not write_node - database is locked");
        return 0;
    }
    unless ($self->{STHS}{ADD_NODE}) {
        $self->{STHS}{ADD_NODE} = $self->dbh->prepare
            ( -name => "Write a node to the database",
              -sql  => "INSERT INTO node (graph_id, node_id) VALUES (?,?)",
              -ignore => 'duplicate key value',
              -level => 3);
    }
    $self->{STHS}{ADD_NODE}->execute($gid, $nid);
    if ($@) {
        return 0;
    } else {
        return 1;
    }
}

sub delete_node {
    my $self = shift;
    my ($gid, $nid) = @_;
    return 0 unless ($gid && $nid);
    if ($self->is_locked) {
        $self->single_message("Can not delete_node - database is locked");
        return 0;
    }
    unless ($self->{STHS}{DEL_NODE}) {
        $self->{STHS}{DEL_NODE} = $self->dbh->prepare
            ( -name => "Delete node from the database",
              -sql  => "DELETE FROM node WHERE graph_id = ? AND node_id = ?",
              -level => 3);
    }
    $self->{STHS}{DEL_NODE}->execute($gid, $nid);
    if ($@) {
        return 0;
    } else {
        return 1;
    }
}

sub read_nodes {
    my $self = shift;
    my ($gid, $limit) = @_;
    my @rv;
    if ($gid) {
        my $sth = $self->dbh->prepare
            ( -name => "Read nodes for a graph",
              -sql  => "SELECT t1.name FROM node n, text2id t1 WHERE a.graph_id = ? AND a.obj_id = t1.id AND t1.name = ? AND a.key_id = t2.id",
                      -level => 3);

        $self->death("Need to actually do something");

    }
    return wantarray ? @rv : \@rv;
}

sub full_query {
    my $self = shift;
    my $args    = $self->parseparams( @_ );
    my $graph   = $args->{GRAPH};
    my $rules   = $args->{RULES} || $args->{RULE};
    my $limit   = 0;
    return undef unless ($rules);

    $graph ||= $self->graph;
    unless ($graph) {
        $self->err("No -graph defined for full_query()");
        return undef;
    }

    # Allow the user to define some variables
    # We need to pre-parse the rules line-by-line
    my %vars;
    my @lines = split(/[\n\r]+/, $rules);
    $rules = "";
    foreach my $line (@lines) {
        if ($line =~ /^\$([A-Za-z0-9_]+)\=(.*)$/) {
            # A variable value is being set
            # $FILTER1=level > 3 AND color = lime
            # $TARG_SET=|HiggledyPiggledy|
            $vars{uc($1)} = $2;
            next;
        }
        if ($line =~ /^\#/) {
            # Comment line, skip
            next;
        }
        while ($line =~ /(\$([A-Za-z0-9_]+))/) {
            # variable interpolation is needed
            my ($ori, $tag) = ($1, uc($2));
            my $rep = $vars{$tag} || "";
            $line =~ s/\Q$ori\E/$rep/g;
        }
        # Add information to the re-built rules set.
        $rules .= "$line\n";
    }

    my @gids    = $self->extract_graph_ids( $args );
    my @gnames  = map { $self->id2text($_) } @gids;
    my %nameAli = ( '*' => 'GRAPH',
                    '+' => 'ACTIVE', );

    # DO NOT FORGET
    # RegExps have problems capturing (.+) across new lines!
    
    my @actions;
    while ($rules) {
        if ($rules =~ /^([\s\n\r\t]+)/) {
            # Bunch of leading whitespace
            my $rep = $1;
            $rules =~ s/^\Q$rep\E//;
            next;
        }
        if ($rules =~ /^(CMD:([^\:]+):(.*?)\n)/) {
            my ($rep, $cmd, $val) = ($1, uc($2), $3);
            $rules =~ s/^\Q$rep\E//;
            push @actions, {
                type => 'CMD',
                cmd  => $cmd,
                val  => $val,
            };
            next;
        }
        if ($rules =~ /^((\>*)\|(\S*?)(\{([^\}]*)\})?\|(\>*))/) {
            # >|foo|  Store to foo
            # >>|foo| Add to foo
            # |foo|>  Read from foo
            # |foo{color = blue and score > 10}|> Filtered read from foo
            # |*| Special set referencing the current graph
            # |+| Special set referencing the active set
            my $rep  = $1;
            my $name = uc($3 || 'ACTIVE');
            my $dat  = {
                type     => 'NODE',
                savetok  => $2,
                name     => $nameAli{$name} || $name,
                readtok  => $6,
            };
            if (my $freq = $5) {
                if (my $f = $self->extract_parenthetical_test( $freq )) {
                    $dat->{nodeattr} = $f;
                    $dat->{nodetext} = $freq;
                } else {
                    $self->err("Failed to interpret node filter", $freq);
                    return undef;
                }
            }
            # If no tokens were provided, assume a read operation:
            $dat->{readtok} = '>' unless ($dat->{readtok} || $dat->{savetok});
            if ($dat->{nodetext} && !$dat->{readtok}) {
                $self->err("Useless filter '$dat->{nodetext}' on set $dat->{name}",
                     "Please note that filters are only used on read, not write");
            }
            push @actions, $dat;
            $rules =~ s/^\Q$rep\E//;
            next;
        }

        if ($rules =~ /^((\>+)\s+)/) {
            # > Store to active set
            # >> Append to active set
            my $rep = $1;
            push @actions, {
                type     => 'NODE',
                savetok  => $2,
                name     => 'ACTIVE',
            };
            $rules =~ s/^\Q$rep\E//;
            next;
        }
        if ($rules =~ /^(([\-\<]\@?[\-\*])\s*(\{([^\}]*)\})?\s*(\{([^\}]+)\})?\s*(\[([^\]]+)\])?\s*([\-\*]\@?[\-\>]))/) {
            # ---- All edges in database linked from ACTIVE set
            # --{rank > 3}-- filtered edge request
            # --{color = green}{size > 9}-- edge + node filter
            # --{}{size > 9}-- node filter only
            # ---> All 'forward' edges from live set
            # <--- All 'reverse' edges from live set
            # -**- All edges in current graph
            # -*{color = green}*- edge filtered graph query
            
            my ($rep, $lft, $ar1, $att1, $ar2, $att2, $nr, $nrt, $rgt) = ( $1,$2,$3,$4,$5,$6,$7,$8,$9 );
            # die join("\n", map { defined $_ ? "'$_'" : '-UNDEF-'} ('',$rep, $lft, $ar1, $att1, $ar2, $att2, $nr, $nrt, $rgt,'') );
            unless ($ar1 && $ar2) {
                # One or no specification was passed - only the edge is being
                # filtered, not the target node.
                $att1 ||= $att2;
                $att2 = undef;
            }
            my $dat = {
                type     => 'EDGE',
                left     => $lft || "",
                right    => $rgt || "",
                edgeattr => $att1,
                nodeattr => $att2,
                nodes    => $nrt,
            };
            
            push @actions, $dat;
            $rules =~ s/^\Q$rep\E//;
            next;
        }
        # warn $self->branch(\@actions);
        $self->err("Could not parse part of your rules:", $rules);
        return undef;
    }
    
    my %init = (QUERY => $args->{QUERY});
    if (my $sReq = $args->{SETS} || $args->{SET}) {
        if (ref($sReq) && ref($sReq) eq 'HASH') {
            map { $init{uc($_)} = $sReq->{$_} } keys %{$sReq};
        } else {
            $self->err("You have specified starting -sets in non-hash format",$sReq);
            return undef;
        }
    }
    my %sets;
    while (my ($name, $dat) = each %init) {
        next unless ($dat);
        $name = uc($name);
        if (my $r = ref($dat)) {
            if ($r eq 'ARRAY') {
                $sets{$name} = [ map { $graph->node($_) } @{$dat} ];
            } elsif ($r =~ /BMS::FriendlyGraph::Node/) {
                # Node object - make sure it is within this graph:
                $sets{$name} = [ $graph->node($dat) ];
            } else {
                $self->err("Improper set deffinition", $dat);
                return undef;
            }
        } else {
            # Single text node specification
            $sets{$name} = [ $graph->node($dat) ];
        }
    }
    return \%sets if ($#actions == -1);

    unless ($sets{QUERY}) {
        my $act1 = $actions[0];
        if ($act1->{type} eq 'NODE' && $act1->{readtok}) {
            # The first action is a set read, we will use that set
            if ($act1->{name} eq '*') {
                $sets{QUERY} = [ $graph->each_node ];
            } else {
                $sets{QUERY} = $sets{ $act1->{name} } ||= [];
            }
        } else {
            # No explicit query defined, use the whole graph
            $sets{QUERY} = [ $graph->each_node ];
        }
    }
    # map { $_->read() } @{$sets{QUERY}};
    $sets{ACTIVE} = $sets{QUERY} || [];
    $sets{RESULT} = [];

    my $splitter = '\s*\,\s*';
    my @log;
    for my $i (0..$#actions) {
        my $dat    = $actions[$i];
        my $type   = $dat->{type};
        my $resSet = $sets{RESULT};
        if ($type eq 'NODE') {
            my $name = $dat->{name};
            if (my $tok = $dat->{savetok}) {
                if (length($tok) == 1) {
                    # '>' Complete replace of set with result set
                    if ($name eq 'GRAPH') {
                        # This is the full graph
                        # I guess we should delete all nodes not in active set?
                        my @nodes = map { $_->name } $graph->each_node;
                        my %keep = map { $_->name => 1 } @{$resSet};
                        my $rm = 0;
                        foreach my $node (@nodes) {
                            next if ($keep{$node});
                            $graph->remove_node($node);
                            $rm--;
                        }
                        if ($rm == 0) {
                            push @log, "Graph update request results in no changes, as current result set is equivalent to graph content";
                        } else {
                            push @log, sprintf("Graph is updated to only contain current set, resulting in %d node%s being removed", $rm, $rm == 1 ? '' : 's');
                        }
                    } else {
                        $sets{$name} = [ @{$resSet} ];
                        push @log, sprintf("Reset set '%s' to contain the %d current result node%s", $name, $#{$resSet} + 1, $#{$resSet} == 0 ? '' : 's');
                    }
                } else {
                    # '>>' Append result set to this set
                    if ($name eq '*') {
                        # We really don't need to do anything, already in graph
                    } else {
                        my %cur = map { $_->name => $_ } @{$resSet};
                        map { delete $cur{ $_->name } } @{$sets{$name} ||= []};
                        my @nov = values %cur;
                        push @log, sprintf("%d novel node%s in result set are appended to set '%s'", $#nov + 1, $#nov == 0 ? '' : 's', $name);
                    }
                }
            }
            if (my $tok = $dat->{readtok}) {
                # Updating the active set
                my $src;
                if ($name eq 'GRAPH') {
                    $src = [ $graph->each_node() ];
                } else {
                    $src = $sets{$name} ||= [];
                }
                my @found;
                my $msg = sprintf("%d node%s in '%s'", $#{$src} + 1, $#{$src} == 0 ? '' : 's', $name);
                if (my $filter = $dat->{nodeattr} ) {
                    # We are only going to keep some of the nodes
                    foreach my $node (@{$src}) {
                        push @found, $node if (&{$filter}($node));
                    }
                    $msg .= sprintf(" (filtered by '%s' to %d)", $dat->{nodetext}, $#found + 1);
                } else {
                    # Taking the whole set
                    @found = @{$src};
                }
                if (length($tok) == 1) {
                    # We want to replace the entire set
                    $msg .= " used to replace active set";
                    $sets{ACTIVE} = [ @found ]
                } else {
                    # We are appending to the active set
                    my %fnd = map { $_->name => $_ } @found;
                    map { delete $fnd{$_->name} } @{$sets{ACTIVE}};
                    my @nov = values %fnd;
                    $msg .= sprintf(" were added to active set (%d new)",
                                    $#nov + 1);
                }
                push @log, $msg;
                $sets{RESULT} = [ @found ];
            }
        } elsif ($type eq 'EDGE') {
            my ($lft, $rgt) = ($dat->{left}, $dat->{right});
            # Determine directionality
            my $nodeParam = 'node';
            my $nodeIndex = 0;
            if ($lft =~ /\</ && $rgt !~ /\>/) {
                $nodeParam = 'node2';
                $nodeIndex = 2;
            } elsif ($lft !~ /\</ && $rgt =~ /\>/) {
                $nodeParam = 'node1';
                $nodeIndex = 1;
            }
            my @params;
            my @filtMsg;
            if (my $f = $dat->{edgeattr}) {
                push @params, (-edgefilter => $f);
                push @filtMsg, "edges matching '$f'";
            }
            if (my $f = $dat->{nodeattr}) {
                push @params, (-nodefilter => $f);
                push @filtMsg, "nodes matching '$f'";
            }
            if (my $targNodes = $dat->{nodes}) {
                my %uniq;
                while ($targNodes =~ /(\|([^\|]*)\|)/) {
                    my ($rep, $setN) = ($1, uc($2));
                    map { $uniq{$_} = undef } @{$sets{$setN} ||[]};
                    $targNodes =~ s/\Q$rep\E/\,/;
                }
                map { $uniq{$_} = undef } $self->extract_strings
                    ( split($splitter, $targNodes ) );
                my @nodeList = keys %uniq;
                if ($nodeParam eq 'node1') {
                    push @params, ( -node2 => \@nodeList );
                } elsif ($nodeParam eq 'node2') {
                    push @params, ( -node1 =>  );
                } else {
                    $nodeParam = 'node1';
                    push @params, ( -node2  => \@nodeList,
                                    -mirror => 1 );
                }
                my $lsz = $#nodeList + 1;
                push @filtMsg, "nodes in ".( $lsz < 5 ? "(".join(',', @nodeList).")" : sprintf("user set of %d node%s", $lsz, $lsz == 1 ? '' : 's'));
            }
            
            if (my $d = $args->{DUMPSQL}) { push @params, (-dumpsql => $d) }
            my $recurse = ($lft =~ /\@/ || $rgt =~ /\@/) ? 1 : 0;
            my $gOnly   = ($lft =~ /\*/ || $rgt =~ /\*/) ? 1 : 0;

            my @nodes = map { $_->name } @{$sets{ACTIVE}};
            my (%recovered, %searched);
            while ($#nodes != -1) {
                my %found;
                if ($gOnly) {
                    # Work only with the current graph
                    foreach my $name (@nodes) {
                        my $node = $graph->node($name);
                        foreach my $edge ($node->each_edge( @params )) {
                            if ($nodeIndex && $edge->directed) {
                                # Directed edge on a directed request
                                my @ns = $edge->nodes;
                                next unless ($ns[$nodeIndex-1]->name eq $name);
                            }
                            my $onode = $edge->other_node( $node );
                            my $oname = $onode->name;
                            $recovered{$oname} ||= $onode;
                            $found{$oname} ||= $onode;
                        }
                    }
                } else {
                    # Query against database
                    my %queries = map { $_ => 1 } @nodes;
                    my $search =  $self->query_edges
                        ( -desc      => "Recover edges from full query",
                          -graphname => \@gnames,
                          $nodeParam => \@nodes,
                          @params );
                    foreach my $row (@{$search}) {
                        my ($nm1, $nm2, $gn) = @{$row};
                        my $edge;
                        if ($nodeIndex) {
                            # Directed edge
                            $edge = $graph->edge( $nm1, $nm2, 1 );
                        } else {
                            $edge = $graph->edge( $nm1, $nm2 );
                        }
                        $edge->read_attributes( $gn );
                        if ($queries{$nm1}) {
                            my $onode = $graph->node($nm2);
                            $onode->read_attributes( $gn );
                            $recovered{$nm2} ||= $onode;
                            $found{$nm2} ||= $onode;
                        }
                        if ($queries{$nm2}) {
                            my $onode = $graph->node($nm1);
                            $onode->read_attributes( $gn );
                            $recovered{$nm1} ||= $onode;
                            $found{$nm1} ||= $onode;
                        }
                    }
                }
                last unless ($recurse);
                # Recurse using new nodes found in this iteration
                map { $searched{$_} = 1 } @nodes;
                # warn "Recursing on ".join(' + ', @nodes);
                @nodes = ();
                map { push @nodes, $_ unless ($searched{$_}) } keys %found;
            }
            $sets{RESULT} = [ values %recovered ];
            unless ($#{$sets{ACTIVE}} == -1) {
                my $msg = sprintf
                    ("%d active node%s searched %s",
                     $#{$sets{ACTIVE}} + 1, $#{$sets{ACTIVE}} == 0 ? '' : 's',
                     $recurse ? 'recursively ' : '');
                $msg .= $gOnly ? 'within the current network' : $#gnames == -1 ?
                    'against all graphs in DB' : 
                    sprintf("against graph%s %s",
                            $#gnames == 0 ? '' : 's', join(' + ', @gnames));
                $msg .= ", keeping ".join(" and ", @filtMsg) unless ($#filtMsg == -1);
                $msg .= sprintf(": %d node%s recovered", $#{$sets{RESULT}} + 1,
                                $#{$sets{RESULT}} == 0 ? '' : 's');
                push @log, $msg;
            }

        } elsif ($type eq 'CMD') {
            my ($cmd, $val) = ($dat->{cmd}, $dat->{val});
            if ($cmd eq 'GRAPHID' || $cmd eq 'GRAPHNAME') {
                @gids    = $self->extract_graph_ids( { $cmd => $val } );
                @gnames  = map { $self->id2text($_) } @gids;
                if ($#gids == -1) {
                    if ($val) {
                        push @log, "Tried to set graph IDs to '$val' but failed";
                    } else {
                        push @log, "All graphs in database available for query";
                    }
                } else {
                    push @log, sprintf("Contstraining queries to graph%s %s",
                                       $#gids == 0 ? '' : 's', join
                                       (' + ', @gnames));
                }
            } elsif ($cmd eq 'SPLIT' || $cmd eq 'SPLITTER') {
                $splitter = $val || '';
                if ($val) {
                    push @log, "Split token set to '$splitter'";
                } else {
                    push @log, "Splitting suppressed";
                }
            } elsif ($cmd eq 'READ') {
                my @objs;
                push @objs, $graph->each_node if ($val =~ /node/i);
                push @objs, $graph->each_edge if ($val =~ /edge/i);
                my @names = $#gnames == -1 ? ('') : @gnames;
                foreach my $gn (@names) {
                    map { $_->read_attributes($gn) } @objs;
                }
                push @log, sprintf("Attributes read from %s for %d object%s currently in graph", $#gnames == -1 ? ' default graph' : join(' + ', @gnames), $#objs + 1, $#objs == 0 ? '' : 's');
            } elsif ($cmd eq 'ADDNODE') {
                my $node = $graph->node($val);
                if ($#gids == -1) {
                    $node->read_attributes();
                } else {
                    map { $node->read_attributes($_) } @gnames;
                }
                push @log, sprintf("Manually added node: ". $node->name);
            } else {
                $self->err("Unrecognized command", "$cmd : $val");
            }
        }
    }
    return wantarray ? (\%sets, \@log, \@actions) : \%sets;
}

my $atN = 0;
sub query_nodes {
    my $self = shift;
    my $args = $self->parseparams( @_ );

    my %tables = map { $_ => 1 } ("node n", "text2id t1", "text2id t2");
    my @where  = ("t1.id = n.node_id", "t2.id = n.graph_id");
    my @binds;
    my @gids = $self->extract_graph_ids( $args );
    unless ($#gids == -1) {
        push @where, "n.graph_id ".($#gids == 0 ? "= $gids[0]" :
                                    " IN (".join(',', @gids).")");
    }
    if (my $nreq = $args->{NODE} || $args->{NODES}) {
        my %matchers;
        foreach my $string ($self->extract_strings( $nreq )) {
            if ($string =~ /[\%\?]/) {
                $string =~ s/_/\\_/g;
                $string =~ s/\?/_/g;
                $matchers{LIKE}{$string}++;
            } else {
                if (my $id = $self->text2id($string)) {
                    $matchers{IN}{$id}++;
                }
            }
        }
        if (my $idh = $matchers{IN}) {
            my @ids = sort {$a <=> $b} keys %{$idh};
            push @where, "n.node_id ". ($#ids == 0 ? "= $ids[0]" : "IN (".join(',', ).")");
        }
        if (my $lk = $matchers{LIKE}) {
            my @lb = sort keys %{$lk};
            push @where, map { "upper(t1.name) LIKE upper(?)" } @lb;
            push @binds, @lb;
        }
        
    }
    if (my $ef = $args->{ATTRFILTER} || $args->{NODEFILTER} || $args->{NODEATTR}) {
        my ($sq, $binds) = $self->extract_parenthetical_sql
            ( $ef, 'n', 'node_id');
        if ($sq) {
            push @where, "($sq)";
            push @binds, @{$binds};
        }
    }
    if (my $ef = $args->{EDGEFILTER} || $args->{EDGEATTR}) {
        my ($sq, $binds) = $self->extract_parenthetical_sql
            ( $ef, 'e', 'edge_id');
        if ($sq) {
            $tables{"edge e"} = 1;
            push @where, "((e.node1 = n.node_id OR e.node2 = n.node_id) AND $sq)";
            push @binds, @{$binds};
        }
    }

    my $sql = "SELECT DISTINCT t1.name AS Node, t2.name AS Graph";
    $sql .= " FROM ".join(', ', sort keys %tables);
    $sql .= " WHERE ".join(" AND ", @where);
    my $sth = $self->dbh->prepare
        ( -name => "Query nodes from the database",
          -sql  => $sql,
          -limit => $args->{LIMIT},
          -level => 1);
    warn $sth->pretty_print( @binds ) if ($args->{DUMPSQL});
    $sth->execute( @binds );
    my $rv = $sth->fetchall_arrayref();
    return wantarray ? ($rv, $sql) : $rv;
}

my $otherNodeCol = { node1 => 'node2', node2 => 'node1' };
sub query_edges {
    my $self = shift;
    my $args = $self->parseparams( @_ );

    my $sql = "SELECT DISTINCT t1.name AS Node1, t2.name AS Node2, t3.name AS Graph";
    my %tables = map { $_ => 1 } ("edge e", "text2id t1", 
                                  "text2id t2", "text2id t3");
    my @where  = ("t1.id = e.node1", "t2.id = e.node2", "t3.id = e.graph_id");
    my @binds;
    my @gids = $self->extract_graph_ids( $args );
    unless ($#gids == -1) {
        push @where, "e.graph_id ".($#gids == 0 ? "= $gids[0]" :
                                    " IN (".join(',', @gids).")");
    }

    my %recover;
    my @sorting;
    if (my $req = $args->{SORT}) {
        my @recs = ref($req) ? @{$req} : split(/\s*[\t\r\n\,]\s*/, $req);
        foreach my $rec (@recs) {
            my $xtra = "";
            my $type = "float";
            if ($rec =~ /(.+)\s+desc$/i) {
                $rec = $1;
                $xtra = " DESC";
            }
            my ($func, $pureKey) = $self->extract_sql_function($rec);
            if (my $id = $self->text2id($pureKey)) {
                my $cname = $pureKey;
                $cname =~ s/[^a-z0-9]/_/ig;
                unless ($recover{$cname}) {
                    my $char = substr($type, 0, 1);
                    my $tok  = join('', substr($type, 0, 1), 't', ++$atN);
                    $recover{$cname} = [$tok, $id, $func];
                    push @sorting, sprintf("%s%s", $cname, $xtra);
                }
            }
        }
    }
    if (my $req = $args->{RECOVER}) {
        my @recs = ref($req) ? @{$req} : split(/\s*[\t\r\n\,]\s*/, $req);
    }
    foreach my $cname (sort keys %recover) {
        my ($tt, $id, $func) = @{$recover{$cname}};
        $sql .= sprintf(", $func AS %s", "$tt.val", $cname);
        my $type = $tt =~ /^f/ ? "float" : "string";
        $tables{"attribute_$type $tt"} = 1;
        push @where, ("$tt.obj_id = e.edge_id", "$tt.key_id = $id");
    }

    my %matchers;
    foreach my $arg qw(NODE NODES NODE1 NODE2) {
        my $nreq = $args->{$arg};
        next unless ($nreq);
        my @targs = qw(node1 node2);
        if ($arg =~ /1$/) {
            @targs = ('node1');
        } elsif ($arg =~ /2$/) {
            @targs = ('node2');
        }
        foreach my $string ($self->extract_strings( $nreq )) {
            if ($string =~ /[\%\?]/) {
                $string =~ s/_/\\_/g;
                $string =~ s/\?/_/g;
                map { $matchers{$_}{LIKE}{$string}++ } @targs;
            } else {
                if (my $id = $self->text2id($string)) {
                    map { $matchers{$_}{IN}{$id}++ } @targs;
                }
            }
        }
    }
    my $nodeFilter  = $args->{NODEFILTER} || $args->{NODEATTR} || $args->{NODEPARAM};
    my $lrSpecified = $args->{NODE1} || $args->{NODE2} ? 1 : 0;
    my $invert      = $args->{MIRROR} && $lrSpecified ? 1 : 0;
    my @doInv       = (0);
    push @doInv, 1 if ($invert);
    my @nodeBits;
    foreach my $di (@doInv) {
        my @nbits;
        foreach my $targ (sort keys %matchers) {
            my $targCol = $di ? $otherNodeCol->{$targ} : $targ;
            my @wbits;

            if (my $idh = $matchers{$targ}{IN}) {
                my @ids = sort {$a <=> $b} keys %{$idh};
                push @wbits, "e.$targCol ". ($#ids == 0 ? "= $ids[0]" : "IN (".join(', ', @ids).")");
            }
            if (my $lk = $matchers{$targ}{LIKE}) {
                my $jtok = $targCol; $jtok =~ s/node/t/;
                my @lb = sort keys %{$lk};
                push @wbits, map { "upper($jtok.name) LIKE upper(?)" } @lb;
                push @binds, @lb;
            }
            my $nsql;
            if ($#wbits == 0) {
                $nsql = $wbits[0];
            } elsif ($#wbits != -1) {
                $nsql = "(".join(" OR ", @wbits).")";
            }
            if ($nodeFilter) {
                my $t2 = $otherNodeCol->{$targCol};
                my ($nfSql, $nfBinds) = $self->extract_parenthetical_sql
                    ( $nodeFilter, 'e', $t2);
                if ($nfSql) {
                    $nsql = "($nsql AND $nfSql)";
                    push @binds, @{$nfBinds};
                }
            }
            push @nbits, $nsql;
        }

        if ($#nbits == -1) {
            # There were no node-specific requests
            if ($nodeFilter) {
                # We need to add attribute filters for the nodes, however
                die "FAILED TO CODE NODE FILTER IN ABSENCE OF EXPLICIT NODES!";
            }
        } elsif ($#nbits == 0) {
            # Single node clause
            push @nodeBits, $nbits[0];
        } else {
            if ($lrSpecified || $args->{PAIRED}) {
                # Explict left/right specification, or a request that only
                # paired matches be found
                push @nodeBits, "(".join(" AND ", @nbits).")";
            } else {
                # If either left or right nodes match it is ok
                push @nodeBits, "(".join(" OR ", @nbits).")";
            }
        }

    }


    if ($#nodeBits == 0) {
        push @where, $nodeBits[0];
    } elsif ($#nodeBits != -1) {
        push @where, "(".join(" OR ", @nodeBits).")";
    }
    
    if (my $ef = $args->{ATTRFILTER} || $args->{EDGEFILTER} || 
        $args->{EDGEATTR}) {
        my ($sq, $binds) = $self->extract_parenthetical_sql
            ( $ef, 'e', 'edge_id');
        if ($sq) {
            push @where, "($sq)";
            push @binds, @{$binds};
        }
    }
    $sql .= " FROM ".join(', ', sort keys %tables);
    $sql .= " WHERE ".join(" AND ", @where);
    $sql .= " ORDER BY ".join(", ", @sorting) unless ($#sorting == -1);

    my $sth = $self->dbh->prepare
        ( -name => $args->{DESC} || "Query edges from the database",
          -sql  => $sql,
          -limit => $args->{LIMIT},
          -level => 1);
    if ($args->{EXPLAIN}) {
        warn $sth->explain_text( \@binds );
    } elsif ($args->{DUMPSQL}) {
        warn $sth->pretty_print( @binds );
    }
    $sth->execute( @binds );
    my $rv = $sth->fetchall_arrayref();
    return wantarray ? ($rv, $sql) : $rv;
}

sub extract_graph_ids {
    my $self = shift;
    my $args = shift;
    my @gids;
    if (my $greq = $self->{GRAPHID} || $self->{GRAPH_ID}) {
        if (my $r = ref($greq)) {
            if ($r eq 'ARRAY') {
                # Stringify it for parsing below
                $greq = join("\t", @{$greq});
            } else {
                $self->err("Unknown object provided as Graph ID", $greq);
                $greq = "";
            }
        }
        foreach my $gid (split(/[\n\r\t\s\,]/, $greq)) {
            if ($gid && $gid =~ /^\d+$/) {
                push @gids, $gid;
            } else {
                $self->err("Non-numeric Graph ID provided", $gid) if ($gid);
            }
        }
    }
    if (my $greq = $args->{GRAPH} || $args->{GRAPHNAME}) {
        my ($ids, $errs) = $self->extract_ids( $greq );
        push @gids, @{$ids};
        $self->err("Failed to recover ID for graph", @{$errs})
            unless ($#{$errs} == -1);
    }
    return wantarray ? @gids : \@gids;
}

my %funcs = (ABS => 'abs(%s)',
             UPPER => 'upper(%s)');
sub extract_sql_function {
    my $self = shift;
    my $key  = shift;
    my $func = '%s';
    while (1) {
        my $repeat = 0;
        while (my ($tag, $frm) = each %funcs) {
            if ($key =~ /^$tag\:(.+)$/) {
                $key = $1;
                $func = sprintf($frm, $func);
                $repeat = 1;
                last;
            }
        }
        last unless ($repeat);
    }
    return ($func, $key);
}

my @validOps = ('!=','==','>=','<=','=','>','<',);
my $opMap = {
    SQL => {
        '==' => '=',
    },
    NUM => {
        '=' => '==',
    },
    TXT => {
        '==' => 'eq',
        '!=' => 'ne',
        '>=' => 'ge',
        '<=' => 'le',
        '>'  => 'gt',
        '<'  => 'lt',
        '='  => 'eq',
    },
};
my $boolMap = {
    SQL => {
        '&&' => 'AND',
        '||' => 'OR',
    },
    CODE => {
        'AND' => '&&',
        'OR'  => '||',
    },
};
sub extract_operation {
    my $self = shift;
    my ($text, $isCode) = @_;
    my ($nonOp, $key, $op, $val) = ("");
    while ($text) {
        if ($text =~ /^([\s\(\)]+)/) {
            # Parentheses
            my $rep = $1;
            $text   =~ s/\Q$rep\E//;
            $rep    =~ s/\s+/ /g;
            $nonOp .= $rep;
            next;
        }
        if ($text =~ /^([\s]*(AND|OR|\&\&|\|\|)[\s]+)/i) {
            # Boolean
            my $rep = $1;
            my $boo = uc($2);
            $text   =~ s/\Q$rep\E//;
            $nonOp .= " ".
                ($boolMap->{ $isCode ? 'CODE' : 'SQL' }{$boo} || $boo) . "";
            next;
        }
        my $phrase = $text;
        if ($text =~ /(.+?)([\s]*\)|[\s]+AND[\s]+|[\s]+OR[\s]+)/i) {
            # Take the phrase only until the next paren/and/or
            $phrase = $1;
        }
        $text =~ s/\Q$phrase\E//;

        for my $o (0..$#validOps) {
            my $opt = $validOps[$o];
            if ($phrase =~ /^[\s]*(.+?)[\s]*(\Q$opt\E)[\s]*(.+?)[\s]*$/) {
                ($key, $op, $val) = ($1, $2, $3);
                last;
            }
        }
        last;
    }
    return ($nonOp, $key, $op, $val, $text);
}

sub extract_parenthetical_test {
    my $self = shift;
    my ($text) = @_;
    return undef unless ($text);
    my $code = "";
    # Non-whitelisted user values are going to be stored in an array to
    # avoid tainting when eval() is called
    my %strings;
    my $strCount = 0;

    while ($text) {
        my ($nonOp, $key, $op, $val, $res) = $self->extract_operation($text,1);
        $text = $res;
        return undef unless ($op || $nonOp);
        $code .= $nonOp;
        next unless ($op);
        my $keyInd = $strings{$key} ||= ++$strCount;
        my $valInd = $strings{$val} ||= ++$strCount;

        if (&is_num($val)) {
            $op = $opMap->{NUM}{$op} || $op
        } else {
            $op = $opMap->{TXT}{$op} || $op
        }
        $code .= sprintf('(defined $obj->attribute($ss[%d]) && $obj->attribute($ss[%d]) %s $ss[%d])'."\n", $keyInd, $keyInd, $op, $valInd);
    }
    my @ss;
    while (my ($string, $ind) = each %strings) {
        $ss[$ind] = $string;
    }
    return undef if ($code =~ /^\s*$/);
    # warn "\nCODE:\n$code\n\n";
    return eval
        ("sub {\n  my \$obj = shift \@_;\n  return ($code) ? 1 : 0;\n}\n");
}

sub extract_parenthetical_sql {
    my $self = shift;
    my ($text, $tabTok, $tabCol) = @_;
    my $sql  = "";
    my @binds;
    my %tables;
    while ($text) {
        my ($nonOp, $key, $op, $val, $resid) = $self->extract_operation($text);
        # warn "----\n$text\n---- [$nonOp] ($key, $op, $val)\n$resid\n****\n\n";
        $text = $resid;
        return ("") unless ($op || $nonOp);
        $sql .= $nonOp;
        next unless ($op);

        my ($func, $pureKey) = $self->extract_sql_function($key);

        my $id = $self->text2id($pureKey);
        my $tab = "attribute_float";
        unless (&is_num($val)) {
            $tab = "attribute_string";
            $val = "" if ($val eq "''" || $val eq '""');
        }
        my $tt = "at".++$atN;
        $tables{"$tab $tt"} = $tt;
        my $test = sprintf($func, "$tt.val");
        $sql .= sprintf(" EXISTS (SELECT val FROM %s %s WHERE %s.obj_id = %s.%s AND %s.graph_id = %s.graph_id AND %s.key_id = %d AND %s %s ?) ", $tab, $tt, $tt, $tabTok, $tabCol, $tabTok, $tt, $tt, $id, $test, $opMap->{SQL}{$op} || $op);
        push @binds, $val;
    }
    return ($sql, \@binds );
}

sub extract_strings {
    my $self = shift;
    my @strings;
    return @strings if ($#_ == -1);
    my @reqs;
    if (my $r = ref($_[0])) {
        if ($r eq 'ARRAY') {
            # Array reference
            @reqs = @{$_[0]};
        } else {
            @reqs = @_;
        }
    } else {
        @reqs = @_;
    }
    foreach my $text (@reqs) {
        next unless ($text);
        if (my $r = ref($text)) {
            if ($r =~ /BMS::FriendlyGraph/ && $text->can('name')) {
                push @strings, $text->name;
            }
            next;
        }
        foreach my $qt ('"',"'") {
            while ($text =~ /($qt\s*([^\"]*?)\s*$qt)/) {
                my $rep = $1;
                push @strings, $2 if ($2);
                $text =~ s/\Q$rep\E//;
            }
        }
        $text =~ s/^\s+//; $text =~ s/\s+$//;
        push @strings, $text if ($text);
    }
    return @strings;
}

sub extract_ids {
    my $self = shift;
    my (@ids, @errs);
    foreach my $string ($self->extract_strings(@_)) {
        if (my $id = $self->text2id($string)) {
            push @ids, $id;
        } else {
            push @errs, $string;
        }
    }
    return wantarray ? (\@ids, \@errs) : \@ids;
}

sub delete_named_edge {
    my $self = shift;
    my ($gn, $nn1, $nn2, $dir) = @_;
    my ($gid, $nid1, $nid2) = map { $self->text2id( $_ ) } ($gn, $nn1, $nn2);
    my $rv = $self->delete_edge( $gid, $nid1, $nid2 );
    unless ($dir) {
        my $rv2 = $self->delete_edge( $gid, $nid2, $nid1 );
        $rv = $rv2 if ($rv);
    }
    return $rv;
}

sub delete_edge {
    my $self = shift;
    my ($gid, $nid1, $nid2) = @_;
    return 0 unless ($gid && $nid1 && $nid2);
    if ($self->is_locked) {
        $self->single_message("Can not delete_edge - database is locked");
        return 0;
    }
    unless ($self->{STHS}{DEL_EDGE}) {
        $self->{STHS}{DEL_EDGE} = $self->dbh->prepare
            ( -name => "Delete edge from the database",
              -sql  => "DELETE FROM edge WHERE graph_id = ? AND node1 = ? AND node2 = ?",
              -level => 3);
    }
    $self->{STHS}{DEL_EDGE}->execute($gid, $nid1, $nid2);
    if ($@) {
        return 0;
    } else {
        return 1;
    }
}

our @attrTables = qw(float string);

sub is_num {
    return shift =~ /^[\+\-]?(\d+|\d*\.?\d+([eE][\+\-]?\d+)?)$/ ? 1 : 0;
}

sub write_attribute {
    my $self = shift;
    my ($gid, $oid, $key, $val) = @_;
    return 0 unless ($gid && $oid && $key);
    return 0 unless (defined $val);
    if ($self->is_locked) {
        $self->single_message("Can not write_attribute - database is locked");
        return 0;
    }
    my $type = &is_num($val) ? 'float' : 'string';

    my $add = $self->{STHS}{"ADD_ATTR_$type"} ||= $self->dbh->prepare
        ( -name => "Add $type attribute for object",
          -sql  => 
          "INSERT INTO attribute_$type (graph_id, obj_id, key_id, val)".
          " VALUES (?,?,?,?)",
          -ignore => 'duplicate key value',
          -level => 3);

    my $chk = $self->{STHS}{"CHECK_ATTR_$type"} ||= $self->dbh->prepare
        ( -name => "Check $type attribute for object parameter",
          -sql  => "SELECT val FROM attribute_$type WHERE".
          " graph_id = ? AND obj_id = ? AND key_id = ?",
          -level => 3);

    my $upd = $self->{STHS}{"UPDATE_ATTR_$type"} = $self->dbh->prepare
        ( -name => "Change $type attribute for object parameter",
          -sql  => "UPDATE attribute_$type SET val = ?".
          " WHERE graph_id = ? AND obj_id = ? AND key_id = ?",
          -level => 3);

    my $rv = 0;
    if (my $kid = $self->text2id($key)) {
        my $prior = $chk->get_single_value($gid, $oid, $kid);
        if (defined $prior) {
            # A value already exists.
            if ($prior eq $val) {
                # Requested value is already stored, nothing needs to be done
                $rv = -1;
            } else {
                # We need to change the value the current value.
                $upd->execute($val, $gid, $oid, $kid);
                $rv = $@ ? 0 : 1;
            }
        } else {
            # No values in DB, just insert
            $add->execute($gid, $oid, $kid, $val);
            $rv = $@ ? 0 : 1;
        }
    }
    return $rv;
}

sub read_attributes_by_name {
    my $self = shift;
    my ($gid, $name) = @_;
    my %rv;
    if ($gid && $name) {
        unless ($gid =~ /^\d+$/) {
            $gid = $self->text2id($gid);
        }
        foreach my $type (@attrTables) {
            my $stn  = "READ_ATTRNAME_$type";
            unless ($self->{STHS}{$stn}) {
                $self->{STHS}{$stn} = $self->dbh->prepare
                    ( -name => "Read $type attributes for object",
                      -sql  => "SELECT t2.name, a.val FROM attribute_$type a, text2id t1, text2id t2 WHERE a.graph_id = ? AND a.obj_id = t1.id AND t1.name = ? AND a.key_id = t2.id",
                      -level => 3);
            }
            my $rows = $self->{STHS}{$stn}->selectall_arrayref($gid, $name);
            map { $rv{ $_->[0] } = $_->[1] } @{$rows};
        }
    }
    return wantarray ? %rv : \%rv;
}

sub read_attributes {
    my $self = shift;
    my ($gid, $id) = @_;
    my %rv;
    if ($gid && $id) {
        my $sth = $self->{STHS}{"READ_FULL_ATTR"} ||= $self->dbh->prepare
            ( -name  => "Read all attributes for object",
              -level => 3,
              -sql   => "
SELECT t2.name, a.val
  FROM attribute_string a, text2id t2 
 WHERE a.graph_id = ? AND a.obj_id = ? AND t2.id = a.key_id
UNION
SELECT t2.name, a.val::Text
  FROM attribute_float a, text2id t2 
 WHERE a.graph_id = ? AND a.obj_id = ? AND t2.id = a.key_id
");

        # warn $sth->pretty_print($gid, $id, $gid, $id);
        my $rows = $sth->selectall_arrayref($gid, $id, $gid, $id);
        map { $rv{ $_->[0] } = $_->[1] } @{$rows};
    }
    return wantarray ? %rv : \%rv;
}

sub delete_attribute {
    my $self = shift;
    my ($gid, $oid, $key) = @_;
    return 0 unless ($gid && $oid && $key);
    my $kid = $self->text2id($key);
    return 0 unless ($kid);
    if ($self->is_locked) {
        $self->single_message("Can not delete_attribute - database is locked");
        return 0;
    }
    my $rv = 1;
    foreach my $type (@attrTables) {
        my $stn  = "DEL_ATTR_$type";
        unless ($self->{STHS}{$stn}) {
            $self->{STHS}{$stn} = $self->dbh->prepare
                ( -name => "Delete $type attribute from object",
                  -sql  => "DELETE FROM attribute_$type WHERE graph_id = ? AND obj_id = ? AND key_id = ?",
                  -level => 3);
        }
        $self->{STHS}{$stn}->execute($gid, $oid, $kid);
        if ($@) {
            $rv = 0;
        }
    }
    return $rv;
}

sub delete_attributes_deep {
    my $self = shift;
    my ($gid, $oid) = @_;
    return 0 unless ($gid && $oid);
    my $rv = 1;
    if ($self->is_locked) {
        $self->single_message("Can not delete_attributes_deep - database is locked");
        return 0;
    }
    foreach my $type (@attrTables) {
        my $stn  = "DEL_ALL_ATTR_$type";
        unless ($self->{STHS}{$stn}) {
            $self->{STHS}{$stn} = $self->dbh->prepare
                ( -name => "Delete all $type attributes from object",
                  -sql  => "DELETE FROM attribute_$type WHERE graph_id = ? AND obj_id = ?",
                  -level => 3);
        }
        $self->{STHS}{$stn}->execute($gid, $oid);
        $rv = 0 if ($@);
    }
    # die;
    return $rv;
}

sub create {
    my $self = shift;
    my $args = shift;
    my $dbn  = $args->{DATABASE} || 'friendlygraph';
    if ($self->is_locked) {
        $self->single_message("Can not create database - database is locked");
        return 0;
    }
    my $cmd  = "createdb $dbn 'Persistent storage of graph data for BMS::FriendlyGraph'";
    system($cmd);
    if (my $dbh = $self->connect($args)) {
        $dbh->schema( $self->schema );
        $dbh->make_all();
    } else {
        $fdbLastError .= "\nAn attempt was made to create the database, but this appears to have failed.\n";
    }
}

sub connect {
    my $self = shift;
    my $args   = shift;
    my $dbn    = $args->{DATABASE};
    my $host   = $args->{HOST};
    my $port   = $args->{PORT};
    my $user   = $args->{USER};
    my $erFile = $args->{ERRFILE};
    my $mail   = $args->{MAIL};
    my $dbh;
    eval {
        $dbh = BMS::FriendlyDBI->connect
            ("dbi:Pg:dbname=$dbn; host=$host; port=$port", 'tilfordc',
             undef, { PrintError => 0 },
             -noenv     => 1,
             -errorfile => $erFile,
             -adminmail => $mail,  );
    };
    return $dbh;
}

sub schema {
    my $self = shift;
    my %tables;
    $tables{"text2id"} = {
        name  => "text2id",
        com   => "Maps all string-based objects to integer IDs",
        sequence => { 
            't2i_seq' => 1,
        },
        index => {
            ti_n_ind   => {
                cols => [ 'name' ],
                unique => 1,
            },
            ti_upper_ind   => {
                cols => [ 'upper(name)' ],
            },
            ti_i_ind => { 
                cols => ['id'],
            },
        },
        cols  => [ ['id', 'integer', 
                   'The primary integer key assigned to the object. Generated automatically by the program'],
                  ['name', 'text',
                   'The text representation of the object. Supplied by the user.' ],
                   ],
    };
    $tables{"node"} = {
        name => 'node',
        com  => "Iterates all nodes in a graph.",
        index => {
            node_ng_ind   => {
                cols => [ 'node_id', 'graph_id' ],
                unique => 1,
            },
            node_g_ind   => {
                cols => [ 'graph_id' ],
            },
        },
        cols  => [ ['graph_id', 'integer', 
                   'The graph that the node is part of (FKEY to text2id)'],
                   ['node_id', 'integer', 
                   'The node itself (FKEY to text2id)'],
                   ],
    };
    $tables{"edge"} = {
        name => 'edge',
        com  => "Connects two nodes together within a specific graph",
        index => {
            edge_n1n2g_ind   => {
                cols => [ 'node1', 'node2', 'graph_id' ],
                unique => 1,
            },
            edge_n2n1_ind   => {
                cols => [ 'node2', 'graph_id' ],
            },
            edge_g_ind   => {
                cols => [ 'graph_id' ],
            },
            edge_e_ind   => {
                cols => [ 'edge_id' ],
            },
        },
        cols  => [ ['graph_id', 'integer', 
                    'The graph that the edge is part of (FKEY to text2id)'],
                   ['node1', 'integer', 
                    'The source node (FKEY to text2id)'],
                   ['node2', 'integer', 
                    'The target node (FKEY to text2id)'],
                   ['edge_id', 'integer', 
                    'An ID representing this edge (PKEY used for attributes)'],
                   ],
    };
    $tables{"attribute_string"} = {
        name => 'attribute_string',
        com  => "Records a string (text) attribute for an object",
        index => {
            as_ogk_ind   => {
                cols => [ 'obj_id','graph_id', 'key_id' ],
           },
            as_kg_ind => {
                cols => [ 'key_id', 'graph_id'  ],
            },
        },
        cols  => [ ['graph_id', 'integer', 
                   'The graph that the data is part of (FKEY to text2id)'],
                   ['obj_id', 'integer', 
                   'The identifier for the object being qualified (text2id)'],
                   ['key_id', 'integer', 
                    'The attribute key (name, via text2id)'],
                   ['val', 'text', 
                    'The value associated with the object, text format'],
                  ],
    };
    $tables{"attribute_float"} = {
        name => 'attribute_float',
        com  => "Records a floating point (numeric) attribute for an object",
        index => {
            as_ofk_ind   => {
                cols => [ 'obj_id','graph_id', 'key_id' ],
            },
            af_kg_ind => {
                cols => [ 'key_id', 'graph_id'  ],
            },
        },
        cols  => [ ['graph_id', 'integer', 
                   'The graph that the data is part of (FKEY to text2id)'],
                   ['obj_id', 'integer', 
                   'The identifier for the object being qualified (text2id)'],
                   ['key_id', 'integer', 
                    'The attribute key (name, via text2id)'],
                   ['val', 'float', 
                    'The value associated with the object, float format'],
                   ],
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
    $tables{ "v_ind" } =
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

    return \%tables
}

return 1;
