#!/usr/bin/perl -w

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

my $isBeta;

use strict;
use BMS::ArgumentParser;
use BMS::SnpTracker::MapLoc;
use Bio::SeqIO;

my ($mapLoc);

my $args = BMS::ArgumentParser->new
    ( -limit    => 0,
      -nocgi    => 1,
      -fork     => 1,
      -verbose  => 1,
      -progress => 300,
      );

my $limit    = $args->val(qw(limit));
my $mlInt     = $args->val(qw(instance));


&load_file($args->val(qw(input file)));


sub load_file {
    my $file = shift;
    return unless ($file);
    $args->msg("Reading Fasta file", $file);
    my $ml = &location_dbi();
    $args->death("Failed to locate file", $file) unless (-e $file);
    my $reader = Bio::SeqIO->new( -file   => "<$file" );
    my $num = 0;
    while ( my $bs = $reader->next_seq() ) {
        my $id = $bs->display_id();
        my ($spec, $type, $chr, $build);
        if ($id =~ /^([a-z_]+)\.([^\.]+)\.([^\.]+)\.([^\.]+)$/i) {
            ($spec, $type, $chr, $build) = ($1, $2, $3, $4);
        }
        next unless ($chr && $build);
        my $txt = BMS::SnpTracker::MapLoc::Text->new($ml, $id);
        $txt->tag( 'Build', $build);
        $txt->tag( 'Description', $bs->desc() );
        $txt->tag( 'Length', $bs->length() );
        $txt->tag( 'Name', $chr );
        $txt->tag( 'Assembly Type', $type );
        if ($spec) {
            $spec =~ s/_/ /g;
            substr($spec, 0, 1) = uc(substr($spec, 0, 1));
            $txt->tag( 'Species', $spec );
        }
        # warn $txt->to_text();
        $txt->write_tags();
        $ml->define_acc_for_loc( $id, $build, $chr );
        $num++;
        last if ($limit && $num >= $limit);
    }
    $args->msg("Loaded $num sequences");
}

sub location_dbi {
    return $mapLoc if ($mapLoc);
    $mapLoc ||= BMS::SnpTracker::MapLoc->new
        ( -build    => $args->val(qw(build)),
          -noenv    => $args->val(qw(noenvironment noenv)),
          -instance => $mlInt,
          -makedb   => $args->val(qw(makedatabase makedb)) );
    return $mapLoc;
}

