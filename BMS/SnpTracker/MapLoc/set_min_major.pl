#!/usr/bin/perl -w

my $isBeta;

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
use BMS::ArgumentParser;
use BMS::SnpTracker::MapLoc;
use BMS::ForkCritter;

my ($mapLoc, $fc, $findFreq, $makeRow, $setRow);

my $args = BMS::ArgumentParser->new
    ( -limit    => 0,
      -nocgi    => 1,
      -fork     => 1,
      -verbose  => 1,
      -base     => "MapLoc-MinMajor",
      -ageall   => '17 sep 2011',
      -progress => 300,
      );

my $vb       = $args->val(qw(vb verbose)) || 0;
my $basefile = $args->val(qw(basefile base)) || "MLMM";

my $forkNum  = $args->val(qw(forknum fork)) || 1;
my $prog     = $args->val(qw(prog progress)) || 300;
my $limit    = $args->val(qw(limit));

my $cage     = $args->val(qw(cloudage ageall));
my $age      = $args->val(qw(age ageall));
my $newOnly  = $args->val(qw(new newonly));
my $clobber  = $newOnly || $args->val(qw(clobber));

&main;

sub main {
    &set_min_major();
}

sub set_min_major {
    undef $mapLoc;
    my @input;
    if (my $req = $args->val(qw(id ids))) {
        @input = ( -input => [map { [$_] } split(/\s*[\n\r\t\,]\s*/, $req)],
                   -input_type => 'array' );
    } else {
        my $idFile = $args->val(qw(idfile)) || &find_loc_ids();
        @input = ( -input       => $idFile,
                   -input_type  => 'tsv', );

    }
    $fc = BMS::ForkCritter->new
        ( -init_meth   => \&location_dbi,
          -finish_meth => \&finish,
          -method      => \&set_freq,
          @input,
          -limit       => $limit,
          -progress    => $prog,
          -verbose     => $vb );
    my $errFile = $basefile."-Errors.err";
    $fc->output_file('Errors', $errFile);
    if (my $failed = $fc->execute( $forkNum )) {
        $args->err("$failed children failed to finish execution");
    }
}

sub set_freq {
    my ($lid) = @{$_[0]};
    my %freqs;
    $findFreq->execute($lid);
    my $found = $findFreq->fetchall_arrayref();
    map { push @{$freqs{$_->[0]}}, $_->[1] } @{$found};
    # Get the major frequency for each source
    my @major;
    while (my ($sid, $fArr) = each %freqs) {
        my ($maj) = sort { $b <=> $a } @{$fArr};
        if ($#{$fArr} < 1) {
            # There is only one allele
            unless ($maj == 1) {
                # Ooo. It's not 100%.
                &fork_err("Single allele less than 100%", "loc_id = $lid AND pop_id = $sid");
                next;
            }
        }
        push @major, [$sid, $maj];
    }
    # Find the source with the smallest major frequency
    my ($best) = sort { $a->[1] <=> $b->[1] } @major;
    # warn $args->branch({found => $found, freq => \%freqs}) if ($best);
    $makeRow->execute($lid);
    $setRow->execute(@{$best || [undef, undef]}, $lid);
}


sub fork_err {
    return if ($#_ == -1);
    if ($_[0] =~ /^\[/) {
        # Also report error to screen
        $args->msg(shift @_, @_);
    }
    $fc->write_output('Errors', join("\t",@_)."\n");
}


sub find_loc_ids {
    my $idFile = $basefile."-LocIDs";
    $idFile   .= "-New" if ($newOnly);
    $idFile   .= ".txt";
    my @ids;
    unless (!$clobber && -s $idFile) {
        # We need to make a new file
        $args->msg("Generating Location ID File", $idFile);
        my $sql = "SELECT l.loc_id FROM location l";
        $sql   .= " WHERE NOT EXISTS (SELECT mm.loc_id FROM min_major_allele mm WHERE mm.loc_id = l.loc_id)" if ($newOnly);
        my $ml  = &location_dbi();
        my $get = $ml->dbh->prepare
            ( -name => "Get all location IDs",
              -sql  => $sql, );
        $get->execute();
        @ids = $get->get_array_for_field();
        open(IDFILE, ">$idFile") || $args->death
            ("Failed to write ID file", $idFile, $!);
        print IDFILE join("\n", @ids)."\n";
        close IDFILE;
    }
    return $idFile;
}

sub location_dbi {
    return $mapLoc if ($mapLoc);
    $mapLoc ||= BMS::SnpTracker::MapLoc->new
        ( -build => $args->val(qw(build)),
          -noenv  => $args->val(qw(noenvironment noenv)),
          -makedb => $args->val(qw(makedatabase makedb)) );

    $findFreq = $mapLoc->dbh->prepare
        ( -name => "Find frequencies for location",
          -sql  => "SELECT pop_id, freq FROM allele WHERE loc_id = ? AND freq IS NOT NULL" );
    $makeRow = $mapLoc->dbh->prepare
        ( -name => "Make frequency row",
          -sql  => "INSERT INTO min_major_allele (loc_id) VALUES (?)",
          -ignore => "duplicate key value" );
    $setRow = $mapLoc->dbh->prepare
        ( -name => "Set frequency value",
          -sql  => "UPDATE min_major_allele SET pop_id = ?, freq = ? WHERE loc_id = ?");
    return $mapLoc;
}

sub finish {
}
