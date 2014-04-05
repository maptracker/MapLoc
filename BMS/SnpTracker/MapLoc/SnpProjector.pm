# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::SnpTracker::MapLoc::SnpProjector;
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
use BMS::Utilities::Benchmark;
use BMS::Utilities::SequenceUtilities;
use BMS::SnpTracker::MapLoc;

use vars qw(@ISA);
@ISA      = qw(BMS::Utilities::SequenceUtilities
               BMS::Utilities::Benchmark);

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self = {
        PARAMS => {},
    };
    bless ($self, $class);
    my $args = $self->parseparams( @_ );
    $self->{BUILD}    = $args->{BUILD};
    $self->{INSTANCE} = $args->{INSTANCE};
    return $self;
}

sub maploc {
    my $self = shift;
    if (my $nv = shift) {
        $self->{MAPLOC} = $nv;
    }
    return $self->{MAPLOC} ||= BMS::SnpTracker::MapLoc->new
        ( -build => $self->{BUILD},
          -instance => $self->{INSTANCE}, );
}

sub param {
    my $self = shift;
    my $key  = uc(shift || "");
    my $val  = shift;
    if (defined $val) {
        if ($val eq '') {
            delete $self->{PARAMS}{$key};
        } else {
            $self->{PARAMS}{$key} = $val;
        }
    }
    return $self->{PARAMS}{$key};
}

sub param_or_argument {
    my $self = shift;
    my $key  = uc(shift || "");
    my $args = shift;
    return $args->{$key} if (defined $args->{$key});
    return $self->param($key);
}

sub ambiguity_for_accession {
    my $self = shift;
    my $args = $self->parseparams( @_ );
    my %rv;
    my $acc  = $args->{ID}  || $args->{ACC};
    return \%rv unless ($acc);
    $self->bench_start();
    $rv{acc} = $acc;
    my $src  = $args->{SRC} || $args->{SOURCE};
    my $ml   = $self->maploc();
    my $su   = $ml->sequtils();
    my @rnas = $ml->get_rnas( $acc, $src );
    if ($#rnas == -1) {
        $self->bench_end();
        return \%rv;
    } elsif ($#rnas != 0) {
        $self->bench_end();
        return \%rv;
    }
    my $rna     = $rnas[0];
    my $bld     = $self->{BUILD};
    my $minFrac = $rv{minFrac} = $self->param_or_argument('minfrac', $args);
    my @alns    = $bld ? $rna->each_alignment_for_build( $bld ):
        $rna->each_alignment();
    my $popids;
    if (my $preq = $self->param_or_argument('population', $args)) {
        my @reqs = ref($preq) ? @{$preq} : ($preq);
        foreach my $req (@reqs) {
            if (my $pop = $ml->get_population($req)) {
                $popids ||= $rv{pop} ||= [];
                push @{$popids}, $pop->pkey();
            }
        }
    }
    my %rcCache;
    foreach my $aln (@alns) {
        my $str  = $aln->strand() || 0;
        my @locs = $ml->locations_for_HSPs( [$aln], 0 );
        foreach my $loc (@locs) {
            # Only keep SNPs
            next unless ($loc->length == 1);
            my $locHsp   = $loc->hsp();
            my $anc      = $locHsp->anchor_id();
            my $int      = $aln->intersect
                ( -hsp => $locHsp, -anchor => $anc);
            next unless ($int && $int->length == 1);
            my ($l, $r) = $int->bounds_for_seq( $acc );
            my $pos = $l + 1;
            my $alls = $loc->each_allele();
            my @bases;
            foreach my $base (keys %{$alls}) {
                my $ub = uc($base);
                next unless ($ub =~ /^[A-Z]$/);
                my $fdat = $alls->{$base};
                my @pops = keys %{$fdat};
                if ($popids) {
                    my @keep;
                    foreach my $pid (@pops) {
                        push @keep, $pid if 
                            ($ml->popid_match($pid, @{$popids}));
                    }
                    @pops = @keep;
                }
                my @freqs;
                foreach my $pid (@pops) {
                    push @freqs, defined $fdat->{$pid}[0] ?
                        $fdat->{$pid}[0] : -1;
                }
                next if ($#freqs == -1);
                my ($max) = sort { $b <=> $a } @freqs;
                next if ($minFrac && $max < $minFrac);
                my $targ = $rv{pos}{$pos} ||= {};
                $ub = $rcCache{$ub} ||= $su->revcom($ub) if ($str < 0);
                $targ->{$ub} = $max if (!$targ->{$ub} ||
                                        $targ->{$ub} < $max);
            }
        }
    }
    $self->bench_end();
    return \%rv;
}


sub add_ambiguity_to_bioseq {
    my $self = shift;
    my $bs   = shift;
    return undef unless ($bs);
    $self->bench_start();
    my $acc;
    if ($bs->can('accession_number')) {
        $acc = $bs->accession_number();
        if (my $v = $bs->seq_version()) {
            $acc .= ".$v";
        }
    }
    $acc ||= $bs->display_id();
    my $ambig = $self->ambiguity_for_accession( @_, -acc => $acc );
    my @poses = keys %{$ambig->{pos}};
    my %rv;
    foreach my $key (keys %{$ambig}) {
        $rv{$key} = $ambig->{$key} unless ($key eq 'pos');
    }
    unless ($#poses == -1) {
        my $args = $self->parseparams( @_ );
        my $seq = $bs->seq();
        my $ml  = $self->maploc;
        my $su  = $ml->sequtils();
        my $altered = 0;
        foreach my $pos (@poses) {
            my $ref     = uc(substr($seq, $pos -1, 1));
            my $fdat    = $ambig->{pos}{$pos};
            $fdat->{$ref} = 1 if (! defined $fdat->{$ref});# || $fdat->{$ref} < 1);
            my @alleles = keys %{$fdat};
            # If there is only one allele here do nothing:
            next if ($#alleles < 1);
            my $amCode = $su->ambiguous_code( \@alleles );
            substr($seq, $pos -1, 1) = lc($amCode);
            $altered++;
            $rv{pos}{$pos} = $fdat;
        }
        if ($altered) {
            $bs->seq($seq);
            my $desc = $bs->desc() || "";
            my $note = sprintf
                ("%d ambiguit%s added from polymorphism",
                 $altered, $altered == 1 ? 'y' : 'ies');
            if (my $minFrac = $ambig->{minFrac}) {
                $note .= sprintf(" at %d+%%", $minFrac * 100);
            }
            if (my $pids = $ambig->{pop}) {
                my @pname;
                foreach my $pid (@{$pids}) {
                    if (my $pop = $ml->get_population($pid)) {
                        push @pname, $pop->name;
                    } else {
                        push @pname , "pop_id = $pid";
                    }
                }
                $note .= " in population ".join(' / ', @pname);
            }
            if ($desc) {
                $desc =~ s/\s+$//;
                $desc .= '.' unless ($desc =~ /\.$/);
                $desc .= ' ';
            }
            $desc .= $note;
            $bs->desc($desc);
        }
    }
    $self->bench_end();
    # my @dbp = sort { $a <=> $b } keys %{$rv{pos}}; my @bits; foreach my $p (@dbp) { push @bits, "$p : ".join(', ', map { "$_ = $rv{pos}{$p}{$_}"} sort keys %{$rv{pos}{$p}}) }; $self->msg("Surviving variants:", @bits);
    return \%rv;
}
