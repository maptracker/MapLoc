# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::SnpTracker::MapLoc::PopulationFilter;
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
use BMS::SnpTracker::MapLoc::Common;

use vars qw(@ISA);
@ISA      = qw(BMS::SnpTracker::MapLoc::Common);

our $handleCounter = 0;
sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self = {
        MINMAF   => 0,
        W        => 1,
        CATS     => ['Unspecified Variants'],
    };
    bless ($self, $class);
    $self->bench_start();
    my $args = $self->parseparams( @_ );

    unless ($self->{MAPLOC} = $args->{MAPLOC}) {
        $self->death("Failed to make population filter",
                     "MapLoc object must be provided with -maploc");
    }
    
    $self->_set_freq_func();
    $self->_set_imp_func();
    $self->bench_end();
    return $self;
}

sub obj_type   { return "PopFilter"; }

sub handle {
    return shift->{HANDLE} ||= "PopFilt-". ++ $handleCounter;
}

sub name {
    my $self = shift;
    if ($_[0]) { 
        # Setting the filter name
        $self->{NAME} = $_[0];
    }
    return $self->{NAME} || "Un-named Filter";
}

sub categories {
    my $self = shift;
    if (my $req = $_[0]) { 
        # Setting the categories
        $req = [$req] unless (ref($req));
        $self->{CATS} = $_[0];
    }
    return wantarray ? @{$self->{CATS}} : join(', ', @{$self->{CATS}});
}

sub weight {
    # Relative weight assigned to the population. Used when determining
    # outcome from two or more conflicting populations at one location
    my $self = shift;
    $self->{W} = $_[0] if (defined $_[0]);
    return $self->{W};
}

sub to_text {
    my $self = shift;
    my $txt = sprintf
        ("[%s] %s\n   Weight %s | MAF %s%% | Nulls %s\n   %s\n",
         $self->obj_type(), $self->name(),
         $self->weight(), 100 * ($self->min_maf() || 0), $self->null_ok() ?
         "allowed" : "rejected", $self->{IMPTEXT});
    my @cats = $self->categories;
    if (my $cnum = $#cats + 1) {
        $txt .= sprintf("   Categor%s: %s\n", $cnum == 1 ? 'y' : 'ies',
                        join(', ', @cats));
    }

    return $txt;
}

sub null_ok {
    # Flag to indicate if non-defined frequencies are considered
    # passing
    my $self = shift;
    if (defined $_[0]) {
        $self->{NULLOK} = $_[0] ? 1 : 0;
        # Recalculate the closure used to test frequencies:
        $self->_set_freq_func();
    }
    return $self->{NULLOK};
}

sub min_maf {
    # Minimum minor allele frequency
    my $self = shift;
    if (defined $_[0]) {
        my $minMaf = $_[0];
        if ($minMaf) {
            $self->{MINMAF} = $minMaf;
        } else {
            delete $self->{MINMAF};
        }
        # Recalculate the closure used to test frequencies:
        $self->_set_freq_func();
    }
    return $self->{MINMAF};
}

sub test_freq {
    # Test a single frequency value, return empty string if it passes filter,
    # otherwise a string describing the failure
    my $self = shift;
    return &{$self->{FREQFUNC}}( shift );
}

sub test_impact {
    # Test a single impact code, return empty string if it passes filter,
    # otherwise a string describing the failure
    my $self = shift;
    return &{$self->{STRICTFUNC}}( shift );
}

sub test_impact_loose {
    # Test a single impact code with loosened criteria, return empty
    # string if it passes, otherwise a string describing the failure
    my $self = shift;
    return &{$self->{LOOSEFUNC}}( shift );
}

sub test {
    # Test a frequency plus an impact code
    my $self = shift;
    # Test the frequency first - if it fails, fail immediately, return reason
    if (my $reason = $self->test_freq( shift )) { return $reason; }
    # Then test the impact:
    return $self->test_impact( shift );
}

sub test_loose {
    # Test a frequency plus an impact code, with loosened criteria
    my $self = shift;
    # Test the frequency first - if it fails, fail immediately, return reason
    if (my $reason = $self->test_freq( shift )) { return $reason; }
    # Then test the impact:
    return $self->test_impact_loose( shift );
}

sub _set_freq_func {
    my $self = shift;
    my $method;
    if (my $maf = $self->{MINMAF}) {
        # A minimum frequency is requested ...
        my $perc = 100 * $maf;
        if ($self->{NULLOK}) {
            # ... or a non-defined frequency
            $method = sub {
                my $f = shift;
                return (!defined $f || $f == -1) ? "" : 
                    $f < $maf ? "MAF less than $perc%" : "";
            };
        } else {
            # ... and null frequencies are not OK
            $method = sub {
                my $f = shift;
                return 
                    (!defined $f || $f == -1) ? "Frequency not defined" :
                    $f < $maf   ? "MAF less than $perc%"  : "";
            };
        }
    } elsif (!$self->{NULLOK}) {
        # No MAF, nulls not explicitly ok
        # Test is true as long as frequency is defined
        $method = sub {
            my $f = shift;
            return (!defined $f || $f == -1) ? "Frequency not defined" : "";
        };
    } else {
        # Anything goes!
        $method = sub { return ""; }
    }
    return $self->{FREQFUNC} = $method;
}

*keep_impact = \&require_impact;
sub require_impact {
    # Impacts we will require of the population in order to keep it
    my $self = shift;
    $self->_update_imp_hash( $_[0], 'KEEP' ) if (defined $_[0]);
    return wantarray ? sort keys %{$self->{KEEPIMP} || {}} : $self->{KEEPIMP};
}

*toss_impact = \&forbid_impact;
sub forbid_impact {
    # We will reject this population if it has these impacts
    my $self = shift;
    $self->_update_imp_hash( $_[0], 'TOSS' ) if (defined $_[0]);
    return wantarray ? sort keys %{$self->{TOSSIMP} || {}} : $self->{TOSSIMP};
}

sub _update_imp_hash {
    my $self = shift;
    $self->bench_start();
    my ($reqs, $key) = @_;
    if (!$reqs) {
        # Clear filter
        delete $self->{$key};
    } else {
        my @chk  = ref($reqs) ? @{$reqs} : split(/\s*[\n\r\,]+\s*/, $reqs);
        my $targ = $self->{$key.'IMP'} = {};
        my $ml   = $self->maploc();
        foreach my $req (@chk) {
            if (my $det = $ml->impact_details($req)) {
                $targ->{$det->{token}} = 1;
            } else {
                $self->msg_once("[!]", "Could not set $key impact filter for unrecognized impact '$req'");
            }
        }
        # Recalculate the closure used to test frequencies:
        $self->_set_imp_func();
    }
    $self->bench_end();
}

sub _set_imp_func {
    my $self = shift;
    my @keep = $self->require_impact();
    my @toss = $self->forbid_impact();
    my $ml   = $self->maploc();
    my @all  = $ml->impact_list();
    my @seed = $#keep == -1 ? @all : @keep;
    # warn $self->name()." Keep=$#seed, Toss=$#toss\n";
    foreach my $key ('STRICT', 'LOOSE') {
        # We will 'subtract' the unwanted codes from %full. We will
        # assume if the user wants to exclude 'COD', they wish to just
        # exclude the GENERIC impact code, and not the children. Such
        # filters are useful for excluding 'variants' for which no
        # alleles are provided.
        my %rem  = map { $_ => 1 } $ml->expand_impacts( \@toss, 0, 1 );
        my $kpm  = 0; # Keep-list Parent Mode
        if ($key eq 'LOOSE') {
            # We want to make sure high-level general impact codes are
            # kept if *any* of the more specific children are
            # requested:
            $kpm = 1;
            # We also want to still keep certain impact codes, even if
            # we have been requested to discard them. This is because
            # loose impact calculations will report specific impacts
            # generically
            map { delete $rem{$_} } qw(SPL RNA);
        }
        my %full = map { $_ => 1 } $ml->expand_impacts( \@seed, $kpm );
        # warn "Full = [".join(' ', sort keys %full)."] - [".join(' ', sort keys %rem)."] via{".join(" ", @toss)."}\n";
        map { delete $full{$_} } keys %rem;
        $self->{$key.'LIST'} = [ sort keys %full ];
    }
    my @strictImps = @{$self->{STRICTLIST}};
    my $strictSz   = $#strictImps;
    my $allSz      = $#all;
    if ($strictSz == -1) {
        # This filter removes everything
        $self->{IMPTEXT}    = "All impacts rejected!";
        $self->{STRICTFUNC} =  $self->{LOOSEFUNC} = 
            sub { return "Filter rejects every impact!"; };
        return;
    } elsif ($strictSz == $allSz) {
        # No impact criteria, always keep
        $self->{IMPTEXT}    = "All impacts kept";
        $self->{STRICTFUNC} =  $self->{LOOSEFUNC} = sub { return ""; };
        return;
    }
    my $tossReason;
    if ($strictSz < $allSz / 2) {
        # Easier to describe what we keep
        my $ltxt = join(' ', @strictImps);
        $self->{IMPTEXT} = "Keeping $ltxt";
        $tossReason = "Not one of [$ltxt]";
    } else {
        # Easier to describe what we discard
        my %ah = map { $_ => 1 } @all;
        map { delete $ah{$_} } @strictImps;
        my @disc = sort keys %ah;
        my $ltxt = join(' ', @disc);
        $self->{IMPTEXT} = "Discarding $ltxt";
        $tossReason = "One of [$ltxt]";
    }
    foreach my $key ('STRICT', 'LOOSE') {
        my %hash = map { $_ => 1 } @{$self->{$key.'LIST'}};
        # warn $self->name(). ":$key = ".join(' ',@{$self->{$key.'LIST'}})."\n";
        $self->{$key.'FUNC'} = sub {
            my $i = uc(shift || "");
            return $hash{$i} ? "" : $tossReason;
        };
    }
}
