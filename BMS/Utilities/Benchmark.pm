# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::Utilities::Benchmark;
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
use Scalar::Util qw(weaken);
use Time::HiRes;

use vars qw(@ISA);
@ISA    = qw(BMS::Utilities);

sub hitime {
    return Time::HiRes::time;
    # If you don't have Time::HiRes, use:
    # return time;
}

# Cumulative benchmark information:
our $benchmarks      = {};
our $sharedInfo      = {};
# Stack for tracking nested benchmark calls:
our $startPoints     = [];
our $recursionTimes  = {};
our $progStartTime   = &hitime();

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self = {
    };
    bless ($self, $class);
    return $self;
}

sub clear_benchmarks {
    $benchmarks      = {};
    $sharedInfo      = {};
    $startPoints     = [];
    $recursionTimes  = {};
    $progStartTime   = &hitime();
}

# Setting a benchmark (pair of benchstart + benchend) costs about 35us
# The inner void time (the recorded time for a benchstart immediately
# followed by a benchend) for a benchmark is about 13us

*bench_start = \&benchstart;
sub benchstart {
    my $self    = shift;
    my $nowTime = &hitime();
    my $func    = $self->calling_package( @_ );
    return unless ($func);
    push @{$startPoints}, [ $func, $nowTime, 0 ];
    $recursionTimes->{$func}[0]++;
}

sub count_stack {
    my $self    = shift;
    if (my $func    = $self->calling_package( @_ )) {
        my $stack = $self->stack_trace();
        $self->{STACK_COUNT}{$func}{$stack}++;
    }
}

sub calling_package {
    my $self = shift;
    # Get the package making the call:
    my $func = (caller(2))[3];
    return unless ($func);
    # Optional sub-level grouping:
    $func .= "-$_[0]" if ($_[0]);
    return $func;
}

*bench_end = \&benchend;
*bench_stop = \&benchend;
*benchstop = \&benchend;
sub benchend {
    my $self = shift;
    my $nowTime   = &hitime();
    my $lastPoint = pop @{$startPoints};
    my $func      = $self->calling_package( @_ );
    return unless ($func);
    unless ($lastPoint) {
        # No start points on the stack
        $self->err("Error setting benchmark for $func",
                   "It looks like benchstart() was never called");
        return;
    }
    my ($funcCheck, $priorTime, $childTime) = @{$lastPoint};
    unless ($funcCheck eq $func) {
        # Hmm the checkpoint stack got messed up.
        # benchstart probably called without a benchend
        $self->err("Error setting benchmark for $func",
                   "It looks like $funcCheck never called benchend()");
        $startPoints = [];
        return;
    }
    # The entire time required for the function to return:
    my $elapsedTime = $nowTime - $priorTime;
    # The part of that time not used by child functions:
    my $selfTime    = $elapsedTime - $childTime;
    if (--$recursionTimes->{$func}[0]) {
        # This function was directly or indirectly called by itself.
        # We need to note how much time was spent recursively in the
        # function to accurately account for total time
        $recursionTimes->{$func}[1] += $elapsedTime;
    } else {
        # This is the first (or only) parent in the recursion chain
        $childTime -= $recursionTimes->{$func}[1] || 0;
        $recursionTimes->{$func}[1] = 0;
    }

    $benchmarks->{ $func }{func} ||= $func;
    $benchmarks->{ $func }{n}++;
    $benchmarks->{ $func }{t}    += $selfTime;
    $benchmarks->{ $func }{tt}   += $selfTime * $selfTime;
    $benchmarks->{ $func }{kids} += $childTime;

    $sharedInfo->{LASTBENCH} = $selfTime;
    $sharedInfo->{LASTFULL}  = $elapsedTime;

    if ($#{$startPoints} != -1) {
        # This benchmark is the 'child' of another - note the child's time
        # in the parent:
        $startPoints->[-1][2] += $elapsedTime;
    }
    if (1) {
        # Also record the amount of time spent benchmarking
        # For a two day load that called this function 28 million times
        # The average was 33 usec (about 0.6% of total time)
        my $bdat = $benchmarks->{'BMS::Utilities::Benchmark::benchend'} ||=
        { func => 'BMS::Utilities::Benchmark::benchend',
          kids => 0 };
        my $benchTime  = &hitime() - $nowTime;
        $bdat->{t}    += $benchTime;
        $bdat->{tt}   += $benchTime * $benchTime;
        $bdat->{n}++;
        # $bdat->{kids} = 0;
    }
}

sub lastbench {
    my $self = shift;
    my $key  = $_[0] ? 'LASTFULL' : 'LASTBENCH';
    return $sharedInfo->{$key};
}

our @textCols = qw(func n t tt kids);

sub benchmarks_to_text {
    my $self  = shift;
    my @funcs = sort { $b->{t} <=> $a->{t} } values %{$benchmarks};
    my $txt   = "";
    foreach my $data ( @funcs ) {
        $txt .= join("\t", map { $data->{$_} } @textCols)."\n";
    }
    return $txt;
}

sub benchmarks_from_text {
    my $self = shift;
    my ($text, $denom) = @_;
    return unless ($text);
    $denom ||=1;
    foreach my $line (split(/\n/, $text)) {
        my ($func, $n, $t, $tt, $kids) = split(/\t/, $line);
        next unless ($n && $n =~ /^\d+$/);
        $benchmarks->{ $func }{func} ||= $func;
        $benchmarks->{ $func }{n}    += $n;
        $benchmarks->{ $func }{t}    += $t / $denom;
        $benchmarks->{ $func }{tt}   += $tt / ($denom * 2);
        $benchmarks->{ $func }{kids} += $kids / $denom;
    }
}

sub benchmarks_to_file {
    my $self = shift;
    my $file = shift;
    if (!$file) {
    } elsif (open(BENCHFILE, ">$file")) {
        print BENCHFILE $self->benchmarks_to_text();
        close BENCHFILE;
    } else { 
        $self->err("Failed to write benchmark file", $file, $!);
    }
}

sub benchmarks_from_file {
    my $self = shift;
    my $file = shift;
    return unless ($file);
    if (open(BENCHFILE, "<$file")) {
        my $txt = "";
        while (<BENCHFILE>) {
            $txt .= $_;
        }
        close BENCHFILE;
        return $self->benchmarks_from_text( $txt, @_);
    } else { 
        $self->err("Failed to read benchmark file", $file, $!);
    }
}

*showbench = \&show_benchmarks;
*show_bench = \&show_benchmarks;
sub show_benchmarks {
    my $self  = shift;
    my $args  = $self->parseparams( @_ );
    my $isH   = $args->{HTML} || $args->{ASHTML};
    my @funcs = sort { $b->{t} <=> $a->{t} } values %{$benchmarks};
    my $total = 0;
    my $frac  = $args->{MINFRAC} || 0;
    my $elaps = &hitime();
    $elaps -= $progStartTime;
    map { $total += $_->{t} } @funcs;
    my $msg = "";
    my @timeUnit  = ('s','ms','us','ns');
    my @elapTime  = (['',1],['second',60],['minute',60*60],['hour',60*60*24],['day']);
    my @tSty      = ('font-weight:bold; color:red',
                     'color:#f90',
                     'color:#39f',
                     'color:#c9f' );
    if ($isH) {
        $timeUnit[0]  = 'sec';
        $timeUnit[2]  = '&micro;s';
        $msg .= "<table";
        if (my $c = $args->{CLASS}) { $msg .= " class='$c'"; }
        $msg .= "><caption>";
    }
    $msg .= "Benchmarks";
    if ($frac && $frac =~ /^0?\.\d+$/ && $frac > 0 && $frac < 1) {
        $frac *= 100;
        $msg .= sprintf(" (showing %.2f%% or larger)", $frac);
    } elsif ($frac) {
        $msg .= " (ignoring unusual -minfrac '$frac')";
        $frac = 0;
    }
    if ($isH) {
        $msg .= "<br />\n";
    } else {
        $msg .= ":\n";
    }
    my $st = 1;
    while ($st < $#elapTime && $elaps / $elapTime[$st][1] > 1.5) { $st++; }
    my ($eUnit, $eMod) = ($elapTime[$st][0], $elapTime[$st-1][1]);
    $msg .= sprintf("Elapsed time %.3f %ss, %.3f are benchmarked\n", 
                    $elaps / $eMod, $eUnit, $total / $eMod);
    my @head = ('%Elaps','w/Kids', 'N','Avg','Avg+Kids','Method');
    my @cols = (6, 6, 5, 9, 9, -40);
    my @frm  = map {'%'.$_.'s' } @cols;
    my $frm = "  " .join(' ', map { '%'.$_.'s' } @cols) ."\n";
    if ($isH) {
        $msg .= "</caption><tbody>\n";
        $msg .= "<tr><th>".join("</th><th>", @head)."</th></tr>\n";
    } else {
        $msg .= sprintf($frm, @head);
        $msg .= sprintf($frm, map { '-' x abs($_) } @cols);
    }
    my @countUnit = ('','k','M','G');
    foreach my $data ( @funcs ) {
        my $t  = $data->{t};
        my $n  = $data->{n};
        my $nm = $data->{func};
        my $ft = $t + $data->{kids};
        my ($f, $of) = map { 100 * $_ / $elaps } ($t, $ft);
        # Technically we should just be able to check $of, but I do not
        # fully trust the child summation code
        next if ($frac && $f < $frac && $of < $frac);
        my $sf  = sprintf("%.1f%%", $f);
        my $sof = sprintf("%.1f%%", $of);

        my $av = $t / $n;
        my $ai = 0;
        while ($av < 1 && $ai < $#timeUnit) {
            $av *= 1000;
            $ai++;
        }
        my $tu  = $timeUnit[$ai];
        $tu = sprintf("<span style='%s'>%s</span>", $tSty[$ai], $tu) if ($isH);
        my $at = sprintf("%.2f%s", $av, $tu);

        my $akv = $ft / $n;
        my $aki = 0;
        while ($akv < 1 && $aki < $#timeUnit) {
            $akv *= 1000;
            $aki++;
        }
        $tu  = $timeUnit[$aki];
        $tu = sprintf("<span style='%s'>%s</span>", $tSty[$aki],$tu) if ($isH);
        my $akt = sprintf("%.2f%s", $akv, $tu);

        my $ni = 0;
        while ($n >= 10000 && $ni < $#countUnit) {
            $n /= 1000;
            $ni++;
        }
        my $nt = sprintf("%d%s", int(0.5 + $n), $countUnit[$ni]);
        my @row = ($sf, $sof, $nt, $at, $akt, $nm);
        if ($isH) {
            my ($fred, $ofred) = map { 255 - int(0.5 + 255 * $_ /100) } ($f, $of);
            my @sty = map {
                " style='background-color:rgb(255,$_,$_)'" } ($fred, $ofred);
            $sty[2] = " style='text-align:center'";
            $sty[$#row] = $sty[1];
            $msg .= " <tr>".join('', map { 
                "<td".($sty[$_] || "").">$row[$_]</td>"
                } (0..$#row)). "</tr>\n";
        } else {
            my @txt = map { sprintf($frm[$_] || "%s", $row[$_]) } (0..$#row);
            if ($args->{SHELL}) {
                # Apply shell coloring
                my @chk = ($f, $of);
                for my $i (0..1) {
                    my $cf;
                    if ($chk[$i] > 20) {
                        $cf = "\033[31m\%s\033[0;49m";
                    } elsif ($chk[$i] > 10) {
                        $cf = "\033[33m\%s\033[0;49m";
                    }
                    if ($cf) {
                        $txt[$i] = sprintf($cf, $txt[$i]);
                        $txt[5]  = sprintf($cf, $txt[5]) if ($i);
                    }
                }
                my $Nf;
                if ($txt[2] =~ /[MG]$/) {
                    $Nf = "\033[31m\%s\033[0;49m";
                } elsif ($txt[2] =~ /[k]$/) {
                    $Nf = "\033[33m\%s\033[0;49m";
                }
                $txt[2] = sprintf($Nf, $txt[2]) if ($Nf);
            }
            $msg .= "  ".join(' ', @txt)."\n";
        }
    }
    $msg .= "</tbody></table>\n" if ($isH);
    if (my $cs = $self->{STACK_COUNT}) {
        if ($isH) {
            $msg .= "<h3>Function entry path tallies:</h3>\n";
        } else {
            $msg .= "Function entry path tallies:\n";
        }
        foreach my $func (sort keys %{$cs}) {
            my $csf = $cs->{$func};
            my @points = sort { $csf->{$b} <=> $csf->{$a}} keys %{$csf};
            my $pn = $#points + 1;
            my $f  = sprintf("%s() = %d entry point%s", $func, $pn, $pn == 1 ? '' : 's');
            if ($isH) {
                $msg .= "<h4>$f</h4>\n";
            } else {
                $msg .= "$f:\n";
            }
            foreach my $stack (@points) {
                my $num = $csf->{$stack};
                $num = sprintf("%d event%s", $num, $num == 1 ? '' : 's');
                if ($isH) {
                    $msg .= "<div><b>$num</b><pre>$stack</pre></div>\n";
                } else {
                    $msg .= "[$num]\n$stack";
                }
            }
        }
    }
    return $msg;
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

return 1;
