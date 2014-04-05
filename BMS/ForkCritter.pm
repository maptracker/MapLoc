# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::ForkCritter;
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
use Bio::SeqIO;
use BMS::FriendlySAX;
use POSIX;
use BMS::ErrorInterceptor;
use BMS::Utilities::Benchmark;

use vars qw(@ISA);
@ISA   = qw(BMS::ErrorInterceptor BMS::Utilities::Benchmark);

my $VERSION = ' $Id: ForkCritter.pm,v 1.37 2014/03/20 18:19:19 tilfordc Exp $ ';

our $scriptStartTime = time;
our $sortTmpDir      = "/tmp/";

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self = {
        METHOD   => \&_bogus_method,
        PIDS     => [], # pid array of child process IDs
        PIDNUM   => {}, # keyed to pid, points to child number
        INPUT    => '',
        INPUTS   => [],
        TYPE     => '', # The type of the input file (sequence, csv, etc)
        INPTARGS => [], # Arguments to pass to the input method (used by SAX)
        LASTREC  => '', # The last item recovered by next_task
        RECORD   => 0, # Total number of records observed
        COUNT    => 0, # Number of records accepted by this child
        LIMIT    => 0, # Limit the number of records read
        PROG     => 0, # Progress notification interval in seconds
        DOALL    => 0, # Flag to parse all records
        REFORK   => 0, # Child number to allow re-running of failed child
        REFORKOK => 0, # Can program re-try crashed children?

        FHNUM    => 0, # Increment for file handles

        CHILD    => 0, # The child count for this process
        LASTKID  => 0, # The last child spawned
        LOWKID   => 1, # The lowest child count for the set
        FORKNUM  => 1, # The total number of children being spawned
        COUNTTAG  => '',      # What are we using to measure progress?
        COUNTWHAT => 'lines', # What are the things we are counting?
        HISTORY   => {},      # Number of records previously read
        VERBOSE   => 0,
        EXITCODE  => 0,
        TERMINATED => 0,
        FORCE_EXIT => 0,
    };
    bless ($self, $class);
    my $args = $self->parseparams
        ( -method => undef,
          -type   => 'ARRAY',
          @_ );

    while (my ($key, $val) = each %{$args}) {
        # Standardize arguments to no underscore
        next unless ($key =~ /_/);
        $key =~ s/_//g;
        $args->{$key} = $val;
    }

    $self->intercept_errors();
    $self->ignore_error("Bio/Root/IO.pm line 505");
    $self->method($args->{METHOD});
    $self->unstable_input( $args->{UNSTABLE} );
    $self->init_method( $args->{INITMETH}  || $args->{INITMETHOD} );
    $self->finish_method( $args->{FINISHMETH}  || $args->{FINISHMETHOD} ||
                          $args->{FINALMETH}  || $args->{FINALMETHOD} );
    $self->group_method($args->{GROUPMETH} || $args->{GROUPMETHOD});
    $self->skip_record_method( $args->{SKIPMETH}  || $args->{SKIP_METH} || 
                               $args->{SKIPMETHOD} );
    $self->last_item_method( $args->{LAST_ITEM_METH} || 
                             $args->{LAST_ITEM_METHOD} );
    $self->input_type(   $args->{INPUTTYPE} || $args->{TYPE} );
    $self->input_args(   $args->{INPUTARGS} || $args->{INPUT_ARGS} );
    $self->exit_code(    $args->{EXITCODE}  || $args->{EXIT_CODE} );
    $self->input(        $args->{INPUT} );
    $self->progress(     $args->{PROGRESS}  || $args->{PROG});
    $self->verbose(      $args->{VERBOSE}   || $args->{DEBUG});
    $self->limit(        $args->{LIMIT});
    $self->allow_refork( $args->{REFORK} );
    $self->set_colmap($args->{COLMAP} || $args->{REMAP});
    foreach my $key (qw(QUIETLIMIT LIMITMSG)) {
        $self->{uc($key)} = $args->{uc($key)};
    }
    # Initialize the task directory:
    $self->_task_directory();
    return $self;
}

sub reset {
    my $self = shift;
    # Reset the number for the 'first' kid
    $self->close_all();
    $self->{LOWKID}      = $self->fork_count + 1;
    $self->{FILES}       = {};
    $self->{PIDS}        = [];
    $self->{PIDNUM}      = {};
    $self->{REFORK}      = 0;
    $self->{INPUTS}      = [];
    $self->{FORCE_EXIT}  = 0;
    $self->{TERMINATED}  = 0;
    $self->{INITIALIZED} = 0; # bad idea? trying to get handles working

    map { $self->{$_} = undef; } qw(TOTAL_COUNT TASKDIR COL_MAP IOSEP IOSTRIP);

    foreach my $extra (map { uc($_) } @_) {
        delete $self->{ $extra };
    }
    $self->_clear_task;
    $scriptStartTime = time;
}

sub set_colmap {
    my $self = shift;
    my $cm = shift || {};
    while (my ($in, $out) = each %{$cm}) {
        $self->remap_header_name( $in, $out );
    }
}

sub set_column_separator {
    my $self = shift;
    my ($nv, $strip) = @_;
    if (defined $nv) {
        if ($nv eq '') {
            # Clear setting
            delete $self->{IOSEP};
            delete $self->{IOSTRIP};
        } else {
            $self->{IOSEP}   = $nv;
            $self->{IOSTRIP} = $strip if (defined $strip);            
        }
    }
    return wantarray ? ($self->{IOSEP}, $self->{IOSTRIP} ) : $self->{IOSEP};
}

sub remap_header_name {
    my $self = shift;
    my ($nameReq, $newval) = @_;
    my $name = $nameReq;
    $name = '' unless (defined $name);
    $name = uc($name); # unless ($self->column_case_matters());
    $name =~ s/[_\-\s]+/ /g;# unless ($self->column_whitespace_matters());
    if (defined $newval) {
        if ($newval ne '') {
            $self->{COL_MAP}{$name} = $newval;
        } else {
            delete $self->{COL_MAP}{$name};
        }
    }
    my $rv = $self->{COL_MAP}{$name};
    return defined $rv ? $rv : $nameReq;
}

sub header {
    my $self = shift;
    return @{$self->{SUBTYPES}{HEADER} || []};
}

sub DESTROY {
    my $self = shift;
    if (my $child = $self->child) {
        # This is a child
        $self->err("[+]", "Child $child DESTROY - " . &nice_date())
            if ($self->verbose > 2);
        $self->_finish;
        $self->err("[+]", "Child $child finished - " . &nice_date())
            if ($self->verbose > 1);
    } else {
        # This is the parent
        $self->_clear_task;
        my $dir = $self->_task_directory;
        if ($dir && -d $dir) {
            #print ">>>>CLEARING $dir\n". $self->stack_trace()."\n";
            delete $self->{TASKDIR};
            rmdir($dir);
        }
        $self->msg("[+]", "Parent ForkCritter finished - " . &nice_date())
            if ($self->verbose);
    }
    if (my $ec = $self->{EXITCODE}) { exit($ec) }
}

sub graceful_death {
    my $self = shift;
    my @msg  = $#_ == -1 ? ("ForkCritter fatal error!") : @_;
    if (my $child = $self->child) {
        push @msg, "Child $child";
    } else {
        push @msg, "Parent object";
    }
    if (my $li = $self->last_item) {
        unshift @msg, "Last Item: $li";
    }
    $self->err( @msg );
    $self->graceful_shutdown();
    die;
}

sub graceful_shutdown {
    my $self = shift;
    return unless ($self);
    $self->close_all();
    
}

sub close_all {
    my $self = shift;
    return unless ($self);
    while (my ($tag, $data) = each %{$self->{FILES} || {}}) {
        my $fh = $data->{FH};
        close $fh if ($fh);
        delete $data->{FH};
        my $file = $data->{FRAG};
        next unless ($file);
        if (-e $file && -s $file) {
            # The file exists and has data
            if ($data->{OPTS} =~ /sort(\d+)?/) {
                # The user wants the file sorted
                my $mem = $1 || 0;
                my $cmd = "sort ";
                $cmd .= "-S ${mem}G " if ($mem);
                $cmd .= "-T $sortTmpDir " if ($sortTmpDir);
                $cmd .= " $file -o $file";
                system($cmd);
            }
        } else {
            # The file exists but is empty
            unlink($file);
        }
    }
    if (my $io = $self->{SEQIO}) {
        $io->close();
    }
}

sub method {
    my $self = shift;
    if (my $meth = $_[0]) {
        if (ref($meth) eq 'CODE') {
            $self->{METHOD} = $meth;
        } else {
            $self->err("'$meth' is not a code reference",
                       "method() NOT set");
            $self->{METHOD} = \&_bogus_method;
        }
    }
    return $self->{METHOD};
}

sub init_method {
    my $self = shift;
    if (defined $_[0]) {
        my $meth = $_[0];
        if (!$meth) {
            $self->{INIT_METHOD} = undef;
        } elsif (ref($meth) eq 'CODE') {
            $self->{INIT_METHOD} = $meth;
        } else {
            $self->err("'$meth' is not a code reference",
                       "init_method() NOT set");
        }
    }
    return $self->{INIT_METHOD};
}

sub finish_method {
    my $self = shift;
    if (defined $_[0]) {
        my $meth = $_[0];
        if (!$meth) {
            $self->{FINISH_METHOD} = undef;
        } elsif (ref($meth) eq 'CODE') {
            $self->{FINISH_METHOD} = $meth;
        } else {
            $self->err("'$meth' is not a code reference",
                       "finish_method() NOT set");
        }
    }
    return $self->{FINISH_METHOD};
}

sub next_record_method {
    my $self = shift;
    if (my $meth = $_[0]) {
        if (ref($meth) eq 'CODE') {
            $self->{NEXTFUNC} = $meth;
        } else {
            $self->err("'$meth' is not a code reference",
                       "next_record_method() NOT set");
        }
    }
    return $self->{NEXTFUNC};
}

sub skip_record_method {
    my $self = shift;
    if (defined $_[0]) {
        my $meth = $_[0];
        if (!$meth) {
            $self->{SKIPFUNC} = undef;
        } elsif (ref($meth) eq 'CODE') {
            $self->{SKIPFUNC} = $meth;
        } else {
            $self->err("'$meth' is not a code reference",
                       "skip_record_method() NOT set");
        }
    }
    return $self->{SKIPFUNC};
}

=head2 group_method

 Title   : group_method
 Usage   : my $codeRef = $fc->group_method( $newval )
 Function: Sets / gets a code callback for grouping records
 Returns : The current callback
 Args    : Optional new value, should be a code reference

A group method is optional logic that will collect multiple records
together. This is useful when your data source should be analyzed in
logical sets, rather than one by one. The method will ONLY work if all
members of a group are in uninterupted sequential order.

The code reference will be provided with a single argument,
representing the currently encountered record. It is up to your code
to maintain the growing 'current' group, and to determine when a group
is complete. When you know you have a complete group (presumably
because the passed record is NOT in the growing group), return the
group object. Otherwise, return a false value (eg undef).

If the passed argument is itself undef, it indicates that ForkCritter
has reached the end of the data stream. Typically, your method will
then still have the 'current' group, which it should return at this
time (otherwise return undef).

=cut

sub group_method {
    my $self = shift;
    if (my $meth = $_[0]) {
        if (ref($meth) eq 'CODE') {
            $self->{GROUPFUNC} = $meth;
        } else {
            $self->err("'$meth' is not a code reference",
                       "group_method() NOT set");
        }
    }
    return $self->{GROUPFUNC};
}

sub count_method {
    my $self = shift;
    if (my $meth = $_[0]) {
        if (ref($meth) eq 'CODE') {
            $self->{COUNTFUNC} = $meth;
            $self->{MANUALCOUNTFUNC} = 1;
        } else {
            $self->err("'$meth' is not a code reference",
                       "count_method() NOT set");
        }
    }
    return $self->{COUNTFUNC};
}

sub input {
    my $self = shift;
    if (my $input = $_[0]) {
        my $type  = $self->input_type( $_[1] );
        $self->{INPUT}       = $input;
        $self->{TYPE}        = $type;
        $self->{TOTAL_COUNT} = undef;
        push @{$self->{INPUTS}}, [$input, $type];
    }
    return $self->{INPUT};
}

my $typeCounters = {
    basic   => \&_count_basic,
    array   => \&_count_array,
    fasta   => \&_count_fasta,
    fastq   => \&_count_fastq,
    genbank => \&_count_genbank,
    xml     => \&_count_xml,
    tsv     => \&_count_tab,
    maf     => \&_count_tab,
    csv     => \&_count_tab,
};
*type = \&input_type;
sub input_type {
    my $self = shift;
    if (my $req = lc($_[0] || "")) {
        my $stnd = $self->{TYPE} = $self->_standard_types($req);
        my $st   = $self->{SUBTYPES} = {};
        if ($req =~ /(head|maf)/) {
            $st->{HEADER} = 1;
            if ($req =~ /hash/) {
                $st->{HASH} = 1;
            }
            if ($req =~ /maf/) {
                $st->{LOOKFOR} = 'NCBI';
                $st->{MAF}     = 1;
            }
        }
        $st->{CLEANASCII} = 1 if ($req =~ /clean/);
        if ($req =~ /lookfor_(\S+)/) {
            $st->{LOOKFOR} = $1;
        }
        if ($req =~ /groupby_(\S+)/) {
            $st->{GROUPBY} = $1;
        }
        if ($stnd eq 'seq') {
            if ($req =~ /(gbk|genbank|gb)/) {
                $st->{SEQFORMAT} = 'genbank';
            } elsif ($req =~ /fastq|fq/) {
                $st->{SEQFORMAT} = 'fastq';
            } elsif ($req =~ /uniprot|swiss|sp/) {
                $st->{SEQFORMAT} = 'swiss';
            } elsif ($req =~ /ipi/) {
                $st->{SEQFORMAT} = 'ipi';
            } else {
                $st->{SEQFORMAT} = 'fasta';
            }
            $stnd = $st->{SEQFORMAT};
            if ($req =~ /pair/) {
                $st->{PAIRED} = 1;
            }
        }
        unless ($self->{MANUALCOUNTFUNC}) {
            $self->{COUNTFUNC} = $typeCounters->{$stnd};
        }
    }
    return $self->{TYPE};
}

sub _standard_types {
    my $self = shift;
    if (my $type  = lc(shift || "")) {
        if ($type =~ /(seq|fasta|gbk|genbank|uniprot|swiss|sp)/) {
            return 'seq';
        } elsif ($type =~ /basic/) {
            return 'basic';
        } elsif ($type =~ /user/) {
            return 'user';
        } elsif ($type =~ /arr/) {
            return 'array';
        } elsif ($type =~ /tab/ || $type =~ /[ct]sv/) {
            return $type =~ /csv/ ? 'csv' : 'tsv';
        } elsif ($type =~ /maf/) {
            return 'maf';
        } elsif ($type =~ /sax/ || $type =~ /xml/) {
            return 'xml';
        }
        warn "Unrecognized input type '$type'";
        return $type;
    }
    return "";
}

sub input_args {
    my $self = shift;
    if ($_[0]) {
        $self->{INPTARGS} = $_[0];
    }
    return $self->{INPTARGS};
}

sub unstable_input {
    my $self = shift;
    if (defined $_[0]) {
        $self->{UNSTABLE_INPUT} = $_[0] ? 1 : 0;
    }
    return $self->{UNSTABLE_INPUT};
}

sub allow_refork {
    my $self = shift;
    if (defined $_[0]) {
        $self->{REFORKOK} = $_[0] ? 1 : 0;
    }
    return $self->{REFORKOK};
}

sub force_exit {
    return shift->{FORCE_EXIT}++;
}

sub exit_code {
    my $self = shift;
    if (defined $_[0]) {
        $self->{EXITCODE} = $_[0];
    }
    return $self->{EXITCODE};
}

sub progress {
    my $self = shift;
    if (my $num = $_[0]) {
        if ($num =~ /^\d+$/) {
            $self->{PROG} = $num;
        } else {
            $self->err("prog() requires an integer argument, not '$num'");
        }
    }
    return $self->{PROG};
}

sub last_item_method {
    # Allows customization of the object stored by last_item()
    my $self = shift;
    if (defined $_[0]) {
        $self->{LAST_ITEM_METH} = $_[0] ? $_[0] : undef;
    }
    return $self->{LAST_ITEM_METH};
}

sub last_item {
    # Sets / Gets the last item recovered for this fork
    my $self = shift;
    if ($_[0]) {
        if (my $meth = $self->{LAST_ITEM_METH}) {
            $self->{LASTREC} = &{$meth}($_[0]) || '';
        } else {
            $self->{LASTREC} = $_[0];
        }
    }
    return $self->{LASTREC};
}

sub total_fork {
    # Set via execute(), or manually by user
    # The total number of child processes spawned in a batch
    my $self = shift;
    if (my $num = $_[0]) {
        if ($num =~ /^\d+$/) {
            $self->{FORKNUM} = $num;
        } else {
            $self->err("total_fork() requires an integer argument, not '$num'");
        }
    }
    return $self->{FORKNUM};
}

sub fork_count {
    # Always incremented via fork()
    # The total number of children forked so far
    my $self = shift;
    return $self->{LASTKID};
}

sub child {
    # The child number of the current batch
    my $self = shift;
    return $self->{CHILD};
}

sub match_modulus {
    my $self  = shift;
    return 1 if ($self->doall);
    my ($num) = @_;
    # Child will be a value 0 to (totalFork -1)
    my $child = $self->child - $self->{LOWKID};
    my $mod   = ($num - 1) % $self->total_fork;
    return ($child == $mod) ? 1 : 0;
}

sub match_block {
    my $self  = shift;
    my ($num) = @_;
    # Child will be a value 0 to (totalFork -1)
    my $child = $self->child - $self->{LOWKID};
    
}

sub doall {
    my $self = shift;
    if (defined $_[0]) {
        $self->{DOALL} = $_[0] ? 1 : 0;
    }
    return $self->{DOALL};
}

sub limit {
    my $self = shift;
    if (defined $_[0]) {
        my $num = $_[0];
        if ($num =~ /^\d+$/) {
            $self->{LIMIT} = $num;
        } else {
            $self->err("limit() requires an integer argument, not '$num'");
        }
    }
    return $self->{LIMIT};
}

sub extend_limit {
    my $self = shift;
    my $amnt = shift || 1;
    if ($self->{LIMIT}) {
        my $nl = $self->{LIMIT} += $amnt;
        return $nl;
    }
    return 0;
}

sub verbose {
    my $self = shift;
    if (defined $_[0]) {
        $self->{VERBOSE} = $_[0];
    }
    return $self->{VERBOSE};
}

sub history {
    my $self = shift;
    my $input = $self->input;
    return 0 unless ($input);
    return $self->{HISTORY}{$input} ? $self->{HISTORY}{$input} : 0;
}

sub history_file {
    my $self = shift;
    if (my $num = $self->child) {
        return 0 unless ($self->{HFILE});
        return $self->{HFILE} . "_$num";
    }
    my ($file) = @_;
    if ($file) {
        $self->{HFILE}   = $file;
        $self->{HISTORY} = {};
        if (-e $file) {
            # Read in prior history
            open(HIST, "<$file") || $self->
                graceful_death("Failed to read history file", $file, $!);
            my %nums;
            while (<HIST>) {
                chomp;
                my ($hfile, $recnum, $child) = split("\t", $_);
                $nums{$hfile}{$recnum}++;
            }
            close HIST;
            foreach my $hfile (keys %nums) {
                my ($lowest) = sort { $a <=> $b } keys %{$nums{$hfile}};
                $self->{HISTORY}{$hfile} = $lowest;
            }
        }
    }
    return $self->{HFILE};    
}

sub output_file {
    my $self = shift;
    my $tag  = uc(shift);
    if ($tag) {
        if (my $path = shift) {
            my $opts = lc(shift || '');
            if ($path =~ /^\*/) {
                # User is passing a file handle
                $self->{FILES}{$tag} = {
                    FH   => $path,
                    PATH => '',
                    FRAG => '',
                };
            } elsif ($self->{INITIALIZED}) {
                $self->err("You can not set output_file() once forking has begun");
            } else {
                $self->{FILES}{ $tag } = {
                    PATH => $path,
                };
            }
            $self->{FILES}{ $tag }{OPTS} = $opts;
        }
        return $self->{FILES}{ $tag }{PATH} || '';
    }
    return '';
}

sub output_fh {
    my $self = shift;
    if ($_[0]) {
        return $self->{FILES}{ uc($_[0]) }{FH};
    }
    return undef;
}

sub write_output {
    my $self = shift;
    my $tag  = uc(shift);
    my $txt  = shift;
    return if (!defined $txt || $txt eq "");
    my $fh   = $self->{FILES}{$tag}{FH};
    if ($fh) {
        print $fh $txt;
    } else {
        $self->err("Attempt to write to closed file handle on tag '$tag'", 
                   "Perhaps you passed the wrong tag, or called reset()?",
                   $tag, $txt, $self->branch($self->{FILES}));
    }
}

sub execute {
    my $self = shift;
    my ($num) = @_;
    $num ||= 1;

    $self->total_fork($num);
    if ($self->verbose) {
        my @bits = (sprintf("Forking %d child%s - %s", $num, $num == 1 ?
                            '':'ren', &nice_date()));
        if (my $lim = $self->limit) {
            push @bits, "User request to process only $lim records";
        }
        $self->msg(@bits);
    }
    $self->{REFORK} = 0;
    my @pids;
    for my $i (1..$num) {
        my $pid = $self->fork( 'quiet' );
        push @pids, $pid;
    }
    $self->msg(sprintf("  Spawned %d Child%s: %s\n", $#pids + 1, $#pids == 0 ?
                       '' : 'ren', join(', ', @pids))) if ($self->verbose > 1);
    my $failed = $self->wait;
    $self->join_files();
    $self->msg(sprintf("All tasks completed - %s", &nice_date()))
        if ($self->verbose);
    $self->reset();
    return $failed;
}

sub fork {
    my $self = shift;
    my ($bequiet) = @_;
    my $num  = $self->{REFORK} || ++$self->{LASTKID};
    my $pid;
    if ($pid = CORE::fork) {
        # parent $pid = pid of child...
        push @{$self->{PIDS}}, $pid;
        $self->{PIDNUM}{$pid} = $num;
        $self->msg(sprintf("%spawning Child %d PID %d\n", 
                           $self->{REFORK} ? 'Res' : 'S', $num, $pid))
            if (! $bequiet && $self->verbose > 1);
    } elsif (defined $pid) {
        # $pid is zero but defined - this is the child
        # Each child calls one of the 'read' methods, eg read_info()
        $self->{CHILD} = $self->{REFORK} || $self->{LASTKID};
        my $lt   = $self->{TIME}{START} = time;
        $self->_init;
        #if ($@) {
        #    $self->graceful_death
        #        ("[!!]", "Child $num failed to initialize prior to fork",
        #         $@);
        #    exit 1;
        #}
        my $meth = $self->method;
        my $prog = $self->progress;
        while ( my $rec = $self->_next_task ) {
            &{$meth}($rec);
            if ($prog && time - $lt > $prog) {
                $self->_show_progress;
                $lt = time;
            }
        }
        $self->_finish;
        if ($prog && $self->verbose > 1) {
            my $elapsed = (time - $self->{TIME}{START}) / 60;
            my $rate    = $self->{COUNT} / ($elapsed || .01);
            $self->msg(sprintf("[%2d]", $self->{CHILD}), sprintf
                       ("Finished %d %s in %.1f min %.1f/min\n",
                         $self->{COUNT}, $self->{COUNTWHAT}, 
                         $elapsed, $rate));
        }
        if (my $ec = $self->exit_code()) {
            exit $ec;
        } else {
            # Fastest way to exit child, without time-consuming clean-up:
            kill('KILL', $$);
            CORE::dump();
            exit 0;
        }
    } else {
        $self->graceful_death("Failure to fork process for iteration $num");
    }
    return wantarray ? ($pid, $num) : $pid;
}

sub wait {
    my $self = shift;
    my @pidarray = @{$self->{PIDS}};
    $self->_create_count();
    my $failed = 0;
    my $expected = $self->exit_code();
    foreach my $pid (@pidarray) {
        waitpid($pid, 0);
        my $err = 0;
        my $exit_value = $? >> 8;
        # print "EXIT VALUE $exit_value\n";
        if (defined $expected) {
            unless ($exit_value == $expected) {
                $self->err("STACK_+2 Child $pid exits with exit value $exit_value, expected $expected");
                $err++;
            }
        } elsif ($exit_value) {
            $self->err("STACK_+2 Child $pid exits with exit value $exit_value");
            $err++;
        }

        if ($!) {
            # Bogus "Inappropriate ioctl for device" errors
            # print "    Child $pid throws error [!] $!\n";
            # $err++;
        }
        if ($@) {
           # print "    Child $pid throws error [@] $!\n";
        }
        if ($? & 128) {
            $self->err("Child $pid exits with core dump");
            # Remove core file - obviously this is not always desirable
            my $core_file = "core.$pid";
            unlink($core_file) if (-e $core_file);
            $err++;
        }
        if ($err) {
            $failed++;
            if ($self->allow_refork) {
                # The user wants to re-try the child when it fails
                $self->{REFORK} = $self->{PIDNUM}{$pid};
                $self->fork;
            }
        }
    }
    return $failed;
}

sub join_files {
    my $self = shift;
    my $fork = $self->fork_count;
    my %tobuild;
    my %opts;
    while (my ($tag, $data) = each %{$self->{FILES}}) {
        my $path = $data->{PATH};
        next unless ($path);
        $path = '>' . $path unless ($path =~ /^\>/);
        next unless ($path);
        my @files;
        for my $i (1..$fork) {
            # Iteratively check for all likely files
            my $file = sprintf("%s_%03d", $path, $i);
            $file =~ s/^\>+//;
            if (-e $file) {
                if (-s $file) {
                    push @files, $file;
                } else {
                    unlink($file);
                }
            }
        }
        next if ($#files < 0);
        $tobuild{$path} = \@files;
        $opts{$path}    = $data->{OPTS};
    }
    my @paths = keys %tobuild;
    return if ($#paths <0);
    $self->msg(sprintf("Assembling %d file%s - %s", $#paths + 1, $#paths == 0 ?
                       '' : 's', &nice_date())) if ($self->verbose > 1);

    while (my ($path, $files) = each %tobuild) {
        my ($ftok,$fname);
        if ($path =~ /^(\>+)(.+)$/) {
            ($ftok,$fname) = ($1, $2);
        } else {
            $self->err("The path '$path' is not set for write operations");
            next;
        }
        my (%errors, @to_kill);
        if ($#{$files} == 0) {
            # We just have a single file, we should be able to just rename it
            my $src = $files->[0];
            $src =~ s/ /\\ /g;
            my $trg = $fname;
            $trg =~ s/ /\\ /g;
            if ($ftok eq '>') {
                system("mv $src $trg");
            } else {
                system("cat $src >> $trg");
            }
            push @to_kill, $files->[0];
        } elsif ($opts{$path} =~ /sort/) {
            # We have sorted the files, now we need to merge them
            $self->bench_start('Merge sort');
            my $cmd = "sort "; 
            $cmd .= "-T $sortTmpDir " if ($sortTmpDir);
            $cmd .= "-m ".join(' ', @{$files})." $path";
            system($cmd);
            push @to_kill, @{$files};
            $self->bench_end('Merge sort');
        } else {
            # We have multiple files, we will need to concatenate
            $self->msg(sprintf("%s = %d fragment%s\n", $fname, $#{$files} + 1,
                               $#{$files} == 0 ? '' : 's'))
                if ($self->verbose > 1);
            if ($ftok eq '>') {
                sysopen(TOFILE, $fname, O_WRONLY | O_TRUNC | O_CREAT) ||
                    $self->graceful_death("Failed to generate output file",
                                          $fname, $!);
            } else {
                sysopen(TOFILE, $fname, O_WRONLY | O_APPEND | O_CREAT) ||
                    $self->graceful_death("Failed to concatenate output",
                                          $fname, $!);
            }
            
            foreach my $file (@{$files}) {
                sysopen(FROMFILE, $file, O_RDONLY) ||
                    $self->graceful_death
                    ("Failed to read from output fragment",
                     $file, $!);
                my $blksize = (stat FROMFILE)[11] || 16384;
                my $buffer;
                while (my $len = sysread(FROMFILE, $buffer, $blksize)) {
                    $self->graceful_death("sysread error", $file, $!)
                        if (!defined $len);
                    my $offset = 0;
                    while ($len) {
                        my $written = syswrite
                            (TOFILE, $buffer, $len, $offset);
                        if ($written) {
                            $offset += $written;
                            $len    -= $written;
                        } else {
                            $errors{$file} += $len;
                        }
                    }
                }
                push @to_kill, $file;
            }
            close TOFILE;
        }
        my @errs = sort keys %errors;
        unless ($#errs == -1) {
            $self->err("OUTPUT ERROR: $path",
                       "Failed to join information from component temp files",
                       (map {sprintf("%s : %d bytes", $_, $errors{$_})} @errs),
                       "Temporary files have NOT been removed:", @to_kill);
        } else {
            foreach my $kf (@to_kill) {
                unlink($kf);
            }
        }
    }
}

sub get_fh {
    my $self = shift;
    my ($file) = @_;
    if ($file !~ /^\>/) {
        # The user plans to read this file
        my $base = $file;
        $base =~ s/^[\<\>]+//;
        unless (-e $base) {
            $self->graceful_death("Could not read '$base'", "Does not exist");
        }
    }
    my ($fh, $ftype);
    $self->ignore_error('Inappropriate ioctl for device');
    undef $!;
    undef $@;
    if ($file =~ /\.gz$/) {
        $file =~ s/^[\<\>]+//;
        $ftype = 'gz';
        $self->ignore_error('Illegal seek');
        open($fh, "gunzip -c $file |");
        $self->ignore_error('Illegal seek', 'StopIgnoring');
    } else {
        unless ($file =~ /^[\<\>]/) {
            $file = "<$file";
        }
        $ftype  = '';
        open($fh, $file);
    }
    if (!$fh || ($! && $! ne 'Illegal seek')) {
        if ($fh) {
            $self->err
                ("Failed to recover file handle glob", $file,
                 $! ? '$! = '.$! : undef, $@ ? '$@ = '.$@ : undef);
        } else {
            $self->graceful_death
                ("Failed to open file handle", $file,
                 $! ? '$! = '.$! : undef, $@ ? '$@ = '.$@ : undef);
        }
    }
    return $fh;
}

sub _create_count {
    my $self = shift;
    my @cdat;
    my $cf        = $self->{COUNT_FILE};
    my $limit     = $self->limit;
    my $countFunc = $self->count_method();
    open(CF, ">$cf") || $self->graceful_death
        ("Failed to write count file", $cf, $!);
    foreach my $idat (@{$self->{INPUTS}}) {
        # It is possible that we have forked analysis of multiple files,
        # We need to count each seperately
        my ($input, $type) = @{$idat};
        my $num = $countFunc ? &{$countFunc}( $self, $input, $type ) : 0;
        print CF join("\t", $input, $num || 0)."\n";
    }
    close CF;
}


sub _show_progress {
    my $self = shift;
    my $child  = $self->child;
    my $hf     = $self->{PROG_HIST};

    # Update the history file to note how many records were parsed:
    open(HF, ">$hf") || $self->graceful_death
        ("Failed to write history file", $hf, $!);
    print HF $self->{COUNT} . "\n";
    close HF;

    my @tfs = $self->_task_files;
    foreach my $file (@tfs) {
        if ($file =~ /^Child_(\d+)$/) {
            # There is a lower numbered task still running - let *it*
            # report progress
            return if ($1 < $child);
        }
    }

    # Calculate the total number of tasks done:
    my $dir  = $self->_task_directory;
    my $count = 0;
    foreach my $file (@tfs) {
        if ($file =~ /^History/) {
            # Read the history
            open(HF, "<$dir/$file") || $self->graceful_death
                ("Failed to read history file", "$dir/$file", $!);
            while (<HF>) {
                chomp;
                $count += $_;
            }
            close HF;
        }
    }

    if (!defined $self->{TOTAL_COUNT} && -e $self->{COUNT_FILE}) {
        my $cf = $self->{COUNT_FILE};
        my $input = $self->input;
        my %counts;
        open(CF, "<$cf") || $self->graceful_death
            ("Failed to read count file", $cf, $!);
        while (<CF>) {
            chomp;
            my ($file, $count) = split(/\t/, $_);
            $counts{$file} = $count;
        }
        $self->{TOTAL_COUNT} = $counts{$input};
        close CF;
    }
    my $total   = $self->{TOTAL_COUNT};
    my $elapsed = (time - $self->{TIME}{START});
    my $rate    = $count / ($elapsed || .01);
    my $remain  = '';
    if ($total && $rate) {
        my ($r, $u) = (($total - $count) / $rate, 'sec');
        if ($r > 60) {
            $r /= 60;
            $u = 'min';
            if ($r > 60) {
                $r /= 60;
                $u = 'hr';
                if ($r > 24) {
                    $r /= 24;
                    $u = 'day';
                }
            }
        }
        $remain = sprintf(", %.1f %s remain", $r, $u);
    }
    my $li      = substr($self->last_item, 0, 50);
    $self->msg(sprintf("[%2d]",$child),sprintf
               ("%4d %s, %.1f min, %.1f per min%s - %s", $count, 
                $self->{COUNTWHAT}, $elapsed / 60, $rate * 60, $remain, $li));
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

*_task_dir = \&_task_directory;
sub _task_directory {
    my $self = shift;
    unless ($self->{TASKDIR}) {
        my $uniqTag = $$;
        if ($self =~ /x([\da-f]+)/) {
            # Use the hash key as a unique tag
            # $uniqTag = $1; # I forget why I was doing this...
            # $uniqTag =~ s/0+/0/g;
        }
        my $pdir = "/tmp/ForkCritter";
        unless (-d $pdir) {
            mkdir($pdir);
            chmod(0777, $pdir);
        }
        my $dir = $self->{TASKDIR} = "$pdir/PID_$uniqTag";
        if (-d $dir) {
            $self->_clear_task;
            # print "<<<<EXISTS: $dir\n". $self->stack_trace()."\n";
        } else {
            # print "<<<<MAKING: $dir\n". $self->stack_trace()."\n";
            mkdir( $dir, 0777 );
            chmod( 0777, $dir );
        }
        $self->{COUNT_FILE} = "$dir/RecordCount";
    }
    return $self->{TASKDIR};
}

sub _task_files {
    my $self = shift;
    my $dir  = $self->_task_directory;
    my @files;
    if (-d $dir) {
        opendir(TMPDIR, $dir) || $self->graceful_death
            ("Failed to read task file directory", $dir, $!);
        foreach my $file (readdir TMPDIR) {
            push @files, $file;
        }
        closedir TMPDIR;
    }
    return @files;
}

sub _clear_task {
    my $self = shift;
    my $dir  = $self->_task_directory;
    foreach my $file ($self->_task_files) {
        unlink("$dir/$file");
    }
}

sub _bogus_method {
    shift->graceful_death
        ("You never set the callback method()!",
         "Be sure to specify -method when calling new()",
         "Or explicitly call \$fc->method()");
}

sub _init {
    my $self = shift;
    return 0 if ($self->{INITIALIZED});
    $self->{FHNUM} = 0;
    my $type = $self->input_type;
    if ($type eq 'user') {
        # The user is initializing themselves
    } elsif ($type eq 'seq') {
        $self->_init_seq;
    } elsif ($type eq 'basic') {
        $self->_init_basic;
    } elsif ($type eq 'array') {
        $self->_init_array;
    } elsif ($type eq 'csv' || $type eq 'tsv' || $type eq 'maf') {
        $self->_init_tab;
    } elsif ($type eq 'xml') {
        # Wait - do it last
    } else {
        $self->graceful_death
            ( "Unknown input type", $type || '-NOT SPECIFIED-');
    }
    my $ctag  = sprintf("_%03d", $self->child);
    my $fork  = $self->fork_count;
    my $dir   = $self->_task_directory();
    my $cf    = $self->{TASK_FILE} = "$dir/Child$ctag";
    my $hf    = $self->{PROG_HIST} = "$dir/History$ctag";
    open(HF, ">$hf") || $self->graceful_death
        ("Failed to write history file", $hf, $!);
    print HF "0\n";
    close HF;
    open(TF, ">$cf") || $self->graceful_death
        ("Failed to write task file", $cf, $!);
    print TF "$$\n";
    close TF;
    while (my ($tag, $data) = each %{$self->{FILES}}) {
        my $path = $data->{PATH};
        next unless ($path);
        my $file = $path;
        $file .= $ctag unless ($file eq '/dev/null');
        $file  = '>' . $file  unless ($file =~ /^\>/);
        my $fh = $self->get_fh($file);
        $data->{FH}   = $fh;
        $data->{FRAG} = $file;
        $data->{FRAG} =~ s/^\>+//;
    }
    if (my $initmeth = $self->init_method) {
        &{$initmeth};
    }
    $self->{INITIALIZED} = 1;
    if ($type eq 'xml') {
        $self->_init_sax;
    }
}

# If this method is changed, _sax_wrapper() should also be changed
my $secPerDay = 60 * 60 * 24;
sub _next_task {
    my $self = shift;
    return undef if ($self->{TERMINATED} || $self->{FORCE_EXIT});
    if ($self->unstable_input) {
        # Verify that the input has not changed since execution started
        my $file = $self->input;
        $file =~ s/^\<+//;
        $self->graceful_death("Unstable Input - file no longer exists!",
                              $file) unless (-e $file);
        $self->graceful_death("Unstable Input - file is now zero size!",
                              $file) unless (-s $file);
        my $runTime = $self->{TIME}{START};
        my $modTime = $scriptStartTime - int((-M $file) * $secPerDay);
        my $diff    = $modTime - $runTime;
        $self->graceful_death("Unstable Input - file modified $diff seconds ".
                              "after analysis started!", $file) if ($diff > 30)
    }
    my $func  = $self->next_record_method;
    my $skip  = $self->skip_record_method;
    my $grpm  = $self->group_method();
    my $limit = $self->limit;
    my $retval;
    while (1) {
        # Return a null entry if we have exceeded the number of requested recs
        return undef if ($limit && $self->{RECORD} >= $limit);
        # Get the next record in the stream
        $retval = &{$func}( $self );
        if ($grpm) {
            # The user is grouping records
            my $record = $retval;
            $retval    = &{$grpm}( $record );
            # An undefined return value indicates an incomplete group
            next if ($record && !$retval);
        }
        # undef always indicates the end of the task
        return undef if (!defined $retval);
        next if ( $skip && &{$skip}( $retval ) );
        my $num = ++$self->{RECORD};
        # Keep cycling until we reach the appropriate modulus for this child,
        # unless this child is tasked to doall ($da) records in the input
        if ($self->{DO_BLOCK}) {
            $self->graceful_death("Charles needs to write match_block()!");
            last if ($self->match_block($num));
        } else {
            last if ($self->match_modulus($num));
        }
    }
    # printf("%s = %d %% %d\n", $retval->[0][0], $self->{RECORD},$self->child);
    $self->{COUNT}++;
    return $retval;
}

sub _finish {
    my $self = shift;
    return 0 if ($self->{TERMINATED}++);
    if (my $finishmeth = $self->finish_method()) {
        &{$finishmeth};
    }
    if (my $func = $self->{TERMFUNC}) {
        &{$func}( $self ) if ($func);
    }
    my $fork = $self->fork_count;
    $self->close_all();
    # Remove the task file
    unlink($self->{TASK_FILE});
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

sub _init_basic {
    my $self  = shift;
    my $input = $self->input;
    my $fh    = $self->get_fh($input);
    $self->{IO} = $fh;
    $self->{COUNTTAG}    = 'LINECOUNT';
    $self->{TERMFUNC}  ||= \&_finish_basic;
}

sub _finish_basic {
    my $self = shift;
    $self->{IO} = undef;
}

sub _count_basic {
    my $self = shift;
    my ($input, $type) = @_;
    # User defined methods, can not count
    return 0;
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

sub _init_array {
    my $self = shift;
    my $array = $self->input;
    if (!$array || !ref($array) || ref($array) ne 'ARRAY') {
        $array ||= '-UNDEF-';
        $self->graceful_death("Input() '$array' is not an array reference");
    }
    $self->{INDEX}       = 0;
    $self->{COUNTTAG}    = 'INDEX';
    $self->{COUNTWHAT}   = 'records';
    $self->{NEXTFUNC}  ||= \&_next_array;
    $self->{TERMFUNC}  ||= \&_finish_basic;
}

sub _next_array {
    my $self = shift;
    my $array = $self->input;
    my $index = $self->{INDEX}++;
    return undef if ($index > $#{$array});
    my $retval = $array->[ $index ];
    $self->last_item( $retval );
    return $retval;
}

sub _count_array {
    my $self  = shift;
    my $array = shift;
    return $#{$array} + 1;
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

sub _init_seq {
    my $self  = shift;
    my $input = $self->input;
    # Use a filehandle to allow analysis of .gz fasta files
    my $fh    = $self->get_fh("<$input");
    my $st    = $self->{SUBTYPES} || {};
    my $fmt   = $st->{SEQFORMAT};
    $self->{IO}          = $fh;
    $self->{SEQIO}       = Bio::SeqIO->new( -fh => $fh, -format => $fmt );
    if ($st->{PAIRED}) {
        if ($input =~ /(.+)1.([^\.]{2,7})$/) {
            my $inp2 = $1.'2.'.$2;
            if (-s $inp2) {
                my $fh    = $self->get_fh("<$inp2");
                $self->{SEQIO2} = Bio::SeqIO->new
                    ( -fh => $fh, -format => $fmt );
            } else {
                $self->death("Paired SeqIO request, but failed to find pair",
                             $inp2);
            }
        } else {
            $self->death("Paired SeqIO request, but can not interpret primary",
                         $input);
        }
    }
    $self->{NEXTFUNC}  ||= \&_next_seq;
    $self->{TERMFUNC}  ||= \&_finish_seq;
    $self->{COUNTTAG}    = 'RECORD';
    $self->{COUNTWHAT}   = 'records';
    unless ($self->last_item_method) {
        $self->last_item_method( sub {
            my $seq = shift;
            return $seq->display_id;
        });
    }
}

sub _next_seq {
    my $self   = shift;
    my $retval = $self->{SEQIO}->next_seq;
    return undef unless ($retval);
    $self->last_item( $retval );
    if (my $sio2 = $self->{SEQIO2}) {
        if (my $bs = $sio2->next_seq()) {
            $retval = [ $retval, $bs ];
        }
    }
    return $retval;
}

sub _finish_seq {
    my $self = shift;
    foreach my $key ('SEQIO', 'SEQIO2') {
        if (my $reader = $self->{$key}) {
            # We are going to handle filehandles on our own
            # BioPerl spews out a bunch of annoying errors via Bio::Root::IO
            $reader->{_filehandle} = undef;
            #$self->{IO}    = undef;
        }
    }
}

sub _count_fasta {
    my $self = shift;
    my ($input) = @_;
    my $num   = 0;
    my $fh    = $self->get_fh($input);
    my $limit = $self->limit;
    while (<$fh>) {
        $num++ if (/^\>/);
        last if ($limit && $num >= $limit);
    }
    close $fh;
    return $num;
}

sub _count_fastq {
    my $self = shift;
    my ($input) = @_;
    my $num   = 0;
    my $fh    = $self->get_fh($input);
    my $limit = $self->limit;
    while (<$fh>) {
        $num++ if (/^\@/);
        last if ($limit && $num >= $limit);
    }
    close $fh;
    return $num;
}

sub _count_genbank {
    my $self = shift;
    my ($input) = @_;
    my $num   = 0;
    my $fh    = $self->get_fh($input);
    my $limit = $self->limit;
    while (<$fh>) {
        $num++ if (/^LOCUS /);
        last if ($limit && $num >= $limit);
    }
    close $fh;
    return $num;
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

my $globalSelf;
sub _init_sax {
    my $self  = shift;
    my $input   = $self->input;
    my $args    = $self->input_args() || [];
    $globalSelf = $self;
    $self->{NEXTFUNC}  ||= \&_next_sax;
    $self->{TERMFUNC}  ||= \&_finish_sax;
    $self->{LASTTIME}    = time;
    $self->{COUNTTAG}    = 'RECORD';
    $self->{COUNTWHAT}   = 'records';
    unless ($self->last_item_method) {
        if ($self->group_method) {
            $self->last_item_method( sub {
                return "Node Group " . $globalSelf->{RECORD};
            });
        } else {
            $self->last_item_method( sub {
                my $rec = shift;
                return sprintf("<%s> %d", $rec->{NAME}, $globalSelf->{RECORD});
            });
        }
    }
    eval {
        # Now get all the nodes
        my $fs = BMS::FriendlySAX->new
            ( -file    => $input,
              # -limit   => $self->limit,
              -verbose => 0,
              -method  => \&_sax_wrapper,
              -limitmsg => $self->{LIMITMSG},
              @{$args},
              );
    };
    if (my $grpm = $self->group_method) {
        # Grouping was being used, and we need to deal with the final group
        my $residual = &{$grpm}();
        # Do not want to trigger grouping logic when we pass the group object
        # (rather than a FriendlySAX node)
        $self->{NOGROUP} = 1;
        &_sax_wrapper($residual) if ($residual);
        $self->{NOGROUP} = 0;
    }
    if ($@) {
        my $expected = $self->{LIMITMSG} || 'user limit';
        unless ($@ =~ /\Q$expected\E|\[IGNORE\]/i) {
            $self->{EXITCODE} = 1;
            $self->graceful_death("FriendlySAX error", $@);
        }
    }
}

sub _sax_wrapper {
    my ($node) = @_;
    my $self   = $globalSelf;
    my $limit  = $self->limit;
    my $grpm  = $self->group_method;
    if ($limit && $self->{RECORD} >= $limit) {
        my $num  = $self->{RECORD};
        my $what = $grpm ? "Node Group $num" : sprintf
            ("<%s> %d", $node->{NAME} || "??", $num );
        my $msg = $self->{LIMITMSG} ||
            "User limit ($limit) halts processing on $what";
        if ($self->{NOGROUP}) {
            $self->msg($msg) if ($self->{CHILD} == $self->{LOWKID}
                                 && !$self->{QUIETLIMIT});
            return 0;
        } else {
            if ($self->{QUIETLIMIT}) {
                # $self->graceful_shutdown();
                die;
            } else {
                # $self->graceful_shutdown();
                die $msg;
            }
        }
    }
    if ($grpm) {
        # The user is grouping records
        $node = &{$grpm}( $node ) unless ($self->{NOGROUP});
        # An undefined return value indicates an incomplete group
        return unless ($node);
    }
    if (my $skip = $self->skip_record_method) {
        return if (&{$skip}( $node ));
    }
    my $num    = ++$self->{RECORD};
    return unless ($self->match_modulus($num));
    $self->last_item( $node );
    $self->{COUNT}++;
    my $meth   = $self->method;
    my $prog   = $self->progress;
    &{$meth}($node);
    if ($prog && time - $self->{LASTTIME} > $prog) {
        $self->_show_progress;
        $self->{LASTTIME} = time;
    }
}

sub _next_sax {
    # Cycling will be accomplished by FriendlySAX, just return null
    return undef;
}

sub _finish_sax {
    my $self = shift;
    $self->{IO}    = undef;
}

sub _count_xml {
    my $self = shift;
    my ($input) = @_;
    my $num     = 0;
    my $args    = $self->input_args() || [];
    my $pargs   = $self->parseparams( @{$args} );
    my $tags    = $pargs->{TAG} || $pargs->{TAGS};
    my $limit   = $self->limit;
    if ($tags) {
        $tags = [ $tags ] unless (ref($tags));
        my %thash = map { $_ => 1 } @{$tags};
        my $fh = $self->get_fh($input);
        while (<$fh>) {
            if ( /\<\s*([^\/]\S+)/) {
                my $tag = $1; $tag =~ s/\>//;
                if ($thash{$tag}) {
                    $num++;
                    last if ($limit && $num >= $limit);
                }
            }
        }
        close $fh;
    } else {
        $num = 1;
    }
    return $num;
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

sub _init_tab {
    my $self  = shift;
    my $input = $self->input;
    my $type  = $self->input_type;
    my $fh    = $self->get_fh($input);
    $self->{IO} = $fh;
    if ( $self->{IOSEP}) {
        # Separator has been manually defined
    } elsif ($type eq 'csv') {
        $self->{IOSEP}   = '[\\\'\\"],[\\\'\\"]';
        $self->{IOSTRIP} = '[\\\'\\"]';
    } else {
        $self->{IOSEP}   = "\t";
        $self->{IOSTRIP} = "";
    }
    my $func = $self->{NEXTFUNC} ||= \&_next_tab;
    $self->{TERMFUNC} ||= \&_finish_basic;
    if ($self->{SUBTYPES}{HEADER}) {
        # There is a header attached to the file
        my $header;
        my $lookFor = $self->{SUBTYPES}{LOOKFOR};
        while (1) {
            $header = &{$func}( $self );
            if (!$header) {
                last;
            } elsif ($lookFor) {
                my $chk = join('\t', @{$header});
                last if ($chk =~ /$lookFor/i);
            } else {
                last;
            }
        }
        $self->death("Failed to recover header") unless ($header);
        my @mapped = map { $self->remap_header_name( $_ ) } @{$header};
        $self->{HEAD_ARRAY} = \@mapped;
        $self->{SUBTYPES}{HEADER} = \@mapped;
        $self->{HEADER}     = \@mapped if ($self->{SUBTYPES}{HASH});
        #for my $i (0..$#{$header}) { printf("[%3d] %s\n", $i, $header->[$i]);}
    }
    $self->{COUNTTAG}    = 'LINECOUNT';
    if (my $gb = $self->{SUBTYPES}{GROUPBY}) {
        # The user wants results grouped by a column
        $self->{GROUP}    = $gb;
        $self->{STACK}    = [];
        $self->{PUSHBACK} = [];
        $self->{NEXTFUNC} = \&_next_tab_group;
        # $self->{COUNTTAG} = 'RECORD';
    }
    unless ($self->last_item_method) {
        $self->last_item_method( sub {
            my $line = shift;
            $line = "" unless (defined $line);
            $line =~ s/[\n\r]+$//;
            return $line;
        });
    }
}

sub _next_tab {
    my $self = shift;
    my $fh   = $self->{IO};
    my $line = <$fh>;
    return undef unless ($line);
    $self->{LINECOUNT}++;
    $self->last_item( $line );
    chomp $line;
    if (my $strip = $self->{IOSTRIP}) {
        # Clean edges:
        $line =~ s/^$strip//;
        $line =~ s/$strip$//;
    }
    # Remove non-ascii
    $line =~ s/\P{IsASCII}//g if ($self->{SUBTYPES}{CLEANASCII});
    my $sep  = $self->{IOSEP};
    my @list = split($sep, $line);
    my $retval = \@list;
    if (my $head = $self->{HEADER}) {
        # The user wants the data back as a hash
        my %hash;
        for my $i (0..$#{$head}) {
            $hash{ $head->[$i] } = $list[$i];
        }
        $retval = \%hash;
    }
    return $retval;
}

sub _next_tab_group {
    my $self   = shift;
    my $in     = $self->{GROUP};
    my $ishash = $self->{HEADER};
    my @group;
    while (1) {
        my $rec = shift @{$self->{PUSHBACK}} || $self->_next_tab;
        last unless ($rec);
        if ($#group < 0) {
            # First record - we need to seed the group
            push @group, $rec;
            next;
        }
        if ($ishash) {
            # We are grouping by hash key
            if ($group[-1]{$in} eq $rec->{$in}) {
                push @group, $rec;
                next;
            }
        } elsif ($group[-1][$in] eq $rec->[$in]) {
            # We are grouping by array index
            push @group, $rec;
            next;
        }
        # This record should be assigned to the next group
        push @{$self->{PUSHBACK}}, $rec;
        last;
    }
    return undef if ($#group < 0);
    $self->last_item('');
    return \@group;
}

sub _count_tab {
    my $self = shift;
    my ($input) = @_;
    my $num   = 0;
    my $fh    = $self->get_fh($input);
    my $limit = $self->limit;
    while (<$fh>) {
        $num++;
        last if ($limit && $num >= $limit);
    }
    close $fh;
    return $num;
}

sub nice_date {
    my $dt = `date +'%d %b %H:%M:%S'`;
    $dt =~ s/\s*[\n\r]+//;
    return $dt;
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
1;
