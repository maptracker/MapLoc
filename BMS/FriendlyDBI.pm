# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::FriendlyDBI;
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

# EXEC DBMS_UTILITY.analyze_schema('SCHEMA_NAME','COMPUTE');
# ANALYZE TABLE foo COMPUTE STATISTICS;
# ANALYZE INDEX foo_ind COMPUTE STATISTICS;

#BEGIN {
#    # use lib '/stf/sys/perl/newlib/';
#    use lib '/perl/newlib';
#}

use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK);
use Scalar::Util qw(weaken);
require Exporter;

use DBI;

@ISA    = qw(DBI);

our $fdbLastError = "";
our $colTypeAliases = {
    postgres => {
        number => 'numeric',
        date   => 'timestamp',
        
    },
    oracle => {
        numeric => 'number',
    },
};


$BMS::FriendlyDBI::VERSION = 
    ' $Id: FriendlyDBI.pm,v 1.66 2014/03/20 18:18:56 tilfordc Exp $ ';
@BMS::FriendlyDBI::EXPORT = qw($fdbLastError);

sub connect {
    my $self = shift;
    my @fdArgs = ( -db => $_[1], -passwd => $_[2], -opts => $_[3] );
    push @fdArgs, splice(@_, 4);
    my $args = {};
    for (my $i = 0; $i < $#fdArgs; $i += 2) {
        my ($key, $val) = (uc($fdArgs[$i]), $fdArgs[$i+1]);
        $key =~ s/^\-//;
        $args->{$key} = $val;
    }
    
    my $name  = $args->{DB};
    if (defined $name) {
        ($name) = split('/', $name);
    } else {
        $name = '-UNDEFINED-';
    }
    my $retry  = $args->{RETRY} || 1;
    my $reWait = $args->{WAIT};
    my $dbh;
    my @errs;
    for my $r (1..$retry) {
        eval { $dbh = $self->SUPER::connect( @_ ); };
        last if ($dbh);
        my $err = $@ || 'Unknown error';
        if (my $de = $DBI::errstr) { $err .= "; $de" };
        push @errs, $err;
        warn "Failed to connect ($err), retrying...\n" 
            if ($args->{RETRYVB} && $r != $retry);
        sleep($reWait) if ($reWait);
    }
    if (!$dbh || $#errs + 1 >= $retry) {
        my $html  = $ENV{'HTTP_HOST'} ? 1 : 0;
        my $msg   = "Unable to connect to database ";
        $msg .= sprintf("%s%s%s (%s)", $html ? "<font size='+2'><b>" : '',
                        $name, $html ? '</b></font>' : '', $args->{DB});
        my $evErr = join(' // ', @errs);
        $msg .= "\nEval Error:\n   $evErr" if ($evErr);
        $fdbLastError = $msg;
        $fdbLastError =~ s/<[^>]+>//g;

        $msg .= "\nArguments:\n";
        my @pargs = (['Connect String', $_[0]]);
        push @pargs, ['DBH Acquired', $dbh] if ($dbh);
        push @pargs, map { [ $_, $args->{$_} ] } sort keys %{$args};
        while (my $data = shift @pargs) {
            my ($key, $val) = @{$data};
            if ($key eq 'OPTS') {
                push @pargs, 
                map { ["Option: $_",$val->{$_}] } sort keys %{$val};
                next;
            }
            $val = '-UNDEF-' unless (defined $val);
            $msg .= "  $key: $val\n";
        }
        &BMS::FriendlyDBI::db::verbose_error($msg, $args);
        die;
    }
    # Set some defaults for the database handle
    $dbh->{private_DEBUG}   = 0;
    $dbh->{private_DB_NAME} = $name;
    $dbh->{private_DB_FULL} = $args->{DB};
    my @bits = split(/\//, uc( $dbh->{Username}) );
    # Suddenly getting very irritating errors from DBI such as:
    # DBI db handle 0x21c0dcd0 has uncleared implementors data
    # These appear to be thrown during object destruction, and it
    # appears they can be safely ignored
    $dbh->ignore_error('uncleared implementors data');
    $dbh->schema_owner( $bits[0] );
    $dbh->error_file( $args->{ERRORFILE} );
    $dbh->admin_mail( $args->{ADMINMAIL} );
    $dbh->no_environment( $args->{NOENV} || $args->{NOENVIRONMENT} );
    return $dbh;
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::FriendlyDBI::db;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

use strict;
use BMS::Utilities::Benchmark;
use BMS::ErrorInterceptor;

use vars qw(@ISA);
@ISA = qw(BMS::Utilities::Benchmark BMS::ErrorInterceptor DBI::db );

# Both BMS::Utilities and DBI:db define err()
# Typeglob a new name for DBI::db to avoid unpleasant collisions:

*dbi_err = \&DBD::_::common::err;

sub DESTROY {
    my $self = shift;
    return unless ($self);
#    eval {
#        if (my $leftover = $self->{private_NESTED_BEGIN}) {
#            $self->error("WARNING: It appears that commit() was not called - ".
#                         "needed to call at least $leftover more times");
#            $self->{private_NESTED_BEGIN} = 1;
#            $self->commit();
#        }
#    };
    # Adding to stop dbih_clearcom warnings:
    my $rv = $self->disconnect();
    # warn "Failed to disconnect:\n". &stack_trace() unless ($rv);
    $self->SUPER::DESTROY( @_ );
}

sub fork_safe   { 
    my $self = shift;
    my $aggressive = shift;
    $self->{InactiveDestroy} = 1;
    if (my $nsths = $self->{private_NAMED_STHS}) {
        if ($aggressive) {
            $self->{private_NAMED_STHS} = {};
        } else {
            foreach my $sth (values %{$nsths}) {
                $sth->{InactiveDestroy} = 1;
            }
        }
    }
}
sub fork_unsafe {
    my $self = shift;
    $self->{InactiveDestroy} = 0;
    foreach my $sth (values %{$self->{private_NAMED_STHS}}) {
        $sth->{InactiveDestroy} = 0;
    }
}

*sql_dump = \&debug_level;
*dump_sql = \&debug_level;
*dumpsql  = \&debug_level;
*sqldump  = \&debug_level;
sub debug_level {
    my $self = shift;
    $self->{private_DEBUG} = $_[0] if (defined $_[0]);
    return $self->{private_DEBUG};
}

sub admin_mail {
    my $self = shift;
    $self->{private_ADMIN_MAIL} = $_[0] if (defined $_[0]);
    return $self->{private_ADMIN_MAIL};
}

sub error_file {
    my $self = shift;
    if (my $newval = $_[0]) {
        $self->{private_ERROR_FILE} = $newval;
        unless (-e $newval) {
            system("touch $newval");
            chmod(0777, $newval);
        }
    }
    return $self->{private_ERROR_FILE};
}

sub driver {
    my $self = shift;
    unless ($self->{private_DRIVER_TYPE}) {
        # $self->dbh_msg("Determining database driver");
        my $driver = lc($self->get_info( 17 ) || '');
        if ($driver =~ /mysql/) {
            $driver = 'mysql';
        } elsif ($driver =~ /oracle/) {
            $driver = 'oracle';
        } elsif ($driver =~ /(postgres|pg)/) {
            $driver = 'postgres';
        } else {
            $self->verbose_error("Unknown database driver '$driver'");
        }
        $self->{private_DRIVER_TYPE} = $driver;
    }
    return $self->{private_DRIVER_TYPE};
}

sub limit_syntax {
    my $self = shift;
    unless ($self->{private_LIMIT_SYNTAX}) {
        my $driver = $self->driver();
        if ($driver =~ /^(mysql|postgres)$/) {
            $self->{private_LIMIT_SYNTAX} = "LIMIT";
        } else {
            $self->{private_LIMIT_SYNTAX} =  "AND rownum <=";
        }
    }
    return $self->{private_LIMIT_SYNTAX};
}

sub _eval_can {
    my $self = shift;
    my $name   = shift;
    $name      =~ s/[^a-z_0-9]//gi;
    my $driver = $self->driver();
    # For some reason can() will not work properly
    # Instead, use eval() to try to find the method for this driver
    my $sname  = "_${name}_$driver";
    my $cb;
    eval('$cb = \\&'.$sname);
    return $cb;
}

sub insert_array {
    my $self = shift;
    unless ($self->{private_BULK_INSERT}) {
        # I can't get can() to work here...
        my $cb = $self->{private_BULK_INSERT} = 
            $self->_eval_can('insert_array');
        unless ($cb) {
            $self->dbh_die("No bulk insert method defined for ".
                           $self->driver());
        }
    }
    return &{$self->{private_BULK_INSERT}}($self,  @_ );
}

sub _insert_array_postgres {
    my $self  = shift;
    $self->bench_start();
    my ($table, $array) = @_;
    # $self->execute("COPY $table FROM STDIN"); # Old syntax
    $self->do("COPY $table FROM STDIN");
    for my $i (0..$#{$array}) {
        my $row = $array->[$i];
        next unless ($row);
        my @esc;
        foreach my $v (@{$row}) {
            if (!defined $v || $v eq '') {
                push @esc, '\N';
            } else {
                # Must escape the escape characters!
                $v =~ s/\\/\\\\/g;
                push @esc, $v;
            }
        }
        my $line = join("\t", @esc);
        # $self->func("$line\n", 'putline'); # Old syntax
        $self->pg_putcopydata($line."\n");
    }
    #$self->func("\\.\n", 'putline'); # Old syntax
    # If there are constraint violations, they will occur here:
    # my $rv = $self->func('endcopy'); # Old syntax
    my $rv =$self->pg_putcopyend();
    $self->bench_end();
    return $rv;
}

sub _insert_array_oracle {
    my $self  = shift;
    $self->bench_start();
    my ($table, $array) = @_;

    # NOT TESTED YET

    # Just loop over array and insert
    # Need a more efficient method for this
    my $ins = $self->{private_ORACLE_INSERT_STH}{$table};
    unless ($ins) {
        my @cols = $self->column_order( $table );
        
        $ins = $self->{private_ORACLE_INSERT_STH}{$table} =
            $self->prepare
            ( -sql    => "INSERT INTO TABLE $table (".join
              (', ', @cols).") VALUES (".join(', ', map {'?'} @cols).")",
              -name   => "Insert a row into table $table",
              -level  => 3, );
    }
    foreach my $row (@{$array}) {
        $ins->execute( @{$row} );
    }
    $self->bench_end();
}

sub nextval {
    my $self  = shift;
    my ($seq) = @_;
    my $sth   = $self->{private_NEXTVAL}{uc($seq)};
    unless ($sth) {
        my $frm   = $self->{private_NEXTVAL_SYNTAX};
        unless ($frm) {
            my $driver  = $self->driver();
            if ($driver =~ /^(postgres)$/) {
                $frm    = $self->{private_NEXTVAL_SYNTAX} =
                    "SELECT nextval('\%s')";
            } else {
                $frm = $self->{private_NEXTVAL_SYNTAX} =
                    "SELECT \%s.nextval FROM DUAL";
            }
        }
        $seq = uc($seq);
        $sth = $self->{private_NEXTVAL}{$seq} = $self->prepare
            ( -sql    => sprintf($frm, $seq),
              -name   => "Next sequence value for $seq",
              -level  => 3, );
    }
    return $sth->get_single_value();
}

*currval = \&lastval;
sub lastval {
    my $self  = shift;
    my ($seq) = @_;
    my $sth   = $self->{private_LASTVAL}{uc($seq)};
    unless ($sth) {
        my $frm   = $self->{private_LASTVAL_SYNTAX};
        unless ($frm) {
            my $driver  = $self->driver();
            if ($driver =~ /^(postgres)$/) {
                $frm    = $self->{private_LASTVAL_SYNTAX} =
                    "SELECT last_value FROM \%s";
            } else {
                $frm = $self->{private_LASTVAL_SYNTAX} =
                    "SELECT \%s.currval FROM DUAL";
            }
        }
        $seq = uc($seq);
        $sth = $self->{private_LASTVAL}{$seq} = $self->prepare
            ( -sql    => sprintf($frm, $seq),
              -name   => "Last sequence value for $seq",
              -level  => 3, );
    }
    return $sth->get_single_value();
}

sub no_environment {
    my $self = shift;
    my $nv   = shift;
    if (defined $nv) {
        $self->{private_NOENV} = $nv ? 1 : 0;
    }
    return $self->{private_NOENV};
}

sub verbose_error {
    my $self;
    $self = shift if (ref($_[0]));
    my ($txt, $args) = @_;
    my $name  = $args->{DB};
    my $owner = $self->admin_mail() if ($self);
    my $efile = $self->error_file() if ($self);
    $efile  ||= $args->{ERRORFILE};
    $owner  ||= $args->{ADMINMAIL} || 'tilfordc@bms.com';
    $name   ||= $self ? $self->{private_DB_NAME} : 'uninitialized';
    ($name)   = split('/', lc($name));
    my $full  = $self ? $self->{private_DB_FULL} : 'uninitialized';
    my $html  = $ENV{'HTTP_HOST'} ? 1 : 0;
    my $nl    = $html ? "<br />" : "\n";
    my $dt    = `date`; chomp($dt);
    my $ignoreAll = {
        '00604' => "error occurred at recursive SQL level 1 (Ctrl-C)",
        '01013' => "user requested cancel of current operation (Ctrl-C)",
        '03106' => "fatal two-task communication protocol error (Ctrl-C)",
    };
    my $admin = {
        '00257' => "Disk space exhausted during archiving",
        '00600' => "Serious fatal error on Oracle backend",
        '01014' => "The database is being brought down, possibly by the DBAs",
        '01033' => "The database is going down (or up), possibly by the DBAs",
        '01034' => "Database is totally inaccessible",
        '01035' => "Generally indicates database maintenance",
        '01077' => "Required Oracle program failed to start",
        '01078' => "An Oracle startup file is improperly formatted",
        '01089' => "The database is being brought down, possibly by the DBAs",
        '01090' => "The database is being brought down, possibly by the DBAs",
        '01092' => "Bad Things have happened to the database",
        '01109' => "The database is not available yet",
        '01110' => "An Oracle support file has a problem",
        '01113' => "Something is wrong with an Oracle file",
        '01119' => "Something is wrong with an Oracle file",
        '01219' => "The database is not available yet",
        '01254' => "The database is being backed up and not available",
        '01340' => "The National Language Support function is broken",
        '01547' => "Database recovery is incomplete",
        '01562' => "Generally means insufficient disk space is available",
        '01650' => "Generally means insufficient disk space is available",
        '01652' => "Generally means insufficient disk space is available",
        '01653' => "Generally means insufficient disk space is available",
        '01654' => "Generally means insufficient disk space is available",
        '01658' => "Generally means insufficient disk space is available",
        '01683' => "Generally means insufficient disk space is available",
        '01685' => "DB is not configured to allow for more data (?)",
        '02068' => "Something is Very Wrong",
        '03135' => 'The database abruptly terminated, maybe timeout',
        '04030' => "The computer running Oracle has run out of RAM",
        '04031' => "The computer running Oracle has run out of RAM",
        '07445' => "Something Very Bad happened to Oracle",
        '08103' => "Mysterious failure on Oracle backend",
        '12154' => "TNS (resource defining database) is unhappy",
        '12162' => "TNS (resource defining database) is unhappy",
        '12170' => "TNS (resource defining database) timeout is too short",
        '12203' => "TNS (resource defining database) is unhappy",
        '12222' => "TNS (resource defining database) is unhappy",
        '12224' => "TNS (resource defining database) is unhappy",
        '12500' => "TNS (resource defining database) could not be started",
        '12504' => "TNS (resource defining database) is unhappy",
        '12505' => "TNS (resource defining database) is unhappy",
        '12514' => "TNS (resource defining database) is unhappy",
        '12516' => "TNS (resource defining database) is unhappy",
        '12519' => "TNS (resource defining database) is unhappy",
        '12520' => "TNS (resource defining database) is unhappy",
        '12523' => "TNS (resource defining database) is unhappy",
        '12528' => "TNS (resource defining database) is unhappy",
        '12533' => "TNS (resource defining database) is unhappy",
        '12534' => "TNS (resource defining database) is unhappy",
        '12535' => "TNS (resource defining database) is unhappy",
        '12537' => "TNS (resource defining database) is unhappy",
        '12538' => "TNS (resource defining database) is unhappy",
        '12540' => "TNS (resource defining database) is unhappy",
        '12541' => "TNS (resource defining database) is not listening",
        '12542' => "TNS (resource defining database) is unhappy",
        '12543' => "TNS (resource defining database) is unhappy",
        '12547' => "TNS (resource defining database) is unhappy",
        '12552' => "TNS (resource defining database) is unhappy",
        '12554' => "TNS (resource defining database) is unhappy",
        '12557' => "TNS (resource defining database) is unhappy",
        '12560' => "TNS (resource defining database) is unhappy",
        '12564' => "TNS (resource defining database) is unhappy",
        '12570' => "TNS (resource defining database) is unhappy",
        '12571' => "TNS (resource defining database) is unhappy",
        '27101' => "Database is totally inaccessible",
        # '01035' => "Database access permissions have been restricted",
    };

    if ($efile) {
        if (open(EFILE, ">>$efile")) {
            print EFILE "[x] $dt\n$txt\n";
            close EFILE;
            $txt .= "\nError log at $efile";
            chmod(0666, $efile);
        } else {
            $txt .= "\nFailed to append error text to $efile\n  $!\n";
        }
    }

    my @admins = ('mg-pri-dtcsdba@bms.com');
    my $msg    = ($html) ? "<h3>$dt</h3>" : $dt;
    $msg      .= "\n$txt\n";
    my $erst = $DBI::errstr || '';
    $msg .= sprintf($html ? '<pre>%s</pre>' : '%s', &stack_trace());
    if ($erst) {
        $msg .= $html ? 
            "\n<b>DBI Error:</b>\n$erst\n" : "\nDBI Error:\n$erst\n";
    }
    my $to     = $owner;
    my $chkErr = 1; # Check first line for error code
    my $code;
    foreach my $line (split(/[\n\r]+/, $erst)) {
        if ($chkErr && $line =~ /ORA\-(\d+)/) {
            $code = $1;
            last;
        } elsif ($line =~ /Eval Error/) {
            # Error code can be on next line
            $chkErr = 1;
        } elsif ($line =~ /^\s*ERROR\s*\:\s+(.+)/) {
            $code = $1;
            $code =~ s/\"[^\"]+\"/ /g;
            $code =~ s/at character \d+$//;
            $code =~ s/for table/ /;
            $code =~ s/in column/ /;
            $code =~ s/\s+/ /g;
            $code =~ s/^ //; $code =~ s/ $//;
            $code = "PG-ERR $code";
            last;
        } else {
            $chkErr = 0;
        }
    }
    my $adminMsg = "";
    if ($code && $code =~ /^\d+$/) {
        if (my $igAll  = $ignoreAll->{$code}) {
            # This is an "expected" error, do not get verbose about it
            die "Expected error: ORA-$code : $igAll NO_STACK\n";
        }
        if (my $reason = $admin->{$code}) {
            $to   = join(",", @admins);
            $adminMsg .=  "This is a known issue : $reason$nl".
                "Database administrators have been alerted: $to\n";
            $to  .= ",$owner" if ($owner);
        }
        $code = "ORA-$code";
    }
    $code ||= 'UNKNOWN ERROR';
    my $subj = "$code for $name";
    my $cmd  = qq(| Mail -s '$subj' $to);
    my $sdir = "/tmp/$name";
    unless (-e $sdir) {
        mkdir($sdir);
        chmod(0777, $sdir);
    }
    $code =~ s/[^a-z0-9\-]+/_/gi;
    my $safety = "$sdir/$code.history";
    if (-e $safety) {
        my $days = (-M $safety);
        my $mins = $days * 24 * 60;
        # Only send a warning once an hour
        if ($mins < 60) {
            $cmd = "";
            $msg .= sprintf("\nAlert mail last sent %.1f minutes ago, no message sent now\n   Alert File: $safety\n", $mins);
        }
    }
    if (open(SAFETY, ">>$safety")) {
        print SAFETY "$0\t" . `date`;
        close SAFETY;
        chmod(0666, $safety);
        $txt .= "\nError History: $safety";
    } else {
        $txt .= "\nFailed to write safety file $safety\n  $!\n";
        $cmd  = "";
    }
    my $ev = $html ?
        "<table border='1'><caption>Environmental Variables</caption><tbody>".
        "\n  <tr><th>Variable</th><th>Value</th></tr>\n" :
        "Environmental Variables\n";
    my %skipenv = map { $_ => 1 } qw(COLORTERM LS_OPTIONS MINICOM HTTP_COOKIE GTK_PATH LS_COLORS MANPATH XFILESEARCHPATH CVSROOT CTIME_ROOT APOLLO_ROOT INFODIR INFOPATH INTEL_LICENSE_FILE);
    my %edat;
    foreach my $var (sort keys %ENV) {
        next if ($skipenv{$var});
        my $val = $ENV{$var} || '';
        $val =~ s/\n/ /g;
        if ($var =~ /(PATH|LIB)$/) {
            $edat{$var} = [ split(/\:/, $val) ];
        } else {
            $edat{$var} = $val;
        }
    }
    foreach my $var (sort keys %edat) {
        my $val = $edat{$var} || '';
        my @vals = ref($val) ? @{$val} : ($val);
        map { s/\&/\&amp;/g;
              s/\>/\&gt;/g;
              s/\</\&lt;/g; } @vals;
        if ($html) {
            $ev .= sprintf("  <tr><th>%s</th><td>%s</td></tr>\n",
                           $var, join("<br />", @vals));
        } else {
            for my $i (0..$#vals) {
                $ev .= sprintf("  %20s : %s\n", $i ? '' : $var, $vals[$i]);
            }
        }
    }
    $ev .= "</tbody></table>\n" if ($html);
    if ($cmd) {
        my $shortProg = $0;
        if ($shortProg =~ /\/([^\/]+)$/) { $shortProg = $1 };
        my $user = 
            $ENV{'HTTP_MAIL'} || $ENV{'REMOTE_USER'} ||
            $ENV{'LDAP_USER'} || $ENV{'HTTP_CN'} ||
            $ENV{'USER'}      || $ENV{'LOGNAME'} || $ENV{'REMOTE_ADDR'};
        my $link = $ENV{REQUEST_URI};
        my $serv = $ENV{HTTP_HOST};
        open (MAIL, $cmd) || die "Could not send mail to Admins!\n  $!";
        print MAIL "$dt\n** AUTOMATED MESSAGE - generated by $shortProg **\n\n$txt";
        print MAIL "\n\nError: $erst" if ($erst);
        print MAIL "\n\nProcess: $$ $0";
        print MAIL "\n\nDB: $full"    if ($full);
        print MAIL "\n\nUser: $user"  if ($user);
        print MAIL "\n\nRequest: http://$serv$link" if ($serv && $link);
        print MAIL "\n\n" . &stack_trace();
        print MAIL "\nEnvironment:\n";
        foreach my $var (sort keys %edat) {
            my $val = $edat{$var} || '';
            if (ref($val)) {
                foreach my $subval (@{$val}) {
                    printf(MAIL "%30s : %s\n", $var, $subval);
                    $var = "";
                }
            } else {
                printf(MAIL "%30s : %s\n", $var, $val);
            }
        }
        close MAIL;
        if ($adminMsg) {
            $adminMsg = "<div style='background-color:yellow'>$adminMsg</div>"
                if ($html);
            $msg .= $adminMsg;
        }
    }
    warn $msg;
    unless ($self && $self->no_environment()) {
        if ($html) {
            print $ev;
        } else {
            warn $ev;
        }
    }
    die;
}

sub verbose {
    my $self = shift;
    $self->{private_VERBOSE} = $_[0] if (defined $_[0]);
    return $self->{private_VERBOSE};
}

sub default_progress_delay {
    my $self = shift;
    $self->{private_PROG_DEFAULT_DELAY} = $_[0] if ($_[0]);
    return $self->{private_PROG_DEFAULT_DELAY};
}

sub initiate_progress {
    my $self = shift;
    return unless ($self->verbose);
    my $hash = {
        start => time,
        last  => time,
        total => $_[0],
        title => $_[1] || '',
        form  => $_[2] || "  %7d done, %4.1f min %s\n",
        wait  => $_[3] || $self->{private_PROG_DEFAULT_DELAY} || 180,
        done  => 0,
    };
    my $key = ++$self->{private_PROG_KEYS};
    $self->{private_PROG_DATA}{$key} = $hash;
    return $key;
}

sub track_progress {
    my $self = shift;
    return unless ($self->verbose);
    my ($key) = @_;
    my $hash  = $self->{private_PROG_DATA}{$key};
    $hash->{done}++;
    my $delta = time - $hash->{last};
    if ($delta > $hash->{wait}) {
        $hash->{last} = time;
        my $elapsed = $hash->{last} - $hash->{start};
        my $done    = $hash->{done};
        my $total   = $hash->{total};
        my $remain = '';
        if ($total) {
            $remain = $elapsed * ( $total - $done ) / $done;
            my $u = 'sec';
            if ($remain > 100) { $remain /= 60; $u = 'min' }
            if ($remain > 100) { $remain /= 60; $u = 'hr' }
            if ($remain >  50) { $remain /= 24; $u = 'day' }
            $remain = sprintf("[ %4.1f %s remain]", $remain, $u);
        }
        if ($hash->{title}) {
            my $msg = $hash->{title};
            $msg .= " - total $total" if ($total);
            warn "$msg\n";
            $hash->{title} = '';
        }
        warn sprintf( $hash->{form}, $done, $elapsed / 60, $remain);
    }
}

sub finish_progress {
    my $self = shift;
    return unless ($self->verbose);
    my ($key, $force) = @_;
    my $hash  = $self->{private_PROG_DATA}{$key};
    if ($force || $hash->{last} != $hash->{start} ) {
        if ($hash->{title}) {
            my $msg = $hash->{title};
            my $total = $hash->{total};
            $msg .= " - total $total" if ($total);
            warn "$msg\n";
            $hash->{title} = '';
        }
        my $elapsed = time - $hash->{start};
        warn sprintf("  Complete %d in %3.1f min\n", 
                     $hash->{done}, $elapsed / 60 );
    }
    delete $self->{private_PROG_DATA}{$key};
}

sub schema {
    my $self = shift;
    if (my $val = $_[0]) {
        $self->{private_SCHEMA} = $val;
    }
    return $self->{private_SCHEMA};
}

sub schema_owner {
    my $self = shift;
    if (defined $_[0]) {
        $self->{private_SCHEMA_OWNER} = $_[0];
    }
    return $self->{private_SCHEMA_OWNER};
}

sub report_owner {
    my $self = shift;
    unless (defined $self->{private_OWNER_VIA_DB}) {
        $self->benchstart;
        if ($self->driver() eq 'postgres') {
            $self->death("Have not implemented report_owner() for postgres");
        } else {
            my $sth = $self->prepare
                ( -name  => "Get Oracle owner",
                  -sql   => "SELECT user FROM dual",
                  -level => 2,);
            $self->dbh_msg("Getting owner information");
            my @own = $sth->get_array_for_field;
            $self->{private_OWNER_VIA_DB} = $own[0] || "";
        }
    }
    return $self->{private_OWNER_VIA_DB};
}

*build_all = *make_all;
sub make_all {
    my $self = shift;
    my $schema = $self->schema;
    unless ($schema) {
        $self->dbh_warn("You can not make_all() unless you define schema()");
        return;
    }
    $self->dbh_msg("Completely building database");
    my $cons = $self->{private_CONSTRAINTS} = {};
    my @views;
    my %needed = map { $_ => 1 } keys %{$schema};
    my %done;
    my $didSomething = 1;
    my $doingTables  = 1;
    while ($didSomething) {
        $didSomething = 0;
        my @need = keys %needed;
        last if ($#need == -1);
        my $tabCount = 0;
        map { $tabCount++ unless 
                  ($schema->{$_}{view} || $schema->{$_}{function}) } @need;
        $doingTables = 0 unless ($tabCount);
        foreach my $table (@need) {
            my $tabdat = $schema->{$table};
            my $delay  = 0;
            $delay++ if (($tabdat->{view} || $tabdat->{function})
                         && $doingTables);
            foreach my $req (@{$tabdat->{requires} || []}) {
                $delay++ unless ($done{uc($req)});
            }
            next if ($delay);
            $self->make_table( $table, $tabdat );
            $done{uc($table)}++;
            $didSomething++;
            delete $needed{$table};
        }
        if ($doingTables && !$didSomething) {
            # Nothing done, switch to processing views and functions
            $doingTables = 0;
            $didSomething = 1;
        }
    }

    my @notMade = keys %needed;
    $self->err("Failed to create some tables",
               @notMade) unless ($#notMade == -1);


    my %constraints;
    my $cOrder = 0;
    my $iinfo = $self->index_info( );
    while (my ($tab, $col) = each %{$cons->{PKEYS}}) {
        my $cname = lc(join('_', 'pk', $tab, $col));
        $cname  =~ s/[aeiou]//g if (length($cname) > 30);
        $cname  = substr($cname, 0, 30) if (length($cname) > 30);
        my $status = $self->index_status($cname, [$col], $tab );
        if ($status) {
            # Index already there
            $self->dbh_warn
                ("Attempt to create primary key '$cname' on $tab ($col) fails: ". $status) unless ($status eq '1');
            next;
        }
        my $ct  = $constraints{$cname} ||= 
        { o => ++$cOrder, n => $cname,  sql => [] };
        push @{$ct->{sql}}, "ALTER TABLE $tab ADD CONSTRAINT $cname ".
            "PRIMARY KEY ( $col )";
    }
    while (my ($srctab, $srccols) = each %{$cons->{FKEYS}}) {
        while (my ($srccol, $targdat) = each %{$srccols}) {
            my ($targtab, $targcol) = @{$targdat};
            my $cname = lc(join('_', 'fk', $srctab, $srccol,
                                $targtab, $targcol));
            $cname =~ s/[aeiou]//g if (length($cname) > 30);
            $cname = substr($cname, 0, 30) if (length($cname) > 30);
            my $ct  = $constraints{$cname} ||= 
            { o => ++$cOrder, n => $cname, sql => [] };
            push @{$ct->{sql}}, "ALTER TABLE $srctab ADD CONSTRAINT $cname ".
                "FOREIGN KEY ($srccol) REFERENCES $targtab ($targcol)";
        }
    }
    my @toDo = sort { $a->{o} <=> $b->{o} } values %constraints;
    foreach my $dat (@toDo) {
        my $cname = $dat->{n};
        my @sql = @{$dat->{sql}};
        if ($#sql == -1) {
            $self->err("Weird. No SQL found for constraint $cname");
            next;
        } elsif ($#sql != 0) {
            $self->err("Multiple constraints generated the same name ($cname)",
                       @sql);
            next;
        }
        if ($cname =~ /[^a-z0-9_]/) {
            $self->err("Skipping column constraint '$cname' with unusual characters");
            next;
        }
        my $sth = $self->prepare
            ( -sql => $sql[0],
              -ignore => [ 'already exists' ],
              -name => "Add DB constraint $cname" );
        $sth->execute();
    }
}

sub apply_constraints {
    my $self = shift;
    my $cons = $self->{CONSTRAINTS};
    return unless ($cons);

    # Not implemented currently ...
    
}

sub update_comments {
    my $self = shift;
    my $schema = $self->schema;
    unless ($schema) {
        $self->dbh_warn
            ("You can not update_comments() unless you define schema()");
        return;
    }
    while ( my ($table, $tabdat) = each %{$schema}) {
        $self->update_table_comments($tabdat);
    }
}

sub update_table_comments {
    my $self = shift;
    my $tabSchema = shift;
    return unless ($tabSchema);
    my %coms;
    my $table = $tabSchema->{name};
    $coms{"TABLE $table"} = $tabSchema->{com};
    if (my $cols = $tabSchema->{cols}) {
        foreach my $cdat (@{$cols}) {
            my ($name, $type, $com) = @{$cdat};
            $coms{"COLUMN $table.$name"} = $com;
        }
    }
    my @clauses = keys %coms;
    return if ($#clauses == -1);
    $self->dbh_msg("Updating ".scalar(@clauses)." comments for $table");
    foreach my $clause (@clauses) {
        if (my $com = $coms{$clause}) {
            $com =~ s/\'/\'\'/g;
            $self->execute( -sql => "COMMENT ON $clause IS '$com'",
                            -ignore => ' ');
        }
    }
}

sub make_view {
    my $self = shift;
    my ($name, $sql) = @_;
    
    $sql =~ s/^\s*[\n\r]*//;
    $sql =~ s/\s*[\n\r]*$//;

    $self->dbh_msg("Create / update view '$name'");
    my $sth = $self->prepare
        ( -name  => "Create view $name",
          -sql   => "CREATE OR REPLACE VIEW $name AS $sql",
          -level => 1,);
    $sth->execute();
}

sub make_function {
    my $self = shift;
    my ($name, $args, $sql, $rv, $lang) = @_;
    $args ||= [];
    my @abits;
    for (my $i = 0; $i < $#{$args}; $i += 2) {
        push @abits, $args->[$i] .' '. $args->[$i+1];
    }
    my $tag    = '$SQLTXT$';
    my $argTxt = join(', ', @abits) || "";
    my $build  = sprintf("CREATE OR REPLACE FUNCTION %s(%s)\n",
                         $name, $argTxt);
    if ($rv) {
        $build .= "  RETURNS $rv\n";
    }
    $build .= " AS\n";
    $sql =~ s/^\s*[\n\r]*//;
    $sql =~ s/\s*[\n\r]*$//;
    $sql .= ";" unless ($sql =~ /\;$/);
    $sql .= "\n" unless ($sql =~ /\n$/);
    $build .= join("\n", $tag, $sql, $tag)."\n";
    if ($lang) {
        $build .= "LANGUAGE $lang";
    }

    $self->dbh_msg("Create / update function '$name'");
    my $sth = $self->prepare
        ( -name  => "Create function $name",
          -sql   => $build,
          -level => 1,);
    return $sth->execute();
}

our $tempTabCounter = 0;
our $ownedTempTabs = {};
our $okTempColumn  = {};

sub _register_temp_table {
    my ($self, $tab, $cols) = @_;
    if ($tab) {
        $tab = lc($tab);
        my $ott = $ownedTempTabs->{$$} ||= {};
        if (my $already = $ott->{$tab}) {
            my $prior = join("\t", @{$already});
            my $curr  = join("\t", @{$cols});
            unless ($prior eq $curr) {
                $self->err("Temporary table '$tab' registered with different column order", "Now: ($already), Denied: ($curr)");
            }
        } else {
            # $self->err("[DEBUG]","TEMP TAB $tab = (".join("\t", @{$cols}).") PID=$$");
            $ott->{$tab} = $cols || [];
        }
    }
}

sub clear_temp_table {
    my $self = shift;
    my $tab = shift;
    return unless ($tab);
    $tab = lc($tab);
    my $ott = $ownedTempTabs->{$$} ||= {};
    my $cols = $ott->{$tab};
    $self->death("Refusing to drop temporary table '$tab' (case senistive)", 
                 "That table was not created by this API")
        unless ( $cols );
    # $self->err("[DEBUG]","CLEAR TEMP TAB $tab (".join("\t", @{$cols}).") PID=$$");
    my $drop = $self->prepare
        ( -name  => "Drop temporary table",
          -sql   => "DROP TABLE $tab",
          -level => 2,);
    $drop->execute();
    delete $ott->{$tab};
}

sub list_to_temp_table {
    # [0] = The table, as 1D or 2D array
    # [1] = Optional column data types
    # [2] = Optional column names
    # [3] = Optional table name
    # [4] = Optional index
    my $self = shift;
    $self->benchstart();
    my $list = shift;
    my $cReq = shift; # lc(shift || "char");
    my @creqs = ref($cReq) ? @{$cReq} : ($cReq);
    $creqs[0] ||= "char";

    # How many columns do we need to make?
    my $firstRow = $list->[0];
    my $colN = 0;
    if (!ref($firstRow)) {
        # This is just a simple 1D array of values
        # Normalize it to 2D
        $list = [ map { [$_] } @{$list} ];
        $colN = 1;
    } else {
        # Count up the number of columns in the first row
        $colN = $#{$firstRow} + 1;
    }
    my @cols;
    if (my $cnms = shift) {
        foreach my $col (@{$cnms}) {
            unless (defined $okTempColumn->{$col}) {
                # We may be doing this a lot, so cache the check
                $okTempColumn->{$col} = $col =~ /^[a-z0-9_]+$/i ? 1 : 0;
            }
            if ($okTempColumn->{$col}) {
                push @cols, $col;
            } else {
                $self->death("Temp table column names must be alphanumeric/underscore, not '$col'"); 
            }
        }
    } else {
        @cols = map { "col$_" } (1..$colN);
    }
    my $tab  = shift;
    if ($tab) {
        # For named tables, clear any prior ones
        my $ti = $self->table_info();
        if (exists $ti->{$tab}) {
            $self->death("Table '$tab' already exists in the database, you can not make a temporary table with that name"); 
        } elsif ($ownedTempTabs->{$$}{lc($tab)}) {
            # We already used this table
            $self->clear_temp_table( $tab );
        }
    } else {
        $tab = "FD_TEMP_TAB_".++$tempTabCounter;
    }

    # What types should the columns be?
    my $dbType = $self->driver();
    my @types;
    foreach my $cr (map { lc($_) } @creqs) {
        my $type = 
            $cr =~ /char/ ? 'varchar(4000)' : 
            $cr =~ /text/ ? ($dbType eq 'postgres' ? 'text' : 'varchar(4000)'):
            $cr =~ /int/ ? 'integer' :
            $cr =~ /real/ ? 'real' :
            $cr =~ /date/ ? 'date' :
            $self->death
            ("I can not make a temporary table with unrecognized type '$cr'");
        push @types, $type;
    }
    while ($#types != -1 && !$types[-1]) { pop @types; }
    for my $c (0..$#cols) { 
        if ($c > $#types || !$types[$c]) {
            $types[$c] = $types[$c-1];
        }
    }
    
    # Make the table
    my $make = $self->prepare
        ( -name  => "Create temporary table for list",
          -sql   => "CREATE GLOBAL TEMPORARY TABLE $tab ( ".
          join(', ', map { "$cols[$_] $types[$_]" } (0..$#cols)) .")",
          -level => 2,);
    $make->execute();
    $self->_register_temp_table( $tab, \@cols );
    # Add indices. I will admit to being a bit baffled by this, since I thought
    # these tables were held in RAM, but it seems they still want indices
    if (my $indReq = shift) {
        $self->_fast_add_index( $tab, $indReq );
    } else {
        $self->_fast_add_index( $tab, \@cols );
        # $self->add_index( $tab, { name => "${tab}_all", cols => \@cols });
        for my $c (1..$#cols) {
            my $col = $cols[$c];
            $self->_fast_add_index( $tab, [ $col ] );
            # $self->add_index( $tab, { name => "${tab}_$col", cols => [ $col ] });
        }
    }
    $self->insert_array( $tab, $list);
    $self->bench_end();
    return wantarray ? ($tab, @cols) : $tab;
}

sub select_to_temp_table {
    my $self = shift;
    $self->bench_start;
    my $colReq  = shift;
    my $sqlBody = shift;
    my $driver  = $self->driver();
    my $tab     = "FD_TEMP_TAB_".++$tempTabCounter;
    
    my @cols = map { lc($_) }
    (ref($colReq) ? @{$colReq} : split(/[\s\,]+/, $colReq));
    my $cTxt = join(', ', @cols);
    my $sql;
    if ($driver eq 'oracle') {
        $sql = "CREATE TEMPORARY TABLE $tab ON COMMIT PRESERVE ROWS AS ".
            "SELECT $cTxt FROM $sqlBody";
    } else {
        $sql = "SELECT $cTxt INTO TEMPORARY TABLE $tab FROM $sqlBody";
    }
    my $make = $self->prepare
        ( -name  => "Create temporary table for SELECT statement",
          -sql   => $sql,
          -level => 2,);
    $make->execute( @_ );
    $self->_register_temp_table( $tab, \@cols );
    $self->bench_end;
    return wantarray ? ($tab, @cols) : $tab;
}

sub make_table {
    my $self = shift;
    my ($table, $schema, $ign) = @_;
    my $driver = $self->driver();
    $table = lc($table);
    unless ($schema) {
        my $sdat = $self->schema;
        if ($sdat) {
            if (exists $sdat->{$table}) {
                $schema = $sdat->{$table};
            } else {
                $self->dbh_die
                    ("make_table('$table') could not find stored schema info");
            }
        } else {
            $self->dbh_die
                ("make_table('$table') requires schema information");
        }
    }

    if (my $vsql = $schema->{view}) {
        return $self->make_view( $table, $vsql );
    } elsif (my $fsql = $schema->{function}) {
        return $self->make_function( $table, $schema->{args}, $fsql,
                                     $schema->{retval}, $schema->{language});
    }
    my $dat  = $self->table_info( );
    # $self->table_info('tilfordc'); my @foo; while (my ($own, $th) = each %{$self->{private_TAB_INFO}}) { map { push @foo, "$own: $_" } sort keys %{$th} }; $self->death("Known tables:", @foo);
    my $cols = $schema->{cols};

    unless ($dat->{$table}) {
        # We need to create the table
        my $obj = uc($schema->{object} || 'TABLE');
        if (!$cols || $#{$cols} < 0) {
            $self->dbh_die
                ("make_table('$table') did not provide any columns");
            
        }
        # We will make the table blank:
        my $tsql = "CREATE $obj $table (";
        if ($driver eq 'oracle') {
            # Oracle can not tolerate an empty table, we need to add columns
            # here
            my @ctxt;
            foreach my $cdat (@{$cols || []}) {
                my ($name, $type, $com) = @{$cdat};
                if (my $typeAlias = $colTypeAliases->{$self->driver()}) {
                    $type = $typeAlias->{lc($type)} || $type;
                }
                push @ctxt, sprintf("  %10s %s", $name, $type);
            }
            $tsql .= "\n". join(",\n", @ctxt) . "\n";
        }
        $tsql .= ")";
        $tsql .= " " . $schema->{modifier} if ($schema->{modifier});
        $self->dbh_msg("Creating $obj $table");
        $self->begin_work;
        $self->execute( -sql    => $tsql, 
                        -name   => "Create $obj $table",
                        -ignore => $ign || $schema->{ignore}{make} || 0, );
        $self->commit;
        if ($driver eq 'oracle') {
            # Easiest just to rebuild the info hash
            $self->clear_table_info();
            $dat  = $self->table_info();
        } else {
            # Recursion occurs if you do not update the hash:
            $dat->{$table} = {};
        }
    }
    # Now create or update the columns:
    foreach my $coldat (@{$cols || []}) {
        $self->add_column( $table, $coldat );
    }
    
    if (my $inds = $schema->{index}) {
        while (my ($name, $inddat) = each %{$inds}) {
            $inddat->{name} = $name;
            $self->add_index( $table, $inddat );
        }
    }
    if (my $seqs = $schema->{sequence}) {
        while (my ($name, $dat) = each %{$seqs}) {
            $self->add_sequence( $name, $dat );
        }
    }
    if (my $fkeys = $schema->{fkey}) {
        while (my ($col, $targ) = each %{$fkeys}) {
            my @bits = split(/\./, lc($targ));
            $self->death("Foreign key constraint for $table.$col is malformed",
                         "It should be of form TARGETTABLE.TARGETCOL",
                         "Instead it is '$targ'") unless ($#bits == 1);
            $self->{private_CONSTRAINTS}{FKEYS}{$table}{lc($col)} = \@bits;
        }
    }
    if (my $pkey = $schema->{pkey}) {
        $self->{private_CONSTRAINTS}{PKEYS}{$table} = lc($pkey);
    }

    if ($driver eq 'oracle') {
        # Have to add in column comments now
        $self->update_table_comments( $schema );
    } elsif (my $com = $schema->{com}) {
        # Note the table comment:
        $com =~ s/\'/\'\'/g;
        $self->execute
            ("COMMENT ON TABLE $table IS '$com'") if ($com);
    }


    # Refresh information in DB
    $self->{private_TAB_INFO} = undef;
}


sub add_column {
    my $self = shift;
    my ($table, $column) = @_;
    $table     = lc($table);
    my $coldat = $column;
    unless (ref($coldat)) {
        $column = lc($column);
        my $sdat = $self->schema;
        if ($sdat) {
            if (exists $sdat->{$table}) {
                my $schema = $sdat->{$table};
                if (exists $schema->{cols}) {
                    foreach my $dat (@{$schema->{cols}}) {
                        if (lc($dat->[0]) eq $column) {
                            $coldat = $dat;
                            last;
                        }
                    }
                    unless ($coldat) {
                        $self->dbh_die
                            ("add_column('$table', '$column') could not find ".
                             "stored schema info column '$column'");
                    }
                } else {
                    $self->dbh_die
                        ("add_column('$table', '$column') could not find ".
                         "stored schema info for any columns");
                }
            } else {
                $self->dbh_die
                    ("add_column('$table', '$column') could not find stored ".
                     "schema info for '$table'");
            }
        } else {
            $self->dbh_die
                ("add_column('$table','$column') requires schema information");
        }
    }
    my ($name, $type, $com, $cons) = @{$coldat};
    $self->dbh_die( "add_column('$table') - column name is not defined")
        unless ($name);
    $name = lc($name);
    $self->dbh_die
        ( "add_column('$table','$name') - column type is not defined")
        unless ($type);

    my $dat = $self->table_info();
    unless ($dat->{$table}) {
        $self->make_table( $table );
    }
    $dat->{$table}{order} ||= [];
    if (!defined $dat->{$table}{cols}{$name}) {
        $self->dbh_msg("Adding column $table.$name");
        $self->begin_work;
        $self->execute("ALTER TABLE $table ADD $name $type");
        if ($com) {
            $com =~ s/\'/\'\'/g;
            $self->execute("COMMENT ON COLUMN $table.$name IS '$com'" );
        }
        $self->commit;
        my $index = $#{$dat->{$table}{order}} + 1;
        push @{$dat->{$table}{order}}, $name;
        $dat->{$table}{cols}{$name}  = $index;
        $dat->{$table}{cols}{$index} = $name;
    }
    # warn "$name, $type".$self->branch($cons);
    if ($cons) {
        my $driver = $self->driver();
        if ($cons->{SEQUENCE}) {
            # Want a default sequence on this column
            my $seq = $cons->{SEQUENCE};
            $self->add_sequence( $seq );
            my $seqSql;
            if ($driver =~ /^(postgres)$/) {
                $seqSql = "ALTER TABLE $table ALTER COLUMN $name SET DEFAULT nextval('$seq')";
            } else {
                my $t   = join('_', $table, $name, 'seq');
                my $tin = $self->trigger_info();
                unless ($tin->{$t}) {
                    $seqSql = <<EOF;
CREATE OR REPLACE TRIGGER $t
 BEFORE INSERT ON $table FOR EACH ROW WHEN (new.$name IS NULL) BEGIN
    SELECT $seq.nextval INTO :new.$name FROM dual;
 end;
EOF
                }
            }
            if ($seqSql) {
                eval {
                    $self->dbh_msg("Adding trigger $name on $table");
                    $self->execute( $seqSql );
                };
                if (my $err = $!) {
                    $self->dbh_msg($err);
                }
            }
        }

        if (my $def = $cons->{DEFAULT}) {
            my $qd = $self->quote($def);
            my $defSql;
            if ($driver =~ /^(postgres)$/) {
                $defSql = "ALTER TABLE $table ALTER COLUMN $name SET DEFAULT $qd";
            } else {
                my $t   = join('_', $table, $name, 'def');
                my $tin = $self->trigger_info();
                unless ($tin->{$t}) {
                    $defSql = <<EOF;
CREATE OR REPLACE TRIGGER $t
 BEFORE INSERT ON $table FOR EACH ROW WHEN (new.$name IS NULL) BEGIN
    SELECT $qd INTO :new.$name FROM dual;
 end;
EOF
                }
            }
            $self->execute( $defSql ) if ($defSql);
        }
    }

}

sub sequence_info {
    my $self = shift;
    unless ($self->{private_SEQ_INFO}) {
        $self->benchstart;
        $self->dbh_msg("Getting all database sequence information");
        if ($self->driver() eq 'postgres') {
            my $sth = $self->prepare
                ( -name  => "Get Postgres sequences",
                  -sql   => "SELECT relname from pg_statio_user_sequences",
                  -level => 2,);
            map {  $self->{private_SEQ_INFO}{lc($_)} = 1 }
            $sth->get_array_for_field;
        } else {
            my $sth;
            eval {
                $sth = $self->SUPER::table_info
                    (undef, $self->schema_owner, undef, 'SEQUENCE');
            };
            my $rows = $sth->fetchall_arrayref;
            foreach my $row (@{$rows}) {
                my ($qual, $owner, $table, $type) = @{$row};
                $table = lc($table);
                $self->{private_SEQ_INFO}{$table} = 1;
            }
        }
        $self->benchend;
    }
    return $self->{private_SEQ_INFO};
}

sub trigger_info {
    my $self = shift;
    unless ($self->{private_TRIGGER_INFO}) {
        $self->benchstart;
        $self->dbh_msg("Getting all database trigger information");
        my $info = $self->{private_TRIGGER_INFO} = {};
        if ($self->driver() eq 'postgres') {
            $self->death("Have not implemented trigger_info() for postgres");
        } else {
            # die $self->report_owner();
            my $sth = $self->prepare
                ( -name  => "Get Oracle triggers",
                  -sql   => "SELECT trigger_name, table_name, trigger_type, trigger_body FROM all_triggers",
                  -level => 2,);
            my $rows = $sth->selectall_arrayref;
            foreach my $row (@{$rows}) {
                my $tname = lc(shift @{$row});
                $info->{$tname} = $row;
            }
        }
        $self->benchend;
    }
    return $self->{private_TRIGGER_INFO};
}

sub add_sequence {
    my $self = shift;
    my ($name, $dat) = @_;
    $name   = lc($name);
    my $inf = $self->sequence_info();
    return 0 if ($inf->{$name});
    my $sql = "CREATE SEQUENCE $name";
    my $sth = $self->prepare( -sql => $sql,
                              -ignore => "name is already used",
                              -name => "Generate new sequence $name" );
    $self->dbh_msg("Adding sequence $name");
    $sth->execute;
    return 1;
}

sub add_index {
    my $self = shift;
    my ($table, $index ) = @_;
    $table     = lc($table);
    my $inddat = $index;
    unless (ref($inddat)) {
        $index   = lc($index);
        my $sdat = $self->schema;
        if ($sdat) {
            if (exists $sdat->{$table}) {
                my $schema = $sdat->{$table};
                if (exists $schema->{index} && 
                    exists $schema->{index}{$index}) {
                    $inddat = $schema->{index}{$index};
                    $inddat->{name} ||= $index;
                } else {
                    $self->dbh_die
                        ("add_index('$table', '$index') could not find ".
                         "stored schema info for '$index'");
                }
            } else {
                $self->dbh_die
                    ( "add_index('$table', '$index') could not find stored ".
                      "schema info for '$table'");
            }
        } else {
            $self->dbh_die
                ( "add_index('$table','$index') requires schema information");
        }
    }
    my $name = $inddat->{name};
    $self->dbh_die( "add_index('$table') - index name is not defined")
        unless ($name);
    $name = lc($name);
    my $cols = $inddat->{cols};
    $self->dbh_die("add_index('$table','$name') - index columns not defined")
        unless ($cols);
    my @list   = map { lc($_) } (ref($cols) ? @{$cols} : ($cols));
    my @cleaned;
    foreach my $term (@list) {
        my $val = $term;
        if ($val =~ /\((\S+)\)/) {
            $val = $1;
        }
        push @cleaned, $val;
    }
    my $passed = join(",", @list);
    my $status = $self->index_status($name, \@list, $table );
    if ($status) {
        # Index already there
        $self->dbh_warn
            ("Attempt to create index '$name' on $table ($passed) fails: ".
             $status) unless ($status eq '1');
        return;
    }

    # Create the index
    my $sql = sprintf
        ("CREATE %sINDEX %s ON %s (%s)", $inddat->{unique} ? 'UNIQUE ' : '',
         $name, $table, join(', ', @list));

    my $sth = $self->prepare( -sql => $sql,
                              -ignore => 'such column list already indexed',
                              -name => "Generate index for $table" );
    $self->dbh_msg("Adding index $name ON $table ($passed)");
    eval {
        $sth->execute;
    };
    $self->dbh_msg($!) if ($!);
    my $iinfo = $self->index_info();
    $iinfo->{$name} = {
        table => $table,
        cols  => \@cleaned,
    };
    if (my $com = $inddat->{com}) {
        # no comments on indices?
    }
}

sub _fast_add_index {
    my $self = shift;
    $self->bench_start();
    my ($table, $cols) = @_;
    my $name = "QUIK_IND_".(++$self->{private_FAST_IND_COUNTER});
    my $sql  = sprintf
        ("CREATE INDEX %s ON %s (%s)", $name, $table, join(', ', @{$cols}));
    my $sth = $self->prepare( -sql => $sql,
                              -ignore => 'such column list already indexed',
                              -name => "Generate index for $table" );
    eval { $sth->execute; };
    $self->bench_end();
    return $!;
}

sub index_status {
    my $self = shift;
    my ($name, $cols, $table) = @_;
    return undef unless ($name);
    my $iinfo = $self->index_info( );
    unless (exists $iinfo->{$name}) {
        return 0;
    }
    # There is already an index with this name
    my $nowTab = $iinfo->{$name}{table};
    return "Index $name exists but is assigned to table $nowTab"
        if ($table && uc($nowTab) ne uc($table));
    if ($cols) {
        my $ctxt   = lc(join(",", @{$cols}) || '');
        my $ntxt   = lc(join(",", @{$iinfo->{$name}{cols}}));
        unless ($ctxt eq $ntxt || $ntxt =~ /function\(\)/ ||
                ($ctxt =~ /\(/ && $ntxt =~ /^sys_.+\$$/)) {
            return "Index exists with column order ($ntxt), not ($ctxt)";
        }
    }
    return 1;
}

sub update_lock_method {
    my $self = shift;
    if (defined $_[0]) {
        $self->{private_LOCK_METH} = $_[0];
    }
    return $self->{private_LOCK_METH};
}

sub update_unlock_method {
    my $self = shift;
    if (defined $_[0]) {
        $self->{private_UNLOCK_METH} = $_[0];
    }
    return $self->{private_UNLOCK_METH};    
}

sub update_via_array_insert {
    my $self = shift;
    if (defined $_[0]) {
        my $nv     = $_[0] =~ /^\d+$/ ? $_[0] : 1;
        my $cb     = $self->_eval_can('insert_array');
        $nv        = 0 unless ($cb);
        $self->{private_UPDATE_BY_INSERT} = $nv;
    }
    return $self->{private_UPDATE_BY_INSERT};
}

sub update_rows {
    my $self = shift;
    my ($table, $rows, $dcols) = @_;
    return if (!$rows || $#{$rows} < 0);
    my ($insSth, $cols, $stndDel, $lenCheck) = $self->get_update_info($table);
    my $onow    = $self->oracle_now();
    my $uai     = $self->update_via_array_insert();
    my $lockSth = $self->prepare
        (-sql   => "LOCK TABLE $table IN EXCLUSIVE MODE",
         -level => 4,
         -name  => "Lock table for deletes" );
    $dcols ||= $stndDel;
    my %deleteGroups;

    my $progKey = $self->initiate_progress
        ( $#{$rows} + 1, "Parsing rows for '$table'");
    $self->benchstart('Parse Rows');
    my %existingUpdate;
    foreach my $row (@{$rows}) {
        if (my $eu = $row->{updated}) {
            $existingUpdate{$eu} = 1;
        } else {
            # Automatically set the date to now
            $row->{updated} = $onow;
            $existingUpdate{$onow} = 1;
        }
    }
    my ($lowDate, $secondDate) = sort keys %existingUpdate;
    if ($lowDate && $secondDate) {
        # There are at least two different dates in the data
        # We want all dates to be identical; use the oldest one
        map { $_->{updated} = $lowDate } @{$rows};
    }
    foreach my $row (@{$rows}) {
        while (my ($cname, $clen) = each %{$lenCheck}) {
            # Auto-truncate some columns (if designated)
            $row->{$cname} = substr($row->{$cname}, 0, $clen)
                if ($row->{$cname});
        }

        my ($key, $delbinds, $delSth) = $dcols ? 
            $self->_delete_with_nulls($table, $dcols, $row) : ('');
        my $data = $deleteGroups{ $key } ||= {
            delsth  => $delSth,
            delbind => $delbinds,
            insbind => [],
        };
        push @{$data->{insbind}}, [ map { $row->{$_} } @{$cols} ];
        $self->track_progress($progKey);
    }
    $self->benchend('Parse Rows');
    $self->finish_progress($progKey);
    my ($custLock, $custUnlock) = ( $self->update_lock_method(),
                                    $self->update_unlock_method() );
    my $isLocked = 0;
    $self->begin_work;
    foreach my $key (keys %deleteGroups) {
        my $data =  $deleteGroups{ $key };
        my $locKey;
        if (my $delSth = $data->{delsth}) {
            $self->benchstart('Delete Rows');
            # We need to delete any pre-existing rows - lock the table
            if ($custLock) {
                $self->bench_start('Get Lock');
                $locKey = &{$custLock}( $self, $table, $key);
                $self->bench_end('Get Lock');
            } elsif (!$isLocked++) {
                $lockSth->execute();
            }
            $delSth->execute( @{$data->{delbind}} );
            $self->benchend('Delete Rows');
        }
        $self->benchstart('Load Rows');
        if ($uai && $#{$data->{insbind}} + 1 >= $uai) {
            # printf("<COPY> Rows for (%s) = %d\n", $key, $#{$data->{insbind}} + 1) if ($#{$data->{insbind}} > 500);
            $self->insert_array( $table, $data->{insbind} );
        } else {
            foreach my $binds (@{$data->{insbind}}) {
                $insSth->execute( @{$binds} );
            }
        }
        $self->benchend('Load Rows');
        if ($locKey && $custUnlock) {
            $self->bench_start('Release Lock');
            &{$custUnlock}( $self, $locKey );
            $self->bench_end('Release Lock');
        }
    }
    $self->commit;
}

sub oracle_now {
    my $self = shift;
    my ($s,$m,$h, $day, $mon, $year) = localtime(time);
    return sprintf("%04d-%02d-%02d %02d:%02d:%02d",
                   $year + 1900, $mon+1, $day, $h, $m, $s);
}

sub date_format {
    return 'yyyy-mm-dd hh24:mi:ss';
}

sub oracle_start {
    my $self = shift;
    if ($_[0] || !$self->{private_ORA_START}) {
        my $nv = $_[0];
        unless ($nv) {
            my $dfm = $self->date_format;
            my $sql = $self->driver() eq 'postgres' ?
                "SELECT to_char( now(), '$dfm')" :
                "SELECT to_char(sysdate, '$dfm') FROM dual";
            my $sth = $self->prepare
                ( -sql    => $sql,
                  -level  => 4,
                  -name   => "Get a static Oracle 'start time'");
            $nv = $sth->get_single_value();
        }
        $self->{private_ORA_START} = $nv;
    }
    return $self->{private_ORA_START};
}

sub get_update_info {
    my $self = shift;
    my ($table) = @_;
    $table = lc($table);
    unless ($self->{private_UPDATE_STH}{$table}) {
        my $dat    = $self->table_info();
        my $schema = $self->schema;
        $self->dbh_die("Can not get_update_info('$table') - no such table")
            unless ($dat->{$table});
        my $sinfo = $schema->{$table};
        my @cols  = @{$dat->{$table}{order}};
        my %upfrm = %{$sinfo->{upform} || {}};
        $upfrm{updated} = "to_timestamp(?, '".$self->date_format()."')";
        my @binds =  map { $upfrm{$_} || '?' } @cols;
        my $isql  = sprintf
            ("INSERT INTO %s ( %s ) VALUES ( %s )", $table, join(", ", @cols),
             join(", ", @binds));
        my $isth = $self->prepare
            ( -sql    => $isql,
              -level  => 4,
              -ignore => $sinfo->{ignore}{insert},
              -name   => "Insert function for $table");
        my @dcols = ($schema && $sinfo)  ? @{$sinfo->{update} || []} : ();
        my %lenCheck;
        if (my $carr = $sinfo->{checklength}) {
            my $ti   = $self->table_info() || {};
            $ti      = $ti->{$table} || {};
            foreach my $column (map { lc($_) } @{$carr}) {
                my $type = $ti->{coltype}{$column} || '-unknown-';
                unless ($type =~ /char|text/i) {
                    $self->err("Length restriction placed on non-character ".
                               "column $table.$column ($type)");
                    next;
                }
                my $width = $ti->{colsize}{$column};
                unless ($width) {
                    $self->err("Unable to establish column width for ".
                               "table.$column ($type)");
                    next;
                }
                $lenCheck{$column} = $width;
            }
        }
        $self->{private_UPDATE_STH}{$table} =
            [ $isth, \@cols, $#dcols > -1 ? \@dcols : undef, \%lenCheck ];
    }
    return @{$self->{private_UPDATE_STH}{$table}};
}

sub _delete_with_nulls {
    my $self = shift;
    my ($table, $dcols, $row) = @_;
    return ('') unless ($dcols && $#{$dcols} > -1);
    $table = lc($table);
    my ($key, $nullkey, @binds) = ("","");
    my $colKey = join(',', @{$dcols});
    foreach my $col (@{$dcols}) {
        my $val = $row->{$col};
        if (defined $val && $val ne '') {
            $key     .= "$val\t";
            $nullkey .= "+";
            push @binds, $val;
        } else {
            $key     .= "\t";
            $nullkey .= "-";
        }
    }
    my $tabKey = "$table $colKey = $nullkey";
    unless ($self->{private_DELETE_STH}{$tabKey}) {
        my $schema = $self->schema;
        my @dfrmt;
        my @nk = split('', $nullkey);
        for my $c (0..$#nk) {
            my $col = $dcols->[$c];
            if ($nk[$c] eq '-') {
                # Null position
                push @dfrmt, "$col IS NULL";
            } else {
                my $frm = $schema->{$table}{upform}{$col} || '?';
                push @dfrmt, "$col = $frm";
            }
        }
        my $dsql = sprintf("DELETE FROM %s WHERE %s", $table,
                           join(" AND ", @dfrmt));
        my $sth = $self->prepare
            ( -sql   => $dsql,
              -level => 4,
              -name  => "Delete from $table ($colKey) with Nulls '$nullkey'");
        $self->{private_DELETE_STH}{$tabKey} = $sth;
    }
    return ($key, \@binds, $self->{private_DELETE_STH}{$tabKey});
}


sub begin_work {
    my $self = shift;
    unless ($self->{private_NESTED_BEGIN}++) {
        $self->{private_BEGIN_BLOCK_ERRS} = [];
        eval {
            $self->SUPER::begin_work( @_ );
        };
        $self->verbose_error( $self->error( ) ) if ($self->dbi_err);
    }
}

*end_work = \&commit;
sub commit {
    my $self = shift;
    $self->{private_NESTED_BEGIN}--;
    if ($self->{private_NESTED_BEGIN} < 0) {
        warn "You have called commit more times than begin_work";
        $self->{private_NESTED_BEGIN} = 0;
    } elsif ($self->{private_NESTED_BEGIN} == 0) {
        if ($self->{private_BEGIN_WARN}) {
            warn "/* COMMIT */\n\n";
            $self->{private_BEGIN_WARN} = 0;
        }
        eval {
            $self->SUPER::commit( @_ );
        };
        $self->verbose_error( $self->error( ) ) if ($self->dbi_err);
        $self->{private_BEGIN_BLOCK_ERRS} = undef;
    }
}

sub execute {
    my $self = shift;
    my ($sth, $args)  = $self->prepare( @_ );
    return $sth->execute( $args->{BIND} ? @{$args->{BIND}} : () );
}

sub prepare {
    my $self = shift;
    my ($sql, $args) = $self->standardize_sql( @_ );
    $self->verbose_error("/* prepare() : No SQL Provided */\n") unless ($sql);
    my $sth;
    if (my $at = $args->{ATTR}) { while (my ($tag,$val) = each %{$at}) { warn " $tag = $val\n" } }
    eval {
        $sth = $self->SUPER::prepare( $sql, $args->{ATTR} );
    };
    $self->verbose_error( $self->error( ) ) if ($self->dbi_err &&
                                                !$args->{NOPREPERR});
    $sth->name( $args->{NAME} );
    $sth->level( $args->{LEVEL} );
    $sth->ignore_errors( $args->{IGNORE} ) if ($args->{IGNORE});
    return wantarray ? ($sth, $args) : $sth;
}

sub prepare_cached {
    my $self = shift;
    my ($sql, $args) = $self->standardize_sql( @_ );
    $self->{private_CACHED_STH}{$sql} = $self->prepare( %{$args} )
        unless ($self->{private_CACHED_STH}{$sql});
    my $sth = $self->{private_CACHED_STH}{$sql};
    return wantarray ? ($sth, $args) : $sth;
}

sub error {
    my $self = shift;
    my $sql  = $self->{Statement} || '';
    my $txt  = $sql ? $self->pretty_print($sql) : "";
    my $estr = $self->errstr;
    $txt .= sprintf("\nDBH Failure: ERR %d : %s\n\n%s  ", $self->dbi_err,
                    $estr, &stack_trace());
   # print "-- $txt\nNESTING: ".($self->{private_NESTED_BEGIN} || 'no')."\n\n";
    return $txt;
}

sub stack_trace {
    my @history;
    my $hist = 2;
    while (1) {
        my ($pack, $file, $j4, $subname) = caller($hist);
        last unless ($subname);
        my ($j1, $j2, $line) = caller($hist-1);
        push @history, sprintf("  %50s : %d\n", $subname, $line);
        $hist++;
    }
    return join('', @history) || '-No stack trace-';
}

# Will die if the name does not exist
sub named_sth {
    my $self = shift;
    my $name = $_[0];
    my $ucname = uc($name);
    if (!$self->{private_NAMED_STHS}{$ucname}) {
        my $sql   = $_[1];
        my $level = $_[2];
        if (my $pend = $self->{private_PENDING_STHS}{$ucname}) {
            # The STH has not been prepared yet, but data is pending:
            ($name, $sql, $level) = @{$self->{private_PENDING_STHS}{$ucname}};
        } elsif (!$sql) {
            $self->dbh_die("No SQL provided for named STH '$name'");
        }
        # The user is defining the statment handle
        my $sth   = $self->prepare( -sql    => $sql,
                                    -name   => $name,
                                    -level  => $level, );
        $self->verbose_error( $self->error( ) ) if ($self->dbi_err);

        $self->{private_NAMED_STHS}{$ucname} = $sth;
        if (defined $level && $self->debug_level >= $level) {
            $sth->pretty_print( );
        }
    }
    unless ($self->{private_NAMED_STHS}{$ucname}) {
        $self->dbh_die("Attempt to recover non-existant statement handle ".
                     "for name '$name'");
    }
    return $self->{private_NAMED_STHS}{$ucname};
}

sub note_named_sth {
    my $self = shift;
    my $name = uc($_[0]);
    die "No SQL provided for pending STH '$_[0]'\n  " unless ($_[1]);
    $self->{private_PENDING_STHS}{$name} = [ $_[0], $_[1], $_[2] ];
}

sub purge_sths {
    my $self = shift;
    foreach my $name (keys %{$self->{private_PENDING_STHS}}) {
        delete $self->{private_NAMED_STHS}{ uc($name) };
    }
}

sub clear_table_info {
    my $self = shift;
    $self->{private_TAB_INFO} = undef;
}

sub table_info {
    my $self = shift;
    my $owner = shift;
    # $owner    = $self->schema_owner() unless (defined $owner);
    $owner    = uc($owner || "");
    my $tinfo = $self->{private_TAB_INFO}{$owner};
    unless ($tinfo) {
        $self->benchstart;
        # $self->dbh_msg("Recover overall table information");
        $tinfo = $self->{private_TAB_INFO}{$owner} = {};
        my $isPg = $self->driver() eq 'postgres' ? 1 : 0;
        my @sths;
        eval {
            my $sch   = 'public';
            my @types = qw(TABLE VIEW);
            unless ($isPg) {
                $sch = $owner ? $owner : $self->schema_owner;
                push @types, qw(ALIAS SYNONYM);
            }
            foreach my $type (@types) {
                if (my $tsth = $self->SUPER::table_info
                    (undef, $sch, undef, $type)) {
                    push @sths, $tsth;
                }
                $self->verbose_error( $self->error( ) ) if ($self->dbi_err);
            }
        };
        $self->death("No queries found to determine table information for this database") if ($#sths == -1);
        my $schema = $self->schema;
        foreach my $sth (@sths) {
            my $rows = $sth->fetchall_arrayref;
            foreach my $row (@{$rows}) {
                my ($qual, $own, $table, $type) = @{$row};
                $table = lc($table);
                # Ignore wierd temp tables:
                # warn $table if ($table =~ /calc/);
                next if ($table =~ /\$/);
                next if ($schema && ! exists $schema->{$table});
                $tinfo->{$table}{table} = $table;
                $tinfo->{$table}{owner} = $own;
                $tinfo->{$table}{type}  = $type;
                # warn "TABLE DEBUG [$own]: ".join(' | ', map { defined $_ ? $_ : '---'} $qual, $own, $table, $type)."\n";
                my $sqT = $own ? "$own.$table" : $table;
                my $colsth;
                eval {
                    $colsth = $self->prepare
                        ( -name  => "Bogus SQL, used to make STH to get column names",
                          -sql   => "SELECT * FROM $sqT where 1=0",
                          -level => 4,
                          -nopreperr => 1);
                };
                if ($self->dbi_err) {
                    $self->err("Failed to recover columns for '$own:$table [$type]") if ($self->verbose());
                    next;
                }
                $colsth->execute();
                # my %foo; map { $foo{$_} = $colsth->{$_} } qw(NUM_OF_FIELDS NUM_OF_PARAMS NAME NAME_hash TYPE PRECISION SCALE NULLABLE CursorName Database ParamValues ParamArrays ParamTypes Statement RowsInCache); die $self->branch(\%foo);
                my @cols  = map { lc($_) } @{$colsth->{NAME}};
                my @tNums = @{$colsth->{TYPE} || []};
                my @ctyps;
                for my $c (0..$#cols) {
                    my $tyNum = $tNums[$c];
                    my $h     = defined $tyNum ? $self->type_info($tyNum) : undef;
                    my $ty = "";
                    if ($h) {
                        unless ($ty = lc(scalar $h->{TYPE_NAME} || "")) {
                            $self->err("Column type number does not map to name",
                                       "$table.$cols[$c] [$tyNum]");
                        }
                    } else {
                        # It looks like Postgres 'text' is not recognized,
                        # Mapping to $tyNum == -1
                        # Try to scavenge values from the schema
                        my $schema = $self->schema();
                        my $found;
                        if (exists $schema->{$table}) {
                            if (exists $schema->{$table}{cols}) {
                                foreach my $dat (@{$schema->{$table}{cols}}) {
                                    if (lc($dat->[0]) eq $cols[$c]) {
                                        $found = lc($dat->[1]);
                                    }
                                }
                            }
                        }
                        if (!$found) {
                            if ($tyNum == 40) {
                                # Appears to be Oracle clob?
                                $found = 'clob';
                            } elsif ($tyNum == 30) {
                                # Appears to be Oracle blob?
                                $found = 'blob';
                            }
                        }
                        
                        $ty = $found || "";
                        $self->err("Failed to find column type info for $type column",
                                   "$table.$cols[$c] [$tyNum]") unless ($ty || lc($type) eq 'view');
                    }
                    $ctyps[$c] = $ty;
                }
                my @csz   = @{$colsth->{PRECISION}};
                # 4 bytes added to the 'precision'. Maybe in Oracle too?
                map { $_ -= 4 if ($_) } @csz if ($isPg);
                $tinfo->{$table}{colsize} = { map {
                    $cols[$_] => $csz[$_] } (0..$#cols) };
                my %coldat;
                for my $index (0..$#cols) {
                    my $cname       = $cols[$index];
                    $coldat{$cname} = $index;
                    $coldat{$index} = $cname;
                }
                # Index -1 stores a list of the columns in proper order
                $tinfo->{$table}{order} = \@cols;
                $tinfo->{$table}{cols}  = \%coldat;
                $tinfo->{$table}{coltype} = { map { $cols[$_] => $ctyps[$_] } (0..$#cols)};
            }
        }
        unless ($isPg) {
            my $tabCom = $self->prepare
                ( -name  => "Get Oracle table comments",
                  -sql   => "SELECT table_name, comments FROM all_tab_comments WHERE upper(owner) = ?",
                  -level => 4,);
            $tabCom->execute(uc($owner));
            my $rows = $tabCom->selectall_arrayref;
            foreach my $row (@{$rows}) {
                if (my $com = $row->[1]) {
                    $tinfo->{lc($row->[0])}{comment} = $com;
                }
            }
            my $colCom =  $self->prepare
                ( -name  => "Get Oracle column comments",
                  -sql   => "SELECT table_name, column_name, comments FROM all_col_comments WHERE upper(owner) = ?",
                  -level => 4,);
            $colCom->execute(uc($owner));
            $rows = $colCom->selectall_arrayref;
            foreach my $row (@{$rows}) {
                if (my $com = $row->[2]) {
                    $tinfo->{lc($row->[0])}{colcom}{lc($row->[1])} = $com;
                }
            }
        }
        my @allTabs = keys %{$tinfo};
        $self->err("No tables were detected for owner '$owner'")
            if ($#allTabs == -1);
        $self->benchend;
    }
    return $tinfo;
}

sub table_info_text {
    my $self = shift;
    my $tinfo = $self->table_info(@_);
    my $txt = "";
    foreach my $tab ($self->all_tables(@_)) {
        $txt .= "*** $tab ***";
        my $ti = $tinfo->{$tab};
        if (my $com = $ti->{comment}) {
            $txt .= " $com";
        }
        $txt .= "\n";
        my @coldat;
        foreach my $col (@{$ti->{order} || []}) {
            push @coldat, [$col, $ti->{coltype}{$col}, $ti->{colsize}{$col}, $ti->{colcom}{$col}];
        }
        $txt .= $self->array_to_text( \@coldat, ["Column", "Type", "Size","Comment"]);
    }
    return $txt;
}

sub array_to_text {
    my $self = shift;
    my $arr  = shift || [];
    my $head = shift;
    # Calculate column widths
    my @w;
    foreach my $row (@{$arr}, $head || []) {
        map { push @{$w[$_]}, length( defined $row->[$_] ? $row->[$_] : "") } (0..$#{$row});
    }
    my $frm = "|";
    my $hr  = "+";
    for my $i (0..$#w) {
        my ($wid) = sort { $b <=> $a } @{$w[$i]};
        $frm .= " %-${wid}s |";
        $hr  .= ('-' x ($wid + 2))."+";
    }
    $frm .= "\n";
    $hr  .= "\n";
    my $txt = $hr;
    if ($head) {
        $txt .= sprintf($frm, map { defined $_ ? $_ : "" } map { $head->[$_] } (0..$#w));
        $txt .= $hr;
    }
    foreach my $row (@{$arr}) {
        $txt .= sprintf($frm, map { defined $_ ? $_ : "" } map { $row->[$_] } (0..$#w));
    }
    $txt .= $hr;
    return $txt;
}

sub all_tables {
    my $self  = shift;
    my $tinfo = $self->table_info(@_);
    return sort keys %{$tinfo};
}

sub column_order {
    my $self  = shift;
    my ($tab) = @_;
    die "You must supply table name for column_order()"
        unless ($tab);
    $tab = lc($tab);
    if (my $co = $ownedTempTabs->{$$}{$tab}) {
        return @{$co};
    }
    my $owner;
    if ($tab =~ /^([^\.]+)\.(.+)$/) {
        ($owner, $tab) = ($1, $2);
    }
    my $tinfo = $self->table_info( $owner );
    #warn "($owner, $tab)"; $debug->branch($tinfo);
    die "Failed to find table for column_order('$tab')"
        unless (exists $tinfo->{$tab});
    return @{$tinfo->{$tab}{order}};
}

sub index_info {
    my $self = shift;
    unless ($self->{private_IDX_INFO}) {
        $self->dbh_msg("Recover index information");
        $self->benchstart;
        my %data;
        if ($self->driver() eq 'postgres') {
            my $dat  = $self->table_info();
            my $sth = $self->prepare
                ( -name => "Get index information for Postgres database",
                  -sql  => 
                  "SELECT s.indexrelname, s.relname, i.indkey".
                  "  FROM pg_index i, pg_stat_all_indexes s".
                  " WHERE i.indexrelid = s.indexrelid".
                  "   AND s.schemaname = 'public'",
                  );
            $sth->execute(  );
            my $rows = $sth->fetchall_arrayref(  );
            foreach my $row (@{$rows}) {
                my ($iname, $tname, $cpos) = @{$row};
                $iname = lc($iname);
                $data{$iname} ||= {
                    table => lc($tname),
                    cols  => [],
                };
                my $corder = $dat->{lc($tname)}{order};
                next unless ($corder);
                foreach my $ind (split(' ', $cpos)) {
                    my $cname = 'UNKNOWN';
                    if ($ind) {
                        $ind--;
                        $cname = $corder->[$ind] unless ($ind > $#{$corder});
                    } else {
                        $cname = 'FUNCTION()';
                    }
                    push @{$data{$iname}{cols}}, $cname;
                }
            }
        } else {
            my $sth = $self->prepare
                ( -name => "Get index information for database",
                  -sql  => 
                  "SELECT index_name, table_name, column_name, column_position".
                  "  FROM USER_IND_COLUMNS ORDER BY index_name, column_position",
                  );
            $sth->execute(  );
            my $rows = $sth->fetchall_arrayref(  );
            foreach my $row (@{$rows}) {
                my ($iname, $tname, $cname, $pos) = @{$row};
                $iname = lc($iname);
                $data{$iname} ||= {
                    table => lc($tname),
                    cols  => [],
                };
                $data{$iname}{cols}[ $pos -1 ] = lc($cname);
            }
        }
        $self->{private_IDX_INFO} = \%data;
        $self->benchend;
    }
    return $self->{private_IDX_INFO};
}

sub dbh_msg {
    my $self = shift;
    return unless ($self->verbose);
    $self->msg("[DB]", @_);
}

sub dbh_warn {
    my $self = shift;
    my ($msg) = @_;
    warn "$msg\n" . &stack_trace();
}

sub dbh_die {
    my $self = shift;
    my ($msg) = @_;
    die "$msg\n" . &stack_trace()."\n  FATAL!\n";
}

sub standardize_sql {
    my $self = shift;
    # DBI::_do_selectrow was passing undef at the end of a prepare call
    while (! defined $_[-1]) { pop @_ };
    unshift @_, '-sql' if ($#_ == 0);
    # print "standardize_sql(".join(', ', @_).");\n";
    my $args  = $self->parseparams( @_ );
    my $sql   = $args->{SQL};
    my $debug = $args->{DUMPSQL} || $self->debug_level;
    if (my $limit = $args->{LIMIT} || $args->{ROWNUM}) {
        my $lsyn = $self->limit_syntax();
        $lsyn =~ s/AND /WHERE / if ($sql !~ /WHERE/);
        my $lsql = sprintf(" %s %d", $lsyn, $limit);
        if ($sql =~ /order by/i && $self->driver eq 'oracle') {
            $sql =~ s/order by/$lsql ORDER BY/i;
        } else {
            $sql .= $lsql;
        }
    }
    $self->pretty_print( $sql, $args) 
        if ($debug && (!$args->{LEVEL} || $debug >= $args->{LEVEL}));
    return wantarray ? ($sql, $args) : $sql;   
}

sub pretty_print {
    my $self = shift;
    my ($sql, $args) = @_;
    $args ||= {};

    my @history;
    for my $hist (2..4) {
	my @f = split "::", (caller($hist))[3] || ""; # Calling funciton
	push @history, $f[$#f] if ($f[$#f]);
    }
    my $callHist = (join " < ", @history) || "";

    return "/* NULL SQL QUERY : $callHist */\n" unless ($sql);
    my $sqlcom = "/* %s : %s */\n";
    $sql =~ s/[\s\n\t]+/ /g;
    $sql = " $sql ";
    my $maxtag = 12;
    my @tags = ("CREATE TABLE", "CREATE OR REPLACE FUNCTION", "SELECT", 
                "AS SELECT", "INSERT INTO",
                "UPDATE", "LOCK", 'LEFT JOIN', "COMMENT ON", "DECLARE",
                'THEN CASE WHEN', 'ELSE CASE WHEN', 'CASE WHEN',
                'THEN', 'ELSE', 'END', "ANY",
		["UNION", "\nUNION\t\n"], "LIMIT", "VALUES", "COPY",
		"DELETE FROM", "FROM", "WHERE", "OR", "AND", "ORDER BY",
                "GROUP BY", "HAVING", 'ON');
    foreach my $set (@tags) {
	# Case sensitive - the SQL should have keywords in caps
	my ($tag, $out);
	if (ref($set)) {
	    ($tag, $out) = @{$set};
	} else {
	    $tag = $set;
	    $out = "\n$tag\t";
	}
	$sql =~ s/[\n ]+$tag[\n ]+/$out/g;
    }
    $sql =~ s/\([\n ]*/\(/g;

    
    $sql =~ s/^ //; $sql =~ s/ $//;
    $sql .= ';' unless ($sql =~ /\;$/);
    $sql =~ s/\;/\;\n/g;
    my @lines = split("\n", $sql);
    my @newlines;
    my $maxline = 60;
    my $indent = 0; my @pad = ("");
    while ($#lines > -1) {
	my $line = shift @lines;
	next if ($line =~ /^\s*$/);
	# Wrap long lines:
	if (length($line) > $maxline) {
            # Temporarily mask spaces inside quotes:
            my $iloop = 0;
            while ($line =~ /(\'[^\']*?) ([^\']*?\')/) {
                $line =~ s/(\'[^\']*?) ([^\']*?\')/$1\n$2/;
                last if (++$iloop > 500);
            }
	    my $pos = rindex($line, " ", $maxline - 2);
            # Reset spaces that were inside quotes:
            $line =~ s/\n/ /g;
	    if ($pos > 0) {
		my $tail = " \t" . substr($line, $pos+1);
		unshift @lines, $tail;
		$line = substr($line, 0, $pos);
	    }
	}
        my @bits = split("\t", $line);
        my $pre = shift @bits;
        my $pro = join(" ", @bits);
	# Manage parentheses indenting:
	if ($indent > 0) {
	    ($pre, $pro) = ("", $pad[-1] . "$pre $pro");
	}
        # print "<pre>($pre, $pro)</pre>";
	my $newindent = $indent;
	$newindent += ( $line =~ tr/\(/\(/ );
	$newindent -= ( $line =~ tr/\)/\)/ );
	if ($newindent > $indent) {
	    my $ppos = index($pro, "(");
	    push @pad, " " x ($ppos+1);
	} elsif ($newindent < $indent) {
	    pop @pad;
	}
	$indent = $newindent;
	my $pline = sprintf("%".$maxtag."s %s", $pre, $pro);
	push @newlines, $pline;
    }
    $sql = join("\n", @newlines);
    my $header = "";
    if (my $bvs = $args->{BIND}) {
        # $header .= sprintf($sqlcom, 'Bind variables' ,join(", ", @{$bvs}));
        $sql =~ s/\?/BIND_VARIABLE/g;
        my @binds = @{$bvs};
        while ($#binds > -1 && $sql =~ /BIND_VARIABLE/) {
            my $bv = shift @binds;
            if (!defined $bv) {
                $bv = "NULL";
            } elsif ($bv !~ /^\-?(\d+|\d*\.\d+)$/) {
                $bv =~ s/\\/\\\\/g;
                $bv =~ s/\'/\\\'/g;
                $bv = "'$bv'";
            }
            $sql =~ s/BIND_VARIABLE/$bv/;
        }
        if ($sql =~ /BIND_VARIABLE/) {
            $sql =~ s/BIND_VARIABLE/\?/g;
            $header .= sprintf($sqlcom, 'SQL ERROR', "Too few bind variables");
        } elsif ($#binds > -1 ) {
            $header .= sprintf($sqlcom, 'SQL ERROR', "Extra bind variables: ".
                               join(', ', @binds));
        }
    }
    $header .= sprintf($sqlcom, $args->{NAME} || "Un-named SQL", $callHist);
    my $text = $header . $sql . "\n";
    return $text;
}

sub explain {
    my $self = shift;
    my $driver = $self->driver();
    if (my $cb = $self->{private_EXPLAIN_METHOD} ||=
        $self->_eval_can('explain')) {
        return &{$cb}($self, @_);
    }
    $self->err("Module does not have SQL to explain queries for $driver");
    return [];
}

sub explain_text {
    my $self = shift;
    my $rows = $self->explain( @_ );
    return join("\n", map { join("|", @{$_}) } @{$rows}) || "";
}

sub show_and_explain {
    my $self = shift;
    return $self->pretty_print(@{$_[0] || []}).
        "[-- Query Plan --]\n".$self->explain_text(@_)."\n";
}

sub statistics {
    my $self = shift;
    my $args   = $self->parseparams( -degree => 8, @_ );
    my @params = ( "ownname => '" . $self->schema_owner() . "'" );
    my $meth   = 'gather_schema_stats';
    if (my $tab = $args->{TABLE}) {
        $meth = 'gather_table_stats';
        push @params, "tabname => '".uc($tab)."'";
    }
    if (my $deg = $args->{DEGREE}) {
        push @params, "degree => $deg";
    }
    if (my $perc = $args->{PERCENT} || $args->{ESTIMATE}) {
        push @params, "estimate_percent => $perc";
    }
    my $sql = sprintf
        ("BEGIN DBMS_STATS.%s( %s ); END;", $meth, join(', ', @params));
    my $sth = $self->prepare
        ( -sql   => $sql, 
          -name  => "Calculate statistics", 
          -level => 1 );
    $sth->execute(  );
}

sub _explain_postgres {
    my $self = shift;
    my $sql  = shift;
    my $bnd  = shift || [];
    my $ana  = shift() ? ' ANALYZE' : '';
    my $sth = $self->prepare
        ( -sql   => "EXPLAIN$ana $sql", 
          -name  => "Explain postgres query", 
          -level => 1 );
    $sth->execute( @{$bnd} );
    my $rows = $sth->fetchall_arrayref;
    return $rows;
}

# http://docs.oracle.com/cd/E11882_01/server.112/e16638/ex_plan.htm#i16938
sub _explain_oracle {
    my $self = shift;
    $self->bench_start();
    my $sql  = shift;
    my $bnd  = shift || [];
    my $ana  = shift() ? ' ANALYZE' : '';
    my $sid  = $$;
    my $clr = $self->prepare
        ( -sql   => "DELETE FROM plan_table WHERE  statement_id = '$sid'", 
          -name  => "Clear explain table",
          -level => 1 );
    $clr->execute();

    my $sth = $self->prepare
        ( -sql   => "EXPLAIN PLAN set statement_id = '$sid' FOR $sql", 
          -name  => "Explain oracle query", 
          -level => 1 );
    $sth->execute(  ); # Apparently does not want bind values provided

    my $get = $self->prepare
        ( -sql => "
 SELECT lpad('  ',level-1) || operation || ' ' || options || ' ' ||
        object_name \"Plan\"
   FROM plan_table
CONNECT BY prior id = parent_id AND prior statement_id = statement_id
  START WITH id = 0 AND statement_id = '$sid'
  ORDER BY id",
          -name => "Get details from PLAN_TABLE",
          -level => 2 );
    $get->execute();
    my $rows = $get->fetchall_arrayref;
    $self->bench_end();
    return $rows;
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::FriendlyDBI::st;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

use strict;
use vars qw(@ISA);
@ISA = qw(DBI::st);

sub name {
    my $self = shift;
    $self->{private_NAME} = $_[0] if (defined $_[0]);
    return $self->{private_NAME};
}

sub level {
    my $self = shift;
    $self->{private_LEVEL} = $_[0] if (defined $_[0]);
    return $self->{private_LEVEL};
}

sub ignore_errors {
    my $self = shift;
    if ($#_ > -1) {
        my @errs;
        if ($#_ == 0 && ref($_[0])) {
            # List passed as array ref
            @errs = @{$_[0]};
        } else {
            @errs = @_;
        }
        $self->{private_IGNORE_ERRORS} = \@errs;
    }
    return $self->{private_IGNORE_ERRORS};
}

sub explain {
    my $self = shift;
    return $self->{Database}->explain
        ( $self->{Statement}, @_ );
}

sub explain_text {
    my $self = shift;
    return $self->{Database}->explain_text
        ( $self->{Statement}, @_ );
}

sub pretty_print {
    my $self = shift;
    return $self->{Database}->pretty_print
        ( $self->{Statement}, 
          { NAME  => $self->name, 
            BIND  => ($#_ > -1) ? [ @_ ] : undef, 
            LEVEL => $self->level } );
}

sub show_and_explain {
    my $self = shift;
    return $self->pretty_print(@{$_[0] || []}).
        "[-- Query Plan --]\n".$self->explain_text(@_);
}

sub substitute_binds {
    my $self = shift;
    my $sql  = $self->{Statement} . " ";
    # The extra space is in case a ? is at the end of the SQL.
    my $dbh  = $self->{Database};
    my @bits = split(/\?/, $sql);
    my $newSql = shift @bits;
    if ($#bits != $#_) {
        $dbh->death("SQL statement has different number of placeholders than provided bind values:\n".$dbh->pretty_print($sql)."\n".join('+', map { !defined $_ ? '-NULL-' : "'$_'" } @_));
    }
    for my $b (0..$#bits) {
        my $v = $_[$b];
        if (defined $v) {
            $newSql .= $dbh->quote($v);
        } else {
            # This should be polymorphic for Oracle / Postgres
            # The line below will not work for PG
            # However, nulls need to be handled in the SQL construction
            # anyway, since foo = ? will never evaluate true for NULL
            # (need foo IS NULL instead)
            $newSql .= $dbh->quote("");
            # $newSql .= 'NULL';
        }
        $newSql .= $bits[$b];
    }
    $newSql =~ s/\s+$//;
    # warn $dbh->pretty_print($newSql);
    return $newSql;
}

sub flex_execute {
    my $self = shift;
    my ($binds, $notPrepared) = @_;
    my $using = $self;
    if (!$binds || $#{$binds} == -1) {
        # No bind values, always execute using existing STH
        $binds ||= [];
    } elsif ($notPrepared) {
        # Request to run the query using a statement handle without bind values
        $using = $self->{Database}->prepare
            ( -level  => $self->level(),
              -name   => $self->name(),
              -ignore => $self->ignore_errors(),
              -sql    => $self->substitute_binds( @{$binds} ) );
        $binds = [];
    }
    # warn "BINDS:\n".join(' + ', @{$binds})."\n";
    $using->execute(@{$binds});
    return $using;
}

sub execute {
    my $self = shift;
    my $dbh = $self->{Database};
    if (my $dbg = $dbh->{private_DEBUG}) {
        if ( $dbg == 9 || 
             ($self->{private_LEVEL} && $self->{private_LEVEL} <= $dbg)) {
            if ($dbh->{private_NESTED_BEGIN} && ! $dbh->{private_BEGIN_WARN}) {
                warn "\n/* BEGIN WORK */\n";
                $dbh->{private_BEGIN_WARN} = 1;
            }
            warn $self->pretty_print( @_ );
        }
    }
    my $rv;
    eval {
        $rv = $self->SUPER::execute( @_ );
    };
    if ($self->err) {
        my $throwit = 1;
        if (my $arr = $self->ignore_errors) {
            # The user wants to ignore some errors
            my $str = $self->errstr;
            foreach my $match (@{$arr}) {
                $throwit = 0 if ($str =~ /\Q$match\E/);
            }
            # warn "[!] Check '$str'" if ($throwit);
        } else {
            # warn "NO CHECK AVAILABLE";
        }
        if ($throwit) {
            $dbh->verbose_error( $self->error( @_ ) );
        } elsif (my $arr = $dbh->{private_BEGIN_BLOCK_ERRS}) {
            push @{$arr}, $self->error( @_ );
        }
    }
    return $rv;
}

*select_single_value = \&get_single_value;
sub get_single_value {
    my $self = shift;
    my ($rv) = $self->selectrow_array( @_ );
    return $rv;
}

sub selectrow_array {
    my $self = shift;
    $self->execute( @_ );
    my $rows = $self->fetchall_arrayref();
    return $rows->[0] ? @{$rows->[0]} : ();
}

sub selectcol_array {
    my $self = shift;
    $self->execute( @_ );
    my $rows = $self->fetchall_arrayref();
    return map { $_->[0] } @{$rows};
}

sub selectall_arrayref {
    my $self = shift;
    $self->execute( @_ );
    return $self->fetchall_arrayref();
}

sub get_array_for_field {
    my $self = shift;
    my $arr  = $self->selectall_arrayref( @_ );
    my @rv   = map { $_->[0] } @{$arr};
    return @rv;
}

sub error {
    my $self = shift;
    my $txt  = $self->pretty_print( @_ ) || "";
    my $estr = $self->errstr;
    $txt .= sprintf("\nSTH Failure: ERR %d : %s\n\n%s\n  ", $self->err,
                    $estr, &BMS::FriendlyDBI::db::stack_trace());
    if ($estr =~ /current transaction is aborted/) {
        if (my $arr = $self->{Database}{private_BEGIN_BLOCK_ERRS}) {
            $txt .= "\nErrors occuring within the current transaction block:";
            for my $i (0..$#{$arr}) {
                $txt .= "\n[BLOCK ERROR $i]\n$arr->[$i]";
            }
        } else {
            $txt .= "\n!! No Transaction block error messages found!\n";
        }
    }
    return $txt;
}

return 1;
