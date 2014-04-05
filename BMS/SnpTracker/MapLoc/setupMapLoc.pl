#!/usr/bin/perl - w

# You may need to modify the above line to point to your version of
# Perl, if it does not reside at the location shown above.

use strict;

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

my $shellPreCol = "\033["; #1;33m
my $shellProCol = "\033[0;49m";
my $shellTermCols = {
    'x'           => ['37',   'ClearMsg'],
    '*'           => ['32',   'GenericMsg'],
    '#'           => ['36',   ''],
    '^'           => ['34',   ''],
    '!'           => ['31;7', 'StrongErrorMsg'],
    '!!'          => ['31;43', 'StrongErrorMsg'],
    '?'           => ['33;7', 'WarnMsg'],
    '+'           => ['37',   ''],
    '-'           => ['33',   ''],
    'FILE'        => ['37', ''],
    ''            => ['', ''],
};

my (@fixes, %stuff);
my $status = "";

my $xtras = {
    # http://stackoverflow.com/questions/6040583/cant-find-the-libpq-fe-h-header-when-trying-to-install-pg-gem-on-ubuntu
    'DBD::Pg' => [ \&dbd_pg_details, 'sudo apt-get install libpq-dev' ],
};

my $baseDir = $0;
$baseDir    =~ s/\/[^\/]+$//;
$baseDir    = "." unless ($baseDir);

&welcome();
&check_exec();
&check_modules();
&check_postgres();

&finish();

sub halt_if_needed {
    return if ($#fixes == -1);
    my $txt = <<INFO;

# At least some dependencies were not found. The commands below should
# install the required components. Please note that 'apt-get' is used
# by Debian-based systems; If you are using a different distribution
# you should substitute the relevant package management code.

INFO

    &markup_notes($txt );
    foreach my $fix (@fixes) {
        if (my $r = ref($fix)) {
            &{$fix}();
        } else {
            warn "$fix\n";
        }
    }
    &try_again();
}

sub check_modules {
    my $file = &progress_file("RequiredPerlModules");
    if (-s $file) {
        &msg("[FILE]", "All required Perl modules are available:", $file);
        return;
    }
    &msg("[^]","Checking Perl Modules");
    my @basicML = map { "BMS::SnpTracker::MapLoc::$_" }
        qw(Accessioned Alignment Common Featured GeneAlignment
           Gene HasAlignments HSPed IsolatedHSP Location LockedLocation
           Locked PopulationFilter Population Range Reporter RNA
           RnaToGeneAlignment SnpProjector TaggedHSP Tagged
           TaggedRange Text TransientAlignment TransientLocation Transient
           TransientText );
    map { &find_module($_) }
    ('CGI', 'DBD::Pg', 'DBI',
     ['Time::HiRes',1], 'LWP::UserAgent', 'Scalar::Util', 'Math::Trig',
     'Spreadsheet::ParseExcel', 'Spreadsheet::XLSX', 
     'Spreadsheet::WriteExcel', 'Excel::Writer::XLSX',

     'Bio::SeqIO', 'Bio::PrimarySeq',

     'XML::Parser::PerlSAX',

     'BMS::Utilities', 'BMS::Utilities::Debug', 'BMS::Utilities::Serialize',
     'BMS::Utilities::FileUtilities', 'BMS::ExcelHelper'

     'BMS::SnpTracker::MapLoc', @basicML );


    &halt_if_needed();
    &write_progress($file);
    
}

sub check_exec {
    my $file = &progress_file("RequiredExecutables");
    if (-s $file) {
        &msg("[FILE]", "All required executables are available:", $file);
        return;
    }
    
    $status = "";
    &msg("[^]","Checking executables");
    
    &find_exec("psql", 0, "Postgres is the database engine that underlies MapLoc. Installation will be impossible without it.", 'postgresql');
    
    &find_exec("cpan", 1, "cpan is a convienent mechansim to install Perl .");
    
    &halt_if_needed();
    &write_progress($file);
}

sub progress_file {
    my $base = shift;
    my $file = $0;
    $file    =~ s/\/[^\/]+$//;
    $file   .= "/$base.txt";
    return $file;
}

sub write_progress {
    my $file = shift;
    open(PFILE, ">$file") || die "Failed to write progress file\n  $file\n  $!\n  ";
    print PFILE $status;
    close PFILE;
}

sub msg {
    my @msg = @_;
    return if ($#msg == -1);
    my $firstTok = '[!]';
    my $tca;
    if ($msg[0] && $msg[0] =~ /^\s*\[([^\]]+)\]$/) {
        $firstTok = shift @msg;
        my $tok = $1;
        $tca = $shellTermCols->{$tok} ? $shellTermCols->{$tok}[0] : undef;
    }
    my $pad = " " x length($firstTok);
    my @bits;
    for my $i (0..$#msg) {
        push @bits, sprintf("%s %s", $i ? $pad : $firstTok, 
                            defined $msg[$i] ? $msg[$i] : "-UNDEF-");
    }
    my $msg  = join("\n", @bits);
    $status .= $msg . "\n";
    $msg     = sprintf("%s%sm%s%s", $shellPreCol, $tca,
                       $msg, $shellProCol) if ($tca);
    $msg    .= "\n";
    warn $msg;
}


sub find_exec {
    my ($cmd, $isOpt, $msg, $pkg) = @_;
    exit &msg("[!!]", "'$cmd' does not look kosher") 
        unless ($cmd =~ /^[a-z]+$/i);
    my $check;
    if ( $check = `which "$cmd"`) {
        $check =~ s/[\n\r]+$//;
        &msg("[+]","$cmd : $check");
    } else {
        &msg("[!]", "Failed to locate executable for '$cmd'", $msg);
        if ($pkg) {
            push @fixes, "sudo apt-get install $pkg";
        } else {
            &msg("[?]", "No package provided, uncertain of remediation steps");
        }
    }
    return $check;
}

sub find_module_cpan {
    my ($mod, $isOpt) = @_;
    my @lines = split(/\n/, `cpan -D $mod`);
    my $tok;
    my @msg = ($mod);
    my ($iv, $cv);
    my $needFix;
    foreach my $line (@lines) {
        if ($line =~ /^\s+CPAN\:\s*(.+?) ((not )?up to date)/i) {
            $cv = $1;
            $cv =~ s/^\s+//; $cv =~ s/\s+$//;
            if ($3) {
                $tok = "[!]";
                $needFix = 1;
            } else {
                $tok ||= "[+]";
            }
        } elsif ($line =~ /^\s+Installed\:\s*$/) {
            $tok = "[!]";
            $needFix = 1;
        } elsif ($line =~ /^\s+Installed\:\s*(\S.*\S)/) {
            $iv = $1;
        } elsif ($line =~ /\.pm$/) {
            $line =~ s/^\s+//;
            push @msg, $line;
        }
    }
    if ($needFix) {
        if (my $xtra = $xtras->{$mod}) {
            push @fixes, @{$xtra};
        }
        push @fixes, "cpan -i $mod";
    }
    $tok ||= "[!!]";
    $iv ||= "";
    push @msg, 
    !$cv ? "Module not located in CPAN" : 
        !$iv ? "Uninstalled (CPAN $cv)" : 
        $iv ne $cv ? "Installed: $iv (CPAN $cv)" :
        "Installed: $iv (up to date)";
    &msg($tok, @msg);
}

sub find_module {

    # This does not work for Spreadsheet::ParseExcel
    # If I 'use' the module in the code body, all is well.
    # If I 'require' the module below, then the require fails with:
    # Modification of a read-only value attempted at /usr/local/share/perl/5.14.2/Digest/Perl/MD5.pm line
    # Going to use cpan to check for non-BMS module availability

    my $mod = shift;
    my $isOpt;
    if (ref($mod)) {
        ($mod, $isOpt) = @{$mod};
    }
    return &find_module_cpan( $mod, $isOpt ) unless ($mod =~ /^BMS/);

    my $path = $mod;
    $path .= '.pm' unless ($path =~ /\.p[ml]$/);
    $path =~ s/\:\:/\//g;
    eval {
        require $path;
        1;
    };
    if (my $err = $@) {
        if ($err =~ /Can\'t locate/) {
            unless ($stuff{INC}++) {
                &msg("[#]","Perl \@INC = ", @INC);
            }
            $err = "Not found in \@INC";
        }
        my $tok = "[!]";
        if ($isOpt) {
            $err = "<OPTIONAL> $err";
            $tok = "[?]";
        } elsif ($mod =~ /^BMS/) {
            $tok = "[!!]";
        } else {
            push @fixes, "cpan -i $mod";
        }
        &msg($tok, "Failed to locate BMS module '$mod'", $err);
    } else {
        my $found = $INC{$path};
        &msg("[+]", $mod, $found);
    }
}


sub try_again {
    &msg("[?]", "Please follow instructions above, then re-run this script");
    exit;
}

sub markup_notes {
    foreach my $txt (@_) {
        foreach my $line (split(/\n/, $txt)) {
            my $tca;
            if ($line =~ /^\#/) {
                $tca = 31;
            } elsif ($line =~ /^http/) {
                $tca = 34;
            }
            my $msg = $line;
            if ($tca) {
                $msg = sprintf("%s%sm%s%s", $shellPreCol, $tca,
                               $msg, $shellProCol);
            }
            warn "$msg\n";
        }
    }
}


sub dbd_pg_details {    
    
    my $txt = <<INFO;

# To install the DBD::pg driver, you will be asked to specify the
# PG version, which can be found with:

psql --version

# You will also be asked for several "valid directories". To be honest
# I am not sure what these directories 'should' be; I just created a
# directory and made it world-accessible. This is maybe a bad idea,
# but it appears to allow the installation:

sudo mkdir /usr/psql
sudo chmod 0777 /usr/psql

# I used the above directory in every query requesting a "valid" directory
# I also had to install libpq-dev to get DBD::pg to pass tests:

http://stackoverflow.com/questions/6040583/cant-find-the-libpq-fe-h-header-when-trying-to-install-pg-gem-on-ubuntu
 
INFO

    &markup_notes($txt );
}

sub check_postgres {
    my $pgList = `psql -d postgres -c '\\dS'`;
    my @lines  = split(/\n/, $pgList);
    if ($#lines < 5) {
        &markup_notes( &postgres_setup() );
        &try_again();
    }
    &msg("[+]","Postgres online and ready");

    my $args;
    eval {
        require BMS::ArgumentParser;
        $args = BMS::ArgumentParser->new();
    };
    unless ($args) {
        &msg("[!!]", "Failed to instantiate ArgumentParser object",
             $@,
             "Uncertain what to do at this point");
        exit;
    }

    my $instance = $args->val(qw(instance)) || 'maploc';

    my $isThere = "psql -d $instance -c '\\l'";
    my $chk = `$isThere`;
    if (!$chk || $chk =~ /FATAL/ || $chk !~ /List of databases/) {
        system("createdb maploc");
        $chk = `$isThere`;
        if (!$chk || $chk =~ /FATAL/ || $chk !~ /List of databases/) {
            &msg("[!!]", "Failed to auto-create the maploc database");
            &database_creation();
            &try_again();
        }
        &msg("[+]","maploc database '$instance' succesfully created");
    } else {
        &msg("[+]","maploc database '$instance' already exists");
    }

    if (my $incs = $args->val(qw(inc))) {
        push @INC, ref($incs) ? @{$incs} : $incs;
    }
    
    my $ml;
    eval {
        require BMS::SnpTracker::MapLoc;
        $ml = BMS::SnpTracker::MapLoc->new();
    };
    unless ($ml) {
        &msg("[!!]", "Failed to instantiate MapLoc object",
             $@,
             "Uncertain what to do at this point");
        exit;
    }
    my $dbh = $ml->dbh;
    &msg("[-]", "Creating / updating database tables");
    $dbh->make_all();
    &msg("[+]", "Database tables created");
    $dbh->disconnect();
}

sub welcome {
    my $txt = <<INFO;

# This script will attempt to configure the MapLoc variation warehouse
# on this computer. It will verify that the required executables and
# Perl libraries are installed, and will make sure that the
# appropriate schema exists for holding data.
 
INFO

    &markup_notes( $txt );
}

sub postgres_setup {
    return <<INFO;

# Generally installation of postgres should set up the underlying
# system databses and set aside a location for the data files to
# reside in. You should see at least one cluster by using:

pg_lsclusters

# If you do not, then you will need to make one manually. The steps
# below should create one. I think.

# First you need a loacation which will hold the data files. This can
# be anywhere you choose, but it should be owned by postgres:

sudo mkdir -p /usr/local/pgsql/data
sudo chown -R postgres /usr/local/pgsql/
sudo chgrp -R postgres /usr/local/pgsql/

# Become postgres for the following operations:

sudo su - postgres

# You want to use initdb, but it is not always obvious where it
# resides; use locate to find it:

locate initdb

# On my install, it was in /usr/lib/postgresql/9.1/bin/, so:

/usr/lib/postgresql/9.1/bin/initdb -D /usr/local/pgsql/data

# This should create a cluster stored within the directory that you
# specified.

# Alternatively, if you know that the database is configured, you may
# be seeing this message simply because the server is currently not
# running.

INFO

}

sub database_creation {

    return <<INFO;

# This installer will try to create the maploc database for you. If it
# fails, you will need to do so manually.  To add youself as a user,
# use CREATE ROLE. The command below allows user tilfordc to create
# databases and login to the system

sudo su - postgres

psql -c 'CREATE ROLE tilfordc WITH CREATEDB LOGIN'

exit

# You may wish to try re-running the installer at this point, as not
# being registered as a role in postgres with create priveleges is
# sufficient to prevent database creation. Alteratively, you can just
# create the database now:

createdb maploc
 
INFO


}

sub finish {
    my $txt = <<INFO;

# The database should now be created and ready for data. You can run
# the following scripts to begin populating it with information. The
# scripts can also serve as examples of how to load your own data into
# the system.


# Load the Ensembl transcriptome into the system:
$baseDir/load_ensembl.pl

# Please ignore any "uncleared implementors" errors. This appears to
# be an issue with some versions of Perl? It does not appear to be a
# problem, and I have only observed it when a database handle is
# freed.

# --NOTE--

# When you first begin loading data into the database, Postgres will
# not have statistics on the tables being used. This will cause
# performance to gradually degrade as the load progresses. In Postgres
# statistics are generated with the ANALYZE command; it is generally a
# good idea to run an analyze after the row composition has changed
# significantly. During the first loads for both RNAs and variants, I
# recommend after the first 30-60 minutes of loading to run a full
# analysis of the DB:

psql -c 'VACUUM ANALYZE VERBOSE' maploc

# Also, run the same command once a large load has finished. The above
# command is structured to run from a system shell, but you could just
# run ANALYZE VERBOSE from within a psql session. VACUUM is not
# neccesary for updating indices, but is not a bad idea. VERBOSE can
# be left out as well. If your database is being frequently altered,
# you can cron the command to keep index statistics fresh.

INFO

    &markup_notes( $txt );
    exit();
}

