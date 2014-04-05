# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::ExcelHelper;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

=head1 DESCRIPTION

  BMS::ExcelHelper - Easily generate Excel workbooks

This module is a wrapper around Excel::Writer::XLSX. It is
designed to be easier to use by automatically performing workbook
configuration tasks. It also allows creation of XLS documents

Note: Even though WriteExcel references multiple child classes (Format,
Formula, Workbook, Worksheet), all documentation for those classes are
centralized in Excel::Writer::XLSX.

=head1 SYNOPSIS

    use strict;
    use BMS::ExcelHelper;

    my $pathBase = "/tmp";
    my $urlBase  = "http://www.example.com/tmp";
    my $fileName = "MasterGuide.xlsx";

    my $eh  = BMS::ExcelHelper->new( "$pathBase/$fileName" );
    # Setting the URL allows html_summary() to link to the file:
    $eh->url("$urlBase/$fileName");

    $eh->sheet( -name    => "Vegetables",
                -freeze  => 1,
                -columns => [ 'Name', 'Color', 'Triva' ],
                -width   => [ 12, 20, 80 ], );

    $eh->sheet( -name    => "Fruits",
                -columns => [ 'Name', 'Color' ], );

    # Red background, bold yellow text
    $eh->format( -name     => 'alert',
                 -bold     => 1,
                 -color    => 'yellow',
                 -bg_color => 'red',
                 -pattern  => 1 );

    # Add some veggies
    $eh->add_row('Vegetables', ['Spinach', 'Green', 'Not so tasty']);
    $eh->add_row('Vegetables', ['Red Pepper', 'Red', 'Quite tasty']);
    # Let's also format this row with the 'alert' style:
    $eh->add_row('Vegetables', ['Jalepeno', 'Green', 'Hot!'],
                 [undef,undef,'alert']);

    my %fruits = ( Blueberry => 'Blue',
                   Apple     => 'Red',
                   Grape     => 'Purple' );

    while (my ($fruit, $color) = each %fruits) {
        $eh->add_row('Fruits', [ $fruit, $color ] );
    }

    $eh->close;
    if ($ENV{'HTTP_HOST'}) {
        print $eh->html_summary;
    } else {
        warn "Workbook saved to " . $eh->file_path( );
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

BEGIN {
    use lib '/stf/biocgi/tilfordc/patch_lib';
}

use vars qw(@ISA);
use strict;
use BMS::ErrorInterceptor;

@ISA = qw(BMS::ErrorInterceptor);


my $VERSION = 
    ' $Id: ExcelHelper.pm,v 1.35 2014/03/20 18:18:12 tilfordc Exp $ ';

=head2 new

 Title   : new
 Usage   : my $eh = BMS::ExcelHelper->new( $filePath );
 Function: Creates a new Perl object
 Returns : A blessed BMS::ExcelHelper object
    Args : [0] The path to the file (required)

This function will also call set_basic_formats(). ExcelHelper objects
are subclasses of Excel::Writer::XLSX::Workbook (but all
documentation for that module is in Excel::Writer::XLSX).

If a single value is passed, it is presumed to be the file
path. Otherwise, it is treated as param/value hash that recognizes the
following keys:

      -file The path for the newly created Excel file

 -forcexlsx If true, then files with 'xls' suffices will be changed to
            XLSX files. Otherwise, if the provide file path has an xls
            suffix, then Spreadsheet::WriteExcel will be used to make
            an XLS document.

If an XLS document is created, the program will warn the user when
close() is called if more than 65536 rows were added to any given
sheet.

=cut


our $globals = {};
sub new {
    my ($class, @arglist) = @_;


    
    unshift @arglist, '-file' if ($#arglist == 0);
    my $otherSelf = {};
    bless ($otherSelf, $class);
    my $args = $otherSelf->parseparams( @arglist );
    $otherSelf->intercept_errors();
    my $file = $args->{FILE};
    foreach my $safeToIgnore ("ConfigLocal",
                           "wrapped in pack",
                           "Can't locate Encode/ConfigLocal.pm",
                           "Can't locate Digest/MD4.pm",
                           "Can't locate Digest/Perl/MD4.pm",
                           "Malformed UTF-8 character") {
        $otherSelf->ignore_error( $safeToIgnore);
    }

    my $self;
    if ($file =~ /\.(xls)$/i) {
        if ($args->{FORCEXLSX}) {
            $file .= ($1 eq 'XLS') ? 'X' : 'x';
        } else {
            $self = BMS::ExcelHelper::XLS->new( $file );
            $self->{isXLS} = 1;
        }
    }
    $self ||= BMS::ExcelHelper::XLSX->new( $file );

    unless ($self) {
        $otherSelf->err("Failed to create excel spreadsheet object",
                        $file, $!);
        return $self;
    }
    # warn $file;
    # This tag is a parasite that will hold the meta-data I want to
    # associate with the workbook:
    $self->{_CAT_} = {
        SHEETS => {},
        WBFORM => {},
        PFORM  => {
            warn => $ENV{'HTTP_HOST'} ? '<pre><font color=red>%s</pre>':'%s\n',
        },
        PARAMS => {
            allow_dup => 0,
        },
        PATHS => {
            file => $file,
        },
        URL => "",
        DEF_WID => {},
    };
    $self->set_basic_formats;
    return $self;
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::ExcelHelper::XLS;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

=head1 DESCRIPTION

This package is a subclass of BMS::ExcelHelper::Common used when
generating XLS workbooks. It includes safety checking that warn when
too many rows (over 65536) are added, or too many hyper links (65530).

=cut

use vars qw(@ISA);
use strict;
use Spreadsheet::WriteExcel;
# use BMS::ExcelHelper::Common;
@ISA = qw(Spreadsheet::WriteExcel BMS::ExcelHelper::Common);

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my ($file) = @_;
    return $class->SUPER::new( $file );
}

sub close {
    my $self = shift;
    return unless ($self);
    my $count = $self->row_count;
    my @probs;
    my $maxRows = 65536;
    foreach my $sheet (sort keys %{$count}) {
        my $num = $count->{$sheet} - 1;
        push @probs, "$sheet : $num" if ($num > $maxRows);
    }
    $self->msg("[DATA LOSS]","Excel 97 can only support $maxRows rows per sheet",
               "Some of your sheets exceed this limit", @probs)
        unless ($#probs == -1);
    if (my $uw = $globals->{URLS_WRITTEN}) {
        my $maxUrl = 65530;
        my $over   = $uw - $maxUrl;
        $self->msg("[CORRUPT]","Your workbook will likely have issues loading",
                   "Excel supports at most $maxUrl hyperlinks, you have $uw (+$over)") if ($over > 0);
    }
    
    return $self->SUPER::close();
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::ExcelHelper::XLSX;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

=head1 DESCRIPTION

This package is a subclass of BMS::ExcelHelper::Common used when
generating XLSX workbooks.

=cut

use vars qw(@ISA);
use strict;
use Excel::Writer::XLSX;
# use BMS::ExcelHelper::Common;
@ISA = qw(Excel::Writer::XLSX BMS::ExcelHelper::Common);

sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my ($file) = @_;
    return $class->SUPER::new( $file );
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
package BMS::ExcelHelper::Common;
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

=head1 DESCRIPTION

This package provides the primary functionality, and is used by both BMS::ExcelHelper::XLSX and BMS::ExcelHelper::XLS .

=cut

use vars qw(@ISA);
use strict;
use BMS::ErrorInterceptor;
@ISA = qw(BMS::ErrorInterceptor);

sub DESTROY {
    my $self = shift;
    $self->close if ($self);
}

=head2 close

 Title   : close 
 Usage   : $eh->close( );
 Function: Finalize creation of the workbook
 Returns : Nothing
    Args : None

It is neccesary to call this function in order to finish writing the
workbook to file! If the object goes out of scope and is garbage
collected, a DESTROY( ) method is supposed to call close( ), but you
should not rely on that.

=cut

sub close {
    my $self = shift;
    # print "<pre>CLOSING [$! / $@]\n".$self->stack_trace()."</pre>";
    unless ($self->{_CAT}{ALREADY_CLOSED}++) {
        if ($self->{_worksheets}) {
            return $self->SUPER::close;
        } else {
            # $self->msg("[!!]", "Excel workbook is empty");
        }
    }
    return undef;
}

=head2 allow_duplicate_rows

 Title   : allow_duplicate_rows
 Usage   : my $val = $eh->allow_duplicate_rows( $newVal );
 Function: Gets / sets flag indicating if duplicate rows are allowed in a sheet
 Returns : The current value. Default is false
    Args : [0] Optional new value

If this value is false (zero or empty string), then a worksheet will
only contain one instance of any given row (defined by the
tab-concatenated values of that row). Additional attempts to add the
same row to a given worksheet using add_row() or add_row_explicit()
will be ignored.

If you wish to be able to add duplicate rows, call the method with a
true value (might we recommend 1?).

=cut

sub allow_duplicate_rows {
    my $self = shift;
    if (defined $_[0]) {
        $self->{_CAT_}{PARAMS}{allow_dup} = $_[0];
    }
    return $self->{_CAT_}{PARAMS}{allow_dup};
}

=head2 default_width

 Title   : default_width
 Usage   : my $width = $eh->default_width( $columnName, $newWidth );
 Function: Gets / Sets the default width for a column name
 Returns : The width
    Args : [0] The name of the column (case insensitive)
           [1] Optional new width value

Allows columns to have a width set based on their name. This is useful
if you have many worksheets with the same column name repeated.

If a value is not set, the default width will be twice the number of
characters in a column name.

This function is called by sheet() when setting up a new worksheet.

=cut

sub default_width {
    my $self = shift;
    my ($col, $val) = @_;
    return undef unless ($col);
    if ($val) {
        $self->{_CAT_}{DEF_WID}{uc($col)} = $val;
    }
    $val ||= $self->{_CAT_}{DEF_WID}{uc($col)} || length($col) * 2;
    return $val;
}

=head2 default_format

 Title   : default_format
 Usage   : my $format = $eh->default_format( $columnName, $newFormat );
 Function: Sets a default format for a column name
 Returns : A Excel::Writer::XLSX::Format object
    Args : [0] The name of the column (case insensitive)
           [1] Optional new format (name or object)

Allows columns to have a format set based on their name. This is useful
if you have many worksheets with the same column name repeated.

Before calling this method the relevant format should be defined using
format(). Formats can be provided as names or format objects.

=cut

sub default_format {
    my $self = shift;
    my ($col, $val) = @_;
    return undef unless (defined $col);
    if ($val) {
        $self->{_CAT_}{DEF_FRM}{uc($col)} = $val;
    }
    my $req = $self->{_CAT_}{DEF_FRM}{uc($col)};
    return $req ? $self->format($req) : undef;
}

=head2 sheet

 Title   : sheet 
 Usage   : my $ws = $eh->sheet( @params );
           my $ws = $eh->sheet( $sheetName );
 Function: Defines a new sheet, or retrieves an existing one
 Returns : A Excel::Writer::XLSX::Worksheet object
    Args : Associative array of arguments.

Recognized keys [Default]:

     -name Required. The name of the worksheet; what is displayed in the
           sheet tab

     -cols An array reference of column names, in the order they
           should appear. Can also use the argument -columns. Will
           cause the array to be added as a row via add_row(). 

 -colformat An optional array of formats to be applied to the column
            heading.

    -width Optional array reference of numbers, defining the column
           widths. 

   -format Optional array reference of formats. These will be used as
           the default format for any cell in the column (other than
           the header itself)

   -freeze Optional, if true then the first row will be frozen (ie
           will not scroll with the rest of the document; useful for
           locking a header row at the top of the page)

If -width or -format are defined, Excel::Writer::XLSX::Worksheet::set_column()
will be used to apply them to the entire column.

The primary use of this function is to establish a new
worksheet. However, it can also be used to retrieve a worksheet by its
name. If the method is called with a single argument - $eh->sheet(
'foobar' ) - then the argument is assumed to be a worksheet name. The
case of the name is preserved for display purposes, but is
case-insensitive for creating and retrieving worksheets.

If a worksheet does not exist with the provided name, it will be
created automatically. See has_sheet() for a method that tests for
sheet existance without triggering creation.

Note that there are disallowed characters for worksheet names. Also,
at most 31 characters will be used for the sheet name.

=cut

sub sheet {
    my $self = shift;
    my ($name, $args, $ws);
    if ($#_ == 0) {
        # Single argument assumed to be the name
        $name = $_[0];
    } else {
        $args = $self->parseparams
            ( -name   => undef,
              -cols   => undef,
              -freeze => 0,
              @_ );
        $name = $args->{NAME};
    }
    unless ($name) {
        $self->err("Request for sheet() with no name");
        return undef;
    }
    if (my $r = ref($name)) {
        if ($r eq 'Excel::Writer::XLSX::Worksheet' ||
            $r eq 'Spreadsheet::WriteExcel::Worksheet') {
            # Already a sheet object
            $ws = $name;
        } else {
            $self->death
                ("Can not interpret unknown object '$name' as worksheet");
        }
    } else {
        my $trim = $name;
        $trim =~ s/[\[\]\:\*\?\\\/]+/_/g;
        $trim = substr($trim, 0, 31);
        my $ucname = uc($trim);
        $ws = $self->{_CAT_}{SHEETS}{$ucname};
        unless ($ws) {
            # Need to make a new sheet
            $self->{_CAT_}{SHEETS}{$ucname} = $ws = $self->add_worksheet($trim);
            $ws->{_CAT_} = {
                DONE        => {},
                CURRENT_ROW => 0,
                NAME        => $name,
                TRIM        => $trim,
            };
        }
    }
    return $ws unless ($args);

    if (my $freeze = $args->{FREEZE}) {
        my ($x, $y) = ref($freeze) ? @{$freeze} : ($freeze);
        $ws->freeze_panes($x || 0, $y || 0);
    }

    my $cols    = $args->{COLUMNS} || $args->{COLS};
    my $widths  = $args->{WIDTH}   || $args->{WIDTHS};
    my $formats = $args->{FORMAT}  || $args->{FORMATS};
    if ($cols) {
        my $frmt = $args->{COLFORMAT};
        $self->add_row($ws, $cols, $frmt || 'header');
        $widths  ||= [ map { $self->default_width($_)  } @{$cols} ];
        $formats ||= [ map { $self->default_format($_) } @{$cols} ];
    }
    if ($widths || $formats) {
        # Proxy is used to cycle through both widths and formats
        my $proxy = $widths || $formats;
        $proxy = $formats 
            if ($widths && $formats && $#{$widths} < $#{$formats});
        for my $i (0..$#{$proxy}) {
            # warn "$name $i : ".($formats->[$i] || 'n/a')."\n";
            my $format = ($formats && $formats->[$i]) ? 
                $self->format($formats->[$i]) : undef;
            $ws->set_column($i, $i, 
                            $widths ? $widths->[$i] : undef, $format);
        }
    }

    return $ws;
}

=head2 has_sheet

 Title   : has_sheet
 Usage   : my $ws = $eh->has_sheet( $sheetName );
 Function: Checks to see if the workbook has a sheet of that name
 Returns : The worksheet, if it exists, or undef if not
    Args : [0] The worksheet name

Note that sheet() will create a worksheet if it does not exist; this
function allows you to test for sheet existance without triggering
creation of the sheet.

=cut

sub has_sheet {
    my $self = shift;
    my ($name) = @_;
    return $self->{_CAT_}{SHEETS}{uc($name)};
}

=head2 row_count

 Title   : row_count
 Usage   : my $num = $eh->row_count( $sheetName);
 Function: Gets the number of rows for a worksheet, or for entire workbook
 Returns : An integer
    Args : [0] Optional sheet name

If a sheet name is provided, sheet() is used to get (or create) the
worksheet, and the number of rows on that sheet will be calculated. If
no arguments are provided, then the total count of all rows in the
workbook will be returned.

=cut

sub row_count {
    my $self = shift;
    my ($sheetname) = @_;
    my $ws = $sheetname ? $self->sheet($sheetname) : undef;
    return $ws->{_CAT_}{CURRENT_ROW} if ($ws);
    my $count = {};
    while (my ($key, $ws) = each %{$self->{_CAT_}{SHEETS}}) {
        my ($row, $name) = ($ws->{_CAT_}{CURRENT_ROW}, $ws->{_CAT_}{NAME});
        $count->{$name} = $row + 1;
    }
    return $count;
}

*add_blank_row = \&blank_row;
sub blank_row {
    my $self = shift;
    my ($sheetname) = @_;
    my $ws = $sheetname ? $self->sheet($sheetname) : undef;
    return 0 unless $ws;
    return ++$ws->{_CAT_}{CURRENT_ROW};
}

=head2 format

 Title   : format
 Usage   : $eh->format( @params );
 Function: Get / create a format
 Returns : A Excel::Writer::XLSX::Format object
    Args : Associative array of arguments.

This function allows organization of formats by name. If you just want
to recover a format, you can pass a single argument representing the
name you chose for it:

  my $bb = $eh->format('bold blue').

If you want to define the properties of a format, you should pass an
associative array defining them. The array should always include the
name you are associating with the format. For example:

  $eh->format( -name     => 'emphasis',
               -bold     => 1,
               -bg_color => 'silver', 
               -fg_color => 'silver', 
               -pattern  => 1,                
               -color    => 'red' );

See the 'cell formatting' section of Excel::Writer::XLSX for more
information on format parameters. Coloring cells can be a little
confusing - in particular, for background coloring you need to
remember to set 'pattern' in addition to the color. As a convienence,
if you pass the parameter -background (eg -background => 'blue'), the
prorgram will automatically set -bg_color and -fg_color to that value,
and set -pattern => 1.

=cut

sub format {
    my $self = shift;
    my ($name, $args, $format);
    if ($#_ == 0) {
        # Single argument assumed to be the name
        $name = $_[0];
    } else {
        $args = $self->parseparams( -name => undef,
                                   @_ );
        $name = $args->{NAME};
    }
   unless (defined $name) {
        $self->err("Request for format() with no name");
        return undef;
    }
    if (ref($name)) {
        # Assume this is a format object
        $format = $name;
    } else {
        # Lookup by name
        my $ucname = uc($name);
        $format = $self->{_CAT_}{WBFORM}{$ucname} ||= $self->add_format();
    }
    return $format unless ($args);

    my %prop = %{$args};
    map { delete $prop{$_} } qw(NAME NOTE COMMENT);
    if (my $bg = $prop{BACKGROUND}) {
        $prop{BG_COLOR} = $prop{FG_COLOR} = $bg;
        $prop{PATTERN}  = 1;
        delete $prop{BACKGROUND};
    }
    my @pargs = map { lc($_) } ( %prop );
    $format->set_properties( @pargs ) if ($#pargs > -1);
    return $format;
}

=head2 set_basic_formats

 Title   : set_basic_formats
 Usage   : $eh->set_basic_formats( );
 Function: Sets a couple 'standard' formats
 Returns : Nothing
    Args : None

This function is called by new(), and sets two formats:

    'header' = centered bold
    'center' = centered

These are used, not surprisingly, to format headers specified by -cols
when the sheet() method is called.

=cut

sub set_basic_formats {
    my $self = shift;
    $self->format( -name  => 'header',
                   -bold  => 1,
                   -align => 'center' );
    $self->format( -name  => 'center',
                   -align => 'center' );
    $self->format( -name  => 'hyperlink',
                   -color => 'blue',
                   -underline => 1, );
}

=head2 add_row

 Title   : add_row
 Usage   : $eh->add_row( $sheetname, $arrayRef, $optionalFormat );
 Function: Adds a single row to the spreadsheet
 Returns : The current row count of the sheet, or zero if no row added
    Args : [0] The worksheet to add the row to
           [1] An array reference of cell contents
           [2] Optional format specification

If the sheet does not already exist, it will be created with
sheet(). Your array reference of formats is typically just going to be
a bunch of strings. The values will be added with write_row( ) from
Excel::Writer::XLSX. Note that if a 2D array is passed, write_row( )
assumes that each entry is a column - this may cause some strange
behavior. See add_row_explicit() if you have complex data you want to
add to the worksheet.

If a format is specified, the cells in the row will be formatted
appropriately. If the format is an array reference, it is assumed to
be matched to the row, such that you can format each cell in the row
independently (use undef for cells that should have no
format). Otherwise, a single format will be applied to every
cell. Formats can be provided by name or as Format objects (see the
format() function).

=cut

#my $cLim = 255;
#my $exNL = "\"\n\""; # 'CHAR(10)';
sub add_row {
    my $self = shift;
    my ($sheetname, $array, $format) = @_;
    my $ws = $self->sheet($sheetname);
    unless ($ws) {
        $self->err("No worksheet with name '$sheetname' in add_row()");
        return 0;
    }

    my $key = join("\t", map { defined $_ ? $_ : '' } @{$array});
    return 0 if ($ws->{_CAT_}{DONE}{$key} && !$self->allow_duplicate_rows);
    $ws->{_CAT_}{DONE}{$key}++;
    my $crow = $ws->{_CAT_}{CURRENT_ROW};
    if ($format && ref($format) eq 'ARRAY') {
        # Individual formats per cell
        my @formats = $self->_map_formats(@{$format});
        for my $col (0..$#{$array}) {
            my $val = $array->[$col];
            $ws->write($crow, $col, $val, $formats[$col]);
        }
    } else {
        ($format) = $self->_map_formats($format);
        $ws->write_row( $crow, 0, $array, $format);
    }
    # print "<pre>$ws->{_CAT_}{NAME} $crow = $key</pre>";
    return ++$ws->{_CAT_}{CURRENT_ROW};
}

sub apply_format {
    my $self = shift;
    my ($sheetname, $format, $crow, $col, $val) = @_;
    my $ws = $self->sheet($sheetname);
    unless ($ws) {
        $self->warn("No worksheet with name '$sheetname' in add_row()");
        return 0;
    }
    ($format) = $self->_map_formats($format);
    if ($format) {
        $ws->write($crow, $col, $val, $format);
    } else {
        $self->warn("Unknown format '$_[3]' in apply_format()");
        return 0;        
    }
    return $format;
}

=head2 add_row_explicit

 Title   : add_row_explicit
 Usage   : $eh->add_row_explicit( $sheet, $cells, $format );
 Function: Adds a single row to the spreadsheet
 Returns : The current row count of the sheet, or zero if no row added
    Args : [0] The worksheet to add the row to
           [1] An array reference of cell contents
           [2] Optional format specification

This function is similar to add_row(), but allows more complex cell
contents to be added to the worksheet. Each member (cell) of the row
array reference is considered seperately. If the member is a scalar,
it is added as a string - the same behavior as add_row().

If the cell is an array reference, it is assumed to have 2 or 3 parts:

  [0] The cell type (a string)
  [1] The primary cell value
  [2] Optional second value, used by some types

The recognized cell types are:

  'string' - The primary cell value is added as a forced string
  'number' - The primary cell value is added as an explicit number
     'url' - The primary value is the URL, the secondary one is used as
             an optional name. They cell will be a hyperlink.
     'bmp' - The primary value is a path to a windows bitmap file. The cell
             will contain an image

=cut

sub add_row_explicit {
    my $self = shift;
    my ($sheetname, $array, $format) = @_;
    my $ws = $self->sheet($sheetname);
    unless ($ws) {
        $self->warn("No worksheet with name '$sheetname' in add_row()");
        return 0;
    }
    my @keybits;
    foreach my $val (@{$array}) {
        if ( ref($val) && ref($val) eq 'ARRAY' ) {
            push @keybits, join('[ARR]', @{$val}) || "";
        } else {
            push @keybits, defined $val ? $val : "";
        }
    }

    my $key = join("\t", @keybits);
    return 0 if ($ws->{_CAT_}{DONE}{$key} && !$self->allow_duplicate_rows);
    $ws->{_CAT_}{DONE}{$key}++;
    my $crow = $ws->{_CAT_}{CURRENT_ROW};
    my @formats;
    
    if ($format && ref($format) eq 'ARRAY') {
        @formats = $self->_map_formats(@{$format});
    } else {
        # Same format in each cell
        ($format) = $self->_map_formats($format);
        @formats = map { $format } @{$array};
    }

    for my $col (0..$#{$array}) {
        my $val = $array->[$col];
        my $fmt = $formats[$col];
        if ( ref($val) && ref($val) eq 'ARRAY' ) {
            my @bits = @{$val};
            my ($type, $val1, $val2) = @bits;
            $type = lc($type);
            if ($type eq 'string') {
                $ws->write_string($crow, $col, $val1, $fmt);
            } elsif ($type eq 'number') {
                $ws->write_number($crow, $col, $val1, $fmt);
            } elsif ($type eq 'url') {
                $self->_write_url($ws, $crow, $col, $val1, $val2, $fmt);
            } elsif ($type eq 'bmp') {
                if (-e $val1) {
                    $ws->insert_bitmap($crow, $col, $val1);
                } else {
                    $val1 ||= "-UNDEF-";
                    $self->warn("Bitmap file '$val1' does not exist");
                }
            } else {
                $ws->write($crow, $col, $val1, $fmt);
            }
        } else {
            # Normal entry
            $ws->write($crow, $col, $val, $fmt);
        }
    }
    # print "<pre>$ws->{_CAT_}{NAME} $crow = $key</pre>";
    return ++$ws->{_CAT_}{CURRENT_ROW};
}

=head2 set_cell

 Title   : set_cell
 Usage   : $eh->set_cell( @params );
 Function: Adds a single row to the spreadsheet
 Returns : The current row count of the sheet, or zero if no row added
    Args : Associative array of arguments. Recognized keys:

    -sheet Required. The name of the worksheet to add the cell to

    -value The primary value to be used

      -col The column number of the cell, where 0 is the first column

      -row The row number of the cell, where 0 is the first row. If
           this value is left out, the most recently written row will
           be used.

     -val2 Used as the secondary value for type 'URL'

     -type If not specified, then the value is added as a generic
           value. Otherwise, can be one of:


This function allows a single value to be 'surgically' placed in the
spreadsheet. It is useful when one cell in a row has special
formatting and you do not wish to call add_row_explicit(), or if you
need to alter previously entered cells.

=cut

sub set_cell {
    my $self = shift;
    my $args = $self->parseparams( @_ );
    my $sheetname = $args->{SHEET};
    my $ws = $self->sheet($sheetname);
    unless ($ws) {
        $sheetname ||= '-UNDEF-';
        $self->warn("No worksheet with name '$sheetname' in set_cell()");
        return 0;
    }
    my $val1      = $args->{VALUE};
    my $format    = $self->_map_formats($args->{FORMAT});
    my $crow      = $args->{ROW};
    my $col       = $args->{COL};
    my $type      = lc($args->{TYPE} || "");
    unless (defined $crow) {
        $crow = $ws->{_CAT_}{CURRENT_ROW};
        $crow-- if ($crow);
    }
    unless (defined $col) {
        $self->warn("set_cell() did not define column. Using Col 0");
        $col = 0;
    }
    if (!$type) {
        # Normal entry
        $ws->write($crow, $col, $val1, $format);
    } elsif ($type eq 'string') {
        $ws->write_string($crow, $col, $val1, $format);
    } elsif ($type eq 'url') {
        my $val2 = $args->{VAL2} || $args->{LABEL};
        $self->_write_url($ws, $crow, $col, $val1, $val2, $format);
    } elsif ($type eq 'bmp') {
        if (-e $val1) {
            $ws->insert_bitmap($crow, $col, $val1);
        } else {
            $val1 ||= "-UNDEF-";
            $self->warn("Bitmap file '$val1' does not exist");
        }
    } else {
        $self->warn("Failed to understand cell type '$type'");
        $ws->write($crow, $col, $val1, $format);
    }
}

sub _write_url {
    my $self = shift;
    my ($ws, $crow, $col, $val1, $val2, $format) = @_;
    $format ||= $self->format('hyperlink');
    # Wow. They changed the param order...
    ($val2, $format) = ($format, $val2); # unless ($self->{isXLS});
    $ws->write_url($crow, $col, $val1, $val2, $format);
    $globals->{URLS_WRITTEN}++;
}

sub _map_formats {
    my $self = shift;
    my @forms;
    foreach my $req (@_) {
        if (defined $req) {
            push @forms, $self->format($req);
        } else {
            push @forms, undef;
        }
    }
    return @forms;
}

#sub err {
#    my $self  = shift;
#    my ($msg, $kill) = @_;
#    printf($self->{_CAT_}{PFORM}{warn}, $msg || "Undefined Error");
#    die " " if ($kill);
#}

=head2 file_path

 Title   : file_path
 Usage   : my $path = $eh->file_path( );
 Function: Gets the path of the created workbook
 Returns : The path (string)
    Args : None

Returns the file path, as initially set when new() was called.

=cut

sub file_path {
    return shift->{_CAT_}{PATHS}{file};
}

=head2 url

 Title   : url
 Usage   : my $link = $eh->url($newValue );
 Function: Gets / Sets a URL to the file
 Returns : The URL (string), if set
    Args : [0] Optional new value

If you are going to be using html_summary() to describe the new
workbook on a webpage, you should use this method to define the URL to
the file so the user will be able to open it over the web.

At BMS, web-accessible paths will correspond to URLs as follows:

    my $url = $eh->file_path( );
    $url =~ s/\/stf/http:\/\/bioinformatics.bms.com/;
    $eh->url($url);

Remember that if you are generating a file from a web page, it will be
owned by httpd. If the file is being accessed from the web, this
should be fine, if not you might want to chmod / chgrp / chown the
file as appropriate.

=cut

sub url {
    my $self = shift;
    if ($_[0]) {
        $self->{_CAT_}{URL} = $_[0];
    }
    return $self->{_CAT_}{URL};
}

=head2 html_summary

 Title   : html_summary
 Usage   : my $html = $eh->html_summary( );
 Function: Generates a little snippet of HTML summarizing the workbook
 Returns : A string of HTML text
    Args : None

Generates a little HTML table that summarizes the document, listing
the name of each worksheet and the number of rows in each. If you have
specified url(), a hyperlink to the file will also be shown.

=cut

sub html_summary {
    my $self = shift;
    my $string = "<table border='1'><tr><td bgcolor='#00ff00'>";
    my $msg = "Excel file summary";
    if (my $url = $self->url) {
        $msg = "<a href='$url'>Click to view your Excel file</a>";
    }
    $string .= "<b>$msg</b></td></tr>\n";
    $string .= "<tr><td align='center'><table><tr><th>Sheet</th><th>Rows</th><tr>";
    my $count = $self->row_count;
    my $total = 0;
    foreach my $sheet (sort keys %{$count}) {
        my $num = $count->{$sheet} - 1;
        $total += $num;
        next unless ($num);
        $string .= sprintf("<tr><td align='right'>%s:</td>".
                           "<td align='center'>%d</td></tr>\n",$sheet,$num);
    }
    $string .= sprintf("<tr><td align='right'><b>%s:</b></td>".
           "<td align='center'>%d</td></tr>\n",'Total',$total);
    
    $string .= "</table></td></tr></table>";
    return $string;
}

our $hexMap = {
    AQUA    => '#00FFFF',
    BLACK   => '#000000',
    BLUE    => '#0000FF',
    BROWN   => '#804000',
    CYAN    => '#00FFFF',
    FUCHSIA => '#FF00FF',
    GOLD    => '#D4A017',
    GRAY    => '#808080',
    GREY    => '#808080',
    GREEN   => '#008000',
    LIME    => '#00FF00',
    MAROON  => '#800000',
    NAVY    => '#000080',
    OLIVE   => '#808000',
    ORANGE  => '#FF8040',
    PINK    => '#FF00FF',
    PURPLE  => '#800080',
    RED     => '#FF0000',
    SILVER  => '#C0C0C0',
    TEAL    => '#008080',
    WHITE   => '#FFFFFF',
    YELLOW  => '#FFFF00',
};

our $hexRE = '[A-F0-9]';
sub hex_color {
    my $self = shift;
    my $col  = uc(shift || "");
    return undef unless ($col);
    # Known color name:
    $col = $hexMap->{$col} if ($hexMap->{$col});
    # Standard 6 digit hex code:
    if ($col =~ /^\#?($hexRE{6})$/) { return "#$1" }
    # 3 digit hex code:
    if ($col =~ /^\#?($hexRE)($hexRE)($hexRE)$/) { return "#$1$1$2$2$3$3" }
    if ($col =~ /rgb\(\s*(\d{1,3})\s*\,\s*(\d{1,3})\s*\,\s*(\d{1,3})\s*\)/) {
        # RGB decimal notation
        return sprintf("#%02X%02X%02X", $1, $2, $3);
    }
    # Hmm...
    return undef;
}

sub rgb_color {
    my $self = shift;
    my @rv;
    if (my $hex = $self->hex_color(@_)) {
        if ($hex =~ /^\#($hexRE{2})($hexRE{2})($hexRE{2})$/) {
            @rv = ( hex($1), hex($2), hex($3) );
        }
    }
    return @rv;
}

sub integer_to_alphabet {
    my $self = shift;
    my $pos  = shift || 0;
    my $rv = "";
    do {
        $pos--;
        my $mod = $pos % 26;
        $rv = chr(65 + $mod) . $rv;
        $pos -= $mod;
        $pos /= 26;
    } while ($pos > 0);
    return $rv;
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
1;

