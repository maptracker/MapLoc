### MapLoc Variant Database

Genetic variation database that unifies polymorphisms and
mutations. Captures information from populations, individuals and
samples, and provides interfaces for queries with a variety of report
formats. The system is species agnostic, using reference sequences
(generally a specific genome build) to scaffold and segregate data.

![Graphical report][MapLoc]

### Utility modules

This release also includes several useful utility Perl modules that
are used by many of my other packages:

* [Utilities.pm][Utilities] - Core module functions, like argument
  parsing and messaging
* [ErrorInterceptor.pm][ErrorInterceptor] - Advanced error management
  via `%SIG`. All my web-served code was tied into this framework,
  which allowed errors (including fatal ones) to send me an email with
  stack trace (via `Debug.pm`) and user context.
* [ArgumentParser.pm][ArgumentParser] - Generalizes argument parsing,
  allowing the same code to function as a web CGI (pulling arguments
  from POST or GET via `CGI`) or from the command line (pulling them
  from `@ARGV`). The switch is made automatically based on inspection
  of the environment. Also handles setting mimetype and parsing of
  parameter files.
* [ExcelHelper.pm][ExcelHelper] - Aids in generation of `xls` and
  `xlsx` spreadsheets, including workbook/worksheet layout, column
  formating, and cell styling. [Example usage][ExEH].
* [ForkCritter.pm][ForkCritter] - Wrapper to abstract handling of
  forked jobs and collection / aggregation of results. [Example][ExFork]
* [FriendlyDBI.pm][FriendlyDBI] - A heavy-weight database helper
  module based around `DBI`. Allows creation of database schemae based
  on hash structures, as well as easy definition of statement handles
  and results recovery. Examples:
  * [Defining MapLoc DB schema][ExDbiSchema], including views and
    `plpgsql` proceedures
  * [Defining a cached statement handle][ExDbiSTH], including
    (commented out) pretty-printing and database EXPLAIN
  * [Hand-building a statement][ExDbiPlan] (no placeholders) when the
    query planner is unhappy (sometimes placeholders, particularly in
    conjunction with `>` or '<', cause a SEQSCAN rather than use of an
    index)
* [FriendlyGraph.pm][FriendlyGraph] - A lightweight graph/network database
* [FriendlySAX.pm][FriendlySAX] - Wrapper around
  `XML::Parser::PerlSAX` that manages callbacks and node
  selection. Here's [an example][ExSAX] that's consuming SAX nodes via
  callback from `ForkCritter`
* [TableReader.pm][TableReader] - Abstracted table reader from TSV,
  CSV and Excel. Includes utilities to manage headers and
  columns. [Example][ExTR] reading GWAS files
* [Benchmark.pm][Benchmark] - Uses `Time::HiRes` to provide
  high-resolution benchmarking utilities. [Example][ExBench] timing a
  SQL query.
* [ColorUtilities.pm][ColorUtilities] - Color management, including a
  function to [hash a string to a color][ExColor], used to colorize
  some output.
* [Debug.pm][Debug] - Used to serialize nested structures (Hash/Array)
  for debugging purposes. Code includes an example with output. Has
  options to prevent explosive recursion, or to exclude/focus on
  certain elements.
* [FileUtilities.pm][FileUtilities] - Utilities to manage files and
  directories, including Perl module locations.
* [SequenceUtilities.pm][SequenceUtilities] - A whole bunch of utility
  functions. Leverages BioPerl's `Bio::` utilities, but also has code
  to manage frustrating [API changes or inconsistencies][ExSeq]
* [Escape.pm][Escape] - Just escaping code for URLs, XML, JSON,
  etc. There's probably better packages available now.
* [Serialize.pm][Serialize] - Used to generate JSON. There are
  *definitely* better options available now.

Installation
------------

To install MapLoc, from this directory run the following script to
determine any public modules that need to be installed on your system:

./BMS/SnpTracker/MapLoc/setupMapLoc.pl

If all dependencies are available, the script will also create an
empty MapLoc schema in Postgres.

Notes
-----

*   2014 Apr 5

    Initial commit. This system is nicely functional within the
    sheltered corporate computing environment from whence it came. It
    will require additional work to make it fucnction on other
    systems. I will be installing it at home to identify presumptions
    that need to be changed, and committing fixes as I go.

Copyright
---------

&copy; 2014 Charles Tilford

 http://mit-license.org/

* [BMS Public Disclosure approval](PubD-Disclosure-Approval.md)

[MapLoc]: img/MapLoc.png

[ArgumentParser]: BMS/ArgumentParser.pm
[ErrorInterceptor]: BMS/ErrorInterceptor.pm
[ExcelHelper]: BMS/ExcelHelper.pm
[ForkCritter]: BMS/ForkCritter.pm
[FriendlyDBI]: BMS/FriendlyDBI.pm
[FriendlyGraph]: BMS/FriendlyGraph.pm
[FriendlySAX]: BMS/FriendlySAX.pm
[TableReader]: BMS/TableReader.pm
[Utilities]: BMS/Utilities.pm
[Benchmark]: BMS/Utilities/Benchmark.pm
[ColorUtilities]: BMS/Utilities/ColorUtilities.pm
[Debug]: BMS/Utilities/Debug.pm
[Escape]: BMS/Utilities/Escape.pm
[FileUtilities]: BMS/Utilities/FileUtilities.pm
[SequenceUtilities]: BMS/Utilities/SequenceUtilities.pm
[Serialize]: BMS/Utilities/Serialize.pm

[ExEH]: https://github.com/maptracker/MapLoc/blob/2f12b24705596ba5a197d6571666d8fff5fdf02d/BMS/SnpTracker/MapLoc/Reporter.pm#L1331
[ExFork]: https://github.com/maptracker/MapLoc/blob/2f12b24705596ba5a197d6571666d8fff5fdf02d/BMS/SnpTracker/MapLoc/load_dbSNP.pl#L241
[ExSAX]: https://github.com/maptracker/MapLoc/blob/2f12b24705596ba5a197d6571666d8fff5fdf02d/BMS/SnpTracker/MapLoc/load_dbSNP.pl#L578
[ExTR]: https://github.com/maptracker/MapLoc/blob/2f12b24705596ba5a197d6571666d8fff5fdf02d/BMS/SnpTracker/MapLoc/load_NHGRI_GWAS.pl#L113-L120
[ExDbiSchema]: https://github.com/maptracker/MapLoc/blob/2f12b24705596ba5a197d6571666d8fff5fdf02d/BMS/SnpTracker/MapLoc.pm#L2901
[ExDbiSTH]: https://github.com/maptracker/MapLoc/blob/2f12b24705596ba5a197d6571666d8fff5fdf02d/BMS/SnpTracker/MapLoc.pm#L4132-L4135
[ExDbiPlan]: https://github.com/maptracker/MapLoc/blob/2f12b24705596ba5a197d6571666d8fff5fdf02d/BMS/SnpTracker/MapLoc.pm#L6636-L6662
[ExBench]: https://github.com/maptracker/MapLoc/blob/2f12b24705596ba5a197d6571666d8fff5fdf02d/BMS/SnpTracker/MapLoc.pm#L5865-L5881
[ExColor]: https://github.com/maptracker/MapLoc/blob/master/BMS/Utilities/ColorUtilities.pm#L320
[ExSeq]: https://github.com/maptracker/MapLoc/blob/2f12b24705596ba5a197d6571666d8fff5fdf02d/BMS/Utilities/SequenceUtilities.pm#L801-L814
