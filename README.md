MapLoc
======

Genetic variation database that unifies polymorphisms and
mutations. Captures information from populations, individuals and
samples, and provides interfaces for queries with a variety of report
formats. The system is species agnostic, using reference sequences
(generally a specific genome build) to scaffold and segregate data.

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

