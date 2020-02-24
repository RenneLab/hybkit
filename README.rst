
hybkit
==================================

Welcome to *hybkit*, a toolkit for analysis of ".hyb" format genomic sequence data.

This project contains multiple major components:
    #. A toolkit of command-line utilities for manipulating,
       analyzing, and plotting data contained within hyb-format files. (ToDo)
    #. Several template analysis piplelines providing example analyses for hybrid sequence data.
    #. The hybkit python package, which contains an extendable codebase and documented API
       for creation of custom analyses of hyb-format data.


Full project documentation is available at LINKHERE


Installation
------------

Hybkit requires Python 3.6+ and the use of the matplotlib package.

The recommended installation method is via 
`python pip <https://pip.pypa.io/en/stable/>`_, which will 
automatically handle version control and dependency installation::
    
    $ pip install hybkit

Acquisition of the package can also be performed by cloning the project's Github repository::

    $ git clone git:://github.com/RenneLab/hybkit

Or by downloading the zipped package::

    $ curl -OL https://github.com/dstrib/hybkit/archive/master.zip
    $ unzip master.zip

Followed by installation using python's setuptools::

    $ python setup.py install


