
.. image:: https://img.shields.io/appveyor/ci/laurentnoe/yass/master.svg?style=flat-square&label=Build%20Status
    :target: https://ci.appveyor.com/project/laurentnoe/yass/
    :alt: Build Status

.. image:: https://img.shields.io/website.svg?style=flat-square&label=Website&url=http%3A%2F%2F193.54.101.81%2Fyass%2F
    :target: https://bioinfo.univ-lille.fr/yass/
    :alt: Website

yass
====

(more at  http://bioinfo.univ-lille.fr/yass/)

``yass`` is a genomic similarity seach tool for nucleic (and only
nucleic) sequences in (multi)fasta or plain text format. ``yass``
produces local pairwise alignments in yass format, blast tabular
format, or PSL format.

The associated tool ``yass2blast.pl``  may be used to convert the
default yass output into blast full output, into fasta alignments, or
into AXT format.

The associated tool ``yass2dotplot.php`` could also be used to
convert the default yass output into  ``png`` dotplots.

Installation
------------

(more at  http://bioinfo.univ-lille.fr/yass/download.php)

You need a C compiler and the autotools. On Linux, you can install
``gcc``, ``autoconf``, ``automake``. On Mac, you can install
``xcode``, or the command line developer tools (or you can use
``macports`` to install ``gcc`` for example).


Using the command line, type::

  git clone https://github.com/laurentnoe/yass.git
  cd yass
  ./configure --with-threads
  make

or::

  git clone https://github.com/laurentnoe/yass.git
  cd yass
  autoreconf
  ./configure --with-threads
  automake
  make

you can install  ``yass`` to a standard ``/local/bin`` directory::

  sudo make install

or copy the binary directly to your homedir::
   
  cp src/yass ~/.

Command-line
------------

(more at  http://bioinfo.univ-lille.fr/yass/help.php)


common usage
~~~~~~~~~~~~

-d <N>
  where *N = [0..5]*, to select the output format (default is 1)

-r <N>
  where *N = [0..2]* to select the *forward*, *reverse*, or *both*
  sense on the first sequence (default is *both*)

-S <N>
  to select *only one* sequence in the first multifasta file (give a
  number between 1 and nbparts).
  
  By default *all the sequences* are processed.


scoring system
~~~~~~~~~~~~~~

-C <N,...>
  with 2,3,4 or 16 parameters to give the:
  
  - Match/Mismatch scores,
  - Match/Transition/Transversion scores,
  - Match/Transition/Transversion/Other IUPAC scores,
  - 4x4 ACGT matrix (and disable scoring correction algorithm).
  

-G <No,Ne>
  with two parameters to change the cost for:

  - the very first gap opening *No*,
  - the subsequent extension costs *Ne*.


-E <N>  to set the E-value threshold *N* (default 10).


-X <N>  to set  the X-drop threshold score *N* (default 25).

and

-L <Nl,Nk>
  to possibly change the Lambda *Nl* and K *Nk* values
  if the one computed do not correspond to your needs.
  (Note that the ALP tool can do the work:
  https://www.ncbi.nlm.nih.gov/CBBresearch/Spouge/html_ncbi/html/software/program.html?uid=6
  )


search parameters
~~~~~~~~~~~~~~~~~

-p <"seedpattern">
    where the seed pattern is one, or several seeds separated by
    comma, where each seed  is a word on the "#@-" alphabet
    
    (Note that the Iedera tool can do the design:
    https://github.com/laurentnoe/iedera
    or
    http://bioinfo.univ-lille.fr/yass/iedera.php
    )

    another possibility is to use "Minimally overlapping words"
    such as the pattern   "RYNNNNNnnnNNNN"   to speed-up
    the search, but at a lower sensitivity.
 
-c <N>
   where *N = [1..2]* for single or double hit criterion


   
  
Example
-------

A very small example::

  yass                   file1.fa  file2.mfa >  yass-output.yop
  yass2blast.pl  -blast  yass-output.yop     >  blastlike-output.blk
  yass2dotplot.php       yass-output.yop  filename1=""  filename2="" ; open dp.png


A second example where the scoring system is modified, the E-value changed::

  yass    file1.fa  file2.mfa    -C 2,-2,-3   -G -5,-2   -E 1e-3   -o yass-output.yop



  

References
----------

how to cite this tool:

    Noe L., Kucherov G., YASS: enhancing the sensitivity of DNA similarity search, 2005, Nucleic Acids Research, 33(2):W540-W543. <http://doi.org/10.1093/nar/gki478>

