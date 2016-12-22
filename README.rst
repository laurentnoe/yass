yass
====

(more at  http://bioinfo.lifl.fr/yass/)

``yass`` is a genomic similarity seach tool for nucleic (and only
nucleic) sequences in (multi)fasta or plain text format. ``yass``
produces local pairwise alignments in yass format, blast tabular
format, or PSL format.

The associated tool ``yass2blast.pl``  may be used to convert the
default yass output into blast full output, into fasta alignments, or
into AXT format.


Installation
------------

(more at  http://bioinfo.lifl.fr/yass/download.php)

You need a C compiler and the autotools. On Linux, you can install
``gcc``, ``autoconf``, ``automake``. On Mac, you can install
``xcode``, or the command line developer tools (or you can use
``macports`` to install ``gcc-mp-5`` for example).


Using the command line, type::

  git clone https://github.com/laurentnoe/yass.git
  cd yass
  ./configure
  make
  
you can install  ``yass`` to a standard ``/local/bin`` directory::

  sudo make install

or copy the binary directly to your homedir::
   
  cp src/yass ~/.

Command-line
------------

(more at  http://bioinfo.lifl.fr/yass/help.php)

The commonly used command-line parameters are :

-d N
  where N = [0..5], to select the output format (default is 1)

-r N
  where N = [0..2] to select the *forward*, *reverse*, or *both*
  sense on the first sequence (default is *both*)

 
For the scoring system:

-C N(,N...)
   with 2,3,4 or 16 parameters to give the
   * Match/Mismacth scores,
   * Match/Transition/Transversion scores,
   * Match/Transition/Transversion/Other IUPAC scores,
   * 4x4 ACGT matrix (and disable scoring correction algorithm).

 
-G NO,NE
    with two parameters to change the cost for the first gap opening NO,
    and subsequent extension costs NE.


-E N  to set the E-value threshold N (default 10).


-X N  to set  the X-drop threshold score N (default 25).


Example
-------

A very small example::

  yass           file1.fa  file2.mfa       >  yass-ouput.yop
  yass2blast.pl  -blast    yass-ouput.yop  >  blastlike-output.blk


A second example where the scoring system is modified, the E-value changed::

  yass           file1.fa  file2.mfa       -C 2,-2,-3   -G -5,-2   -E 1e-3 -o  yass-ouput.yop



  

