
Python based Time Domain Moment Tensor Inversion
------------------------------------------------

Current Version
---------------
V2.4 Roma 08.08.2013
Fully tested on Mac OS 10.7.5

Download
--------
http://webservices.rm.ingv.it/pyTDMT

A. Install
----------

Requirements: . Obspy        see http://obspy.org  
              . rdseed-5.0   see http://www.iris.edu/forms/rdseed_request.htm
              . FKRPROG      into src/fkrprog of this package
                             Generate Frequency waves number for greens
              . earth-models see earth-models/ of this package of this package
                             for PREM,AK135 and Central Italy models
                         

FKRPROG:      cd src/fkrpog
              gfortran -O3 FKRPROG.f -o FKRPROG
              or
              g77 -ff90 -O3 FKRPROG.f -o FKRPROG
              mv FKRPROG ../../bin

Add Into your .profile or .bashrc or bash_profile depending on your OS:    
              PATH=$PATH:/pathTo/tdmt/bin
              PATH=$PATH:/pathTo/tdmt/py

B. Usage and Examples
---------------------

Usage: type
 
    pytdmt.py -h
    pytdmt.py --help

Examples:    cd examples of this package

C. Tutorial
-----------
http://webservices.rm.ingv.it/pyTDMT


D. references
-------------
1. DUNKIN J. W. (1965),BSSA,335-358   
2. HASKELL N. A. (1964),BSSA,337-393   
3. WANG AND HERRMANN (1980),BSSA,1015-1036   
4. WATSON T. H. (1970),BSSA,161-166.
5. http://seismo.berkeley.edu/~dreger/mtindex.html
   (Dreger and Helmberger, 1993; Dreger, 2003)



Contact fabrizio.bernardi@ingv.it
---------------------------------
