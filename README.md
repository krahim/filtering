filtering
==========

R code implementing filtering tools from Percival and Walden (1993). This is based on their LISP code, and it uses FFTW and requires compiling. This is my current working version and there are hard coded path's in some of the R files that will have to be changed. 

I may migrate this to an appropriate R package, and it currently it uses a Makefile. I build this on on Linux, but you can build it in Window, or Mac using either Rtools or Fink if you set up fftw and set the path variables appropriately.

This requires PWLISPR and will overwrite some of the filtering functionality.
