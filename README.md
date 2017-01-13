# VSDP: Verified SemiDefinite Programming (Version 2006 by Christian Jansson)

VSDP is a [MATLAB](https://www.mathworks.com) /
[GNU Octave](https://www.gnu.org/software/octave) software package for
rigorously solving semidefinite programming problems.  It expresses these
problems in a notation closely related to the form given in textbooks and
scientific papers.  Functions for computing verified forward error bounds
of the true optimal value and verified certificates of feasibility and
infeasibility are provided.  All rounding errors due to floating point
arithmetic are taken into account.  Computational results are given,
including results for the
[SDPLIB benchmark problems](http://euler.nmt.edu/~brian/sdplib/sdplib.html).
This package supports interval input data and sparse format.


## Prerequisites

VSDP requires the approximate solver [SDPT3](https://github.com/sqlp/sdpt3)
and the [interval toolbox "INTLAB"](http://www.ti3.tuhh.de/rump/intlab) in
version 9 or higher.


## Some demonstration

Published output from
[demovsdp.m](https://rawgit.com/siko1056/vsdp-2006-ng/master/doc/demovsdp.html)


## References

- [Jansson2005] C. Jansson, Termination and Verification for Ill-posed
  Semidefinite Programming Problems,
  http://www.optimization-online.org/DB_HTML/2005/06/1150.html

- [Jansson2004] C. Jansson, Rigorous Lower and Upper Bounds in Linear
  Programming, SIAM J. OPTIM. Vol.14, No.3, pp. 914-935,
  https://dx.doi.org/10.1137/S1052623402416839
