# VSDP: Verified SemiDefinite Programming (Version 2006 by Christian Jansson)

> The latest version is [VSDP 2018](https://github.com/vsdp/vsdp-2018).

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

- VSDP requires the approximate solver [SDPT3](https://github.com/sqlp/sdpt3)
  and the interval toolbox [INTLAB](http://www.ti3.tuhh.de/rump/intlab).


## Some demonstration

- Published output from
  [demovsdp.m](https://rawgit.com/vsdp/vsdp-2006/master/doc/demovsdp.html)


For more information, please read the Manual
[vsdp-2006-doc.pdf](https://github.com/vsdp/vsdp-2006/raw/master/doc/vsdp-2006-doc.pdf)
or visit <https://vsdp.github.io/>.
