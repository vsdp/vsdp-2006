%DEMOVSDP    A little demonstration of VSDP.
%
% written   11/21/06   Christian Jansson
% modified
% Reference: C. JANSSON, Termination and Verification for
%            Ill-posed Semidefinite Programming Problems,
%            to appear
% http://optimization-online.org/DBHTML/|2005/06/1150.html.

echo off
intvaldisplay = intvalinit('display');
close
clc
format compact
format short
intvalinit('displaymidrad')
echo on
clear
format long e
clc

% Welcome to VSDP: Verified SemiDefinite Programming
%
%     A software package computing verified results
%                  for SDP-problems
%
%         min  sum(j=1:n | <C{j}, X{j}>)
%         s.t. sum(j=1:n | <A{i,j}, X{j}>) = b(i) for i = 1 : m
%              X{j} must be positive semidefinite for j = 1 : n
%              Notation: <C, X> := trace(C * X)
%
% VSDP exploits the block-diagonal structure by cell-arrays:
% The j-th block C{j} and the blocks A{i,j} for i = 1,...,m
% are real symmetric matrices of common size s_j . This
% structure is described with the cell-array
%    blk{j,1} = 's', blk{j,2} = s_j

pause


% EXAMPLE: only one block


DELTA = 1e-4;

C{1}    = [ 0   1/2   0;
  1/2 DELTA 0;
  0    0   DELTA ];

pause

A{1,1} =  [ 0 -1/2 0;
  -1/2 0  0;
  0   0  0];
A{2,1} =  [ 1 0 0 ;
  0 0 0;
  0 0 0];
A{3,1} =  [ 0 0 1;
  0 0 0;
  1 0 0];
A{4,1} =  [ 0 0 0;
  0 0 1;
  0 1 0];

b = [1; 2*DELTA; 0; 0];

pause

blk{1,1} = 's'; blk{1,2} = 3;

pause

% It is easy to prove that for every
% DELTA > 0: the optimal value is -0.5, and strong duality holds,
% DELTA = 0: ill-posed problem with nonzero duality gap,
% DELTA < 0: primal and dual infeasible.


pause

clc

% Check the size and data types of the input data
[m, n] = vsdpcheck(blk,A,C,b)

% All checks are executed correctly
pause
clc

% VSDP can be used for computing approximations only.
% Different solvers can be chosen by setting the global
% variable VSDP_CHOICE_SDP in the file SDP_GLOBALPARAMETER.
% VSDP_CHOICE_SDP = 1:    SDPT3 is selected
% VSDP_CHOICE_SDP = 2:    SDPA is selected
%
% Other solvers can be integrated easily in MYSDPS by the user.

pause



% The call is [objt,Xt,yt,Zt,info] = mysdps(blk,A,C,b);



pause

[objt,Xt,yt,Zt,info] = mysdps(blk,A,C,b);

pause

% The primal and dual approximations


objt, X = Xt{1}, yt

pause
clc


% SDPT3: shows normal termination, weak duality is violated.
%
% SDPA: computes erroneous approximations together with a warning.
info(1)
pause
clc
% Rigorous lower bound with call  [fL, Y, dL] = vsdplow(blk,A,C,b,Xt,yt,Zt);

pause

[fL, Y, dL] = vsdplow(blk,A,C,b,Xt,yt,Zt);

pause

fL, Y, dL

% SDPT3: A close finite rigorous lower bound is computed.
% Y provides a rigorous certificate of dual
% feasibility (i.e. feasibility of LMI's).
%
% SDPA: The rigorous lower bound is -infinity,
% and no certificates are computed, i.e. Y = NaN.

pause
clc

% Now, the upper bound  [fU, X, lb] = vsdpup(blk,A,C,b,Xt,yt,Zt);

pause

[fU, X, lb] = vsdpup(blk,A,C,b,Xt,yt,Zt);


pause

% Output:

fU, X, lb

%  SDPT3: with these approximations we obtain a
%  close rigorous upper bound of the optimal value and
%  an interval enclosure of a primal interior point.
%
%  SDPA: approximations yield the upper bound +infinity.

pause
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Summarizing, we have verified for the above problem with the SDPT3
%  approximations strong duality and the inequality
%
%   -5.000000062262615e-001 <=  f* <=-4.999677693282600e-001,
%
%  Moreover, we have obtained strictly  primal and dual
%  feasible solutions. The upper and lower bounds of the optimal
%  value show a modest overestimation, mainly due to the accuracy of
%  SDPT3. The bounds are rigorous, and thus must satisfy weak
%  duality.
%
%
%  Comparing with the SDPA results,
%  we see that the quality of the rigorous results depends
%  strongly on the quality of the computed approximations
%
pause
% Here, the intention is not to compare these solvers
% (for a comparison see Mittelmann Math. Program. Ser. B, 2003),
% but we want to illustrate that the
% rigorous error bounds depend very much on the quality of
% the computed approximations.

pause
%  VSDP has further features, for example:
%       - Certificates of infeasibility
%       - Other data formats (sparse, interval)
%       - Verified solutions for underdetermined systems
%       - Verified eigenvalue bounds
%       - Conversion routines for SDPA and SDPT3 format
%       - Applications
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pause
clc

% Upper bounds for the primal variables sometimes accellerate
% the routine [fL, Y, dL] = vsdplow(blk,A,C,b,Xt,yt,Zt,xu)

pause

xu = 1e5;

pause

[fL, Y, dL] = vsdplow(blk,A,C,b,Xt,yt,Zt,xu)

% SDPT3: we obtain the same lower bound as before,
%
% SDPA: approximation yield a finite overestimated lower bound

pause
clc

%  Now, we suppose the existence of dual upper bounds

yu = 1e5 * [1 1 1 1]';

pause

[fU, X, lb] = vsdpup(blk,A,C,b,Xt,yt,Zt,yu)

% SDPT3: we obtain a reasonable bound as before,
%
% SDPA: approximation yield a finite overestimated upper bound

pause
clc

% The  following results are very similar for the
% SDPT3 and for the SDPA approximations
pause

% For negative DELTA our problem is primal and dual infeasible

DELTA = -1e-4;

C{1}   = [ 0   1/2   0 ;
  1/2 DELTA 0 ;
  0   0   DELTA ];

b = [1; 2*DELTA; 0; 0];

% Call of MYSDPS: [objt,Xt,yt,Zt,info] = mysdps(blk,A,C,b);

pause
clc

[objt,Xt,yt,Zt,info] = mysdps(blk,A,C,b);

% Lower bound: [fL, Y, dL] = vsdplow(blk,A,C,b,yt)

pause

[fL, Y, dL] = vsdplow(blk,A,C,b,Xt,yt,Zt)

% Upper bound: [fU, X, lb] = vsdpup(blk,A,C,b,Xt,yt,Zt)

pause

[fU, X, lb] = vsdpup(blk,A,C,b,Xt,yt,Zt)

pause
clc
% Notice that fL = -Inf and fU = Inf are the exact optimal values,
% since the problem is primal and dual infeasible.
%
% However, VSDLOW and VSDPUP only indicate infeasibility,
% but they do not proof this property.
%
% This is the task of VSDPINFEAS, which tries to verify
% primal or dual unbounded rays.
pause
% For proving primal infeasibility we use VSDPINFEAS as follows:
%    choose = 'p';
%    [isinfeas, X, Y] = vsdpinfeas(blk,A,C,b,choose,Xt,yt,Zt),

pause

choose = 'p';
[isinfeas, X, Y] = vsdpinfeas(blk,A,C,b,choose,Xt,yt,Zt),

% Since ISINFEAS = 0, VSDPINFEAS has not proved primal infeasibility.

pause
clear
% Reason: each dual improving ray has a zero eigenvalue, and
% positive semidefiniteness cannot be verified for matrices with zero
% eigenvalues. This is out of scope of VSDP (analogous result for
% dual infeasility).

pause
clc
% Another example:
C{1} = [0 0; 0 0];
A{1,1} = [1 0; 0 0];
A{2,1} = [0 1; 1 0.005];
b = [-0.01 1]';
blk{1,1} = 's'; blk{1,2} = 2;

% is primal infeasible, since the first primal equation implies
% X(1,1) = -0.01}. Hence, X cannot be positive semidefinite.

choose = 'p';
% Now we call VSDPINFEAS

pause
[isinfeas, X, Y] = vsdpinfeas(blk,A,C,b,choose),
% Since ISINFEAS = 1, primal infeasibility is verified,
% and Y is a rigorous improving ray
pause
clc
% SDPLIB benchmarks (here from Truss Topology Design)

pause

% Choose the directory containing the SDPLIB

vsdp_path = regexprep(which('demovsdp'),'demovsdp.m','');
sdplib_vsdp_path = strcat(vsdp_path,'sdplib_vsdp/');

sdplibfiles = dir(strcat(vsdp_path,'sdplib_vsdp/*.dat-s'));

pause

% Conversion from SDPA to VSDP format

[blk,A,C,b] = sdpa_to_vsdp...
  ([sdplib_vsdp_path sdplibfiles(3).name]);

Problem = sdplibfiles(3).name,

pause
% Dimensions
m = length(b),

pause
n = length(C),

pause
size(C{1}), size(C{2}),

pause

% Call of SDP solver:    [objt,Xt,yt,Zt,info] = mysdps(blk,A,C,b);
pause

[objt,Xt,yt,Zt,info] = mysdps(blk,A,C,b);

pause
% Lower bound: [fL, Y, dL] = vsdplow(blk,A,C,b,Xt,yt,Zt);

pause

[fL, Y, dL] = vsdplow(blk,A,C,b,Xt,yt,Zt);

pause
fL

pause
size(Y), dL,

pause

% Upper bound: [fU, X, lb] = vsdpup(blk,A,C,b,Xt,yt,Zt)

pause

[fU, X, lb] = vsdpup(blk,A,C,b,Xt,yt,Zt);
fU,

pause
X, dL

pause
clc
% Summarizing: ARCH4 is strictly feasible,
% strong duality holds with rigorous bounds

[fL; fU]

pause

clc


% The same SDPLIB problem
% with interval input data (small radius)

[blk,A,C,b] = sdpa_to_vsdp...
  ([sdplib_vsdp_path sdplibfiles(3).name]);
Problem = sdplibfiles(3).name


r=1e-8;
m = length(b);
n = length(C);
%    for j = 1 : n
%        CI{j} = midrad(C{j},r*abs(C{j}));
%        for i = 1 : m
%            AI{i,j} = midrad(A{i,j},r*abs(A{i,j}));
%        end
%    end
pause
clc

for j = 1 : n
  CI{j} = midrad(C{j},r*abs(C{j}));
  for i = 1 : m
    AI{i,j} = midrad(A{i,j},r*abs(A{i,j}));
  end
end

bI = midrad(b,r*abs(b));

pause
clc

% Call of SDPT3:    [objt,Xt,yt,Zt,info] = mysdps(blk,A,C,b);

pause
[objt,Xt,yt,Zt,info] = mysdps(blk,A,C,b);

%
%
% Lower bound: [fL, Y, dL] = vsdplow(blk,AI,CI,bI,Xt,yt,Zt)

pause

[fL, Y, dL] = vsdplow(blk,AI,CI,bI,Xt,yt,Zt);

%
%
fL, size(Y), dL

%
%
%
% Upper bound [fU, X, lb] = vsdpup(blk,AI,CI,bI,Xt,yt,Zt);

pause

[fU, X, lb] = vsdpup(blk,AI,CI,bI,Xt,yt,Zt);

fU, size(X), lb

pause
clc

% Summarizing: For all data within the intervals
% ARCH4 is strictly feasible,
% strong duality holds with rigorous bounds

[fL; fU]


return
