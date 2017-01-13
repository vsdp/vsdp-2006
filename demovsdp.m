%% Demonstration of VSDP: Verified SemiDefinite Programming
%
% VSDP is a software package computing verified results for SDP-problems
%
% $$
% \begin{equation}
% \begin{array}{lll}
% \text{minimize}
% & \sum\limits_{j=1}^{n} \langle C_{j}, X_{j}\rangle, & \\
% \text{subject to}
% & \sum\limits_{j=1}^{n} \langle A_{i,j}, X_{j}\rangle = b_{i},
% & i = 1, \ldots, m, \\
% & X_{j} \succeq 0,
% & j = 1, \ldots, n,
% \end{array}
% \end{equation}
% $$
%
% where $\langle C, X\rangle := trace(CX)$.
%%

%% Installation
%
% VSDP is completely written in MATLAB and uses
%
% * the the MATLAB-toolbox INTLAB (INTerval LABoratory), which is available
%   from <http://www.ti3.tu-harburg.de/rump/intlab>,
% * the semidefinite solver SDPT3 <http://www.math.cmu.edu/~reha/sdpt3.html>,
%   a more recent version is available from <https://github.com/sqlp/sdpt3>,
% * or alternatively the semidefinite solver SDPA
%   <http://sdpa.sourceforge.net/>.
%
% The VSDP source files can be downloaded from
% <http://www.ti3.tuhh.de/jansson/vsdp2006> or a more recent version from
% <https://github.com/siko1056/vsdp-2006-ng>.
% The installation was successful if this script |demovsdp.m| runs through.
%



%% First example
%
% VSDP exploits the block-diagonal structure by cell-arrays:
%
% The j-th block |C{j}| and the blocks |A{i,j}| for |i = 1,...,m|
% are real symmetric matrices of common size |s_j|.
%
% This structure is described with the cell-array
%
%   blk(j,:) = {'s'; s_j};
%
% For the purpose of illustration, we start with the following semidefinite
% programming problem of dimension |m = 4|, |n = 1|, and |s_1 = 3|, i.e. the
% matrices consist of only one block.  The problem depends on a fixed parameter
% |DELTA|:
%

DELTA = 1e-4;

C{1} = [ 0   1/2    0;
        1/2 DELTA   0;
         0    0   DELTA];

A{1,1} = [  0  -1/2 0;
          -1/2   0  0;
            0    0  0];

A{2,1} = [1 0 0;
          0 0 0;
          0 0 0];

A{3,1} = [0 0 1;
          0 0 0;
          1 0 0];

A{4,1} = [0 0 0;
          0 0 1;
          0 1 0];

b = [1; 2*DELTA; 0; 0];

blk(1,:) = {'s'; 3};


%%
% It is easy to prove that for every
%
% * |DELTA > 0|: the optimal value is -0.5, and strong duality holds,
% * |DELTA = 0|: ill-posed problem with nonzero duality gap,
% * |DELTA < 0|: primal and dual infeasible.
%



%% |vsdpcheck| - check the problem data
%
% The routine |vsdpcheck| may be used to perform a check whether there are
% inconsistencies in the input data:
%

format short
[m, n] = vsdpcheck(blk,A,C,b)

%%
% Since no error messages occur, all checks w.r.t. the sizes and the data
% types of the input data are carried out correctly.
%



%% |mysdps| - computing approximations only
%
% VSDP can be used for computing only approximations with different solvers.
% The user can integrate other solvers very easily in the file |mysdps.m|.
% By default, the function |mysdps| calls the semidefinite programming solver
% SDPT3.  Different solvers can be chosen by setting the global variable
% |VSDP_CHOICE_SDP| in the file |SDP_GLOBALPARAMETER.m|.
%
% * |VSDP_CHOICE_SDP = 1|: SDPT3 is selected (default)
% * |VSDP_CHOICE_SDP = 2|: SDPA is selected
%
% We want to illustrate that the proposed rigorous error bounds depend very
% much on the quality of the computed approximations.
%

[objt,Xt,yt,Zt,info] = mysdps(blk,A,C,b);

%%
% The output consists of approximations of
%
% # the primal and dual optimal value both stored in |objt|,
% # the primal and dual solutions |Xt|, |yt|, |Zt|, and
% # information about termination and performance stored in |info|.
%

format longE
objt, X = Xt{1}, yt

%%
%

format short; info

%%
% The solver terminates without any warning if |info = 0|.  The first four
% decimal digits of the primal and dual optimal value are correct, but week
% duality is not satisfied since the approximate primal optimal value is
% smaller than the dual one.  In other words, the algorithm is not backward
% stable for this example.  Hence, for the nonzero coefficients about four
% decimal digits are correct.
%
% SDPA computes erroneous approximations together with a warning.
%



%% |vsdplow| - rigorous lower bounds
%
% To obtain more reliability we can use the function |vsdplow| which computes
% a verified lower bound of the primal optimal value by using previously
% computed approximations |Xt|, |yt|, and |Zt|.  There are no assumptions about
% the quality of |yt|.  First, we treat the situation where all components of
% |xu| are infinite (i.e. finite bounds of the maximal eigenvalues of a primal
% optimal solution |Xt| are not known or are not available).  The call has the
% form:
%

format longE
[fL, Y, dl] = vsdplow(blk,A,C,b,Xt,yt,Zt)

%%
% The output |fL|, |Y|, and |dl| corresponds to the lower bound, the
% certificate of dual feasibility, and the corresponding vector of eigenvalue
% bounds, respectively.  In the case where |Y = NaN|, all components of the
% vector |dl| are set |NaN|.
%
% With SDPT3 a finite rigorous lower bound close to the optimal value together
% with a certificate of dual feasibility is computed.  Since |dl| is positive,
% the dual problem contains the strictly dual feasible interior point |Y|,
% i.e. the dual Slater condition holds true.  Therefore, strong duality is
% verified, and |fL| is also a lower bound of the dual optimal value.
%
% Secondly, we describe the case, where additional information about finite
% upper bounds |xu| of the maximal eigenvalues of a primal optimal solution
% |Xt| is available.
%
% If we use, as above, the upper bound |xu = 1e5|, then for this example we
% obtain the same output:
%

%%
% *Note:* tiny pertubation of |yt| neccessary to work with the approximation 
% from SDPT3-4.0 for this particular example.
%
%   yt(2) = yt(2) - DELTA;
%

xu = 1e5;
[fL, Y, dl] = vsdplow(blk,A,C,b,Xt,yt,Zt,xu)

%%
% We see that the quality of the rigorous results depends strongly on the
% quality of the computed approximations, which is an immediate consequence
% of the proposed error bound.
%
% Using SDPA, the rigorous lower bound is minus infinity and no certificates
% are computed, i.e. |Y = NaN|.
%



%% |vsdpup| - rigorous upper bounds
%
% Similarly, with function |vsdpup| we can compute a verified upper bound |fU|
% of the dual optimal value using the previously computed approximations.
% There are no assumptions about the quality of the primal approximation |Xt|.
%
% First, we treat the case where all components of |yu| are infinite and
% describe shortly |vsdpup|.  The primal equations
% $\sum_{j=1}^{n} \lbrace A_{ij}, X_{j} \rbrace = b_{i}$, $(i = 1,\ldots,m)$
% are solved with a linear interval solver.  If an enclosure |X| containing an
% exact solution of the primal equations can be computed, it is automatically
% proved that the linear system of equations has full rank.  Then, with Weyl's
% Perturbation Theorem (c.f. routine |veigsym|) verified lower bounds of the
% eigenvalues for the blocks in |X| are computed and stored in a vector |lb|.
% If |lb| is nonnegative, then |X| contains a primal feasible solution.
% If |lb| is positive, then |X| contains a strictly primal feasible solution;
% that is, the primal Slater condition holds true.  If |lb| has negative
% components, then |vsdpup| tries to find iteratively new primal approximations
% by solving appropriate perturbed semidefinite problems, until either |lb| is
% nonnegative, or primal infeasibility is indicated by the semidefinite solver
% for the perturbed problem.  In the latter case |vsdpup| sets |fU = +inf|,
% and the certificate of primal feasibility |X| and the vector |lb| are defined
% as |NaN|.  The call of |vsdpup| has the form:
%

[fU, X, lb] = vsdpup(blk,A,C,b,Xt,yt,Zt)

%%
% Hence, we obtain a lower bound of the optimal value together with the
% interval matrix |X| containing a strictly feasible primal solution.
% The radius of the first component is of order |10e−20|.  Since, |lb > 0|
% Slater's condition is fulfilled and strong duality holds valid.
%
% Using SDPA's approximations, we yield the upper bound plus infinity.
%
% Summarizing, we have verified for the above problem with the SDPT3
% approximations strong duality and the inequality
%
%   -5.000000062262615e-001 <=  f* <= -4.999677693282600e-001,
%
% is fulfilled, and certificates of strictly primal and strictly dual feasible
% solutions are computed.  The Strong Duality Theorem implies that the primal
% and the dual problem have a nonempty compact set of optimal solutions.
% The upper and lower bounds of the optimal value show a modest overestimation,
% mainly due to the accuracy of SDPT3.  The bounds are rigorous, and thus
% satisfy weak duality.
%
% In the case where finite upper bounds |yu| are available, |vsdpup| computes
% rigorously the finite upper bound |fU| without verifying primal feasibility.
% Hence, solving the primal equations with a linear interval solver is not
% necessary, yielding an acceleration in most cases.  The fact that primal
% feasibility is not verified is expressed in |vsdpup| by |NaN| for |X| and
% |lb|.  The call of |vsdpup| with using finite upper bounds has the form
%

yu = 1e5 * [1 1 1 1]';
[fU, X, lb] = vsdpup(blk,A,C,b,Xt,yt,Zt,yu)

%%
% For this example the upper bound is worse than the previous one.
%



%% |vsdpinfeas| - certificates of infeasibility
%
% The following results are very similar for the SDPT3 and for the SDPA
% approximations.  For this example the upper bound is worse than the previous
% one.  For negative |DELTA| our problem is primal and dual infeasible.
%

DELTA = -1e-4;

C{1} = [ 0   1/2    0;
        1/2 DELTA   0;
         0    0   DELTA];

b = [1; 2*DELTA; 0; 0];

%%
% The routine SDPT3
%

[objt,Xt,yt,Zt,info] = mysdps(blk,A,C,b);

%%
%

format short; info

%%
% gives the termination code 1 with the warning primal infeasibility has
% deteriorated too much.  In accordance the rigorous lower bound is
%

format longE
[fL, Y, dL] = vsdplow(blk,A,C,b,Xt,yt,Zt)

%%
% and for the upper bound we obtain
%

[fU, X, lb] = vsdpup(blk,A,C,b,Xt,yt,Zt)

%%
% Notice that these are the exact optimal values, since the problem is primal
% and dual infeasible.  However, |vsdplow| and |vsdpup| only indicate
% infeasibility, but do not proof this property (both certificates are |NaN|).
% This is the task of |vsdpinfeas| which tries to verify primal or dual
% unbounded rays that deliver certificates of dual or primal infeasibility.
% For proving primal infeasibility we can use |vsdpinfeas| as follows:
%

choose = 'p';
[isinfeas, X, Y] = vsdpinfeas(blk,A,C,b,choose,Xt,yt,Zt);

format short; isinfeas
format longE; X, Y

%%
% The first line declares that we want to prove primal infeasibility.  In the
% second line the input of |vsdpinfeas| contains the data of the problem and
% the selection variable choose as well as the previously computed SDPT3
% approximations |Xt|, |yt|, and |Zt|.  Since |isinfeas = 0| the routine
% |vsdpinfeas| has not proved primal infeasibility for this example, and the
% certificates |X|, |Y| are set to |NaN|.
%
% The reason is that each dual improving |Y| ray has a zero eigenvalue, and
% positive semidefiniteness cannot be verified for matrices with zero
% eigenvalues.  This is out of scope of VSDP (analogous results are obtained
% for dual infeasility).
%

%%
% The example
%

clear C A b blk

C{1} = [0 0; 0 0];
A{1,1} = [1 0; 0 0];
A{2,1} = [0 1; 1 0.005];
b = [-0.01 1]';
blk(1,:) = {'s'; 2};

%%
% is primal infeasible, since the first primal equation implies
% |X(1,1) = -0.01|.  Hence, |X| cannot be positive semidefinite.
%

choose = 'p';
[isinfeas, X, Y] = vsdpinfeas(blk,A,C,b,choose);

format short; isinfeas
format longE; X, Y

%%
% Now, |isinfeas = 1| shows that the problem is primal infeasible and a
% rigorous certificate is provided by the dual improving ray |Y|.  Observe
% that vsdpinfeas can be called also without the optional input |Xt|, |yt|,
% and |Zt|.  Then appropriate approximations are generated in this routine,
% with the consequence of more computational effort.  Since this example is
% dual feasible, proving dual infeasibility should not be successful:
%

choose = 'd';
[isinfeas, X, Y] = vsdpinfeas(blk,A,C,b,choose);

format short; isinfeas
format longE; X, Y


%% SDPLIB benchmarks (here from Truss Topology Design)
%
% For computing verified results of the SDPLIB collection it is required that
% these problems are encoded in SDPA sparse format and stored in an appropriate
% directory.  This version of VSDP supplies the problem instances inside the
% directory |sdplib_vsdp|.  In this example we take the third SDPLIB problem
% which is |arch4.dat-s| from truss topology design:
%

arch4_path = ...
  fullfile(fileparts(which('demovsdp')), 'sdplib_vsdp', 'arch4.dat-s');

%%
% The routine |sdpa_to_vsdp| reads such a problem in SDPA sparse format, and
% converts it to VSDP format:
%

[blk,A,C,b] = sdpa_to_vsdp(arch4_path);

%%
% The size of the problem can be obtained with
%

format short
s_j = blk(:,2)
m   = length(b)


%%
% Hence, this is a block diagonal problem consisting of two symmetric blocks
% with dimensions 161 and 174 (yielding 28266 variables), and it has |m = 174|
% constraints.  At first, we solve this problem approximately with SDPT3, and
% then we compute the rigorous error bounds for the optimal value.
%

[objt,Xt,yt,Zt,info] = mysdps(blk,A,C,b);
[fL, Y, dl] = vsdplow(blk,A,C,b,Xt,yt,Zt);
[fU, X, lb] = vsdpup (blk,A,C,b,Xt,yt,Zt);

format longE
objt, fL, fU

%%
% The bounds for the optimal value provide a guaranteed accuracy of about eight
% decimal digits for this example.  The time for computing the lower bound is
% small compared to the times for the approximations and for the upper bound.
% This behavior is very typical.  There are two reasons.  Firstly, solving the
% primal equations with a linear interval solver is time consuming.  Secondly,
% solving perturbed semidefinite programming problems during the iteration is
% also time consuming.
%
% The computed quantities |X| and |Y| are both unequal |NaN| with positive
% minimal eigenvalue bounds
%

min(lb), min(dl)

%%
% Hence, for |arch4.dat-s| strictly primal and dual feasible solutions are
% verified, implying strong duality and compactness of the sets of optimal
% solutions.
%


%% VSDP with interval input data
%
% If the input data of a semidefinite programming problem are intervals then
% the rigorous bounds computed by VSDP hold valid for all problems with real
% input data inside the intervals.  In the following, we solve again the Truss
% Topology Design Problem |arch4.dat-s| with interval input data:
%

[blk,A,C,b] = sdpa_to_vsdp(arch4_path);

%%
%

r = 1e-8;
m = length(b);
n = length(C);

for j = 1:n
  CI{j} = midrad(C{j},r*abs(C{j}));
  for i = 1:m
    AI{i,j} = midrad(A{i,j},r*abs(A{i,j}));
  end
end

bI = midrad(b,r*abs(b));

%%
% Compute an approximate solution using SDPT3 and a rigorous lower and upper
% bound.
%

[objt,Xt,yt,Zt,info] = mysdps(blk,A,C,b);
[fL, Y, dL] = vsdplow(blk,AI,CI,bI,Xt,yt,Zt);
[fU, X, lb] = vsdpup (blk,AI,CI,bI,Xt,yt,Zt);

fL, fU

%%
% For |arch4.dat-s| the optimal value for all problems with real input data
% inside the intervals is between |fL| and |fU|.  Hence, a relative
% perturbation radius of $10^{-8}$ yields a relative radius of
% $4.9 \cdot 10^{-5}$ for the set of optimal values.  Looking at the other
% output data it follows that |Y| is a strictly dual feasible solution for all
% real problems, and the interval quantity |X| contains for each real problem
% a strictly primal feasible solution.  Especially, strong duality holds valid
% for all real problems.
%
% VSDP has further features, for example:
%
% * Certificates of infeasibility
% * Other data formats (sparse, interval)
% * Verified solutions for underdetermined systems
% * Verified eigenvalue bounds
% * Conversion routines for SDPA and SDPT3 format
% * Applications
%

%% Reference
%
% _C. Jansson_, *Termination and Verification for Ill-posed Semidefinite
% Programming Problems*,
% <http://www.optimization-online.org/DB_HTML/2005/06/1150.html>
%

% Copyright 2004-2006 Christian Jansson (jansson@tuhh.de)
