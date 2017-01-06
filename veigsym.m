function lambda = veigsym(A)
%VEIGSYM verified enclosure (i.e. interval vector) LAMBDA containing
%        all eigenvalues of a symmetric real or interval matrix A.
%        A may be real or interval, and full or sparse.

%
% Example:
% A = [1 2 3;
%      2 1 4;
%      3 4 5];
% A=midrad(A, 0.01*ones(3));
% lambda = veigsym(A)
% intval lambda =
% [   -1.5170,   -1.4569]
% [   -0.6226,   -0.5625]
% [    9.0495,    9.1096]

% written   11/15/05   Christian Jansson
% modified  12/28/05
% Reference: C. JANSSON, Termination and Verification for
%            Ill-posed Semidefinite Programming Problems,
%            to appear

if ~min(min( mid(A) == mid(A)'))
  error('VEIGSYM: matrix must be symmetric')
end

% Main routine using Weyl's Perturbation Theorem
% see Lecture Notes on Optimization with Result Verification
[V,D] = eig(full(mid(A)));
E = A - V * intval(D) * V';
r = abss(norm(E,inf));
lambda = midrad(diag(D),r);
