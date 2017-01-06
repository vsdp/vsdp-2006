function [fL, Y, dl] = vsdplow(blk,A,C,b,Xt,yt,Zt,xu)
% VSDPLOW Verified Semidefinite Programming Lower Bound
%         for the minimum value of the block-diagonal problem
%
%         min  sum(j=1:n| <C{j}, X{j}>)
%         s.t. sum(j=1:n| <A{i,j}, X{j}>) = b(i) for i = 1 : m
%              X{j} must be positive semidefinite for j = 1 : n
%
%         Moreover, a certificate of feasibility for LMI's is provided.
%
% The block-diagonal structure is described
% by an n*2-cell-array blk, n-cell-arrays C, X, and an
% m*n-cell-array A as follows:
% The j-th block C{j} and the blocks A{i,j} for i = 1 : m
% are real symmetric matrices of common size s_j, and
%    blk{j,1} = 's', blk{j,2} = s_j
% The blocks C{j} and A{i,j} must be stored as individual
% matrices in dense or sparse format.
%
% The input A, C and the m-vector b
% may be floating-point or interval quantities,
%
% Xt,yt,Zt    approximate floating-point solution or initial
%             starting point for the primal (this is Xt) and
%             the dual problem (this is yt, Zt).
%
% [fL, Y, dl] = VSDPLOW(blk,A,C,b,Xt,yt,Zt) returns
%        fL   a rigorous lower bound of the minimum value
%             for all real input data (C,A,b) within the
%             interval input data.
%             fL = -inf, if no finite rigorous lower bound
%             can be computed.
%        Y    =NAN and dl=repmat(NaN,n,1), if dual feasibility is not verified.
%             Otherwise, Y is a floating-point vector which is for all
%             real input data (C,A,b) within the interval input data a dual
%             feasible solution (i.e. LMI certificate of feasibility), and
%        dl   is a n-vector, where dl(j)is a rigorous lower bound of
%             the smallest eigenvalue of C(j)-sum(i=1:m| Y(i)*A{j,i}).
%             dl > 0 implies the existence of strictly dual feasible solutions
%             and strong duality.
%
% [fL, Y, dl] = VSDPLOW(blk,A,C,b,Xt,yt,Zt,xu) returns a rigorous lower
%             bound of the primal minimum value for the above SDP-problem,
%             where upper bounds xu(j) (j = 1,...,n) for the
%             maximal eigenvalues of some optimal X(j) are known.
%             The upper bounds xu(j) also may be equal to infinity.
%
%
% We recommend to use infinite bounds xu(j) instead of unreasonable
% large bounds xu(j). This improves the quality of the lower bound
% in many cases, but may increase the computational time.
%
% EXAMPLE:
% C{1} = [1 0; 0 1];
% A{1,1} = [0 1; 1 0];
% A{2,1} = [1 1; 1 1];
% b = [1;2.0001];
% blk{1,1} = 's'; blk{1,2} = 2;
% [objt,Xt,yt,Zt,info] = mysdps(blk,A,C,b); % Computes approximations
% [fL, Y, dl] = vsdplow(blk,A,C,b,Xt,yt,Zt)
% fL =
%     1.0001
% Y =
%    -1.0000
%     1.0000
% dl =
%   9.3912e-011

% written   01/25/05   Christian Jansson
% modified  09/18/06
% modified  10/30/06
% Reference: C. JANSSON, Termination and Verification for
%            Ill-posed Semidefinite Programming Problems,
%            to appear
% http://optimization-online.org/DBHTML/|2005/06/1150.html.

% Calling the m-file of the global variables
SDP_GLOBALPARAMETER;

% Choice of the maximal number of iterations in VSDP
global VSDP_ITER_MAX;
% Choice of the perturbation parameter; DEFAULT = 2;
global VSDP_ALPHA ;

% Dimension, checks;
yt = yt(:);
b = b(:);
m = length(b);
n = length(C);

fL = -inf;
Y = NaN;
dl = NaN(n,1);

dim = size(A);
if ((dim(1) ~= m) || (dim(2) ~= n)) || (length(yt) ~= m)
  disp('SDP has wrong dimension.');
  return;
end

% NaN check
if max(isnan(yt)) || max(isinf(yt))
  disp('VSDPINFEAS: SDP-solver in MYSDPS computes NaN components.');
  return;
end

% Initialization
if nargin <= 7
  xu = inf * ones(n,1);
elseif isempty(xu)
  xu = inf * ones(n,1);
elseif length(xu) ~= n
  disp('VSDPLOW: Dimension of xu not compatible.');
  return;
elseif min(xu) < 0
  disp('VSDPLOW: primal upper bounds must be nonnegative.');
  return;
end
xu = xu(:);
% Possibly some further checks

% Preallocations
midA = cell(m,n);
midC = cell(1,n);
Ceps = cell(1,n);
D = cell(1,n);
Dlow = cell(1,n);
Dup = cell(1,n);
l = zeros(n,1);

iter = 0;
k = zeros(n,1);
epsj = zeros(n,1);

stop = 0;

% Intval input check
intvalinput = 0;
for j = 1 : n
  for i = 1 : m
    if isintval(A{i,j})
      intvalinput = 1;
      break;
    end
  end
  if (isintval(C{j})) || (intvalinput == 1)
    intvalinput = 1;
    break;
  end
end
if isintval(b)
  intvalinput = 1;
end


% Generating the midpoint problem using sparse format
if intvalinput == 1    % All data are transformed to sparse intval
  for j = 1 : n
    for i = 1 : m
      A{i,j} = sparse(A{i,j});
      if isintval(A{i,j})
        midA{i,j} = mid(A{i,j});
      else
        midA{i,j} = A{i,j};
        A{i,j} = intval(A{i,j});
      end
    end
    C{j} = sparse(C{j});
    if isintval(C{j})
      midC{j} = mid(C{j});
    else
      midC{j} = C{j};
      C{j} = intval(C{j});
    end
  end
  if isintval(b)
    midb = mid(b);
  else
    midb = b;
    b = intval(b);
  end
else % non-interval case: data transformed to sparse
  for j = 1 : n
    for i = 1 : m
      A{i,j} = sparse(A{i,j});
    end
    C{j} = sparse(C{j});
  end
end


% Algorithm
while (~stop) && (iter <= VSDP_ITER_MAX)
  % 1.step: efficient defect computation using monotonic
  %         roundings for real input
  if intvalinput == 1 % Using interval arithmetic for the D{j}
    for j = 1 : n
      D{j} = C{j};
      for i = 1 : m
        D{j} = D{j} - yt(i) * A{i,j};
      end
    end
  else      % Using monotonic roundings for the D{j} if noninterval input
    ytneg = -yt;
    for j = 1 : n
      setround(-1);
      Dlow{j} = C{j};
      for i = 1 : m
        Dlow{j} = Dlow{j} + ytneg(i) * A{i,j};
      end
      setround(1);
      Dup{j} = C{j};
      for i = 1 : m
        Dup{j} = Dup{j} + ytneg(i) * A{i,j};
      end
      D{j} = infsup(Dlow{j},Dup{j});
    end
  end
  setround(0);
  % 2.step
  for j = 1 : n
    Dj = D{j};
    lowbounds = inf_(veigsym(Dj));
    dl(j) = min(lowbounds);
    l(j) = sum(lowbounds < 0);
    % Implement alternatively another Eigenvalue solver
  end
  dlminus = min(0,dl);
  % 3.step
  if max((dl < 0)&(xu==inf)) == 0
    fL = inf_(b'*intval(yt));
    setround(-1);
    for j = 1 : n
      if xu(j) < inf
        fL = fL + l(j) * dlminus(j) * xu(j);
      end
    end
    setround(0);
    if max(dl < 0) == 0
      Y = yt;
    else
      Y = NaN;
      dl = NaN(n,1);
    end
    return;
  end
  % 4.step:  perturbed problem
  for j = 1 : n
    if (dl(j) < 0) && (xu(j) == inf)
      k(j) = k(j) +1;
      epsj(j) = -(VSDP_ALPHA^k(j)) * dl(j) + epsj(j);
    end
    if intvalinput == 1
      Ceps{j} = midC{j} - epsj(j) * speye(size(C{j},1));
    else
      Ceps{j} = C{j} - epsj(j) * speye(size(C{j},1));
    end
  end
  % 5.step: Call of the SDP-solver
  if intvalinput == 1
    [~,Xt,yt,Zt,info] = mysdps(blk,midA,Ceps,midb,Xt,yt,Zt);
  else
    [~,Xt,yt,Zt,info] = mysdps(blk,A,Ceps,b,Xt,yt,Zt);
  end
  
  % NaN check
  if max(isnan(yt)) || max(isinf(yt))
    disp('VSDPINFEAS: SDP-solver in MYSDPS computes NaN components.');
    return;
  end
  
  if ((info(1) == 2) || (info(1) == 3))
    stop = 1;
  end
  % 6.step
  iter = iter+1;
  dl = NaN(n,1);
end

end
