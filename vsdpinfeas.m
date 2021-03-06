function [isinfeas,X,Y] = vsdpinfeas(blk,A,C,b,choose,Xt,yt,~)
% VSDPINFEAS  Infeasibility-check for the block-diagonal problem:
%
%         min  sum(j=1:n| <  C{j}, X{j}>)
%         s.t. sum(j=1:n| <A{i,j}, X{j}>) = b(i)  for i = 1:m
%              X{j} must be positive semidefinite for j = 1:n
%
%   as well as its dual:
%
%         max  b'*y
%         s.t. C{j} - sum(i=1:m| y{i}*A{i,j}) must be positive semidefinite
%                                             for j = 1:n.
%
%   [isinfeas,X,Y] = VSDPINFEAS(blk,A,C,b,choose)
%      The block-diagonal format (blk,A,C,b) is explained in 'mysdps.m'.
%
%         'choose'    If the character is 'p', primal infeasibility should be
%                     verified.  If 'd', dual infeasibility should be verified.
%
%      The output is:
%
%         'isinfeas'  Returns 1 if the primal or dual problem is proved to
%                     be infeasible and 0 if infeasibility cannot be verified.
%
%         'X'         Contains a rigorous certificate (improving ray) of dual
%                     infeasibility, if it is not equal to NaN.
%
%         'Y'         Is a rigorous certificate (improving ray) of primal
%                     infeasibility, if it is not equal to NaN.
%
%   VSDPINFEAS(...,Xt,yt,Zt) optionally provide already by 'mysdps' computed
%      approximate solutions (Xt,yt,Zt).  This avoids calling 'mysdps' from
%      within this function, if approximate solutions are already present.
%
%   Example:
%
%       EPS = -0.01;
%       DELTA = 0.1;
%       blk(1,:) = {'s'; 2};
%       A{1,1} = [1 0; 0 0];
%       A{2,1} = [0 1; 1 DELTA];
%         C{1} = [0 0; 0 0];
%            b = [EPS; 1];
%       [objt,Xt,yt,Zt,info] = mysdps(blk,A,C,b);
%       choose = 'p';
%       [isinfeas,X,Y] = vsdpinfeas(blk,A,C,b,choose,Xt,yt,Zt);
%
%   See also mysdps.

% Copyright 2004-2006 Christian Jansson (jansson@tuhh.de)

% Setting default output
isinfeas = 0;
X = NaN;
Y = NaN;

% Dimension, checks, initialization
b = b(:);
m = length(b);
n = length(C);
dim = size(A);
if ((dim(1) ~= m) || (dim(2) ~= n))
  disp('VSDPINFEAS: SDP has wrong dimension.');
  return;
end

% Basic and nonbasic Indizes
I = [];
N = [];

% Preallocations
Amid = cell(m,n);
Ai = cell(m,n);
Cmid = cell(1,n);
D = cell(1,n);
Dlow = cell(1,n);
Dup = cell(1,n);
dup = zeros(n,1);
blkvec = zeros(n,1);
lb = zeros(n,1);

if choose == 'p'
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
  
  if nargin <= 5
    if intvalinput == 0
      % Solving phase I problem approximately
      % A stable phase I problem should be incorporated !!!!!!
      [~,Xt,yt,~,~] = mysdps(blk,A,C,b);
    else
      % Transformation to Phase I midpoint-problem
      bmid = mid(b);
      for j = 1 : n
        for i = 1 : m
          Amid{i,j} = mid(A{i,j});
        end
        Cmid{j} = mid(C{j});
      end
      % Solving phase I problem approximately
      [~,Xt,yt,~,~] = mysdps(blk,Amid,Cmid,bmid);
    end
  end
  
  if max(isnan(yt)) || max(isinf(yt))
    disp('VSDPINFEAS: SDP-solver in MYSDPS computes NaN components.');
    return;
  end
  
  % Verification of a primal improving ray using approximation yt
  if intvalinput == 1      % then interval arithmetic
    % Compuation of the LMI's
    for j = 1 : n
      D{j} = intval(yt(1)) * A{1,j};
      for i = 2 : m
        D{j} = D{j} + intval(yt(i)) * A{i,j};
      end
    end
  else        % then monotonic roundings for point data
    for j = 1 : n
      setround(-1);
      Dlow{j} = yt(1) * A{1,j};
      for i = 2 : m
        Dlow{j} = Dlow{j} + yt(i) * A{i,j};
      end
      setround(1);
      Dup{j} = yt(1) * A{1,j};
      for i = 2 : m
        Dup{j} = Dup{j} + yt(i) * A{i,j};
      end
      D{j} = infsup(Dlow{j},Dup{j});
    end
  end
  setround(0);
  % The maximal eigenvalue of D{j}
  for j = 1 : n
    upbounds = sup(veigsym(D{j}));
    dup(j) = max(upbounds);
  end
  % Check certificate
  if ( all(dup <= 0) ) && ( (intval(b)'*yt) > 0 )
    isinfeas = 1;
    X = NaN;
    Y = yt;
    return;
  end
end


if choose == 'd'
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
  
  if nargin <= 5
    if intvalinput == 0
      % Solving phase I problem approximately
      % A stable phase I problem should be incorporated !!!!!!!!!
      [~,Xt,~,~,~] = mysdps(blk,A,C,b);
    else
      % Transformation to Phase I midpoint-problem
      bmid = mid(b);
      for j = 1 : n
        Cmid{j} = mid(C{j});
        for i = 1 : m
          Amid{i,j} = mid(A{i,j});
        end
      end
      % Solving phase I problem approximately
      [~,Xt,~,~,~] = mysdps(blk,Amid,Cmid,bmid);
    end
  end
  
  
  % Transformations to the m * (n*(n+1)/2) linear system
  % using sparse format
  if intvalinput == 0      % Case of non-interval input data
    vC = sparse(vsvec(C,0,2));
    nls = length(vC);
    vXt = sparse(vsvec(Xt,0,1));
    if max(isnan(vXt))
      disp('VSDPINFEAS: SDP-solver in MYSDPS computes NaN components.');
      return;
    end
    vX = vXt;
    Xbounds = infsup(-inf*ones(nls,1), inf*ones(nls,1));
    b = zeros(m+1,1);
    b(m+1) = vC' * vXt;
    if b(m+1) >= 0
      isinfeas = 0;
      X = NaN;
      return;
    end
    Amat = speye(m+1, nls);
    for i = 1:m
      for j = 1:n
        Ai{j} = A{i,j};
        blkvec(j) = blk{j,2};
      end
      Amat(i,:) = sparse(vsvec(Ai,0,2));
    end
    Amat(m+1,:) = vC;
    % Verified Solution of the linear system
    [vX,~,~,~] = vuls([], [], Amat, b, inf_(Xbounds), sup(Xbounds),...
      vX,I,N);
    if isnan(vX)
      disp('VSDINFEAS: system matrix may have no full rank');
      return;
    end
    Xwork = vsmat(vX,blkvec,0,1);
    for j = 1 : n
      Xj = Xwork{j};
      lowbounds = inf_(veigsym(Xj));
      lb(j) = min(lowbounds);
    end
    % Check certificate
    if min(lb) >= 0
      isinfeas = 1;
      X = Xwork;
      return;
    end
  else    % Case of interval input data
    vC = intval(sparse(vsvec(C,0,2)));
    nls = length(vC);
    vXt = sparse(vsvec(Xt,0,1));
    if max(isnan(vXt))
      disp('VSDPINFEAS: SDP-solver in MYSDP has NaN components.');
      return;
    end
    vX = vXt;
    Xbounds = infsup(-inf*ones(nls,1), inf*ones(nls,1));
    b = intval(zeros(m+1,1));
    b(m+1) = intval(inf_(vC' * vXt));
    if b(m+1) >= 0
      isinfeas = 0;
      X = NaN;
      return;
    end
    Amat = intval(speye(m+1, nls));
    for i = 1:m
      for j = 1:n
        Ai{j} = intval(A{i,j});
        blkvec(j) = blk{j,2};
      end
      Amat(i,:) = intval(sparse(vsvec(Ai,0,2)));
    end
    Amat(m+1,:) = vC;
    % Verified Solution of the linear system
    [vX,~,~,~] = vuls([], [], Amat, b, inf_(Xbounds), sup(Xbounds),...
      vX,I,N);
    if isnan(vX)
      disp('VSDINFEAS: system matrix may have no full rank');
      return;
    end
    Xwork = vsmat(vX,blkvec,0,1);
    for j = 1 : n
      Xj = Xwork{j};
      lowbounds = inf_(veigsym(Xj));
      lb(j) = min(lowbounds);
    end
    % Check certificate
    if min(lb) >= 0
      isinfeas = 1;
      X = Xwork;
      return;
    end
  end
end

end
