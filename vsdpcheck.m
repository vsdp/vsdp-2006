function [m, n] = vsdpcheck(blk,A,C,b,X0,y0,Z0)
% VSDPCHECK: routine for checking the VSDP format
%        of the block-diagonal sdp problem:
%        min  sum(j=1:n| <C{j}, X{j}>),
%        s.t. sum(j=1:n| <A{i,j}, X{j}> = b(i)) for i = 1 : m,
%             X{j} positive semidefinite for j = 1 : n,
%        A, C, b can be real  or interval quantities.
%
% The block-diagonal structure of VSDP is
% described by an n*2 cell array named blk as follows:
% If the j-th block of C and the blocks A{i,j} for i = 1 : m
% are single real symmetric matrices of common size s_j, then
%         blk{j,1} = 's', blk{j,2} = [s_j]
%
% A is a m*n  cell array containing the matrices A{i,j}.
% C is a n-cell containing the block matrices C{j} for
% j = 1,...,n.

% written   09/10/06   Christian Jansson


m = length(b);
n = length(C);
dim = size(A);
if ((dim(1) ~= m) || (dim(2) ~= n))
  disp('VSDP_CHECK: SDP has wrong dimension.');
  return;
end

% Conversion to noninterval data
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

if intvalinput == 1
  for j = 1 : n
    for i = 1 : m
      if isintval(A{i,j})
        A{i,j} = mid(A{i,j});
      end
    end
    if isintval(C{j})
      C{j} = mid(C{j});
    end
  end
  if isintval(b)
    b = mid(b);
  end
end

[At, Ct] = vsdp_to_sdpt3(blk,A,C,b);

if nargin <= 4
  X0 = Ct;
  Z0 = Ct;
  y0 = b;
end

check(blk,At,Ct,b,X0,y0,Z0);

end



function [A,C,X0,Z0] = check(blk,A,C,b,X0,y0,Z0)
spdensity = 0.5;

if ~iscell(blk)
  error('VSDP_CHECK: blk must be a cell array'); end
if size(blk,2)~=2
  error('VSDP_CHECK: blk must be a cell array with 2 columns');
end
if ~iscell(A) || ~iscell(C)
  error('VSDP_CHECK: A, C must be cell arrays'); end
if (min(size(A))~=1 || min(size(C))~=1)
  error('VSDP_CHECK: cell array A, C can only have 1 column or row'); end
if (size(A,2) > size(A,1))
  A = A'; end
if (size(C,2) > size(C,1))
  C = C'; end
if (nargin == 7)
  if ~iscell(X0) || ~iscell(Z0)
    error('VSDP_CHECK: X0, Z0 must be cell arrays');
  end
  if (min(size(X0))~=1 || min(size(Z0))~=1)
    error('VSDP_CHECK: cell array X, Z can only have 1 column or row');
  end
  if (size(X0,2) > size(X0,1)); X0 = X0'; end
  if (size(Z0,2) > size(Z0,1)); Z0 = Z0'; end
end


m = length(b);
for p=1:size(blk,1)
  pblk = blk(p,:);
  n = sum(pblk{2});
  numblk = length(pblk{2});
  if strcmp(pblk{1},'s')
    n2 = sum(pblk{2}.*pblk{2});  n22 = sum(pblk{2}.*(pblk{2}+1))/2;
    if ~all(size(C{p}) == n)
      error('VSDP_CHECK: blk and C are not compatible'); end
    if (norm(C{p}-C{p}',inf) > 1e-13)
      error('VSDP_CHECK: C is not symmetric'); end
    if (all(size(A{p}) == [m, n22]) && m~=n22)
      A{p} = A{p}'; end
    if ~all(size(A{p}) == [n22, m])
      error('VSDP_CHECK: blk and A not compatible'); end
    if (nnz(A{p}) < spdensity*n22*m)
      if ~issparse(A{p}); A{p} = sparse(A{p}); end
      %else
      %   if issparse(A{p}); A{p} = full(A{p}); end
    end
    if (nnz(C{p}) < spdensity*n2) || (numblk > 1)
      if ~issparse(C{p}); C{p} = sparse(C{p}); end
    else
      if issparse(C{p}); C{p} = full(C{p}); end
    end
    if (nargin == 7)
      if ~all(size(X0{p}) == n) || ~all(size(Z0{p}) == n)
        error('VSDP_CHECK: blk and X0,Z0 are not compatible'); end
      if (length(y0) ~= m)
        error('VSDP_CHECK: length of b and y0 not compatible'); end
      if (norm([X0{p}-X0{p}' Z0{p}-Z0{p}'],inf) > 2e-13)
        error('VSDP_CHECK: X0,Z0 not symmetric'); end
      if (nnz(X0{p}) < spdensity*n2) || (numblk > 1)
        if ~issparse(X0{p}); X0{p} = sparse(X0{p}); end
      else
        if issparse(X0{p}); X0{p} = full(X0{p}); end
      end
      if (nnz(Z0{p}) < spdensity*n2) || (numblk > 1)
        if ~issparse(Z0{p}); Z0{p} = sparse(Z0{p}); end
      else
        if issparse(Z0{p}); Z0{p} = full(Z0{p}); end
      end
    end
  elseif strcmp(pblk{1},'q') || strcmp(pblk{1},'l') || strcmp(pblk{1},'u')
    if (size(C{p},2) ~= 1)
      error(['VSDP_CHECK: ',num2str(p),'-th block of C must be column vectors']);
    end
    if (size(C{p},1) ~= n)
      error('VSDP_CHECK: blk and C are not compatible');
    end
    if (all(size(A{p}) == [m,n]) && m~=n)
      A{p} = A{p}'; end
    if ~all(size(A{p}) == [n,m])
      error('VSDP_CHECK: blk and A not compatible'); end
    if ~issparse(A{p})
      A{p} = sparse(A{p}); end
    if (nnz(C{p}) < spdensity*n)
      if ~issparse(C{p}); C{p} = sparse(C{p}); end
    else
      if issparse(C{p}); C{p} = full(C{p}); end
    end
    if (nargin == 7)
      if ~all([size(X0{p},2) size(Z0{p},2)]==1)
        error(['VSDP_CHECK: ',num2str(p),'-th block of X0,Z0 must be column vectors']);
      end
      if ~all([size(X0{p},1) size(Z0{p},1)]==n)
        error('VSDP_CHECK: blk, and X0,Z0, are not compatible');
      end
      if (nnz(X0{p}) < spdensity*n)
        if ~issparse(X0{p}); X0{p} = sparse(X0{p}); end
      else
        if issparse(X0{p}); X0{p} = full(X0{p}); end
      end
      if (nnz(Z0{p}) < spdensity*n)
        if ~issparse(Z0{p}); Z0{p} = sparse(Z0{p}); end
      else
        if issparse(Z0{p}); Z0{p} = full(Z0{p}); end
      end
      if strcmp(pblk{1},'u')
        Z0{p} = sparse(n,1);
      end
    end
  else
    error(' blk: some fields are not specified correctly');
  end
end

end
