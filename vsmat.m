function A = vsmat(vA,blk,sparseflag,mult)
% VSMAT  Inverse operation of 'vsvec'.
%
%   For a comprehensive explanation, see 'vsvec.m'.
%
%   A = VSMAT(vA,blk) is the inverse operation of 'vsvec', such that
%
%         isequal(A,vsmat(vsvec(A),blk))
%
%      evaluates true, given the correct block dimensions 'blk' of 'A'.  The
%      returned matrix is a block diagonal matrix with n blocks.
%
%      'vA'   A full or sparse vector of lengh sum(j=1:n |blk(j)*(blk(j)+1)/2),
%             for example created by 'vsvec(A)'.
%
%      'blk'  An n-vector where each element blk(j) represents the block size
%             of the resulting block diagonal matrix block A{j}.
%
%   VSMAT(...,sparseflag) optionally decide whether to return only full matrix
%      blocks using 'sparseflag = 0' (default) or only sparse matrix blocks
%      using 'sparseflag = 1'.
%
%   VSMAT(...,sparseflag,mult) optionally use 'sparseflag' as before and
%      specifiy another scaling factor for the off-diagonal elements of the
%      resulting block-diagonal matrix.  The default scaling factor is
%      1/sqrt(2).
%
%   Example:
%
%       A = {[10  2  3;
%             2 11  4;
%             3  4 12];
%           ones(3)};
%      vA = vsvec(A,0,2);
%     blk = [3; 3];
%      AA = vsmat(vA,blk,0,1/2);
%
%   See also vsvec.

% Copyright 2004-2006 Christian Jansson (jansson@tuhh.de)

if nargin < 4
  mult = 1/sqrt(2);
  if nargin < 3
    sparseflag = 0;
  end
end
l = length(blk);
dim = size(vA);
if dim(2) > dim(1)
  vA = vA';
end

start = 1;
ende = 0;
for j = 1:l
  blocksize = blk(j);
  blocklength = blocksize*(blocksize+1)/2;
  ende = ende + blocklength;
  %The j-th block
  if isintval(vA)
    block = intval(zeros(blocksize));
  else
    block = zeros(blocksize);
  end
  Index = repmat((1:blocksize),blocksize,1);
  Jndex = Index';
  block(Index<=Jndex) = vA(start:ende);
  start = start + blocklength;
  %Multiplication of lower and upper part,
  blocklow= mult * tril(block,-1);
  blockup = blocklow';
  block = blocklow + blockup + diag(diag(block));
  
  if sparseflag ~= 0
    A{j} = sparse(block);
  else
    A{j} = block;
  end
end

end
