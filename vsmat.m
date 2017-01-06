function A = vsmat(vA,blk,sparseflag,mult)
%A = VSMAT(vA,blk,sparseflag,mult) is the inverse operation of
%    vA = VSVEC(A,sparseflag,1/mult), and returns the
%    blocks of the diagonal matrix A = diag(A{1}, ..., A{l}).
%
%The blocks A{j} are real or interval matrices of size blk(j) for j = 1 : l,
%and they are stored as a cell array.
%If sparseflag = 0, VSMAT returns full blocks and sparse blocks otherwise.
%The default is sparseflag = 0. Default for the multiplier is mult = 1/sqrt(2).
%
%The concatenated vector vA is of length sum{blk(j)*(blk(j)+1)/2: j=1, .. ,l},
%and the blocks A{j} are stored in vA as follows:
%        vA = [A{1}(1,1),mult*A{1}(2,1), ... ,mult*A{1}(blk(1),1)
%              A{1}(2,2),mult*A{1}(3,2), ... ,mult*A{1}(blk(1),2)
%              A{1}(3,3), ... , A{1}(blk(1),blk(1)), A{2}(1,1), ...].
%In the case of mult = sqrt(2), the inner product of symmetric matrices A, B
%satisfies the identity
%              <A,B> = VSVEC(A,blk)' * VSVEC(B,blk)
%If one wishes to avoid rounding errors, then the use of
%              <A,B> = VSVEC(A,blk,0,2)' * VSVEC(B,blk,0,1)
%is preferable.
%
%EXAMPLE:
%A1 = ones(3,3); A1(3,1) = 3; A1(1,3) = 3;
%A{1} = A1; blk(1) = 3;
%vA = vsvec(A,0,sqrt(2));
%AA = vsmat(vA,blk,0,1/sqrt(2));
%AA{1}
%ans =
%     1     1     3
%     1     1     1
%     3     1     1
%
%Author: Christian Jansson      Last Update: June 04, Mar 05

if nargin < 4
  mult = 1/sqrt(2);
  if nargin < 3
    sparseflag = 0;
  end
end
l = length(blk);
dim = size(vA);
if dim(2) > dim(1),
  vA = vA';
end;

start = 1;
ende = 0;
for j = 1:l,
  blocksize = blk(j);
  blocklength = blocksize*(blocksize+1)/2;
  ende = ende + blocklength;
  %The j-th block
  if isintval(vA)
    block = intval(zeros(blocksize));
  else
    block = zeros(blocksize);
  end
  Index =  repmat((1:blocksize),blocksize,1);
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
  end;
end;

%________________________________END OF VSMAT___________________________

