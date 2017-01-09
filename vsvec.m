function vA = vsvec(A,sparseflag,mult)
%VSVEC:  returns the concatenated vector vA which corresponds to the
%        block diagonal matrix A = diag(A{1}, ..., A{l}).
%
%The blocks A{i} may be a real or interval matrices, and must be stored
%as a cell array. If sparseflag = 0, VSVEC returns a full vector and
%a sparse vector otherwise. The default is sparseflag = 0.
%Default for the real multiplier is mult = sqrt(2).
%
%The concatenated vector vA of length sum{blk(j)*(blk(j)+1)/2 : j=1, .. ,l},
%where blk(j) denotes the size of the j-th bock, is defined as
%        vA = [A{1}(1,1),mult*A{1}(2,1), ... ,mult*A{1}(blk(1),1)
%              A{1}(2,2),mult*A{1}(3,2), ... ,mult*A{1}(blk(1),2)
%              A{1}(3,3), ... , A{1}(blk(1),blk(1)), A{2}(1,1), ...],
%and satisfies in the case of mult = sqrt(2) for the
%inner product of symmetric matrices A, B the identity
%              <A,B> = VSVEC(A,blk)' * VSVEC(B,blk)
%If one wishes to avoid rounding errors then the use of
%              <A,B> = VSVEC(A,blk,0,2)' * VSVEC(B,blk,0,1)
%is preferable.
%
%EXAMPLE:
%A{1} = ones(3,3); vA = vsvec(A); vA'
%ans =
%    1.0000    1.4142    1.4142    1.0000    1.4142    1.0000
%

% Copyright 2004-2006 Christian Jansson (jansson@tuhh.de)

if nargin < 3
  mult = sqrt(2);
  if nargin < 2
    sparseflag = 0;
  end
end

l = length(A);
vA = zeros(l,1);
vAstart = 1;
vAende = 0;

for j = 1 : l
  %Multiplication of lower and upper part
  Aj = A{j};
  Ajl= mult * tril(Aj,-1);
  Aju = Ajl';
  Aj = Ajl + Aju + diag(diag(Aj));
  if isintval(Aj)
    vA = intval(vA);
  end
  blocksize = size(Aj,1);
  blocklength = blocksize*(blocksize+1)/2;
  vAende = vAende + blocklength;
  Index =  repmat((1:blocksize),blocksize,1);
  Jndex = Index';
  vA(vAstart:vAende,1) = Aj(Index<=Jndex);
  %     vA(vAstart:vAende,1) = Aj(find(tril(ones(blocksize))));
  vAstart = vAstart + blocklength;
end

if sparseflag ~= 0
  vA = sparse(vA);
else
  vA = full(vA);
end

end
