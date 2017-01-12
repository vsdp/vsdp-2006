function vA = vsvec(A,sparseflag,mult)
% VSVEC  Vectorize a symmetric block diagonal matrix.
%
%   Given a block diagonal matrix A with n blocks:
%
%           [A{1}         ]
%       A = [     ...     ]
%           [         A{n}]
%
%   stored as cell(n,1)-array.  Each of the j = 1:n blocks represents a
%   symmetric matrix of dimension blk(j).  VSVEC vectorizes the triangular
%   lower matrix of each of the n blocks in the following way:
%
%                                            [a     ]
%              [a b c]                       [b*mult]
%       A{j} = [b d e]   ==>   vsvec(A(j)) = [c*mult]
%              [c e f]                       [d     ]
%                                            [e*mult]
%                                            [f     ]
%
%   The off diagnoal elements are scaled with 'mult'.  The resulting column
%   vector of all n blocks has a length of sum(j=1:n |blk(j)*(blk(j)+1)/2).
%
%   vA = VSVEC(A) returns a full vector with off-diagonal elements scaled with
%      'mult = sqrt(2)'.
%
%   VSVEC(...,sparseflag) optionally decide whether to return a full vector
%      using 'sparseflag = 0' (default) or a sparse vector using
%      'sparseflag = 1'.
%
%   VSVEC(...,sparseflag,mult) optionally use 'sparseflag' as before and
%      specifiy another scaling factor for the off-diagonal elements.  The
%      default scaling factor is 'mult = sqrt(2)'.
%
%   The function VSVEC is useful for efficiently computing the inner product
%   of symmetric block diagonal matrices A, B, using the identity
%
%              <A,B> = VSVEC(A)'     * VSVEC(B)
%                    = VSVEC(A,0,2)' * VSVEC(B,0,1).
%
%   The latter identity is preferable, if one wishes to avoid rounding errors.
%
%   Example:
%
%      A = {[10  2  3;
%             2 11  4;
%             3  4 12];
%           ones(3)};
%      vA = vsvec(A,0,2);
%
%   See also vsmat.

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
