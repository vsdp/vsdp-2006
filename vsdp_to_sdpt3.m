function [At,Ct] = vsdp_to_sdpt3(blk,A,C,b)
% VSDP_TO_SDPT3  Convert VSDP to SDPT3 format.
%
%   [At,Ct] = vsdp_to_sdpt3(blk,A,C,b)
%      The block-diagonal problem is:
%
%         min  sum(j=1:n| <  C{j}, X{j}>)
%         s.t. sum(j=1:n| <A{i,j}, X{j}>) = b(i)  for i = 1:m
%              X{j} must be positive semidefinite for j = 1:n
%
%      VSDP and SDPT3 share the storage format for 'blk', 'C', and 'b', see
%      'mysdps.m' for details.
%
%      The constraint matrices of the block-diagonal problem 'A' are stored in
%      VSDP as 'cell(m,n)' containing the matrices in the cells A{i,j}.
%      In SDPT3 the blocks A{i,j} are stored in n cells using the 'svec'
%      operation from SDPT3:
%
%         At{j} = [svec(A{1,j) ... svec(A{m,j})];
%
%      That is, each cell At{j} contains the m column vectorized contraint
%      matrices for block j.
%
%   This routine requires the SDPT3 file 'svec'.
%
%   Example:
%
%       blk(1,:) = {'s'; 2};
%       A{1,1} = [0 1; 1 0];
%       A{2,1} = [1 1; 1 1];
%         C{1} = [1 0; 0 1];
%            b = [1; 2.0001];
%       [At,Ct] = vsdp_to_sdpt3(blk,A,C,b);
%
%   See also mysdps, sdpt3_to_vsdp.

% Copyright 2004-2006 Christian Jansson (jansson@tuhh.de)


% Transformation to cells, and using the transposed
At = A'; Ct = C;
if ~iscell(A); At = {A}; end
if ~iscell(C);  Ct = {C}; end

%Problem size
m = length(b);
n = size(blk,1);
%  m,n ,size(At)
%Main routine
if all(size(At) == [n, m])
  transyes = zeros(n,1);
  for j = 1:n
    if strcmp(blk{j,1},'s') && all(size(At{j,1}) == sum(blk{j,2}))
      %sum(blk{j,2}) because of the special format for
      %defining many small subblocks as one block
      transyes(j) = 1;
    end
  end
  if any(transyes)
    %fprintf(' VSDP: converting At into SDPT3 format...\n');
    At = svec(blk,At);   %C - routine of SDPT3
  end
end

end
