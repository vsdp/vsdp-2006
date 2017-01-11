function [A,C] = sdpt3_to_vsdp(blk,At,Ct,b)
% SDPT3_TO_VSDP  Convert SDPT3 to VSDP format.
%
%   [A,C] = sdpt3_to_vsdp(blk,At,Ct,b)
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
%   This routine requires the SDPT3 file 'smat'.
%
%   Example:
%
%       blk(1,:) = {'s'; 2};
%       At{1} = [   0,       1;
%                sqrt(2), sqrt(2);
%                   0,       1];
%       Ct{1} = [1 0; 0 1];
%       b = [1; 2.0001];
%       [A,C] = sdpt3_to_vsdp(blk,At,Ct,b);
%
%   See also mysdps, vsdp_to_sdpt3.

% Copyright 2004-2006 Christian Jansson (jansson@tuhh.de)

C = Ct;
m = length(b);
n = length(C);

A = cell(m,n);

for j = 1 : n
  Atj = At{j};
  %full(Atj)
  blkj{1,1} = 's';
  blkj{1,2} = blk{j,2};
  for i = 1 : m
    A{i,j} = smat(blkj,Atj(:,i),ones(n,1));
    %full(A{j,i})
  end
end

end
