function [A, C] = sdpt3_to_vsdp(blk,At,Ct,b)
% SDPT3_TO_VSDP: Conversion routine from SDPT3 to VSDP format
%        for the block-diagonal sdp structure:
%        min  sum(j=1:n| <C{j}, X{j}>),
%        s.t. sum(j=1:n| <A{i,j}, X{j}> = b(i)) for i = 1 : m,
%             X{j} positive semidefinite for j = 1 : n,
%        A, C, b must be real (non-interval) quantities.
%
% The block-diagonal structure of both, VSDP and SDPT3, is
% described by an n*2 cell array named blk as follows:
% If the j-th block of C and the blocks A{i,j} for i = 1 : m
% are single real symmetric matrices of common size s_j, then
%         blk{j,1} = 's', blk{j,2} = [s_j]
% The difference between both formats is:
% In VSDP:
% A is a m*n  cell array containing the matrices A{i,j}.
% In SDPT3:
% The blocks A{i,j} are stored in n cells using the SVEC operation
%      At{j} = [svec(A{1,j) ... svec(A{m,j})];
% that is, each cell At{j} contains in succession the blocks
% corresponding to type j.
% Ct = C, which is a n-cell containing the block matrices C{j} for
% j = 1,...,n.
% This routine requires the SDPT3 file SMAT, which may cause
% conversion errors for the input data.

% Copyright 2004-2006 Christian Jansson (jansson@tuhh.de)

C = Ct;
m = length(b);
n = length(C);

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
