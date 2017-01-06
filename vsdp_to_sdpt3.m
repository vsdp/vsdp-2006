function [At, Ct] = vsdp_to_sdpt3(blk,A,C,b)
% VSDP_TO_SDPT3: Conversion routine from VSDP to SDPT3 format
%         for the block-diagonal sdp structure:
%         min  sum(j=1:n| <C{j}, X{j}>),
%         s.t. sum(j=1:n| <A{i,j}, X{j}> = b(i)) for i = 1 : m,
%              X{j} positive semidefinite for j = 1 : n,
%         A, C, b must be real (non-interval) quantities.
%
% The block-diagonal structure of VSDP  is
% described by an n*2 cell array named blk as follows:
% If the j-th block of C and the blocks A{i,j} for i = 1 : m
% are single real symmetric matrices of common size s_j, then
%          blk{j,1} = 's', blk{j,2} = [s_j]
% The difference between both formats is:
% In VSDP:
% A is a n*m cell array containing the matrices A{i,j}.
% In SDPT3:
% The blocks A{i,j} are stored in n cells using the SVEC operation
%      At{j} = [svec(A{1,j) ... svec(A{m,j})];
% that is, each cell At{j} contains in succession the blocks
% corresponding to type j.
% Ct = C, which is a n-cell containing the block matrices C{j} for
%j = 1,...,n.
%This routine requires the SDPT3 file SVEC, which may cause
%conversion errors for the input data.
%It is written such that the input C, A may be also matrices.


% written   07/25/05   Christian Jansson
% modified  09/10/06

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
