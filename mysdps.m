function [obj,X,y,Z,info] = mysdps(blk,A,C,b,X0,y0,Z0)
% MYSDPS  My Semidefinite Programming solver
%         for the block-diagonal problem
%
%         min  sum(j=1:n| <C{j}, X{j}>)
%         s.t. sum(j=1:n| <A{i,j}, X{j}>) = b(i) for i = 1 : m
%              X{j} must be positive semidefinite for j = 1 : n
%
%         Moreover, a certificate of feasibility for LMI's is provided.
%         The routine uses the semidefinite solver SDPT3.
%
% The block-diagonal structure is described
% by an n*2-cell-array blk, n-cell arrays C, X, and an
% m*n-cell-array A as follows:
% The j-th block C{j} and the blocks A{i,j} for i = 1 : m
% are real symmetric matrices of common size s_j, and
%    blk{j,1} = 's', blk{j,2} = [s_j]
% The blocks C{j} and A{i,j} must be stored as individual
% matrices in dense or sparse format.
%
%Output: obj  = [<C,X> <b,y>].
%        (X0,y0,Z0): an approximately optimal solution or a primal or dual
%                     infeasibility certificate.
%        info = termination-code with
%        info = 0 : indication of optimality (normal termination),
%        info = 1 : indication of primal infeasibility,
%        info = 2 : indication of dual infeasibility,
%        info = 3 : SDPT3: indication of both primal and dual infeasibilities,
%                   SDPA-M: indication of primal or dual infeasibilities,
%        info = -1 : otherwise

% Christian Jansson      Last Update: December 2005

SDP_GLOBALPARAMETER;
%SDP_GLOBALPARAMETER; calls the desired SDP-solver in MYSDPS
%VSDP_CHOICE_SDP = 1: SDPT3 solver
%VSDP_CHOICE_SDP = 2: SDPA solver
global VSDP_CHOICE_SDP;
global VSDP_USE_STARTING_POINT;

if VSDP_CHOICE_SDP == 1
  %Choose is SDPT3 Solver
  OPTIONS.vers           = 1;
  OPTIONS.gam            = 0;
  OPTIONS.predcorr       = 1;
  OPTIONS.expon          = [1 1];
  OPTIONS.gaptol         = 1e-8;    %1e-6;
  OPTIONS.inftol         = 1e-8;    %1e-6;
  OPTIONS.steptol        = 1e-6;
  OPTIONS.maxit          = 50;
  OPTIONS.printyes       = 1;
  OPTIONS.scale_data     = 0;
  OPTIONS.randnstate     = 0;
  OPTIONS.spdensity      = 0.5;
  OPTIONS.rmdepconstr    = 0;
  OPTIONS.cachesize      = 256;
  OPTIONS.smallblkdim    = 15;
  A = A';
  
  if (nargin <= 5) || (VSDP_USE_STARTING_POINT == 0)
    [obj,X,y,Z,info] = sqlp(blk,A,C,b,OPTIONS);
  else
    [obj,X,y,Z,info] = sqlp(blk,A,C,b,OPTIONS, X0,y0,Z0);
  end
  
  % for SDPT3-4.0
  if (isstruct(info) && isfield(info, 'termcode'))
    info = info.termcode; 
  end
  
elseif VSDP_CHOICE_SDP == 2
  %Choose is SDPA Solver
  [obj,X,y,Z,info] = sqlp_buf(blk,A',C,b);
  
end

end



function [obj,X,y,Z,infos] = sqlp_buf(blk,Avec,C,b)
mDIM = size(b,1);
nBLOCK = max(size(C,1),size(C,2));
if size(C,1) < size(C,2)
  C=C';
end
for i = 1 : nBLOCK
  buf = blk(i,2);     %check in get.sdpa.m
  bLOCKsTRUCT(i) = buf{1,1};
end
clear buf;
F(:,1) = C(1:nBLOCK,1);
F(:,2:mDIM+1) = Avec;
for i = 1:size(F,1)
  for j = 1:size(F,2)
    F{i,j} = -F{i,j};
  end
end

OPTION.maxIteration   = 40;
OPTION.epsilonStar    = 1.0e-7; %1.0e-07;
OPTION.lambdaStar     = 100;
OPTION.omegaStar      = 2;
OPTION.lowerBound     = -100000;
OPTION.upperBound     = 100000;
OPTION.betaStar       = 0.1;
OPTION.betaBar        = 0.2;
OPTION.gammaStar      = 0.9;
OPTION.epsilonDash    = 1.0e-7; %1.0e-07;
OPTION.isSymmetric    = 1;
%OPTION.print          = 'display';
OPTION.print          = 'no';

[obj,y,Z,X,INFO]=sdpam(mDIM,nBLOCK,bLOCKsTRUCT,-b,F,OPTION);
INFO.phasevalue,
buf = obj(1);
obj(1) = -obj(2); obj(2) = -buf;
switch INFO.phasevalue
  case 'pdOPT'
    infos(1) = 0;
  case 'pdINF'
    infos(1) = 3;
  case 'pFEAS_dINF'
    infos(1) = 1;
  case 'pINF_dFEAS'
    infos(1) = 2;
  otherwise
    infos(1) = -1;
end

end
