function [fU, X, lb] = vsdpup(blk,A,C,b,Xt,yt,Zt,yu)
% VSDPUP  Verified Semidefinite Programming Upper Bound   
%         for the minimum value of the block-diagonal problem
%
%         min  sum(j=1:n| <C{j}, X{j}>)
%         s.t. sum(j=1:n| <A{i,j}, X{j}>) = b(i) for i = 1 : m
%              X{j} must be positive semidefinite for j = 1 : n
%
%         Moreover, a rigorous certificate of primal feasibility 
%         is provided.
%
% The block-diagonal structure is described
% by an n*2-cell-array blk, n-cell-arrays C, X, and an 
% m*n-cell-array A as follows: 
% The j-th block C{j} and the blocks A{i,j} for i = 1 : m
% are real symmetric matrices of common size s_j, and
%    blk{j,1} = 's', blk{j,2} = s_j
% The blocks C{j} and A{i,j} must be stored as individual 
% matrices in dense or sparse format.
%
% The input A, C and the m-vector b 
% may be floating-point or interval quantities. 
%
% Xt,yt,Zt    approximate floating-point solution or initial 
%             starting point for the primal (this is Xt) and 
%             the dual problem (this is yt, Zt). 
%
%
% [fU, X, lb] = VSDPUP(blk,A,C,b,Xt,yt,Zt) returns 
%       fU    a rigorous upper bound of the minimum value
%             for all real input data (C,A,b) within the 
%             interval input data. 
%             fU = inf, if no finite rigorous upper bound
%             can be computed.
%        X    =NAN and lb=repmat(NaN,n,1), if primal feasibility is not verified. 
%             Otherwise, X is an interval quantity which contains a primal
%             feasible solution (certificate) of the block-diagonal
%             problem for all real input data within the interval data, and           
%        lb   is a n-vector, where lb(j) is a rigorous lower bound of 
%             the smallest eigenvalue of block X{j}. lb > 0 means 
%             that all symmetric matrices within the interval
%             matrix X are rigorously certified as positiv definite.
%             In this case, the existence of strictly feasible solutions
%             and strong duality is proved.
%             
% [fU, X, lb] = VSDPUP(blk,A,C,b,Xt,yt,Zt, yu) returns a rigorous 
%             upper bound fU of the dual optimal value, where
%             the following dual boundedness assumption is assumed:
%             an optimal dual solution y satisfies
%                     -yu <= y <= yu,  yu nonnegative m-vector.
%
% We recommend to use infinite bounds yu instead of unreasonable 
% large bounds yu. This improves the quality of the lower bound 
% in many cases, but may increase the computational time.
%          
% EXAMPLE:
% C{1} = [1 0; 0 1];
% A{1,1} = [0 1; 1 0]; 
% A{2,1} = [1 1; 1 1];
% b = [1;2.0001];
% blk{1,1} = 's'; blk{1,2} = 2;
% [objt,Xt,yt,Zt,info] = mysdps(blk,A,C,b); % Computes approximations
% [fU, X, lb] = vsdpup(blk,A,C,b,Xt,yt,Zt);
% intvalinit('displaymidrad'); format short e;
% fU, X{1},lb
% fU =
%   1.0001e+000
% intval ans = 
% <  5.0005e-001, 2.1337e-010> <  5.0000e-001, 0.0000e+000> 
% <  5.0000e-001, 0.0000e+000> <  5.0005e-001, 2.1337e-010> 
% lb =
%   5.0000e-005
% 

% written   04/21/05   Christian Jansson      
% modified  04/11/06   
% Reference: C. JANSSON, Termination and Verification for 
%            Ill-posed Semidefinite Programming Problems,
%            to appear
% http://optimization-online.org/DBHTML/|2005/06/1150.html. 

% Calling the m-file of the global variables
SDP_GLOBALPARAMETER;

% Choice of the maximal number of iterations in VSDP:
global VSDP_ITER_MAX;
% Choice of the perturbation parameter; DEFAULT = 2;
global VSDP_ALPHA ;

                                           
% Dimension, checks 
b = b(:);                                    
m = length(b);                                           
n = length(C);

fU = inf;                                           
X = NaN;
lb = repmat(NaN,n,1);

dim = size(A);
if (dim(1) ~= m) | (dim(2) ~= n)
   disp('VSDPUP: SDP has wrong dimension.'); 
   return; 
end;

% Initialization
if nargin <= 7   
   yu = inf * ones(m,1); 
elseif isempty(yu) 
   yu = inf * ones(m,1);   
elseif length(yu) ~= m
   disp('VSDPUP: Dimension of yu not compatible.'); 
   return;
elseif min(yu) < 0
   disp('VSDUP: dual upper bounds must be nonnegative.'); 
   return;    
end;
yu = yu(:);

iter = 0;
k = zeros(n,1);
rho = zeros(n,1);

% Preallocations
midA = cell(m,n);
Ai = cell(1,n);
midC = cell(1,n);
Ceps = cell(1,n);
Eps = cell(1,n);
II = cell(1,n);
U = cell(1,n);
D = cell(1,n);
Dlow = cell(1,n);
Dup = cell(1,n);
l = zeros(n,1);
blkvec = zeros(n,1);

% Basic and nonbasic Indizes
I = [];
N = [];

stop = 0;

% Intval input check
intvalinput = 0;   
for j = 1 : n 
   for i = 1 : m
      if isintval(A{i,j})
         intvalinput = 1; 
         break;
      end   
   end
   if (isintval(C{j})) | (intvalinput == 1)
         intvalinput = 1; 
         break;
   end
end
if isintval(b)
   intvalinput = 1; 
end  

%Transformations to the m * (n*(n+1)/2) linear system using sparse format
if intvalinput == 0      % Case of non-interval input data
   vC = sparse(vsvec(C,0,2));
   nls = length(vC);   
   vXt = sparse(vsvec(Xt,0,1));                           
   vX = vXt;     
   Xbounds = infsup(-inf*ones(nls,1), inf*ones(nls,1));
   Amat = speye(m, nls);
   for i = 1:m 
      for j = 1:n
         Ai{j} = A{i,j};
      end   
      Amat(i,:) = sparse(vsvec(Ai,0,2));
   end 
   setround(1);
   for j = 1 : n
      blkvec(j) = blk{j,2};
      II{j} = speye(blkvec(j));
      Eps{j} = sparse(zeros(blkvec(j)));
      % Computation of rho can be avoided if 
      % Xt is shifted appropriately
      if ~isempty(yu) & (max(yu) < inf)  
         U{j} = abs(C{j});
         for i = 1 : m
             U{j} = U{j} + yu(i) * abs(A{i,j});
         end
         rho(j) =  norm(U{j},inf);
      end 
   end 
   setround(0);
else    % Case of interval input data
   vC = intval(sparse(vsvec(C,0,2)));
   nls = length(vC);   
   vXt = sparse(vsvec(Xt,0,1));                           
   vX = vXt;                                                      
   Xbounds = infsup(-inf*ones(nls,1), inf*ones(nls,1));
   Amat = intval(speye(m, nls));
   for i = 1:m 
      for j = 1:n
         if isintval(A{i,j})
            midA{i,j} = mid(A{i,j});
         else
            midA{i,j} = A{i,j};
            A{i,j} = intval(A{i,j});% A should be intval for all i,j
         end
         Ai{j} = A{i,j};
      end   
      Amat(i,:) = intval(sparse(vsvec(Ai,0,2)));
   end 
   midAmat = sparse(mid(Amat));
   if isintval(b)
      midb = mid(b);
   else 
      midb = b;
      b = intval(b);
   end   
   for j = 1:n
      if isintval(C{j})
         midC{j} = mid(C{j});
      else
         midC{j} = C{j};
         C{j} = intval(C{j});% C should be intval for all j
      end             
      blkvec(j) = blk{j,2};
      II{j} = speye(blkvec(j));
      Eps{j} = sparse(zeros(blkvec(j)));
      setround(1);
      if ~isempty(yu) & (max(yu) < inf)   % Computation of rho
         U{j} = mag(C{j});
         for i = 1 : m
                U{j} = U{j} + yu(i) * mag(A{i,j});
         end
         rho(j) = norm(U{j},inf);     
      end 
      setround(0);
   end
end   
                                                     
if intvalinput == 1
   bperturb = midb;
else
   bperturb = b;
end   
vEps = zeros(length(vXt),1);                          
                                                  
% Algorithm with finite dual bounds yu
if max(yu) < inf
   if intvalinput == 0
      setround(1);
      rup = mid(b) + rad(b) + (- Amat) * vXt;
      setround(-1);
      rlow = mid(b) + (- rad(b)) + (- Amat) * vXt;
      setround(0);
      rabs = mag(infsup(rlow,rup));
   else   
      r = b - Amat * vXt; 
      rabs = mag(r);
   end;                                               
   for j = 1 : n
      % Implement other eigenvalue method !!
      Xj = Xt{j};
      lowbounds = inf_(veigsym(Xj));
      lbwork(j) = min(lowbounds);
      l(j) = sum(lowbounds < 0); 
   end                                           
   l = l(:); lbwork = lbwork(:); 
   lbminus = min(0,lbwork);                            
   if (min(lbminus) == 0) & (max(rabs) == 0)
      X = Xt;
      lb = lbwork;
   end   
   if intvalinput == 0   % non-interval data
      setround(1);   % full avoids index representation
      fU = full(vC'*vXt + (- rho)'*(lbminus .* l) + rabs'*yu); 
      setround(0);
      return;
   else    % full avoids index representation
      fU = sup(full(vC' * vXt - intval(rho)' * (lbminus .* intval(l)) ...
               + rabs' *intval(yu))); 
      return;  
   end   
end



%Algorithm with infinite dual bounds yu
if max(yu) == inf                                          
   while (~stop) & (iter <= VSDP_ITER_MAX)
      %1.step                                                                    
      [vX,J,I,N] = vuls([], [], Amat, b, inf_(Xbounds), sup(Xbounds),...
                     mid(vX),I,N); 
      if isnan(vX)
         disp('VSDPUP: system matrix may have no full rank');
         return;            
      end  
      Xwork = vsmat(vX,blkvec,0,1);                    
      %2.step                       
      for j = 1 : n
         % Implement other eigenvalue method 
         Xj = Xwork{j};
         lowbounds = inf_(veigsym(Xj));
         lb(j) = min(lowbounds);
      end  
      %3.step
      if min(lb) >= 0
         X = Xwork;
         fU = sup(full(vC' * intval(vX)));  
         return;   
      end
      %4.step 
      for j = 1 : n
         if lb(j) < 0
            k(j) = k(j) +1;
            Eps{j} = -VSDP_ALPHA^k(j) * lb(j) * II{j} + Eps{j}; 
         end   
      end 
      lb = repmat(NaN,n,1);
      %Perturbed problem
      vEps = vsvec(Eps,0,1);
      if intvalinput == 0
         bperturb = b - Amat * vEps;
      else   
         bperturb = midb - midAmat * vEps; 
      end   
      %5.step
      if intvalinput == 0
         [objt,Xt,yt,Zt,info] = mysdps(blk,A,C,bperturb,Xt,yt,Zt);
         % NaN check 
         if max(isnan(yt)) | max(isinf(yt)) 
            disp('VSDPINFEAS: SDP-solver in MYSDPS computes NaN components.'); 
            return;
         end 
      else   
         [objt,Xt,yt,Zt,info] = mysdps(blk,midA,midC,bperturb,Xt,yt,Zt);
         if max(isnan(yt)) | max(isinf(yt)) 
            disp('VSDPINFEAS: SDP-solver in MYSDPS computes NaN components.'); 
            return;
         end 
      end   
      if ((info(1) == 1) | (info(1) == 3)) 
         stop = 1;              
      end
      %6.step
      vX = sparse(vsvec(Xt,0,1)) + vEps;
      iter = iter+1;
   end   
end    

%________________________________END OF VSDPUP___________________________

