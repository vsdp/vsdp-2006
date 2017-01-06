function [X, J, I, N] = vuls(A, a, B, b, xl, xu, x0, I, N)
% VULS: Verification for Underdetermined Linear Systems
%       of inequalities and equations:
%              A * x <= a,
%              B * x  = b,
%              xl <= x <= xu
% The input A (m*n matrix), a (m-vector), B (p*n matrix), and  b (p-vector)
% can be real or interval quantities; the simple bounds xl, xu  must
% be real, but may be infinite; the approximate solution x0 must be real.
% The optional input vector I must contain p indices out of {1,...,n} such
% that the submatrix B(:,I) is nonsingular; N contains the remaining n-p
% nonbasic indices.
%
% Output:
%      X   a box (n-interval vector), containing for every real
%          input (A,a,B,b)  within the interval input data a  solution
%          x of the above system, provided J is empty. Especially,
%          existence of solutions is verified, and moreover
%          X is computed close to x0 in a specified manner; for details see
%          C. Jansson, Rigorous Lower and Upper Bounds in Linear
%          Programming,  SIAM J. OPTIM. Vol.14, No.3, pp. 914-935.
%
%          X := intval(repmat(NaN,n,1)), J = NaN; I = NaN; N = NaN;
%          if existence of solutions cannot be proved, and verified finite
%          bounds cannot be computed. This is the case, if
%          (i) B has no full rank, or (ii) the linear interval solver
%          VERIFYLSS cannot compute rigorous bounds, or (iii) the box
%          of simple bounds [xl,xu] has no appropriate interior. Otherwise,
%      J   returns the vector of indices of violated inequalities
%          for X:
%             J.ineqlin:  violated row indices of A * X <= b,
%             J.lower: violated indices of  xl <= X ,
%             J.upper: violated indices of  X <= xu .
%      I   index vector such that the p*p submatrix B(:,I) is
%          nonsingular,
%      N   index vector containig the n-p indices of 1:n not in I.
%
% EXAMPLE:
% A = [1 1 1 1]; a =infsup(2.9,3.1) ;
% B = [0 1 0 infsup(0.9,1.1)]; b = 2;
% xl = [0 1 0 0]'; xu = [4 1 1 2]'; x0 = [0 1 0 1];
% xl = [0 1 0 0]'; xu = [4 1 1 2]'; x0 = [0 1 0 1];
% [X,J,I,N] = vuls(A, a, B, b, xl, xu, x0)
% gives the output:
% intval X =
% [    0.0000,    0.0001]
% [    1.0000,    1.0000]
% [    0.0000,    0.0001]
% [    0.8878,    1.1122]
% J =
%     ineqlin: []
%       lower: [0x1 double]
%       upper: [0x1 double]
% I =
%      4
% N =
%      1
%      2
%      3
%
%Used Files: verifylss (INTLAB)
% written   01/12/05   Christian Jansson
% modified  03/21/06
% Reference: C. JANSSON, Termination and Verification for
%            Ill-posed Semidefinite Programming Problems,
%            to appear
% http://optimization-online.org/DBHTML/|2005/06/1150.html.

% Calling the m-file of the global variables
SDP_GLOBALPARAMETER;

%Selection variable of the precision for computing the defect;
global VSDP_HIGHER_PREC;

%Choice, whether in VULS the interval system is solved as a full
%system, or as a sparse system by computing the psd matrix B*B';
%VSDP_CHOICE_FULL = 0 in the latter case.
global VSDP_CHOICE_FULL;

A = sparse(A);
B = sparse(B);

% Transformation to intervals and midpoints
if isa(A,'intval')
  mA = mid(A);
else
  mA = A; A = intval(mA);
end
a = a(:);
if isa(a,'intval')
  ma = mid(a);
else
  ma = a; a = intval(ma);
end
if isa(B,'intval')
  mB = mid(B);
else
  mB = B; B = intval(mB);
end
b = b(:);
if isa(b,'intval')
  mb = mid(b);
else
  mb = b; b = intval(mb);
end


n = length(x0);
p = length(b);
J.ineqlin = [];
J.lower = [];
J.upper = [];
if nargin < 9 || (length(I) ~= p) || (length(N) ~= n-p)
  I = []; N =[];
end

x0 = x0(:);
xl = xl(:);
xu = xu(:);
if any((xu-xl)<0)
  disp('VULS: simple bounds are not feasible');
  X = NaN(n,1);
  J.ineqlin = NaN;
  J.lower = NaN;
  J.upper = NaN;
  I = NaN; N = NaN;
  return;
end

% Projection of x0 into the interior of [xl,xu]
Iwork=find((xu-xl)> 20*eps);
xlint = xl;
xuint = xu;
xlint(Iwork) = xl(Iwork) + 10*eps*abs(xl(Iwork)) + 5*eps;
xuint(Iwork) = xu(Iwork) - 10*eps*abs(xl(Iwork)) - 5*eps;
if any((xuint < xlint) < 0) || (length(Iwork) < p)
  X = NaN(n,1);
  J.ineqlin = NaN;
  J.lower = NaN;
  J.upper = NaN;
  I = NaN; N = NaN;
  disp('VULS: simple bounds are degenerate');
  return;
end
x0(x0<xlint) = xlint(x0<xlint);
x0(x0>xuint) = xuint(x0>xuint);
X = intval(x0);                                                 %X
%size(mB)
% Enclosure X
if ~isempty(mB)
  % Determine the basis with lu decomposition
  if 1
    if isempty(I)
      [~,~,P] = lu(mB(:,Iwork)'); %P*mB(:,Iwork)' - L*U, Iwork
      piv = P * Iwork;            %piv
      I = piv(1:p)';              %I, full(mB(:,I)),  RANK=rank(full(mB(:,I)))
      if p == n
        N = [];
      else   % Finding the nonbasic indices N
        % This way yields out of memory for lrge problems
        %     I = I(:);
        %     Iv=repmat(I',n,1);
        %     Index=repmat((1:n)',1,length(I));
        %     Ifind = sum((Iv == Index)',1);     I,Iv,Index,Ifind
        %     Ifind = Ifind(:);
        Ifind = zeros(n,1);
        for i = 1 : n
          for j = 1 : length(I)
            if i == I(j)
              Ifind(i) = 1;
              break;
            end
          end
        end
        N = find(~Ifind);              % P,Iwork,piv, Iv,Index,I,N
      end
    end
  else   % Determinate the basis with qr decomposition
    if isempty(I)
      [Q,R,P] = qr(full(mB(:,Iwork)));  % P*mB(:,Iwork)' - (Q*R)',Iwork, full(P)
      piv = P' * Iwork;                 % full(mB), P, piv
      I = piv(1:p);                     % I, full(mB(:,I)),  RANK=rank(full(mB(:,I)))
      if p == n
        N = [];
      else   % Finding the nonbasic indices N
        I = I(:);
        Iv=repmat(I',n,1);
        Index=repmat((1:n)',1,length(I));
        Ifind = sum((Iv == Index)',1);
        Ifind = Ifind(:);
        N = find(~Ifind);             % P, I ,Iwork,piv, Iv,Index,N
      end
    end
  end
  % Solving underdetermined interval system
  BI = B(:,I);                                   %B,BI,N,VSDP_HIGHER_PREC
  if isempty(N)
    bI = b;
  else
    if VSDP_HIGHER_PREC == 0
      bI = b - B(:, N) * x0(N);                  % maxrad=max(rad(bI))
    else          % Higher precision for bI NOCH IMPLEMENTIEREN * TESTEN!!!!
      bI = -dot_(mB(:,N),x0(N),-1,mb,-1);
    end
  end
  if VSDP_CHOICE_FULL == 0
    bI = BI'*bI;
    BI = BI'*BI;
    XI = verifylss(BI,bI);
  else
    XI = verifylss(full(BI),full(bI));
  end
  if isnan(XI)
    X = NaN(n,1);
    J.ineqlin = NaN;
    J.lower = NaN;
    J.upper = NaN;
    I = NaN; N = NaN;
  else
    X(I) = XI;
  end
end                                      %maxradX=max(rad(X)), maxabsX=max(abs(mid(X)))
%X-x0
%Test of inequalities
if ~isempty(a)
  J.ineqlin = find(~(A * X <= a));
end
J.lower = find(~(xl <= X));
J.upper = find(~(X <= xu));

end
