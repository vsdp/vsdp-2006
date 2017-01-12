function tests = vsdpTest
% VSDPTEST  Runs a Matlab testsuite for VSDP (version 2006).
%
%   About the tests:
%
%    - example1      contained in the original documentation of vsdplow and
%                    vsdpup
%    - example2      contained in demovsdp
%    - SDPLIB_ARCH4  contained in demovsdp
%
%   Example:
%
%       clc; table (runtests ('vsdpTest'))
%
%   See also vsdpcheck.

% Copyright 2016-2017 Kai T. Ohlhus (kai.ohlhus@tuhh.de)

tests = functiontests(localfunctions);
end

function testVSVEC(testCase)
vsvec_len = @(blk) sum(blk .* (blk + 1) ./ 2);
A = {ones(3)};
vAref = [1 sqrt(2) sqrt(2) 1 sqrt(2) 1]';
vA = vsvec(A);
verifyEqual(testCase, vA, vAref)
verifyEqual(testCase, length(vA), vsvec_len(3))
vA = vsvec(A,0);
verifyEqual(testCase, vA, vAref)
verifyEqual(testCase, length(vA), vsvec_len(3))
vAref = sparse(vAref);
verifyEqual(testCase, vsvec(A,1), vAref)
A = {[10 2 3; 2 11  4; 3 4 12]; ones(3)};
vAref = [10 4 6 11 8 12 1 2 2 1 2 1]';
vA = vsvec(A,0,2);
verifyEqual(testCase, vA, vAref)
verifyEqual(testCase, length(vA), 2*vsvec_len(3))
A = {[10 2 3; 2 11  4; 3 4 12]; sparse(ones(3))};
vA = vsvec(A,0,2);
verifyEqual(testCase, vA, vAref)
verifyEqual(testCase, length(vA), 2*vsvec_len(3))
end

function testVSMAT(testCase)
AA = {{ones(3)};
  {[10 2 3; 2 11  4; 3 4 12]; ones(3)};};
for i = 1:length(AA)
  A = AA{i};
  blk = zeros(length(A),1);
  for j = 1:length(A)
    blk(j) = length(A{j});
  end
  % test full input + full output
  % default call
  verifyEqual(testCase, isequal(A,vsmat(vsvec(A),blk)), true)
  % sparseflag = 0
  verifyEqual(testCase, isequal(A,vsmat(vsvec(A),blk,0)), true)
  verifyEqual(testCase, isequal(A,vsmat(vsvec(A,0),blk)), true)
  verifyEqual(testCase, isequal(A,vsmat(vsvec(A,0),blk,0)), true)
  % intermediate sparse vector
  verifyEqual(testCase, isequal(A,vsmat(vsvec(A,1),blk,0)), true)
  % sparseflag = 0 + different mults
  mult = [sqrt(2), 2, 1];
  for j = 1:length(mult)
    verifyEqual(testCase, ...
      isequal(A,vsmat(vsvec(A,0,mult(j)),blk,0,1/mult(j))), true)
  end
  % test sparse input + sparse output
  for j = 1:length(A)
    A{j} = sparse(A{j});
  end
  % intermediate full vector
  verifyEqual(testCase, isequal(A,vsmat(vsvec(A,0),blk,1)), true)
  % intermediate sparse vector
  verifyEqual(testCase, isequal(A,vsmat(vsvec(A,1),blk,1)), true)
end
end

function testVEIGSYM(testCase)
A = [1 2 3; 2 1 4; 3 4 5];
% test real input
lambda = veigsym(A);
verifyEqual(testCase, mid(lambda), eig(A), 'RelTol', 2*eps())
verifyEqual(testCase, rad(lambda), zeros(3,1), 'AbsTol', 1e-14)
% test interval input
lambda_inf = [-1.5170; -0.6226; 9.0495];
lambda_sup = [-1.4569; -0.5625; 9.1096];
A = midrad(A, 1e-2*ones(3));
lambda = veigsym(A);
verifyEqual(testCase, inf(lambda), lambda_inf, 'RelTol', 1e-4)
verifyEqual(testCase, sup(lambda), lambda_sup, 'RelTol', 1e-4)
verifyEqual(testCase, sup(lambda), lambda_sup, 'RelTol', 1e-4)
end

function [blk,A,C,b] = example1()
blk(1,:) = {'s'; 2};
C{1}   = [1 0; 0 1];
A{1,1} = [0 1; 1 0];
A{2,1} = [1 1; 1 1];
b = [1; 2.0001];
end

function [blk,A,C,b,DELTA] = example2()
DELTA = 1e-4;
C{1} = [ 0   1/2    0;
  1/2 DELTA   0;
  0    0   DELTA];
A{1,1} = [  0  -1/2 0;
  -1/2   0  0;
  0    0  0];
A{2,1} = [1 0 0;
  0 0 0;
  0 0 0];
A{3,1} = [0 0 1;
  0 0 0;
  1 0 0];
A{4,1} = [0 0 0;
  0 0 1;
  0 1 0];
b = [1; 2*DELTA; 0; 0];
blk(1,:) = {'s'; 3};
end

function testVSDPCHECK_example1(testCase)
[blk,A,C,b] = example1();
[m,n] = vsdpcheck(blk,A,C,b);
verifyEqual(testCase, m, 2)
verifyEqual(testCase, n, 1)
end

function testVSDPCHECK_example2(testCase)
[blk,A,C,b,~] = example2();
[m,n] = vsdpcheck(blk,A,C,b);
verifyEqual(testCase, m, 4)
verifyEqual(testCase, n, 1)
end

function testSDPT3_conversion_example1(testCase)
[blk,A,C,b] = example1();
[At,Ct] = vsdp_to_sdpt3(blk,A,C,b);
[AA,CC] = sdpt3_to_vsdp(blk,At,Ct,b);
for i = 1:length(b)
  verifyEqual(testCase, full(AA{i,1}), A{i,1})
end
verifyEqual(testCase, CC, C)
end

function testSDPT3_conversion_example2(testCase)
[blk,A,C,b,~] = example2();
[At,Ct] = vsdp_to_sdpt3(blk,A,C,b);
[AA,CC] = sdpt3_to_vsdp(blk,At,Ct,b);
for i = 1:length(b)
  verifyEqual(testCase, full(AA{i,1}), A{i,1})
end
verifyEqual(testCase, CC, C)
end

function [blk,A,C,b,X,y,Z] = testMYSDPS_example1(testCase)
[blk,A,C,b] = example1();
[objt,X,y,Z,info] = mysdps(blk,A,C,b);
verifyEqual(testCase, objt, [1, 1], 'RelTol', 1e-4)
verifyEqual(testCase, X{1}, ones(2)/2, 'RelTol', 1e-4)
verifyEqual(testCase, objt(1), trace(C{1} * X{1}), 'RelTol', 1e-4)
verifyEqual(testCase, objt(2), b'*y, 'RelTol', 1e-4)
verifyEqual(testCase, Z{1}, C{1} - A{1,1}*y(1) - A{2,1}*y(2), 'RelTol', 1e-4)
verifyEqual(testCase, info, 0)
end

function [blk,A,C,b,X,y,Z,DELTA] = testMYSDPS_example2(testCase)
[blk,A,C,b,DELTA] = example2();
[objt,X,y,Z,info] = mysdps(blk,A,C,b);
verifyEqual(testCase, objt, [-0.5 -0.5], 'RelTol', DELTA)
testX = X{1};
verifyEqual(testCase, testX(1:2,1:2), [2*DELTA, -1; -1, 1/(2*DELTA)], ...
  'RelTol', DELTA)
verifyEqual(testCase, testX(3,3) >= 0, true)
verifyEqual(testCase, y(2), -1/(4*DELTA), 'RelTol', DELTA)
verifyEqual(testCase, y([1,3,4]), [0;0;0], 'AbsTol', DELTA)
verifyEqual(testCase, objt(1), trace(C{1} * X{1}), 'RelTol', DELTA)
verifyEqual(testCase, objt(2), b'*y, 'RelTol', DELTA)
D = C{1};
for j = 1:4
  D = D - y(j) * A{j,1};
end
verifyEqual(testCase, Z{1}, D, 'RelTol', DELTA)
verifyEqual(testCase, info, 0)
end

function testVSDPUP_example1(testCase)
[blk,A,C,b,Xt,yt,Zt] = testMYSDPS_example1(testCase);
[fU,X,lb] = vsdpup(blk,A,C,b,Xt,yt,Zt);
verifyEqual(testCase, fU, 1.0001, 'RelTol', 1e-4)
verifyEqual(testCase, all (all (in (X{1}, midrad(ones(2)/2,5e-4)))), true)
verifyEqual(testCase, lb, 5e-5, 'RelTol', 1e-4)
end

function testVSDPUP_example2(testCase)
[blk,A,C,b,Xt,yt,Zt,DELTA] = testMYSDPS_example2(testCase);
[fU,X,lb] = vsdpup(blk,A,C,b,Xt,yt,Zt);
verifyEqual(testCase, fU, -0.5, 'RelTol', DELTA)
midX = [2*DELTA, -1, 0; -1, 1/(2*DELTA), 0; 0, 0, DELTA];
verifyEqual(testCase, all (all (in (X{1}, midrad(midX,DELTA)))), true)
verifyEqual(testCase, lb, 1.289076539560735e-8, 'AbsTol', DELTA)
end

function testVSDPUP_example2_finite_bnd(testCase)
[blk,A,C,b,Xt,yt,Zt,DELTA] = testMYSDPS_example2(testCase);
yu = 1e5 * [1 1 1 1]';
[fU,X,lb] = vsdpup(blk,A,C,b,Xt,yt,Zt,yu);
% In this case: VSDPUP computes rigorously the finite upper bound without
% verifying primal feasibility.
verifyEqual(testCase, fU, -0.5, 'RelTol', DELTA)
verifyEqual(testCase, isnan(X), true)
verifyEqual(testCase, isnan(lb), true)
end

function testVSDPLOW_example1(testCase)
[blk,A,C,b,Xt,yt,Zt] = testMYSDPS_example1(testCase);
[fL,Y,dl] = vsdplow(blk,A,C,b,Xt,yt,Zt);
verifyEqual(testCase, fL, 1.0001, 'RelTol', 1e-4)
verifyEqual(testCase, Y, [-1; 1], 'RelTol', 1e-4)
verifyEqual(testCase, dl, 9.3912e-11, 'AbsTol', 1e-9)
end

function testVSDPLOW_example2(testCase)
[blk,A,C,b,Xt,yt,Zt,DELTA] = testMYSDPS_example2(testCase);
[fL,Y,dl] = vsdplow(blk,A,C,b,Xt,yt,Zt);
verifyEqual(testCase, fL, -0.5, 'RelTol', DELTA)
verifyEqual(testCase, Y(2), -1/(4*DELTA), 'RelTol', DELTA)
verifyEqual(testCase, Y([1,3,4]), [0;0;0], 'AbsTol', DELTA)
verifyEqual(testCase, dl, 4.319659051028185e-13, 'AbsTol', 1e-9)
end

function testVSDPLOW_example2_finite_bnd(testCase)
[blk,A,C,b,Xt,yt,Zt,DELTA] = testMYSDPS_example2(testCase);
xu = 1e5;
% FIXME: tiny pertubation neccessary to work with approximation from SDPT3-4.0
yt(2) = yt(2) - DELTA/4;
[fL,Y,dl] = vsdplow(blk,A,C,b,Xt,yt,Zt,xu);
verifyEqual(testCase, fL, -0.5, 'RelTol', DELTA)
verifyEqual(testCase, Y(2), -1/(4*DELTA), 'RelTol', DELTA)
verifyEqual(testCase, Y([1,3,4]), [0;0;0], 'AbsTol', DELTA)
verifyEqual(testCase, dl, 4.319659051028185e-13, 'AbsTol', 1e-9)
end

function testSDPLIB_ARCH4(testCase)
n = 2;
m = 174;
s = [161; 174];
% test sdpa_to_vsdp
[blk,A,C,b] = sdpa_to_vsdp(fullfile('sdplib_vsdp', 'arch4.dat-s'));
verifyEqual(testCase, blk, {'s', s(1); 's', s(2)})
verifyEqual(testCase, length(b), m)
for j = 1:n
  verifyEqual(testCase, size(C{j}), [s(j), s(j)])
  for i = 1:m
    verifyEqual(testCase, size(A{i,j}), [s(j), s(j)])
  end
end
% test vsdpcheck
[m,n] = vsdpcheck(blk,A,C,b);
verifyEqual(testCase, m, 174)
verifyEqual(testCase, n, 2)
% test SDPT3 conversion
[At,Ct] = vsdp_to_sdpt3(blk,A,C,b);
[AA,CC] = sdpt3_to_vsdp(blk,At,Ct,b);
verifyEqual(testCase, CC, C)
verifyEqual(testCase, AA, A, 'RelTol', 1e-8)
% test SDPT3 mysdps
[objt,Xt,yt,Zt,info] = mysdps(blk,A,C,b);
verifyEqual(testCase, size(objt), [1, 2])
verifyEqual(testCase, size(Xt), [2, 1])
verifyEqual(testCase, size(yt), size(b))
verifyEqual(testCase, size(Zt), [2, 1])
verifyEqual(testCase, objt(1), trace(C{1} * Xt{1}), 'RelTol', 1e-4)
verifyEqual(testCase, objt(2), b'*yt, 'RelTol', 1e-4)
verifyEqual(testCase, info, 0)
% test SDPT3 vsdpup
[fU,X,lb] = vsdpup(blk,A,C,b,Xt,yt,Zt);
verifyEqual(testCase, objt(2) <= fU, true)
verifyEqual(testCase, isscalar(X) && isnan(X), false);
verifyEqual(testCase, any(isnan(lb)), false);
% test SDPT3 vsdplow
[fL,Y,dl] = vsdplow(blk,A,C,b,Xt,yt,Zt);
verifyEqual(testCase, fL <= objt(1), true)
verifyEqual(testCase, isscalar(Y) && isnan(Y), false);
verifyEqual(testCase, any(isnan(dl)), false);
end
