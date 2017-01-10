function tests = vsdpTest
% VSDPTEST  Runs a Matlab testsuite for VSDP (version 2006).
%
%   About the tests:
%
%      example1  contained in the original documentation of vsdplow and vsdpup
%      example2  contained in demovsdp
%
%   Example:
%
%       clc; table (runtests ('vsdpTest'))
%
%   See also vsdpcheck.

% Copyright 2016-2017 Kai T. Ohlhus (kai.ohlhus@tuhh.de)

tests = functiontests(localfunctions);
end

function [blk,A,C,b,X,y,Z] = testMYSDPS_example1(testCase)
blk(1,:) = {'s'; 2};
C{1}   = [1 0; 0 1];
A{1,1} = [0 1; 1 0];
A{2,1} = [1 1; 1 1];
b = [1; 2.0001];
[objt,X,y,Z,info] = mysdps(blk,A,C,b);
verifyEqual(testCase, objt, [1, 1], 'RelTol', 1e-4)
verifyEqual(testCase, X{1}, ones(2)/2, 'RelTol', 1e-4)
verifyEqual(testCase, objt(1), trace(C{1} * X{1}), 'RelTol', 1e-4)
verifyEqual(testCase, objt(2), b'*y, 'RelTol', 1e-4)
verifyEqual(testCase, Z{1}, C{1} - A{1,1}*y(1) - A{2,1}*y(2), 'RelTol', 1e-4)
verifyEqual(testCase, info, 0)
end

function [blk,A,C,b,X,y,Z] = testMYSDPS_example2(testCase)
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
DELTA = 1e-4;
[blk,A,C,b,Xt,yt,Zt] = testMYSDPS_example2(testCase);
[fU,X,lb] = vsdpup(blk,A,C,b,Xt,yt,Zt);
verifyEqual(testCase, fU, -0.5, 'RelTol', DELTA)
midX = [2*DELTA, -1, 0; -1, 1/(2*DELTA), 0; 0, 0, DELTA];
verifyEqual(testCase, all (all (in (X{1}, midrad(midX,DELTA)))), true)
verifyEqual(testCase, lb, 1.289076539560735e-8, 'AbsTol', DELTA)
end

function testVSDPUP_example2_finite_bnd(testCase)
DELTA = 1e-4;
[blk,A,C,b,Xt,yt,Zt] = testMYSDPS_example2(testCase);
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
DELTA = 1e-4;
[blk,A,C,b,Xt,yt,Zt] = testMYSDPS_example2(testCase);
[fL,Y,dl] = vsdplow(blk,A,C,b,Xt,yt,Zt);
verifyEqual(testCase, fL, -0.5, 'RelTol', DELTA)
verifyEqual(testCase, Y(2), -1/(4*DELTA), 'RelTol', DELTA)
verifyEqual(testCase, Y([1,3,4]), [0;0;0], 'AbsTol', DELTA)
verifyEqual(testCase, dl, 4.319659051028185e-13, 'AbsTol', 1e-9)
end

function testVSDPLOW_example2_finite_bnd(testCase)
DELTA = 1e-4;
[blk,A,C,b,Xt,yt,Zt] = testMYSDPS_example2(testCase);
xu = 1e5;
% FIXME: tiny pertubation neccessary to work with approximation from SDPT3-4.0
yt(2) = yt(2) - DELTA/4;
[fL,Y,dl] = vsdplow(blk,A,C,b,Xt,yt,Zt,xu);
verifyEqual(testCase, fL, -0.5, 'RelTol', DELTA)
verifyEqual(testCase, Y(2), -1/(4*DELTA), 'RelTol', DELTA)
verifyEqual(testCase, Y([1,3,4]), [0;0;0], 'AbsTol', DELTA)
verifyEqual(testCase, dl, 4.319659051028185e-13, 'AbsTol', 1e-9)
end