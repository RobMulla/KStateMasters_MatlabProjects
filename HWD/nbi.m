%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NORMALIZED BOUNDARY INTERSECTION ALGORITHM %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xpf] = nbi(N, x0,lb, ub, options)

global z f1 f2 nv
% Finds the pareto optimal solutions
% of bi-objective optimization problems
%
% Argument List:
% --------------
%
% f1  : Handle to objective function 1 | Calling program should pass
% f2  : Handle to objective function 2 | @fi and @f2 as arguments to nbi
%
% N   : Desired number of Pareto optimal solutions (must be > 2)
%
% x0  : D x 1 vector of starting value of x to be used by fmincon().
%
% Ane : Matrix for inequality constraint | Inequality onstraint:
% bne : Vector for inequality constraint |   Ane*x <= bne
%
% Aeq : Matrix for equality constraint   | Equality onstraint:
% beq : Vector for equality constraint   |   Aeq*x = beq
%
% lb  : D x 1 vector of lower bounds     | Bounds:
% ub  : D x 1 vector of lower bounds     | lb <= x <= ub
%
% optimset : parameters for fmincon()
%
% Output list:
% ------------
%
% xpf : D x N Matrix of Pareto optimal solutions
% fpf : N x 2 Matrix of values of f1 and f2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the two anchor points 'a' and 'b'
[a, f1a, eflag1] = fmincon(f1, [0 x0], [], [], [], [], [-100 lb], [100 ub], @constraints, options);
[b, f2b, eflag2] = fmincon(f2, [0 x0], [], [], [], [], [-100 lb], [100 ub], @constraints, options);
f2a = f2(a);
f1b = f1(b);

% Define NBI parameters phi and nv (normalized vector)
phi = [f1a f1b; f2a f2b];
nv = [1 1]./sqrt(2); nv = nv';

% Find other optimal solutions
N = N - 2;
for i = 1 : N
    w = i*[1 0]' + (N-i)*[0 1]';
    w = w./N; %sum(w) should be 1
    z = phi*w;
    x = fmincon(@nbif, [0 x0], [], [], [], [], [-1000 lb], [1000 ub], @nbic, options);
    xpf(:, i) = x(2:numel(x));
    % INSERT CODE HERE
end

% INSERT CODE HERE
end

function f = nbif(x)
global z f1 f2 nv
% lamda is assumed to be x(1) such that z + lamda*nv = [f1 f2]'
f = -abs(x(1));
end

function [c, ceq] = nbic(x)
global z f1 f2 nv
[c ceq] = constraints(x);
% lamda is assumed to be x(1) such that z + lamda*nv = [f1 f2]'
ceq1 = z + x(1).*nv - [f1(x); f2(x)];
ceq = [ceq; ceq1];
end
