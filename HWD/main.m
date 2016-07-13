function [xpf, fpf] = main()
global f1 f2
f1 = @myObj1;
f2 = @myObj2;
options=optimset('TolCon', 1e-5, 'TolX', 1e-5, 'MaxFunEvals', 10000, 'Algorithm', 'active-set');
x0 = 0*ones(1, 2);
lb = -20*ones(1,2);
ub = 20*ones(1,2);
N = 20;
[xpf, fpf] = nbi(N, x0, lb, ub, options);
plot(fpf(:,1), fpf(:,2), '.');