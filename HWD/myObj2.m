function f = myObj2(x)
n = numel(x);
x = x(2:n);
f = 9*x(1) - (x(2) - 1).^2;
end