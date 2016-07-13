function f = myObj1(x)
x = x(2:numel(x));
f = (x(1) - 2).^2 + (x(2) - 1).^2 + 2;
end