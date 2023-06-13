function y = iseven(x)
y = csmu.isint(x) & ~mod(x, 2);
end