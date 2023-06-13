function y = isodd(x)
y = csmu.isint(x) & logical(mod(x, 2));
end