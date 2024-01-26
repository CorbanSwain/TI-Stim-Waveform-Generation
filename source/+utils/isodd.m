function y = isodd(x)
y = utils.isint(x) & logical(mod(x, 2));
end