function y = iseven(x)
y = utils.isint(x) & ~mod(x, 2);
end