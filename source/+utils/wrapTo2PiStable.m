function y = wrapTo2PiStable(x, tol)
y = mod(x, 2 * pi);
y(y > ((2 * pi) - tol)) = 0;
end