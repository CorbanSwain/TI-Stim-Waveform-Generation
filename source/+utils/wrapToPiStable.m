function y = wrapToPiStable(x, tol)
y = wrapTo2PiStable(x + pi, tol) - pi;
end