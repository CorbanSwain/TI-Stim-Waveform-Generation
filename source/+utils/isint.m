%ISINT Returns true if the value is numerically an integer.

% Corban Swain , 2019

function output = isint(x)
output = isinteger(x) | mod(x, 1) <= eps(x);
end

