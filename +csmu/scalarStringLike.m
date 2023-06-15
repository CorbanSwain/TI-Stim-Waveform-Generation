function isvalid = scalarStringLike(x)
if isstring(x)
   isvalid = isscalar(x);
elseif ischar(x)
   isvalid = isequal(size(x), [0, 0]) ...
      || (isvector(x) && ismatrix(x) && size(x, 1) == 1);
else
   isvalid = false;
end
end