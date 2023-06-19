function output = quadraticSigmoid(x, xStart, xEnd, yStart, yEnd)

% handle inputs
p.addRequired('x', @isnumeric);
p.addRequired('xStart', @(x) isnumeric(x) && isscalar(x));
p.addRequired('xEnd', @(x) isnumeric(x) && isscalar(x));
p.addRequired('yStart', @(x) isnumeric(x) && isscalar(x));
p.addRequired('yEnd', @(x) isnumeric(x) && isscalar(x));
p.parse(x, xStart, xEnd, yStart, yEnd);

if ~(xStart < xEnd)
    error('xEnd must be greater than xStart.');
end

% renaming
x1 = xStart;
x2 = xEnd;
y1 = yStart;
y2 = yEnd;

% output
output = zeros(size(x), 'like', x);

% clipped values
isPreTeepee = x <= x1;
output(isPreTeepee) = 0;
isPostTeepee = x >= x2;
output(isPostTeepee) = 1;
isConstRegion = isPreTeepee | isPostTeepee;

% main definition
if ~all(isConstRegion)
    tMidpoint = (x1 + x2) / 2;
    dx = -x2 + x1;
    
    isXPrePeak = (x < tMidpoint) & ~isConstRegion;
    isXPostPeak = ~isXPrePeak & ~isConstRegion;
    
    output(isXPrePeak) = ...
        1.0 ./ dx .^ 2 .* (-x(isXPrePeak) + x1) .^ 2 .* 2.0;
    output(isXPostPeak) = ...
        (1.0 ./ dx .^ 2 .* (-x(isXPostPeak) + x2) .^ 2 .* -2.0) + 1;
end

output = y1 + (output * (y2 - y1));
end