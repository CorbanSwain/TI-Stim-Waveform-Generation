function output = halfSinSigmoid(x, xStart, xEnd, yStart, yEnd)

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
isClip1 = x <= x1;
isClip2 = x >= x2;
output(isClip1) = 0;
output(isClip2) = 1;

% main definition
isFxDefined = ~isClip1 & ~isClip2;
if any(isFxDefined)
    deltaX = x2 - x1;
    halfSinFx = @(x) (sin((2 * pi * (x - x1) / 2 / deltaX) ...
        - (pi / 2)) / 2) + (1 / 2);
    output(isFxDefined) = halfSinFx(x(isFxDefined));
end

output = y1 + (output * (y2 - y1));
end