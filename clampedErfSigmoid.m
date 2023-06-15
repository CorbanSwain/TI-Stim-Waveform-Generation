function output = clampedErfSigmoid(x, xStart, xEnd, yStart, yEnd, ...
    varargin)

% handle inputs
p = inputParser();
p.addRequired('x', @isnumeric);
p.addRequired('xStart', @(x) isnumeric(x) && isscalar(x));
p.addRequired('xEnd', @(x) isnumeric(x) && isscalar(x));
p.addRequired('yStart', @(x) isnumeric(x) && isscalar(x));
p.addRequired('yEnd', @(x) isnumeric(x) && isscalar(x));
p.addParameter('NumSigmas', 2, @(x) isnumeric(x) && (x > 0));
p.parse(x, xStart, xEnd, yStart, yEnd, varargin{:})
numSigmas = p.Results.NumSigmas;

if ~(xStart < xEnd)
    error('xStart must be greater than xEnd.');
end

% renaming
x1 = xStart;
x2 = xEnd;
y1 = yStart;
y2 = yEnd;

% output
output = zeros(size(x), 'like', x);

% handle clipped values
isClip1 = x < x1;
isClip2 = x > x2;
output(isClip1) = y1;
output(isClip2) = y2;

isFxDefined = ~isClip1 & ~isClip2;
if ~any(isFxDefined)
    return
end

% main implementation
deltaX = x2 - x1;
sigma = deltaX / (numSigmas * 2);
mu = (x1 + x2) / 2;

cdfFx = @(x) normcdf(x, mu, sigma);
pdfFx = @(x) normpdf(x, mu, sigma);

cdf_x1 = cdfFx(x1);
cdf_x2 = cdfFx(x2);
pdf_x1 = pdfFx(x1);

aucBase = cdf_x2 - cdf_x1 - (pdf_x1 * deltaX);
sigmoidFx = @(x) y1 + ((y2 - y1) / aucBase ...
    * (cdfFx(x) - cdf_x1 - (pdf_x1 * (x - x1))));

output(isFxDefined) = sigmoidFx(x(isFxDefined));
end