function scaleFactor = getScaleFactorFromAmp(a)
% returns a factor with which to scale a plot's limits, for example,
%  based on a wave amplitude.

a = abs(a);
if a < eps()
    % if the amplitude is near zero, set the scale factor to 1
    scaleFactor = 1;
    return
end
scaleFactor = a;
end