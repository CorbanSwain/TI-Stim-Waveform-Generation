function phaseFcn = generateFMPhaseFcn(timeSteps, freqSteps, ...
    modulationTime, varargin)
%GENERATEFMPHASEFCN  Creates a function for the phase of an FM signal.
%
%   phaseFcn = generateFMPhaseFcn(timeSteps, freqSteps, modulationTime)
%   returns the function_handle phaseFcn which will return the phase of a
%   frequency modulated signal that linearly modulates from frequency
%   freqSteps(i) Hz to freqSteps(i+1) Hz over ther time interval
%   [timeSteps(i) s, (timeSteps(i) + modulationTime) s]. The phase function
%   will begin with frequency freqSteps(1) at time 0 s and is defined for
%   any time value between 0 s and timeSteps(end) s, inclusive. The
%   function will perform (length(timeSteps) - 1) modulations), beginning
%   at freqSteps(1) Hz at time 0 s and ending at freqSteps(end) Hz at time
%   timeSteps(end) s. timeSteps must be a monotonically increasing numeric
%   vector with timeSteps(1) > 0; freqSteps must be a positive-valued
%   numeric vector with the same length as timeSteps.
%
%   The phase will be chosen for each component such that it is piecewise
%   continuous with no discontinuitues in phaseFcn.
%
%   If a frequency of 0 is provided the phase will remain constant. All
%   values of freqSteps must be non-negative.
%
%   modulationTime can be a scalar or a vector with length
%   (length(timeSteps) - 1); in the case of supplying a vector, the
%   function will perform the modulation beginning at timeSteps(i) over
%   modulationTime(i) seconds. If a modulationTime of 0 is chosen, no
%   smooth ramp in frequency will be applied (i.e. there may be a
%   discontinuity in the first derivative of the phaseFcn).
%
%   generateFMPhaseFcn(__, 'FMScale', fmScale) The frequency modulations
%   will be perfomed with either a ramp in frequency-space or a ramp in
%   log(frequency) space. Valid values for fmScale are:
%      - 'logf'        : ramp in log(frequency) space
%      - 'f' (default) : ramp in frequency space
%
%   If fmScale is 'log' than all values for freqSteps must be > 0. 
%
%   generateFMPhaseFcn(__, 'FMShape', fmShape) The frequency modulations
%   will be perfomed with a smooth ramp having either a sinusoidal or a 
%   linear shape. Valid values for fmShape are:
%      - 'sin' (default) : sinusoidal ramp in fmScale space
%      - 'lin'           : inear ramp in fmScale space
%  
%   generateFMPhaseFcn(__, 'PhaseOffsetSteps', poSteps, 'RefPhaseFcn',
%   refPhaseFcn) If specified, the phase will also be modulated during each
%   frequency modulation to an offset from a reference phase, given by
%   refPhaseFcn, by poSteps(i + 1) radians at timeSteps(i) +
%   modulationTime. The first signal, beginning at time 0 s will begin at 0
%   s with phase offset poSteps(1) with no modulation. refPhaseFcn must be
%   provided if poSteps is set and must be defined on the closed interval
%   [0, timeSteps(end)]. If phaseOffsetSteps is not specified the phase
%   will simply be continuous with no additional adjustment.
%
%   If modulationTime=0 when poSteps is provided, there may be a
%   discontiuity in phaseFcn.
% 
%   generateFMPhaseFcn(__, 'Sigmoid', sigmoidStr, 'SigmoidArgs',
%   sigmoidArgs) The phase modulation will be applied using sigmoid
%   function corresponding to sigmoidStr supplied with parameters from
%   sigmoidArgs. Valid values for sigmoidStr are:
%      - 'erf' (default) : clamped error function (the gaussian CDF)
%      - 'quad'          : piecewise quadratic function
%      - 'sin'           : half of a sinusiod waveform
%
%   generateFMPhaseFcn(__, 'Debug', doDebug) will run the generation code
%   in debug mode if doDebug is true.
%
%   Essentially the sequence can be understood as follows:
%      1. The function begins at time t=0, having frequency f=fq(1), and
%         phase offset po(1). 
%      2. Beginning at t=tm(1) the function will modulate to
%         fq(2) and po(2). modTime(1) seconds will be taken to perform 
%         the modulation (i.e. the modulation occurs over the interval 
%         [t(1), t(1) + modTime(1)]).
%      3. Likewise, beginning at t=tm(i) the function will modulate to 
%         fq(i + 1) and po(i + 1) over modTime(i).
%      4. At t=tm(n-1) the function makes its final modulation to fq(n) and
%         po(n) over modTime(n - 1).
%      5. Finally, at t=tm(n) the function ends and is no longer defined 
%         for greater time inputs.
%   Where:
%      tm      : timeSteps
%      fq      : freqSteps
%      po      : PhaseOffsetSteps
%      modTime : modulationTime
%      n       : length(timeSteps)
%
%   See also CLAMPEDERFSIGMOID, HALFSINSIGMOID, QUADRATICSIGMOID.

% Corban Swain, June 2023

%% Input Handling
p = inputParser();
p.addRequired('timeSteps', @(x) isnumeric(x) && isvector(x));
p.addRequired('freqSteps', @(x) isnumeric(x) && isvector(x));
p.addRequired('modulationTime', @(x) isnumeric(x) && isvector(x));
p.addParameter('FMScale', 'f', @(x) any(strcmpi(x, {'logf', 'f'})));
p.addParameter('FMShape', 'sin', @(x) any(strcmpi(x, {'lin', 'sin'})));
p.addParameter('RefPhaseFcn', [], @(x) isa(x, 'function_handle'));
p.addParameter('PhaseOffsetSteps', [], @(x) isnumeric(x) && isvector(x));
p.addParameter('Sigmoid', 'erf', ...
    @(x) any(strcmpi(x, {'erf', 'quad', 'sin'})));
p.addParameter('SigmoidArgs', {}, @(x) iscell(x) && isvector(x));
p.addParameter('Debug', false, @(x) isscalar(x) && islogical(x));
p.parse(timeSteps, freqSteps, modulationTime, varargin{:});

fmScale = p.Results.FMScale;
fmShape = p.Results.FMShape;
refPhaseFcn = p.Results.RefPhaseFcn;
poSteps = wrapToPi(p.Results.PhaseOffsetSteps);
hasPO = ~isempty(poSteps);
sigmoidStr = p.Results.Sigmoid;
sigmoidArgs = p.Results.SigmoidArgs;
doDebug = p.Results.Debug;

numSteps = length(timeSteps);
numMods = numSteps - 1;
lengthCompare = [length(freqSteps)];
if hasPO
    lengthCompare = [lengthCompare, length(poSteps)];
end

if ~all(numSteps == lengthCompare)
    error(['Inconsistent vector lengths between timeSteps, freqSteps', ...
        ' and phaseOffsetSteps; vectors must have the same length.']);
end

if ~isscalar(modulationTime) && ~(length(modulationTime) == numMods)
    error(['modulationTime must be scalar or a vector with length', ...
        ' (length(timeSteps) - 1).'])
end
modTimes = ...
    modulationTime(:)' .* ones(1, numMods, 'like', modulationTime);

timeDeltas = diff([0, timeSteps(:)']);
if ~all(timeDeltas > 0)
    error('The vector [0, timeSteps] must be monotonically increassing.')
end

croppedTimeDeltas = timeDeltas(2:end);
isModTimeTooLong = modTimes > croppedTimeDeltas;
if any(isModTimeTooLong)
    warning(['Detected modulation times longer than the', ...
        ' corresponding time step for %d of %d modulation steps.', ...
        ' Clipping these values to the duration of the', ...
        ' correspomding time step.'], ...
        sum(isModTimeTooLong), numMods)    
    modTimes(isModTimeTooLong) = croppedTimeDeltas(isModTimeTooLong);
end

if hasPO && isempty(refPhaseFcn)
    error('refPhaseFcn must be passed if phaseOffsetSteps is specified.')
end

if strcmpi(fmScale, 'log') && ~all(freqSteps > 0)
    error('If fmScale is ''log'', all values of freqSteps must be > 0.')
end

%% Define Unique Frequency Moulation Functions
freqMods = [freqSteps(1:(end - 1))', freqSteps(2:end)'];
freqModsWithTime = [num2cell(freqMods, 2), num2cell(modTimes(:))];

unqFreqModsWithTime = cell(0, 2);
unqFreqModIds = zeros(1, numMods);
for iMod = 1:numMods
    [fm, modTime] = freqModsWithTime{iMod, :};
    
    % check if current fm, modTime pair is already saved
    didMatch = false;
    for iUnqMod = 1:size(unqFreqModsWithTime, 1)
        [uFM, uModTime] = unqFreqModsWithTime{iUnqMod, :};        
        if all(fm == uFM) && modTime == uModTime
            unqFreqModIds(iMod) = iUnqMod;
            didMatch = true;
            break
        end
    end

    if didMatch
        continue
    end

    % store the new unique entry
    unqFreqModsWithTime = [unqFreqModsWithTime; {fm}, {modTime}];
    unqFreqModIds(iMod) = size(unqFreqModsWithTime, 1);
end
nUnqFreqMods = size(unqFreqModsWithTime, 1);

switch lower(fmShape)
    case 'lin'
        doSinFM = false;
    case 'sin'
        doSinFM = true;
    otherwise
        error('Unexpected value for fmShape encountered, ''%s''', ...
            fmShape)
end

switch lower(fmScale)
    case 'logf'
        if doSinFM
            baseFMPhaseFcn = @logfSinFMPhaseRamp;
        else
            % function returning phase as a function of time for performing
            % a LINEAR ramp in LOG(frequency) space from frequency f1 Hz to
            % f2 Hz between time t1 s and t2 s
            baseFMPhaseFcn = @(t, t1, t2, f1, f2) ...
                (f1 * pi * ((f1 / f2) .^ ((t - t1) / (t1 - t2))) ...
                 * (t1 - t2) * 2) / log(f1 / f2);
        end

    case 'f'
        if doSinFM
            % function returning phase as a function of time for performing
            % a SINUSOIDAL ramp from frequency f1 Hz to f2 Hz between time
            % t1 s and t2 s
            baseFMPhaseFcn = @(t, t1, t2, f1, f2) 2 * pi ...
                * ((t * (f1 + f2) / 2) ...
                   + (sin((pi * (t - t1)) / (t1 - t2)) ...
                      * (f1 - f2) * (t1 - t2)) / (2 * pi));
        else
            % function returning phase as a function of time for performing
            % a LINEAR ramp from frequency f1 Hz to f2 Hz between time t1 s
            % and t2 s
            baseFMPhaseFcn = @(t, t1, t2, f1, f2) ...
                2 * pi * (f1 * t + ((f2 - f1) * ((t - t1) .^ 2) ...
                / (2 * (t2 - t1))));
        end

    otherwise
        error('Unexpected value for fmScale encountered, ''%s''', ...
            fmScale)
end

fmFcnByUnqFreqMod = cell(1, nUnqFreqMods);
for iUnqMod = 1:nUnqFreqMods
    [fm, modTime] = unqFreqModsWithTime{iUnqMod, :};   
    if modTime == 0
        fmPhaseFcn = @(t, t1) 2 * pi * ((fm(1) + fm(2)) / 2) * (t - t1);
    elseif diff(fm) == 0
        % frequency does not change
        fmPhaseFcn = @(t, t1) 2 * pi * fm(1) * (t - t1);
    else
        fmPhaseFcn = @(t, t1) ...
            baseFMPhaseFcn(t, t1, t1 + modTime, fm(1), fm(2));           
    end
     fmFcnByUnqFreqMod{iUnqMod} = fmPhaseFcn;
end

%% Determine Sigmoid Function
if hasPO
    switch lower(sigmoidStr)
        case 'erf'
            sigmoidFuncHand = @clampedErfSigmoid;
        case 'sin'
            sigmoidFuncHand = @halfSinSigmoid;
        case 'quad'
            sigmoidFuncHand = @quadraticSigmoid;
        otherwise
            error('Unexpected value for sigmoidStr encounterd, ''%s''', ...
                sigmoidStr);
    end
    sigmoidFcn = @(t, tStart, tEnd, phiStart, phiEnd) ...
        sigmoidFuncHand(t, tStart, tEnd, phiStart, phiEnd, ...
        sigmoidArgs{:});
else
    sigmoidFcn = [];
end


%% Build and Compose Piecewise Phase Function
numPieces = (numSteps * 2) - 1;
isModPiece = csmu.iseven(1:numPieces);

stepIds = zeros(1, numPieces);
stepIds(~isModPiece) = 1:numSteps;
stepIds(isModPiece) = NaN;
modIds = zeros(1, numPieces);
modIds(isModPiece) = 1:numMods;
modIds(~isModPiece) = NaN;

pieceEndTimes = repelem(timeSteps, 2);
pieceEndTimes = pieceEndTimes(1:(end - 1));
pieceEndTimes(isModPiece) = pieceEndTimes(isModPiece) + modTimes;
pieceStartTimes = [0, pieceEndTimes(1:(end - 1))];

if hasPO
    priorPhase = wrapTo2Pi(refPhaseFcn(pieceStartTimes(1)) + poSteps(1));
else
    priorPhase = 0;
end

if doDebug
    figure();
    if hasPO
        ax1 = subplot(3, 1, 1);
    else
        ax1 = subplot(2, 1, 1);
    end
    yticks(0:45:360);
    ylim([-10, 370]);
    ylabel('Phase (deg)');
    xlabel('Time (ms)');
    hold on
    grid on
    box off
    if hasPO
        ax2 = subplot(3, 1, 2);
        yticks(-180:45:180);
        ylim([-190, 190]);
        ylabel('Phase Offset (deg)');
        xlabel('Time (ms)')
        hold on
        grid on
        box off
    end
    if hasPO
        ax3 = subplot(3, 1, 3);
    else
        ax3 = subplot(2, 1, 2);
    end
    ax3.YScale = 'log';
    ylabel('Instant. Frequency (Hz)');
    xlabel('Time (ms)')
    hold on
    grid on
    box off
    if hasPO
        linkaxes([ax1, ax2, ax3], 'x');
    else
        linkaxes([ax1, ax3], 'x');
    end
end

phaseFcnPieces = cell(1, numPieces);
for iPiece = 1:numPieces   
    startTime = pieceStartTimes(iPiece);
    endTime = pieceEndTimes(iPiece);

    if isModPiece(iPiece)
        modId = modIds(iPiece);
        unqFMId = unqFreqModIds(modId);
        
        modTime = modTimes(modId);
                 
        fmPhaseFcn = @(t) fmFcnByUnqFreqMod{unqFMId}(t, startTime);
               
        if hasPO
            nextStepId = stepIds(iPiece + 1);
            targetPO = poSteps(nextStepId);
                                    
            if modTime == 0
                phaseModShift = wrapTo2Pi(targetPO ...
                    + refPhaseFcn(startTime) - fmPhaseFcn(startTime));
                phaseModFcn = @(t) phaseModShift;
                
                doContAdj = false;
            else
                startPO = wrapTo2Pi(priorPhase - refPhaseFcn(startTime));
                fmDPhase = ...
                    wrapTo2Pi(fmPhaseFcn(endTime) ...
                    - fmPhaseFcn(startTime));
                refDPhase = ...
                    wrapTo2Pi(refPhaseFcn(endTime) ...
                    - refPhaseFcn(startTime));
                ddPhase = fmDPhase - refDPhase;
                phaseModShift = wrapToPi(targetPO - startPO - ddPhase);
                phaseModFcn = @(t) sigmoidFcn(...
                    t, startTime, endTime, 0, phaseModShift);     

                doContAdj = true;
            end
            pieceBasePhaseFcn = @(t) fmPhaseFcn(t) + phaseModFcn(t);
        else
            pieceBasePhaseFcn = fmPhaseFcn;

            doContAdj = true;
        end        

        if doDebug
            % debug strings
            pieceIdStr = sprintf('Jump %2d', modId);   
        end
    else
        stepId = stepIds(iPiece);        
        stepFreq = freqSteps(stepId);
               
        pieceBasePhaseFcn = @(x) 2 * pi * stepFreq * x;

        doContAdj = true;

        if doDebug
            % debug strings
            pieceIdStr = sprintf('Step %2d', stepId);
        end
    end   

    if ~doContAdj
        contAdj = 0;
    else
        % enforce phase continuity    
        contAdj = wrapTo2Pi(priorPhase - pieceBasePhaseFcn(startTime));
    end
    piecePhaseFnc = @(t) wrapTo2Pi(pieceBasePhaseFcn(t) + contAdj);
    phaseFcnPieces{iPiece} = piecePhaseFnc;
    priorPhase = piecePhaseFnc(endTime);

    if ~doDebug
        continue
    end

    % debug plots

    if abs(startTime - endTime) < eps(startTime)
        % e.g. for modTime == 0
        testT = startTime;
        plotArgs = {'o', 'LineWidth', 3, 'MarkerSize', 10};
    else
        testT = linspace(startTime, endTime, 500);
        plotArgs = {'.', 'LineWidth', 3, 'MarkerSize', 10};
    end
    testPhi = piecePhaseFnc(testT);
    p = plot(ax1, testT * 1e3, rad2deg(testPhi), ...
        plotArgs{:}, 'DisplayName', sprintf('Piece %2d', iPiece));
    xlim(ax1, 'auto')
    if hasPO
        plot(ax2, testT * 1e3, ...
            rad2deg(wrapToPi(testPhi - refPhaseFcn(testT))), ...
            plotArgs{:}, 'Color', p.Color)
        xlim(ax2, 'auto')
    end
    instFreq = diff(unwrap(testPhi)) ./ diff(testT) / 2/ pi;
    plot(ax3, testT(2:end) * 1e3, instFreq, plotArgs{:}, 'Color', p.Color); 
    xlim(ax3, 'auto')
    ylim(ax3, 'auto')
end

% ensure that the function will return NaN for value < 0
composeFcnPieces = [{@(t) NaN}, phaseFcnPieces];
composeKeytimes = [(0 - eps(0)), pieceEndTimes];
phaseFcn = composePiecewiseFcn(composeFcnPieces, composeKeytimes);
end

function phi = logfSinFMPhaseRamp(t, t1, t2, f1, f2)
% function defining a sinusoidal ramp from f1 to f2 in log(f)-space
freqFcn = @(delt) f1 ...
    * exp((log(f2) - log(f1)) ...
    * ((sin((2 * pi * delt / 2 / (t2 - t1)) ...
    - (pi / 2)) / 2) + (1 / 2)));

% the following lines perform an integral of (2 * pi * freqFcn) from t1 to
% t for all values in the vector t

% rescale integral to be from zero to one
tDelta = t(:) - t1;
reScaledFcn = @(tt, ~) freqFcn(tt * tDelta) .* tDelta;

% perform integration with ode45 solver
[~, result] = ode45(reScaledFcn, [0, 0.5, 1], zeros(size(t)));
phi = 2 * pi * result(end, :)';
end