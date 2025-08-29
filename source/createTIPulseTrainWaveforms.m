function [waveforms, varargout] = createTIPulseTrainWaveforms( ...
    pulseDuration, rampDuration, pulseFreq, pulsesPerTrain, ...
    carrierFreq, interferenceBeatFreq, varargin)
%CREATETIPULSETRAINWAVEFORMS  Produces pulse trains for neural TI stim.
%
%   waveforms = createTIPulseTrainWaveforms(pulseDuration, rampDuration,
%   pulseFreq, pulsesPerTrain, carrierFreq, interferenceBeatFreq) Produces
%   a pair (2 column matrix) of stimulation waveforms of that consist of a
%   train of pulsesPerTrain pulses with each pulse lasting pulseDuration
%   seconds plus a ramp up and ramp down time as specified by rampDuration.
%   Each pulse modulates the each waveform to a frequency near carrierFreq
%   Hz to interfere to produce beats at interferenceBeatFreq Hz when the
%   waveforms are summed.
%
%   The ramp specified by rampDuration linearly ramps the stimulation 
%   amplitudes from 0 to amp1 and amp2 over rampT(1) seconds at the 
%   beginning of each pulse and from amp1 and amp2 to 0 over rampT(2) 
%   seconds at the end of each pulse. Where 
%   rampT = rampDuration(:)' .* [1, 1]. To specify different values for 
%   the ramp up time (1) and ramp down time (2), supply a 2-element vector 
%   for this parameter.
%
%   createTIPulseTrainWaveforms(__, 'NumTrains', numTrains)
%   {numTrains=1} Produces waveforms which repeats the stimulation 
%   numTrains time(s). Must be >= 1.
%
%   createTIPulseTrainWaveforms(__, 'InterTrainInterval', interTrainT)
%   {interTrainT=1} Defines the amount of time in seconds between
%   subsequent pulse trains. Has no effect if numTrains==1.
%
%   The total strict stimulation time is given by:
%       totalStimTime = ((((1 / pulseFreq) * pulsesPerTrain) * numTrains)
%                        + (interTrainT * (numTrains - 1))) 
%
%   createTIPulseTrainWaveforms(__, 'SamplingRate', sampRate) 
%   {sampRate=100e3} Defines the sampling rate for the waveforms in Hz.
%
%   createTIPulseTrainWaveforms(__, 'A1', amp1, 'A2', amp2) 
%   {amp1=1, amp2=1} Sets the amplitude of each waveform in arbitrairy 
%   units. This value can be supplied as a scalar or a vector, the value
%   for ampVecMode determine how vector values are handled and applied to
%   modulate each pulse train.
%       Scalar Example: 
%           (__, 'A1', 1, 'A2', 2, ...)
%       Vector Example with ampVecMode = 'exact':
%           (__, 'NumTrains', 3, 'A1', [1, 3, 1], 'A2', [3, 1, 3], 
%            'AmpVecMode', 'exact')
%       Vector Example with ampVecMode = 'repeat':
%           (__, 'NumTrains', 4, 'A1', [1, 2], 'A2', [2, 1, 0],
%            'AmpVecMode', 'repeat')
%           amp1 will be expanded to [1, 2, 1, 2]
%           amp2 will be expanded to [2, 1, 0, 2]
%
%   createTIPulseTrainWaveforms(__, 'AmpVecMode', ampVecMode) 
%   {ampVecMode='exact'} Determines the behavior when a vector is passed
%   for either amp1 or amp2 as follows:
%       - 'exact'  : If non-scalar vectors, amp1 and amp2 must have length 
%                    numTrains.
%       - 'repeat' : amp1 and amp2 can be vectors of any length and the 
%                    amplidude modulation pattern will be repeated
%                    (and/or truncated) until the last train is reached.
%
%   ampVecMode must be a valid member or initializer of the
%   UTILS.AMPLITUDEVECTORMODE enumeration class.
%
%   createTIPulseTrainWaveforms(__, 'ModulationSymetry', symetry)
%   {symetry='signal1'} symetry determines if:
%       - 'signal1'  : only the first waveform will be modulated while the
%                      second remains at the carrierFreq.
%                      signal 1 will be modulated to a higher frequency
%                      than the carrier frequency during the pulse.
%       - 'signal2'  : only the second waveform will be modulated. signal 2
%                      will be modulated to a higher frequency than than
%                      the carrier frequency during the pulse.
%       - 'symetric' : both waveforms are modulated symetrically to produce
%                      the desired summation waveform. The first waveform
%                      will be modulated to a frequency higher than the
%                      carrier frequency, while the second to a lower.
%   All values of symetry will produce the same summation waveform, however
%   the frequency and phase modulations to the two waveforms differently as
%   specified above. 
%   
%   symetry must be a valid member or initializer of the
%   UTILS.MODULATIONSYMETRY enumeration class.
%
%   createTIPulseTrainWaveforms(__, 'Flip', doFlip) {doFlip=false} If 
%   doFlip is true, then waveform 1 will be inverted (i.e. multiplied by 
%   -1 or shifted in phase 180 deg). This setting might be useful to 
%   correct for hardware setup polarity errors and/or implementation 
%   details.
%
%   createTIPulseTrainWaveforms(__, 'RampTime', rampT) {rampT=0.5} Linearly 
%   ramps the stimulation amplitudes from 0 to amp1 and amp2 over rampT 
%   seconds. To specify different values for the ramp up time (1) and ramp 
%   down time (2), supply a 2-element vector for this parameter.
% 
%   createTIPulseTrainWaveforms(__, 'WaitTime', waitSpec) {waitSpec=0.1}
%   Adds a period of no stimulation before waitT(1) and after waitT(2) the
%   stimulation period in seconds; where waitT = waitSpec(:)' * [1, 1]. To
%   specify different values for the pre-stimulation wait and the
%   post-stimulation wait, supply a 2-element vector for the waitSpec
%   parameter (similar to rampDuration).
%
%   The ramp occurs at the beginning and end of every pulse WITHIN the 
%   strict stimulation time, while the wait occurs at the beginning and end 
%   OUTSIDE of the strict stimualtin time.
%
%   The total waveform time is on the interval:
%       waveformDomain = [-waitT(1), totalStimTime + waitT(2)] 
%
%   createTIPulseTrainWaveforms(__, 'Control', doControl) {doControl=false} 
%   If doControl is true no stimulation will be performed and reference 
%   waveforms will be produced.
%
%   createTIPulseTrainWaveforms(__, 'Plot', doPlot) {doPlot=true} If doPlot 
%   is true a preview of the stimulation waveforms and their sum will be
%   displayed in a plot.
%
%   createTIPulseTrainWaveforms(__, 'OptimizePlot', doOptimizePlot) 
%   {doOptimizePlot=true} If doOptimizePlot is true then the plots, 
%   if generated, will be resampled to ~16x the max frequency or sampRate, 
%   whichever is lower.
%
%   createTIPulseTrainWaveforms(__, 'Debug', doDebug) {doDebug=false} will 
%   run the waveform creation code in debug mode if doDebug is true.
%
%   [waveforms, T] = createTIPulseTrainWaveforms(__) also returns the
%   time column vector T having the same number of rows as waveforms. The
%   stimulation begins at time 0s, so if waitT(1) > 0, there will be
%   negative time points.
%
%   [waveforms, T, parameters] = createTIPulseTrainWaveforms(__) also 
%   returns the parameters used for running the function as a struct.
%
%   createTIPulseTrainWaveforms(parameters) will run the creation code 
%   using the provided parameters struct.
%
%   See also CREATETISTIMWAVEFORMS, UTILS.GENERATEFMPHASEFCN, 
%       UTILS.MODULATIONSYMETRY, UTILS.AMPLITUDEVECTORMODE.

% Corban Swain, August 2025

VERSION = 'v0.4.1';

%% Input Handling
isNumericScalar = @(x) isnumeric(x) && isscalar(x);
isPositiveScalar = @(x) isNumericScalar(x) && (x > 0);
isPositiveIntegerScalar = @(x) isPositiveScalar(x) && utils.isint(x);
isNonNegativeScalar = @(x) isNumericScalar(x) && (x >= 0);
isLogicalScalar = @(x) islogical(x) && isscalar(x);
isNonNegativeVecWithLen = @(x, lengths) isnumeric(x) && (all(x >=0)) ...
    && isvector(x) && any(length(x) == lengths, 'all');
isNonNegativeVector = @(x) isnumeric(x) && (all(x >=0)) && isvector(x);

p = inputParser();
PULSE_DUR_KEY = 'PulseDuration';
p.addRequired(PULSE_DUR_KEY, isNonNegativeScalar); % seconds
RAMP_DUR_KEY = 'RampDuration';
p.addRequired(RAMP_DUR_KEY, @(x) isNonNegativeVecWithLen(x, [1, 2]));
PULSE_FQ_KEY = 'PulseFreq';
p.addRequired(PULSE_FQ_KEY, isPositiveScalar); % hertz
PPT_KEY = 'PulsesPerTrain';
p.addRequired(PPT_KEY, isPositiveIntegerScalar);
CARRIER_FQ_KEY = 'CarrierFreq';
p.addRequired(CARRIER_FQ_KEY, isPositiveScalar); % hertz
INTERF_BEAT_FREQ_KEY = 'InterferenceBeatFreq';
p.addRequired(INTERF_BEAT_FREQ_KEY, isNonNegativeScalar); % hertz

NUM_TRAIN_KEY = 'NumTrains';
p.addParameter(NUM_TRAIN_KEY, 1, isPositiveIntegerScalar);
ITI_KEY = 'InterTrainInterval';
p.addParameter(ITI_KEY, 1, isNonNegativeScalar); % seconds
SAMP_RATE_KEY = 'SamplingRate';
p.addParameter(SAMP_RATE_KEY, 100e3, isPositiveScalar); % hertz
A1_KEY = 'A1';
p.addParameter(A1_KEY, 1, isNonNegativeVector);
A2_KEY = 'A2';
p.addParameter(A2_KEY, 1, isNonNegativeVector);
AMP_VEC_MODE_KEY = 'AmpVecMode';
p.addParameter(AMP_VEC_MODE_KEY, 'exact', ...
    @utils.AmplitudeVectorMode.isMember);
MOD_SYM_KEY = 'ModulationSymetry';
p.addParameter(MOD_SYM_KEY, 'signal1', ...
    @utils.ModulationSymetry.isMember);
DO_FLIP_KEY = 'Flip';
p.addParameter(DO_FLIP_KEY, false, isLogicalScalar);
WAIT_TIME_KEY = 'WaitTime';
p.addParameter(WAIT_TIME_KEY, 0.5, ... % seconds
    @(x) isNonNegativeVecWithLen(x, [1, 2]));
DO_BRUTE_FORCE_KEY = 'auto';
p.addParameter(DO_BRUTE_FORCE_KEY, false, ...
    @(x) isLogicalScalar(x) ...
    || (utils.scalarStringLike(x) && strcmpi(x, 'auto')));
DO_CONTROL_KEY = 'Control';
p.addParameter(DO_CONTROL_KEY, false, isLogicalScalar);
DO_PLOT_KEY = 'Plot';
p.addParameter(DO_PLOT_KEY, true, isLogicalScalar);
DO_OPT_PLOT_KEY = 'OptimizePlot';
p.addParameter(DO_OPT_PLOT_KEY, true, isLogicalScalar);
DO_DEBUG_KEY = 'Debug';
p.addParameter(DO_DEBUG_KEY, false, isLogicalScalar);

% handle the situation if a single input argument is supplied 
% it must be a struct of parameter values
if nargin() == 1
    if ~isstruct(pulseDuration)
        error(['If a single argument is supplied, it must be a', ...
            ' parameter struct.']);
    end

    inputStruct = pulseDuration;
    fnames = fieldnames(inputStruct);
    reqParamNames = {
        PULSE_DUR_KEY
        RAMP_DUR_KEY 
        PULSE_FQ_KEY
        PPT_KEY
        CARRIER_FQ_KEY
        INTERF_BEAT_FREQ_KEY};
    numReqParams = length(reqParamNames);
    reqParams = cell(1, numReqParams);
    for iParam = 1:numReqParams
        paramName = reqParamNames{iParam};

        matchIdx = find(strcmpi(paramName, fnames));        

        if isempty(matchIdx)
            error(['The parameter "%s" is required; it must', ...
                ' be supplied if running this function with a' ...
                ' single argument parameter struct.'], paramName);
        elseif numel(matchIdx) > 1
            error(['Found multiple parameters with' ...
                ' case-insensitive name "%s", only one should be' ...
                ' passed.'], paramName);
        end

        actualParamName = fnames{matchIdx};
        reqParams{iParam} = inputStruct.(actualParamName);
        inputStruct = rmfield(inputStruct, actualParamName);
    end

    [pulseDuration, rampDuration, pulseFreq, ...
        pulsesPerTrain, carrierFreq, interferenceBeatFreq] ...
        = reqParams{:};
    varargin = {inputStruct};
end

p.parse(pulseDuration, rampDuration, pulseFreq, pulsesPerTrain, ...
    carrierFreq, interferenceBeatFreq, varargin{:});

nargoutchk(0, 3);

inputs = p.Results;
parameters = inputs;

rampT = rampDuration(:)' .* [1, 1];
parameters.(RAMP_DUR_KEY) = rampT;

numTrains = inputs.(NUM_TRAIN_KEY);
interTrainInterval = inputs.(ITI_KEY);
sampRate = inputs.(SAMP_RATE_KEY);
amp1 = inputs.(A1_KEY);
amp2 = inputs.(A2_KEY);
ampVecMode = utils.AmplitudeVectorMode(inputs.(AMP_VEC_MODE_KEY));

symetry = utils.ModulationSymetry(inputs.(MOD_SYM_KEY));
parameters.(MOD_SYM_KEY) = symetry;

doFlip = inputs.(DO_FLIP_KEY);

waitT = inputs.WaitTime(:)' .* [1, 1];
parameters.WaitTime = waitT;

doControl = inputs.(DO_CONTROL_KEY);
doPlot = inputs.(DO_PLOT_KEY);
doOptimizePlot = inputs.(DO_OPT_PLOT_KEY);
doDebug = inputs.(DO_DEBUG_KEY);

fmArgs = {};
if doDebug
    debugFM = true;
    for iArg = 1:length(fmArgs)
        % if debug has already been pass in the fmArgs array then dont add
        % a duplicate and potentially conflicting entry
        if strcmpi('debug', fmArgs{iArg})
            debugFM = false;
        end
    end

    % by default if debugging this func, also debug the fm func
    if debugFM
        fmArgs = [fmArgs, {'Debug', true}];
    end
end

if interferenceBeatFreq >= carrierFreq
    error('`%s` (%.3f) must be less than `%s` (%.3f).', ...
        INTERF_BEAT_FREQ_KEY, interferenceBeatFreq, ...
        CARRIER_FQ_KEY, carrierFreq);
end

if interferenceBeatFreq > (carrierFreq / 5)
    warning(['The interf. beat frequency (%.3f Hz) should be much' ...
        ' lower than the carrier frequency (%.3f Hz), consider' ...
        ' increasing the carrier frequency.'], interferenceBeatFreq, ...
        carrierFreq)
end

fullPulseDuration = sum(rampT, 'all') + pulseDuration;
pulsePeriod = 1 / pulseFreq;

if fullPulseDuration > pulsePeriod
    error(['The pulse + ramp duration ' ...
        '(%.3f + sum([%s], ''all'') = %.3f s) ' ...
        'must not be longer ' ...
        'than 1 / pulseFreq (pulsePeriod = %.3f s)'], ...
        pulseDuration, num2str(rampT), fullPulseDuration, pulsePeriod)
end

if ampVecMode == utils.AmplitudeVectorMode.EXACT
    % issue error if amp1 or amp2 vectors have improper length; valid
    % values must have a length of 1 (scalar) or numTrains
    if ~any(length(amp1) == [1, numTrains])
        error(['`%s` must be a scalar or a vector with length `%s`' ...
            ' (%d); instead recieved a vector of length %d. Update the' ...
            ' value of either `%s` or `%s` to fix.'], ...
            A1_KEY, NUM_TRAIN_KEY, numTrains, length(amp1), A1_KEY, ...
            NUM_TRAIN_KEY);
    end
    if ~any(length(amp2) == [1, numTrains])
        error(['`%s` must be a scalar or a vector with length `%s`' ...
            ' (%d); instead recieved a vector of length %d. Update the' ...
            ' value of either `%s` or `%s` to fix.'], ...
            A2_KEY, NUM_TRAIN_KEY, numTrains, length(amp2), A2_KEY, ...
            NUM_TRAIN_KEY);
    end
end


maxRecommendCarrierFreq = (sampRate / 2.3) - interferenceBeatFreq;
if carrierFreq > maxRecommendCarrierFreq
    warning(['With the given parameter configuration, the maximum' ...
        ' recommended carrier frequency is %.3f Hz. However %.3f Hz was'...
        ' passed for the carrier frequency. Consider decreasing the' ...
        ' carrier frequency, increasing the sampling rate, or' ...
        ' decreasing the interference beat frequency.'], ...
        maxRecommendCarrierFreq, carrierFreq);
end

%% Computed Values
% fullPulseDuration = sum(rampT, 'all') + pulseDuration;
% pulsePeriod = 1 / pulseFreq;
sampPeriod = 1 / sampRate;
breakT = pulsePeriod - fullPulseDuration;
noBreak = breakT < eps();
breakCarrierFreq = 0;
if noBreak
    modT = 0;
else
    modT = breakT / 3;
end
totalTrainDuration = pulsePeriod * pulsesPerTrain;

numSignals = 2;

numFreqStepsPerPulse = 2;
numStepsPerTrain = (pulsesPerTrain * numFreqStepsPerPulse) + 1;
cycleNumsHelper = [0, repelem(0:(pulsesPerTrain - 1), ...
    numFreqStepsPerPulse)];
stimTimeStepsHelper = [0, ...
    repmat([fullPulseDuration + modT, pulsePeriod], ...
    [1, pulsesPerTrain])];
stimTimeStepsHelper = stimTimeStepsHelper - modT;

stimTimeSteps = (breakT / 2) ...
    + (stimTimeStepsHelper + (cycleNumsHelper .* pulsePeriod));
stimTimeSteps(end) = stimTimeSteps(end) - (breakT / 2) + modT;   

refFreqSteps = [repmat([breakCarrierFreq, carrierFreq], ...
    [1, pulsesPerTrain]), breakCarrierFreq];   
pulseSel = 2:2:numStepsPerTrain;

% make sure amp1 and amp2 are row vectors
amp1 = amp1(:)';
amp2 = amp2(:)';

% handle 'repeat' case for ampVecMode; repeat and resize amplitude vectors
if ampVecMode == utils.AmplitudeVectorMode.REPEAT
    a1Reps = ceil(numTrains / length(amp1));
    a2Reps = ceil(numTrains / length(amp2));
    amp1 = repmat(amp1, 1, a1Reps);
    amp2 = repmat(amp2, 1, a2Reps);
    % trim vectore since original length may not have been a multiple of 
    % numTrains
    amp1 = amp1(:, 1:numTrains); 
    amp2 = amp2(:, 1:numTrains);
end

% If amp1 or amp2 are scalars in 'exact' mode, multiply out to vector;
% otherwise this step leavs amp1 and amp2 unchanged.
amp1 = ones([1, numTrains]) .* amp1;
amp2 = ones([1, numTrains]) .* amp2;

if ~doControl
    %% Define Frequency Modulation Steps
    stimFreqSteps = repmat(refFreqSteps, [numSignals, 1]);
    switch symetry
        case utils.ModulationSymetry.SIGNAL1
            stimFreqSteps(1, pulseSel) = ...
                stimFreqSteps(1, pulseSel) + interferenceBeatFreq;
        case utils.ModulationSymetry.SIGNAL2
            stimFreqSteps(2, pulseSel) = ...
                stimFreqSteps(2, pulseSel) + interferenceBeatFreq;
        case utils.ModulationSymetry.SYMETRIC
            stimFreqSteps(:, pulseSel) = ...
                stimFreqSteps(:, pulseSel)...
                + ([+1; -1] * interferenceBeatFreq / 2);
    end 
end

if ~doControl
    % during the break the signals must be perfectly out of phase;
    % that is 180 deg or pi rad of phase offset.
    breakPO = 0;

    pulsePO = pi;

    basePOSteps = [repmat([breakPO, pulsePO], [1, pulsesPerTrain]), ... 
        breakPO];
end

if ~doControl
    %% Define Frequency Modulation Steps
    stimFreqSteps = repmat(refFreqSteps, [numSignals, 1]);
    switch symetry
        case utils.ModulationSymetry.SIGNAL1
            stimFreqSteps(1, pulseSel) = ...
                stimFreqSteps(1, pulseSel) + interferenceBeatFreq;
        case utils.ModulationSymetry.SIGNAL2
            stimFreqSteps(2, pulseSel) = ...
                stimFreqSteps(2, pulseSel) + interferenceBeatFreq;
        case utils.ModulationSymetry.SYMETRIC
            stimFreqSteps(:, pulseSel) = ...
                stimFreqSteps(:, pulseSel)...
                + ([+1; -1] * interferenceBeatFreq / 2);
    end

    %% Define Phase Modulation Steps
    switch symetry
        case utils.ModulationSymetry.SIGNAL1
            % NaN phase offset for signal 2 since it will be the reference
            % phase. It's phase will be determind based on continuity.
            stimPOSteps = [basePOSteps; repelem(NaN, numStepsPerTrain)];
        case utils.ModulationSymetry.SIGNAL2
            % NaN phase offset for signal 1 since it will be the reference
            % phase. It's phase will be determind based on continuity.
            stimPOSteps = [repelem(NaN, numStepsPerTrain); basePOSteps];
        case utils.ModulationSymetry.SYMETRIC
            if noBreak || noPulse
                stimPOSteps = [basePOSteps; -basePOSteps] / 2;
            else
                halfPODelta = (pulsePODelta / 2);
                symetricPOStepCycle1 = [
                    (     breakPO / 2) + [0,      halfPODelta], ...
                    (-1 * breakPO / 2) + [0,      halfPODelta]];
                symetricPOStepCycle2 = [
                    (-1 * breakPO / 2) + [0, -1 * halfPODelta], ...
                    (     breakPO / 2) + [0, -1 * halfPODelta]];
                nReps = ceil(pulsesPerTrain/2) + 1;
                stimPOStepsLong = [
                    repmat(symetricPOStepCycle1, [1, nReps]);
                    repmat(symetricPOStepCycle2, [1, nReps])];
                stimPOSteps = stimPOStepsLong(:, 1:numStepsPerTrain);
            end
    end
end

%% Generate Phase Functions
refPhaseFcn = utils.generateFMPhaseFcn(...
    stimTimeSteps, refFreqSteps, modT, ...
    'PhaseOffsetSteps', zeros(1, numStepsPerTrain), ...
    'RefPhaseFcn', @(~) 0, ...
    fmArgs{:});

if doControl
    phaseFcns = {refPhaseFcn, @(t) wrapTo2Pi(refPhaseFcn(t) + pi)};
else
    phaseFcns = cell(1, numSignals);
    for iSignal = 1:numSignals
        if symetry == utils.ModulationSymetry.SIGNAL1 && iSignal == 2
            phaseFcns{iSignal} = refPhaseFcn;
            continue
        elseif symetry == utils.ModulationSymetry.SIGNAL2 && iSignal == 1
            phaseFcns{iSignal} = refPhaseFcn;
            continue
        end

        phaseFcns{iSignal} = utils.generateFMPhaseFcn(...
            stimTimeSteps, ...
            stimFreqSteps(iSignal, :), ...
            modT, ...
            'PhaseOffsetSteps', stimPOSteps(iSignal, :), ...
            'RefPhaseFcn', refPhaseFcn, ...
            fmArgs{:});
    end
end

%% Handle Ramping and Waiting
numRampStepsPerPulse = 4;
numRampStepsPerTrain = (numRampStepsPerPulse * pulsesPerTrain) + 1;
rampCycleNumsHelper = [0, repelem(0:(pulsesPerTrain - 1), ...
    numRampStepsPerPulse)];
rampTimeHelper = [0, repmat(...
    cumsum([rampT(1), pulseDuration, rampT(2), breakT]), ...
    [1, pulsesPerTrain])];
rampTimes = (breakT / 2) + rampTimeHelper ...
    + (rampCycleNumsHelper * pulsePeriod);
rampTimes(end) = rampTimes(end) - (breakT / 2);
rampStep = [0, repmat([1, 2, 3, 0], [1, pulsesPerTrain])];

rampFcns = cell(numRampStepsPerPulse * pulsesPerTrain, 1);
for iFcn = 1:numRampStepsPerTrain
    switch rampStep(iFcn)
        case 0
            stepFcn = @(~) 0;
        case 1
            stepFcn = @(t) (t - rampTimes(iFcn) + rampT(1)) ...
                / rampT(1);
        case 2
            stepFcn = @(~) 1;
        case 3
            stepFcn = @(t) -1 * (t - rampTimes(iFcn)) / rampT(2);
    end
    rampFcns{iFcn} = stepFcn;
end

ampFcn = utils.composePiecewiseFcn(rampFcns, rampTimes);

%% Generate Waveforms
timeArray = colon(0, sampPeriod, totalTrainDuration)';
numSamples = length(timeArray);
stimSelection = (timeArray > 0) & (timeArray < totalTrainDuration);
stimTimeArray = timeArray(stimSelection);
phaseforms = zeros(sum(stimSelection), numSignals);
waveforms = zeros(numSamples, numSignals);
for iSignal = 1:numSignals
    phaseforms(:, iSignal) = phaseFcns{iSignal}(stimTimeArray);
    waveforms(stimSelection, iSignal) = ...
        sin(phaseforms(:, iSignal)) .* ampFcn(stimTimeArray); 
end

%% Generate Sequential Pulse Trains
numInterTrainIntervals = numTrains - 1;
if numInterTrainIntervals > 0
    interTrainIntervalLength = ceil(interTrainInterval / sampPeriod);
    interTrainMatrix = zeros([interTrainIntervalLength, numSignals], ...
        'like', waveforms);
    interTrainPhaseMatrix = interTrainMatrix * nan();
    
    originalPhaseforms = phaseforms;
    phaseforms = zeros(numSamples, numSignals) * nan();
    phaseforms(stimSelection, :) = originalPhaseforms;
    phaseforms = cat(1, ...
        repmat(cat(1, phaseforms, interTrainPhaseMatrix), ...
        [numInterTrainIntervals, 1]), ...
        phaseforms);

    % create the final waveforms matrix folded with size
    % [numSamples, numSignals, numTrains]
    waveforms = repmat(waveforms, [1, 1, numTrains]);
    
    % create a matrix of the amplitudes with size 
    % [1, numSignals, numTrains]
    ampMatrix = cat(1, amp1, amp2); % original vectors are row vectors
    ampMatrix = permute(ampMatrix, [3, 1, 2]);
    
    % use element wise multiplication to scale waveforms to the amplitude
    waveforms = waveforms .* ampMatrix;

    if interTrainIntervalLength > 0
        % create the interval matrix folded with size
        % [interTrainIntervalLength, numSignals, numTrains]. There will be 
        % one extra interval at the end which we will remove later.
        allIntervalMatrix = repmat(interTrainMatrix, [1, 1, numTrains]);
    
        % update all train matrix by concatenating the intervals along the 
        % sample/time dimension (1)
        waveforms = cat(1, waveforms, allIntervalMatrix);
    end

    % reshape the 3D train matrix into a 2D waveform matrix by unfolding 
    % the train dimension (3)
    waveforms = permute(waveforms, [1, 3, 2]);
    waveforms = reshape(waveforms, [], numSignals);

    % drop the extra interval at the end, if needed
    if interTrainIntervalLength > 0
        startTrimIndex = size(waveforms, 1) + 1 ...
            - interTrainIntervalLength;
        waveforms(startTrimIndex:end, :) = [];
    end

    numSamples = size(waveforms, 1);
    timeArray = colon(0, (numSamples - 1))' * sampPeriod;
    stimTimeArray = timeArray;
else
    % In the case where there is only one train (numInterTrainIntervals 
    % is 0) we can know amp1 and amp2 are scalars and apply them to the two
    % columns of waveforms, respectively.
    waveforms = waveforms .* [amp1, amp2];
end

%% Add Wait Before and After Stimulus
preWaitLength = ceil(waitT(1) / sampPeriod);
if preWaitLength > 0
    preWaitMatrix = zeros([preWaitLength, numSignals], 'like', ...
        waveforms);

    waveforms = cat(1, preWaitMatrix, waveforms);
    timeSubArray = colon(-preWaitLength, -1)' * sampPeriod;
    timeArray = cat(1, timeSubArray, timeArray);
    numSamples = numSamples + preWaitLength;
end

postWaitLength = ceil(waitT(2) / sampPeriod);
if postWaitLength > 0
    postWaitMatrix = zeros([postWaitLength, numSignals], 'like', ...
        waveforms);

    waveforms = cat(1, waveforms, postWaitMatrix);
    timeSubArray = timeArray(end) ...
        + (colon(1, postWaitLength)' * sampPeriod);
    timeArray = cat(1, timeArray, timeSubArray);
    numSamples = numSamples + postWaitLength;
end

%% Display Plots
if doPlot
    if doControl
        freqs = refFreqSteps;
    else
        freqs = stimFreqSteps;
    end

    if doOptimizePlot
        maxFreq = max(freqs, [], 'all');
        nyquistBufferFreq = maxFreq * 16;
        if nyquistBufferFreq < sampRate
            rateRatio = round(sampRate / nyquistBufferFreq);
            newPeriod = sampPeriod * rateRatio;
            plotSampRate = 1 / newPeriod;
            plotWaveforms = waveforms(1:rateRatio:end, :);
            plotTimeArray = timeArray(1:rateRatio:end, :);
            plotPhaseforms = phaseforms(1:rateRatio:end, :);
            plotStimTimeArray = stimTimeArray(1:rateRatio:end, :);
        else
            % no optimization to perform
            plotWaveforms = waveforms;
            plotTimeArray = timeArray;
            plotStimTimeArray = stimTimeArray;
            plotPhaseforms = phaseforms;
            plotSampRate = sampRate;
        end
    end
    makePlots(...
        plotWaveforms, ...
        plotTimeArray, ...
        doDebug, ...
        doFlip, ...
        freqs, ...
        plotStimTimeArray, ...
        plotPhaseforms, ...
        totalTrainDuration, ...
        plotSampRate, ...
        amp1, ...
        amp2, ...
        parameters, ...
        VERSION)
end    

%% Apply Flipping
if doFlip
    waveforms = waveforms .* [-1, 1];  
end

%% Check Output
if any(isnan(waveforms), 'all')
    warning(['NaN value(s) found in output waveforms; an unexpected' ...
        ' error has occured.'])
end

%% Format Outputs
switch nargout
    case 2
        varargout = {timeArray};
    case 3
        varargout = {timeArray, parameters};
    otherwise
        varargout = {};
end
end

function makePlots(...
    waveforms, ...
    timeArray, ...
    doDebug, ...
    doFlip, ...
    freqs, ...
    stimTimeArray, ...
    phaseforms, ...
    totalTrainDuration, ...
    sampRate, ...
    amp1, ...
    amp2, ...
    parameters, ...
    VERSION)

    useSubplot = verLessThan('matlab', '9.7');
    compWaveform = sum(waveforms, 2);

    maxF = max(freqs, [], 'all');
    minF = min(freqs, [], 'all');

    if doDebug
        figure('Name', 'TI Pulse Train DEBUG Plots');
        if doFlip
            title('Phase Offset vs. Time Between (-1 * S1) and S2');
        else
            title('Phase Offset vs Time Between S1 and S2');
        end
        plot(stimTimeArray, rad2deg(wrapToPi(diff(phaseforms, 1, 2))), ...
            'k.', 'MarkerSize', 10);
        ylabel('Phase Offset (deg)');
        xlabel('Time (s)');
        ylim([-180, 180]);
        xlim([0, totalTrainDuration]);

        figure();
        
        spectArgs = {sampRate, 'spectrogram',...
            'TimeResolution', 1/ sampRate * 1000 ...
            'OverlapPercent', 90, ...
            'Leakage', 0.5, ...
            'FrequencyLimits', [minF / 5, maxF * 2.5]};

        if useSubplot
            ax1 = subplot(3, 1, 1);
        else
            t = tiledlayout('vertical');
            title(t, ...
                'Spectrograms for Individual and Composite Waveforms');
            t.TileSpacing = 'compact';
            t.Padding = 'compact';

            ax1 = nexttile();
        end
        try
            pspectrum(waveforms(:, 1), spectArgs{:});
        catch ME
            warning('pspectrum plot 1 failed.\n\t%s\n\t', ...
                ME.identifier, ME.message)
        end

        if useSubplot
            ax2 = subplot(3, 1, 1);
        else
            ax2 = nexttile();
        end
        try
            pspectrum(waveforms(:, 2), spectArgs{:});
        catch ME
            warning('pspectrum plot 2 failed.\n\t%s\n\t%s', ...
                ME.identifier, ME.message)
        end


        if useSubplot
            ax3 = subplot(3, 1, 1);
        else
            ax3 = nexttile();
        end
        try
            pspectrum(compWaveform, spectArgs{:});
        catch ME
            warning('pspectrum plot 3 failed.\n\t%s\n\t%s', ...
                ME.identifier, ME.message)
        end

        linkaxes([ax1, ax2, ax3], 'x'); 
    end
        
    phaseTol = maxF * 2 * pi / sampRate * 1.05;
    instantFreqFcn = @(x) diff(unwrap(x, (2 * pi) - phaseTol, 1), 1, 1) ...
        * sampRate / 2 / pi;
    instantFreq = instantFreqFcn(phaseforms);

    s1Color = [228, 26, 28] / 255;
    s2Color = [55, 126, 184] / 255;
    sumColor = [77, 175, 74] / 255;    

    f = figure('Name', 'TI Pulse Train Summary Plots');    
    plotArgs = {'-', 'LineWidth', 2};
    if ~useSubplot
        t = tiledlayout('vertical');
        title(t, 'TI Pulse Train Waveforms and Frequency Summary');
        t.TileSpacing = 'compact';
        t.Padding = 'compact';
    end

    if useSubplot
        ax1 = subplot(5, 1, 1);
    else
        ax1 = nexttile();
    end
    if doFlip
        title('FLIPPED Signal 1, (-1 * S1)');
    else
        title('Signal 1, S1');
    end
    hold(ax1, 'on');
    plot(timeArray, waveforms(:, 1), plotArgs{:}, 'Color', s1Color);
    ylabel('Amplitude (a.u.)');
    referenceAmp = max(cat(1, amp1, amp2), [], 'all');
    referenceYLim = (...
        [-1, 1] * utils.getScaleFactorFromAmp(referenceAmp) * 1.1);
    ylim(referenceYLim);
    xlim([timeArray(1), timeArray(end)]);
    box('off');
    grid('on');

    if useSubplot
        ax2 = subplot(5, 1, 2);
    else
        ax2 = nexttile();
    end
    title('Signal 2, S2');
    hold(ax2, 'on');
    plot(timeArray, waveforms(:, 2), plotArgs{:}, 'Color', s2Color);
    ylabel('Amplitude (a.u.)');
    ylim(referenceYLim);
    xlim([timeArray(1), timeArray(end)]);
    box('off');
    grid('on');

    if useSubplot
        ax3 = subplot(5, 1, 3);
    else
        ax3 = nexttile();
    end
    if doFlip
        title('Summed Signal, (-1 * S1) + S2');
    else
        title('Summed Signal, S1 + S2');
    end
    hold(ax3, 'on');
    plot(timeArray, compWaveform, plotArgs{:}, 'Color', sumColor);
    ylabel('Amplitude (a.u.)');
    maxSummedAmp = max(amp1 + amp2, [], 'all');
    ylim([-1, 1] * utils.getScaleFactorFromAmp(maxSummedAmp) * 1.1);
    xlim([timeArray(1), timeArray(end)]);
    box('off');
    grid('on');

    if useSubplot
        ax4 = subplot(5, 1, 4);
    else
        ax4 = nexttile();
    end
    title('Instantaneous Frequency');
    hold(ax4, 'on');
    plot(stimTimeArray(2:end), instantFreq(:, 1), plotArgs{:}, ...
        'Color', s1Color, ...
        'DisplayName', 'S1');
    hold on;
    plot(stimTimeArray(2:end), instantFreq(:, 2), plotArgs{:}, ...
        'Color', s2Color, ...
        'DisplayName', 'S2');
    ylabel('Frequency (Hz)');
    xlabel('Time (s)');    
    legend(ax4);    
    xlim([timeArray(1), timeArray(end)]); 
    
    allFreqs = [instantFreq(:)];
    maxIF = max(allFreqs, [], 'all');
    minIF = min(allFreqs, [], 'all');
    deltaLogF = log(maxIF) - log(minIF);
    if exp(deltaLogF) <= (1 + 1e-6)
        ylm = ([-1, 1] + [minIF, maxIF]);
    elseif isinf(deltaLogF)
        rnge = maxIF - minIF;
        ylm = [minIF, maxIF] + ([-1, 1] * rnge * 0.05);
    else
        lgylm = [log(minIF), log(maxIF)] + (deltaLogF * 0.1 * [-1, 1]);
        ylm = exp(lgylm); 
    end

    ylim(ylm);
    if minIF > 0
        ax4.YScale = 'log';
    end
    box('off');
    grid('on');

    if useSubplot
        ax5 = subplot(5, 1, 5);
    else
        ax5 = nexttile();
    end
    title(ax5, 'Parameters');
    hold(ax5, 'on');
    box(ax5, 'on');
    str = '';
    plot(0, NaN);
    yticks([]);
    xticks([]);
    paramNames = fieldnames(parameters);
    for iParam = 1:length(paramNames)
        paramName = paramNames{iParam};
        val = parameters.(paramName);
        if isempty(str)
            joinStr = '';
        else
            joinStr = ', ';
        end
        str = [str, joinStr, sprintf('%s=%s', paramName, ...
            utils.toStr(val))];
        % str = strrep(str, '{', '\{');
        % str = strrep(str, '}', '\}');
    end

    if verLessThan('matlab', '25.1')
        datetext = datestr(now());
    else
        datetext = datetime('now');
    end

    str = sprintf('%s (Generated on %s, Version %s)', str, ...
        datetext, VERSION);
    curPos = f.Position;
    newW = curPos(3) * 1.2;
    newH = curPos(4) * 2.2;
    f.Position = [curPos(1:2) - ([newW, newH] - curPos(3:4)), ...
        newW, newH];
         
    annotation('textbox', ax5.Position, 'String', str, 'LineStyle', ...
        'none', 'FontName', 'FixedWidth', ...
        'VerticalAlignment', 'middle', 'Interpreter', 'none');

    linkaxes([ax1, ax2, ax3, ax4], 'x');       
end