function [waveforms, varargout] = createTIStimWaveformsPulseTrain( ...
    pulseDuration, rampDuration, pulseFreq, pulsesPerTrain, ...
    carrierFreq, interferenceBeatFreq, varargin)
%CREATETISTIMWAVEFORMSPULSETRAIN  Produces pulse trains for neural TI stim.
%
%   waveforms = createTIStimWaveforms(duration, carrierFreq, pulseFreq)
%   Produces a pair of stimulation waveforms of duration seconds using
%   carrier frequency carrierFreq Hz and interferes to produce a pulse
%   frequency of pulseFreq Hz.
%
%   waveforms will have size = [(duration + sum(waitT .* [1, 1])) *
%   sampRate, 2]
%
%   createTIStimWaveforms(__, 'SamplingRate', sampRate) {sampRate=100e3}
%   Defines the sampling rate for the waveforms in Hz.
%
%   createTIStimWaveforms(__, 'A1', amp1, 'A2', amp2) {amp1=1, amp2=1}
%   Sets the amplitude of the wavforms in arbitrairy units.
%
%   createTIStimWaveforms(__, 'DutyCycle', dutyCycle) {dutyCycle=1}
%   Modulates the signals such that the stimulation only occurs for
%   dutyCycle * (1 / pulseFreq) seconds (i.e. dutyCycle fraction) of the
%   pulse period.
%
%   createTIStimWaveforms(__, 'BreakCarrierFreq', breakFreq) 
%   {breakFreq='same'} If dutyCycle < 1, the carrier frequency will be
%   modulated to breakFreq during the no stimulation period. If unset,
%   breakFreq will be set to match the carrierFreq.
%
%   createTIStimWaveforms(__, 'FMTime', modT, 'FMArgs', fmArgs)
%   {modT='auto', fmArgs=cell()} modT sets the modulation time during which
%   the wavforms will smoothly transition from an initial frequency and
%   relative phase to a new frequency and relative phase. By default
%   modT=((1/pulseFreq)*dutyCycle*0.2), meaning 20% of the pulse time is
%   used to modulate into the pulse and another 20% of the pulse time is
%   used to modulate out of the pulse. The frequency modulation occurs only
%   after the pulse is set to begin and completes before the pulse is set
%   to end. The value for modT is limited to the interval [0,
%   ((1/pulseFreq)*dutyCycle/2)]; values greater than the upper limit will
%   be clipped and a warning will be raised.
% 
%   No modulation is performd and modT is ignored if dutyCycle==1.
%
%   Parameters for controlling the generation of the modulated signal
%   (UTILS.GENERATEFMPHASEFCN) can be set using fmArgs.
%
%   createTIStimWaveforms(__, 'ModulationSymetry', symetry)
%   {symetry='signal1'} symetry determines if:
%       - 'signal1'  : only the first waveform will be modulated while the
%                      second remains at the carrierFreq (or breakFreq).
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
%   createTIStimWaveforms(__, 'Flip', doFlip) {doFlip=false} If doFlip is
%   true, then waveform 1 will be inverted (i.e. multiplied by -1 or
%   shifted in phase 180 deg). This setting might be useful to correct for
%   hardware setup errors and/or implementation details.
%
%   createTIStimWaveforms(__, 'RampTime', rampT) {rampT=0.5} Linearly ramps
%   the stimulation amplitudes from 0 to amp1 and amp2 over rampT seconds.
%   To specify different values for the ramp up time (1) and ramp down time
%   (2), supply a 2-element vector for this parameter.
% 
%   createTIStimWaveforms(__, 'WaitTime', waitT) {waitT=0.1} Adds a period
%   of no stimulation before (1) and after (2) the stimulation period in
%   seconds. To specify different values for the pre-stimulation wait and
%   the post-stimulation wait, supply a 2-element vector for this
%   parameter.
%
%   The ramp occurs at the beginning and end WITHIN the stimulation time,
%   the wait occurs at the beginning and end OUTSIDE of the stimualtin
%   time. The total waveform time will be duration + sum([1, 1] .* waitT).
%
%   createTIStimWaveforms(__, 'Control', doControl) {doControl=false} If
%   doControl is true no stimulation will be performed, only unmodulated
%   waveformes at breakFreq will be applied.
%
%   createTIStimWaveforms(__, 'Plot', doPlot) {doPlot=true} If doPlot is
%   true a preview of the stimulation waveforms and their sum will be
%   displayed on a plot.
%
%   createTIStimWaveforms(__, 'Debug', doDebug) {doDebug=false} will run
%   the waveform creation code in debug mode if doDebug is true.
%
%   [waveforms, T] = createTIStimWaveforms(__) also returns the
%   time column vector T having the same number of rows as waveforms. The
%   stimulation begins at time 0s, so if waitT(1) > 0, there will be
%   negative time points.
%
%   [waveforms, T, parameters] = createTIStimWaveforms(__) also returns the
%   parameters used for running the function as a struct.
%
%   createTIStimWaveforms(parameters) will run the creation code using the
%   provided parameters struct.
%
%   See also UTILS.GENERATEFMPHASEFCN, UTILS.MODULATIONSYMETRY.

% Corban Swain, June 2023

VERSION = 'v0.3.0';

%% Input Handling
isNumericScalar = @(x) isnumeric(x) && isscalar(x);
isPositiveScalar = @(x) isNumericScalar(x) && (x > 0);
isPositiveIntegerScalar = @(x) isPositiveScalar(x) && utils.isint(x);
isNonNegativeScalar = @(x) isNumericScalar(x) && (x >= 0);
isLogicalScalar = @(x) islogical(x) && isscalar(x);
isNonNegativeVecWithLen = @(x, lengths) isnumeric(x) && (all(x >=0)) ...
    && isvector(x) && any(length(x) == lengths, 'all');

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
p.addRequired(INTERF_BEAT_FREQ_KEY, isPositiveScalar);

NUM_TRAIN_KEY = 'NumTrains';
p.addParameter(NUM_TRAIN_KEY, 1, isPositiveIntegerScalar);
ITI_KEY = 'InterTrainInterval';
p.addParameter(ITI_KEY, 1, isNonNegativeScalar); % seconds
SAMP_RATE_KEY = 'SamplingRate';
p.addParameter(SAMP_RATE_KEY, 100e3, isPositiveScalar); % hertz
A1_KEY = 'A1';
p.addParameter(A1_KEY, 1, isNonNegativeScalar);
A2_KEY = 'A2';
p.addParameter(A2_KEY, 1, isNonNegativeScalar);
MOD_SYM_KEY = 'ModulationSymetry';
p.addParameter(MOD_SYM_KEY, 'signal1', ...
    @utils.ModulationSymetry.isMember);
DO_FLIP_KEY = 'Flip';
p.addParameter(DO_FLIP_KEY, false, isLogicalScalar);
WAIT_TIME_KEY = 'WaitTime';
p.addParameter(WAIT_TIME_KEY, 0.5, ... % seconds
    @(x) isNonNegativeVecWithLen(x, [1, 2]));
DO_CONTROL_KEY = 'Control';
p.addParameter(DO_CONTROL_KEY, false, isLogicalScalar);
DO_PLOT_KEY = 'Plot';
p.addParameter(DO_PLOT_KEY, true, isLogicalScalar);
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

symetry = utils.ModulationSymetry(inputs.(MOD_SYM_KEY));
parameters.(MOD_SYM_KEY) = symetry;

doFlip = inputs.(DO_FLIP_KEY);

waitT = inputs.WaitTime(:)' .* [1, 1];
parameters.WaitTime = waitT;

doControl = inputs.(DO_CONTROL_KEY);
doPlot = inputs.(DO_PLOT_KEY);
doDebug = inputs.(DO_DEBUG_KEY);

if interferenceBeatFreq >= carrierFreq
    error('`%s` (%.1f) must be less than `%s` (%.1f).', ...
        INTERF_BEAT_FREQ_KEY, interferenceBeatFreq, ...
        CARRIER_FQ_KEY, carrierFreq);
end


if interferenceBeatFreq > (carrierFreq / 5)
    warning(['The interf. beat frequency (%.1f Hz) should be much' ...
        ' lower than the carrier frequency (%.1f Hz), consider' ...
        ' increasing the carrier frequency.'], interferenceBeatFreq, ...
        carrierFreq)
end

fullPulseDuration = sum(rampT, 'all') + pulseDuration;
pulsePeriod = 1 / pulseFreq;

if fullPulseDuration > pulsePeriod
    error(['The pulse + ramp duration ' ...
        '(%.1f + sum([%s], ''all'') = %.1f s) ' ...
        'must not be longer ' ...
        'than 1 / pulseFreq (pulsePeriod = %.1f s)'], ...
        pulseDuration, rampT, fullPulseDuration, pulsePeriod)
end

interfCyclePeriod = 1 / interferenceBeatFreq;


%% Computed Values
% fullPulseDuration = sum(rampT, 'all') + pulseDuration;
% pulsePeriod = 1 / pulseFreq;
sampPeriod = 1 / sampRate;
breakT = pulsePeriod - fullPulseDuration;
% noBreak = breakT < eps();
% noPulse = fullPulseDuration < eps();
breakCarrierFreq = 0;
modT = 0;
totalTrainDuration = pulsePeriod * pulsesPerTrain;

numSignals = 2;

numFreqStepsPerPulse = 2;
numStepsPerTrain = (pulsesPerTrain * numFreqStepsPerPulse) + 1;
cycleNumsHelper = [0, repelem(0:(pulsesPerTrain - 1), ...
    numFreqStepsPerPulse)];
stimTimeStepsHelper = [0, ...
    repmat([fullPulseDuration, pulsePeriod], [1, pulsesPerTrain])];

stimTimeSteps = (breakT / 2) ...
    + (stimTimeStepsHelper + (cycleNumsHelper .* pulsePeriod));
stimTimeSteps(end) = stimTimeSteps(end) - (breakT / 2);
% cycleNums = repelem(1:numPulses, 2);    

refFreqSteps = [repmat([breakCarrierFreq, carrierFreq], ...
    [1, pulsesPerTrain]), breakCarrierFreq];   
pulseSel = 2:2:numStepsPerTrain;

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

%% Generate Phase Functions
refPhaseFcn = utils.generateFMPhaseFcn(...
    stimTimeSteps, refFreqSteps, modT);

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
            'RefPhaseFcn', refPhaseFcn);
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
waveforms = waveforms .* [amp1, amp2];


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

    waveforms = cat(1, ...
        repmat(cat(1, waveforms, interTrainMatrix), ...
        [numInterTrainIntervals, 1]), ...
        waveforms);

    numSamples = size(waveforms, 1);
    timeArray = colon(0, (numSamples - 1))' * sampPeriod;
    stimTimeArray = timeArray;
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
    useSubplot = verLessThan('matlab', '9.7');
    compWaveform = sum(waveforms, 2);

    if doControl
        freqs = refFreqSteps;
    else
        freqs = stimFreqSteps;
    end
    maxF = max(freqs, [], 'all');
    minF = min(freqs, [], 'all');

    if doDebug
        figure();
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
            'TimeResolution', sampPeriod * 1000 ...
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
        pspectrum(waveforms(:, 1), spectArgs{:});

        if useSubplot
            ax2 = subplot(3, 1, 1);
        else
            ax2 = nexttile();
        end
        pspectrum(waveforms(:, 2), spectArgs{:});


        if useSubplot
            ax3 = subplot(3, 1, 1);
        else
            ax3 = nexttile();
        end
        pspectrum(compWaveform, spectArgs{:});

        linkaxes([ax1, ax2, ax3], 'x'); 
    end
        
    phaseTol = maxF * 2 * pi * sampPeriod * 1.05;
    instantFreqFcn = @(x) diff(unwrap(x, (2 * pi) - phaseTol, 1), 1, 1) ...
        / sampPeriod / 2 / pi;
    instantFreq = instantFreqFcn(phaseforms);

    s1Color = [228, 26, 28] / 255;
    s2Color = [55, 126, 184] / 255;
    sumColor = [77, 175, 74] / 255;    

    f = figure('Name', 'TI Stimulation Summary Plots');    
    plotArgs = {'-', 'LineWidth', 2};
    if ~useSubplot
        t = tiledlayout('vertical');
        title(t, 'TI Stimulation Waveforms and Frequency Summary');
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
    ylim([-1, 1] * amp1 * 1.1);
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
    ylim([-1, 1] * amp2 * 1.1);
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
    ylim([-1, 1] * (amp1 + amp2) * 1.1);
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
    str = sprintf('%s (Generated on %s, Version %s)', str, ...
        datestr(now()), VERSION);
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

%% Apply Flipping
if doFlip
    waveforms = waveforms .* [-1, 1];
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


