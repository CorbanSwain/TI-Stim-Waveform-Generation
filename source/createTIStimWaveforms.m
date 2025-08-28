function [waveforms, varargout] = createTIStimWaveforms(duration, ...
    carrierFreq, pulseFreq, varargin)
%CREATETISTIMWAVEFORMS  Produces waveforms for neural TI stimulation.
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

VERSION = 'v0.2.2';

%% Input Handling
isNumericScalar = @(x) isnumeric(x) && isscalar(x);
isPositiveScalar = @(x) isNumericScalar(x) && (x > 0);
isNonNegativeScalar = @(x) isNumericScalar(x) && (x >= 0);
isLogicalScalar = @(x) islogical(x) && isscalar(x);
isNonNegativeVecWithLen = @(x, lengths) isnumeric(x) && (all(x >=0)) ...
    && isvector(x) && any(length(x) == lengths, 'all');

p = inputParser();
p.addRequired('Duration', isNonNegativeScalar);
p.addRequired('CarrierFreq', isPositiveScalar);
% TODO - refactor function to use interferenceBeatFreq rather than
%   pulseFreq ... simply renaming and updating the documentation and
%   function calls.
p.addRequired('PulseFreq', isNonNegativeScalar);
p.addParameter('SamplingRate', 100e3, isPositiveScalar);
p.addParameter('DutyCycle', 1, @(x) isNonNegativeScalar(x) && (x <= 1));
p.addParameter('A1', 1, isNonNegativeScalar);
p.addParameter('A2', 1, isNonNegativeScalar);
p.addParameter('BreakCarrierFreq', [], ...
    @(x) (utils.scalarStringLike(x) && isequal(lower(char(x)), 'same')) ...
    || isNonNegativeScalar(x) || isempty(x));
p.addParameter('FMTime', [], ...
    @(x) (utils.scalarStringLike(x) && isequal(lower(char(x)), 'auto')) ...
    || isNonNegativeScalar(x) || isempty(x));
p.addParameter('FMArgs', {}, ...
    @(x) iscell(x) && (isvector(x) || isempty(x)));
p.addParameter('ModulationSymetry', 'signal1', ...
    @utils.ModulationSymetry.isMember);
p.addParameter('Flip', false, isLogicalScalar);
p.addParameter('RampTime', 0.5, @(x) isNonNegativeVecWithLen(x, [1, 2]));
p.addParameter('WaitTime', 0.1, @(x) isNonNegativeVecWithLen(x, [1, 2]));
p.addParameter('Control', false, isLogicalScalar);
p.addParameter('Plot', true, isLogicalScalar);
p.addParameter('Debug', false, isLogicalScalar);

% handle the situation if a single input argument is supplied 
% it must be a struct of parameter values
if nargin() == 1
    if ~isstruct(duration)
        error(['If a single argument is supplied, it must be a', ...
            ' parameter struct.']);
    end

    inputStruct = duration;
    fnames = fieldnames(inputStruct);
    reqParamNames = {'Duration', 'CarrierFreq', 'PulseFreq'};
    numReqParams = length(reqParamNames);
    reqParams = cell(1, numReqParams);
    for iParam = 1:numReqParams
        paramName = reqParamNames{iParam};

        matchIdx = find(strcmpi(paramName, fnames));        

        if isempty(matchIdx)
            error(['The parameter "%s" is required; it must', ...
                ' be supplied if running this function with a single', ...
                ' argument parameter struct.'], paramName);
        elseif numel(matchIdx) > 1
            error(['Found multiple parameters with case-insensitive', ...
                ' name "%s", only one should be passed.'], paramName);
        end

        actualParamName = fnames{matchIdx};
        reqParams{iParam} = inputStruct.(actualParamName);
        inputStruct = rmfield(inputStruct, actualParamName);
    end

    [duration, carrierFreq, pulseFreq] = reqParams{:};
    varargin = {inputStruct};
end

p.parse(duration, carrierFreq, pulseFreq, varargin{:});

nargoutchk(0, 3);

inputs = p.Results;
parameters = inputs;
sampRate = inputs.SamplingRate;
dutyCycle = inputs.DutyCycle;
amp1 = inputs.A1;
amp2 = inputs.A2;
breakCFreq = inputs.BreakCarrierFreq;
modT = inputs.FMTime;
fmArgs = inputs.FMArgs;
doFlip = inputs.Flip;
doControl = inputs.Control;
doPlot = inputs.Plot;
doDebug = inputs.Debug;

symetry = utils.ModulationSymetry(inputs.ModulationSymetry);
parameters.ModulationSymetry = symetry;

rampT = inputs.RampTime(:)' .* [1, 1];
parameters.RampTime = rampT;

waitT = inputs.WaitTime(:)' .* [1, 1];
parameters.WaitTime = waitT;

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
        parameters.FMArgs = fmArgs;
    end
end

if pulseFreq >= carrierFreq
    error('pulseFreq must be less than carrierFreq.')
end

if pulseFreq > (carrierFreq / 5)
    warning(['The pulse frequency (%.1f Hz) should be much lower', ...
        ' than the carrier frequency (%.1f Hz), consider increasing', ... 
        ' the carrier frequency.'], pulseFreq, carrierFreq)
end

cyclePeriod = 1 / pulseFreq;
pulseT = cyclePeriod * dutyCycle;

% set computed default values for breakCFreq, modT
doAutoSet = isempty(breakCFreq);
if ~doAutoSet && utils.scalarStringLike(breakCFreq)
    if strcmpi(breakCFreq, 'same')
        doAutoSet = true;
    else
        error(['Encountered unexpected value for breakFreq, ''%s'';', ...
            ' should be ''same'' or a positve scalar value.'], ...
            breakCFreq)
    end
end
if doAutoSet
    breakCFreq = carrierFreq;
    parameters.BreakCarrierFreq = breakCFreq;
end

doAutoSet = isempty(modT);
if ~doAutoSet && utils.scalarStringLike(modT)
    if strcmpi(modT, 'auto')
       doAutoSet = true;
    else
        error(['Encountered unexpected value for modT, ''%s'';',...
            ' should be ''auto'' or a non-negative scalar value.'], modT);
    end
end
if doAutoSet
    modT = pulseT * 0.2;
    parameters.FMTime = modT;
end

if modT > (pulseT / 2)
    warning(['The modulation time cannot be more than half the pulse', ...
        ' period (pulsePeriod/2=%f s); clipping the modulation time', ...
        ' to this max value.'], pulseT/2);
    modT = pulseT / 2;
    parameters.FMTime = modT;
end

totalRampTime = sum(rampT);
if totalRampTime > duration
    warning(['The total ramp time (sum(rampT(:)'' .* [1, 1])) cannot', ... 
        ' be greater than the stimulus duration. Scalling total ramp', ...
        ' time to equal the stimulus duration (%f s).'], duration)
    rampT = (rampT / totalRampTime) * duration;
    parameters.RampTime = rampT;
end

%% Computed Values
% cyclePeriod = 1 / pulseFreq;
% pulseT = cyclePeriod * dutyCycle;
sampPeriod = 1 / sampRate;
breakT = cyclePeriod - pulseT;
pulseCFDelta = (1 / pulseT);
noBreak = dutyCycle == 1;
noPulse = dutyCycle == 0;
numSignals = 2;

if noBreak || noPulse
    stimTimeSteps = duration;
    numSteps = 1;
    refFreqSteps = carrierFreq;  

    if noBreak 
        pulseSel = 1;
    else 
        % noPulse == true
        pulseSel = [];
    end

    basePOSteps = pi;    
else    
    numPulses = ceil(duration / cyclePeriod) + 1;
    numSteps = numPulses * 2;
    cycleNumsHelper = [0, repelem(0:(numPulses - 1), 2)];
    stimTimeStepsHelper = [0, ...
        repmat([pulseT - modT, cyclePeriod], [1, numPulses])];

    stimTimeSteps = (breakT / 2) ...
        + (stimTimeStepsHelper + (cycleNumsHelper .* cyclePeriod));
    stimTimeSteps = stimTimeSteps(1:(end - 1));
    % cycleNums = repelem(1:numPulses, 2);    

    refFreqSteps = repmat([breakCFreq, carrierFreq], [1, numPulses]);   
    pulseSel = 2:2:numSteps;

    if ~doControl
        % during the break the signals must be perfectly out of phase;
        % that is 180 deg or pi rad of phase offset.
        breakPO = pi;
        
        % pulsePODelta is chosen such that at the end of the frequency
        % modulation the modulated signal will be at the phase offset
        % location corresponding to where the phase offset would be for a
        % pulse beginning at pi phase offset with no frequency modulation
        % (i.e. an abrupt frequency and phase step). This allows the "duty
        % cycle" modified pulse to mimic a non duty cycle modified pulse
        % (given the same pulse-carrierFrequency difference, which--as an
        % aside--is not the same as having the same pulse frequency).
        pulsePODelta = wrapToPi(pulseCFDelta * 2 * pi * modT);
        
        pulsePO = wrapToPi(pulsePODelta + breakPO);

        basePOSteps = repmat([breakPO, pulsePO], [1, numPulses]);
    end
end

if ~doControl
    %% Define Frequence Modulation Steps
    stimFreqSteps = repmat(refFreqSteps, [numSignals, 1]);
    switch symetry
        case utils.ModulationSymetry.SIGNAL1
            stimFreqSteps(1, pulseSel) = ...
                stimFreqSteps(1, pulseSel) + pulseCFDelta;
        case utils.ModulationSymetry.SIGNAL2
            stimFreqSteps(2, pulseSel) = ...
                stimFreqSteps(2, pulseSel) + pulseCFDelta;
        case utils.ModulationSymetry.SYMETRIC
            stimFreqSteps(:, pulseSel) = ...
                stimFreqSteps(:, pulseSel) + ([+1; -1] * pulseCFDelta / 2);
    end

    %% Define Phase Modulation Steps
    switch symetry
        case utils.ModulationSymetry.SIGNAL1
            % NaN phase offset for signal 2 since it will be the reference
            % phase. It's phase will be determind based on continuity.
            stimPOSteps = [basePOSteps; repelem(NaN, numSteps)];
        case utils.ModulationSymetry.SIGNAL2
            % NaN phase offset for signal 1 since it will be the reference
            % phase. It's phase will be determind based on continuity.
            stimPOSteps = [repelem(NaN, numSteps); basePOSteps];
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
                nReps = ceil(numPulses/2) + 1;
                stimPOStepsLong = [
                    repmat(symetricPOStepCycle1, [1, nReps]);
                    repmat(symetricPOStepCycle2, [1, nReps])];
                stimPOSteps = stimPOStepsLong(:, 1:numSteps);
            end
    end
end

%% Generate Phase Functions
refPhaseFcn = utils.generateFMPhaseFcn(stimTimeSteps, refFreqSteps, modT, ...
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
rampWaitTimes = [rampT(1), duration - rampT(2), duration];
rampWaitFcns = {
    @(t) t / rampT(1)
    @(~) 1
    @(t) -1 * (t - duration) / rampT(2)};
ampFcn = utils.composePiecewiseFcn(rampWaitFcns, rampWaitTimes);

%% Generate Waveforms
timeArray = colon(-waitT(1), sampPeriod, duration + waitT(2))';
numSamples = length(timeArray);
stimSelection = (timeArray > 0) & (timeArray < duration);
stimTimeArray = timeArray(stimSelection);
phaseforms = zeros(sum(stimSelection), numSignals);
waveforms = zeros(numSamples, numSignals);
for iSignal = 1:numSignals
    phaseforms(:, iSignal) = phaseFcns{iSignal}(stimTimeArray);
    waveforms(stimSelection, iSignal) = ...
        sin(phaseforms(:, iSignal)) .* ampFcn(stimTimeArray); 
end
waveforms = waveforms .* [amp1, amp2];

%% Display Plots
% TDOD - code duplication between this and the pulse train script consider
%    consilidating into a utility function
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
        xlim([0, duration]);

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
    ylim([-1, 1] * utils.getScaleFactorFromAmp(amp1) * 1.1);
    xlim([-waitT(1), duration + waitT(2)]);
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
    ylim([-1, 1] * utils.getScaleFactorFromAmp(amp2) * 1.1);
    xlim([-waitT(1), duration + waitT(2)]);
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
    ylim([-1, 1] * utils.getScaleFactorFromAmp(amp1 + amp2) * 1.1);
    xlim([-waitT(1), duration + waitT(2)]);
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
    xlim([-waitT(1), duration + waitT(2)]); 
    
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
    f.Position = [curPos(1:2) - ([newW, newH] - curPos(3:4)), newW, newH];
         
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


