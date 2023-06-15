function [waveforms, varargout] = createTIStimWaveforms(duration, ...
    carrierFreq, pulseFreq, varargin)
%CREATETISTIMWAVEFORMS  Produces waveforms for neural TI stimulation.
%
%   waveforms = createTIStimWaveforms(duration, carrierFreq, pulseFreq)
%   Produces a pair of stimulation waveforms of duration seconds using
%   carrier frequency carrierFreq Hz and interferes to produce a pulse
%   frequency of pulseFreq Hz.
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
%   (GENERATEFMPHASEFCN) can be set using fmArgs.
%
%   createTIStimWaveforms(__, 'ModulationSymetry', symetry)
%   {symetry='signal1'} symetry determines if:
%       - 'signal1'  : only the first waveform will be modulated while the
%                      second remains at the carrierFreq (or breakFreq).
%       - 'signal2'  : only the second waveform will be modulated.
%       - 'symetric' : both waveforms are modulated symetrically to produce
%                      the desired summation waveform.
%   All values of symetry will produce the same summation waveform, however
%   the frequency and phase modulations to the two waveforms differently as
%   specified above. 
%   
%   symetry must be a valid member or initializer of the ModulationSymetry
%   enumeration class.
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
%   [waveforms, parameters] = createTIStimWaveforms(__) also returns the
%   parameters used for running the function as a struct.
%
%   createTIStimWaveforms(parameters) will run the creation code using the
%   provided parameters struct.
%
%   See also GENERATEFMPHASEFCN, MODULATIONSYMETRY.

% Corban Swain, June 2023

%% Input Handling
isNumericScalar = @(x) isnumeric(x) && isscalar(x);
isPositiveScalar = @(x) isNumericScalar(x) && (x > 0);
isNonNegativeScalar = @(x) isNumericScalar(x) && (x >= 0);
isLogicalScalar = @(x) islogical(x) && isscalar(x);

p = inputParser();
p.addRequired('duration', isNonNegativeScalar);
p.addRequired('carrierFreq', isPositiveScalar);
p.addRequired('pulseFreq', isPositiveScalar);
p.addParameter('SamplingRate', 100e3, isPositiveScalar);
p.addParameter('DutyCycle', 1, @(x) isNonNegativeScalar(x) && (x <= 1));
p.addParameter('A1', 1, isNonNegativeScalar);
p.addParameter('A2', 1, isNonNegativeScalar);
p.addParameter('BreakCarrierFreq', [], ...
    @(x) (csmu.scalarStringLike(x) && isequal(lower(char(x)), 'same')) ...
    || isPositiveScalar(x) || isempty(x));
p.addParameter('FMTime', [], ...
    @(x) (csmu.scalarStringLike(x) && isequal(lower(char(x)), 'auto')) ...
    || isNonNegativeScalar(x) || isempty(x));
p.addParameter('FMArgs', {}, @(x) iscell(x) && isvector(x));
p.addParameter('ModulationSymetry', 'signal1', ...
    @ModulationSymetry.isMember);
p.addParameter('RampTime', 0.5, isNonNegativeScalar);
p.addParameter('WaitTime', 0.1, isNonNegativeScalar);
p.addParameter('Control', false, isLogicalScalar);
p.addParameter('Plot', true, isLogicalScalar);
p.addParameter('Debug', false, isLogicalScalar);

if nargin() == 1
    p.parse(duration);
else
    p.parse(duration, carrierFreq, pulseFreq, varargin{:});
end

inputs = p.Results;
sampRate = inputs.SamplingRate;
dutyCycle = inputs.DutyCycle;
amp1 = inputs.A1;
amp2 = inputs.A2;
breakCFreq = inputs.BreakCarrierFreq;
modT = inputs.FMTime;
fmArgs = inputs.FMArgs;
symetry = ModulationSymetry(inputs.ModulationSymetry);
rampT = inputs.RampTime;
waitT = inputs.WaitTime;
doControl = inputs.Control;
doPlot = inputs.Plot;
doDebug = inputs.Debug;

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
if ~doAutoSet && csmu.scalarStringLike(breakCFreq)
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
end

doAutoSet = isempty(modT);
if ~doAutoSet && csmu.scalarStringLike(modT)
    if strcmpi(modT, 'auto')
       doAutoSet = true;
    else
        error(['Encountered unexpected value for modT, ''%s'';',...
            ' should be ''auto'' or a non-negative scalar value.'], modT);
    end
end
if doAutoSet
    modT = pulseT * 0.2;
end

if modT > (pulseT / 2)
    warning(['The modulation time cannot be more than half the pulse', ...
        ' period (pulsePeriod/2=%f s); clipping the modulation time', ...
        ' to this max value.'], pulseT/2);
    modT = pulseT / 2;
end

%% Computed Values
% cyclePeriod = 1 / pulseFreq;
% pulseT = cyclePeriod * dutyCycle;
breakT = cyclePeriod - pulseT;
pulseCFDelta = (1 / pulseT);
noBreak = dutyCycle == 1;

if noBreak
    stimTimeSteps = duration;
    numSteps = 1;
    refFreqSteps = carrierFreq;  
    pulseSel = 1;
    basePOSteps = pi;    
else    
    numPulses = ceil(duration / cyclePeriod) + 1;
    numSteps = numPulses * 2;
    cycleNumsHelper = [0, repelem(0:(numPulses - 1), 2)];
    stimTimeStepsHelper = [0, ...
        repmat([pulseT - modT, cyclePeriod], [1, numPulses])];

    stimTimeSteps = breakT ...
        + (stimTimeStepsHelper + (cycleNumsHelper .* cyclePeriod));
    stimTimeSteps = stimTimeSteps(1:(end - 1));
    cycleNums = repelem(1:numPulses, 2);    

    refFreqSteps = repmat([breakCFreq, carrierFreq], [1, numPulses]);   
    pulseSel = 2:2:numSteps;

    breakPO = pi;
    pulsePO = wrapToPi((pulseCFDelta * 2 * pi * modT) + pi);
    basePOSteps = repmat([breakPO, pulsePO], [1, numPulses]);
end

%% Define Frequence Modulation Steps
stimFreqSteps = repmat(refFreqSteps, [2, 1]);  
switch symetry
    case ModulationSymetry.SIGNAL1
        stimFreqSteps(1, pulseSel) = ...
            stimFreqSteps(1, pulseSel) + pulseCFDelta;
    case ModulationSymetry.SIGNAL2
        stimFreqSteps(2, pulseSel) = ...
            stimFreqSteps(2, pulseSel) + pulseCFDelta;
    case ModulationSymetry.SYMETRIC
         stimFreqSteps(:, pulseSel) = ...
            stimFreqSteps(:, pulseSel) + ([+1; -1] * pulseCFDelta / 2);
end

%% Define Phase Modulation Steps
switch symetry
    case ModulationSymetry.SIGNAL1
        % NaN phase offset for signal 2 since it will be the reference
        % phase. It's phase will be determind based on continuity.
        stimPOSteps = [basePOSteps; repelem(NaN, numSteps)];
    case ModulationSymetry.SIGNAL2
        % NaN phase offset for signal 1 since it will be the reference
        % phase. It's phase will be determind based on continuity.
        stimPOSteps = [repelem(NaN, numSteps); basePOSteps];
    case ModulationSymetry.SYMETRIC
        stimPOSteps = [basePOSteps; -basePOSteps] / 2;
end

%% Generate Reference Phase Function
refPhaseFcn = generateFMPhaseFcn(stimTimeSteps, refFreqSteps, modT, ...
    fmArgs{:});
end

