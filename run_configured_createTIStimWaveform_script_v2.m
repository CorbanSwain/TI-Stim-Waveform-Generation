%% Get Inputs for Pulse Generation

% template for inputs
% input name | value
inputSpecification = {
     % 'The duration of the stimulation in s.'
    'Duration', 3 % seconds

    % The carrier frequency of the signals in Hz.
    'CarrierFreq', 2000 % hertz

    % The carrier frequency of the signals in Hz.
    'BreakCarrierFreq', 2000 % hertz

    % The frequency of the pulses in Hz.
    'PulseFreq', 1 % hertz   
    
    % The duty cycle of the pulses. A value of 1 means the pulses are 
    % continuous, while a value of 0 means no pulses are present.
    'DutyCycle', 0.2 % (0 <= DutyCycle <= 1)
    
    % The amplitude of the first signal, I1.
    'A1', 1 % volts
    
    % The amplitude of the second signal, I2.
    'A2', 1 % volts

    % The sampling rate in Hz used for generating the signals.
    'SamplingRate', 200e3 % hertz
    
    % The duration of the ramp-[up, down] period of the signals in s,
    % i.e., the time for the signals to gradually increase from 0 to their
    % full amplitude.
    'RampTime', 0.5 % seconds    
    
    % the modulation time during which the wavforms will smoothly
    % transition from an initial frequency and relative phase to a new
    % frequency and relative phase in s.
    'FMTime', 'auto' % seconds (or 'auto')  

    % The duration of time [before, after] the stimulation starts 
    % pre-stimulation time) in s.
    'WaitTime', 0.1 % seconds    
    
    % Controls which signals are modulated to produce the stimulation
    % pulses. Valid values:
    %   - 'signal1'  : only the first waveform will be modulated while the
    %                  second remains at the carrierFreq (or breakFreq).
    %   - 'signal2'  : only the second waveform will be modulated.
    %   - 'symetric' : both waveforms are modulated symetrically to produce
    %                  the desired summation waveform.
    % All values of symetry will produce the same summation waveform,
    % however the frequency and phase modulations to the two waveforms
    % differently as specified above.
    'ModulationSymetry', 'signal1'  

    % A flag to generate an unmodulated high-frequency signal.
    'Control', false % (true/false)
    
    % A flag to plot a summary of the output signals.
    'Plot', true % (true/false)

    % A cell array of Name, Value arguments sent to generateFMPhaseFcn
    'FMArgs', [{}]   
    
    % A flag to run the code in debug mode.
    'Debug', false % (true/false)
    };

% let's create a loop for printing values for each of these inputs
fprintf('Set Parameters:\n')
numArgs = size(inputSpecification, 1);
for iArg = 1:numArgs
    [name, value] = inputSpecification{iArg, :};            
    fprintf('%20s = %-14s\n', name, csmu.toStr(value));
end

paramStruct = ...
    cell2struct(inputSpecification(:, 2), inputSpecification(:, 1));

%% running function with arguments
fprintf('Running function...\n');
[waveforms, timeArray, usedParamStruct] = ...
    createTIStimWaveforms(paramStruct);
fprintf('Function run complete.\n');

%% Running Stimulation Sequence on NIDAQ
fprintf('Preparing to run system control sequence for stimulation...\n');
%%% Variables
% delay between laser on and device trigger (s)
delayTimeBefore= 1;
% delay between device trigger and stimulation start (s)
delay_time = 0.5;
% delay between stimulation stop and laser off (s)
delayTimeAfter= 1;
samplingRate = 1 / (argumentList{11} / 1e3);

%%% Setup
% create device triggeriong daq session
deviceTriggerSession = daq.createSession('ni');
% this controls; [1.0: Not Connected, 1.1: Camera, 1.2: Axon CNS Digidata]
addDigitalChannel(...
    deviceTriggerSession,'Dev2','Port1/Line0:2', 'Outputonly');

% create laser triggering daq session
laserSession = daq.createSession('ni');
% this controls the imaging laser
addDigitalChannel(laserSession,'Dev2', 'Port0/Line0', 'OutputOnly');

% create stimulation daq session
stimSession = daq.createSession('ni');
% this controls the stimulation signals
addAnalogOutputChannel(stimSession,'Dev2',[0 1],'Voltage');
stimSession.Rate = samplingRate;
queueOutputData(stimSession, waveforms);

% make sure devices (camera, digidata) are off
outputSingleScan(deviceTriggerSession, [0, 0, 0]);
% make sure stimulation is off
outputSingleScan(stimSession, [0, 0]);

%%% control sequence
% > begin timer
t0 = tic();
% > start laser
outputSingleScan(laserSession, 1);
fprintf('Laser Started at %6.3f s.\n', toc(t0));
% > pause
pause(delayTimeBefore);
fprintf('Beginning Stimulation at %6.3f s.\n', toc(t0));
% > trigger devices
outputSingleScan(deviceTriggerSession, [1, 1, 1])
% > pause
pause(delay_time);
% > begin stimulation
startForeground(stimSession);
% > stop stimulation
outputSingleScan(stimSession, [0, 0]);
fprintf('Stimultion Ended at %6.3f s.\n', toc(t0));
% > pause
pause(delayTimeAfter);
% > turn off lasers
outputSingleScan(laserSession, 0);
fprintf('Laser Stopped at %6.3f s.\n', toc(t0));
fprintf('\n---\nRun complete!\n');