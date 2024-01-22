% Script for Configuring Stimulation, Generating Waveforms, and Running
% TI Stimulation Experiments
% Corban Swain , June 2023

%% 1. Get Inputs for Pulse Generation
clear();

% template for inputs
% input name | value
inputSpecification = {
    % 'The duration of the stimulation in s.'
    'Duration', 5 % seconds

    % The carrier frequency of the signals in Hz.
    'CarrierFreq', 2000 % hertz

    % The carrier frequency during the break of the signals in Hz.
    'BreakCarrierFreq', 5000 % hertz (or 'same')

    % The frequency of the pulses in Hz.
    'PulseFreq', 10 % hertz   
    
    % The duty cycle of the pulses. A value of 1 means the pulses are 
    % continuous, while a value of 0 means no pulses are present.
    'DutyCycle', 0.5 % (0 <= DutyCycle <= 1)
    
    % The amplitude of the first signal, S1.
    'A1', 1 % volts
    
    % The amplitude of the second signal, S2.
    'A2', 1 % volts

    % The sampling rate in Hz used for generating the signals.
    'SamplingRate', 200e3 % hertz
    
    % The duration of the ramp-[up, down] period of the signals in s,
    % i.e., the time for the signals to gradually increase from 0 to their
    % full amplitude.
    'RampTime', 0.5 % seconds   
    
    % Whether or not to flip one of the signals (shift in phase 180 deg).
    % This might be useful to correct for hardware setup errors/details.
    'Flip', false % (true/false)
    
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
    % differ as specified above.
    'ModulationSymetry', 'signal1'  

    % A flag to generate an unmodulated high-frequency signal.
    'Control', false % (true/false)
    
    % A flag to plot a summary of the output signals.
    'Plot', true % (true/false)

    % A cell array of Name, Value arguments sent to generateFMPhaseFcn
    'FMArgs', {} 
    
    % A flag to run the code in debug mode.
    'Debug', false % (true/false)    
    };

% let's create a loop for printing values for each of these inputs
fprintf('Set Parameters for stimulation:\n')
numArgs = size(inputSpecification, 1);
for iArg = 1:numArgs
    [name, value] = inputSpecification{iArg, :};            
    fprintf('%20s = %-14s\n', name, utils.toStr(value));
end
paramStruct = ...
    cell2struct(inputSpecification(:, 2), inputSpecification(:, 1));

%% 2. Use Configuration to Generate Waveforms
fprintf('---\n')
fprintf('Running waveform generation function...\n');
[waveforms, ~, usedParamStruct] = createTIStimWaveforms(paramStruct);
fprintf('Waveforms sucessfully generated.\n');

%% 3. Run Stimulation Experiment Sequence on NIDAQ
fprintf('---\n')
fprintf('Preparing to run system control sequence for stimulation...\n');
%%% Variables
% delay between laser on and device trigger (s)
preDeviceDelay = 1;

% delay between device trigger and stimulation start (s)
preStimulationDelay = 0.1;

% delay between stimulation stop and laser off (s)
postStimulaionDelay = 1;

% DAQ device name
deviceName = 'Dev2';

samplingRate = usedParamStruct.SamplingRate; % must be same as configured

%%% NIDAQ Setup
vendorID = 'ni';
newDaqApi = ~verLessThan('matlab', '9.8');

% detect devices
if newDaqApi
    devs = daqlist();
else
    devs = daq.getDevices();
end
if isempty(devs)
    error(['No DAQ device is detected, cannot continue with', ...
        ' stimulation. Connect to or check the connection to the', ...
        ' appropriate DAQ device to resolve this error.'])
end

% create device triggeriong daq session
% this controls; [1.0: Not Connected, 1.1: Camera, 1.2: Axon CNS Digidata]
devTriggerLines = 'Port1/Line0:2';
if newDaqApi
    deviceTriggerDAQ = daq(vendorID);
    deviceTriggerDAQ.addoutput(deviceName, devTriggerLines, 'Digital');
else
    deviceTriggerSession = daq.createSession(vendorID);
    addDigitalChannel(...
        deviceTriggerSession,'Dev2', devTriggerLines, 'Outputonly');
end

% create laser triggering daq session
% this controls the imaging laser
laserLine = 'Port0/Line0';
if newDaqApi
    laserDAQ = daq(vendorID);
    laserDAQ.addoutput(deviceName, laserLine, 'Digital');
else
    laserSession = daq.createSession(vendorID);
    addDigitalChannel(laserSession,'Dev2', laserLine, 'OutputOnly');
end

% create stimulation daq session
% this controls the stimulation signals
stimLines = [0, 1];
if newDaqApi
    stimDAQ = daq(vendorID);
    stimDAQ.addoutput(deviceName, stimLines, 'Voltage');
    stimDAQ.Rate = samplingRate;
else
    stimSession = daq.createSession(vendorID);
    addAnalogOutputChannel(stimSession, deviceName, stimLines, 'Voltage');
    stimSession.Rate = samplingRate;
    queueOutputData(stimSession, waveforms);
end

if newDaqApi
    daqSingleOutput = @(daq, data) daq.write(data);
    daqForegroundOutput = daqSingleOutput;

    stimIO = stimDAQ;
    laserIO = laserDAQ;
    devTrigIO = deviceTriggerDAQ;
else
    daqSingleOutput = @(sess, data) outputSingleScan(sess, data);
    daqForegroundOutput = @(sess, ~) startForeground(sess);

    stimIO = stimSession;
    laserIO = laserSession;
    devTrigIO = deviceTriggerSession;
end

% make sure devices (camera, digidata) are off
% make sure stimulation is off
daqSingleOutput(devTrigIO, [0, 0, 0]);
daqSingleOutput(stimIO, [0, 0]);

%%% Experiment Control Sequence
experimentStatusFmt = 'NIDAQ : %25s at %6.3f s. \n';
fprintf('---\n');
t0 = tic(); % begin timer

% > start laser
fprintf(experimentStatusFmt, 'Starting Laser', toc(t0));
daqSingleOutput(laserIO, 1);

% > pause
pause(preDeviceDelay);

% > trigger devices
fprintf(experimentStatusFmt, 'Triggering Devices', toc(t0));
daqSingleOutput(devTrigIO, [1, 1, 1]);

% > pause
pause(preStimulationDelay);

% > begin stimulation
fprintf(experimentStatusFmt, 'Beginning Stimulation', toc(t0));
daqForegroundOutput(stimIO, waveforms);

% > stop stimulation
fprintf(experimentStatusFmt, 'Ending Stimulation', toc(t0));
daqSingleOutput(stimIO, [0, 0]);

% > pause
pause(postStimulaionDelay);

% > turn off lasers and devices
fprintf(experimentStatusFmt, 'Stopping Laser', toc(t0));
daqSingleOutput(laserIO, 0);
daqSingleOutput(devTrigIO, [0, 0, 0]);

fprintf('---\n\nDONE, Experiment Run Complete!\n\n');

%%% Cleanup
% make sure session variables are destroyed
clear('laserIO', 'stimIO', 'deviceTriggerIO');
if newDaqApi
    clear('laserDAQ', 'stimDAQ', 'deviceTriggerDAQ');
else
    clear('laserSession', 'stimSession', 'deviceTriggerSession');
end