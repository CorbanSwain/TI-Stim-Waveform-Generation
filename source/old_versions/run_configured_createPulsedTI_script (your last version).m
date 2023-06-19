%% Get Inputs for Pulse Generation

% template for inputs
% input name | required? | default value | description
inputSpecification = {
    % The carrier frequency of the signals in Hz.
    'carrier_f', 2000, 'Hz'
    
    % The frequency of the pulses in Hz.
    'pulse_f', 1, 'Hz'
    
    % 'The duration of the stimulation in ms.'
    'stim_t', 3000, 'ms'
    
    % The duty cycle of the pulses. A value of 1 means the pulses are 
    % continuous, while a value of 0 means no pulses are present.
    'duty_cycle', 0.02, 'a.u.'
    
    % The amplitude of the first signal, I1.
    'A1', 1, 'V'
    
    % The amplitude of the second signal, I2.
    'A2', 1, 'V'
    
    % The duration of time before the stimulation starts 
    % pre-stimulation time) in ms.
    'pre_t', 0, 'ms'
    
    % The duration of time after the stimulation ends 
    % (post-stimulation time) in ms.
    'post_t', 0, 'ms'
    
    % The duration of the ramp-up period of the signals in ms, i.e., the
    % time for the signals to gradually increase from 0 to their full 
    % amplitude.
    'ramp_up_t', 500, 'ms'
    
    % The duration of the ramp-down period of the signals in ms, i.e., 
    % the time for the signals to gradually decrease from their full 
    % amplitude to 0.
    'ramp_down_t', 500, 'ms'
    
    % The time step in ms used for generating the signals. It corresponds 
    % to the time resolution of the signals.
    'dt', 0.000001, 'ms'
    
    % A flag to generate an unmodulated high-frequency signal.
    'control', false, ''
    
    % A flag to plot the output signal, I. If set to 1, the signal is 
    % plotted.
    'plot_on', true, ''
    
    % A flag to plot the Hilbert transform of the output signal on the 
    % same plot as I.
    'h_on', true, ''
    
    % A flag to save the output signals (I, I1, and I2) as text files. If 
    % set to 1/true, the signals are saved.
    'save_ascii', false, ''};

% let's create a loop for prompting values for each of these inputs
fprintf('Set Parameters:\n')
numArgs = size(inputSpecification, 1);
argumentList = cell(1, numArgs);

for iArg = 1:numArgs
    [name, value, unit] = inputSpecification{iArg, :};
    if islogical(value)
        if value
            valStr = 'true';
        else
            valStr = 'false';
        end
    elseif isnumeric(value)
        valStr = sprintf('%f', value);
    else
        error(['Unexpeced non-numeric, non-logical value given for', ...
           ' parameter "%s"'], name)
    end              
    fprintf('%12s = %14s %s\n', name, valStr, unit);
    argumentList{iArg} = value;
end

%% running function with arguments
fprintf('Running function...\n');
[I, I1, I2] = createPulsedTI(argumentList{:});
fprintf('Function run complete.\n');

%% Plotting output results
stimTime = argumentList{3};

numSamples = size(I, 2);
timeArr = linspace(0, stimTime, numSamples);

f = figure();
ax0 = subplot(3, 1, 1);
plot(timeArr, I);
ylabel('I');

ax1 = subplot(3, 1, 2);
plot(timeArr, I1);
ylabel('I1');

ax2 = subplot(3, 1, 3);
plot(timeArr, I2);
ylabel('I2');
xlabel('Time (ms)');
linkaxes([ax0, ax1, ax2], 'x');

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
queueOutputData(stimSession,[I1', I2']);

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