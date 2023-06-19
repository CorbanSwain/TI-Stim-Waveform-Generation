%% Get Inputs for Pulse Generation

% template for inputs
% input name | required? | default value | description
inputSpecification = {
    'carrier_f', true, [], 'The carrier frequency of the signals in Hz.'
    'pulse_f', true, [] 'The frequency of the pulses in Hz.'
    'stim_t', false, 1000, 'The duration of the stimulation in ms.'
    'duty_cycle', false, 1, 'The duty cycle of the pulses. A value of 1 means the pulses are continuous, while a value of 0 means no pulses are present.'
    'A1', false, 0.5, 'The amplitude of the first signal, I1.'
    'A2', false, 0.5, 'The amplitude of the second signal, I2.'
    'pre_t', false, 0, 'The duration of time before the stimulation starts (pre-stimulation time) in ms.'
    'post_t', false, 0, 'The duration of time after the stimulation ends (post-stimulation time) in ms.'
    'ramp_up_t', false, 0, 'The duration of the ramp-up period of the signals in ms, i.e., the time for the signals to gradually increase from 0 to their full amplitude.'
    'ramp_down_t',false, 0, 'The duration of the ramp-down period of the signals in ms, i.e., the time for the signals to gradually decrease from their full amplitude to 0.'
    'dt', false, 0.01, 'The time step in ms used for generating the signals. It corresponds to the time resolution of the signals.'
    'control', false, false, ' A flag to generate an unmodulated high-frequency signal.'
    'plot_on', false, true, ' A flag to plot the output signal, I. If set to 1, the signal is plotted.'
    'h_on', false, true, 'A flag to plot the Hilbert transform of the output signal on the same plot as I.'
    'save_ascii', false, false, 'A flag to save the output signals (I, I1, and I2) as text files. If set to 1, the signals are saved.'};


% let's create a loop for prompting values for each of these inputs
numArgs = size(inputSpecification, 1);
argumentList = cell(1, numArgs);

for iArg = 1:numArgs
    fprintf('---\n');
    [name, isReq, defaultVal, descrip] = inputSpecification{iArg, :};
   
    if isReq
        isLogVal = false;  % FIXME - this should not necessarily be hard-coded
        reqStr = 'REQUIRED.';
        defaultStr = '';
        dValStr = '[]';
    else
        isLogVal = islogical(defaultVal);
        reqStr = 'NOT required, press ENTER to skip.';
        defaultStrFmt = 'deinputfault value = %s\n';
        if isLogVal
            if defaultVal
                dValStr = 'true';
            else 
                dValStr = 'false';
            end
        else
            dValStr = sprintf('%f', defaultVal);
        end
        defaultStr = sprintf(defaultStrFmt, dValStr);
    end
    
    promptStr = sprintf(...
        'Set the value for >> %s <<\nValue is %s\n%sDescription: %s\n? ', ...
        name, reqStr, defaultStr, descrip);
    
    inputVal = input(promptStr);
    
    if (isempty(inputVal))
        if isReq
            error('The input "%s" is required, but a blank value was recieved.', ...
                name);
        else
            fprintf('Using default value; %s = %s\n', name, dValStr)
            argVal = defaultVal;
        end
    else
        if isLogVal
            if isLogical(inuptVal)
                argVal = inputVal;
            elseif isnumeric(inputVal)
                argVal = logical(inputVal);
            else
                error('Got an unexpected input for a true/false parameter.');
            end
            
            if argVal
                valStr = 'true';
            else
                valStr = 'false';
            end
        else
            argVal = inputVal; 
            valStr = sprintf('%f', argVal);   
        end
              
        fprintf('Value set; %s = %s\n', name, valStr);
    end
    argumentList{iArg} = argVal;
end

fprintf('---\n\nArguments List:\n');
disp(inputSpecification(:, 1))
disp(argumentList);

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
% fprintf('Press ENTER to continue\n ');
input('');
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