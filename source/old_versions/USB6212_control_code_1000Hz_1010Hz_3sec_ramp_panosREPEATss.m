%% Define text of questions to be asked ot the user.
% % prompt1 = 'What is the stimulation frequency you want to use for F1? ';
% % prompt2 = 'What is the stimulation frequency you want to use for F2? ';
% % prompt3 = 'What is the stimulation duration you want to use (in seconds)? ';
prompt4 = 'What is the stimulation intensity you want to use for 1kHz (in Volts; min: 0, max: 10)? ';
prompt5 = 'What is the stimulation intensity you want to use for 1.01kHz (in Volts; min: 0, max: 10)? ';
prompt6 = 'What is the delay you want to use before starting the stimulation (in seconds)? ';
prompt7 = 'What is the delay you want to use AFTER starting the stimulation (in seconds)? ';
prompt8 = 'Time that the laser is on BEFORE the experiments (in seconds)? ';
prompt9 = 'How many repetitions of stimulation protocol?';
prompt10= 'How many SEC between repetitions?';


% prompt7 = 'What is the ramp up time you want to use for the signals (in seconds)?';2
% prompt8 = 'What is the ramp down time you want to use for the signals (in seconds)?';
%%---------------------------------------------
%% DEFINE Variables (either directly or through user-questions(input) command)
% % f1 = input(prompt1);
% % f2 = input(prompt2);
% % duration = input(prompt3);  
f1 =1000;
f2 =1100;
% Intensity Variables
intensity1 = input(prompt4);
intensity2 = input(prompt5);

%Variables for delaying start_of_ramping,start_of_camera, end_of_camera 
delay_time = input(prompt6); 
delay_time_after=input(prompt7);
delay_time_before=input(prompt8); 

% Variables for repeating ramps:
num_of_repetions=input(prompt9); %Or directly put number after the = 
delay_between_repetions=input(prompt10); %Or directly put number after the = 

% Variables for defining ramping_up and down
% ramp_up_t = input(prompt7);
% ramp_down_t = input(prompt8);
ramp_up_t = 0.003;
ramp_down_t = 0.003;
duration = 0.003 + ramp_up_t + ramp_down_t;

sampling_rate = 100000;  % 100 kHz (maximum for USB6212 is 250kHz)
dt = 1/sampling_rate;

total_duration=delay_time_before+delay_time+num_of_repetions.*(duration+delay_between_repetions)+delay_time_after
display(['Expected total time in seconds:' num2str(total_duration)]);
%%--------------------------------------------


%% Definition of the ramping -Do not play with it directly(only throught the varibales above)
% t=0:1/sampling_rate:(duration - 1/sampling_rate);
t = 0:dt:(duration-dt);

y1=intensity1*sin(2*pi*f1*t+pi);
y2=intensity2*sin(2*pi*f2*t);

if ramp_up_t
    ramp_up_size = round(ramp_up_t/dt);
    ramp_vec = 0:1/ramp_up_size:1;

    ramp_up_tt = dt:dt:ramp_up_t;
    ramp_vec = ramp_vec(1:length(ramp_up_tt));
 
    y1(1:length(ramp_vec)) = ramp_vec.*y1(1:length(ramp_vec));
    y2(1:length(ramp_vec)) = ramp_vec.*y2(1:length(ramp_vec));
end

if ramp_down_t
    ramp_down_size = round(ramp_down_t/dt);
    ramp_vec = fliplr(0:1/ramp_down_size:1);

    ramp_down_tt = dt:dt:ramp_down_t;
    ramp_vec = ramp_vec(1:length(ramp_down_tt));

    y1(end-length(ramp_vec)+1:end) = ramp_vec.*y1(end-length(ramp_vec)+1:end);
    y2(end-length(ramp_vec)+1:end) = ramp_vec.*y2(end-length(ramp_vec)+1:end);
end

    y1(end) = 0;
    y2(end) = 0;
% End of the defination of ramping.
%% Initialization of DAQs
%Initialize Camera Trigger control
s1 = daq.createSession('ni');
addDigitalChannel(s1,'Dev1','Port1/Line0:2', 'Outputonly');

%Initialize Analog-voltage control
s2 = daq.createSession('ni');
addAnalogOutputChannel(s2,'Dev1',[0 1],'Voltage');
s2.Rate = sampling_rate;
queueOutputData(s2,[y1' y2']);

%Initialize Analog-voltage control
s3_laserPower = daq.createSession('ni');
addDigitalChannel(s3_laserPower,'Dev1', 'Port0/Line0', 'OutputOnly');

%--so far nothing has been send to the DAQs---

%% ACTUAL EXPERIMENT.

tic %command to start elapsed time measurments.

% Lines bellow send a trigger to the laser.
s3_laserPower.outputSingleScan([1]);
outputSingleScan(s3_laserPower,[1])
display(['LaserStarted' num2str(toc)]);
time_to_laser_start=toc; 
outputSingleScan(s1,[0 0 0]) %this just initialized DAQ

% We wait for some time before sending main trigger to CAMERA
pause(delay_time_before);

% Command send to main DAQ-->CAMERA
display(['StimultionStarted' num2str(toc)]);
outputSingleScan(s1,[1 1 1])
time_to_camera_start=toc; 

% We wait for some time before starting the ramp up/down process.
pause(delay_time);
time_to_ramp_start=toc; 

for i_rep=1:num_of_repetions
    display(['Repetition' num2str(i_rep) 'of' num2str(num_of_repetions) 'starts...'])
    % Command send to start ramping up & down.
    startForeground(s2);
    outputSingleScan(s2,[0 0]);
    display(['StimultionStoped' num2str(toc)]);
    time_to_ramp_stop=toc;
    queueOutputData(s2,[y1' y2']);
    pause(delay_between_repetions)
end

% We keep recording for some more time.
pause(delay_time_after);
time_to_camera_stop=toc;

% We turn of the laser.
outputSingleScan(s3_laserPower,[0])
display(['LaserStopped' num2str(toc)]);
time_to_laser_stop=toc;

