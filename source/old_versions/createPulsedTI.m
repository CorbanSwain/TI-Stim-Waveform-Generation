function [I, I1, I2] = createPulsedTI(carrier_f, pulse_f, stim_t, duty_cycle,...
    A1, A2, pre_t, post_t, ramp_up_t, ramp_down_t, dt, ...
    control, plot_on, h_on, save_ascii)
    % @author: pdzialecka
    
    % Pulsed TI simulator function
    % Carrier freq and pulse_f must be defined; other inputs are optional
    
    % Control mode generates unmodulated high frequency 'pulsed' signal -
    % doesn't work for anything else apart from continuous high f?
    
    % Time unit: ms
    % Frequency unit: Hz
    % phi1, phi2 - removed from inputs for now
    
    
    %% Default simulation values

    if ~exist('duty_cycle', 'var'); duty_cycle = 1; end

    if ~exist('A1', 'var'); A1 = 0.5; end
    if ~exist('A2', 'var'); A2 = 0.5; end

    % if ~exist('phi1', 'var'); phi1 = 0; end
    % if ~exist('phi2', 'var'); phi2 = pi; end

    if ~exist('stim_t', 'var'); stim_t = 1000; end         % 1 for s
    if ~exist('dt', 'var'); dt = 0.004; end                 % 0.01 = 10e-6 for s

    if ~exist('pre_t', 'var'); pre_t = 0; end
    if ~exist('post_t', 'var'); post_t = 0; end

    if ~exist('ramp_up_t', 'var'); ramp_up_t = 0; end     % 0.05 for s
    if ~exist('ramp_down_t', 'var'); ramp_down_t = 0; end

    if ~exist('control', 'var'); control = 0; end  % control = unmodulated high f

    if ~exist('plot_on', 'var'); plot_on = 1; end
    if ~exist('h_on', 'var'); h_on = 1; end

    if ~exist('save_ascii','var'); save_ascii = 0; end

    %% Testing as script
    % carrier_f = 1000;
    % pulse_f = 10;
    % stim_t = 0.4;
    % duty_cycle = 0.6;
    % phi1 = 0;
    % phi2 = pi;
    % A1 = 1;
    % A2 = 1;
    % dt = 10e-6;
    % h_on = 0;
    % 
    % pre_t = 0.1;
    % post_t = pre_t;
    % ramp_up_t = 0.05;
    % ramp_down_t = ramp_up_t;

    %% Set up stim parameters
    % account for break periods so intended f achieved

    if pulse_f == 0 % if no pulsation, treat like continuous high frequency
        control = 1;
        duty_cycle = 1;
    end

    if control
        phi1 = -pi/2;
        phi2 = -pi/2;
    else
        phi1 = 0;
        phi2 = pi;
    end

    % convert to ms
    pulse_f = pulse_f/1000;
    carrier_f = carrier_f/1000;

    f1 = carrier_f;
    f_diff = pulse_f/duty_cycle;  % f_diff not the same as pulse f!

    if control
        f2 = carrier_f;
    else
        f2 = carrier_f + f_diff;
    end

    % f2 = carrier_f + f_diff;

    pulse_t = 1/f_diff;
    break_t = 1/duty_cycle*pulse_t*(1-duty_cycle);
    cycle_t = pulse_t + break_t;

    stim_t = stim_t+ramp_up_t+ramp_down_t; % ramp times not included in stim_t var

    cycles = stim_t/cycle_t;
    half_cycles = round(cycles*2); % number of pulses

    tot_t = 0:dt:stim_t; % dt:dt:stim_t

    % **indexing**
    % tot_t: starting from 0 gives one extra point.

    % Creating signals: 1-5001, 5002-10001, 10002-15001 etc.
    % - first cycle one timestep longer (0 time)

    switch_A = 0; % hardcoded for now

    %% Create signals I1 and I2
    if duty_cycle == 0.5  % easier case - switch f2 to f1 during break

        I1 = A1*cos(2*pi*f1*tot_t+phi1);
        I2 = A2*cos(2*pi*f2*tot_t+phi2);

        if control
            phi2 = pi/2; % for breaks only
        end

        for i=1:2:half_cycles
            idx = round(i*pulse_t/dt)+2;
            break_idxs = idx:1:idx+break_t/dt-1;

            if break_idxs(end) > length(tot_t)
                break_idxs = idx:1:length(tot_t);
            end

            break_ts = tot_t(break_idxs);
            if switch_A
                I2(break_idxs) = A1*cos(2*pi*f1*break_ts+phi2); % I2 switch + A2 -> A1
            else
                I2(break_idxs) = A2*cos(2*pi*f1*break_ts+phi2);
            end
        end

    elseif duty_cycle == 1

        I1 = A1*cos(2*pi*f1*tot_t+phi1);
        I2 = A2*cos(2*pi*f2*tot_t+phi2);

    elseif duty_cycle == 0
    %     f1 = 0;
    %     f2 = 0;
        A1 = 0;
        A2 = 0;
        I1 = A1*cos(2*pi*f1*tot_t+phi1);
        I2 = A2*cos(2*pi*f2*tot_t+phi2);

    else
        I1 = zeros(1,length(tot_t));
        I2 = zeros(1,length(tot_t));

        start_idx = 1;
        phase_shift = 1; % use more realistic phase shift

        for i=1:half_cycles
            
            % PULSE TIME
            if mod(i,2)
    %             start_idx

                if control
                    phi2 = phi1; % keep same phase during pulse
                end

                end_idx = start_idx + round(pulse_t/dt) - (i~=1); % don't subtract when i==1
                % -1 for all for dt start

                pulse_idxs = start_idx:1:end_idx;

                if pulse_idxs(end) > length(tot_t)
                    pulse_idxs = start_idx:1:length(tot_t);
                end


                if phase_shift % change the phase for offset (dc ~= 0.5)
                    pulse_ts = tot_t(pulse_idxs);

                    if i == 1
                        phase1 = 0;
                        phase2 = 0;
                    else
                        % offset forwards - smaller phases
                        time_shift = (i-1)/2*(pulse_t-break_t); %*10e3; % in ms
                        phase1 = 2*pi*f1*time_shift;
                        phase2 = 2*pi*f2*time_shift;

                        % scaling factor
                        c1 = floor(phase1/(2*pi));
                        c2 = floor(phase2/(2*pi));

                        % phase in the range -2pi to +2pi
                        phase1 = phase1 - c1*2*pi;
                        phase2 = phase2 - c2*2*pi;


                        % offset backwards - easier
    %                     time_offset = cycle_t*(i-1)/2;  % = pulse_ts(1);
    %                     phase1 = -2*pi*f1*time_offset;
    %                     phase2 = -2*pi*f2*time_offset;
                    end

                    % add phase offset 
                    I1(pulse_idxs) = A1*cos(2*pi*f1*pulse_ts+phi1+phase1);
                    I2(pulse_idxs) = A2*cos(2*pi*f2*pulse_ts+phi2+phase2);

                    start_idx = end_idx+1;


                % NOT USED ANYMORE
                else  % move indeces of time for offset (not realistic)

                    if i == 1
                        pulse_idxs_off = pulse_idxs;
                    else
                        % 'restarting' the time to 0 also works
                        pulse_idxs_off = 1:length(pulse_idxs);
    %                     pulse_idxs_off = pulse_idxs - round((i-1)/2* break_t/dt); % offset the pulse
                    end

                    pulse_ts = tot_t(pulse_idxs_off);

                    I1(pulse_idxs) = A1*cos(2*pi*f1*pulse_ts+phi1);
                    I2(pulse_idxs) = A2*cos(2*pi*f2*pulse_ts+phi2);

                    start_idx = end_idx+1;
                end

            % BREAK TIME
            else
            % set f2 to f1 and A2 to A1 (if different) during break

                if control
                    phi2 = pi/2; % change phase for break
                end

                end_idx = start_idx + round(break_t/dt) - 1;
                break_idxs = start_idx:1:end_idx;

                if break_idxs(end) > length(tot_t)
                    break_idxs = start_idx:1:length(tot_t);
                end

                break_ts = tot_t(break_idxs);
                I1(break_idxs) = A1*cos(2*pi*f1*break_ts+phi1);
                if switch_A
                    I2(break_idxs) = A1*cos(2*pi*f1*break_ts+phi2); % use f1
                else
                    I2(break_idxs) = A2*cos(2*pi*f1*break_ts+phi2); % use f1
                end
                start_idx = end_idx+1;
            end
            
        end
    end


    %% Add ramping to the signals
    if ramp_up_t
        ramp_up_size = round(ramp_up_t/dt);
        ramp_vec = 0:1/ramp_up_size:1;

        I1(1:length(ramp_vec)) = ramp_vec.*I1(1:length(ramp_vec));
        I2(1:length(ramp_vec)) = ramp_vec.*I2(1:length(ramp_vec));
    end

    if ramp_down_t
        ramp_down_size = round(ramp_down_t/dt);
        ramp_vec = fliplr(0:1/ramp_down_size:1);

        I1(end-length(ramp_vec)+1:end) = ramp_vec.*I1(end-length(ramp_vec)+1:end);
        I2(end-length(ramp_vec)+1:end) = ramp_vec.*I2(end-length(ramp_vec)+1:end);
    end


    %% Add pre and post times
    if pre_t || post_t
        I1 = cat(2,zeros(1,pre_t/dt),I1,zeros(1,post_t/dt));
        I2 = cat(2,zeros(1,pre_t/dt),I2,zeros(1,post_t/dt));
        tot_t = 0:dt:stim_t+post_t+pre_t;%+ramp_up_t+ramp_down_t;
    end

    %% Display created signal
    I = I1+I2;

    if plot_on
        fig = figure, plot(tot_t,I);  % tot_t in ms i think
        xlim([0,tot_t(end)]);

        if h_on
            y_h = hilbert(I);
            y_env = abs(y_h);
            hold on
            plot(tot_t, y_env, 'r');
        end

        % Save the figure
        filename = fullfile('.', 'signal_plot.png');
        saveas(fig, filename);
    end

    if save_ascii
        name1 = sprintf('dc_%1.1f_I1_cf_%d_pf_%d_A1_%1.1d_A2_%1.1d.txt',...
                        duty_cycle,carrier_f*1000,pulse_f*1000,A1,A2); % 'I1.txt';
        t_I1 = [tot_t;I];
        save(name1,'t_I1','-ascii')

        name2 = sprintf('dc_%1.1f_I2_cf_%d_pf_%d_A1_%1.1d_A2_%1.1d.txt',...
                        duty_cycle,carrier_f*1000,pulse_f*1000,A1,A2); % 'I2.txt';
        t_I2 = [tot_t;I];
        save(name2,'t_I2','-ascii')
    end

    %%
    %%%%%%%%%%%% older/unused for now %%%%%%%%%%%%%%%%%%%%


    % extend time if needed to ensure full pulses only are send? necessary?
    % fix this if want to use it

    % % if mod(sim_t, pulse_t) && mod(half_cycles, 2)
    % if mod(sim_t, cycle_t) || sim_t/cycle_t - round(sim_t/cycle_t) ~= 0 % && mod(half_cycles, 2) - keep full cycle for now
    % %     sim_t = half_cycles*pulse_t;
    %     sim_t = round(cycles)*cycle_t; %(ceil(cycles)-1)*cycle_t + pulse_t;
    % 
    %     display(['Extending simulation time to ', num2str(sim_t), 's to ensure full pulses only are applied'])
    %     cycles = sim_t/cycle_t;
    % end


    % if mod(cycles,2)
    %     if round(cycles) < cycles %round(cycles,1) < round(cycles)
    %         cycles = round(cycles)-0.5
    %     elseif round(cycles,1) > round(cycles)
    %         cycles = round(cycles);
    %     end
    % end
