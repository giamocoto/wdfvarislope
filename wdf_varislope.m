%% wdf_varislope.m
% Wave Digital Filter (WDF) simulation of the Serge Variable Slope
% (Varislope) filter.
%
% Read the article on vesia.it
%
% For each slope setting the script:
%   1. builds the R-node scattering matrix,
%   2. runs the WDF sample by sample on the chosen input signal,
%   3. identifies a 2-pole transfer function from the simulated response
%      and plots its magnitude Bode diagram.
%
% Requirements: Control System Toolbox (bodeoptions, bodeplot, tf) and
%               System Identification Toolbox (iddata, tfest).

clear; close all; clc;

%% Simulation settings
Fs  = 192000;               % sampling frequency [Hz]
Ts  = 1/Fs;                 % sampling period [s]
sim_duration = 1;           % simulation length [s]
Nsamp = round(Fs*sim_duration);

% Input signal: 'step' or 'impulse'
input_type = 'step';
step_delay = 0.01;          % time at which the step starts [s]

Vin = zeros(Nsamp, 1);
switch input_type
    case 'step'
        Vin(round(step_delay*Fs)+1:end) = 1;
    case 'impulse'
        Vin(10) = 1;
    otherwise
        error('Unknown input_type "%s". Use ''step'' or ''impulse''.', input_type);
end

%% Slope settings to simulate
% User-facing slope value: 1 = maximum slope, 50 = minimum slope.
% Internally gm1 = gm_m*(51-slope) and gm2 = gm_m/(51-slope).
slope_setting = [1 35 45 50];
slope_array   = 51 - slope_setting;
n_iter        = length(slope_array);

%% Circuit values (independent of slope)
% External resistors
Zf = 10e3;                  % R1
Zb = 20e3;                  % R2
Za =  1e3;                  % R3
Zc = 20e3;                  % R4
% OTA transconductance (reference value) and load capacitors
gm_m = 10e-6;               % [S]
C    = 800e-12;             % [F]

% Gain factor of the resistive network: expected DC gain of the filter
% output after the step transient settles
beta = Zc/((Za*Zb)/(Za+Zb)) + 1;

%% Bode plot options
figure;
hold on;
bode_opts = bodeoptions('cstprefs');
bode_opts.FreqUnits    = 'Hz';
bode_opts.PhaseVisible = 'off';
bode_opts.XLim         = {[20 2e4]};

tic
for jj = 1:n_iter

    %% Slope-dependent OTA transconductances
    slope = slope_array(jj);
    gm1 = gm_m*slope;
    gm2 = gm_m/slope;

    %% Theoretical poles (reference only)
    % Pole locations from the analytical 2nd-order transfer function.
    % Damping ratio < 0.5 means the filter is resonating.
    s1 = (gm_m*slope/(2*beta*C)) * (-1 + sqrt(1 - 4*beta/slope^2));
    s2 = (gm_m*slope/(2*beta*C)) * (-1 - sqrt(1 - 4*beta/slope^2));
    fn1   = abs(s1)/(2*pi);     % natural frequency of pole 1 [Hz]
    fn2   = abs(s2)/(2*pi);     % natural frequency of pole 2 [Hz]
    zeta1 = -real(s1)/abs(s1);  % damping ratio of pole 1
    zeta2 = -real(s2)/abs(s2);  % damping ratio of pole 2
    fprintf('slope %2d: pole 1 = %.1f Hz (zeta = %.3f), pole 2 = %.1f Hz (zeta = %.3f)\n', ...
            slope_setting(jj), fn1, zeta1, fn2, zeta2);

    %% WDF port impedances / conductances
    % Capacitor port resistance given by the bilinear transform: Z = Ts/(2C)
    Zd = Ts/(2*C);
    Ze = Ts/(2*C);

    Ga = 1/Za;  Gb = 1/Zb;  Gc = 1/Zc;
    Gd = 1/Zd;  Ge = 1/Ze;  Gf = 1/Zf;

    %% MNA matrix of the R-node (built by inspection / stamp method)
    % Unknowns (columns): node voltages 1..12 | norator currents nor1, nor2 |
    %                     port currents eA..eF
    % Equations (rows):   KCL at nodes 1..12 | nullator constraints nul1, nul2 |
    %                     port constraints eA..eF
    %
    %   Ports: A = Ra, B = Rb, C = Rc, D = Cd, E = Ce, F = Rf (input source)
    %   Nullator/norator pairs model the ideal op-amps (virtual short).
    %   gm1, gm2 stamp the two OTAs as voltage-controlled current sources.
    %
    %        1    2    3    4    5    6    7    8    9   10   11   12  | nor1 nor2  eA   eB   eC   eD   eE   eF
    X = [   Ga, -Ga,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,     0,   0,  +1,   0,   0,   0,   0,   0;  % KCL 1
           -Ga,  Ga,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,     0,   0,   0,  -1,  -1,   0,   0,   0;  % KCL 2
             0,   0,  Gb, -Gb,   0,   0,   0,   0,   0,   0,   0,   0,     0,   0,   0,  +1,   0,   0,   0,   0;  % KCL 3
             0,   0, -Gb,  Gb,   0,   0,   0,   0,   0,   0,   0,   0,    +1,   0,   0,   0,   0,   0,   0,   0;  % KCL 4
             0,   0,   0,   0,  Gc, -Gc,   0,   0,   0,   0,   0,   0,     0,   0,   0,   0,  +1,   0,   0,   0;  % KCL 5
             0,   0,   0,   0, -Gc,  Gc,   0,   0,   0,   0,   0,   0,     0,  +1,   0,   0,   0,   0,   0,   0;  % KCL 6
             0,   0,   0,-gm2,   0,   0,  Gd,   0,   0, -Gd,   0,   0,     0,   0,   0,   0,   0,   0,   0,   0;  % KCL 7
             0, gm1,   0,   0,   0,   0,   0,  Ge, -Ge,   0,   0,-gm1,     0,   0,   0,   0,   0,   0,   0,   0;  % KCL 8
             0,   0,   0,   0,   0,   0,   0, -Ge,  Ge,   0,   0,   0,     0,   0,   0,   0,   0,   0,  +1,   0;  % KCL 9
             0,   0,   0,   0,   0,   0, -Gd,   0,   0,  Gd,   0,   0,     0,   0,   0,   0,   0,  +1,   0,   0;  % KCL 10
             0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  Gf,   0,     0,   0,   0,   0,   0,   0,   0,  +1;  % KCL 11
             0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,     0,   0,   0,   0,   0,   0,   0,  -1;  % KCL 12
        %---------------------------------------------------------------------------------------------------------
             0,   0,   0,  -1,   0,   0,   0,  +1,   0,   0,   0,   0,     0,   0,   0,   0,   0,   0,   0,   0;  % nullator 1
             0,   0,   0,   0,   0,  -1,  +1,   0,   0,   0,   0,   0,     0,   0,   0,   0,   0,   0,   0,   0;  % nullator 2
            +1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,     0,   0,   0,   0,   0,   0,   0,   0;  % port eA
             0,  -1,  +1,   0,   0,   0,   0,   0,   0,   0,   0,   0,     0,   0,   0,   0,   0,   0,   0,   0;  % port eB
             0,  -1,   0,   0,  +1,   0,   0,   0,   0,   0,   0,   0,     0,   0,   0,   0,   0,   0,   0,   0;  % port eC
             0,   0,   0,   0,   0,   0,   0,   0,   0,  +1,   0,   0,     0,   0,   0,   0,   0,   0,   0,   0;  % port eD
             0,   0,   0,   0,   0,   0,   0,   0,  +1,   0,   0,   0,     0,   0,   0,   0,   0,   0,   0,   0;  % port eE
             0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  +1,  -1,     0,   0,   0,   0,   0,   0,   0,   0]; % port eF

    %% R-node scattering matrix
    % S = I + 2*A1*X^-1*A2', where A1 selects the port-current unknowns
    % scaled by the port impedances and A2 selects them unscaled.
    port_Z = [Za Zb Zc Zd Ze Zf];
    n_cols = size(X, 1);                     % 20 unknowns
    A1 = [zeros(6, n_cols-6), diag(port_Z)];
    A2 = [zeros(6, n_cols-6), eye(6)];

    S_rnode = eye(6) + 2*A1*(X\A2');

    %% Wave Digital Filter processing (sample by sample)
    a = zeros(6, 1);                         % incident waves at the R-node
    vLpf = zeros(Nsamp, 1);                  % output across capacitor D (low-pass)
    vBpf = zeros(Nsamp, 1);                  % output across capacitor E (band-pass)

    for ii = 1:Nsamp
        % Forward scan: apply the input source at port F
        a(6) = Vin(ii);

        % Backward scan: reflected waves from the R-node
        b = S_rnode*a;

        % Capacitor ports: the reflected wave becomes the next incident
        % wave (unit delay)
        a(4) = b(4);
        a(5) = b(5);

        % Port voltages v = (a+b)/2; sign flipped to match the circuit output
        vLpf(ii) = -(a(4) + b(4))/2;
        vBpf(ii) = -(a(5) + b(5))/2;
    end
    %% System identification of the band-pass transfer function
    % Fit a 2-pole transfer function to the simulated step/impulse response
    % and plot its magnitude response. To identify the low-pass output
    % instead, replace vBpf with vLpf.
    np = 2;                                  % number of poles
    sys_id  = tfest(iddata(vLpf, Vin, Ts), np);
    LPF_tf  = tf(sys_id.Numerator, sys_id.Denominator);
    bodeplot(LPF_tf, bode_opts);
    LPF_tf_all{jj} = LPF_tf;                 % keep for the LTspice comparison
    %% System identification of the band-pass transfer function
    % % Fit a 2-pole transfer function to the simulated step/impulse response
    % % and plot its magnitude response. To identify the low-pass output
    % % instead, replace vBpf with vLpf.
    % np = 2;                                  % number of poles
    % sys_id  = tfest(iddata(vBpf, Vin, Ts), np);
    % BPF_tf  = tf(sys_id.Numerator, sys_id.Denominator);
    % bodeplot(BPF_tf, bode_opts);
end
toc

grid on;
legend(arrayfun(@(s) sprintf('slope = %d', s), slope_setting, 'UniformOutput', false), ...
       'Location', 'best');

%% Comparison with LTspice
% Raw LTspice export: one block per "Slope=" step, each line "f  (mag dB, phase deg)"
ltspice_file = 'ltspice_raw_tf.txt';
txt    = fileread(ltspice_file, 'Encoding', 'ISO-8859-1');   % file contains the degree symbol
blocks = regexp(txt, 'Slope=(\d+)[^\n]*\n(.*?)(?=Step Information|$)', 'tokens');

figure;
hold on;
% colors = lines(n_iter);
h_sim  = gobjects(1, n_iter);
h_ltsp = gobjects(1, n_iter);
for jj = 1:n_iter
    % LTspice data of the matching slope step
    ltspice_slope = cellfun(@(t) str2double(t{1}), blocks);
    kk  = find(ltspice_slope == slope_setting(jj), 1);
    pts = regexp(blocks{kk}{2}, '(\S+)\t\((\S+)dB', 'tokens');
    pts = str2double(vertcat(pts{:}));
    f_ltspice   = pts(:,1);                  % [Hz]
    mag_ltspice = pts(:,2);                  % [dB]

    % Simulated (identified) low-pass response on the same frequency axis
    H_sim = squeeze(freqresp(LPF_tf_all{jj}, 2*pi*f_ltspice));

    
    h_sim(jj)  = semilogx(f_ltspice, 20*log10(abs(H_sim)), 'r', 'LineWidth', 1.5);
    h_ltsp(jj) = semilogx(f_ltspice, mag_ltspice,'--b', 'LineWidth', 1.5);
end
set(gca, 'XScale', 'log');
grid on;
xlim([20 2e4]);
ylim([-45 45]);
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
title('WDF simulation (solid)  vs LTspice (dashed)');
% legend(h_sim, arrayfun(@(s) sprintf('slope = %d', s), slope_setting, 'UniformOutput', false), ...
%        'Location', 'best');
legend("LTspice","WDF")
