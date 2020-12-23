%% Description
% The code of WBM system:
% 1) The script is for the case that the communication system includes
% the central AP and one node.
% 2) The central AP, as the receiver, utilize two beams to aim at two
% trasnsmitting beams of a sigle node.


%%
clc;
close all;
clear;

%% Set random
seed = 1;
rng(seed);

%% Set Noise
noise_switch = 1;

%% Set SNR
SNR_dB = -30:1:0; 

%% Generate the packet with known preamble bits
Ns = 1;
sym_bit = 1;
N_sym_per_frame = Ns * sym_bit;
N_frame = 10000;
N_packet = 10;
N_bits_per_packet = N_sym_per_frame * N_frame;


pre_bits = [1 0 0 1 0 0 1 1 0 1 1 0 0 1 0 0 1 1 0 1 1 0 0 1 0 0 1 1 0 1 1 0 0 1 0 0 1 1 0 1 ...
    1 0 0 1 0 0 1 1 0 1 1 0 0 1 0 0 1 1 0 1 1 0 0 1 0 0 1 1 0 1 1 0 0 1 0 0 1 1 0 1];
% pre_bits = [1 0 0 1 0 0 1 1 0 1];


%% Set the communication frequency
Ts = 1 / (200 * 1e2);  
f = 1 / Ts; 
fc = 26 * 1e9; 
fs = 4 * fc; 
Ts_reso = fs / f; 

%% Generate time-domain wave of the one sym
A_t = ones(1, Ts_reso);
S_t = modulate(A_t,fc,fs, 'am');

% figure();
% plot(S_t);

%% Set the system parameter
N_beam = 2; 
Nt_antenna_group = 4;  
Nr_antenna_group = 16;

%% =============== I: Generate the channel matrix =======================
% Step1: Generate the channel parameter
N_c = N_beam + unidrnd(3); 
N_ray = 1; 

aod_range_max = 120;
aod_deg_H = rand(N_c, N_ray) * aod_range_max - aod_range_max / 2;
aoa_deg_H = rand(N_c, N_ray) * aod_range_max - aod_range_max / 2;
aod_rad_H = aod_deg_H / 180 * pi;
aoa_rad_H = aoa_deg_H / 180 * pi;
alpha = rand(N_c, N_ray);

% Step2: Generate the channel
[H_node, Frf_node, Qrf_node, theta_aod, theta_aoa] = channel_generation_angle(...
Nt_antenna_group, Nr_antenna_group, N_c, N_ray, aod_rad_H, aoa_rad_H, alpha);

% Step3: Generate the ARV
theta_aod = [30; -30] / 180 * pi;
Frf_vec = steer_vector(Nt_antenna_group, theta_aod);
Qrf_vec = steer_vector(Nr_antenna_group, theta_aoa);

% %  Ploting the beampattern
% beampattern(Frf_node);
% beampattern(Qrf_node);
% beampattern(Frf_vec);
% beampattern(Qrf_vec);


%% ======================== Ergodic SNR ===========================
N_node = 1;
ber_num = zeros(length(SNR_dB), N_packet);   
for index_snr = 1: length(SNR_dB)
    SNR = SNR_dB(index_snr);
    [Nt_antenna_group, SNR]

	%% =================== beam alignment =========================
    % Step1: Get SNR and generate the normalized noise
    Pn = 1 / 10^(SNR / 10); 
    u_data = sqrt(Pn/2).*(randn(Nr_antenna_group,Ts_reso * length(pre_bits))+...
        1i*randn(Nr_antenna_group,Ts_reso * length(pre_bits)));
    u_prebits_P = mean(mean(abs(u_data).^2, 2));
    
    % Step2: Generate the normalized signal
    bit0_t = kron(1 - pre_bits, S_t);
    beam0_s_t = Frf_vec(:, 1) * bit0_t;
    beam0_s_t_P = sum(2 * mean(abs(beam0_s_t).^2, 2));
    beam0_s_t_norm = beam0_s_t / sqrt(beam0_s_t_P);
    beam0_s_t_norm_P = sum(mean(abs(beam0_s_t_norm).^2, 2));

    bit1_t = kron(pre_bits, S_t);
    beam1_s_t = Frf_vec(:, 2) * bit1_t;
    beam1_s_t_P = sum(2 * mean(abs(beam1_s_t).^2, 2));
    beam1_s_t_norm = beam1_s_t / sqrt(beam1_s_t_P);
    beam1_s_t_norm_P = sum(mean(abs(beam1_s_t_norm).^2, 2));

    beam0_dB = 10 * log10(beam0_s_t_norm_P / u_prebits_P);
    beam1_dB = 10 * log10(beam1_s_t_norm_P / u_prebits_P);
    beam_s_t_norm_P = beam0_s_t_norm_P + beam1_s_t_norm_P;
    beam_dB = 10 * log10(beam_s_t_norm_P / u_prebits_P);
    
    % Step3: Beam alignment
    recv_sig_beam0 = H_node * beam0_s_t_norm + noise_switch * u_data;
    recv_sig_beam1 = H_node * beam1_s_t_norm + noise_switch * u_data;
    Wrf_beam0 = beam_alignment(recv_sig_beam0, Nr_antenna_group);
    Wrf_beam1 = beam_alignment(recv_sig_beam1, Nr_antenna_group);
    test = 0;
    
%     beampattern(Wrf_beam0);
%     beampattern(Wrf_beam1);
    
    % Step4: Get the received noise and interfernce 
    Wrf = [Wrf_beam0 Wrf_beam1];
	y_tb0 = Wrf' * (H_node * beam0_s_t_norm + noise_switch * u_data);
    y_tb1 = Wrf' * (H_node * beam1_s_t_norm + noise_switch * u_data);
    
    % Demodulate the received signal
    y_F0_W0 = demod(y_tb0(1, :), fc, fs, 'am');
    y_F0_W1 = demod(y_tb0(2, :), fc, fs, 'am');
	y_F1_W0 = demod(y_tb1(1, :), fc, fs, 'am');
    y_F1_W1 = demod(y_tb1(2, :), fc, fs, 'am');
    
    % Average the received signal
    y_F0_W0_energy = 2 * mean(abs(y_F0_W0));
    y_F0_W1_energy = 2 * mean(abs(y_F0_W1));
    y_F1_W0_energy = 2 * mean(abs(y_F1_W0));
    y_F1_W1_energy = 2 * mean(abs(y_F1_W1));
    test = 0;
    
    % Get the decision
    y_beam0 = y_F0_W1_energy - y_F0_W0_energy;
    y_beam1 = y_F1_W1_energy - y_F1_W0_energy;
    threshold = y_beam1 + y_beam0;
    node_gap = y_beam1 - y_beam0;
    test = 0;
    

    %% ================== IV: Transmite Dataframe =================
    for packet_index = 1: N_packet
        % Step1: Generate the random sequence
        data_frame = zeros(1, N_bits_per_packet);
        randam_index = randperm(N_bits_per_packet);
        half_point = round(N_bits_per_packet / 2);
        random_high_index = randam_index(1: half_point);
        data_frame(random_high_index) = 1;
        test = sum(data_frame, 2) / length(data_frame);

        % Step2: Get SNR and Generate the normalized noise
        Pn = 1 / 10^(SNR / 10); 
        u_data = sqrt(Pn/2).*(randn(Nr_antenna_group,Ts_reso * length(data_frame))+...
            1i*randn(Nr_antenna_group,Ts_reso * length(data_frame)));
        u_data_P = mean(mean(abs(u_data).^2, 2));

         % Step3: Generate the normalized signal
        bit0_t = kron(1 - data_frame, S_t);
        beam0_s_t = Frf_vec(:, 1) * bit0_t;
        beam0_s_t_P = sum(mean(abs(beam0_s_t).^2, 2));
        beam0_s_t_norm = beam0_s_t / sqrt(beam0_s_t_P);
        beam0_s_t_norm_P = sum(mean(abs(beam0_s_t_norm).^2, 2));

        bit1_t = kron(data_frame, S_t);
        beam1_s_t = Frf_vec(:, 2) * bit1_t;
        beam1_s_t_P = sum(mean(abs(beam1_s_t).^2, 2));
        beam1_s_t_norm = beam1_s_t / sqrt(beam1_s_t_P);
        beam1_s_t_norm_P = sum(mean(abs(beam1_s_t_norm).^2, 2));

        beam0_dB = 10 * log10(beam0_s_t_norm_P / u_data_P);
        beam1_dB = 10 * log10(beam1_s_t_norm_P / u_data_P);

        % Step4: Get and demodulate the received signal
        beam_s_t = beam1_s_t_norm + beam0_s_t_norm;
        recv_sig_W = Wrf' * (H_node * beam_s_t + noise_switch * u_data);
        recv_sig_W0 = demod(recv_sig_W(1, :), fc, fs, 'am');
        recv_sig_W1 = demod(recv_sig_W(2, :), fc, fs, 'am'); % first minus, second demodulate
        recv_sig_Ts_reso = abs(recv_sig_W1) - abs(recv_sig_W0);
        test = 0;

        % Step5: Judgement
        recv_sig_Ts = reshape(recv_sig_Ts_reso, Ts_reso, []);
        recv_sig = mean(recv_sig_Ts, 1);
        recv_sig_bit_high_index = find(recv_sig >= threshold);
        recv_sig_bit_low_index = find(recv_sig < threshold);
        data_r(recv_sig_bit_high_index) = 1;
        data_r(recv_sig_bit_low_index) = 0;
        % Calculate the ber
        frame_ber_num = biterr(data_r, data_frame);
        ber_num(index_snr, packet_index) = frame_ber_num;
    end
end

ber_rate = sum(ber_num, 2) / (N_packet * length(data_frame));
test = 0;


figure();
semilogy(SNR_dB, ber_rate, 'linewidth', 1.0);
% xlim([min(SNR_dB) max(SNR_dB)]);
ylim([9 * 10^(-5) 1])
grid on;
legend_obj = legend('Nt= 4');
set(legend_obj, 'Fontname', 'Times New Roman','FontSize',10);
xlabel('SNR (dB)', 'Fontname', 'Times New Roman','FontSize',10);
ylabel('BER', 'Fontname', 'Times New Roman','FontSize',10);


