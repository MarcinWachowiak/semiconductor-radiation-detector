%Moving average filter simulation
%with ADC quantization, tested on noisy expotential pulses
clear all; close all;

%ADC parameters
Fs = 100*1e+6;             %[MHz] Sampling frequency
N_bits = 14;               %ADC resolution
ADC_res = 2^N_bits-1;
ADC_V_max = 2;                  
ADC_quant = ADC_V_max/ADC_res; %ADC step

%Simulation preset
T_sim = 3000*1e-6;      %[uS]   %time of simulation
A_pulse =2.5*1e-3;      %[mV]   %amplitude of the expotential pulse
A_noise =10*1e-3;       %[mV]   %amplitude of the noise
SNR=(A_pulse/A_noise)^2;        %signal to noise ratio
T_pulse = 500*1e-6;     %[uS]   %interval between pulses   
N_pulses = 5;                   %pulse count
        
% TIA values
C= 1e-10;                   %TIA feedback capacitance
R= 68*1e+6 ;                %TIA feedback resistance [Mohms]
Corr_coef =0.007;           %Correction coeficient based on real measurments
RC_tau=R*C*Corr_coef;       %RC constant

%Input signal processing
T_ext =(T_sim-(T_pulse*N_pulses))/2 ;
dt = 1/Fs;
Sig_time = (0:dt:T_pulse*N_pulses-dt);
Pre_ext = (0:dt:T_ext-dt);
Post_ext =(T_sim-T_ext:dt:T_sim);
T_sim_base = (0:dt:T_sim);
Exp_pulses = A_pulse/(1-exp(-T_pulse/RC_tau))*(1*exp(-mod(Sig_time,T_pulse)/RC_tau) - exp(-T_pulse/RC_tau));
Pre_ext_sig = zeros(1,length(Pre_ext));
Post_ext_sig = zeros(1,length(Post_ext));
Cln_sig =cat(2,cat(2,Pre_ext_sig,Exp_pulses),Post_ext_sig);

w_g_n = randn(1,length(T_sim_base));
Wgn_noise = A_noise*w_g_n/9;        %correction of randn generator range 
WgnMIN = min(Wgn_noise);
WgnMAX = max(Wgn_noise);
WgnAmpl = WgnMAX-WgnMIN;
Wgn_nois_mean = mean(Wgn_noise);
Combined_signal = Cln_sig+Wgn_noise;
H_factor = 170;
Shift_factor = 0;
Processed_sig = Combined_signal*H_factor+Shift_factor;
Cln_q_sig = floor(Cln_sig/ADC_quant);
Quantized_sig = floor(Processed_sig/ADC_quant);

figure
subplot(3,1,1); grid on;plot(T_sim_base*1e+6,Cln_sig*1e+3);title('Input signal without noise');
xlabel('Time [\mus]');ylabel('Amplitude [mV]');xlim([0 T_sim*1e+6]); 
subplot(3,1,2); grid on;plot(T_sim_base*1e+6,Combined_signal*1e+3);title('Input signal with noise');
xlabel('Time [\mus]');ylabel('Amplitude [mV]');xlim([0 T_sim*1e+6]); 
subplot(3,1,3); grid on;plot(T_sim_base*1e+6,Processed_sig);title('Amplified and shifted signal');
xlabel('Time [\mus]');ylabel('Amplitude [V]');xlim([0 T_sim*1e+6]); 

%Moving average filters parameters
MAV_window_points_1=100;    %Filter buffer length  
MAV_window_points_2=1000;
MAV_window_points_3=10000;

Buffer_1 = ones(1,MAV_window_points_1)/MAV_window_points_1;
Buffer_2= ones(1,MAV_window_points_2)/MAV_window_points_2;
Buffer_3 = ones(1,MAV_window_points_3)/MAV_window_points_3;

MAV_filtered_sig_1=filter(Buffer_1,1,Quantized_sig);
MAV_filtered_sig_2=filter(Buffer_2,1,Quantized_sig);
MAV_filtered_sig_3=filter(Buffer_3,1,Quantized_sig);

MAV_quantized_filtered_sig_1 = floor(MAV_filtered_sig_1);
MAV_quantized_filtered_sig_2 = floor(MAV_filtered_sig_2);
MAV_quantized_filtered_sig_3 = floor(MAV_filtered_sig_3);

figure
G_title = sprintf('MAV filter output signal');
G1_title = sprintf('Sampling freqency: %d MSPS',Fs/1e+6 );
G2_title = sprintf('Buffer size: %d',MAV_window_points_1);
subplot(3,1,1); grid on;plot(T_sim_base*1e+6,MAV_quantized_filtered_sig_1*ADC_quant*10e+2,'.');title( {G_title, G1_title, G2_title});
xlabel('Time [\mus]');ylabel('Amplitude [mV]');xlim([0 T_sim*1e+6]);
subplot(3,1,2); grid on;plot(T_sim_base*1e+6,MAV_quantized_filtered_sig_2*ADC_quant*10e+2,'.');title( ['Buffer size: ' num2str( MAV_window_points_2 )] );
xlabel('Time [\mus]');ylabel('Amplitude [mV]');xlim([0 T_sim*1e+6]);
subplot(3,1,3); grid on;plot(T_sim_base*1e+6,MAV_quantized_filtered_sig_3*ADC_quant*10e+2,'.');title( ['Buffer size: ' num2str( MAV_window_points_3 )] );
xlabel('Time [\mus]');ylabel('Amplitude [mV]');xlim([0 T_sim*1e+6]);

%Low-pass filters parameters
fs_multiplier = 8; %determines local sampling rate after decimation (multiple of fpass)
% 10kHz
f_1_pass = 10*1e+3;
fs_1 = fs_multiplier*f_1_pass;
r_1 =ceil(Fs/fs_1);
Y_Fs_10kHz = decimate(Quantized_sig,r_1,'fir');
Lowpass_10kHz_filtered_sig = lowpass(Y_Fs_10kHz,f_1_pass,fs_1);
% 50kHz
f_2_pass = 50*1e+3;
fs_2 = fs_multiplier*f_2_pass;
r_2 =ceil(Fs/fs_2);
Y_Fs_50kHz = decimate(Quantized_sig,r_2,'fir');
Lowpass_50kHz_filtered_sig = lowpass(Y_Fs_50kHz,f_2_pass,fs_2);
% 100kHz
f_3_pass = 100*1e+3;
fs_3 = fs_multiplier*f_3_pass;
r_3 =ceil(Fs/fs_3);
Y_Fs_100kHz = decimate(Quantized_sig,r_3,'fir');
Lowpass_100kHz_filtered_sig = lowpass(Y_Fs_100kHz,f_3_pass,fs_3);
% 500kHz
f_4_pass = 500*1e+3;
fs_bandpass = fs_multiplier*f_4_pass;
r_4 =ceil(Fs/fs_bandpass);
Y_Fs_Bandpass = decimate(Quantized_sig,r_4,'fir');
Lowpass_500kHz_filtered_sig = lowpass(Y_Fs_Bandpass,f_4_pass,fs_bandpass);

figure
G_title = sprintf('Digital low-pass filter output signal');
G1_title = sprintf('Fpass: %d [kHz]',f_1_pass/(1*1e+3));
subplot(4,1,1); grid on;plot(Lowpass_10kHz_filtered_sig*ADC_quant*1e+3,'.');title( {G_title, G1_title});
xlabel('Sample');ylabel('Amplitude [mV]');xlim([0 length(Lowpass_10kHz_filtered_sig)]); 
subplot(4,1,2); grid on;plot(Lowpass_50kHz_filtered_sig*ADC_quant*1e+3,'.');title( ['Fpass:: ' num2str( f_2_pass/(1*1e+3) ) ' [KHz]'] );
xlabel('Sample');ylabel('Amplitude [mV]');xlim([0 length(Lowpass_50kHz_filtered_sig)]); 
subplot(4,1,3); grid on;plot(Lowpass_100kHz_filtered_sig*ADC_quant*1e+3,'.');title( ['Fpass:: ' num2str( f_3_pass/(1*1e+3) ) ' [KHz]'] );
xlabel('Sample');ylabel('Amplitude [mV]');xlim([0 length(Lowpass_100kHz_filtered_sig)]); 
subplot(4,1,4); grid on;plot(Lowpass_500kHz_filtered_sig*ADC_quant*1e+3,'.');title( ['Fpass:: ' num2str( f_4_pass/(1*1e+3) ) ' [KHz]'] );
xlabel('Sample');ylabel('Amplitude [mV]');xlim([0 length(Lowpass_500kHz_filtered_sig)]); 

% % Band-pass filtering test
% % 10kHz - 200kHz
% f_1_bandpass = 10*1e+3;
% f_2_bandpass = 200*1e+3;
% fs_bandpass = fs_multiplier*f_2_bandpass;
% r_5 =ceil(Fs/fs_bandpass);
% Y_Fs_Bandpass = decimate(Quantized_sig,r_4,'fir');
% Bandpass_filtered_sig = bandpass(Y_Fs_Bandpass,[f_1_bandpass f_2_bandpass],fs_bandpass);
% 
% figure
% G_title = sprintf('Digital band-pass filter output signal');
% G1_title = sprintf('Fpass: %d-%d [kHz]',f_1_bandpass/(1*1e+3),f_2_bandpass/(1*1e+3));
% plot(Bandpass_filtered_sig*ADC_quant*1e+3,'.');title( {G_title, G1_title});grid on;
% xlabel('Sample');ylabel('Amplitude [mV]');xlim([0 length(Bandpass_filtered_sig)]); 
