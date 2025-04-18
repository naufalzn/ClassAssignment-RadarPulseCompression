clc;clear all;close all
%% LFM Input Parameter
fs=1e3 %Sampling rate (Hz)
f_low=0 %Chirp start frequency (Hz)
f_high=50 %Chirp end frequency (Hz)
t_on = 1 %Chirp duration (s)
t_max=15 %Echo received duration in one pulse (s)
A = 1 %Amplitude
fs_low = f_low/fs;
fs_high = f_high/fs;
tt= 0:1/fs:t_max-1/fs 

%% LFM Echo Simulated Parameter
t0 = 4; %Echo first time received 
dt = 1; %Period of received echo (s)
lgs = t0:dt:t_max;
att = 1.1; %Attenuation factor
ref = 0.2; %Echo amplitude

%% Generate LFM
N=fs*t_on
d_fs = fs_high-fs_low;
t=0:N-1
lfm_pw = A*cos((t.*(pi*d_fs)/(N-1)+(2*pi*fs_low)).*t) + j*A*sin((t.*(pi*d_fs)/(N-1)+(2*pi*fs_low)).*t);
t_lfm_pw = 0:1/fs:t_on-1/fs

%% Add AWGN Noise to LFM
SNR = 15 %SNR of AWGN
rpls = pulstran(tt,[lgs;ref*att.^-(lgs-t0)]',lfm_pw,fs); %Make echo received periodic
r=randn(size(tt))*std(lfm_pw)/db2mag(SNR)
rplsnoise=r+rpls %Add noise to lfm

%% Cross Correlation
[c,lags]=xcorr(rplsnoise,lfm_pw) %Cross correlation of received signal with transmitted signal
y_xcorr=c(lags>=0)
t_xcorr=lags(lags>=0)/fs;

%% Convolution
pls_rev = fliplr(lfm_pw) %Pulse flipped
h=conj(pls_rev) %Pulse conjugated
y_conv=conv(rplsnoise, h,'same') %Convolution of received signal with flipped conjugate transmitted signal

%% Plot
figure(1)
plot(tt,rplsnoise)  
hold on;
plot(t_lfm_pw, lfm_pw)
hold on
plot(tt, rpls)
title('Transmitted and Received Signal')
xlabel('Time (s)')
ylabel('Magnitude')
legend('Noise signal', 'Transmitted pulse', 'Noiseless echo received')

figure(2)
plot(t_xcorr,y_xcorr)
hold on
plot(t_xcorr,y_conv)
title('Pulse Compression')
xlabel('Time (s)')
ylabel('Magnitude')
legend('Autocorrelation', 'Convolution')



