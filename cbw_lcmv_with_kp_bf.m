
%% calc kp-lcmv
close all; clc; clear all;

%plot consts
plot_deg = true;  % deg/rad
plot_dB = true;   % pow ratio in dB / abs(amplitude)
% plot_2d = true;  % plot 2d grid (theta, f axes)
plot_2d = false;  % or plot 1d grid (theta axis only)

% model params
theta_cbw = deg2rad(20);
theta_d = deg2rad(0);
% theta_nulls = [deg2rad(40)];
theta_nulls = [deg2rad(40) deg2rad(70)];
% theta_nulls = [deg2rad(40) deg2rad(-30)];

delta=0.035; % spatial sampling distance
M1=11; % num of mics in wider virtual array
M2=1+length(theta_nulls); % num of mics in narrower virtual array
M=M1+M2-1; % num of mics in the global array

f = [4500];
% f = [3000 4500 6000];
% f = linspace(2500,8000,800);
% f = linspace(1,8000,100);
% f = linspace(1,8000,800);
% f = linspace(1700,8000,100);

[w, w1, w2] = calc_kp_lcmv(M1, delta, f, theta_cbw, theta_d, theta_nulls);

plot_beampattern(w1,f,M1,delta,plot_deg,plot_dB,plot_2d);
plot_beampattern(w2,f,M2,delta,plot_deg,plot_dB,plot_2d);

plot_beampattern(w,f,M,delta,plot_deg,plot_dB,plot_2d);


%% calc lcmv

% theta_nulls = [deg2rad(50) deg2rad(30)];
% f = linspace(1,8000,800);
noise_field = "white";

w_lcmv = calc_lcmv(M, delta, f, theta_d, theta_nulls, noise_field);
plot_2d_temp = false;
plot_beampattern(w_lcmv,f,M,delta,plot_deg,plot_dB,plot_2d_temp);


%% compare kp-lcmv to lcmv
f = [f(1) f(1)];

w_temp = [w_lcmv w];
plot_2d_temp = false;
plot_beampattern(w_temp,f,M,delta,plot_deg,plot_dB,plot_2d_temp);

%% WNG

%plot consts
linewd = 1;
hcfontsize = 20;
is_plot = false;

% f = linspace(1,8000,500);
f = linspace(1,8000,100);
noise_field = "white";

w = calc_kp_lcmv(M1, delta, f, theta_cbw, theta_d, theta_nulls);
w_lcmv = calc_lcmv(M, delta, f, theta_d, theta_nulls, noise_field);
w_ds = ones([M length(f)])/M;

snr_gain_dB = calc_snr_gain(w, f, delta, theta_d, noise_field, is_plot);
snr_gain_dB_lcmv = calc_snr_gain(w_lcmv, f, delta, theta_d, noise_field, is_plot);
snr_gain_dB_ds = calc_snr_gain(w_ds, f, delta, theta_d, noise_field, is_plot);

figure
plot(f, snr_gain_dB,'linewidth',linewd);
hold on;

plot(f, snr_gain_dB_lcmv,'linewidth',linewd);
plot(f, snr_gain_dB_ds,'linewidth',linewd);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;


xlabel('f [Hz]');
ylabel('SNR Gain [dB]');
lgd = legend('kp-lcmv', 'lcmv', 'ds');
lgd.FontSize = 28;


%% DF

%plot consts
linewd = 1;
hcfontsize = 20;

% f = linspace(1,8000,500);
noise_field = "diffuse";
% w_mdf = ones([M length(f)])/M;
% f = linspace(1,8000,100);
% f = linspace(1,8000,100);
is_plot = false;

snr_gain_dB = calc_snr_gain(w, f, delta, theta_d, noise_field, is_plot);
snr_gain_dB_lcmv = calc_snr_gain(w_lcmv, f, delta, theta_d, noise_field, is_plot);
% snr_gain_dB_ds = calc_snr_gain(w_ds, f, delta2, theta_d, noise_field, is_plot);

figure
plot(f, snr_gain_dB,'linewidth',linewd);
hold on;

plot(f, snr_gain_dB_lcmv,'linewidth',linewd);
% plot(f, snr_gain_dB_ds,'linewidth',linewd);
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;


xlabel('f [Hz]');
ylabel('SNR Gain [dB]');
lgd = legend('kp-lcmv', 'null steering bf');
lgd.FontSize = 28;


