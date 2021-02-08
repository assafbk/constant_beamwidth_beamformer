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
% theta_nulls = [deg2rad(40) deg2rad(70)];
theta_nulls = [deg2rad(40) deg2rad(-30)];

% f = [4500];
% f = [3000 4500 6000];
% f = linspace(2500,8000,800);
% f = linspace(1,8000,100);
f = linspace(1700,8000,100);

c=340; % speed of sound
M1=11; % num of mics in wider virtual array
M2=1+length(theta_nulls); % num of mics in narrower virtual array
M=M1*M2; % num of mics in the global array
delta1=0.035; % spatial sampling distance of first virtual array
delta2=delta1/M2; % spatial sampling distance [in meters]

w1 = zeros([M1 length(f)]);
w2 = zeros([M2 length(f)]);
w = zeros([M length(f)]);

if rem(M2,2) == 1  % init for virtual array 2
    N=(M2-1)/2;
    m = (-N:N)';
else
    N=M2/2;
    m = (-(N-1):N)';
end

for i=1:length(f)
    [w1(:,i), beta] = calc_kaiser(theta_cbw, f(i), M1, delta1);

    d_d = exp(-1j*(2*pi*f(i)*delta2*sin(theta_d)/c)*m);
    C = [d_d];
    i_c = [1];
    for l=1:length(theta_nulls)
        C = [C exp(-1j*(2*pi*f(i)*delta2*sin(theta_nulls(l))/c)*m)];
        i_c = [i_c; 0];
    end
    w2(:,i) = C*inv(C'*C)*i_c;

    w(:,i) = kron(w1(:,i),w2(:,i));
end

% plot_beampattern(w1,f,M1,delta1,plot_deg,plot_dB,plot_2d);
% plot_beampattern(w2,f,M2,delta2,plot_deg,plot_dB,plot_2d);

% plot_beampattern(w,f,M,delta2,plot_deg,plot_dB,plot_2d);


%% compare to lcmv
f = [f(1) f(1)];

if rem(M,2) == 1  % init for lcmv
    N=(M-1)/2;
    m = (-N:N)';
else
    N=M/2;
    m = (-(N-1):N)';
end

d_d = exp(-1j*(2*pi*f(1)*delta2*sin(theta_d)/c)*m);
C = [d_d];
i_c = [1];
for l=1:length(theta_nulls)
    C = [C exp(-1j*(2*pi*f(1)*delta2*sin(theta_nulls(l))/c)*m)];
    i_c = [i_c; 0];
end

% white noise field
%     sigma_w = 1;
%     Phi_w = sigma_w*eye(M);

% diffuse noise field
[I,J] = meshgrid(1:M,1:M);
Phi_w = sinc(2*pi*f(1)*(J-I)*delta2/c);

% some other noise field
% alpha=0.8;
% Phi_w = alpha^abs(J-I);

Phi_y = Phi_w;

w_lcmv = inv(Phi_y)*C*inv(C'*inv(Phi_y)*C)*i_c;

w_temp = [w_lcmv w];
plot_2d_temp = false;
plot_beampattern(w_temp,f,M,delta2,plot_deg,plot_dB,plot_2d_temp);


%% lcmv

% f = linspace(1,8000,800);
theta_nulls = [deg2rad(40) deg2rad(-30)];
f = linspace(1,8000,800);

w_lcmv = zeros([M length(f)]);

if rem(M,2) == 1  % init for lcmv
    N=(M-1)/2;
    m = (-N:N)';
else
    N=M/2;
    m = (-(N-1):N)';
end

for i=1:length(f)
    d_d = exp(-1j*(2*pi*f(i)*delta2*sin(theta_d)/c)*m);
    C = [d_d];
    i_c = [1];
    for l=1:length(theta_nulls)
        C = [C exp(-1j*(2*pi*f(i)*delta2*sin(theta_nulls(l))/c)*m)];
        i_c = [i_c; 0];
    end
    
    
    % white noise field
    sigma_w = 1;
    Phi_w = sigma_w*eye(M);

    % diffuse noise field
%     [I,J] = meshgrid(1:M,1:M);
%     Phi_w = sinc(2*pi*f(i)*(J-I)*delta2/c);

    % some other noise field
%     alpha=0.8;
%     Phi_w = alpha^abs(J-I);

    Phi_y = Phi_w;

    w_lcmv(:,i) = inv(Phi_y)*C*inv(C'*inv(Phi_y)*C)*i_c;
end

plot_2d_temp = true;
plot_beampattern(w_lcmv,f,M,delta2,plot_deg,plot_dB,plot_2d_temp);


%% WNG

%plot consts
linewd = 1;
hcfontsize = 20;

% f = linspace(1,8000,500);
noise_field = "white";
w_ds = ones([M length(f)])/M;
% f = linspace(1,8000,100);
% f = linspace(1,8000,100);
is_plot = false;

snr_gain_dB = calc_snr_gain(w, f, delta2, theta_d, noise_field, is_plot);
snr_gain_dB_lcmv = calc_snr_gain(w_lcmv, f, delta2, theta_d, noise_field, is_plot);
snr_gain_dB_ds = calc_snr_gain(w_ds, f, delta2, theta_d, noise_field, is_plot);

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

snr_gain_dB = calc_snr_gain(w, f, delta2, theta_d, noise_field, is_plot);
snr_gain_dB_lcmv = calc_snr_gain(w_lcmv, f, delta2, theta_d, noise_field, is_plot);
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
lgd = legend('KP-LCMV', 'Null Steering BF');
lgd.FontSize = 28;


