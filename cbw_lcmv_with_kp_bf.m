close all; clc; clear all;

%plot consts
plot_deg = true;  % deg/rad
plot_dB = true;   % pow ratio in dB / abs(amplitude)
plot_2d = true;  % plot 2d grid (theta, f axes)
% plot_2d = false;  % or plot 1d grid (theta axis only)

% model params
theta_cbw = deg2rad(20);
theta_d = deg2rad(0);
theta_nulls = [deg2rad(40) deg2rad(-30)];

% f = [3000 4500 6000];
% f = linspace(2500,8000,800);
% f = linspace(0,8000,800);
f = linspace(1,8000,800);

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

plot_beampattern(w,f,M,delta2,plot_deg,plot_dB,plot_2d);
