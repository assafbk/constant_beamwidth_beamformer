close all; clc; clear all;

%plot consts
plot_deg = true;  % deg/rad
plot_dB = true;   % pow ratio in dB / abs(amplitude)
% plot_2d = true;  % plot 2d grid (theta, f axes)
plot_2d = false;  % or plot 1d grid (theta axis only)

% model params
delta=0.035; % spatial sampling distance
c=340; % speed of sound
M=11; % num of mics in array
N=(M-1)/2;
theta_cbw = deg2rad(20);

% Mc = 3; % num of constraints
Mc = 2; % num of constraints
theta_d = deg2rad(0);
theta_1 = deg2rad(40);
% theta_2 = deg2rad(70);
% constraint_thetas = [theta_d, theta_1 theta_2];
constraint_thetas = [theta_d, theta_1];

% f = [3000 4500 6000];
% f = [4000 5000 6000];
% f = linspace(3000,8000,60);
% f = [4500 6000 7500];
f = [6000];

m = (-N:N)';
w_gsc = zeros([M length(f)]);

for i=1:length(f)
    
    d_d = exp(-1j*(2*pi*f(i)*delta*sin(theta_d)/c)*m);
    d_1 = exp(-1j*(2*pi*f(i)*delta*sin(theta_1)/c)*m);
%     d_2 = exp(-1j*(2*pi*f(i)*delta*sin(theta_2)/c)*m);

    % assuming the noise field is a diffuse noise field
    % zeroing all constraint corr matrices, this causes Phi_y = Phi_w
    % until now we get the same results
%     sigma_x = 0;    % source energy  % FIXME think about this value
%     sigma_u = 0;    % angular interference of white gaussian noise
    % sigma_w = 0.01;  % thermal noise energy
%     T = 512;        % sampling interval
%     Phi_x = T*sigma_x*d_d*(d_d'); 
%     Phi_u = T*sigma_u*d_1*(d_1');
%     Phi_u = T*sigma_u*d_1*d_1' + T*sigma_u*d_2*d_2';

    % white noise field
    % Phi_w = T*sigma_w*eye(M);

    % diffuse noise field
    [I,J] = meshgrid(1:M,1:M);
    Phi_w = sinc(2*pi*f(i)*(J-I)*delta/c);
    
    % some other noise field
    % alpha=0.8;
    % Phi_w = alpha^abs(J-I);

%     Phi_v = Phi_u + Phi_w;
%     Phi_y = Phi_x + Phi_v;
    Phi_y = Phi_w;


%     i_c = [1 0 0]';
%     C = [d_d d_1 d_2];
    i_c = [1 0]';
    C = [d_d d_1];
    B = (eye(M) - C*inv(C'*C)*C') * eye(M,M-Mc);

%     w_fbf = C*inv(C'*C)*i_c;  % MN
%     w_fbf = calc_kaiser(theta_cbw, f(i), M);
    attn_at_constraint_dB = -70;
    w_fbf = calc_constrained_dpss(constraint_thetas, theta_cbw, attn_at_constraint_dB, f(i), M);
    w_ad = inv(B'*Phi_y*B)*B'*Phi_y*w_fbf;
    w_gsc(:,i) = conj(w_fbf) - (w_ad'*B').';
    
end

% plot_beampattern(w_gsc,f,M,plot_deg,plot_dB,plot_2d);



% plot all parts of gsc
% w_gsc_temp = [conj(w_fbf) w_gsc (w_ad'*B').'];
% f = [6000 6000 6000];
% plot_beampattern(w_gsc_temp,f,M,plot_deg,true,plot_2d);


% plot parts of gsc
w_gsc_temp = [conj(w_fbf) w_gsc];
f = [6000 6000];
plot_beampattern(w_gsc_temp,f,M,plot_deg,plot_dB,plot_2d);



%% compare to lcmv (should be identical)
w_lcmv = conj(inv(Phi_y)*C*inv(C'*inv(Phi_y)*C)*i_c);
f = [4000 4000];
w_temp = [w_lcmv w_gsc];
plot_beampattern(w_temp,f,M,plot_deg,plot_dB,plot_2d);





