close all; clc; clear all;

%plot consts
plot_deg = true;  % deg/rad
plot_dB = true;   % pow ratio in dB / abs(amplitude)
% plot_2d = true;  % plot 2d grid (theta, f axes)
plot_2d = false;  % or plot 1d grid (theta axis only)

% model params
delta=0.035; % spatial sampling distance
c=340; % speed of sound
M=9; % num of mics in array
N=(M-1)/2;
theta_cbw = deg2rad(20);

% Mc = 3; % num of constraints
Mc = 5; % num of constraints
theta_d = deg2rad(0);
theta_1 = deg2rad(40);
theta_2 = deg2rad(70);
constraint_thetas = [theta_d, theta_1 theta_2];
% constraint_thetas = [theta_d, theta_1];

% params for beta search
% delta_beta = 0.2;     % choosed by taking diff of the betas for cbw=20, cbw=19.5
delta_beta = 1;
num_of_angles = 10001; % for theta axis
theta = linspace(0,pi/2,(num_of_angles+1)/2)'; % angles for B(f.theta)
theta_res = theta(2)-theta(1);
theta_tol = 2*theta_res; % allowed tolerance for the mainlobe width around theta_cbw


% f = [3000 4500 6000];
% f = [4000 5000 6000];
% f = linspace(3000,8000,60);
% f = [4500 6000 7500];
f = [4000];

m = (-N:N)';
w_gsc = zeros([M length(f)]);

for i=1:length(f)
    
    d_d = exp(-1j*(2*pi*f(i)*delta*sin(theta_d)/c)*m);
    d_1 = exp(-1j*(2*pi*f(i)*delta*sin(theta_1)/c)*m);
    d_2 = exp(-1j*(2*pi*f(i)*delta*sin(theta_cbw)/c)*m);
    d_3 = exp(-1j*(2*pi*f(i)*delta*sin(-theta_cbw)/c)*m);
    d_4 = exp(-1j*(2*pi*f(i)*delta*sin(theta_2)/c)*m);

%{
% assuming the noise field is a diffuse noise field
      % zeroing all constraint corr matrices, this causes Phi_y = Phi_w
      % until now we get the same results
     sigma_x = 1;    % source energy  % FIXME think about this value
     sigma_u = 1;    % angular interference of white gaussian noise
     sigma_w = 0.01;  % thermal noise energy
     T = 512;        % sampling interval
     Phi_x = T*sigma_x*d_d*(d_d');
     Phi_u = T*sigma_u*d_1*(d_1');
%      Phi_u = T*sigma_u*d_1*d_1' + T*sigma_u*d_2*d_2'; 

    % white noise field
    % Phi_w = T*sigma_w*eye(M);
%}
    % diffuse noise field
    [I,J] = meshgrid(1:M,1:M);
    Phi_w = sinc(2*pi*f(i)*(J-I)*delta/c);
    
    % some other noise field
    % alpha=0.8;
    % Phi_w = alpha^abs(J-I);
%{ 
    Phi_v = Phi_u + Phi_w;
    Phi_y = Phi_x + Phi_v;
%}
%     Phi_cbw = 0.1*abs(J-I);
%     Phi_y = Phi_w+Phi_cbw;
    Phi_y = Phi_w;

%     i_c = [1 0 0]';
%     C = [d_d d_1 d_2];
    i_c = [1 0 0 0 0]';
    C = [d_d d_1 d_2 d_3 d_4];
    B = (eye(M) - C*inv(C'*C)*C') * eye(M,M-Mc);
    
    % initial kaiser FBF
%     [w_fbf, initial_beta] = calc_kaiser(theta_cbw, f(i), M); 
%     beta_ = initial_beta;
    
    % initial cdpss FBF
    attn_at_constraint_dB = -50;
%     attn_at_constraint_dB = -21.14;
    constraint_thetas = [theta_d, theta_1, theta_2, theta_cbw, -theta_cbw];
    w_fbf = calc_constrained_dpss(constraint_thetas, theta_cbw, attn_at_constraint_dB, f(i), M);
    beta_ = 1;  % FIXME remove
    
    % initial D&S FBF
%     w_fbf = ones([M 1])/M;
    
    w_ad = inv(B'*Phi_y*B)*B'*Phi_y*w_fbf;
    w_gsc_initial = conj(w_fbf) - (w_ad'*B').';
    
    m = (-N:N);
    for k=1:100
        % calc gsc with new beta
        w_ad = inv(B'*Phi_y*B)*B'*Phi_y*w_fbf;
        w_gsc(:,i) = conj(w_fbf) - (w_ad'*B').';
        
%         w_gsc_temp = [conj(w_fbf) w_gsc(:,i) (w_ad'*B').'];
%         f_ = [f(1) f(1) f(1)];
        w_gsc_temp = [conj(w_fbf) w_gsc(:,i)];
        f_ = [f(1) f(1)];
        plot_beampattern(w_gsc_temp,f_,M,plot_deg,plot_dB,plot_2d);
        return;
        
        % calc beampattern for 0<theta<pi/2 and find min
        u = 2*pi*f(i)*delta*sin(theta)/c;
        d = exp(-1j*u*m(1,:));
        BP = abs(d*w_gsc(:,i));
        BP(BP<(10^-3)) = 10^-3; % set all low points to same value
        [~, min_BP_idx] = min(BP);  % min_B_idx = first null's index
        theta_first_null = theta(min_BP_idx);
        if abs(theta_first_null-theta_cbw) < theta_tol
            break
        else
            if theta_first_null < theta_cbw
                beta_ = beta_ - delta_beta;
            else
                beta_ = beta_ + delta_beta;
            end
            
            if beta_ < 0
                disp("beta < 0");
                break;
            end
            
            w_fbf = calc_kaiser_with_beta(beta_, M);
            
        end
        
    end
    
end

% plot_beampattern(w_gsc,f,M,plot_deg,plot_dB,plot_2d);



% plot all parts of gsc
w_gsc_temp = [conj(w_fbf) w_gsc (w_ad'*B').'];
f = [f(1) f(1) f(1)];
plot_beampattern(w_gsc_temp,f,M,plot_deg,true,plot_2d);


% plot parts of gsc
% w_gsc_temp = [conj(w_fbf) w_gsc];
% f = [f(1) f(1)];
% plot_beampattern(w_gsc_temp,f,M,plot_deg,plot_dB,plot_2d);



%% compare to lcmv (should be identical)
w_lcmv = conj(inv(Phi_y)*C*inv(C'*inv(Phi_y)*C)*i_c);
f = [4000 4000];
w_temp = [w_lcmv w_gsc];
plot_beampattern(w_temp,f,M,plot_deg,plot_dB,plot_2d);





