close all; clc; clear all;

%consts
% note - we assume theta_d=0
delta=0.035; % spatial sampling distance
c=340; % speed of sound
M=11; % num of mics in array
N=(M-1)/2;

%plot consts
plot_deg = true;  % deg/rad
plot_dB = true;   % dB/pow
plot_2d = true;   % plot 2d grid (theta, f axes)
                  % or plot 1d grid (theta axis only)

% model params
theta_cbw = deg2rad(20); % the angle of the first mainlobe null 
% f = [4000, 6000]; % beampattern freqs
% f = [3000, 4500, 6000]; % beampattern freqs
f = linspace(0,8000,800);
% f=4000;



%% plot delay and sum
 w_ds = (1/M)*ones([M,1]);
 w_ds = w_ds*ones([1 length(f)]); % replicate weights to all freqs
 plot_beampattern(w_ds,f,M,plot_deg,plot_dB,plot_2d);


%% plot modified rectangle

% calc Ks
% for some freq f, m_tilde = 2k-1 is the largest size of the array
% that still has a mainlobe null larger than theta_cbw (will be smoothed to theta_cbw by g)
f_low = c/(M*delta*sin(theta_cbw));
K = zeros([length(f) 1]);
for i=1:length(f)
    k = 1:N;
    k((2*k-1)*f(i) > M*f_low) = -1;
    K(i) = max(k);    
end

% calc windows
w_mod_rect = zeros([M length(f)]);

f_min = c/(delta*M*sin(theta_cbw));
f_max = c/(delta*3*sin(theta_cbw));
for i=1:length(f)
    
    if (f(i) < f_min) | (f(i) > f_max) | (f(i) > c/delta)
        w_mod_rect(:,i)=1/M;
    else
        u = 2*pi*f(i)*delta*sin(theta_cbw)/c;
        m = (-(K(i)-1):K(i)-1);
        exp_sum = sum(exp(-1j*u*m));
        g_denom = -2*cos(2*K(i)*pi*f(i)*delta*sin(theta_cbw)/c);
        g = exp_sum/g_denom;

        array_mid = (M+1)/2;
        w_mod_rect(array_mid-K(i)+1:array_mid+K(i)-1, i)=1;
        w_mod_rect(array_mid-K(i), i)=g;
        w_mod_rect(array_mid+K(i), i)=g;
        w_mod_rect(:,i) = w_mod_rect(:,i)/(2*g+2*K(i)-1);  % normalize
    end
end


plot_beampattern(w_mod_rect,f,M,plot_deg,plot_dB,plot_2d);



%% plot DPSS
% f = [3000 5000 8000]; % beampattern freqs
w_dpss = zeros([M length(f)]);
num_of_angles = 1001; % for theta axis
theta = linspace(0,pi/2,(num_of_angles+1)/2)'; % angles for B(f.theta)
theta_res = theta(2)-theta(1);
theta_tol = 2*theta_res; % allowed tolerance for the mainlobe width around theta_cbw

[m,n] = meshgrid(-N:N);
f_min = c/(delta*M*sin(theta_cbw));
f_max = c/delta;
for i=1:length(f)
    
    if (f(i) < f_min) | (f(i) > f_max)
        w_dpss(:,i)=1/M;
    else
        
        u0 = (2*pi*f(i)*delta/c)*sin(theta_cbw);
        % iterative search for the correct u0
        tmp_theta_cbw = theta_cbw;
        for k=1:100
            A = 2*sin((m-n)*u0)./(m-n);
            A(1:M+1:end) = 2*u0; % set diagonal to lim (m -> n) A[m,n];
            [w_dpss(:,i),D] = eigs(A,1);  % A is symmetric and real, hence always has real eigenvalues
                                              % this returns the eigenvector corr. to
                                              % the largest eigenvalue.
            w_dpss(:,i) = w_dpss(:,i)/sum(w_dpss(:,i));  % normalize

            % calc beampattern for 0<theta<pi/2 and find min
            u = 2*pi*f(i)*delta*sin(theta)/c;
            d = exp(-1j*u*m(1,:));
            B = abs(d*w_dpss(:,i));
            B(B<(10^-3)) = 10^-3; % set all low points to same value
            [min_B, min_B_idx] = min(B);
            if abs(theta(min_B_idx)-theta_cbw) < theta_tol
                break
            else
                if theta(min_B_idx) < theta_cbw
                    tmp_theta_cbw = tmp_theta_cbw + theta_res;
                else
                    tmp_theta_cbw = tmp_theta_cbw - theta_res;
                end
                u0 = (2*pi*f(i)*delta/c)*sin(tmp_theta_cbw);
            end
                
            
        end
    end

end


plot_beampattern(w_dpss,f,M,plot_deg,plot_dB,plot_2d);



