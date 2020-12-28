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

% model params
theta_cbw = deg2rad(20); % the angle of the first mainlobe null 
f = [4000, 5000, 6000]; % beampattern freqs
% f=4000;



%% plot delay and sum
 w_ds = (1/M)*ones([M,1]);
 w_ds = [w_ds w_ds w_ds]; % FIXME make this nicer
 plot_beampattern(w_ds,f,M,plot_deg,plot_dB);


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

for i=1:length(f) 
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


plot_beampattern(w_mod_rect,f,M,plot_deg,plot_dB);



%% plot DPSS

w_dpss = zeros([M length(f)]);

[m,n] = meshgrid(-N:N);
for i=1:length(f) 
    u0 = (2*pi*f(i)*delta/c)*sin(theta_cbw);
    A = 2*sin((m-n)*u0)./(m-n);
    A(1:M+1:end) = 1; % set diagonal to 1;
    A = A*(1/pi);
    [w_dpss(:,i),D] = eigs(A,1);  % A is symmetric and real, hence always has real eigenvalues
                                  % this returns the eigenvector corr. to
                                  % the largest eigenvalue.
    w_dpss(:,i) = w_dpss(:,i)/sum(w_dpss(:,i));  % normalize
end

plot_beampattern(w_dpss,f,M,plot_deg,plot_dB);


%%
% w - spatial filter weights vector
% f - vector of wanted freqs
% M - num of mics in array. choose odd value for M
function res = plot_beampattern (w, f, M, plot_deg, plot_dB)

    %consts 
    N=(M-1)/2;
    delta=0.035; % spatial sampling distance
    c=340; % speed of sound
    num_of_angles = 1001; % for theta axis

    %plot consts
    linewd = 0.8;
    hcfontsize = 20;

    m = (-N:N);
    theta = linspace(-pi/2,pi/2,num_of_angles)'; % angles of f
    
    B=zeros([length(theta), length(f)]);

    for i=1:length(f)

        u = 2*pi*f(i)*delta*sin(theta)/c;
        d = exp(-1j*u*m);
        B(:,i) = d*w(:,i);

    end
    
    B_dB = 20*log10(abs(B));

    figure
%     if plot_deg 
%         plot(rad2deg(theta), 10*log(abs(B(:,1))),'linewidth',linewd);
%         xlabel('theta [deg]');
%         xlim([-90 90]);
%     else
%         plot(theta, 10*log(abs(B(:,1))),'linewidth',linewd);
%         xlabel('theta [rad]');
%         xlim([-pi/2 pi/2]);
%     end
   
    if plot_dB
        plot(rad2deg(theta), B_dB(:,1),'linewidth',linewd);
        hold on;
        plot(rad2deg(theta), B_dB(:,2),'linewidth',linewd);
        plot(rad2deg(theta), B_dB(:,3),'linewidth',linewd);
    else
        plot(rad2deg(theta), abs(B(:,1)),'linewidth',linewd);
        hold on;
        plot(rad2deg(theta), abs(B(:,2)),'linewidth',linewd);
        plot(rad2deg(theta), abs(B(:,3)),'linewidth',linewd);
        hold off;
    end
    
    set(gca, 'Color', [1, 1, 1]); 
    set(gca, 'FontName', 'Times New Roman');
    set(gca, 'FontSize', hcfontsize);
    set(gca, 'LineWidth', linewd); 
    box on; grid on;
    xlabel('theta [deg]');
    ylabel('|B(f_0,\theta)| [dB]');
    xlim([-90 90]);
%     ylim([-50 0]);
    lgd = legend('f=4000', 'f=5000', 'f=6000');
    lgd.FontSize = 10;
    
    res=0;

end


