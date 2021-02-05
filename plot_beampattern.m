%%
% w - spatial filter weights vector
% f - vector of wanted freqs
% M - num of mics in array. choose odd value for M
function res = plot_beampattern (w, f, M, plot_deg, plot_dB, plot_2d)

    %consts 
    N=(M-1)/2;
    delta=0.035; % spatial sampling distance
    c=340; % speed of sound
    num_of_angles = 10001; % for theta axis

    %plot consts
    linewd = 1;
    hcfontsize = 20;

    m = (-N:N);
    theta = linspace(-pi/2,pi/2,num_of_angles)'; % angles of f
    theta(7223)=deg2rad(40);
    theta(1668)=deg2rad(-60);
    theta(8890)=deg2rad(70);
    
    B=zeros([length(theta), length(f)]);

    for i=1:length(f)

        u = 2*pi*f(i)*delta*sin(theta)/c;
        d = exp(-1j*u*m);
        B(:,i) = d*w(:,i);

    end
    
    B_dB = 20*log10(abs(B));

    figure
    if plot_2d
        if plot_deg 
            theta = rad2deg(theta);
        end
            [THETA,F] = meshgrid(theta,f);  
            refline(-1,1)
            hold on
            surf(THETA', F', B_dB, 'edgecolor', 'none');
            caxis([-60 0])
            xlabel('\theta [deg]');
            ylabel('f [Hz]');
%             title("Delay&Sum beampattern")
        if plot_deg 
            xlim([-90 90])
        end
            colorbar()
        
    else
        
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
            plots = B_dB;  % plots power ratio in dB
        else
            plots = abs(B); % plots |amplitude|
        end
        plot(rad2deg(theta), plots(:,1),'linewidth',linewd);
        if length(f)>1
            hold on;
          for i=2:length(f)
              plot(rad2deg(theta), plots(:,i),'linewidth',linewd);
          end
          hold off;
        end
        

        set(gca, 'Color', [1, 1, 1]); 
        set(gca, 'FontName', 'Times New Roman');
        set(gca, 'FontSize', hcfontsize);
        set(gca, 'LineWidth', linewd); 
        box on; grid on;
        
        if plot_deg
            xlim([-90 90]);
            xlabel('\theta [deg]');
        else
            xlim([-pi/2 pi/2]);
            xlabel('\theta [rad]');
        end
        
        if plot_dB
%             ylim([-80 10]);
            ylim([-80 inf]);
            ylabel('|B(f_0,\theta)|^2 [dB]');
        else
            ylabel('|B(f_0,\theta)|');
        end
%         lgd = legend('f=1000', 'f=2000', 'f=3000', 'f=4000', 'f=5000', 'f=6000');
%         lgd = legend('f=3000', 'f=4500', 'f=6000');
%         lgd = legend('f=4500', 'f=6000', 'f=7500');
%           lgd = legend('f=4500');
%         lgd = legend('upper filter', 'gsc', 'adaptive (lower) filter');
%         lgd = legend('fbf', 'gsc', 'ad');
        lgd = legend('fbf', 'gsc');
%         lgd = legend('f=6000');
        lgd.FontSize = 28;

        res=0;
    end
end


