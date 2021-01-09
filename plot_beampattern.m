%%
% w - spatial filter weights vector
% f - vector of wanted freqs
% M - num of mics in array. choose odd value for M
function res = plot_beampattern (w, f, M, plot_deg, plot_dB, plot_2d)

    %consts 
    N=(M-1)/2;
    delta=0.035; % spatial sampling distance
    c=340; % speed of sound
    num_of_angles = 1001; % for theta axis

    %plot consts
    linewd = 1;
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
            plot(rad2deg(theta), B_dB(:,1),'linewidth',linewd);
            hold on;
            plot(rad2deg(theta), B_dB(:,2),'linewidth',linewd);
            plot(rad2deg(theta), B_dB(:,3),'linewidth',linewd);
%             plot(rad2deg(theta), B_dB(:,4),'linewidth',linewd);
%             plot(rad2deg(theta), B_dB(:,5),'linewidth',linewd);
%             plot(rad2deg(theta), B_dB(:,6),'linewidth',linewd);
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
        xlabel('\theta [deg]');
        ylabel('|B(f_0,\theta)| [dB]');
        xlim([-90 90]);
        ylim([-80 0]);
%         lgd = legend('f=1000', 'f=2000', 'f=3000', 'f=4000', 'f=5000', 'f=6000');
        lgd = legend('f=3000', 'f=4500', 'f=6000');
%         lgd = legend('f=6000');
        lgd.FontSize = 28;

        res=0;
    end
end


