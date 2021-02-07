function snr_gain_dB = calc_snr_gain(w, f, delta, theta_d, noise_field, is_plot)

%plot consts
linewd = 1;
hcfontsize = 20;

% model params
c=340; % speed of sound
M = length(w(:,1)); % num of mics
[I,J] = meshgrid(1:M,1:M); % for diffuse noise field

if rem(M,2) == 1  % init for desired angle's steering vector
    N=(M-1)/2;
    m = (-N:N)';
else
    N=M/2;
    m = (-(N-1):N)';
end

snr_gain = zeros(size(f));

    for i=1:length(f)
        
        d_d = exp(-1j*(2*pi*f(i)*delta*sin(theta_d)/c)*m);
        
        if strcmp("white", noise_field)
            Gamma = eye(M);
        elseif strcmp("diffuse", noise_field)
            Gamma = sinc(2*pi*f(i)*(J-I)*delta/c);
        else
            disp("unsupported noise field");
            break;
        end
        
        snr_gain(i) = (abs(w(:,i)'*d_d)^2)/(w(:,i)'*Gamma*w(:,i));
    end
snr_gain_dB = 10*log10(abs(snr_gain)); 

if is_plot == false
    return

figure
plot(f, snr_gain_dB,'linewidth',linewd);

set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;


xlabel('f [Hz]');
ylabel('SNR Gain [dB]');
% xlim([-90 90]);
% ylim([-80 inf]);
% lgd = legend('f=3000', 'f=4500', 'f=6000');
% lgd.FontSize = 28;


end