function w_lcmv = calc_lcmv(M, delta, f, theta_d, theta_nulls, noise_field)
    
    % model params
    c = 340; % speed of sound [m/s]
    w_lcmv = zeros([M length(f)]);

    if rem(M,2) == 1  % init for lcmv
        N=(M-1)/2;
        m = (-N:N)';
    else
        N=M/2;
        m = (-(N-1):N)';
    end

    for i=1:length(f)
        d_d = exp(-1j*(2*pi*f(i)*delta*sin(theta_d)/c)*m);
        C = [d_d];
        i_c = [1];
        for l=1:length(theta_nulls)
            C = [C exp(-1j*(2*pi*f(i)*delta*sin(theta_nulls(l))/c)*m)];
            i_c = [i_c; 0];
        end

    % set the noise corr matrix (Phi_v)
    if strcmp(noise_field,"white")
        sigma_w = 1;
        Phi_v = sigma_w*eye(M);

    elseif strcmp(noise_field,"diffuse")
        [I,J] = meshgrid(1:M,1:M);
        Phi_v = sinc(2*pi*f(i)*(J-I)*delta/c);

    elseif strcmp(noise_field,"power") % some other noise field   
        alpha=0.8;
        [I,J] = meshgrid(1:M,1:M);
        Phi_v = alpha^abs(J-I);
    else
        disp("unsupported noise field");
    end

        w_lcmv(:,i) = inv(Phi_v)*C*inv(C'*inv(Phi_v)*C)*i_c;
    end

end

