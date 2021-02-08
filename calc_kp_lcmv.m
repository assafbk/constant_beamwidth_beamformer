%
% w - the kp-lcmv
% w1, w2 - the weights for the virtual arrays
%
function [w, w1, w2] = calc_kp_lcmv(M1, delta, f, theta_cbw, theta_d, theta_nulls)

    % model params
    c=340; % speed of sound
    M2=1+length(theta_nulls); % num of mics in narrower virtual array
    M=M1+M2-1; % num of mics in the global array

    w1 = zeros([M1 length(f)]);
    w2 = zeros([M2 length(f)]);
    w_prod = zeros([M1*M2 length(f)]);

    if rem(M2,2) == 1  % init for virtual array 2
        N=(M2-1)/2;
        m = (-N:N)';
    else
        N=M2/2;
        m = (-(N-1):N)';
    end

    for i=1:length(f)
        [w1(:,i), beta] = calc_kaiser(theta_cbw, f(i), M1, delta);

        d_d = exp(-1j*(2*pi*f(i)*delta*sin(theta_d)/c)*m);
        C = [d_d];
        i_c = [1];
        for l=1:length(theta_nulls)
            C = [C exp(-1j*(2*pi*f(i)*delta*sin(theta_nulls(l))/c)*m)];
            i_c = [i_c; 0];
        end
        w2(:,i) = C*inv(C'*C)*i_c;

        w_prod(:,i) = kron(w2(:,i),w1(:,i));
    end

    w = produced_array_weights_to_equivalent_physical_array_weights(w_prod,M1,M2);

end

