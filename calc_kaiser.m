function w_kaiser = calc_kaiser(theta_cbw, f, M)

    %consts 
    N=(M-1)/2;
    delta=0.035; % spatial sampling distance
    c=340; % speed of sound

    m = -N:N;
    num_of_angles = 5001; % for theta axis
    theta = linspace(0,pi/2,(num_of_angles+1)/2)'; % angles for B(f.theta)
    theta_res = theta(2)-theta(1);
    theta_tol = 2*theta_res; % allowed tolerance for the mainlobe width around theta_cbw

    w_kaiser = zeros([M length(f)]);
    f_min = c/(delta*M*sin(theta_cbw));
    f_max = c/delta;
    for i=1:length(f)

        if (f(i) < f_min) | (f(i) > f_max)
            w_kaiser(:,i)=1/M;
        else

            A_sl = (26*M*f(i)*delta/c)*sin(theta_cbw)-12;
            beta = 0.76608*(A_sl-13.26)^0.4 + 0.09834*(A_sl-13.26);

            % iterative search for the correct beta
            tmp_theta_cbw = theta_cbw;
            for k=1:100
                w_kaiser(:,i) = besseli(0,beta*sqrt(1-(m/N).^2))/besseli(0,beta);
                w_kaiser(:,i) = w_kaiser(:,i)/sum(w_kaiser(:,i));  % normalize

                % calc beampattern for 0<theta<pi/2 and find min
                u = 2*pi*f(i)*delta*sin(theta)/c;
                d = exp(-1j*u*m(1,:));
                B = abs(d*w_kaiser(:,i));
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
                    A_sl = (26*M*f(i)*delta/c)*sin(tmp_theta_cbw)-12;
                    beta = 0.76608*(A_sl-13.26)^0.4 + 0.09834*(A_sl-13.26);
                end


            end
        end

end
