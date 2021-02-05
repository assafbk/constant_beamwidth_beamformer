function w_cdpss = calc_constrained_dpss(constraint_thetas, theta_cbw, attn_at_constraint_dB, f, M)
    
    addpath(genpath('./CRQPACK'))

    crq_opts.maxit=300;
    crq_opts.tol=1e-10;
    crq_opts.method=1;
%     crq_opts.resopt=1;

    %consts 
    N=(M-1)/2;
    delta=0.035; % spatial sampling distance
    c=340; % speed of sound
    
    w_cdpss = zeros([M length(f)]);
    num_of_angles = 2500; % for theta axis
    theta = linspace(deg2rad(10),pi/2,(num_of_angles+1)/2)'; % angles for B(f.theta)
    theta_res = theta(2)-theta(1);
    theta_tol = 2*theta_res; % allowed tolerance for the mainlobe width around theta_cbw
    attn_at_constraint = sqrt(db2pow(attn_at_constraint_dB)); % the constraint's units are amplitude (they are w'*d, not |w'*d|^2)

    [m,n] = meshgrid(-N:N);
%     f_min = c/(delta*M*sin(theta_cbw));
    f_min = 3000;
    f_max = c/delta;
    for i=1:length(f)
        
        % calc constraints
        theta_d = constraint_thetas(1);
        theta_1 = constraint_thetas(2);
        theta_2 = constraint_thetas(3);
        theta_3 = constraint_thetas(4);
        theta_4 = constraint_thetas(5);

        m_tag = (-N:N).';
        d_d = exp(-1j*(2*pi*f(i)*delta*sin(theta_d)/c)*m_tag);
        d_1 = exp(-1j*(2*pi*f(i)*delta*sin(theta_1)/c)*m_tag);
        d_2 = exp(-1j*(2*pi*f(i)*delta*sin(theta_2)/c)*m_tag);
        d_3 = exp(-1j*(2*pi*f(i)*delta*sin(theta_3)/c)*m_tag);
        d_4 = exp(-1j*(2*pi*f(i)*delta*sin(theta_4)/c)*m_tag);
        i_c = [1 attn_at_constraint attn_at_constraint attn_at_constraint attn_at_constraint]';
        C = [d_d d_1 d_2 d_3 d_4];
        
%         theta_2 = constraint_thetas(3);
%         d_2 = exp(-1j*(2*pi*f(i)*delta*sin(theta_2)/c)*m_tag);
%         i_c = [1 attn_at_constraint attn_at_constraint]';
%         C = [d_d d_1 d_2];


        if (f(i) < f_min) | (f(i) > f_max)
            w_cdpss(:,i)=C*inv(C'*C)*i_c;  % if freq is too low, take MN solution
        else

            u0 = (2*pi*f(i)*delta/c)*sin(theta_cbw);
            % iterative search for the correct u0
            tmp_theta_cbw = theta_cbw;
            for k=1:100
                A = 2*sin((m-n)*u0)./(m-n);
                A(1:M+1:end) = 2*u0; % set diagonal to lim (m -> n) A[m,n];
%                 [w_cdpss(:,i),info] = CRQ_Lanczos(-A, C, i_c, crq_opts);
                [w_cdpss(:,i),info] = CRQ_Lanczos(-A, C, i_c, crq_opts);
%                 [w_cdpss(:,i),info] = CRQ_explicit(-A, C, i_c);
%                 w_cdpss(:,i) = w_cdpss(:,i)/sum(w_cdpss(:,i));  % normalize

                % calc beampattern for 0<theta<pi/2 and find min
                u = 2*pi*f(i)*delta*sin(theta)/c;
                d = exp(-1j*u*m(1,:));
                B = abs(d*w_cdpss(:,i));
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
end




