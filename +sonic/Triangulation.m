% This software is made available under the MIT license. See SONIC/LICENSE
% for details.
classdef Triangulation

    methods (Static)

        function [rI, cov_r] = triangulate(x_C, p_I, attitudes, ...
                method, cov_x, compute_cov)
            %% [rI, cov_r] = triangulate(x_C, p_I, attitudes, method, cov_x, compute_cov)
            %   Triangulate using the specified method.
            %
            %   References:
            %   - [1] S. Henry and J. A. Christian. Absolute Triangulation
            %   Algorithms for Space Exploration. JGCD (2023).
            %   https://doi.org/10.2514/1.G006989
            %
            %   Inputs:
            %       - x_C (sonic.Points2): image_plane coordinates 
            %         (equivalent of K^-1 * pixel coordinate)
            %       - p_I  (sonic.Points3): respective 3D position objects
            %         that are used to triangulate (resection) OR position of
            %         the respective cameras (intersection)
            %       - attitudes
            %           - (nx1 sonic.Attitude): attitudes of the cameras
            %   
            %           OR
            %
            %           - (1x1 sonic.Attitude): then all cameras have 
            %             the same attitude
            %       - method (string): "midpoint", "DLT", or "LOST"
            %       - cov_x (OPTIONAL)
            %           - DEFAULT: unit variance
            %           - (1x1 double): isotropic variance, same for all
            %             measurements
            %           - (1xn double): isotropic variance of points2
            %           - (2x2xn double): NOT YET IMPLEMENTED. Isotropic
            %             variance MUST be assumed.
            %       - compute_cov (1x1 logical) (OPTIONAL): flag to decide 
            %         whether to compute covariance or not. This will make 
            %         the triangulation function slightly slower, 
            %         especially for DLT and midpoint. Defaults to false.
            %           
            %   Outputs:
            %       - rI (sonic.Points3): triangulated 3D point in inertial
            %       - cov_r (3x3 double): 3D covariance of the triangulated
            %         point
            %
            %   Last revised: 4/26/24
            %   Last author: Sebastien Henry

            arguments
                x_C         (1, 1)  sonic.Points2
                p_I         (1, 1)  sonic.Points3
                attitudes   (1, :)  sonic.Attitude
                method      (1, 1)  string
                cov_x       (1, :)  double = 1
                compute_cov (1, 1)  logical = false
            end
            
            % obtain the number of measurements
            n_meas = x_C.n;

            % obtain number of attitudes given
            length_attitudes = length(attitudes);
            
            % some input verifications
            if n_meas < 2
                error('sonic:Triangulation:inputError', ...
                    'the number of measurements should be at least 2');
            end
            
            if length_attitudes > 1 && length_attitudes ~= n_meas
                disp(n_meas)
                disp(length_attitudes)
                error('sonic:Triangulation:inputError', ...
                    ['Attitudes should either be a scalar, or a vector' ...
                    ' of size n_meas']);
            end
            
            if p_I.n ~= n_meas
                error('sonic:Triangulation:inputError', ...
                    ['the number of 3D points (p_I) should be the same' ...
                    ' as the number of measurements (x_C)']);
            end

            if length(cov_x) == 1
                cov_x = repmat(cov_x, 1, n_meas);
            end

            % preallocation to speed-up code
            H = zeros(2*n_meas, 3);
            y = zeros(2*n_meas, 1);
            if compute_cov
                Prr_inv = zeros(3);
            end

            % S is a matrix such that S*A removes the third row of A
            S = [eye(2), [0; 0]];

            % build the system in Eq. 20 from Ref. [1]
            for ii = 1:n_meas
                % fetch data related to the ii-th measurement
                xi = x_C.p2(:, ii);

                if length_attitudes==1
                    T_I_to_Ci = attitudes.dcm;
                else
                    T_I_to_Ci = attitudes(ii).dcm;
                end
                pIi = p_I.r3(:, ii);

                switch method
                    case 'DLT'
                        % DLT sets the weights to 1
                        qi = 1;

                    case 'midpoint'
                        % midpoint is equivalent to DLT, but with unit
                        % normalized measurements
                        qi = 1;
                        xi = xi/norm(xi);

                    case 'LOST'
                        % selection of another measurement jj to estimate
                        % the gamma coefficient (see below). The following
                        % line arbitrarily selects the measurement that
                        % is n_meas/2 indexes away and cycles when
                        % jj > n_meas
                        jj = max([mod( round(ii+n_meas/2), n_meas), 1]);

                        % adress the case jj == ii,
                        % since distinct measurements are needed
                        if ii == jj
                            jj = mod(ii+1,n_meas);
                        end
                        if jj == 0
                            if ii==1
                                jj = 2;
                            else
                                jj = 1;
                            end
                        end

                        % fetch data related to the jj-th measurement
                        xj =  x_C.p2(:, jj);
                        if length_attitudes==1
                            T_I_to_Cj = attitudes.dcm;
                        else
                            T_I_to_Cj = attitudes(jj).dcm;
                        end
                        pIj = p_I.r3(:, jj);

                        sig2 = cov_x(ii);

                        % take the difference of positions associated with
                        % the two measurements as in Eq. 46 from Ref. [1]
                        dIij = pIj-pIi;

                        % compute the square of the gamma coefficient in
                        % Eq. 108 from Ref. [1]
                        gamma2 = (norm(cross(dIij,T_I_to_Cj'*xj)) ...
                            /norm(cross(T_I_to_Ci'*xi,T_I_to_Cj'*xj)))^2;
                        qi = sqrt(1/sig2/gamma2);

                    otherwise
                        error('sonic:Triangulation:invalidMethod', ...
                            ['Invalid triangulation method specified. ' ...
                            'Specify one of `DLT`, `midpoint`, or `LOST`.']);
                end

                % precomputation of a matrix to speedup code.
                % Note that the third row is linearly dependant to the two
                % others. Therefore it is removed by pre-multiplication
                % by the S matrix.
                C = qi * S * sonic.Math.crossmat(xi) * T_I_to_Ci;

                % fill the matrix on the left side of Eq. 20 from Ref. [1]
                H((2*ii-1):2*ii, :) = C;

                % fill the matrix on the right side of Eq. 20 from Ref. [1]
                y((2*ii-1):2*ii) = C * pIi;
                
                if compute_cov
                    % compute the sum inside the parenthesis
                    % of Eq. 117 from Ref. [1]
                    Prr_inv = Prr_inv + C' * C;
                end

            end
            
            % solve the system to find the position
            [~, ~, V ] = svd([H y]);
            rI = sonic.Points3(-V(1:3,end)/V(end,end));

            if compute_cov
                switch method
                    case 'LOST'
                        % invert the sum to obtain the position covariance
                        % of Eq. 117 from Ref. [1]
                        cov_r = inv(Prr_inv);

                    otherwise
                        % compute covariance for DLT or midpoint
                        % from the estimated position

                        sum_eps = zeros(3);
                        for ii = 1:n_meas
                            % fetch data related to the ii-th measurement
                            xi = x_C.p2(:, ii);
                            pIi = p_I.r3(:, ii);
                            sig2 = cov_x(ii);

                            % fetch attitude
                            if length_attitudes==1
                                T_I_to_Ci = attitudes.dcm;
                            else
                                T_I_to_Ci = attitudes(ii).dcm;
                            end  

                            switch method
                                case 'midpoint'
                                    % normalize covariance of measurement
                                    sig2 = sig2/(norm(xi)^2);
                            end
                            
                            % covariance of residual
                            R_eps = - sig2 * S * sonic.Math.crossmat(T_I_to_Ci*(pIi-rI.r3)) * (S'*S) * sonic.Math.crossmat(T_I_to_Ci*(pIi-rI.r3)) * S';
                            
                            % H^T * R_eps * H
                            sum_eps = sum_eps + ... 
                                H((2*ii-1):2*ii,:)' * R_eps * H((2*ii-1):2*ii,:);
                        end
                        
                        % the well known least squares covariance
                        inv_HTH = (H'*H)^-1;
                        cov_r = inv_HTH * sum_eps * inv_HTH;
                end
            end

            
        end

    end
end