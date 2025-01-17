classdef PositionEstimation
    % This software is made available under the MIT license. 
    % See SONIC/LICENSE for details.
    
    methods (Static)

        function [rI, cov_r] = triangulate(x_C, p_I, attitudes, method, cov_x, computeCov)
        %% [rI, cov_r] = triangulate(x_C, p_I, attitudes, method, cov_x, computeCov)
        %  Triangulate using the specified method.
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
        %       - computeCov (1x1 logical) (OPTIONAL): flag to decide 
        %         whether to compute covariance or not. This will make 
        %         the triangulation function slightly slower, 
        %         especially for DLT and midpoint. Defaults to false.
        %           
        %   Outputs:
        %       - rI (sonic.Points3): triangulated 3D point in inertial
        %       - cov_r (3x3 double): 3D covariance of the triangulated
        %         point
        %
        %   References:
        %       - [1] S. Henry and J. A. Christian. Absolute Triangulation
        %       Algorithms for Space Exploration. JGCD (2023).
        %       https://doi.org/10.2514/1.G006989
        %
        %   Last revised: 4/26/24
        %   Last author: Sebastien Henry

            arguments
                x_C         (1, 1)  sonic.Points2
                p_I         (1, 1)  sonic.Points3
                attitudes   (1, :)  sonic.Attitude
                method      (1, 1)  string
                cov_x       (1, :)  double = 1
                computeCov  (1, 1)  logical = false
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

            if isscalar(cov_x)
                cov_x = repmat(cov_x, 1, n_meas);
            end

            % preallocation to speed-up code
            H = zeros(2*n_meas, 3);
            y = zeros(2*n_meas, 1);
            if computeCov
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
                
                if computeCov
                    % compute the sum inside the parenthesis
                    % of Eq. 117 from Ref. [1]
                    Prr_inv = Prr_inv + C' * C;
                end

            end
            
            % solve the system to find the position
            [~, ~, V ] = svd([H y]);
            rI = sonic.Points3(-V(1:3,end)/V(end,end));

            if computeCov
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

        function posI = withConics(conics, matchedConics_E, attItoC, conicPose, method)
        %% posI = withConics(conics, matchedConics_E, attItoC, conicPose, method)
        % Perform least-squares positioning using image-plane conics, 
        % matched with the observed conics in their local frames, E 
        % Section 9.1 of "Lunar Crater ID..." paper
        % This method assumes that conic association has already been
        % performed and the image-plane conics from the image are each
        % associated with the reference conics defined in their local 
        % frames, matchedConics_E
        %
        %   Inputs:
        %   - conics (Nx1 sonic.Conic): array of image-plane conics
        %   - attItoC (1x1 sonic.Attitude): attitude from the reference 
        %       frame of interest I, wrt which we want to position   
        %       ourselves, to the camera frame C. The reference frame can  
        %       be any frame (inertial, body-fixed). 
        %   - matchedConics_E (Nx1 sonic.Conic): array of matched conics in
        %       their local frames, E, (e.g., could be a local ENU frame)
        %   - attConicsEtoI_arr (Nx1 sonic.Attitude): array of attiude
        %       objects which represent the attitude transformation of each 
        %       conic from their respective local frames to the reference 
        %       frame of interest I
        %   - conicCenters_I (1x1 sonic.Points3): center points of each 
        %       conic in the reference frame of interest I
        %   - method (1xn string): String indicating method of position
        %       estimation using conics. Supported methods below:
        %       - "leastsquares" 
        %
        %   Outputs:
        %   - posI (1x1 sonic.Points3): solved least-squares position 
        %       solution for least-squares conic positioning
        %
        %   Last revised: Nov 21, 2024
        %   Last author: Tara Mina
            arguments
                conics              (:,1) sonic.Conic
                matchedConics_E     (:,1) sonic.Conic
                attItoC             (1,1) sonic.Attitude             
                % attConicsEtoI_arr   (:,1) sonic.Attitude
                % conicCenters_I      (1,1) sonic.Points3
                conicPose           (:,1) sonic.Pose
                method              (1,:) string   
            end
            
            % Extract attitudes and positions from Pose object
            attConicsEtoI_arr = [conicPose.att]';
            conicCenters_I = [conicPose.t]';

            nconics = size(conics,1);
            assert(nconics == size(matchedConics_E,1));
            assert(nconics == size(attConicsEtoI_arr,1));
            assert(nconics == size(conicCenters_I,1));
        
            % Extract camera attitude transformation matrix
            TMtoC = attItoC.dcm;
        
            % Check specified method to solve
            switch method
                case "leastsquares"
                    % for each conic, create necessary matrices to setup LS problem
                    S = [eye(2); zeros(1,2)];
                    Als = nan*ones(2*nconics,3);
                    bls = nan*ones(2*nconics,1);
                    kvec = [0;0;1];
                    for i=1:nconics
                        % extract image crater locus for conic i
                        imgConic_i = conics(i,1);
                        Ai = imgConic_i.locus;
                    
                        % extract true crater's conic locus for conic i
                        refConicE_i = matchedConics_E(i,1);
                        Ci = refConicE_i.locus;
                    
                        % define Bi matrix for conic i 
                        Bi = TMtoC' * Ai * TMtoC;
                    
                        % compute the unknown scale for conic i
                        matStCS = S' *  Ci * S;
                        vecStCS = matStCS(:);
                        TEtoI_i = attConicsEtoI_arr(i,1).dcm;
                        matBi_toE = S' * TEtoI_i' * Bi * TEtoI_i * S;
                        si = vecStCS'*matBi_toE(:) / (vecStCS' * vecStCS);
                    
                        % build the system of equations for least squares
                        rowidcs = ((2*i)-1):(2*i);
                        matSTB = S' * TEtoI_i' * Bi;
                        posM_i = conicCenters_I(i).r3;
                        Als(rowidcs,:) = matSTB;
                        bls(rowidcs,:) = matSTB*posM_i - si*(S' * Ci * kvec);
                    end
                    
                    % Now we've built our system of equations, solve in LS sense
                    posI_arr = Als \ bls; 
                    posI = sonic.Points3(posI_arr);

                otherwise
                        error('sonic:PositionEstimation:withConics', ...
                        'Input method type to solve is not supported. Please check the list of supported method types.');
            end
        end

        function [r_C, cov_r] = withHorizonPts(limbPts, radii_P, att_P2C, cam, cov_x, computeCov)
        %% [r_C, cov_r] = withHorizonPts(limbPts, radii_P, att_P2C, cam, cov_x, computeCov)
        %   Performs horizon based optical navigation using the Christian-Robinson
        %   algorithm as detailed in "A Tutorial on Horizon-Based Optical 
        %   Navigation and Attitude Determination With Space Imaging Systems"
        %
        %   Reference DOI: 10.1109/ACCESS.2021.3051914
        %
        %   Inputs:
        %       - limbPts (1x1 sonic.Points2): Image plane coordinates of 
        %         the limb
        %       - radii_P (3x1 double): Celestial body's principle ai radii 
        %         as described in the planet's principle axis frame
        %       - att_P2C (1x1 sonic.Attitude): Attitude transformation from the 
        %         celestial body's principle axis frame to the camera frame
        %       - cam (1x1 sonic.)
        %       - cov_x (OPTIONAL)
        %           - DEFAULT: unit variance
        %           - (1x1 double): isotropic variance, same for all
        %             measurements
        %       - computeCov (1x1 logical) (OPTIONAL): flag to decide 
        %         whether to compute covariance or not. This will make 
        %         the triangulation function slightly slower, 
        %         especially for DLT and midpoint. Defaults to false.
        % 
        %   Outputs:
        %       - r_C (1x1 sonic.Points3): vector from the camera to the 
        %           center of the celestial body expressed in the camera 
        %           frame
        %       - cov_r (3x3 double): 3D covariance of the estimated
        %           celestial body position
        %
        %   Last revised: 3/27/24
        %   Last author: Ava Thrasher
        
            arguments
                limbPts     (1,1) sonic.Points2
                radii_P     (3,1) double
                att_P2C     (1,1) sonic.Attitude
                cam         (1,1) sonic.Camera
                cov_x       (1,1)  double = 1
                computeCov  (1,1) logical = false
            end

            % Get rotation matrix
            T_P2C = att_P2C.dcm;
            T_C2P = T_P2C';
        
            % Compute D eq 98
            Dinv = diag(radii_P);
            D = diag(1./radii_P);
        
            % Compute R = DT eq 102
            R = D*T_C2P;
                    
            % Get homogeneous pixel coordinates (p2) 
            % %%
            xbar = limbPts.p2;
        
            % Initialize storage for H
            H = zeros(length(xbar),3);
            % for i 1 to n
            for i = 1:length(xbar)
                % xbar_i' = R ubari eq 95
                xbari = R*xbar(:,i);
        
                % s_i' = xbar_i'/norm(xbar_i') eq 105
                si = xbari./norm(xbari);
                H(i,:) = si'; % eq 109
            end
            
            % Compute LS solution for n eq 108
            n = H\ones(length(H),1);
        
            % Compute rc eq 119
            r_C = (1/sqrt(n'*n - 1))*T_P2C*Dinv*n;
        
            % Output to Points3 object
            r_C = sonic.Points3(r_C);

            if computeCov
                cov_r = sonic.PositionEstimation.CRACov(limbPts, r_C, radii_P, cam, att_P2C, cov_x);
            end
        end

        function P = parametricCovCRA(rEst_KM,Rp_KM,camObj,esun,n,cov_x,thetaMax_RAD)
        %% P = parametricCovCRA(rEst_KM,Rp_KM,camObj,esun,n,cov_x,thetaMax_RAD)
        %   Computes the parametric covariance for a position estimate from
        %   horizon based position estimation. ONLY VALID FOR SPHERICAL
        %   OBSERVED BODY ASSUMPTION.
        %   Inputs:
        %       - rEst_KM (1x1 sonic.Points3): vector from spacecraft to 
        %           target in the camera frame
        %       - Rp_KM (1x1 double): radius of spherical planet
        %       - camObj (1x1 sonic.Camera): camera object
        %       - esun (1x1 sonic.Points3): vector from target to 
        %           sun expressed in camera frame
        %       - n (1x1 double): number of observed horizon points
        %       - cov_x (1x1 double): variance of geometric 
        %       distance between observed edge points and best fit ellipse
        %       - thetaMax_RAD (1x1 double): limb observation half angle
        %
        %   Outputs:
        %       - P (3x3 double): covariance matrix expressed in the camera 
        %           frame
        
            arguments
                rEst_KM         (1,1) sonic.Points3
                Rp_KM           (1,1) double
                camObj          (1,1) sonic.Camera
                esun            (1,1) sonic.Points3
                n               (1,1) double
                cov_x           (1,1) double
                thetaMax_RAD    (1,1) double
            end

            % Distance from spacecraft to planet
            rho = norm(rEst_KM.r3);
            
            % Define cone principal axis frame
            % See Eqs. 63-66 in [Hikes et al., 2017]
            ez = rEst_KM.r3/rho;
            esun = esun.r3/norm(esun.r3);
            ey = cross(ez, esun);
            ey=ey/norm(ey);
            ex = cross(ey,ez); 
            ex = ex/norm(ex);
            ang = 0*pi/180;
            R = [cos(ang) sin(ang) 0; -sin(ang) cos(ang) 0; 0 0 1];
            T_coneP2cam = R*[ex ey ez]; % T_coneP2cam = [ex ey ez];
            T_cam2coneP = T_coneP2cam';
            
            % Compute scalar parameter D
            % See Eq. 84 in [Hikes et al., 2017]
            D = (thetaMax_RAD/4)*(2*thetaMax_RAD + sin(2*thetaMax_RAD)) - sin(thetaMax_RAD)^2;
            
            % Compute covariance matrix
            % See Eq. 85 in [Hikes et al., 2017]
            dx = camObj.d_x;
            P11 = thetaMax_RAD/D;
            P13 = (sqrt(rho^2-Rp_KM^2)/(D*Rp_KM))*sin(thetaMax_RAD);
            P22 = 4/(2*thetaMax_RAD-sin(2*thetaMax_RAD));
            P33 = ((rho^2-Rp_KM^2)*(2*thetaMax_RAD + sin(2*thetaMax_RAD)))/(4*D*Rp_KM^2);
            P = (cov_x*rho^4*thetaMax_RAD)/(n*dx^2*(rho^2-Rp_KM^2)) * T_coneP2cam * ...
                [P11    0       P13;...
                0       P22     0;...
                P13     0       P33] * T_cam2coneP;
        
        end

        function P = CRACov(limbPts, rc_KM, radii_P, cam, att_P2C, cov_x)
        %% P = CRACov(limbPts, rc_KM, radii_P, cam, att_P2C, cov_x)
        % Calculate the covariance of the maximum likelihood estimate of the camera
        % position from horizon based position estimation
        %
        %   Inputs:
        %       - limbPts (1x1 sonic.Points2): Image plane coordinates of 
        %           the limb
        %       - rc_KM (1x1 sonic.Points3): The estimated position of the
        %           celestial body in the camera frame in kilometers
        %       - radii_P (3x1 double): Celestial body's principle axis radii 
        %           as described in the planet's principle axis frame
        %       - cam (1x1 sonic.Camera): Camera object
        %       - att_P2C (1x1 sonic.Attitude): Attitude from the planetary
        %           axis frame to the camera frame
        %       - cov_x (1x1 double): variance in horizon pixel
        %       localization
        %
        %   Outputs:
        %       - P (3x3 double): Position estimation covariance associated with rc
        %
        %   Last revised: 1/8/25
        %   Last author: Ava Thrasher
        arguments
            limbPts     (1,1) sonic.Points2
            rc_KM       (1,1) sonic.Points3
            radii_P     (3,1) double
            cam         (1,1) sonic.Camera
            att_P2C     (1,1) sonic.Attitude
            cov_x       (1,1) double
        end
        
            % Set up estimate and shape matrix
            rc = rc_KM.r3;
            Ap = diag(radii_P);
            att_C2P = sonic.Attitude(att_P2C.dcm');
            Ac = att_C2P.dcm'*Ap*att_C2P.dcm; 
        
            % Calculate Mc eq. 9.116
            Mc = Ac*(rc*rc')*Ac - (rc'*Ac*rc - 1)*Ac;
            
            % Set the measurement covariance ex. 8.6
            dx = cam.d_x;

            Rxi = cov_x/(dx)^2*diag([1;1;0]);
            
            % Initialize inverted covariance
            PiInv = zeros(3);
            
            % Iterate through each point to construct the covariance
            for i = 1:limbPts.n
            
                xi = limbPts.p2(:,i);
        
                % Calculate variance eq. 9.145
                vari = 4*xi'*Mc*Rxi*Mc*xi;
                
                % Calculate B eq. 9.147
                Hi = 2*rc'*((xi'*Ac*xi)*Ac - Ac*(xi*xi')*Ac);
                
                % Update Pi
                PiInv = PiInv + (Hi'*Hi)/vari; 
            
            end
            
            % Invert to get the position estimate covariance
            P = inv(PiInv);
        end

    end

end