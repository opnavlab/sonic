classdef PnP
    % A class to solve for the PnP

    methods (Static)

        function [pose, covVecP] = PnPCalibrated(x_C, points_I, method, covx, computeCov)
        %% [pose, covVecP] = PnPCalibrated(x_C, points_I, method, covx, computeCov)
        % Estimate the pose of a calibrated camera using the
        % Perspective-n-Point with specified method.
        %
        %   References:
        %   - [1] Henry, S., & Christian, J. A. (2024).
        %      Optimal DLT-based Solutions for the Perspective-n-Point.
        %      https://doi.org/10.48550/arXiv.2410.14164
        %
        %   Inputs:
        %       - x_C (sonic.Points2): image plane coordinates
        %       (equivalent of K^-1 * pixel coordinates)
        %       - points_I  (sonic.Points3): respective 3D position objects
        %         that are used for estimation
        %       - method (string): "DLT" (normalized DLT)
        %                       or "oDLT" (optimal DLT)
        %       - covx (OPTIONAL)
        %           - DEFAULT: unit variance
        %           - (1x1 double): isotropic variance, same for all
        %             measurements
        %           - (1xn double): isotropic variance of points2
        %       - computeCov (1x1 logical) (OPTIONAL): flag to decide
        %         whether to compute covariance or not. This will make
        %         the PnP function slightly slower. Defaults to false.
        %
        %   Outputs:
        %       - pose (sonic.Pose): camera pose
        %       - covVecP (12x12 double) covariance of the vectorized 
        %       camera matrix P
        %
        %   Last revised: 11/24/24
        %   Last author: Sebastien Henry

            arguments
                x_C         (1, 1)  sonic.Points2
                points_I    (1, 1)  sonic.Points3
                method      (1, 1)  string
                covx        (1, :)  double = 1
                computeCov  (1, 1)  logical = false
            end

            [pose, ~, covVecP] = sonic.PnP.PnPDLT(x_C, points_I, method, true, covx, computeCov);

        end

        function [pose, camera, covVecP] = PnPUncalibrated(u_C, p_I, method, covu, computeCov)
        %% [pose, camera, covVecP] = PnPUncalibrated(u_C, p_I, method, covu, computeCov)
        % Estimate the pose of a calibrated camera using the
        % Perspective-n-Point with specified method.
        %
        %   References:
        %   - [1] Henry, S., & Christian, J. A. (2024).
        %      Optimal DLT-based Solutions for the Perspective-n-Point.
        %      https://doi.org/10.48550/arXiv.2410.14164
        %
        %   Inputs:
        %       - u_C (sonic.Points2): pixel coordinates
        %       - points_I  (sonic.Points3): respective 3D position objects
        %         that are used for estimation
        %       - method (string): "DLT" (normalized DLT)
        %                       or "oDLT" (optimal DLT)
        %       - covu (OPTIONAL)
        %           - DEFAULT: unit variance
        %           - (1x1 double): isotropic variance, same for all
        %             measurements
        %           - (1xn double): isotropic variance of points2
        %       - computeCov (1x1 logical) (OPTIONAL): flag to decide
        %         whether to compute covariance or not. This will make
        %         the PnP function slightly slower. Defaults to false.
        %
        %   Outputs:
        %       - pose (sonic.Pose): camera pose
        %       - camera (sonic.Camera): calibrated camera
        %       - covVecP (12x12 double) covariance of the vectorized 
        %       camera matrix P
        %
        %   Last revised: 11/24/24
        %   Last author: Sebastien Henry

            arguments
                u_C         (1, 1)  sonic.Points2
                p_I         (1, 1)  sonic.Points3
                method      (1, 1)  string
                covu        (1, :)  double = 1
                computeCov  (1, 1)  logical = false
            end

            [pose, camera, covVecP] = sonic.PnP.PnPDLT(u_C, p_I, method, false, covu, computeCov);

        end
    end

    methods (Static, Access=?tests.PnPTest) % private access except for test class

        function [pose, camera, covVecP] = PnPDLT(meas_C, points_I, method, isCalib, covMeas, computeCov)
        %% [pose, camera, covVecP] = PnPDLT(meas_C, points_I, method, isCalib, covMeas, computeCov)
        %  Solve the Perspective-n-Point using the specified method
        %
        %   References:
        %   - [1] Henry, S., & Christian, J. A. (2024).
        %      Optimal DLT-based Solutions for the Perspective-n-Point.
        %      https://doi.org/10.48550/arXiv.2410.14164
        %
        %   Inputs:
        %       - meas_C (sonic.Points2): measurement coordinates
        %       - points_I  (sonic.Points3): respective 3D position objects
        %         that are used for estimation
        %       - method (string): "DLT" (normalized DLT)
        %                       or "oDLT" (optimal DLT)
        %       - isCalib (boolean): flag whether the camera is calibrated
        %       - covMeas (OPTIONAL)
        %           - DEFAULT: unit variance
        %           - (1x1 double): isotropic variance, same for all
        %             measurements
        %           - (1xn double): isotropic variance of points2
        %       - computeCov (1x1 logical) (OPTIONAL): flag to decide
        %         whether to compute covariance or not. This will make
        %         the triangulation function slightly slower,
        %         especially for DLT and midpoint. Defaults to false.
        %
        %   Outputs:
        %       - pose (sonic.Pose): camera pose
        %       - camera (sonic.Camera): calibrated camera
        %       - covVecP (12x12 double) covariance of the vectorized 
        %       camera matrix P
        %
        %   Last revised: 11/24/24
        %   Last author: Sebastien Henry
            arguments
                meas_C      (1, 1)  sonic.Points2
                points_I    (1, 1)  sonic.Points3
                method      (1, 1)  string
                isCalib     (1, 1)  logical
                covMeas     (1, :)  double
                computeCov  (1, 1)  logical
            end

            n_meas = meas_C.n;

            if n_meas < 6
                error('sonic:PnP:inputError', ...
                    'The number of measurements should be at least 6');
            end

            if points_I.n ~= n_meas
                error('sonic:PnP:inputError', ...
                    ['The number of 3D points (p_I) should be the same' ...
                    ' as the number of measurements (u_C)']);
            end

            if length(covMeas)~=1 && length(covMeas)~=n_meas
                error('sonic:PnP:inputError', ...
                    ['The number of input covariances should be' ...
                    'either one (same for all measurements) or equal to' ...
                    'the number of measurements (one different' ...
                    ' covariance per measurement)']);
            end


            if points_I.has_inf_points
                error('sonic:PnP:inputError', ...
                    ['The 3D points (p_I) should not have points at' ...
                    ' infinity']);
            end

            if computeCov
                error('sonic:PnP:inputError', ...
                    'Computation of the covariance not yet implemented');
            end

            % Extract measurements
            meas = meas_C.p2;

            % Extract points in 3d
            points = points_I.p3;

            % Scale points for better behaviour
            [scaledMeas, measTransform, invMeasTransform] = sonic.PnP.normalizePoints(meas(1:2,:));
            [scaledPoints, pointTransform, invPointTransform] = sonic.PnP.normalizePoints(points(1:3,:));

            % Build the DLT system
            A = sonic.PnP.buildDLTSystem(scaledMeas', scaledPoints');

            switch method
                case 'oDLT'
                    % Choose a small number of points to pre-solve the
                    % PnP
                    n_small = min(n_meas, 20);

                    % We use linspace over a small number of points
                    idxSelect = linspace(1, 2*double(n_meas)-1, n_small);

                    % idxSelect+1 is the second pixel component of the
                    % selected measurements
                    idxSelect = [idxSelect, idxSelect+1];

                    % Solve for an initial camera matrix P
                    [~,~,V] = svd(A(idxSelect, :), 'econ');
                    vecInitP = V(:,end);
                    initP = reshape(vecInitP, 3, 4);

                    % Compute weights for the DLT system
                    sigu = sqrt(covMeas);
                    scale3 = initP(3,:) * points;
                    weights = 1./ (sigu.*scale3);

                    % Apply weights to the DLT system
                    A = kron(weights', [1; 1]) .* A;
                case 'DLT'
                    % Nothing special to do
            end

            % Solve the system by finding smallest eigen vector of A^T A
            [~, D, V] = svd(A, 'econ');

            % Camera matrix and inverse covariance of camera matrix
            vecP = V(:,end);
            invCovVecP = V * (D.^2) * V';

            % De-normalize P
            [vecP, invCovVecP] = sonic.PnP.denormalizeKronecker(vecP, invCovVecP, ...
                invMeasTransform, pointTransform);

            if isCalib
                camera = [];
                switch method
                    case 'DLT'
                        % Regular Procrustes Problem
                        [dcm, translation] = sonic.PnP.PToSE3(vecP);
                    case 'oDLT'
                        % Weighted Procrustes Problem
                        [dcm, ~] = sonic.PnP.PToSE3(vecP, invCovVecP);

                        % Re-triangulation for optimal translation

                        % Re-normalize the optimal DCM 
                        optNormDcm = measTransform * dcm * invPointTransform(1:3,1:3);

                        % Equivalent to LOST triangulation
                        optNormTranslation = -A(:,10:12)\(A(:,1:9)*optNormDcm(:));

                        % De-normalize the translation
                        translation = invMeasTransform * ...
                            [optNormDcm, optNormTranslation] * ...
                            pointTransform(:,end);
                end
            else
                % Decompose the camera matrix
                [dcm, translation, K] = sonic.PnP.PToCalibSE3(vecP);

                % Store calibration
                d_x = K(1,1);
                d_y = K(2,2);
                u_p = K(1,3);
                v_p = K(2,3);
                alpha = K(1,2);

                % Build SONIC Camera
                camera = sonic.Camera(d_x, d_y, alpha, u_p, v_p, sonic.Pinhole());
            end

            % Store pose
            att = sonic.Attitude(dcm);
            translation = sonic.Points3(translation);
            pose = sonic.Pose(att, translation);
     
            if computeCov
                % Store covariance of camera matrix
                covVecP = pinv(invCovVecP);
            else 
                covVecP = [];
            end
        end

        function A = buildDLTSystem(meas2d, points3d)
        %% A = buildDLTSystem(meas2d, points3d)
        % Efficiently build the DLT system, accounting for the third redundant row.
        %
        % References:
        %   - [1] S. Henry and J. A. Christian. Optimal DLT-based
        %   Solutions for the Perspective-n-Point. (2024).
        %
        % Inputs:
        %   - meas2d (nx2 double): 2D pixel (or image plane) measurements
        %   - points3d (nx3 double): 3D points in the world coordinates
        %
        % Outputs:
        %   - A (2nx12 double): Matrix incorporating the DLT constraint
        %
        % Last revised: 11/24/24
        % Last author: Sebastien Henry

            n = size(points3d, 1);

            % Extract columns
            x = points3d(:,1);
            y = points3d(:,2);
            z = points3d(:,3);
            u = meas2d(:,1);
            v = meas2d(:,2);

            % Pre-compute multiplications
            ux = u .* x;
            uy = u .* y;
            uz = u .* z;
            vx = v .* x;
            vy = v .* y;
            vz = v .* z;

            % Pre-compute matrices of zeros and ones
            zero = zeros(n, 1);
            one = ones(n, 1);

            % Fill A using vectorized operations
            A = zeros(2*n, 12);
            A(1:2:2*n,:) = [zero, -x, vx, zero, -y, vy, zero, -z, vz, zero, -one, v];
            A(2:2:2*n,:) = [x,  zero, -ux, y,  zero, -uy, z, zero, -uz, one, zero, -u];
        end

        function [vecP, invCovVecP] = denormalizeKronecker(vecP, invCovVecP, ...
                invMeasTransform, pointTransform)
        %% [vecP, invCovVecP] = denormalizeKronecker(vecP, invCovVecP, invMeasTransform, pointTransform)
        % Denormalize the vectorized camera matrix and its inverse
        % covariance.
        %
        % References:
        %   - [1] S. Henry and J. A. Christian. Optimal DLT-based
        %   Solutions for the Perspective-n-Point. (2024).
        %
        % Inputs:
        %   - vecP (12x1 double): vectorized camera matrix in the similarity space
        %   - invCovVecP (12x12 double): inverse covariance of h in the similarity space
        %   - invMeasTransform (3x3 double): inverse similarity transformation on the pixel measurements
        %   - pointTransform (4x4 double): similarity transformation on the 3D points
        %
        % Outputs:
        %   - vecP (12x1 double):  vectorized camera matrix in the regular space
        %   - invCovVecP (12x12 double): inverse covariance of h in the regular space
        %
        % Last revised: 11/24/24
        % Last author: Sebastien Henry

            % Step 1: Create the transformation matrix
            M = kron(pointTransform', invMeasTransform);

            % Step 2: Apply the transformation to P
            vecP = M * vecP;

            invM = inv(M);
            % Step 3: Transform the covariance matrix
            invCovVecP = invM' * invCovVecP * invM;
        end


        function [dcm, translation] = PToSE3(vecP, invCovVecP)
            %% [dcm, translation] = PToSE3(vecP, invCovVecP)
            % Extract the rotation matrix and translation vector
            % from an unscaled and unorthogonalized camera matrix.
            % Accounts for covariance information
            %
            % References:
            %   - [1] S. Henry and J. A. Christian. Optimal DLT-based
            %   Solutions for the Perspective-n-Point. (2024).
            %
            % Inputs:
            %   - vecP (12x1 double): unscaled and unorthogonalized vectorized
            %   camera matrix
            %   - invCovVecP (12x12 double): the inverse covariance of h
            %
            % Outputs:
            %   - dcm (3x3 double): rotation matrix from world to camera frame
            %   - translation (3x1 double): translation vector representing the camera
            %   position in the world frame
            %
            % Last revised: 11/24/24
            % Last author: Sebastien Henry

            % Extract the rotation matrix and compute the scaling
            initDcm = reshape(vecP(1:9), 3, 3);
            [~, ~, scale] = sonic.Math.A_toDet1(initDcm);

            % Scale P
            vecP = vecP/scale;
            initDcm = reshape(vecP(1:9), 3, 3);

            if nargin == 1 || rank(invCovVecP) < 12
                [U, ~, V] = svd(initDcm);
                dcm = U * diag([1, 1, det(U * V')]) * V';
            else
                W = reshape(diag(invCovVecP(1:9, 1:9)), 3, 3);

                dcm = sonic.PnP.weightedOrthogonalProcrustes(initDcm, W);

            end

            translation = vecP(10:12);
        end


        function [dcm, translation, K] = PToCalibSE3(vecP)
            %% [dcm, translation, K] = PToCalibSE3(vecP)
            % Extract the rotation matrix, translation vector, and 
            % intrinsic calibration matrix from an unscaled camera matrix.
            %
            % References:
            %   - [1] S. Henry and J. A. Christian. Optimal DLT-based
            %   Solutions for the Perspective-n-Point. (2024).
            %
            % Inputs:
            %   - vecP (12x1 double): unscaled and unorthogonalized vectorized
            %   camera matrix
            %
            % Outputs:
            %   - dcm (3x3 double): rotation matrix from world to camera frame
            %   - translation (3x1 double): translation vector representing the camera
            %   position in the world frame
            %   - K (3x3 double): Camera intrinsic calibration matrix
            %
            % Last revised: 11/24/24
            % Last author: Sebastien Henry

            P = reshape(vecP, 3, 4);
            
            invP33 = inv(P(1:3,1:3));
            Rot = [-1, 0, 0; 0, -1, 0; 0, 0, 1];
            
            % First 3x3 matrix in P is K*dcm
            % Recover the inverse rotation and inverse of K
            [invDcm, invK] = qr(invP33);

            % Recover the calibration matrix
            K = inv(invK);

            % Rescale such that the last element of the calibration matrix is 1
            scale = K(3,3);
            K = K*Rot/scale;

            % Recover extrinsics
            dcm = Rot*invDcm';
            r_I = - invP33 * P(:,4);
            translation = - dcm * r_I;
        end

        function [dcm] = weightedOrthogonalProcrustes(A, W)
        %% dcm = weightedOrthogonalProcrustes(A, W)
        % Find the closest orthogonal matrix "dcm" to "A" in SO(3)
        % with element-wise weights W. Uses linear estimation and
        % rotation vector parametrization to ensure SO(3).
        %
        % References:
        %   - [1] S. Henry and J. A. Christian. Optimal DLT-based
        %   Solutions for the Perspective-n-Point. (2024).
        %
        % Inputs:
        %   - A (3x3 double) input matrix
        %   - W (3x3 double) element-wise weights corresponding to A
        %
        % Outputs:
        %   - dcm (3x3 double): optimal rotation matrix from world to
        % camera frame, weighted closest to A in SO(3)
        %
        % Last revised: 11/24/24
        % Last author: Sebastien Henry

            % Ensure that A and W are 3x3 matrices
            assert(isequal(size(A), [3, 3]), 'A must be a 3x3 matrix.');
            assert(isequal(size(W), [3, 3]), 'W must be a 3x3 matrix.');

            Wvec = sqrt(W(:));

            % Convert the initial guess matrix A into an angle vector for optimization
            [U, ~, V] = svd(A);

            dcm0 = U * diag([1, 1, det(U * V')]) * V';

            %%% Compute linear deviation for optimal procrustes problem
            % Linearize dcm = dcm0 * (eye - crossmat(dphi))
            eps = dcm0(:) - A(:);

            % Compute the partial
            Jphi = - kron(eye(3), dcm0) * [...
                    0, 0, 0;
                    0, 0, 1;
                    0, -1, 0;
                    0, 0, -1;
                    0, 0, 0;
                    1, 0, 0;
                    0, 1, 0;
                    -1, 0, 0;
                    0, 0, 0;
                    ];

            dphi = -(Wvec.*Jphi) \ (Wvec.*eps);
            dcm = dcm0 * sonic.Attitude.rotationVecToDCM(dphi); 

        end

        function [Xtilde, transform, invTransform] = normalizePoints(X)
        %% [Xtilde, transform, invTransform] = normalizePoints(X)
        % Normalize points such that they have zero mean
        % and average distance sqrt(d) from origin, where d is the
        % dimension of the point
        %
        % References:
        %   - [1] Hartley, R., & Zisserman, A. (2003). 
        %   Multiple view geometry in computer vision.
        %   Cambridge university press.
        %
        % Inputs:
        %   - X (dxn double): points to be normalized
        %
        % Outputs:
        %   - Xtilde (dxn double): normalized points
        %   - transform (d+1xd+1 double): transform such that
        %       Xtilde = transform * [X; 1]
        %   - invTransform (d+1xd+1 double): inverse transform
        % 
        %
        % Last revised: 11/24/24
        % Last author: Sebastien Henry
            
            % Extract dimension of points
            [d, ~] = size(X);

            % Compute the centroid of the points
            centroid = mean(X, 2);

            % Translate the points so the centroid is at the origin
            translatedX = X - centroid;

            % compute the average distance of points
            averageSquaredDistance = mean(sum(translatedX.^2, 1));

            % Scale factor
            scale = sqrt(d / averageSquaredDistance);

            % Scaled points
            Xtilde = scale*translatedX;

            % Transformation from regular to scaled
            transform = [scale*eye(d), -scale*centroid;
                         zeros(1,d)  ,               1];

            % Transformation from scaled to regular
            invTransform = [1/scale*eye(d),   centroid;
                            zeros(1,d),              1];
        end
    end
end