classdef PnPTest < matlab.unittest.TestCase
    
    methods(TestClassSetup)
        % Shared setup for the entire test class
    end
    
    methods(TestMethodSetup)
        % Setup for each test
    end
    
    methods(Test)
        % Test methods
        
        function PnPCalibratedTest(testCase)
            
            % Create a random pose
            rotation_vec = 1e-1*[1; 2; 3];
            att = sonic.Attitude(rotation_vec);
            camTrans = sonic.Points3([0.5;0.5;0.5]);
            expPose = sonic.Pose(att, camTrans);

            %%% Create 6 random 3D points
            r3 = [
                1, 1, 1, -1, -1, -1;
                1, -1, 3, 1, -1, 2;
                1, 2, 3, 4, 5, 6];
            p_I = sonic.Points3(r3);

            % Project these 3D points
            x_C = expPose*p_I;
            x_C = sonic.Points2(x_C.r3);
            
            % Make sure we get the exact result, as perfect measurements
            actPose = sonic.PnP.PnPCalibrated(x_C, p_I, 'DLT');
            actPose2 = sonic.PnP.PnPCalibrated(x_C, p_I, 'oDLT');

            testCase.verifyEqual(actPose.T, expPose.T, ...
                AbsTol=sonic.Tolerances.CompZero)
            testCase.verifyEqual(actPose2.T, expPose.T, ...
                AbsTol=sonic.Tolerances.CompZero)


            %%% Noisy Measurements
            % Project these 3D points

% <<<<<<< HEAD
            % dx = 0.01*randn(2,6)
            % dx = [    
            %     0.0109   -0.0086   -0.0121   -0.0001   -0.0077   -0.0023;
            %     0.0111    0.0008   -0.0111    0.0153    0.0037    0.0112;
            %     0         0         0         0         0         0];
            % x_C = sonic.Points2(x_C.p2 + dx);
            % 
            % % Make sure we get OK results with noisy measurements
            % actPose = sonic.PnP.PnPCalibrated(x_C, p_I, 'DLT');
            % actPose2 = sonic.PnP.PnPCalibrated(x_C, p_I, 'oDLT');
% =======
            % dx = 0.001*randn(2,6)
            dx = [    
                0.00109   -0.00086   -0.00121   -0.00001   -0.00077   -0.00023;
                0.00111    0.00008   -0.00111    0.00153    0.00037    0.00112;
                0         0         0         0         0         0];
            xNoisy_C = sonic.Points2(x_C.p2 + dx);

            % Make sure we get OK results with noisy measurements
            actPose = sonic.PnP.PnPCalibrated(xNoisy_C, p_I, 'DLT');
            actPose2 = sonic.PnP.PnPCalibrated(xNoisy_C, p_I, 'oDLT');
% >>>>>>> 873d3fae56460702f1305f0ed7529614f90965e8
            
            errAngle = sonic.Attitude.angleBetween(actPose.att, expPose.att);
            errAngle2 = sonic.Attitude.angleBetween(actPose2.att, expPose.att);
            
            testCase.verifyEqual(errAngle, 0, AbsTol=0.1)
            testCase.verifyEqual(errAngle2, 0, AbsTol=0.1)

            testCase.verifyEqual(actPose.t.r3, expPose.t.r3, RelTol=0.5)
            testCase.verifyEqual(actPose2.t.r3, expPose.t.r3, RelTol=0.5)
% <<<<<<< HEAD
% 
% =======
            % 
            %%% Compute Monte-Carlo samples and make sure oDLT behave
            %%% better
            nMC = 100;
            errAngles = zeros(1, nMC);
            errAngles2 = zeros(1, nMC);
            errTranslations = zeros(3, nMC);
            errTranslations2 = zeros(3, nMC);
            xStd = 0.001;
            rng(1);
            for i = 1:nMC
                dx = [xStd; xStd; 0].*randn(3,6);
                xNoisy_C = sonic.Points2(x_C.p2 + dx);
                actPose = sonic.PnP.PnPCalibrated(xNoisy_C, p_I, 'DLT');
                actPose2 = sonic.PnP.PnPCalibrated(xNoisy_C, p_I, 'oDLT');

                errAngles(i) = sonic.Attitude.angleBetween(actPose.att, expPose.att);
                errAngles2(i) = sonic.Attitude.angleBetween(actPose2.att, expPose.att);
                errTranslations(:,i) = actPose.t.r3 - expPose.t.r3;
                errTranslations2(:,i) = actPose2.t.r3 - expPose.t.r3;
            end
            % Verify angle errors lower for oDLT
            angleRmseDLT = sqrt(mean(errAngles.^2));
            angleRmseoDLT = sqrt(mean(errAngles2.^2));
            testCase.verifyGreaterThanOrEqual(angleRmseDLT,angleRmseoDLT);

            % Verify translation errors lower for oDLT
            translationRmseDLT = sqrt(mean(norm(errTranslations)^2));
            translationRmseoDLT = sqrt(mean(norm(errTranslations2)^2));
            testCase.verifyGreaterThanOrEqual(translationRmseDLT,translationRmseoDLT);

% >>>>>>> 873d3fae56460702f1305f0ed7529614f90965e8
            %%% verify that I can put a vector of pixel variances
            sonic.PnP.PnPCalibrated(x_C, p_I, 'oDLT', ones(1, 6));
            
            testCase.verifyError( ...
                @() sonic.PnP.PnPCalibrated(x_C, p_I, 'oDLT', ones(1, 5)), ...
                'sonic:PnP:inputError')

% <<<<<<< HEAD
            %%% CANNOT YET COMPUTE COVARIANCE
% =======
            %% Cannot yet compute covariance
            testCase.verifyError( ...
                @() sonic.PnP.PnPCalibrated(x_C, p_I, 'DLT', 1, true), ...
                'sonic:PnP:inputError')
% >>>>>>> 873d3fae56460702f1305f0ed7529614f90965e8
            testCase.verifyError( ...
                @() sonic.PnP.PnPCalibrated(x_C, p_I, 'oDLT', 1, true), ...
                'sonic:PnP:inputError')

% <<<<<<< HEAD
% =======

% >>>>>>> 873d3fae56460702f1305f0ed7529614f90965e8
            %%% Test when P_I does not have the same number of points as
            % measurements
            p_I = sonic.Points3(r3(:,1:5));
            testCase.verifyError( ...
                @() sonic.PnP.PnPCalibrated(x_C, p_I, 'DLT'), ...
                'sonic:PnP:inputError')

            %%% Test when P_I has inf points
            p_I = sonic.Points3(Inf*r3);
            testCase.verifyError( ...
                @() sonic.PnP.PnPCalibrated(x_C, p_I, 'DLT'), ...
                'sonic:PnP:inputError')


            %%% Create 3 random 3D points
            p_I = sonic.Points3(r3(:,1:3));

            % Project these 3D points
            x_C = expPose*p_I;
            x_C = sonic.Points2(x_C.r3);

            % Make sure it fails, as PnP with DLT needs at least
            % 6 measurements
            testCase.verifyError( ...
                @() sonic.PnP.PnPCalibrated(x_C, p_I, 'DLT'), ...
                'sonic:PnP:inputError')
            testCase.verifyError( ...
                @() sonic.PnP.PnPCalibrated(x_C, p_I, 'oDLT'), ...
                'sonic:PnP:inputError')

        end

        function PnPUnCalibratedTest(testCase)
            % Create a random pose
            rotation_vec = 1e-1*[1; 2; 3];
            att = sonic.Attitude(rotation_vec);
            camTrans = sonic.Points3([0.5;0.5;0.5]);
            expPose = sonic.Pose(att, camTrans);

            %%% Create 6 random 3D points
            r3 = [
                1, 1, 1, -1, -1, -1;
                1, -1, 3, 1, -1, 2;
                1, 2, 3, 4, 5, 6];

            %%% Assume uncalibrated measurements now
            expCamera = sonic.Camera(480, 320, 0, 240, 160, sonic.Pinhole());
            p_I = sonic.Points3(r3);

            % Project these 3D points
            x_C = expPose*p_I;
            x_C = sonic.Points2(x_C.r3);
            u_C = sonic.Points2(expCamera.K*x_C.p2);

            [actPose, actCamera] = sonic.PnP.PnPUncalibrated(u_C, p_I, 'DLT');
            [actPose2, actCamera2] = sonic.PnP.PnPUncalibrated(u_C, p_I, 'oDLT');
            
            % Verify the pose
            testCase.verifyEqual(actPose.T, expPose.T, ...
                AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actPose2.T, expPose.T, ...
                AbsTol=sonic.Tolerances.CompZero);
            
            % Verify the camera calibration
            % Increase the tolerance slightly around 1e-10
            testCase.verifyEqual(actCamera.K, expCamera.K, ...
                AbsTol=100*sonic.Tolerances.CompZero);
            testCase.verifyEqual(actCamera2.K, expCamera.K, ...
                AbsTol=100*sonic.Tolerances.CompZero);
        end

        function normalizePointsTest(testCase)

            % Some random points
            x = [
             1.1174    0.5525    0.0859   -1.0616    0.7481   -0.7648    0.4882    1.4193    1.5877    0.8351
            -1.0891    1.1006   -1.4916    2.3505   -0.1924   -1.4023   -0.1774    0.2916   -0.8045   -0.2437
             0.0326    1.5442   -0.7423   -0.6156    0.8886   -1.4224   -0.1961    0.1978    0.6966    0.2157
             ];

            [xtilde, T, invT] = sonic.PnP.normalizePoints(x);
            % Verify mean is zero
            testCase.verifyEqual(mean(xtilde, 2), zeros(3, 1), ...
                AbsTol=sonic.Tolerances.CompZero);

            % Verify RMS distance is sqrt(d)
            testCase.verifyEqual(sqrt(mean(sum(xtilde.^2, 1))), sqrt(3), ...
                AbsTol=sonic.Tolerances.CompZero);
            
            % Verify transform
            testCase.verifyEqual([xtilde; ones(1, size(x,2))], ...
                T*[x; ones(1, size(x,2))], ...
                AbsTol=sonic.Tolerances.CompZero);

            % Verify inverse transform
            testCase.verifyEqual([x; ones(1, size(x,2))], ...
                invT*[xtilde; ones(1, size(x,2))], ...
                AbsTol=sonic.Tolerances.CompZero);
        end
    end
    
end