classdef StarIdTest < matlab.unittest.TestCase
    properties
        max_angle = 10*pi/180; % approx FOV of camera
        max_mag = 6; % maximum visual magnitude
        att_truth % attitude of the view
        tol = 0.001*pi/180; % angular tolerance for search
        hip_cat 
        GNC_cat
        
        % 
        bin_width_cos = 1-cos(0.001);
        hipkvec
        GNCkvec
        
        % Hipparcos 
        % filtered by FOV and visual magnitude
        min_n_perfect = 5; % the minimum number of measurement allowed for the test to continue
        hip_measured_perfect
        hip_indices_measured_perfect
        hip_u_measured_world_perfect
        hip_u_measured_cam_perfect
        hip_matches_sol_perfect

        % filtered by FOV but not visual magnitude
        % this way some stars are not represented in the kvector
        hip_measured_all
        hip_indices_measured_all
        hip_u_measured_world_all
        hip_u_measured_cam_all
        hip_matches_sol_all

        % GNC
        % filtered by FOV and visual magnitude
        GNC_measured_perfect
        GNC_indices_measured_perfect
        GNC_u_measured_world_perfect
        GNC_u_measured_cam_perfect
        GNC_matches_sol_perfect

        % filtered by FOV but not visual magnitude
        % this way some stars are not represented in the kvector
        GNC_measured_all
        GNC_indices_measured_all
        GNC_u_measured_world_all
        GNC_u_measured_cam_all
        GNC_matches_sol_all

    end

    methods (TestClassSetup)
        % Shared setup for the entire test class
        function setAttitude(testCase)
            % generate a random DCM
            rotation_vec = 3*pi/180*[1/sqrt(2); 0; 1/sqrt(2)];
            testCase.att_truth = sonic.Attitude(rotation_vec);
        end

        function createHipCat(testCase)
            % evaluate the catalog
            testCase.hip_cat = sonic.Hipparcos();
            hip_ssb = testCase.hip_cat.eval(0, [0;0;0]);

            % extract stars line of sights
            u_world = hip_ssb.u;
            
            % filter all stars that are not within max_angle / 2 from center of image
            filtFov = testCase.att_truth.dcm(3,:)*u_world > cos(testCase.max_angle/2);

            % filter all stars that are within max visual magnitude
            filtVmag = testCase.hip_cat.Vmag < testCase.max_mag;
            
            % the filter that accounts for both conditions above
            filtAll = filtFov & filtVmag';

            % assert that the number of visible measurements is at least
            % min_n_perfect
            % if not, there is no hope of testing the function below
            % -> direct FAILURE OF TEST
            % need to look at another region of the sky because the test
            % may not proceed.
            testCase.assertGreaterThanOrEqual(sum(filtAll), testCase.min_n_perfect)

            % only prune the min_n_perfect first measurements for
            % tracktable test time
            % since we are going to explicitely test all triads
            filtAll = filtAll & cumsum(filtAll) <= testCase.min_n_perfect;

            %%%
            % these are stars that are measured by camera AND have
            % suffificent visual magnitude to be in the kvector
            testCase.hip_measured_perfect = testCase.hip_cat.filter(filtAll);
            testCase.hip_u_measured_world_perfect = sonic.PointsS2(u_world(:,filtAll));
            testCase.hip_indices_measured_perfect = 1:testCase.hip_cat.n;
            testCase.hip_indices_measured_perfect = testCase.hip_indices_measured_perfect(filtAll);
            % measurements are transferred to camera frame
            testCase.hip_u_measured_cam_perfect = testCase.att_truth*testCase.hip_u_measured_world_perfect;
            % solution of the attribution measurement-hipparcos_id
            testCase.hip_matches_sol_perfect = [1:testCase.hip_u_measured_world_perfect.n; testCase.hip_indices_measured_perfect];


            %%%
            % these are stars that are in the field of view BUT
            % some of them are too dim to be included in the kvector table
            testCase.hip_measured_all = testCase.hip_cat.filter(filtFov);
            testCase.hip_u_measured_world_all = sonic.PointsS2(u_world(:,filtFov));
            testCase.hip_indices_measured_all = 1:testCase.hip_cat.n;
            testCase.hip_indices_measured_all = testCase.hip_indices_measured_all(filtFov);
            % measurements are transferred to camera frame
            testCase.hip_u_measured_cam_all = testCase.att_truth*testCase.hip_u_measured_world_perfect;
            % solution of the attribution measurement-hipparcos_id
            testCase.hip_matches_sol_all = [1:testCase.hip_u_measured_world_perfect.n; testCase.hip_indices_measured_perfect];

        end

        function createGNCCat(testCase)
            % evaluate the catalog
            testCase.GNC_cat = sonic.USNOGNC();
            GNC_ssb = testCase.GNC_cat.eval(0, [0;0;0]);

            % extract stars line of sights
            u_world = GNC_ssb.u;
            
            % filter all stars that are not within max_angle / 2 from center of image
            filtFov = testCase.att_truth.dcm(3,:)*u_world > cos(testCase.max_angle/2);

            % filter all stars that are within max magnitude
            filtmag = testCase.GNC_cat.gmag < testCase.max_mag;
            
            % the filter that accounts for both conditions above
            filtAll = filtFov & filtmag';

            % assert that the number of visible measurements is at least
            % min_n_perfect
            % if not, there is no hope of testing the function below
            % -> direct FAILURE OF TEST
            % need to look at another region of the sky because the test
            % may not proceed.
            testCase.assertGreaterThanOrEqual(sum(filtAll), testCase.min_n_perfect)

            % only prune the min_n_perfect first measurements for
            % tracktable test time
            % since we are going to explicitely test all triads
            filtAll = filtAll & cumsum(filtAll) <= testCase.min_n_perfect;

            %%%
            % these are stars that are measured by camera AND have
            % suffificent visual magnitude to be in the kvector
            testCase.GNC_measured_perfect = testCase.GNC_cat.filter(filtAll);
            testCase.GNC_u_measured_world_perfect = sonic.PointsS2(u_world(:,filtAll));
            testCase.GNC_indices_measured_perfect = 1:testCase.GNC_cat.n;
            testCase.GNC_indices_measured_perfect = testCase.GNC_indices_measured_perfect(filtAll);
            % measurements are transferred to camera frame
            testCase.GNC_u_measured_cam_perfect = testCase.att_truth*testCase.GNC_u_measured_world_perfect;
            % solution of the attribution measurement-hipparcos_id
            testCase.GNC_matches_sol_perfect = [1:testCase.GNC_u_measured_world_perfect.n; testCase.GNC_indices_measured_perfect];


            %%%
            % these are stars that are in the field of view BUT
            % some of them are too dim to be included in the kvector table
            testCase.GNC_measured_all = testCase.GNC_cat.filter(filtFov);
            testCase.GNC_u_measured_world_all = sonic.PointsS2(u_world(:,filtFov));
            testCase.GNC_indices_measured_all = 1:testCase.GNC_cat.n;
            testCase.GNC_indices_measured_all = testCase.GNC_indices_measured_all(filtFov);
            % measurements are transferred to camera frame
            testCase.GNC_u_measured_cam_all = testCase.att_truth*testCase.GNC_u_measured_world_perfect;
            % solution of the attribution measurement-hipparcos_id
            testCase.GNC_matches_sol_all = [1:testCase.GNC_u_measured_world_perfect.n; testCase.GNC_indices_measured_perfect];

        end
        
        function createKvecs(testCase)
            % max visual magnitude for Kvector
            testCase.hipkvec = sonic.Kvector(testCase.hip_cat, 0, testCase.max_angle, testCase.max_mag, testCase.bin_width_cos);
            testCase.GNCkvec = sonic.Kvector(testCase.GNC_cat, 0, testCase.max_angle, testCase.max_mag, testCase.bin_width_cos);
        end

    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)
        % Test methods

        function hipInterstarAngleTest(testCase)
    
            %%% test with perfect measurements in FOV: 
            %%% only the ones we could find in the kvector
            [matches, att_est] = ...
                sonic.StarId.interstarAngle(testCase.hipkvec, testCase.hip_cat, ...
                testCase.hip_u_measured_cam_perfect, testCase.tol, testCase.max_angle);
            
            err_att = sonic.Attitude.angleBetween(att_est, testCase.att_truth);
            testCase.verifyEqual(err_att, 0, sonic.Tolerances.SmallAngle);
            % should have recovered all stars in kvector
            testCase.verifyEqual(matches, testCase.hip_matches_sol_perfect, sonic.Tolerances.CompZero);

            %%% test with all measurements in fov, some of which are not 
            %%% included in kvector
            [matches, att_est] = ...
                sonic.StarId.interstarAngle(testCase.hipkvec, testCase.hip_cat, ...
                testCase.hip_u_measured_cam_all, testCase.tol, testCase.max_angle);
            
            err_att = sonic.Attitude.angleBetween(att_est, testCase.att_truth);
            testCase.verifyEqual(err_att, 0, sonic.Tolerances.SmallAngle);
            % should have recovered all stars in kvector, and all stars in
            % FOV
            testCase.verifyEqual(matches, testCase.hip_matches_sol_all, sonic.Tolerances.CompZero);

            %%% test input sanity checks
            testCase.verifyError(...
                @() sonic.StarId.interstarAngle(testCase.hipkvec, testCase.hip_cat, ...
                testCase.hip_u_measured_cam_all, testCase.tol, testCase.max_angle, 3), ...
                'sonic:StarId:minMatchesTooFew')

            testCase.verifyError(...
                @() sonic.StarId.interstarAngle(testCase.hipkvec, testCase.hip_measured_perfect, ...
                testCase.hip_u_measured_cam_all, testCase.tol, testCase.max_angle, 4), ...
                'sonic:StarId:nonBaseCatalogEntered')

        end

        function GNCInterstarAngleTest(testCase)
    
            %%% test with perfect measurements in FOV: 
            %%% only the ones we could find in the kvector
            [matches, att_est] = ...
                sonic.StarId.interstarAngle(testCase.GNCkvec, testCase.GNC_cat, ...
                testCase.GNC_u_measured_cam_perfect, testCase.tol, testCase.max_angle);
            
            err_att = sonic.Attitude.angleBetween(att_est, testCase.att_truth);
            testCase.verifyEqual(err_att, 0, sonic.Tolerances.SmallAngle);
            % should have recovered all stars in kvector
            testCase.verifyEqual(matches, testCase.GNC_matches_sol_perfect, sonic.Tolerances.CompZero);

            %%% test with all measurements in fov, some of which are not 
            %%% included in kvector
            [matches, att_est] = ...
                sonic.StarId.interstarAngle(testCase.GNCkvec, testCase.GNC_cat, ...
                testCase.GNC_u_measured_cam_all, testCase.tol, testCase.max_angle);
            
            err_att = sonic.Attitude.angleBetween(att_est, testCase.att_truth);
            testCase.verifyEqual(err_att, 0, sonic.Tolerances.SmallAngle);
            % should have recovered all stars in kvector, and all stars in
            % FOV
            testCase.verifyEqual(matches, testCase.GNC_matches_sol_all, sonic.Tolerances.CompZero);

            %%% test input sanity checks
            testCase.verifyError(...
                @() sonic.StarId.interstarAngle(testCase.GNCkvec, testCase.GNC_cat, ...
                testCase.GNC_u_measured_cam_all, testCase.tol, testCase.max_angle, 3), ...
                'sonic:StarId:minMatchesTooFew')

            testCase.verifyError(...
                @() sonic.StarId.interstarAngle(testCase.GNCkvec, testCase.GNC_measured_perfect, ...
                testCase.GNC_u_measured_cam_all, testCase.tol, testCase.max_angle, 4), ...
                'sonic:StarId:nonBaseCatalogEntered')

        end

        function hipCheckTriadTest(testCase)            
            % needs to be a fresh hipparcos catalog!!!
            hip_cat_idx = 1:testCase.hip_cat.n;
            pointsS2_cat = testCase.hip_cat.eval(testCase.hipkvec.et, testCase.hipkvec.r_obs_AU);
            
            %%% Make sure of perfect result for any triad since perfect
            % measurements
            n = testCase.hip_u_measured_cam_perfect.n;
            for i = 1:n
                for j = 1:n
                    for k = 1:n
                        if i~=j && i~=k && j~=k
                            
                           [matches, att_est] = sonic.StarId.checkTriad( ...
                           testCase.hipkvec, hip_cat_idx, testCase.hip_u_measured_cam_perfect, pointsS2_cat, ...
                                [i,j,k], testCase.tol, testCase.max_angle, testCase.min_n_perfect);
                            % verify that the matches are perfect
                            testCase.verifyEqual(matches, testCase.hip_matches_sol_perfect);

                            % verify that the attitude is perfect
                            err_att = sonic.Attitude.angleBetween(att_est, testCase.att_truth);
                            testCase.verifyEqual(err_att, 0, sonic.Tolerances.SmallAngle);
                        end
                    end
                end
            end
            
            %%% now only work with 5 measurements

            %%% verify that it works on 5 perfect measurements
            u_measured_cam_5 = sonic.PointsS2(testCase.hip_u_measured_cam_perfect.u(:,1:5));
            
            i = 1; j = 2; k = 3;
            [matches, att_est] = sonic.StarId.checkTriad( ...
                testCase.hipkvec, hip_cat_idx, u_measured_cam_5, pointsS2_cat, ...
                [i,j,k], testCase.tol, testCase.max_angle);

            % verify that the matches are perfect
            testCase.verifyEqual(matches, testCase.hip_matches_sol_perfect);

            % verify that the attitude is perfect
            err_att = sonic.Attitude.angleBetween(att_est, testCase.att_truth);
            testCase.verifyEqual(err_att, 0, sonic.Tolerances.SmallAngle);
            
            %%% verify that it does not work with 4 good measurements and a
            %%% minimum of 5 matches required
      
            % make a random last measurement
            u_measured_cam_5_noisy = sonic.PointsS2(...
                [u_measured_cam_5.u(:,1:4), sonic.PointsS2([pi/4; pi/5]).u]);
            
            i = 1; j = 2; k = 3;
            [~, ~, isValid] = sonic.StarId.checkTriad( ...
                testCase.hipkvec, hip_cat_idx, u_measured_cam_5_noisy, pointsS2_cat, ...
                [i,j,k], testCase.tol, testCase.max_angle, 5);
            testCase.verifyFalse(isValid);

            %%% verify that it does work with 4 good measurements and a
            %%% minimum of 4 matches required
            i = 1; j = 2; k = 3;
            [~, att_est, isValid] = sonic.StarId.checkTriad( ...
                testCase.hipkvec, hip_cat_idx, u_measured_cam_5_noisy, pointsS2_cat, ...
                [i,j,k], testCase.tol, testCase.max_angle, 4);
            
            % verify that the attitude is perfect
            err_att = sonic.Attitude.angleBetween(att_est, testCase.att_truth);
            testCase.verifyEqual(err_att, 0, sonic.Tolerances.SmallAngle);
            testCase.verifyTrue(isValid);

            %%% verify there is no match if I put a bad measurement in the
            %%% triad
            i = 1; j = 2; k = 5; % <- k is the bad measurement
            [matches, att_est] = sonic.StarId.checkTriad( ...
              testCase.hipkvec, hip_cat_idx, u_measured_cam_5_noisy, pointsS2_cat, ...
              [i,j,k], testCase.tol, testCase.max_angle, 4);

            % verify that no match has been made with one wrong index
            % in the triad
            testCase.verifyEmpty(matches);
            testCase.verifyEmpty(att_est);


        end

        function GNCCheckTriadTest(testCase)            
            % needs to be a fresh hipparcos catalog!!!
            GNC_cat_idx = 1:testCase.GNC_cat.n;
            pointsS2_cat = testCase.GNC_cat.eval(testCase.GNCkvec.et, testCase.GNCkvec.r_obs_AU);
            
            %%% Make sure of perfect result for any triad since perfect
            % measurements
            n = testCase.GNC_u_measured_cam_perfect.n;
            for i = 1:n
                for j = 1:n
                    for k = 1:n
                        if i~=j && i~=k && j~=k
                            
                           [matches, att_est] = sonic.StarId.checkTriad( ...
                           testCase.GNCkvec, GNC_cat_idx, testCase.GNC_u_measured_cam_perfect, pointsS2_cat, ...
                                [i,j,k], testCase.tol, testCase.max_angle, testCase.min_n_perfect);
                            % verify that the matches are perfect
                            testCase.verifyEqual(matches, testCase.GNC_matches_sol_perfect);

                            % verify that the attitude is perfect
                            err_att = sonic.Attitude.angleBetween(att_est, testCase.att_truth);
                            testCase.verifyEqual(err_att, 0, sonic.Tolerances.SmallAngle);
                        end
                    end
                end
            end
            
            %%% now only work with 5 measurements

            %%% verify that it works on 5 perfect measurements
            u_measured_cam_5 = sonic.PointsS2(testCase.GNC_u_measured_cam_perfect.u(:,1:5));
            
            i = 1; j = 2; k = 3;
            [matches, att_est] = sonic.StarId.checkTriad( ...
                testCase.GNCkvec, GNC_cat_idx, u_measured_cam_5, pointsS2_cat, ...
                [i,j,k], testCase.tol, testCase.max_angle);

            % verify that the matches are perfect
            testCase.verifyEqual(matches, testCase.GNC_matches_sol_perfect);

            % verify that the attitude is perfect
            err_att = sonic.Attitude.angleBetween(att_est, testCase.att_truth);
            testCase.verifyEqual(err_att, 0, sonic.Tolerances.SmallAngle);
            
            %%% verify that it does not work with 4 good measurements and a
            %%% minimum of 5 matches required
      
            % make a random last measurement
            u_measured_cam_5_noisy = sonic.PointsS2(...
                [u_measured_cam_5.u(:,1:4), sonic.PointsS2([pi/4; pi/5]).u]);
            
            i = 1; j = 2; k = 3;
            [~, ~, isValid] = sonic.StarId.checkTriad( ...
                testCase.GNCkvec, GNC_cat_idx, u_measured_cam_5_noisy, pointsS2_cat, ...
                [i,j,k], testCase.tol, testCase.max_angle, 5);
            testCase.verifyFalse(isValid);

            %%% verify that it does work with 4 good measurements and a
            %%% minimum of 4 matches required
            i = 1; j = 2; k = 3;
            [~, att_est, isValid] = sonic.StarId.checkTriad( ...
                testCase.GNCkvec, GNC_cat_idx, u_measured_cam_5_noisy, pointsS2_cat, ...
                [i,j,k], testCase.tol, testCase.max_angle, 4);
            
            % verify that the attitude is perfect
            err_att = sonic.Attitude.angleBetween(att_est, testCase.att_truth);
            testCase.verifyEqual(err_att, 0, sonic.Tolerances.SmallAngle);
            testCase.verifyTrue(isValid);

            %%% verify there is no match if I put a bad measurement in the
            %%% triad
            i = 1; j = 2; k = 5; % <- k is the bad measurement
            [matches, att_est] = sonic.StarId.checkTriad( ...
              testCase.GNCkvec, GNC_cat_idx, u_measured_cam_5_noisy, pointsS2_cat, ...
              [i,j,k], testCase.tol, testCase.max_angle, 4);

            % verify that no match has been made with one wrong index
            % in the triad
            testCase.verifyEmpty(matches);
            testCase.verifyEmpty(att_est);

        end

    end

end