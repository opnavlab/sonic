classdef KvectorTest < matlab.unittest.TestCase
    properties
        hip_cat
        min_angle = 0.5*pi/180; % 0.5 degrees
        max_angle = 3 * pi / 180; % 3 degrees
        max_Vmag = 6;
        bin_width_cos_angle = 1-cos(0.001);
    end

    methods(TestClassSetup)
        % Shared setup for the entire test class
        function createStarCat(testCase)
            % fresh star catalog
            testCase.hip_cat = sonic.Hipparcos();
        end
    end
    
    methods(TestMethodSetup)
        % Setup for each test
    end
    
    methods(Test)
        % Test methods
        function KvectorInitTest(testCase)
            et = 13.7;
            r_obs_AU = [2.07;1;4];
            
            % if the bin width is provided, then there should be a
            % different number of bins and stars in the kvec
            kvec = sonic.Kvector(testCase.hip_cat, testCase.min_angle, ...
                testCase.max_angle, testCase.max_Vmag, ...
                testCase.bin_width_cos_angle, et, r_obs_AU);
            testCase.verifyNotEqual(kvec.n, kvec.n_bins)
            
            % if the bin width is zero, then should have one bin per item
            kvec = sonic.Kvector(testCase.hip_cat, testCase.min_angle, ...
                testCase.max_angle, testCase.max_Vmag, ...
                0, et, r_obs_AU);
            testCase.verifyEqual(kvec.n, kvec.n_bins)

            % no bin width is the same as bin_width = 0
            kvec = sonic.Kvector(testCase.hip_cat, testCase.min_angle, ...
                testCase.max_angle, testCase.max_Vmag);
            testCase.verifyEqual(kvec.n, kvec.n_bins)


            % if no et and r_obs_AU, then default to 0
            kvec = sonic.Kvector(testCase.hip_cat, testCase.min_angle, ...
                testCase.max_angle, testCase.max_Vmag, ...
                0);
            testCase.verifyEqual(kvec.et, 0)
            testCase.verifyEqual(kvec.r_obs_AU, [0;0;0])

            % a fresh catalog should be included
            testCase.verifyError(...
                @() sonic.Kvector(testCase.hip_cat.filter(1:3), testCase.min_angle, testCase.max_angle, testCase.max_Vmag, 0), ...
                'sonic:StarId:nonBaseCatalogEntered');

        end

        function queryTest(testCase)
            et = 0;
            r_obs_AU = [0;0;0];
            kvec = sonic.Kvector(testCase.hip_cat, testCase.min_angle, ...
                testCase.max_angle, testCase.max_Vmag, ...
                testCase.bin_width_cos_angle, et, r_obs_AU);

            min_query_angle = 1.3*pi/180;
            max_query_angle = 2.3*pi/180;
            [min_idx, max_idx] = kvec.query(min_query_angle, max_query_angle);

            % verify that the angles are approximately right
            %%% minimum angle
            actual = max(min_query_angle, acos(kvec.cos_interstar_angle(end)));
            queried = acos(kvec.cos_interstar_angle(max_idx));

            % queried min angle should be less or equal to actual min
            testCase.verifyLessThanOrEqual(queried, actual)

            % verify the angles are pretty close, should be within one bin
            % width
            testCase.verifyEqual(queried, actual, AbsTol=acos(testCase.bin_width_cos_angle))

            %%% minimum angle
            actual = max_query_angle;
            queried = acos(kvec.cos_interstar_angle(min_idx));

            % queried min angle should be less or equal to actual min
            testCase.verifyGreaterThanOrEqual(queried, actual)

            % verify the angles are pretty close, should be within one bin
            % width
            testCase.verifyEqual(queried, actual, AbsTol=acos(testCase.bin_width_cos_angle))


            %%% some input sanity checks
            % should not have a valid angle bracket
            testCase.verifyError( ...
                @() kvec.query(max_query_angle, min_query_angle), ...
                'sonic:KvectorError:invalidInput')

            % if the min angle is greater than the kvector max angle
            [min_idx, max_idx] = kvec.query(testCase.max_angle+0.001, testCase.max_angle+0.002);
            testCase.verifyEmpty(min_idx)
            testCase.verifyEmpty(max_idx)

            % if the max angle is smaller than the kvector min angle
            [min_idx, max_idx] = kvec.query(testCase.min_angle-0.002, testCase.min_angle-0.001);
            testCase.verifyEmpty(min_idx)
            testCase.verifyEmpty(max_idx)

            % if we put exactly the min and max angle
            [min_idx, max_idx] = kvec.query(testCase.min_angle, testCase.max_angle);
            testCase.verifyEqual(min_idx, 1);
            testCase.verifyEqual(int64(max_idx), kvec.n)
        end
    end
    
end