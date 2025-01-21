classdef PointsS2Test < matlab.unittest.TestCase

    methods (TestClassSetup)
        % Shared setup for the entire test class
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)
        % Test methods

        function PointsS2InitTest(testCase)
            %%% Nominal case
            % Specify points in degrees:
            ra_DEG = [0; 90; 0; -120];
            dec_DEG = [0; 0; 90; -45];
            
            n = length(dec_DEG);
            % Convert to radians
            exp_ra_RAD = sonic.Units.DEGtoRAD(ra_DEG);
            exp_dec_RAD = sonic.Units.DEGtoRAD(dec_DEG);
            exp_ra_dec = [exp_ra_RAD'; exp_dec_RAD'];

            % Vectors
            exp_u = [1, 0, 0, -sqrt(2)/2*1/2;
                    0, 1, 0, -sqrt(2)/2*sqrt(3)/2;
                    0, 0, 1, -sqrt(2)/2];

            pts = sonic.PointsS2(exp_ra_dec);
            pts2 = sonic.PointsS2(exp_u);
            
            testCase.verifyEqual(pts.n, uint64(n));
            testCase.verifyEqual(pts2.n, uint64(n));

            testCase.verifyEqual(pts.u, exp_u, AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(pts2.u, exp_u, AbsTol=sonic.Tolerances.CompZero);

            testCase.verifyEqual(pts.ra_dec , exp_ra_dec, AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(pts2.ra_dec, exp_ra_dec, AbsTol=sonic.Tolerances.CompZero);

            %%% Verify errors in input
            % Wrong dimension
            testCase.verifyError(@() sonic.PointsS2([1; 2; 3; 4]), ...
                'sonic:PointsS2:IncorrectDimension')
            
            % Not a unit vector
            testCase.verifyError(@() sonic.PointsS2([1; 2; 3]), ...
                'sonic:PointsS2:notUnitVectors')
        end
    end

end