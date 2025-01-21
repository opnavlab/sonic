classdef Points3Test < matlab.unittest.TestCase
    
    methods(TestClassSetup)
        % Shared setup for the entire test class
    end
    
    methods(TestMethodSetup)
        % Setup for each test
    end
    
    methods(Test)
        % Test methods
        
        function Points3InitTest(testCase)
            %%% Points not at infinity
            % Create points 3 with r3            
            expPointsR3 = [1, 2, 3, 4;
                           5, 6, 7, 8;
                           9, 10, 11, 12];

            n = length(expPointsR3);
            % Create the same points using p3
            expPointsP3 = [expPointsR3; ones(1, n)];

            sonicPoints = sonic.Points3(expPointsR3);
            sonicPoints_2 = sonic.Points3(expPointsP3);
            
            testCase.verifyEqual(sonicPoints.r3, expPointsR3)
            testCase.verifyEqual(sonicPoints_2.r3, expPointsR3)

            testCase.verifyEqual(sonicPoints.p3, expPointsP3)
            testCase.verifyEqual(sonicPoints_2.p3, expPointsP3)

            testCase.verifyFalse(sonicPoints.has_inf_points)
            testCase.verifyFalse(sonicPoints_2.has_inf_points)

            testCase.verifyEqual(sonicPoints.inf_points, false(1,n))
            testCase.verifyEqual(sonicPoints_2.inf_points, false(1,n))

            testCase.verifyEqual(sonicPoints.n, uint64(n))
            testCase.verifyEqual(sonicPoints_2.n, uint64(n))

            %%% Now a mix of points at infinity
            fourthRow = [0, 1, 0, 1];
            expPointsP3 = [expPointsR3; fourthRow];
            sonicPoints_mixed = sonic.Points3(expPointsP3);
            
            testCase.verifyTrue(~allfinite(sonicPoints_mixed.r3))
            testCase.verifyEqual(sonicPoints_mixed.p3, expPointsP3)
            testCase.verifyTrue(sonicPoints_mixed.has_inf_points)
            testCase.verifyEqual(sonicPoints_mixed.inf_points, ~logical(fourthRow))

            % test the getters now
            testCase.verifyEqual(sonicPoints_mixed.r3_finite, expPointsP3(1:3, logical(fourthRow)))
            testCase.verifyEqual(sonicPoints_mixed.p3_finite, expPointsP3(:, logical(fourthRow)))
            testCase.verifyEqual(sonicPoints_mixed.p3_infinite, expPointsP3(:, ~fourthRow))

            %%% Wrong input
            testCase.verifyError(@() sonic.Points3([1;2]), 'sonic:Points3:IncorrectDimension')
        end
    end
    
end