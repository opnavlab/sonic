classdef PoseTest < matlab.unittest.TestCase
    
    methods(TestClassSetup)
        % Shared setup for the entire test class
    end
    
    methods(TestMethodSetup)
        % Setup for each test
    end
    
    methods(Test)
        % Test methods
        
        function PoseInitTest(testCase)
            % 
            att = sonic.Attitude([1; 2; 3]);
            trans = sonic.Points3([1; 2; 3]);
            expT = [att.dcm, trans.r3; 0, 0, 0, 1];
            
            %%% check that we can create a pose from position and attitude
            pose = sonic.Pose(att, trans);
            
            testCase.verifyEqual(pose.att.dcm, att.dcm, AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(pose.t.r3, trans.r3, AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(pose.T, expT)

            %%% check that we can create a pose from a 4x4 matrix
            pose = sonic.Pose(expT);

            testCase.verifyEqual(pose.att.dcm, att.dcm, AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(pose.t.r3, trans.r3, AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(pose.T, expT)
            
            %%% check that it throws an error if multiple positions
            testCase.verifyError(...
                @() sonic.Pose(att, sonic.Points3([1 2;3 4; 5 6])), ...
                'sonic:Pose:invalidPoints3')

            %%% what about if translation is at infinity?
            testCase.verifyError(...
                @() sonic.Pose(att, sonic.Points3([1; 2; 3; 0])), ...
                'sonic:Pose:invalidPoints3')
            
            %%% what about if T is not valid
            nonValidT1 = expT;
            nonValidT1(4,3) = 1;
            testCase.verifyError(...
                @() sonic.Pose(nonValidT1), ...
                'sonic:Pose:invalidSE3')

            nonValidT2 = expT;
            nonValidT2(4,4) = 1.5;
            testCase.verifyError(...
                @() sonic.Pose(nonValidT2), ...
                'sonic:Pose:invalidSE3')

            testCase.verifyError(...
                @() sonic.Pose([1, 2]), ...
                'sonic:Pose:invalidSize')

            %%% what if wrong objects
            testCase.verifyError(...
                @() sonic.Pose(1, sonic.Points3([1; 2; 3; 0])), ...
                'sonic:Pose:invalidInputType')

            testCase.verifyError(...
                @() sonic.Pose(att, 1), ...
                'sonic:Pose:invalidInputType')
        end

        function compTest(testCase)
            att = sonic.Attitude([1; 2; 3]);
            trans = sonic.Points3([1; 2; 3]);

            pose = sonic.Pose(att, trans);
            newPose = pose.comp(pose);

            testCase.verifyEqual(newPose.T, pose.T*pose.T)
        end

        function invTest(testCase)
            att = sonic.Attitude([1; 2; 3]);
            trans = sonic.Points3([1; 2; 3]);

            pose = sonic.Pose(att, trans);
            
            % expect identity for the compose of the inverse
            expT = eye(4);
            actT = pose.comp(pose.inv()).T;

            testCase.verifyEqual(actT, expT, AbsTol=sonic.Tolerances.CompZero);

        end

        function transformTest(testCase)
            att = sonic.Attitude([1; 2; 3]);
            trans = sonic.Points3([1; 2; 3]);

            pose = sonic.Pose(att, trans);
            
            % test with non-infinite points
            points = sonic.Points3([1 2; 3 4; 5 6; 7 8]);
            transformedPoints = pose.transform(points);

            expr3 = att.rotate(points).r3 + trans.r3;
            testCase.verifyEqual(transformedPoints.r3, expr3, AbsTol=sonic.Tolerances.CompZero)

            % test with mix of finite and infinite points
            points = sonic.Points3([1 2; 3 4; 5 6; 7 0]);
            transformedPoints = pose.transform(points);
            

            expp3 = zeros(4, 2);
            % the point not at infinity should have the same transformation
            expp3(:,1) = sonic.Points3(att.rotate(points).r3 + trans.r3).p3(:,1);

            % the point at infinity should only have the rotation
            % transformation
            expp3(:,2) = att.rotate(points).p3(:,2);

            testCase.verifyEqual(transformedPoints.p3, expp3, AbsTol=sonic.Tolerances.CompZero)
        end


        function mtimesTest(testCase)
            att = sonic.Attitude([1; 2; 3]);
            trans = sonic.Points3([1; 2; 3]);
            pose = sonic.Pose(att, trans);
            
            %%% Test the pose * pose multiplication
            newPose = pose.comp(pose);
            newPose2 = pose * pose;

            testCase.verifyEqual(newPose.T, newPose2.T, AbsTol=sonic.Tolerances.CompZero)

           %%% Test the pose * point multiplication
           points = sonic.Points3([1 2; 3 4; 5 6; 7 0]);
           newPoints = pose.transform(points);
           newPoints2 = pose * points;

           testCase.verifyEqual(newPoints.p3, newPoints2.p3, AbsTol=sonic.Tolerances.CompZero)

           %%% Test with a wrong object
           testCase.verifyError(@() pose * 1, 'sonic:Pose:mtimes:invalidInputType')
        end
    end
    
end