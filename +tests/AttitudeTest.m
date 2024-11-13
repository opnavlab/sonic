classdef AttitudeTest < matlab.unittest.TestCase

    methods (TestClassSetup)
        % Shared setup for the entire test class
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)
        % Test methods
        function AttitudeInitTest(testCase)
            %%% create attitude from rotation angle and axis
            expAngle1 = 27*pi/180;
            expAx1 = [1/sqrt(2); 0; -1/sqrt(2)];
            
            % create attitude from rotation vector
            rotation_vec = sonic.Attitude.axAngToRotationVec(expAx1, expAngle1);
            att1 = sonic.Attitude(rotation_vec);
            
            % create attitude from quaternion
            [qv, qs] = sonic.Attitude.axAngToQuat(expAx1, expAngle1);
            att2 = sonic.Attitude([qs; qv], 'sv');

            % create attitude from DCM
            T = sonic.Attitude.axAngToDCM(expAx1, expAngle1);
            att3 = sonic.Attitude(T);

            % verify DCM
            testCase.verifyEqual(att1.dcm, att2.dcm, AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(att1.dcm, att3.dcm, AbsTol=sonic.Tolerances.CompZero);

            % verify quaternion
            testCase.verifyEqual(att1.quat_v, att2.quat_v, AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(att1.quat_v, att3.quat_v, AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(att1.quat_s, att2.quat_s, AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(att1.quat_s, att3.quat_s, AbsTol=sonic.Tolerances.CompZero);

            % verify rotationvector
            testCase.verifyEqual(att1.rotation_vec, att2.rotation_vec, AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(att1.rotation_vec, att3.rotation_vec, AbsTol=sonic.Tolerances.CompZero);


            %%% verify that invalid inputs yield to errors
            % what if I input a vectorized rotation matrix?
            testCase.verifyError(@() sonic.Attitude(T(:)), 'sonic:Attitude:invalidInput');

            % what if I input a transposed rotation vector?
            testCase.verifyError(@() sonic.Attitude(rotation_vec'), 'sonic:Attitude:invalidInput');

            % what if I input a transposed quaternion?
            testCase.verifyError(@() sonic.Attitude([qs; qv]'), 'sonic:Attitude:invalidInput');

            % what if I input a vector that does not have the right size?
            testCase.verifyError(@() sonic.Attitude(ones(5,1)), 'sonic:Attitude:invalidInput')

            % verify that the quat convention needs to be set when
            % quaternion
            testCase.verifyError(@() sonic.Attitude([qs; qv]), 'sonic:Attitude:enterQuatConvention')

            % verify that the quat convention needs to be set when
            % quaternion
            testCase.verifyError(@() sonic.Attitude([qs; qv], 'svq'), 'sonic:Attitude:incorrectQuatConvention')

            % verify that a non valid DCM throws an error
            testCase.verifyError(@() sonic.Attitude(diag([2, 1, 1])), 'sonic:Attitude:verifyValidDCM:DetNot1');

            % verify that a non valid DCM with det 1 but not orthogonal
            dcm = [1,0,1; 0,1,0; 0,0,1];
            testCase.verifyEqual(det(dcm), 1, AbsTol=sonic.Tolerances.CompZero)
            testCase.verifyError(@() sonic.Attitude(dcm), 'sonic:Attitude:verifyValidDCM:NotOrthog');

            % verify that a quaternion not of norm 1 throws an error
            testCase.verifyError(@() sonic.Attitude([2; 1; 1; 3]), 'sonic:Attitude:verifyValidQuat:NormNot1');

            % verify that a rotation vector with NaN or Inf does not
            % work
            testCase.verifyError(@() sonic.Attitude([1; nan; 2]), 'sonic:Attitude:verifyValidRotationVec:NotFinite');
            testCase.verifyError(@() sonic.Attitude([1; inf; 2]), 'sonic:Attitude:verifyValidRotationVec:NotFinite');
            testCase.verifyError(@() sonic.Attitude([1; 3+1i; 2]), 'sonic:Attitude:verifyValidRotationVec:NotReal');
        end

        function compTest(testCase)
            %%% a 90 deg attitude
            rotation_vec = pi/2 * [1/sqrt(3); 1/sqrt(3); 1/sqrt(3)];
            att = sonic.Attitude(rotation_vec);

            actAtt = att.comp(att);
            expDCM = att.dcm*att.dcm;

            testCase.verifyEqual(actAtt.dcm, expDCM, AbsTol=sonic.Tolerances.CompZero);
            
            % 4 rotations of 90 degrees should lead to the same attitude
            fullTour = att.comp(att).comp(att).comp(att).comp(att);
            testCase.verifyEqual(att.dcm, fullTour.dcm, AbsTol=sonic.Tolerances.CompZero);
        end

        function invTest(testCase)
            %%% the attitude does not really matter here
            rotation_vec = [47; 32; 3];
            att = sonic.Attitude(rotation_vec);
            
            % verify that the compose of the inverse is gives no rotation
            actAtt = att * att.inv();
            expDCM = eye(3);

            testCase.verifyEqual(actAtt.dcm, expDCM, AbsTol=sonic.Tolerances.CompZero);
       end


        function solveWahbasProblemTest(testCase)
            %%% test with a couple of points

            % create an attitude matrix
            expAngle1 = 42*pi/180;
            expAx1 = [1/sqrt(2); 0; -1/sqrt(2)];
            rotation_vec = sonic.Attitude.axAngToRotationVec(expAx1, expAngle1);
            expAtt = sonic.Attitude(rotation_vec);

            % create a couple pointsS2 in world
            RAs = [5, 30] * pi / 180;
            DECs = [24, 21] * pi / 180;
            pointsS2_world = sonic.PointsS2([RAs; DECs]);

            % rotate world points in camera frame
            pointsS2_cam = expAtt*pointsS2_world;

            % solve Wahba's problem
            actAtt = sonic.Attitude.solveWahbasProblem(pointsS2_cam, pointsS2_world);

            % expect the attitudes to be the same
            error = sonic.Attitude.angleBetween(expAtt, actAtt);
            testCase.verifyEqual(error, 0, AbsTol=sonic.Tolerances.CompZero);

            %%% test the input safety checks
            testCase.verifyError(@() ...
                sonic.Attitude.solveWahbasProblem(sonic.PointsS2([0;1]), pointsS2_world), ...
                'sonic:Attitude:solveWahbasProblem:badSizes');

            testCase.verifyError(@() ...
                sonic.Attitude.solveWahbasProblem(pointsS2_world, sonic.PointsS2([0;1])), ...
                'sonic:Attitude:solveWahbasProblem:badSizes');

            RAs = [5, 30, 31] * pi / 180;
            DECs = [24, 21, 24] * pi / 180;
            pointsS2_world_large = sonic.PointsS2([RAs; DECs]);
            testCase.verifyError(@() ...
                sonic.Attitude.solveWahbasProblem(pointsS2_cam, pointsS2_world_large), ...
                'sonic:Attitude:solveWahbasProblem:badSizes');

        end

        function angleBetweenTest(testCase)
            %%% create a random attitude, and another one that is the same
            rotation_vec = 27*pi/180 * [1/sqrt(2); 0; -1/sqrt(2)];
            att1 = sonic.Attitude(rotation_vec);
            att2 = sonic.Attitude(att1.dcm);

            % expect zero angle
            expAngle1 = 0;
            actAngle1 = sonic.Attitude.angleBetween(att1, att2);
            testCase.verifyEqual(expAngle1, actAngle1, AbsTol=sonic.Tolerances.SmallAngle);

            % expect a scalar as an output
            testCase.verifySize(actAngle1, [1,1]);

            %%% create two attitudes that are 90 deg apart
            att1 = sonic.Attitude(eye(3));
            att3 = sonic.Attitude([1, 0, 0; 0, 0, -1; 0, 1, 0]);

            % expect 90 deg angle (pi/2 rad)
            expAngle2 = pi/2;
            actAngle2 = sonic.Attitude.angleBetween(att1, att3);
            testCase.verifyEqual(expAngle2, actAngle2, AbsTol=sonic.Tolerances.SmallAngle);

            % expect the same result if attitudes are switched in input
            actAngle3 = sonic.Attitude.angleBetween(att3, att1);
            testCase.verifyEqual(expAngle2, actAngle3, AbsTol=sonic.Tolerances.SmallAngle);
        end

        function rotationVecToDCMTest(testCase)
            %%% create a random rotation angle and axis
            expAngle = 27*pi/180;
            expAx = [1/sqrt(2); 0; -1/sqrt(2)];
            rotation_vec1 = expAngle*expAx;

            expDCM = expm(-sonic.Math.crossmat(rotation_vec1));
            actDCM = sonic.Attitude.rotationVecToDCM(rotation_vec1);

            testCase.verifyEqual(actDCM, expDCM, AbsTol=sonic.Tolerances.CompZero);

            %%% what if I add 2*pi to the angle with same axis? It should be the same.
            expAngle2 = 27*pi/180 + 2*pi;
            rotation_vec2 = expAngle2*expAx;

            actDCM2 = sonic.Attitude.rotationVecToDCM(rotation_vec2);

            testCase.verifyEqual(actDCM2, expDCM, AbsTol=sonic.Tolerances.CompZero);

            %%% create a zero vector. It should be identity.
            rotation_vec3 = [0;0;0];

            expDCM3 = eye(3);
            actDCM3 = sonic.Attitude.rotationVecToDCM(rotation_vec3);

            testCase.verifyEqual(actDCM3, expDCM3, AbsTol=sonic.Tolerances.CompZero);

        end

        function axAngToDCMTest(testCase)
            %%% create a random rotation angle and axis
            expAngle = 27*pi/180;
            expAx = [1/sqrt(2); 0; -1/sqrt(2)];
            rotation_vec1 = expAngle*expAx;

            expDCM = expm(-sonic.Math.crossmat(rotation_vec1));
            actDCM = sonic.Attitude.axAngToDCM(expAx, expAngle);

            testCase.verifyEqual(actDCM, expDCM, AbsTol=sonic.Tolerances.CompZero);

            %%% what if I add 2*pi to the angle with same axis? It should be the same.
            expAngle2 = 27*pi/180 + 2*pi;
            actDCM2 = sonic.Attitude.axAngToDCM(expAx, expAngle2);

            testCase.verifyEqual(actDCM2, expDCM, AbsTol=sonic.Tolerances.CompZero);

            %%% create a zero vector. It should be identity.
            expAngle3  = 0;

            expDCM3 = eye(3);
            actDCM3 = sonic.Attitude.axAngToDCM(expAx, expAngle3);

            testCase.verifyEqual(actDCM3, expDCM3, AbsTol=sonic.Tolerances.CompZero);

        end

        function rotationVecToAxAngTest(testCase)
            %%% create a random rotation angle and axis
            expAngle1 = 27*pi/180;
            expAx1 = [1/sqrt(2); 0; -1/sqrt(2)];
            rotation_vec1 = expAngle1*expAx1;
            [actAx1, actAngle1] = sonic.Attitude.rotationVecToAxAng(rotation_vec1);

            testCase.verifyEqual(expAngle1, actAngle1, AbsTol=sonic.Tolerances.SmallAngle);
            testCase.verifyEqual(expAx1, actAx1, AbsTol=sonic.Tolerances.CompZero);

            %%% create a zero vector
            expAngle2 = 0*pi/180;
            Ax2 = [1/sqrt(2); 0; -1/sqrt(2)];
            rotation_vec2 = expAngle2*Ax2;
            [actAx2, actAngle2] = sonic.Attitude.rotationVecToAxAng(rotation_vec2);

            % this time, expect zero axes because angle is zero
            expAx2 = [0; 0; 0];
            testCase.verifyEqual(expAngle2, actAngle2, AbsTol=sonic.Tolerances.SmallAngle);
            testCase.verifyEqual(expAx2, actAx2, AbsTol=sonic.Tolerances.CompZero);
        end


        function interpolateTest(testCase)
            % Note: verifying with rotation vector here to test the healthy
            % behaviour of the quaternion SLERP

            % Random rotation
            ang1 = pi/3;
            ax1  = [1; 2; 3];
            ax1  = ax1 / norm(ax1);
            
            att1 = sonic.Attitude(ang1*ax1);
            
            % Other rotation on the same axis but with other angle
            ang2 = 4/3*pi;
            att2 = sonic.Attitude(ang2*ax1);
            
            %%% Propagate zero time the arclength
            att1bis = sonic.Attitude.interpolate(att1, att2, 0);
            testCase.verifyEqual(att1bis.dcm, att1.dcm, AbsTol=sonic.Tolerances.CompZero);

            %%% Propagate the entire arclength
            att2bis = sonic.Attitude.interpolate(att1, att2, 1);
            testCase.verifyEqual(att2bis.dcm, att2.dcm, AbsTol=sonic.Tolerances.CompZero);
            
            %%% Propagate only a fraction of the arclength
            t = 0.4;
            attbetween = sonic.Attitude.interpolate(att1, att2, t);
            
            % expect q SLERP to give linear interpolation of rot_vec
            expAng = (1-t)*ang1 + t*ang2;
            actAng = norm(attbetween.rotation_vec);

            testCase.verifyEqual(actAng, expAng, AbsTol=sonic.Tolerances.CompZero);

            %%% Propagate more than the arclength
            t = 1.4;
            attbetween = sonic.Attitude.interpolate(att1, att2, t);
            
            % expect q SLERP to give linear interpolation of rot_vec
            expAng = (1-t)*ang1 + t*ang2;
            actAng = norm(attbetween.rotation_vec);

            testCase.verifyEqual(actAng, expAng, AbsTol=sonic.Tolerances.CompZero);
        end

        function rotationVecToQuatTest(testCase)
            % just making sure that this works at small angle too
            ang = sonic.Tolerances.SmallAngle * 0.9;
            ax = [1; 2; 3];
            ax = ax/norm(ax);
            
            expq = [ax*sin(ang/2); cos(ang/2)];
            att = sonic.Attitude(ax*ang);
            actq = [att.quat_v; att.quat_s];

            testCase.verifyEqual(actq, expq, AbsTol=sonic.Tolerances.CompZero)
        end

        function rotateTest(testCase)

            %%% test with non-infinity points
            homogeneousPoints = [0 1;0 2;0 3; 1 1];
            att = sonic.Attitude([1; 2; 3]);
            T = [att.dcm, zeros(3,1); zeros(1,3), 1];

            points = sonic.Points3(homogeneousPoints);
            rotatedPoints = att.rotate(points);
            
            expRotatedPoints = T * homogeneousPoints;
            testCase.verifyEqual(rotatedPoints.p3, expRotatedPoints, AbsTol=sonic.Tolerances.CompZero)

            %%% test with points at infinity
            homogeneousPoints = [2; 3; 1; 0];
            att = sonic.Attitude([1; 2; 3]);
            T = [att.dcm, zeros(3,1); zeros(1,3), 1];

            points = sonic.Points3(homogeneousPoints);
            rotatedPoints = att.rotate(points);
            
            expRotatedPoints = T * homogeneousPoints;
            testCase.verifyEqual(rotatedPoints.p3, expRotatedPoints, AbsTol=sonic.Tolerances.CompZero)

            %%% test with a mix of points at, and not at, infinity
            homogeneousPoints = [0 1 4; 
                                 0 2 5;
                                 0 3 6; 
                                 1 1 0];
            att = sonic.Attitude([1; 2; 3]);
            T = [att.dcm, zeros(3,1); zeros(1,3), 1];

            points = sonic.Points3(homogeneousPoints);
            rotatedPoints = att.rotate(points);
            
            expRotatedPoints = T * homogeneousPoints;
            testCase.verifyEqual(rotatedPoints.p3, expRotatedPoints, AbsTol=sonic.Tolerances.CompZero)

            %%% Test with a wrong object
            testCase.verifyError(@() att.rotate(1), ...
                'sonic:Attitude:rotate:invalidType')
        end


        function mtimesTest(testCase)
            %%% Attitude times Attitude
            rotation_vec = [1; 2; 3];
            att1 = sonic.Attitude(rotation_vec);

            expAttComp = att1.comp(att1);
            actAttComp = att1 * att1;

            testCase.verifyEqual(expAttComp.dcm, actAttComp.dcm, AbsTol=sonic.Tolerances.CompZero)

            %%% Attitude times Points3
            p3 = [1 2 3 4; 
                  5 6 7 8;
                  1 2 3 4;
                  1 0 2 0];
            points = sonic.Points3(p3);

            expPoints3 = att1.rotate(points);
            actPoints3 = att1 * points;
           
            testCase.verifyEqual(expPoints3.p3, actPoints3.p3, AbsTol=sonic.Tolerances.CompZero)

            %%% Attitude times PointsS2
            ra_dec = [1, 2, 3, -3;
                      1, -1, 0.5, 0.5];
            points = sonic.PointsS2(ra_dec);

            expPointsS2 = att1.rotate(points);
            actPointsS2 = att1 * points;

            testCase.verifyEqual(expPointsS2.u, actPointsS2.u, AbsTol=sonic.Tolerances.CompZero)

            %%% If wrong object is entered
            testCase.verifyError(@() att1*1, ...
                'sonic:Attitude:mtimes:invalidInputType')
        end

    end



end