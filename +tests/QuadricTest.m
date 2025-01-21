classdef QuadricTest < matlab.unittest.TestCase

    methods (TestClassSetup)
        % Shared setup for the entire test class
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)
        % Test methods

        function QuadricHyperboloidInitTest(testCase)
            % example hyperboloid
            locus =  [13     3     9     6;
                3   -16     9     0;
                9     9   -14     6;
                6     0     6   -19];
            envelope = adjoint(locus);
            expLocus = sonic.Math.A_toDet1(locus);
            expEnvelope = sonic.Math.A_toDet1(envelope);

            quadric1 = sonic.Quadric(locus,"locus");

            quadric2 = sonic.Quadric(envelope,"envelope");

            actLocus1 = quadric1.locus;
            actEnvelope1 = quadric1.envelope;

            actLocus2 = quadric2.locus;
            actEnvelope2 = quadric2.envelope;

            testCase.verifyEqual(actLocus1,expLocus,AbsTol = sonic.Tolerances.CompZero);
            testCase.verifyEqual(actLocus2,expLocus,AbsTol = sonic.Tolerances.CompZero);
            testCase.verifyEqual(actEnvelope1,expEnvelope,AbsTol = sonic.Tolerances.CompZero);
            testCase.verifyEqual(actEnvelope2,expEnvelope,AbsTol = sonic.Tolerances.CompZero)
            testCase.verifyMatches(quadric1.type,"hyperboloid")
            testCase.verifyMatches(quadric2.type,"hyperboloid")
            testCase.verifyTrue(quadric1.proper);
            testCase.verifyTrue(quadric2.proper);
            testCase.verifyTrue(quadric1.stable);
            testCase.verifyTrue(quadric2.stable);

        end

        function QuadricUnstableHyperboloidInitTest(testCase)
            % example hyperboloid
            locus =  [13     3     9     0;
                3   -16     9     0;
                9     9   -14     0;
                0     0     0  -19];
            xc = -140;
            yc = -100;
            zc = -800;

            T = [1,0,0,-xc;
                0,1,0,-yc;
                0,0,1,-zc;
                0,0,0,1];

            locus = inv(T)'*locus*inv(T);

            envelope = adjoint(locus);

            [~,indL] = max(abs(locus(:)));
            [~,indE] = max(abs(envelope(:)));

            expLocus = locus/locus(indL);
            expEnvelope = envelope/envelope(indE);

            quadric1 = sonic.Quadric(locus,"locus");
            quadric2 = sonic.Quadric(envelope,"envelope");

            actLocus1 = quadric1.locus;
            actEnvelope1 = quadric1.envelope;
            [~,indLS] = max(abs(actLocus1(:)));
            actLocus1 = actLocus1/actLocus1(indLS);
            [~,indES] = max(abs(actEnvelope1(:)));
            actEnvelope1 = actEnvelope1/actEnvelope1(indES);

            actLocus2 = quadric2.locus;
            actEnvelope2 = quadric2.envelope;
            [~,indLS] = max(abs(actLocus2(:)));
            actLocus2 = actLocus2/actLocus2(indLS);
            [~,indES] = max(abs(actEnvelope2(:)));
            actEnvelope2 = actEnvelope2/actEnvelope2(indES);

            testCase.verifyEqual(actLocus1,expLocus,AbsTol = sonic.Tolerances.RelaxedCompZero);
            testCase.verifyEqual(actLocus2,expLocus,AbsTol = sonic.Tolerances.RelaxedCompZero);
            testCase.verifyEqual(actEnvelope1,expEnvelope,AbsTol = sonic.Tolerances.RelaxedCompZero);
            testCase.verifyEqual(actEnvelope2,expEnvelope,AbsTol = sonic.Tolerances.RelaxedCompZero)
            testCase.verifyMatches(quadric1.type,"hyperboloid")
            testCase.verifyMatches(quadric2.type,"hyperboloid")
            testCase.verifyTrue(quadric1.proper);
            testCase.verifyTrue(quadric2.proper);
            testCase.verifyFalse(quadric1.stable);
            testCase.verifyFalse(quadric2.stable);

        end

        function QuadricEllipsoidInitTest(testCase)
            % example hyperboloid
            locus =  [13     3     9     1;
                3   16     9     2;
                9     9   14     3;
                1     2     3  15];
            envelope = adjoint(locus);
            expLocus = sonic.Math.A_toDet1(locus);
            expEnvelope = sonic.Math.A_toDet1(envelope);

            quadric1 = sonic.Quadric(locus,"locus");

            quadric2 = sonic.Quadric(envelope,"envelope");

            actLocus1 = quadric1.locus;
            actEnvelope1 = quadric1.envelope;

            actLocus2 = quadric2.locus;
            actEnvelope2 = quadric2.envelope;

            testCase.verifyEqual(actLocus1,expLocus,AbsTol = sonic.Tolerances.CompZero);
            testCase.verifyEqual(actLocus2,expLocus,AbsTol = sonic.Tolerances.CompZero);
            testCase.verifyEqual(actEnvelope1,expEnvelope,AbsTol = sonic.Tolerances.CompZero);
            testCase.verifyEqual(actEnvelope2,expEnvelope,AbsTol = sonic.Tolerances.CompZero)
            testCase.verifyMatches(quadric1.type,"ellipsoid")
            testCase.verifyMatches(quadric2.type,"ellipsoid")
            testCase.verifyTrue(quadric1.proper);
            testCase.verifyTrue(quadric2.proper);
            testCase.verifyTrue(quadric1.stable);
            testCase.verifyTrue(quadric2.stable);

        end

        function QuadricUnstableEllipsoidInitTest(testCase)
            % example hyperboloid
            locus =  [13     3     9     0;
                3   16     9     0;
                9     9   14     0;
                0     0     0  15];
            xc = -100;
            yc = -400;
            zc = -500;

            T = [1,0,0,-xc;
                0,1,0,-yc;
                0,0,1,-zc;
                0,0,0,1];

            locus = inv(T)'*locus*inv(T);

            envelope = adjoint(locus);

            [~,indL] = max(abs(locus(:)));
            [~,indE] = max(abs(envelope(:)));

            expLocus = locus/locus(indL);
            expEnvelope = envelope/envelope(indE);

            quadric1 = sonic.Quadric(locus,"locus");
            quadric2 = sonic.Quadric(envelope,"envelope");

            actLocus1 = quadric1.locus;
            actEnvelope1 = quadric1.envelope;
            [~,indLS] = max(abs(actLocus1(:)));
            actLocus1 = actLocus1/actLocus1(indLS);
            [~,indES] = max(abs(actEnvelope1(:)));
            actEnvelope1 = actEnvelope1/actEnvelope1(indES);

            actLocus2 = quadric2.locus;
            actEnvelope2 = quadric2.envelope;
            [~,indLS] = max(abs(actLocus2(:)));
            actLocus2 = actLocus2/actLocus2(indLS);
            [~,indES] = max(abs(actEnvelope2(:)));
            actEnvelope2 = actEnvelope2/actEnvelope2(indES);

            testCase.verifyEqual(actLocus1,expLocus,AbsTol = sonic.Tolerances.RelaxedCompZero);
            testCase.verifyEqual(actLocus2,expLocus,AbsTol = sonic.Tolerances.RelaxedCompZero);
            testCase.verifyEqual(actEnvelope1,expEnvelope,AbsTol = sonic.Tolerances.RelaxedCompZero);
            testCase.verifyEqual(actEnvelope2,expEnvelope,AbsTol = sonic.Tolerances.RelaxedCompZero)
            testCase.verifyMatches(quadric1.type,"ellipsoid")
            testCase.verifyMatches(quadric2.type,"ellipsoid")
            testCase.verifyTrue(quadric1.proper);
            testCase.verifyTrue(quadric2.proper);
            testCase.verifyFalse(quadric1.stable);
            testCase.verifyFalse(quadric2.stable);

        end

        function QuadricDiskQuadricTest(testCase)
            envelope = [0,0,2,1;
                0,2,0,-2;
                2,0,0,3;
                1,-2,3,5];

            locus = adjoint(envelope);

            [~,indL] = max(abs(locus(:)));
            [~,indE] = max(abs(envelope(:)));

            expLocus = locus/locus(indL);
            expEnvelope = envelope/envelope(indE);

            quadric = sonic.Quadric(envelope,"envelope");
            actLocus = quadric.locus;
            actEnvelope = quadric.envelope;
            [~,indLS] = max(abs(actLocus(:)));
            [~,indES] = max(abs(actEnvelope(:)));

            actLocus = actLocus / actLocus(indLS);
            actEnvelope = actEnvelope / actEnvelope(indES);

            testCase.verifyEqual(actLocus,expLocus,AbsTol = sonic.Tolerances.CompZero);
            testCase.verifyEqual(actEnvelope,expEnvelope,AbsTol = sonic.Tolerances.CompZero);
            testCase.verifyMatches(quadric.type,"diskQuadric");
            testCase.verifyFalse(quadric.proper);
            testCase.verifyTrue(quadric.stable);



        end

        function QuadricGetTypeTest(testCase)

            %disk quadric
            envelope = [0,0,2,1;
                0,2,0,-2;
                2,0,0,3;
                1,-2,3,5];

            % wrong matrix type string
            testCase.verifyError(@()sonic.Quadric.getType(envelope,"envElope"),'sonic:Quadric:getType:invalidInput')

            envelope = [0,0,2,6;
                0,2,0,-2;
                0,-2,0,2;
                0,0,2,6];

            % rank 2 envelopes not accepted
            testCase.verifyError(@()sonic.Quadric.getType(envelope,"envelope"),'sonic:Quadric:getType:invalidInput')

        end

        function QuadricCenterTest(testCase)

            % example hyperboloid
            initialLocus =  [13     3     9     0;
                3   16     9     0;
                9     9   14     0;
                0     0     0  15];

            xc = -100;
            yc = -400;
            zc = -500;

            T = [1,0,0,-xc;
                0,1,0,-yc;
                0,0,1,-zc;
                0,0,0,1];

            locus = inv(T)'*initialLocus*inv(T);
            expLocus = sonic.Math.A_toDet1(initialLocus);

            quadric = sonic.Quadric(locus,"locus");
            translatedQuadric = quadric.center();
            actLocus= translatedQuadric.locus;

            testCase.verifyEqual(actLocus,expLocus,AbsTol = sonic.Tolerances.CompZero);

            % construct a disk quadric with center at the origin
            xc = 0; yc = 0; a = 5; b = 3; psi = 50*pi/180; 

            A = a^2*sin(psi)^2 + b^2*cos(psi)^2;
            B = 2*(b^2 - a^2)*cos(psi)*sin(psi);
            C = a^2*cos(psi)^2 + b^2*sin(psi)^2;
            D = -2*A*xc - B*yc;
            E = -B*xc - 2*C*yc;
            F = A*xc^2 + B*xc*yc + C*yc^2 - a^2*b^2;
            
            % construct the conic locus
            locusUnnorm = [A, B/2, D/2; B/2, C, E/2; D/2, E/2, F];

            % construct the conic envelope
            C_star = adjoint(locusUnnorm);

            % vector from the origin to the center of the disk quadric
            center = [5;10;20]; k = [0;0;1];

            ui = center/norm(center);
            ei = cross(k,ui)/norm(cross(k,ui));
            ni = cross(ui,ei)/norm(cross(ui,ei));
            H_M = [ei,ni, center];
            diskQuadricEnvelope = [H_M*C_star*H_M',  H_M*C_star*k;
                k'*C_star*H_M',    k'*C_star*k];

            originalDiskQuadric = sonic.Quadric(diskQuadricEnvelope,"envelope");

            center = [center;1];
            P = [center(4), 0,0, -center(1);
                    0, center(4),0, -center(2);
                    0, 0,center(4), -center(3);
                    0,0,0,center(4)];

            expDiskQuadric = P*diskQuadricEnvelope*P';

            finalDiskQuadric = originalDiskQuadric.center();
            actDiskQuadric = finalDiskQuadric.envelope;

            actDiskQuadric = actDiskQuadric/actDiskQuadric(4,4);
            expDiskQuadric = expDiskQuadric/expDiskQuadric(4,4);

            testCase.verifyEqual(actDiskQuadric,expDiskQuadric,AbsTol = sonic.Tolerances.CompZero);


        end
    end

end