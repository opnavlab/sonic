classdef ConicTest < matlab.unittest.TestCase

    methods (TestClassSetup)
        % Shared setup for the entire test class
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)

        % Test methods

        function ConicIntersectingLinesTest(testCase)
            % test creation of conic object for degenerate conic (intersecting lines)

            % define the locus
            locus = [-1,0,1;
                0,1,-2;
                1,-2,3];

            % define implicit
            implicit = [-1,0,1,2,-4,3]';

            % define envelope
            envelope =  [1, 2,1;
                2,4, 2;
                1, 2, 1];

            % define P2 coordinates of the two lines
            L = [1,-1;
                1,1;
                -3,-1];

            % normalize expected P2 lines for comparison
            expDegenerate = L./L(3,:);

            % expected locus and envelope are the same as the original ones
            expLocus = locus;
            expEnvelope = envelope;
            expImplicit = implicit/norm(implicit);

            % create conic object using input locus, degenerate and
            % implicit
            conicLocus = sonic.Conic(locus,"locus");
            conicDegenerate = sonic.Conic(sonic.Lines2(L));
            conicImplicit = sonic.Conic(implicit);

            % extract locus
            actLocusLocus = conicLocus.locus;
            actEnvelopeLocus = conicLocus.envelope;
            actImplicitLocus = conicLocus.implicit;
            actExplicitLocus = conicLocus.explicit;
            actDegenerateLocus = conicLocus.degenerateLocus.p2;
            actTypeLocus = conicLocus.type;
            actProperLocus = conicLocus.proper;
            actStableLocus = conicLocus.stable;

            % normalize locus, envelope, lines and implicit
            actEnvelopeLocus = actEnvelopeLocus/actEnvelopeLocus(1,1);
            actLocusLocus = actLocusLocus/actLocusLocus(1,3);
            actDegenerateLocus = actDegenerateLocus./actDegenerateLocus(3,:);

            if actDegenerateLocus(1)>0
                actDegenerateLocus = [actDegenerateLocus(:,2),actDegenerateLocus(:,1)];
            end

            if actImplicitLocus(1)>0
                actImplicitLocus = -actImplicitLocus;
            end

            % extract degenerate
            actLocusDegenerate = conicDegenerate.locus;
            actImplicitDegenerate = conicDegenerate.implicit;
            actEnvelopeDegenerate = conicDegenerate.envelope;
            actExplicitDegenerate = conicDegenerate.explicit;
            actDegenerateDegenerate = conicDegenerate.degenerateLocus.p2;
            actTypeDegenerate = conicDegenerate.type;
            actProperDegenerate = conicDegenerate.proper;
            actStableDegenerate = conicDegenerate.stable;

            % normalize locus, envelope, lines and implicit
            actEnvelopeDegenerate = actEnvelopeDegenerate/actEnvelopeDegenerate(1,1);
            actLocusDegenerate = actLocusDegenerate/actLocusDegenerate(1,3);
            actDegenerateDegenerate = actDegenerateDegenerate./actDegenerateDegenerate(3,:);

            if actDegenerateDegenerate(1,1)>0
                actDegenerateDegenerate = [actDegenerateDegenerate(:,2),actDegenerateDegenerate(:,1)];
            end

            if actImplicitDegenerate(1)>0
                actImplicitDegenerate = -actImplicitDegenerate;
            end


            % extract implicit
            actLocusImplicit = conicImplicit.locus;
            actImplicitImplicit = conicImplicit.implicit;
            actEnvelopeImplicit = conicImplicit.envelope;
            actExplicitImplicit = conicImplicit.explicit;
            actDegenerateImplicit = conicImplicit.degenerateLocus.p2;
            actTypeImplicit = conicImplicit.type;
            actProperImplicit = conicImplicit.proper;
            actStableImplicit = conicImplicit.stable;

            % normalize locus, envelope, lines and implicit
            actEnvelopeImplicit = actEnvelopeImplicit/actEnvelopeImplicit(1,1);
            actLocusImplicit = actLocusImplicit/actLocusImplicit(1,3);
            actDegenerateImplicit = actDegenerateImplicit./actDegenerateImplicit(3,:);

            if actDegenerateImplicit(1,1)>0
                actDegenerateImplicit = [actDegenerateImplicit(:,2),actDegenerateImplicit(:,1)];
            end

            if actImplicitImplicit(1)>0
                actImplicitImplicit = -actImplicitImplicit;
            end

            % test locus
            testCase.verifyEqual(actLocusLocus,expLocus,AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actImplicitLocus,expImplicit,AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actEnvelopeLocus,expEnvelope,AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actExplicitLocus,zeros(5,1),AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actDegenerateLocus, expDegenerate, AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyTrue(actStableLocus);
            testCase.verifyMatches(actTypeLocus,"intersectingLines");
            testCase.verifyFalse(actProperLocus);

            % test degenerate
            testCase.verifyEqual(actLocusDegenerate,expLocus,AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actImplicitDegenerate,expImplicit,AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actEnvelopeDegenerate,expEnvelope,AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actExplicitDegenerate,zeros(5,1),AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actDegenerateDegenerate,expDegenerate, AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyTrue(actStableDegenerate);
            testCase.verifyMatches(actTypeDegenerate,"intersectingLines");
            testCase.verifyFalse(actProperDegenerate)


            % test implicit
            testCase.verifyEqual(actLocusImplicit,expLocus,AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actImplicitImplicit,expImplicit,AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actEnvelopeImplicit,expEnvelope,AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actExplicitImplicit,zeros(5,1),AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actDegenerateImplicit,expDegenerate, AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyTrue(actStableImplicit);
            testCase.verifyMatches(actTypeImplicit,"intersectingLines");
            testCase.verifyFalse(actProperImplicit);
        end

        function ConicParallelLinesTest(testCase)
            % test creation of conic object for degenerate conic (parallel lines)

            locus = [3,6,-7;
                6,12,-14;
                -7,-14,-5];
            implicit = [3,12,12,-14,-28,-5]';
            implicit = implicit/norm(implicit);

            envelope = [-256, 128,0;
                128,-64  0;
                0,0,0];

            L = [3,1;6,2;1,-5];
            expDegenerate = L./L(3,:);
            expLocus = locus/locus(1,1);
            expEnvelope = envelope/envelope(1,1);
            expImplicit = implicit;

            conicLocus = sonic.Conic(locus,"locus");
            conicDegenerate = sonic.Conic(sonic.Lines2(L));
            conicImplicit = sonic.Conic(implicit);

            % extract locus
            actLocusLocus = conicLocus.locus;
            actEnvelopeLocus = conicLocus.envelope;
            actImplicitLocus = conicLocus.implicit;
            actExplicitLocus = conicLocus.explicit;
            actDegenerateLocus = conicLocus.degenerateLocus.p2;
            actTypeLocus = conicLocus.type;
            actProperLocus = conicLocus.proper;
            actStableLocus = conicLocus.stable;

            actEnvelopeLocus = actEnvelopeLocus/actEnvelopeLocus(1,1);
            actLocusLocus = actLocusLocus/actLocusLocus(1,1);
            actDegenerateLocus = actDegenerateLocus./actDegenerateLocus(3,:);

            if actDegenerateLocus(1)<0
                actDegenerateLocus = [actDegenerateLocus(:,2),actDegenerateLocus(:,1)];
            end

            if actImplicitLocus(1)<0
                actImplicitLocus = -actImplicitLocus;
            end

            % extract degenerate
            actLocusDegenerate = conicDegenerate.locus;
            actImplicitDegenerate = conicDegenerate.implicit;
            actEnvelopeDegenerate = conicDegenerate.envelope;
            actExplicitDegenerate = conicDegenerate.explicit;
            actDegenerateDegenerate = conicDegenerate.degenerateLocus.p2;
            actTypeDegenerate = conicDegenerate.type;
            actProperDegenerate = conicDegenerate.proper;
            actStableDegenerate = conicDegenerate.stable;

            actEnvelopeDegenerate = actEnvelopeDegenerate/actEnvelopeDegenerate(1,1);
            actLocusDegenerate = actLocusDegenerate/actLocusDegenerate(1,1);
            actDegenerateDegenerate = actDegenerateDegenerate./actDegenerateDegenerate(3,:);

            if actDegenerateDegenerate(1,1)<0
                actDegenerateDegenerate = [actDegenerateDegenerate(:,2),actDegenerateDegenerate(:,1)];
            end

            if actImplicitDegenerate(1)<0
                actImplicitDegenerate = -actImplicitDegenerate;
            end


            % extract implicit
            actLocusImplicit = conicImplicit.locus;
            actImplicitImplicit = conicImplicit.implicit;
            actEnvelopeImplicit = conicImplicit.envelope;
            actExplicitImplicit = conicImplicit.explicit;
            actDegenerateImplicit = conicImplicit.degenerateLocus.p2;
            actTypeImplicit = conicImplicit.type;
            actProperImplicit = conicImplicit.proper;
            actStableImplicit = conicImplicit.stable;

            actEnvelopeImplicit = actEnvelopeImplicit/actEnvelopeImplicit(1,1);
            actLocusImplicit = actLocusImplicit/actLocusImplicit(1,1);
            actDegenerateImplicit = actDegenerateImplicit./actDegenerateImplicit(3,:);


            if actDegenerateImplicit(1,1)<0
                actDegenerateImplicit = [actDegenerateImplicit(:,2),actDegenerateImplicit(:,1)];
            end

            if actImplicitImplicit(1)<0
                actImplicitImplicit = -actImplicitImplicit;
            end

            % test locus
            testCase.verifyEqual(actLocusLocus,expLocus,AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actImplicitLocus,expImplicit,AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actEnvelopeLocus,expEnvelope,AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actExplicitLocus,zeros(5,1),AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actDegenerateLocus, expDegenerate, AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyTrue(actStableLocus);
            testCase.verifyMatches(actTypeLocus,"parallelLines");
            testCase.verifyFalse(actProperLocus);

            % test degenerate
            testCase.verifyEqual(actLocusDegenerate,expLocus,AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actImplicitDegenerate,expImplicit,AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actEnvelopeDegenerate,expEnvelope,AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actExplicitDegenerate,zeros(5,1),AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actDegenerateDegenerate,expDegenerate, AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyTrue(actStableDegenerate);
            testCase.verifyMatches(actTypeDegenerate,"parallelLines");
            testCase.verifyFalse(actProperDegenerate)


            % test implicit
            testCase.verifyEqual(actLocusImplicit,expLocus,AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actImplicitImplicit,expImplicit,AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actEnvelopeImplicit,expEnvelope,AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actExplicitImplicit,zeros(5,1),AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actDegenerateImplicit,expDegenerate, AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyTrue(actStableImplicit);
            testCase.verifyMatches(actTypeImplicit,"parallelLines");
            testCase.verifyFalse(actProperImplicit);

        end

        function ConicCoincidingLinesTest(testCase)
            % test creation of conic object for degenerate conic (doulbe line)

            locus = [9, 36/2,3;
                36/2,36,6;
                3,6,1];

            implicit = [9, 36, 36, 6, 12, 1]';
            implicit = implicit/norm(implicit);

            L = [3,3;6,6;1,1];
            expDegenerate = L./L(3,:);

            expEnvelope = zeros(3,3);
            expImplicit = implicit;
            expLocus = locus;

            conicLocus = sonic.Conic(locus,"locus");
            conicDegenerate = sonic.Conic(sonic.Lines2(L));
            conicImplicit = sonic.Conic(implicit);

            % extract locus
            actLocusLocus = conicLocus.locus;
            actEnvelopeLocus = conicLocus.envelope;
            actImplicitLocus = conicLocus.implicit;
            actExplicitLocus = conicLocus.explicit;
            actDegenerateLocus = conicLocus.degenerateLocus.p2;
            actTypeLocus = conicLocus.type;
            actProperLocus = conicLocus.proper;
            actStableLocus = conicLocus.stable;

            actLocusLocus = actLocusLocus/actLocusLocus(3,3);
            actDegenerateLocus = actDegenerateLocus./actDegenerateLocus(3,:);

            if actDegenerateLocus(1)<0
                actDegenerateLocus = [actDegenerateLocus(:,2),actDegenerateLocus(:,1)];
            end

            if actImplicitLocus(1)<0
                actImplicitLocus = -actImplicitLocus;
            end

            % extract degenerate
            actLocusDegenerate = conicDegenerate.locus;
            actImplicitDegenerate = conicDegenerate.implicit;
            actEnvelopeDegenerate = conicDegenerate.envelope;
            actExplicitDegenerate = conicDegenerate.explicit;
            actDegenerateDegenerate = conicDegenerate.degenerateLocus.p2;
            actTypeDegenerate = conicDegenerate.type;
            actProperDegenerate = conicDegenerate.proper;
            actStableDegenerate = conicDegenerate.stable;

            actLocusDegenerate = actLocusDegenerate/actLocusDegenerate(3,3);
            actDegenerateDegenerate = actDegenerateDegenerate./actDegenerateDegenerate(3,:);

            if actDegenerateDegenerate(1,1)<0
                actDegenerateDegenerate = [actDegenerateDegenerate(:,2),actDegenerateDegenerate(:,1)];
            end

            if actImplicitDegenerate(1)<0
                actImplicitDegenerate = -actImplicitDegenerate;
            end


            % extract implicit
            actLocusImplicit = conicImplicit.locus;
            actImplicitImplicit = conicImplicit.implicit;
            actEnvelopeImplicit = conicImplicit.envelope;
            actExplicitImplicit = conicImplicit.explicit;
            actDegenerateImplicit = conicImplicit.degenerateLocus.p2;
            actTypeImplicit = conicImplicit.type;
            actProperImplicit = conicImplicit.proper;
            actStableImplicit = conicImplicit.stable;

            actLocusImplicit = actLocusImplicit/actLocusImplicit(3,3);
            actDegenerateImplicit = actDegenerateImplicit./actDegenerateImplicit(3,:);


            if actDegenerateImplicit(1,1)<0
                actDegenerateImplicit = [actDegenerateImplicit(:,2),actDegenerateImplicit(:,1)];
            end

            if actImplicitImplicit(1)<0
                actImplicitImplicit = -actImplicitImplicit;
            end

            % test locus
            testCase.verifyEqual(actLocusLocus,expLocus,AbsTol=sonic.Tolerances.RelaxedCompZero);
            testCase.verifyEqual(actImplicitLocus,expImplicit,AbsTol=sonic.Tolerances.RelaxedCompZero);
            testCase.verifyEqual(actEnvelopeLocus,expEnvelope,AbsTol=sonic.Tolerances.RelaxedCompZero);
            testCase.verifyEqual(actExplicitLocus,zeros(5,1),AbsTol=sonic.Tolerances.RelaxedCompZero);
            testCase.verifyEqual(actDegenerateLocus, expDegenerate, AbsTol=sonic.Tolerances.RelaxedCompZero);
            testCase.verifyTrue(actStableLocus);
            testCase.verifyMatches(actTypeLocus,"parallelLines");
            testCase.verifyFalse(actProperLocus);

            % test degenerate
            testCase.verifyEqual(actLocusDegenerate,expLocus,AbsTol=sonic.Tolerances.RelaxedCompZero);
            testCase.verifyEqual(actImplicitDegenerate,expImplicit,AbsTol=sonic.Tolerances.RelaxedCompZero);
            testCase.verifyEqual(actEnvelopeDegenerate,expEnvelope,AbsTol=sonic.Tolerances.RelaxedCompZero);
            testCase.verifyEqual(actExplicitDegenerate,zeros(5,1),AbsTol=sonic.Tolerances.RelaxedCompZero);
            testCase.verifyEqual(actDegenerateDegenerate,expDegenerate, AbsTol=sonic.Tolerances.RelaxedCompZero);
            testCase.verifyTrue(actStableDegenerate);
            testCase.verifyMatches(actTypeDegenerate,"parallelLines");
            testCase.verifyFalse(actProperDegenerate)


            % test implicit
            testCase.verifyEqual(actLocusImplicit,expLocus,AbsTol=sonic.Tolerances.RelaxedCompZero);
            testCase.verifyEqual(actImplicitImplicit,expImplicit,AbsTol=sonic.Tolerances.RelaxedCompZero);
            testCase.verifyEqual(actEnvelopeImplicit,expEnvelope,AbsTol=sonic.Tolerances.RelaxedCompZero);
            testCase.verifyEqual(actExplicitImplicit,zeros(5,1),AbsTol=sonic.Tolerances.RelaxedCompZero);
            testCase.verifyEqual(actDegenerateImplicit,expDegenerate, AbsTol=sonic.Tolerances.RelaxedCompZero);
            testCase.verifyTrue(actStableImplicit);
            testCase.verifyMatches(actTypeImplicit,"parallelLines");
            testCase.verifyFalse(actProperImplicit);

        end

        function ConicPointTest(testCase)
            % test creation of conic object for degenerate conic (point)

            locus = [13,13,-30;
                13,26,-50;
                -30,-50,100];
            implicit =  [13,26,26,-60,-100,100]';
            L = [3+2*1i,3-2*1i;
                5-1i,5+1i;
                -10,-10];
            envelope= [1.0000,2.0000, 1.3000;
                2.0000,4.0000,2.6000;
                1.3000,2.6000,1.6900];

            implicit = implicit/norm(implicit);
            expDegenerate = L./L(3,:);
            expLocus = locus/locus(1,1);
            expEnvelope = envelope;
            expImplicit = implicit;

            conicLocus = sonic.Conic(locus,"locus");
            conicDegenerate = sonic.Conic(sonic.Lines2(L));
            conicImplicit = sonic.Conic(implicit);

            % extract locus
            actLocusLocus = conicLocus.locus;
            actEnvelopeLocus = conicLocus.envelope;
            actImplicitLocus = conicLocus.implicit;
            actExplicitLocus = conicLocus.explicit;
            actDegenerateLocus = conicLocus.degenerateLocus.p2;
            actTypeLocus = conicLocus.type;
            actProperLocus = conicLocus.proper;
            actStableLocus = conicLocus.stable;

            actEnvelopeLocus = actEnvelopeLocus/actEnvelopeLocus(1,1);
            actLocusLocus = actLocusLocus/actLocusLocus(1,1);
            actDegenerateLocus = actDegenerateLocus./actDegenerateLocus(3,:);

            if imag(actDegenerateLocus(1,1))>0
                actDegenerateLocus = [actDegenerateLocus(:,2),actDegenerateLocus(:,1)];
            end

            if actImplicitLocus(1)<0
                actImplicitLocus = -actImplicitLocus;
            end

            % extract degenerate
            actLocusDegenerate = conicDegenerate.locus;
            actImplicitDegenerate = conicDegenerate.implicit;
            actEnvelopeDegenerate = conicDegenerate.envelope;
            actExplicitDegenerate = conicDegenerate.explicit;
            actDegenerateDegenerate = conicDegenerate.degenerateLocus.p2;
            actTypeDegenerate = conicDegenerate.type;
            actProperDegenerate = conicDegenerate.proper;
            actStableDegenerate = conicDegenerate.stable;

            actEnvelopeDegenerate = actEnvelopeDegenerate/actEnvelopeDegenerate(1,1);
            actLocusDegenerate = actLocusDegenerate/actLocusDegenerate(1,1);
            actDegenerateDegenerate = actDegenerateDegenerate./actDegenerateDegenerate(3,:);

            if imag(actDegenerateDegenerate(1,1))>0
                actDegenerateDegenerate = [actDegenerateDegenerate(:,2),actDegenerateDegenerate(:,1)];
            end


            if actImplicitDegenerate(1)<0
                actImplicitDegenerate = -actImplicitDegenerate;
            end

            % extract implicit
            actLocusImplicit = conicImplicit.locus;
            actImplicitImplicit = conicImplicit.implicit;
            actEnvelopeImplicit = conicImplicit.envelope;
            actExplicitImplicit = conicImplicit.explicit;
            actDegenerateImplicit = conicImplicit.degenerateLocus.p2;
            actTypeImplicit = conicImplicit.type;
            actProperImplicit = conicImplicit.proper;
            actStableImplicit = conicImplicit.stable;

            actEnvelopeImplicit = actEnvelopeImplicit/actEnvelopeImplicit(1,1);
            actLocusImplicit = actLocusImplicit/actLocusImplicit(1,1);
            actDegenerateImplicit = actDegenerateImplicit./actDegenerateImplicit(3,:);


            if imag(actDegenerateImplicit(1,1))>0
                actDegenerateImplicit = [actDegenerateImplicit(:,2),actDegenerateImplicit(:,1)];
            end

            if actImplicitImplicit(1)<0
                actImplicitImplicit = -actImplicitImplicit;
            end

            % test locus
            testCase.verifyEqual(actLocusLocus,expLocus,AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actImplicitLocus,expImplicit,AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actEnvelopeLocus,expEnvelope,AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actExplicitLocus,zeros(5,1),AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actDegenerateLocus, expDegenerate, AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyTrue(actStableLocus);
            testCase.verifyMatches(actTypeLocus,"point");
            testCase.verifyFalse(actProperLocus);

            % test degenerate
            testCase.verifyEqual(actLocusDegenerate,expLocus,AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actImplicitDegenerate,expImplicit,AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actEnvelopeDegenerate,expEnvelope,AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actExplicitDegenerate,zeros(5,1),AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actDegenerateDegenerate,expDegenerate, AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyTrue(actStableDegenerate);
            testCase.verifyMatches(actTypeDegenerate,"point");
            testCase.verifyFalse(actProperDegenerate)


            % test implicit
            testCase.verifyEqual(actLocusImplicit,expLocus,AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actImplicitImplicit,expImplicit,AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actEnvelopeImplicit,expEnvelope,AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actExplicitImplicit,zeros(5,1),AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(actDegenerateImplicit,expDegenerate, AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyTrue(actStableImplicit);
            testCase.verifyMatches(actTypeImplicit,"point");
            testCase.verifyFalse(actProperImplicit);
        end

        function ConicEllipseInitTest(testCase)
            % test creation of conic object for ellipse

            xc = 1;
            yc = 2;
            a = 4;
            b = 3;
            psi = pi/6;

            A = a^2*sin(psi)^2 + b^2*cos(psi)^2;
            B = 2*(b^2 - a^2)*cos(psi)*sin(psi);
            C = a^2*cos(psi)^2 + b^2*sin(psi)^2;
            D = -2*A*xc - B*yc;
            E = -B*xc - 2*C*yc;
            F = A*xc^2 + B*xc*yc + C*yc^2 - a^2*b^2;

            expType = "ellipse";
            expProper = true;
            expStable = true;

            explicit = [xc;yc;a;b;psi];
            implicit = [A,B,C,D,E,F]';
            locus = [A,B/2,D/2;
                B/2,C,E/2;
                D/2,E/2,F];
            envelope = adjoint(locus);
            expExplicit = [xc;yc;a;b;psi];


            tolerance = sonic.Tolerances.CompZero;


            expEnvelope = sonic.Math.A_toDet1(envelope);
            expLocus = sonic.Math.A_toDet1(locus);

            % normalize explicit and implicit for comparison
            expImplicit = implicit/norm(implicit);
            if sign(expImplicit(1))<0
                expImplicit = -expImplicit;
            end

            %% create conic from implicit
            conicImplicit = sonic.Conic(implicit);

            % extract properties
            actImplicitImplicit = conicImplicit.implicit;
            actEnvelopeImplicit = conicImplicit.envelope;
            actLocusImplicit = conicImplicit.locus;
            actExplicitImplicit = conicImplicit.explicit;
            actProperImplicit = conicImplicit.proper;
            actTypeImplicit = conicImplicit.type;
            actStableImplicit=conicImplicit.stable;

            if sign(actImplicitImplicit(1))<0
                actImplicitImplicit = -actImplicitImplicit;
            end

            if abs(psi-actExplicitImplicit(end))>sonic.Tolerances.CompZero
                actExplicitImplicit(end) = pi+actExplicitImplicit(end);
            end

            % test implicit
            testCase.verifyEqual(actImplicitImplicit,expImplicit,AbsTol=tolerance);
            testCase.verifyEqual(actEnvelopeImplicit,expEnvelope,AbsTol=tolerance);
            testCase.verifyEqual(actLocusImplicit,expLocus,AbsTol=tolerance);
            testCase.verifyEqual(actProperImplicit,expProper);
            testCase.verifyEqual(actStableImplicit,expStable);
            testCase.verifyMatches(actTypeImplicit,expType);
            testCase.verifyEqual(actExplicitImplicit,expExplicit,AbsTol=tolerance);

            %% create conic from explicit
            conicExplicit = sonic.Conic(explicit);

            % extract properties
            actImplicitExplicit = conicExplicit.implicit;
            actEnvelopeExplicit = conicExplicit.envelope;
            actLocusExplicit = conicExplicit.locus;
            actExplicitExplicit = conicExplicit.explicit;
            actProperExplicit = conicExplicit.proper;
            actTypeExplicit = conicExplicit.type;
            actStableExplicit=conicExplicit.stable;

            if sign(actImplicitExplicit(1))<0
                actImplicitExplicit = -actImplicitExplicit;
            end

            if abs(psi-actExplicitExplicit(end))>sonic.Tolerances.CompZero
                actExplicitExplicit(end) = pi+actExplicitExplicit(end);
            end

            % test explicit
            testCase.verifyEqual(actImplicitExplicit,expImplicit,AbsTol=tolerance);
            testCase.verifyEqual(actEnvelopeExplicit,expEnvelope,AbsTol=tolerance);
            testCase.verifyEqual(actLocusExplicit,expLocus,AbsTol=tolerance);
            testCase.verifyEqual(actExplicitExplicit,expExplicit,AbsTol=tolerance);
            testCase.verifyEqual(actProperExplicit,expProper);
            testCase.verifyEqual(actStableExplicit,expStable);
            testCase.verifyMatches(actTypeExplicit,expType);
            testCase.verifyEqual(actExplicitImplicit,expExplicit,AbsTol=tolerance);


            %% create conic from locus
            conicLocus = sonic.Conic(locus,"locus");

            actImplicitLocus= conicLocus.implicit;
            actEnvelopeLocus = conicLocus.envelope;
            actLocusLocus = conicLocus.locus;
            actExplicitLocus = conicLocus.explicit;
            actProperLocus = conicLocus.proper;
            actStableLocus=conicLocus.stable;
            actTypeLocus = conicLocus.type;

            % if necessary, normalize
            if sign(actImplicitLocus(1))<0
                actImplicitLocus = -actImplicitLocus;
            end

            if abs(psi-actExplicitLocus(end))>sonic.Tolerances.CompZero
                actExplicitLocus(end) = pi+actExplicitLocus(end);
            end

            % test locus
            testCase.verifyEqual(actImplicitLocus,expImplicit,AbsTol=tolerance);
            testCase.verifyEqual(actEnvelopeLocus,expEnvelope,AbsTol=tolerance);
            testCase.verifyEqual(actLocusLocus,expLocus,AbsTol=tolerance);
            testCase.verifyEqual(actProperLocus,expProper);
            testCase.verifyEqual(actStableLocus,expStable);
            testCase.verifyMatches(actTypeLocus,expType);
            testCase.verifyEqual(actExplicitLocus,expExplicit,AbsTol=tolerance);


            %% create conic from envelope
            conicEnvelope = sonic.Conic(envelope,"envelope");

            % extract properties
            actImplicitEnvelope = conicEnvelope.implicit;
            actLocusEnvelope = conicEnvelope.locus;
            actEnvelopeEnvelope = conicEnvelope.envelope;
            actExplicitEnvelope = conicEnvelope.explicit;
            actTypeEnvelope = conicEnvelope.type;
            actStableEnvelope=conicEnvelope.stable;
            actProperEnvelope = conicEnvelope.proper;

            % if necessary, normalize
            if sign(actImplicitEnvelope(1))<0
                actImplicitEnvelope = -actImplicitEnvelope;
            end

            if abs(psi-actExplicitEnvelope(end))>sonic.Tolerances.CompZero
                actExplicitEnvelope(end) = pi+actExplicitEnvelope(end);
            end

            % test envelope
            testCase.verifyEqual(actImplicitEnvelope,expImplicit,AbsTol=tolerance);
            testCase.verifyEqual(actEnvelopeEnvelope,expEnvelope,AbsTol=tolerance);
            testCase.verifyEqual(actLocusEnvelope,expLocus,AbsTol=tolerance);
            testCase.verifyEqual(actProperEnvelope,expProper);
            testCase.verifyEqual(actStableEnvelope,expStable);
            testCase.verifyMatches(actTypeEnvelope,expType);
            testCase.verifyEqual(actExplicitEnvelope,expExplicit,AbsTol=tolerance);


        end

        function ConicUnstableEllipseInitTest(testCase)
            % test creation of conic object for numerically unstable
            % ellipse

            xc = 750;
            yc = 500;
            a = 4;
            b = 3;
            psi = pi/6;

            A = a^2*sin(psi)^2 + b^2*cos(psi)^2;
            B = 2*(b^2 - a^2)*cos(psi)*sin(psi);
            C = a^2*cos(psi)^2 + b^2*sin(psi)^2;
            D = -2*A*xc - B*yc;
            E = -B*xc - 2*C*yc;
            F = A*xc^2 + B*xc*yc + C*yc^2 - a^2*b^2;

            expType = "ellipse";
            expProper = true;
            expStable = false;

            explicit = [xc;yc;a;b;psi];
            implicit = [A,B,C,D,E,F]';
            locus = [A,B/2,D/2;
                B/2,C,E/2;
                D/2,E/2,F];
            envelope = adjoint(locus);
            expExplicit = [xc;yc;a;b;psi];

            tolerance = sonic.Tolerances.RelaxedCompZero;
            [~,indL] = max(abs(locus(:)));
            expLocus= locus / locus(indL);
            [~,indE] = max(abs(envelope(:)));
            expEnvelope = envelope / envelope(indE);

            % normalize explicit and implicit for comparison
            expImplicit = implicit/norm(implicit);
            if sign(expImplicit(1))<0
                expImplicit = -expImplicit;
            end

            %% create conic from implicit
            conicImplicit = sonic.Conic(implicit);

            % extract properties
            actImplicitImplicit = conicImplicit.implicit;
            actEnvelopeImplicit = conicImplicit.envelope;
            actLocusImplicit = conicImplicit.locus;
            actExplicitImplicit = conicImplicit.explicit;
            actProperImplicit = conicImplicit.proper;
            actTypeImplicit = conicImplicit.type;
            actStableImplicit=conicImplicit.stable;

            [~,indE] = max(abs(actEnvelopeImplicit(:)));
            actEnvelopeImplicit = actEnvelopeImplicit / actEnvelopeImplicit(indE);
            [~,indL] = max(abs(actLocusImplicit(:)));
            actLocusImplicit = actLocusImplicit / actLocusImplicit(indL);

            if sign(actImplicitImplicit(1))<0
                actImplicitImplicit = -actImplicitImplicit;
            end

            if abs(psi-actExplicitImplicit(end))>sonic.Tolerances.CompZero
                actExplicitImplicit(end) = pi+actExplicitImplicit(end);
            end

            % test implicit
            testCase.verifyEqual(actImplicitImplicit,expImplicit,AbsTol=tolerance);
            testCase.verifyEqual(actEnvelopeImplicit,expEnvelope,AbsTol=tolerance);
            testCase.verifyEqual(actLocusImplicit,expLocus,AbsTol=tolerance);
            testCase.verifyEqual(actProperImplicit,expProper);
            testCase.verifyEqual(actStableImplicit,expStable);
            testCase.verifyMatches(actTypeImplicit,expType);
            testCase.verifyEqual(actExplicitImplicit,expExplicit,AbsTol=tolerance);

            %% create conic from explicit
            conicExplicit = sonic.Conic(explicit);

            % extract properties
            actImplicitExplicit = conicExplicit.implicit;
            actEnvelopeExplicit = conicExplicit.envelope;
            actLocusExplicit = conicExplicit.locus;
            actExplicitExplicit = conicExplicit.explicit;
            actProperExplicit = conicExplicit.proper;
            actTypeExplicit = conicExplicit.type;
            actStableExplicit=conicExplicit.stable;

            [~,indE] = max(abs(actEnvelopeExplicit(:)));
            actEnvelopeExplicit = actEnvelopeExplicit / actEnvelopeExplicit(indE);
            [~,indL] = max(abs(actLocusExplicit(:)));
            actLocusExplicit = actLocusExplicit / actLocusExplicit(indL);

            if sign(actImplicitExplicit(1))<0
                actImplicitExplicit = -actImplicitExplicit;
            end

            if abs(psi-actExplicitExplicit(end))>sonic.Tolerances.CompZero
                actExplicitExplicit(end) = pi+actExplicitExplicit(end);
            end

            % test explicit
            testCase.verifyEqual(actImplicitExplicit,expImplicit,AbsTol=tolerance);
            testCase.verifyEqual(actEnvelopeExplicit,expEnvelope,AbsTol=tolerance);
            testCase.verifyEqual(actLocusExplicit,expLocus,AbsTol=tolerance);
            testCase.verifyEqual(actExplicitExplicit,expExplicit,AbsTol=tolerance);
            testCase.verifyEqual(actProperExplicit,expProper);
            testCase.verifyEqual(actStableExplicit,expStable);
            testCase.verifyMatches(actTypeExplicit,expType);
            testCase.verifyEqual(actExplicitImplicit,expExplicit,AbsTol=tolerance);


            %% create conic from locus
            conicLocus = sonic.Conic(locus,"locus");

            actImplicitLocus= conicLocus.implicit;
            actEnvelopeLocus = conicLocus.envelope;
            actLocusLocus = conicLocus.locus;
            actExplicitLocus = conicLocus.explicit;
            actProperLocus = conicLocus.proper;
            actStableLocus=conicLocus.stable;
            actTypeLocus = conicLocus.type;


            [~,indE] = max(abs(actEnvelopeLocus(:)));
            actEnvelopeLocus = actEnvelopeLocus / actEnvelopeLocus(indE);
            [~,indL] = max(abs(actLocusLocus(:)));
            actLocusLocus = actLocusLocus / actLocusLocus(indL);


            % if necessary, normalize
            if sign(actImplicitLocus(1))<0
                actImplicitLocus = -actImplicitLocus;
            end

            if abs(psi-actExplicitLocus(end))>sonic.Tolerances.CompZero
                actExplicitLocus(end) = pi+actExplicitLocus(end);
            end

            % test locus
            testCase.verifyEqual(actImplicitLocus,expImplicit,AbsTol=tolerance);
            testCase.verifyEqual(actEnvelopeLocus,expEnvelope,AbsTol=tolerance);
            testCase.verifyEqual(actLocusLocus,expLocus,AbsTol=tolerance);
            testCase.verifyEqual(actProperLocus,expProper);
            testCase.verifyEqual(actStableLocus,expStable);
            testCase.verifyMatches(actTypeLocus,expType);
            testCase.verifyEqual(actExplicitLocus,expExplicit,AbsTol=tolerance);


            %% create conic from envelope
            conicEnvelope = sonic.Conic(envelope,"envelope");

            % extract properties
            actImplicitEnvelope = conicEnvelope.implicit;
            actLocusEnvelope = conicEnvelope.locus;
            actEnvelopeEnvelope = conicEnvelope.envelope;
            actExplicitEnvelope = conicEnvelope.explicit;
            actTypeEnvelope = conicEnvelope.type;
            actStableEnvelope=conicEnvelope.stable;
            actProperEnvelope = conicEnvelope.proper;

            [~,indE] = max(abs(actEnvelopeEnvelope(:)));
            actEnvelopeEnvelope = actEnvelopeEnvelope / actEnvelopeEnvelope(indE);
            [~,indL] = max(abs(actLocusEnvelope(:)));
            actLocusEnvelope = actLocusEnvelope / actLocusEnvelope(indL);


            % if necessary, normalize
            if sign(actImplicitEnvelope(1))<0
                actImplicitEnvelope = -actImplicitEnvelope;
            end

            if abs(psi-actExplicitEnvelope(end))>sonic.Tolerances.CompZero
                actExplicitEnvelope(end) = pi+actExplicitEnvelope(end);
            end

            % test envelope
            testCase.verifyEqual(actImplicitEnvelope,expImplicit,AbsTol=tolerance);
            testCase.verifyEqual(actEnvelopeEnvelope,expEnvelope,AbsTol=tolerance);
            testCase.verifyEqual(actLocusEnvelope,expLocus,AbsTol=tolerance);
            testCase.verifyEqual(actProperEnvelope,expProper);
            testCase.verifyEqual(actStableEnvelope,expStable);
            testCase.verifyMatches(actTypeEnvelope,expType);
            testCase.verifyEqual(actExplicitEnvelope,expExplicit,AbsTol=tolerance);
        end

        function ConicHyperbolaInitTest(testCase)
            % test creation of conic object for hyperbola

            xc = 5;
            yc = 10;
            a = -20;
            b = -8;
            psi = pi/6;

            A = b^2*cos(psi)^2 - a^2*sin(psi)^2;
            B= 2*(a^2+b^2)*cos(psi)*sin(psi);
            C = b^2*sin(psi)^2 - a^2*cos(psi)^2;
            D = -2*A*xc - B*yc;
            E = -B*xc - 2*C*yc;
            F = A*xc^2 + B*xc*yc + C*yc^2 - a^2*b^2;

            expProper = true;
            expStable = true;
            expType = "hyperbola";


            explicit = [xc;yc;a;b;psi];
            implicit = [A,B,C,D,E,F]';
            locus = [A,B/2,D/2;
                B/2,C,E/2;
                D/2,E/2,F];
            envelope = adjoint(locus);
            expExplicit = [xc;yc;a;b;psi];

            tolerance = sonic.Tolerances.CompZero;



            expEnvelope = sonic.Math.A_toDet1(envelope);
            expLocus = sonic.Math.A_toDet1(locus);

            % normalize explicit and implicit for comparison
            expImplicit = implicit/norm(implicit);
            if sign(expImplicit(1))<0
                expImplicit = -expImplicit;
            end

            %% create conic from implicit
            conicImplicit = sonic.Conic(implicit);

            % extract properties
            actImplicitImplicit = conicImplicit.implicit;
            actEnvelopeImplicit = conicImplicit.envelope;
            actLocusImplicit = conicImplicit.locus;
            actExplicitImplicit = conicImplicit.explicit;
            actProperImplicit = conicImplicit.proper;
            actTypeImplicit = conicImplicit.type;
            actStableImplicit=conicImplicit.stable;

            if sign(actImplicitImplicit(1))<0
                actImplicitImplicit = -actImplicitImplicit;
            end

            if abs(psi-actExplicitImplicit(end))>sonic.Tolerances.CompZero
                actExplicitImplicit(end) = pi+actExplicitImplicit(end);
            end

            % test implicit
            testCase.verifyEqual(actImplicitImplicit,expImplicit,AbsTol=tolerance);
            testCase.verifyEqual(actEnvelopeImplicit,expEnvelope,AbsTol=tolerance);
            testCase.verifyEqual(actLocusImplicit,expLocus,AbsTol=tolerance);
            testCase.verifyEqual(actProperImplicit,expProper);
            testCase.verifyEqual(actStableImplicit,expStable);
            testCase.verifyMatches(actTypeImplicit,expType);
            testCase.verifyEqual(actExplicitImplicit,expExplicit,AbsTol=tolerance);

            %% create conic from explicit
            conicExplicit = sonic.Conic(explicit);

            % extract properties
            actImplicitExplicit = conicExplicit.implicit;
            actEnvelopeExplicit = conicExplicit.envelope;
            actLocusExplicit = conicExplicit.locus;
            actExplicitExplicit = conicExplicit.explicit;
            actProperExplicit = conicExplicit.proper;
            actTypeExplicit = conicExplicit.type;
            actStableExplicit=conicExplicit.stable;

            if sign(actImplicitExplicit(1))<0
                actImplicitExplicit = -actImplicitExplicit;
            end

            if abs(psi-actExplicitExplicit(end))>sonic.Tolerances.CompZero
                actExplicitExplicit(end) = pi+actExplicitExplicit(end);
            end

            % test explicit
            testCase.verifyEqual(actImplicitExplicit,expImplicit,AbsTol=tolerance);
            testCase.verifyEqual(actEnvelopeExplicit,expEnvelope,AbsTol=tolerance);
            testCase.verifyEqual(actLocusExplicit,expLocus,AbsTol=tolerance);
            testCase.verifyEqual(actExplicitExplicit,expExplicit,AbsTol=tolerance);
            testCase.verifyEqual(actProperExplicit,expProper);
            testCase.verifyEqual(actStableExplicit,expStable);
            testCase.verifyMatches(actTypeExplicit,expType);
            testCase.verifyEqual(actExplicitImplicit,expExplicit,AbsTol=tolerance);


            %% create conic from locus
            conicLocus = sonic.Conic(locus,"locus");

            actImplicitLocus= conicLocus.implicit;
            actEnvelopeLocus = conicLocus.envelope;
            actLocusLocus = conicLocus.locus;
            actExplicitLocus = conicLocus.explicit;
            actProperLocus = conicLocus.proper;
            actStableLocus=conicLocus.stable;
            actTypeLocus = conicLocus.type;

            % if necessary, normalize
            if sign(actImplicitLocus(1))<0
                actImplicitLocus = -actImplicitLocus;
            end

            if abs(psi-actExplicitLocus(end))>sonic.Tolerances.CompZero
                actExplicitLocus(end) = pi+actExplicitLocus(end);
            end

            % test locus
            testCase.verifyEqual(actImplicitLocus,expImplicit,AbsTol=tolerance);
            testCase.verifyEqual(actEnvelopeLocus,expEnvelope,AbsTol=tolerance);
            testCase.verifyEqual(actLocusLocus,expLocus,AbsTol=tolerance);
            testCase.verifyEqual(actProperLocus,expProper);
            testCase.verifyEqual(actStableLocus,expStable);
            testCase.verifyMatches(actTypeLocus,expType);
            testCase.verifyEqual(actExplicitLocus,expExplicit,AbsTol=tolerance);


            %% create conic from envelope
            conicEnvelope = sonic.Conic(envelope,"envelope");

            % extract properties
            actImplicitEnvelope = conicEnvelope.implicit;
            actLocusEnvelope = conicEnvelope.locus;
            actEnvelopeEnvelope = conicEnvelope.envelope;
            actExplicitEnvelope = conicEnvelope.explicit;
            actTypeEnvelope = conicEnvelope.type;
            actStableEnvelope=conicEnvelope.stable;
            actProperEnvelope = conicEnvelope.proper;

            % if necessary, normalize
            if sign(actImplicitEnvelope(1))<0
                actImplicitEnvelope = -actImplicitEnvelope;
            end

            if abs(psi-actExplicitEnvelope(end))>sonic.Tolerances.CompZero
                actExplicitEnvelope(end) = pi+actExplicitEnvelope(end);
            end

            % test envelope
            testCase.verifyEqual(actImplicitEnvelope,expImplicit,AbsTol=tolerance);
            testCase.verifyEqual(actEnvelopeEnvelope,expEnvelope,AbsTol=tolerance);
            testCase.verifyEqual(actLocusEnvelope,expLocus,AbsTol=tolerance);
            testCase.verifyEqual(actProperEnvelope,expProper);
            testCase.verifyEqual(actStableEnvelope,expStable);
            testCase.verifyMatches(actTypeEnvelope,expType);
            testCase.verifyEqual(actExplicitEnvelope,expExplicit,AbsTol=tolerance);


        end

        function ConicUnstableHyperbolaInitTest(testCase)
            % test creation of conic object for numerically unstable hyperbola

            xc = 500;
            yc = 1000;
            a = -5;
            b = -3;
            psi = pi/6;

            A = b^2*cos(psi)^2 - a^2*sin(psi)^2;
            B= 2*(a^2+b^2)*cos(psi)*sin(psi);
            C = b^2*sin(psi)^2 - a^2*cos(psi)^2;
            D = -2*A*xc - B*yc;
            E = -B*xc - 2*C*yc;
            F = A*xc^2 + B*xc*yc + C*yc^2 - a^2*b^2;

            expProper = true;
            expStable = false;
            expType = "hyperbola";
            tolerance = sonic.Tolerances.RelaxedCompZero;

            explicit = [xc;yc;a;b;psi];
            implicit = [A,B,C,D,E,F]';
            locus = [A,B/2,D/2;
                B/2,C,E/2;
                D/2,E/2,F];
            envelope = adjoint(locus);
            expExplicit = [xc;yc;a;b;psi];


            [~,indL] = max(abs(locus(:)));
            expLocus= locus / locus(indL);
            [~,indE] = max(abs(envelope(:)));
            expEnvelope = envelope / envelope(indE);

            % normalize explicit and implicit for comparison
            expImplicit = implicit/norm(implicit);
            if sign(expImplicit(1))<0
                expImplicit = -expImplicit;
            end

            %% create conic from implicit
            conicImplicit = sonic.Conic(implicit);

            % extract properties
            actImplicitImplicit = conicImplicit.implicit;
            actEnvelopeImplicit = conicImplicit.envelope;
            actLocusImplicit = conicImplicit.locus;
            actExplicitImplicit = conicImplicit.explicit;
            actProperImplicit = conicImplicit.proper;
            actTypeImplicit = conicImplicit.type;
            actStableImplicit=conicImplicit.stable;

            [~,indE] = max(abs(actEnvelopeImplicit(:)));
            actEnvelopeImplicit = actEnvelopeImplicit / actEnvelopeImplicit(indE);
            [~,indL] = max(abs(actLocusImplicit(:)));
            actLocusImplicit = actLocusImplicit / actLocusImplicit(indL);

            if sign(actImplicitImplicit(1))<0
                actImplicitImplicit = -actImplicitImplicit;
            end

            if abs(psi-actExplicitImplicit(end))>sonic.Tolerances.CompZero
                actExplicitImplicit(end) = pi+actExplicitImplicit(end);
            end

            % test implicit
            testCase.verifyEqual(actImplicitImplicit,expImplicit,AbsTol=tolerance);
            testCase.verifyEqual(actEnvelopeImplicit,expEnvelope,AbsTol=tolerance);
            testCase.verifyEqual(actLocusImplicit,expLocus,AbsTol=tolerance);
            testCase.verifyEqual(actProperImplicit,expProper);
            testCase.verifyEqual(actStableImplicit,expStable);
            testCase.verifyMatches(actTypeImplicit,expType);
            testCase.verifyEqual(actExplicitImplicit,expExplicit,AbsTol=tolerance);

            %% create conic from explicit
            conicExplicit = sonic.Conic(explicit);

            % extract properties
            actImplicitExplicit = conicExplicit.implicit;
            actEnvelopeExplicit = conicExplicit.envelope;
            actLocusExplicit = conicExplicit.locus;
            actExplicitExplicit = conicExplicit.explicit;
            actProperExplicit = conicExplicit.proper;
            actTypeExplicit = conicExplicit.type;
            actStableExplicit=conicExplicit.stable;

            [~,indE] = max(abs(actEnvelopeExplicit(:)));
            actEnvelopeExplicit = actEnvelopeExplicit / actEnvelopeExplicit(indE);
            [~,indL] = max(abs(actLocusExplicit(:)));
            actLocusExplicit = actLocusExplicit / actLocusExplicit(indL);

            if sign(actImplicitExplicit(1))<0
                actImplicitExplicit = -actImplicitExplicit;
            end

            if abs(psi-actExplicitExplicit(end))>sonic.Tolerances.CompZero
                actExplicitExplicit(end) = pi+actExplicitExplicit(end);
            end

            % test explicit
            testCase.verifyEqual(actImplicitExplicit,expImplicit,AbsTol=tolerance);
            testCase.verifyEqual(actEnvelopeExplicit,expEnvelope,AbsTol=tolerance);
            testCase.verifyEqual(actLocusExplicit,expLocus,AbsTol=tolerance);
            testCase.verifyEqual(actExplicitExplicit,expExplicit,AbsTol=tolerance);
            testCase.verifyEqual(actProperExplicit,expProper);
            testCase.verifyEqual(actStableExplicit,expStable);
            testCase.verifyMatches(actTypeExplicit,expType);
            testCase.verifyEqual(actExplicitImplicit,expExplicit,AbsTol=tolerance);


            %% create conic from locus
            conicLocus = sonic.Conic(locus,"locus");

            actImplicitLocus= conicLocus.implicit;
            actEnvelopeLocus = conicLocus.envelope;
            actLocusLocus = conicLocus.locus;
            actExplicitLocus = conicLocus.explicit;
            actProperLocus = conicLocus.proper;
            actStableLocus=conicLocus.stable;
            actTypeLocus = conicLocus.type;


            [~,indE] = max(abs(actEnvelopeLocus(:)));
            actEnvelopeLocus = actEnvelopeLocus / actEnvelopeLocus(indE);
            [~,indL] = max(abs(actLocusLocus(:)));
            actLocusLocus = actLocusLocus / actLocusLocus(indL);


            % if necessary, normalize
            if sign(actImplicitLocus(1))<0
                actImplicitLocus = -actImplicitLocus;
            end

            if abs(psi-actExplicitLocus(end))>sonic.Tolerances.CompZero
                actExplicitLocus(end) = pi+actExplicitLocus(end);
            end

            % test locus
            testCase.verifyEqual(actImplicitLocus,expImplicit,AbsTol=tolerance);
            testCase.verifyEqual(actEnvelopeLocus,expEnvelope,AbsTol=tolerance);
            testCase.verifyEqual(actLocusLocus,expLocus,AbsTol=tolerance);
            testCase.verifyEqual(actProperLocus,expProper);
            testCase.verifyEqual(actStableLocus,expStable);
            testCase.verifyMatches(actTypeLocus,expType);
            testCase.verifyEqual(actExplicitLocus,expExplicit,AbsTol=tolerance);


            %% create conic from envelope
            conicEnvelope = sonic.Conic(envelope,"envelope");

            % extract properties
            actImplicitEnvelope = conicEnvelope.implicit;
            actLocusEnvelope = conicEnvelope.locus;
            actEnvelopeEnvelope = conicEnvelope.envelope;
            actExplicitEnvelope = conicEnvelope.explicit;
            actTypeEnvelope = conicEnvelope.type;
            actStableEnvelope=conicEnvelope.stable;
            actProperEnvelope = conicEnvelope.proper;

            [~,indE] = max(abs(actEnvelopeEnvelope(:)));
            actEnvelopeEnvelope = actEnvelopeEnvelope / actEnvelopeEnvelope(indE);
            [~,indL] = max(abs(actLocusEnvelope(:)));
            actLocusEnvelope = actLocusEnvelope / actLocusEnvelope(indL);


            % if necessary, normalize
            if sign(actImplicitEnvelope(1))<0
                actImplicitEnvelope = -actImplicitEnvelope;
            end

            if abs(psi-actExplicitEnvelope(end))>sonic.Tolerances.CompZero
                actExplicitEnvelope(end) = pi+actExplicitEnvelope(end);
            end

            % test envelope
            testCase.verifyEqual(actImplicitEnvelope,expImplicit,AbsTol=tolerance);
            testCase.verifyEqual(actEnvelopeEnvelope,expEnvelope,AbsTol=tolerance);
            testCase.verifyEqual(actLocusEnvelope,expLocus,AbsTol=tolerance);
            testCase.verifyEqual(actProperEnvelope,expProper);
            testCase.verifyEqual(actStableEnvelope,expStable);
            testCase.verifyMatches(actTypeEnvelope,expType);
            testCase.verifyEqual(actExplicitEnvelope,expExplicit,AbsTol=tolerance);


        end

        function ConicParabolaInitTest(testCase)
            % test creation of conic object for a parabola

            xc = 5;
            yc = 1;
            a = Inf;
            b = Inf;
            psi = 50*pi/180;
            p = 5;

            A = sin(psi)^2;
            C = cos(psi)^2;
            B = -cos(psi)*sin(psi);
            D = -A*xc-B*yc-p*cos(psi);
            E = -C*yc-B*xc-p*sin(psi);
            F = A*xc^2+C*yc^2+2*B*xc*yc+2*p*(xc*cos(psi)+yc*sin(psi));

            B = 2*B;
            D = 2*D;
            E = 2*E;

            expProper = true;
            expStable = true;
            expType = "parabola";
            tolerance = sonic.Tolerances.CompZero;

            explicit = [xc;yc;a;b;psi];
            implicit = [A,B,C,D,E,F]';
            locus = [A,B/2,D/2;
                B/2,C,E/2;
                D/2,E/2,F];
            envelope = adjoint(locus);
            expExplicit = [xc;yc;a;b;psi];



            expEnvelope = sonic.Math.A_toDet1(envelope);
            expLocus = sonic.Math.A_toDet1(locus);

            % normalize explicit and implicit for comparison
            expImplicit = implicit/norm(implicit);
            if sign(expImplicit(1))<0
                expImplicit = -expImplicit;
            end

            %% create conic from implicit
            conicImplicit = sonic.Conic(implicit);

            % extract properties
            actImplicitImplicit = conicImplicit.implicit;
            actEnvelopeImplicit = conicImplicit.envelope;
            actLocusImplicit = conicImplicit.locus;
            actExplicitImplicit = conicImplicit.explicit;
            actProperImplicit = conicImplicit.proper;
            actTypeImplicit = conicImplicit.type;
            actStableImplicit=conicImplicit.stable;

            if sign(actImplicitImplicit(1))<0
                actImplicitImplicit = -actImplicitImplicit;
            end

            if abs(psi-actExplicitImplicit(end))>sonic.Tolerances.CompZero
                actExplicitImplicit(end) = pi+actExplicitImplicit(end);
            end

            % test implicit
            testCase.verifyEqual(actImplicitImplicit,expImplicit,AbsTol=tolerance);
            testCase.verifyEqual(actEnvelopeImplicit,expEnvelope,AbsTol=tolerance);
            testCase.verifyEqual(actLocusImplicit,expLocus,AbsTol=tolerance);
            testCase.verifyEqual(actProperImplicit,expProper);
            testCase.verifyEqual(actStableImplicit,expStable);
            testCase.verifyMatches(actTypeImplicit,expType);

            testCase.verifyEqual(actExplicitImplicit([1,2,5]),expExplicit([1,2,5]),AbsTol=tolerance);
            testCase.verifyTrue(isinf(actExplicitImplicit(3)));
            testCase.verifyTrue(isinf(actExplicitImplicit(4)));

            %% create conic from explicit
            testCase.verifyError(@()sonic.Conic(explicit),'sonic:Conic:explicitToLocus:invalidInput');


            %% create conic from locus
            conicLocus = sonic.Conic(locus,"locus");

            actImplicitLocus= conicLocus.implicit;
            actEnvelopeLocus = conicLocus.envelope;
            actLocusLocus = conicLocus.locus;
            actExplicitLocus = conicLocus.explicit;
            actProperLocus = conicLocus.proper;
            actStableLocus=conicLocus.stable;
            actTypeLocus = conicLocus.type;

            % if necessary, normalize
            if sign(actImplicitLocus(1))<0
                actImplicitLocus = -actImplicitLocus;
            end

            if abs(psi-actExplicitLocus(end))>sonic.Tolerances.CompZero
                actExplicitLocus(end) = pi+actExplicitLocus(end);
            end

            % test locus
            testCase.verifyEqual(actImplicitLocus,expImplicit,AbsTol=tolerance);
            testCase.verifyEqual(actEnvelopeLocus,expEnvelope,AbsTol=tolerance);
            testCase.verifyEqual(actLocusLocus,expLocus,AbsTol=tolerance);
            testCase.verifyEqual(actProperLocus,expProper);
            testCase.verifyEqual(actStableLocus,expStable);
            testCase.verifyMatches(actTypeLocus,expType);

            testCase.verifyEqual(actExplicitLocus([1,2,5]),expExplicit([1,2,5]),AbsTol=tolerance);
            testCase.verifyTrue(isinf(actExplicitLocus(3)));
            testCase.verifyTrue(isinf(actExplicitLocus(4)));


            %% create conic from envelope
            conicEnvelope = sonic.Conic(envelope,"envelope");

            % extract properties
            actImplicitEnvelope = conicEnvelope.implicit;
            actLocusEnvelope = conicEnvelope.locus;
            actEnvelopeEnvelope = conicEnvelope.envelope;
            actExplicitEnvelope = conicEnvelope.explicit;
            actTypeEnvelope = conicEnvelope.type;
            actStableEnvelope=conicEnvelope.stable;
            actProperEnvelope = conicEnvelope.proper;

            % if necessary, normalize
            if sign(actImplicitEnvelope(1))<0
                actImplicitEnvelope = -actImplicitEnvelope;
            end

            if abs(psi-actExplicitEnvelope(end))>sonic.Tolerances.CompZero
                actExplicitEnvelope(end) = pi+actExplicitEnvelope(end);
            end

            % test envelope
            testCase.verifyEqual(actImplicitEnvelope,expImplicit,AbsTol=tolerance);
            testCase.verifyEqual(actEnvelopeEnvelope,expEnvelope,AbsTol=tolerance);
            testCase.verifyEqual(actLocusEnvelope,expLocus,AbsTol=tolerance);
            testCase.verifyEqual(actProperEnvelope,expProper);
            testCase.verifyEqual(actStableEnvelope,expStable);
            testCase.verifyMatches(actTypeEnvelope,expType);

            testCase.verifyEqual(actExplicitEnvelope([1,2,5]),expExplicit([1,2,5]),AbsTol=tolerance);
            testCase.verifyTrue(isinf(actExplicitEnvelope(3)));
            testCase.verifyTrue(isinf(actExplicitEnvelope(4)));
        end

        function ConicUnstableParabolaInitTest(testCase)
            % test construction of a conic object for numerically unstable
            % parabola

            xc = 0;
            yc = 0;

            a = Inf;
            b = Inf;
            psi = 50*pi/180;
            p = 5;

            A = sin(psi)^2;
            C = cos(psi)^2;
            B = -cos(psi)*sin(psi);
            D = -A*xc-B*yc-p*cos(psi);
            E = -C*yc-B*xc-p*sin(psi);
            F = A*xc^2+C*yc^2+2*B*xc*yc+2*p*(xc*cos(psi)+yc*sin(psi));

            locus = [A,B,D;
                B,C,E;
                D,E,F];

            envelope = adjoint(locus);

            vertex = [50000; 77500; 1];
            transl = [-vertex(1);-vertex(2);1];

            % Translate the conic envelope
            P = [transl(3),0,-transl(1);
                0,transl(3),-transl(2);
                0,0,transl(3)];

            envelope = P*envelope*P';
            envelope = sonic.Math.A_toDet1(envelope);


            % translate conic locus
            P_star = [transl(3), 0, transl(1);
                0, transl(3), transl(2);
                0, 0, transl(3)];

            locus = P_star'*locus*P_star;
            locus = sonic.Math.A_toDet1(locus);

            explicit =  [vertex(1);vertex(2);a;b;psi];
            expExplicit = [vertex(1);vertex(2);a;b;psi];

            implicit = [locus(1,1),2*locus(1,2),locus(2,2),locus(1,3)*2,2*locus(2,3),locus(3,3)]';


            expProper = true;
            expStable = false;
            expType = "parabola";

            tolerance = sonic.Tolerances.RelaxedCompZero;

            [~,indL] = max(abs(locus(:)));
            expLocus= locus / locus(indL);
            [~,indE] = max(abs(envelope(:)));
            expEnvelope = envelope / envelope(indE);

            % normalize explicit and implicit for comparison
            expImplicit = implicit/norm(implicit);
            if sign(expImplicit(1))<0
                expImplicit = -expImplicit;
            end

            %% create conic from implicit
            conicImplicit = sonic.Conic(implicit);

            % extract properties
            actImplicitImplicit = conicImplicit.implicit;
            actEnvelopeImplicit = conicImplicit.envelope;
            actLocusImplicit = conicImplicit.locus;
            actExplicitImplicit = conicImplicit.explicit;
            actProperImplicit = conicImplicit.proper;
            actTypeImplicit = conicImplicit.type;
            actStableImplicit=conicImplicit.stable;


            [~,indE] = max(abs(actEnvelopeImplicit(:)));
            actEnvelopeImplicit = actEnvelopeImplicit / actEnvelopeImplicit(indE);
            [~,indL] = max(abs(actLocusImplicit(:)));
            actLocusImplicit = actLocusImplicit / actLocusImplicit(indL);

            if sign(actImplicitImplicit(1))<0
                actImplicitImplicit = -actImplicitImplicit;
            end

            if abs(psi-actExplicitImplicit(end))>sonic.Tolerances.CompZero
                actExplicitImplicit(end) = pi+actExplicitImplicit(end);
            end

            % test implicit
            testCase.verifyEqual(actImplicitImplicit,expImplicit,AbsTol=tolerance);
            testCase.verifyEqual(actEnvelopeImplicit,expEnvelope,AbsTol=tolerance);
            testCase.verifyEqual(actLocusImplicit,expLocus,AbsTol=tolerance);
            testCase.verifyEqual(actProperImplicit,expProper);
            testCase.verifyEqual(actStableImplicit,expStable);
            testCase.verifyMatches(actTypeImplicit,expType);

            testCase.verifyEqual(actExplicitImplicit([1,2,5]),expExplicit([1,2,5]),AbsTol=tolerance);
            testCase.verifyTrue(isinf(actExplicitImplicit(3)));
            testCase.verifyTrue(isinf(actExplicitImplicit(4)));

            %% create conic from explicit
            testCase.verifyError(@()sonic.Conic(explicit),'sonic:Conic:explicitToLocus:invalidInput');


            %% create conic from locus
            conicLocus = sonic.Conic(locus,"locus");

            actImplicitLocus= conicLocus.implicit;
            actEnvelopeLocus = conicLocus.envelope;
            actLocusLocus = conicLocus.locus;
            actExplicitLocus = conicLocus.explicit;
            actProperLocus = conicLocus.proper;
            actStableLocus=conicLocus.stable;
            actTypeLocus = conicLocus.type;


            [~,indE] = max(abs(actEnvelopeLocus(:)));
            actEnvelopeLocus = actEnvelopeLocus / actEnvelopeLocus(indE);
            [~,indL] = max(abs(actLocusLocus(:)));
            actLocusLocus = actLocusLocus / actLocusLocus(indL);


            % if necessary, normalize
            if sign(actImplicitLocus(1))<0
                actImplicitLocus = -actImplicitLocus;
            end

            if abs(psi-actExplicitLocus(end))>sonic.Tolerances.CompZero
                actExplicitLocus(end) = pi+actExplicitLocus(end);
            end

            % test locus
            testCase.verifyEqual(actImplicitLocus,expImplicit,AbsTol=tolerance);
            testCase.verifyEqual(actEnvelopeLocus,expEnvelope,AbsTol=tolerance);
            testCase.verifyEqual(actLocusLocus,expLocus,AbsTol=tolerance);
            testCase.verifyEqual(actProperLocus,expProper);
            testCase.verifyEqual(actStableLocus,expStable);
            testCase.verifyMatches(actTypeLocus,expType);

            testCase.verifyEqual(actExplicitLocus([1,2,5]),expExplicit([1,2,5]),AbsTol=tolerance);
            testCase.verifyTrue(isinf(actExplicitLocus(3)));
            testCase.verifyTrue(isinf(actExplicitLocus(4)));


            %% create conic from envelope
            conicEnvelope = sonic.Conic(envelope,"envelope");

            % extract properties
            actImplicitEnvelope = conicEnvelope.implicit;
            actLocusEnvelope = conicEnvelope.locus;
            actEnvelopeEnvelope = conicEnvelope.envelope;
            actExplicitEnvelope = conicEnvelope.explicit;
            actTypeEnvelope = conicEnvelope.type;
            actStableEnvelope=conicEnvelope.stable;
            actProperEnvelope = conicEnvelope.proper;

            [~,indE] = max(abs(actEnvelopeEnvelope(:)));
            actEnvelopeEnvelope = actEnvelopeEnvelope / actEnvelopeEnvelope(indE);
            [~,indL] = max(abs(actLocusEnvelope(:)));
            actLocusEnvelope = actLocusEnvelope / actLocusEnvelope(indL);

            % if necessary, normalize
            if sign(actImplicitEnvelope(1))<0
                actImplicitEnvelope = -actImplicitEnvelope;
            end

            if abs(psi-actExplicitEnvelope(end))>sonic.Tolerances.CompZero
                actExplicitEnvelope(end) = pi+actExplicitEnvelope(end);
            end

            % test envelope
            testCase.verifyEqual(actImplicitEnvelope,expImplicit,AbsTol=tolerance);
            testCase.verifyEqual(actEnvelopeEnvelope,expEnvelope,AbsTol=tolerance);
            testCase.verifyEqual(actLocusEnvelope,expLocus,AbsTol=tolerance);
            testCase.verifyEqual(actProperEnvelope,expProper);
            testCase.verifyEqual(actStableEnvelope,expStable);
            testCase.verifyMatches(actTypeEnvelope,expType);

            testCase.verifyEqual(actExplicitEnvelope([1,2,5]),expExplicit([1,2,5]),AbsTol=tolerance);
            testCase.verifyTrue(isinf(actExplicitEnvelope(3)));
            testCase.verifyTrue(isinf(actExplicitEnvelope(4)));


        end

        function ConicCircleInitTest(testCase)
            % test creation of conic object for ellipse

            xc = 1;
            yc = 2;
            a = 3;
            b = 3;
            psi = 0;

            A = a^2*sin(psi)^2 + b^2*cos(psi)^2;
            B = 2*(b^2 - a^2)*cos(psi)*sin(psi);
            C = a^2*cos(psi)^2 + b^2*sin(psi)^2;
            D = -2*A*xc - B*yc;
            E = -B*xc - 2*C*yc;
            F = A*xc^2 + B*xc*yc + C*yc^2 - a^2*b^2;

            expType = "circle";
            expProper = true;
            expStable = true;

            explicit = [xc;yc;a;b;psi];
            implicit = [A,B,C,D,E,F]';
            locus = [A,B/2,D/2;
                B/2,C,E/2;
                D/2,E/2,F];
            envelope = adjoint(locus);
            expExplicit = [xc;yc;a;b;psi];


            tolerance = sonic.Tolerances.CompZero;


            expEnvelope = sonic.Math.A_toDet1(envelope);
            expLocus = sonic.Math.A_toDet1(locus);

            % normalize explicit and implicit for comparison
            expImplicit = implicit/norm(implicit);
            if sign(expImplicit(1))<0
                expImplicit = -expImplicit;
            end

            %% create conic from implicit
            conicImplicit = sonic.Conic(implicit);

            % extract properties
            actImplicitImplicit = conicImplicit.implicit;
            actEnvelopeImplicit = conicImplicit.envelope;
            actLocusImplicit = conicImplicit.locus;
            actExplicitImplicit = conicImplicit.explicit;
            actProperImplicit = conicImplicit.proper;
            actTypeImplicit = conicImplicit.type;
            actStableImplicit=conicImplicit.stable;

            if sign(actImplicitImplicit(1))<0
                actImplicitImplicit = -actImplicitImplicit;
            end

            if abs(psi-actExplicitImplicit(end))>sonic.Tolerances.CompZero
                actExplicitImplicit(end) = pi+actExplicitImplicit(end);
            end

            % test implicit
            testCase.verifyEqual(actImplicitImplicit,expImplicit,AbsTol=tolerance);
            testCase.verifyEqual(actEnvelopeImplicit,expEnvelope,AbsTol=tolerance);
            testCase.verifyEqual(actLocusImplicit,expLocus,AbsTol=tolerance);
            testCase.verifyEqual(actProperImplicit,expProper);
            testCase.verifyEqual(actStableImplicit,expStable);
            testCase.verifyMatches(actTypeImplicit,expType);
            testCase.verifyEqual(actExplicitImplicit,expExplicit,AbsTol=tolerance);

            %% create conic from explicit
            conicExplicit = sonic.Conic(explicit);

            % extract properties
            actImplicitExplicit = conicExplicit.implicit;
            actEnvelopeExplicit = conicExplicit.envelope;
            actLocusExplicit = conicExplicit.locus;
            actExplicitExplicit = conicExplicit.explicit;
            actProperExplicit = conicExplicit.proper;
            actTypeExplicit = conicExplicit.type;
            actStableExplicit=conicExplicit.stable;

            if sign(actImplicitExplicit(1))<0
                actImplicitExplicit = -actImplicitExplicit;
            end

            if abs(psi-actExplicitExplicit(end))>sonic.Tolerances.CompZero
                actExplicitExplicit(end) = pi+actExplicitExplicit(end);
            end

            % test explicit
            testCase.verifyEqual(actImplicitExplicit,expImplicit,AbsTol=tolerance);
            testCase.verifyEqual(actEnvelopeExplicit,expEnvelope,AbsTol=tolerance);
            testCase.verifyEqual(actLocusExplicit,expLocus,AbsTol=tolerance);
            testCase.verifyEqual(actExplicitExplicit,expExplicit,AbsTol=tolerance);
            testCase.verifyEqual(actProperExplicit,expProper);
            testCase.verifyEqual(actStableExplicit,expStable);
            testCase.verifyMatches(actTypeExplicit,expType);
            testCase.verifyEqual(actExplicitImplicit,expExplicit,AbsTol=tolerance);


            %% create conic from locus
            conicLocus = sonic.Conic(locus,"locus");

            actImplicitLocus= conicLocus.implicit;
            actEnvelopeLocus = conicLocus.envelope;
            actLocusLocus = conicLocus.locus;
            actExplicitLocus = conicLocus.explicit;
            actProperLocus = conicLocus.proper;
            actStableLocus=conicLocus.stable;
            actTypeLocus = conicLocus.type;

            % if necessary, normalize
            if sign(actImplicitLocus(1))<0
                actImplicitLocus = -actImplicitLocus;
            end

            if abs(psi-actExplicitLocus(end))>sonic.Tolerances.CompZero
                actExplicitLocus(end) = pi+actExplicitLocus(end);
            end

            % test locus
            testCase.verifyEqual(actImplicitLocus,expImplicit,AbsTol=tolerance);
            testCase.verifyEqual(actEnvelopeLocus,expEnvelope,AbsTol=tolerance);
            testCase.verifyEqual(actLocusLocus,expLocus,AbsTol=tolerance);
            testCase.verifyEqual(actProperLocus,expProper);
            testCase.verifyEqual(actStableLocus,expStable);
            testCase.verifyMatches(actTypeLocus,expType);
            testCase.verifyEqual(actExplicitLocus,expExplicit,AbsTol=tolerance);


            %% create conic from envelope
            conicEnvelope = sonic.Conic(envelope,"envelope");

            % extract properties
            actImplicitEnvelope = conicEnvelope.implicit;
            actLocusEnvelope = conicEnvelope.locus;
            actEnvelopeEnvelope = conicEnvelope.envelope;
            actExplicitEnvelope = conicEnvelope.explicit;
            actTypeEnvelope = conicEnvelope.type;
            actStableEnvelope=conicEnvelope.stable;
            actProperEnvelope = conicEnvelope.proper;

            % if necessary, normalize
            if sign(actImplicitEnvelope(1))<0
                actImplicitEnvelope = -actImplicitEnvelope;
            end

            if abs(psi-actExplicitEnvelope(end))>sonic.Tolerances.CompZero
                actExplicitEnvelope(end) = pi+actExplicitEnvelope(end);
            end

            % test envelope
            testCase.verifyEqual(actImplicitEnvelope,expImplicit,AbsTol=tolerance);
            testCase.verifyEqual(actEnvelopeEnvelope,expEnvelope,AbsTol=tolerance);
            testCase.verifyEqual(actLocusEnvelope,expLocus,AbsTol=tolerance);
            testCase.verifyEqual(actProperEnvelope,expProper);
            testCase.verifyEqual(actStableEnvelope,expStable);
            testCase.verifyMatches(actTypeEnvelope,expType);
            testCase.verifyEqual(actExplicitEnvelope,expExplicit,AbsTol=tolerance);


        end

        function ConicInitTest(testCase)

            % test input errors
            mat = [1,2,3;
                2,1,5;
                3,5,8];

            testCase.verifyError(@()sonic.Conic(mat),'sonic:Conic:invalidInput');
            testCase.verifyError(@()sonic.Conic(mat,"ENVELOPE"),'sonic:Conic:invalidInput');

            % single line not accepted
            testCase.verifyError(@()sonic.Conic(sonic.Lines2([3;2;1])),'sonic:Conic:invalidInput');

            % degenerate conic not accepted as envelope
            envelope = [100.0000  200.0000  130.0000;
                200.0000  400.0000  260.0000;
                130.0000  260.0000  169.0000];
            testCase.verifyError(@()sonic.Conic(envelope,"envelope"),'sonic:Conic:invalidInput');

            % wrong dimension of input
            mat = [1,2,3];
            testCase.verifyError(@()sonic.Conic(mat),'sonic:Conic:invalidInput');

            % numerical issues cannot be avoided
            xc = 75000000;
            yc = 50000000;
            a = 1;
            b = 0.2;
            psi = pi/2;

            explicit = [xc;yc;a;b;psi];

            testCase.verifyError(@()sonic.Conic(explicit),'sonic:Conic:invalidInput');


            explicitDegenerate = [1,2,Inf,Inf,pi/3];
            testCase.verifyError(@()sonic.Conic(explicitDegenerate),'sonic:Conic:explicitToLocus:invalidInput');


        end

        function ConicParabolaToLocusTest(testCase)

            xc = 5;
            yc = 1;
            psi = 50*pi/180;
            p = 5;

            parabola = [xc;yc;psi;p];


            A = sin(psi)^2;
            C = cos(psi)^2;
            B = -cos(psi)*sin(psi);
            D = -A*xc-B*yc-p*cos(psi);
            E = -C*yc-B*xc-p*sin(psi);
            F = A*xc^2+C*yc^2+2*B*xc*yc+2*p*(xc*cos(psi)+yc*sin(psi));

            locus = [A,B,D;
                B,C,E;
                D,E,F];

            expLocus = sonic.Math.A_toDet1(locus);

            actLocus = sonic.Conic.parabolaToLocus(parabola);

            testCase.verifyEqual(expLocus,actLocus,AbsTol = sonic.Tolerances.CompZero);

        end

        function ConicLocusToParameterTest(testCase)

            % test parabola
            xc = 10;
            yc = 20;
            psi = 50*pi/180;
            p = 5;

            A = sin(psi)^2;
            C = cos(psi)^2;
            B = -cos(psi)*sin(psi);
            D = -A*xc-B*yc-p*cos(psi);
            E = -C*yc-B*xc-p*sin(psi);
            F = A*xc^2+C*yc^2+2*B*xc*yc+2*p*(xc*cos(psi)+yc*sin(psi));

            locusParabola = [A,B,D;
                B,C,E;
                D,E,F];

            expParameter = p;

            actParameter = sonic.Conic.locusToParameter(locusParabola);

            % test ellipse
            xcE = 1;
            ycE = 2;
            aE = 4;
            bE = 3;
            psiE = pi/6;

            AE = aE^2*sin(psiE)^2 + bE^2*cos(psiE)^2;
            BE = 2*(bE^2 - aE^2)*cos(psiE)*sin(psiE);
            CE = aE^2*cos(psiE)^2 + bE^2*sin(psiE)^2;
            DE = -2*AE*xcE - BE*ycE;
            EE = -BE*xcE - 2*CE*ycE;
            FE = AE*xcE^2 + BE*xcE*ycE + CE*ycE^2 - aE^2*bE^2;

            locusEllipse = [AE,BE/2,DE/2;
                BE/2,CE,EE/2;
                DE/2,EE/2,FE];

            expParameterEllipse = bE^2/aE;
            actParameterEllipse = sonic.Conic.locusToParameter(locusEllipse);


            % test unstable ellipse
            xcEU = 1;
            ycEU = 2;
            aEU = 4;
            bEU = 3;
            psiEU = pi/6;

            AEU = aEU^2*sin(psiEU)^2 + bEU^2*cos(psiEU)^2;
            BEU = 2*(bEU^2 - aEU^2)*cos(psiEU)*sin(psiEU);
            CEU = aEU^2*cos(psiEU)^2 + bEU^2*sin(psiEU)^2;
            DEU = -2*AEU*xcEU - BEU*ycEU;
            EEU = -BEU*xcEU - 2*CEU*ycEU;
            FEU = AEU*xcEU^2 + BEU*xcEU*ycEU + CEU*ycEU^2 - aEU^2*bEU^2;

            locusEllipseUnstable = [AEU,BEU/2,DEU/2;
                BEU/2,CEU,EEU/2;
                DEU/2,EEU/2,FEU];

            expParameterEllipseUnstable = bEU^2/aEU;
            actParameterEllipseUnstable = sonic.Conic.locusToParameter(locusEllipseUnstable);


            % test hyperbola
            xcH = 5;
            ycH = 10;
            aH = -20;
            bH = -8;
            psiH = pi/6;

            AH = bH^2*cos(psiH)^2 - aH^2*sin(psiH)^2;
            BH= 2*(aH^2+bH^2)*cos(psiH)*sin(psiH);
            CH = bH^2*sin(psiH)^2 - aH^2*cos(psiH)^2;
            DH = -2*AH*xcH - BH*ycH;
            EH = -BH*xcH - 2*CH*ycH;
            FH = AH*xcH^2 + BH*xcH*ycH + CH*ycH^2 - aH^2*bH^2;

            locusHyperbola = [AH,BH/2,DH/2;
                BH/2,CH,EH/2;
                DH/2,EH/2,FH];


            expParameterHyperbola = -bH^2/aH;
            actParameterHyperbola = sonic.Conic.locusToParameter(locusHyperbola);

            locusDegenerate = [-1,0,1;
                0,1,-2;
                1,-2,3];

            testCase.verifyEqual(actParameter,expParameter,AbsTol = sonic.Tolerances.CompZero);
            testCase.verifyEqual(actParameterEllipse,expParameterEllipse,AbsTol = sonic.Tolerances.CompZero);
            testCase.verifyEqual(actParameterEllipseUnstable,expParameterEllipseUnstable,AbsTol = sonic.Tolerances.CompZero);
            testCase.verifyEqual(actParameterHyperbola,expParameterHyperbola,AbsTol = sonic.Tolerances.CompZero);
            testCase.verifyError(@()sonic.Conic.locusToParameter(locusDegenerate),'sonic:Conic:locusToParameter:invalidInput')
        end

        function ConicCenterTest(testCase)
            % test function for conic object translation

            % describe initial conic in explicit form
            firstExplicit = [-5;-4;3;2;pi/3];

            % build the initial conic object
            firstConic = sonic.Conic(firstExplicit);

            % the expected explicit conic has center at the origin
            secondExplicit = [0;0;3;2;pi/3];

            % build the expected conic object
            expConic = sonic.Conic(secondExplicit);

            % extract expected properties
            expImpl = expConic.implicit;
            expExpl = expConic.explicit;
            expLocus = expConic.locus;
            expEnvelope = expConic.envelope;
            expType = expConic.type;
            expProper = expConic.proper;

            % build actual conic object
            % actConic = sonic.Conic.center(firstConic);
            actConic = firstConic.center();

            % extract actual properties
            actImpl = actConic.implicit;
            actExpl = actConic.explicit;
            actLocus = actConic.locus;
            actEnvelope = actConic.envelope;
            actType = actConic.type;
            actProper = actConic.proper;
            actStable = actConic.stable;

            % remove sign ambiguity
            if expImpl(1)*actImpl(1)<0
                expImpl = -expImpl;
            end

            % remove pi ambiguity
            if actExpl(end)*expExpl(end)<0
                expExpl(end) = expExpl(end)-pi;
            end

            testCase.verifyEqual(actImpl,expImpl,AbsTol = sonic.Tolerances.CompZero)
            testCase.verifyEqual(actExpl,expExpl,AbsTol = sonic.Tolerances.CompZero)
            testCase.verifyEqual(actLocus,expLocus,AbsTol = sonic.Tolerances.CompZero)
            testCase.verifyEqual(actEnvelope,expEnvelope,AbsTol = sonic.Tolerances.CompZero)
            testCase.verifyMatches(actType,expType);
            testCase.verifyEqual(expProper,actProper);
            testCase.verifyTrue(actStable);
        end

        function ConicRotateTest(testCase)
            % test function for conic object translation

            % describe initial conic in explicit form
            firstExplicit = [0;0;3;2;0];

            % the expected rotated version
            secondExplicit = [0;0;3;2;pi/3];

            % build the initial conic object
            firstConic = sonic.Conic(firstExplicit);

            % build the expected conic object
            expConic = sonic.Conic(secondExplicit);

            % extract expected properties
            expImpl = expConic.implicit;
            expExpl = expConic.explicit;
            expLocus = expConic.locus;
            expEnvelope = expConic.envelope;
            expType = expConic.type;
            expProper = expConic.proper;

            % build actual conic object
            % actConic = sonic.Conic.center(firstConic);
            actConic = firstConic.rotate(pi/3);

            % extract actual properties
            actImpl = actConic.implicit;
            actExpl = actConic.explicit;
            actLocus = actConic.locus;
            actEnvelope = actConic.envelope;
            actType = actConic.type;
            actProper = actConic.proper;
            actStable = actConic.stable;

            % remove sign ambiguity
            if expImpl(1)*actImpl(1)<0
                expImpl = -expImpl;
            end

            if actExpl(end)*expExpl(end)<0 
                expExpl(end) = expExpl(end)-pi;
            end

            testCase.verifyEqual(actImpl,expImpl,AbsTol = sonic.Tolerances.CompZero)
            testCase.verifyEqual(actExpl,expExpl,AbsTol = sonic.Tolerances.CompZero)
            testCase.verifyEqual(actLocus,expLocus,AbsTol = sonic.Tolerances.CompZero)
            testCase.verifyEqual(actEnvelope,expEnvelope,AbsTol = sonic.Tolerances.CompZero)
            testCase.verifyMatches(actType,expType);
            testCase.verifyEqual(expProper,actProper);
            testCase.verifyTrue(actStable);
        end

        function ConicMeetConicConicTest(testCase)

            % two intersecting conics
            explicit1 = [1,2,4,3,pi/6];

            explicit2 = [-2,-4,7,5,pi/2];

            conic1 = sonic.Conic(explicit1);
            locus1 = conic1.locus;

            conic2 = sonic.Conic(explicit2);
            locus2 = conic2.locus;


            intPoints = sonic.Conic.meetConicConic(locus1,locus2);

            intPointsCoords = intPoints.p2;

            % two conics tangent in one point
            implicit3 = [0,0,1,-12,-8,52]';
            implicit4 = [1,0,1,-14,-18,128];

            conic3 = sonic.Conic(implicit3);
            locus3 = conic3.locus;

            conic4 = sonic.Conic(implicit4);
            locus4 = conic4.locus;

            intTangPoints = sonic.Conic.meetConicConic(locus3,locus4);

            intTangPointsCoords = intTangPoints.p2;

            % two conics tangent in two points

            implicit5 = [1,0,3,-4,-2,3]';
            implicit6 = [1,0,11,-4,-2,3]';

            conic5 = sonic.Conic(implicit5);
            locus5 = conic5.locus;

            conic6 = sonic.Conic(implicit6);
            locus6 = conic6.locus;

            tangPoints = sonic.Conic.meetConicConic(locus5,locus6);

            tangPointsCoords = tangPoints.p2;

            % check if they lie on both the conics
            liesOnConic1 = zeros(4,1);
            liesOnConic2 = zeros(4,1);
            liesOnConic3 = zeros(4,1);
            liesOnConic4 = zeros(4,1);
            liesOnConic5 = zeros(4,1);
            liesOnConic6 = zeros(4,1);

            for i = 1 : 4
                liesOnConic1(i) = transpose(intPointsCoords(:,i))*conic1.locus*intPointsCoords(:,i);
                liesOnConic2(i) = transpose(intPointsCoords(:,i))*conic2.locus*intPointsCoords(:,i);
                liesOnConic3(i) = transpose(intTangPointsCoords(:,i))*locus3*intTangPointsCoords(:,i);
                liesOnConic4(i) = transpose(intTangPointsCoords(:,i))*locus4*intTangPointsCoords(:,i);
                liesOnConic5(i) = transpose(tangPointsCoords(:,i))*locus5*tangPointsCoords(:,i);
                liesOnConic6(i) = transpose(tangPointsCoords(:,i))*locus6*tangPointsCoords(:,i);
            end

            testCase.verifyTrue(all(abs(liesOnConic1)<sonic.Tolerances.RelaxedCompZero));
            testCase.verifyTrue(all(abs(liesOnConic2)<sonic.Tolerances.RelaxedCompZero));
            testCase.verifyTrue(all(abs(liesOnConic3)<sonic.Tolerances.RelaxedCompZero));
            testCase.verifyTrue(all(abs(liesOnConic4)<sonic.Tolerances.RelaxedCompZero));

            degLocus = [9, 36/2,3;
                36/2,36,6;
                3,6,1];
            testCase.verifyError(@() sonic.Conic.meetConicConic(degLocus,locus1),'sonic:Conic:meetConicConic:degenerateConics')
        end

        function ConicMeetLineConicTest(testCase)

            % line and conic in generic position

            explicit1 = [1,2,4,3,pi/6];
            conic1 = sonic.Conic(explicit1);

            locus1 = conic1.locus;

            line1 = [1;1;1];
            line1 = sonic.Lines2(line1);

            line2 = [0;0;1];

            intPoints1 = sonic.Conic.meetLineConic(line1,locus1);
            intPointsCoords1 = intPoints1.p2;

            % line at infinity and conic
            line2 = sonic.Lines2(line2);

            intPoints2 = sonic.Conic.meetLineConic(line2,locus1);
            intPointsCoords2 = intPoints2.p2;

            % second conic and line at infinity for complete coverage

            explicit3 = [-4,1,4,3,3*pi/2];
            conic3 = sonic.Conic(explicit3);

            locus2 = conic3.locus;
            intPoints3 = sonic.Conic.meetLineConic(line2,locus2);
            intPointsCoords3 = intPoints3.p2;


            % line tangent to conic
            implicit2 = [144,0,288,-1728,-1152,3744]';
            line3 = sonic.Lines2([-576,288,288]');
            conic4 = sonic.Conic(implicit2);
            locus3 = conic4.locus;

            intPoints4 = sonic.Conic.meetLineConic(line3,locus3);
            intPointsCoords4 = intPoints4.p2;


            liesOnConic1 = zeros(2,1);
            liesOnLine1 = zeros(2,1);
            liesOnConic2 = zeros(2,1);
            liesOnLine2 = zeros(2,1);
            liesOnConic3 = zeros(2,1);
            liesOnLine3 = zeros(2,1);
            liesOnConic4 = zeros(2,1);
            liesOnLine4= zeros(2,1);

            for i = 1 : 2
                liesOnConic1(i) = transpose(intPointsCoords1(:,i))*locus1*intPointsCoords1(:,i);
                liesOnLine1(i) = transpose(intPointsCoords1(:,i))*line1.p2;
                liesOnConic2(i) = transpose(intPointsCoords2(:,i))*locus1*intPointsCoords2(:,i);
                liesOnLine2(i) = transpose(intPointsCoords2(:,i))*line2.p2;
                liesOnConic3(i) = transpose(intPointsCoords3(:,i))*locus2*intPointsCoords3(:,i);
                liesOnLine3(i) = transpose(intPointsCoords3(:,i))*line2.p2;

                liesOnConic4(i) = transpose(intPointsCoords4(:,i))*locus3*intPointsCoords4(:,i);
                liesOnLine4(i) = transpose(intPointsCoords4(:,i))*line3.p2;
            end



            testCase.verifyTrue(all(abs(liesOnConic1)<sonic.Tolerances.IsOnConic));
            testCase.verifyTrue(all(abs(liesOnLine1)<sonic.Tolerances.CompZero));
            testCase.verifyTrue(all(abs(liesOnConic2)<sonic.Tolerances.IsOnConic));
            testCase.verifyTrue(all(abs(liesOnLine2)<sonic.Tolerances.CompZero));



        end


    end

end