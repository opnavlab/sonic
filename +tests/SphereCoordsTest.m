classdef SphereCoordsTest < matlab.unittest.TestCase
    properties

    end
    methods(TestClassSetup)
        % Shared setup for the entire test class
    end
    
    methods(TestMethodSetup)
        % Setup for each test
    end
    
    methods(Test)
        % Test methods
        
        function raDecToCartTest(testCase)
            %%% test the simple case of angle conversion
            % Specify points in degrees:
            ra_DEG = [0; 90; 0; 240; -120];
            dec_DEG = [0; 0; 90; -45; -45];

            % Convert to radians:
            ra_RAD = sonic.Units.DEGtoRAD(ra_DEG);
            dec_RAD = sonic.Units.DEGtoRAD(dec_DEG);

            % And convert to cartesian:
            u = sonic.SphereCoords.raDecToCart(ra_RAD, dec_RAD);
            
            u_exp = [1, 0, 0, -sqrt(2)/2*1/2, -sqrt(2)/2*1/2;
                     0, 1, 0, -sqrt(2)/2*sqrt(3)/2, -sqrt(2)/2*sqrt(3)/2;
                     0, 0, 1, -sqrt(2)/2, -sqrt(2)/2];

            testCase.verifyEqual(u, u_exp, AbsTol=sonic.Tolerances.CompZero);

            %%% Make sure the weighting works
            % for vector of weights
            w = [1; 2; 3; 4; 5];
            u = sonic.SphereCoords.raDecToCart(ra_RAD, dec_RAD, w);
            testCase.verifyEqual(u, w'.*u_exp, AbsTol=sonic.Tolerances.CompZero);
            
            % also for scalar weights
            w = 3;
            u = sonic.SphereCoords.raDecToCart(ra_RAD, dec_RAD, w);
            testCase.verifyEqual(u, w.*u_exp, AbsTol=sonic.Tolerances.CompZero);


            %%% Throws exception if mismatched dimensions
            testCase.verifyError( ...
                @() sonic.SphereCoords.raDecToCart([1; 2], [1; 2; 3]), ...
                'sonic:SphereCoords:raDecToCart:mismatchedValues')

            testCase.verifyError( ...
                @() sonic.SphereCoords.raDecToCart([1; 2; 3], [1; 2]), ...
                'sonic:SphereCoords:raDecToCart:mismatchedValues')

            testCase.verifyError( ...
                @() sonic.SphereCoords.raDecToCart([1; 2; 3], [1; 2; 3], [1; 2]), ...
                'sonic:SphereCoords:raDecToCart:mismatchedValues')
        end

        function cartToRaDecTest(testCase)
            % Expected angles
            ra_DEG = [0; 90; 0; -120];
            dec_DEG = [0; 0; 90; -45];

            % Convert to radians
            exp_ra_RAD = sonic.Units.DEGtoRAD(ra_DEG);
            exp_dec_RAD = sonic.Units.DEGtoRAD(dec_DEG);

            % Vectors
            % u = sonic.SphereCoords.raDecToCart(exp_ra_RAD, exp_dec_RAD);
            u = [1, 0, 0, -sqrt(2)/2*1/2;
                 0, 1, 0, -sqrt(2)/2*sqrt(3)/2;
                 0, 0, 1, -sqrt(2)/2];

            % And convert to cartesian:
            [act_ra_RAD, act_dec_RAD] = sonic.SphereCoords.cartToRaDec(u);
            testCase.verifyEqual(act_ra_RAD, exp_ra_RAD, AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(act_dec_RAD, exp_dec_RAD, AbsTol=sonic.Tolerances.CompZero);
            
            % Weighted version
            w = [1, 2, 3, 4];
            [act_ra_RAD, act_dec_RAD] = sonic.SphereCoords.cartToRaDec(w.*u);
            testCase.verifyEqual(act_ra_RAD, exp_ra_RAD, AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(act_dec_RAD, exp_dec_RAD, AbsTol=sonic.Tolerances.CompZero);

            % Weighted version
            w = 3;
            [act_ra_RAD, act_dec_RAD] = sonic.SphereCoords.cartToRaDec(w*u);
            testCase.verifyEqual(act_ra_RAD, exp_ra_RAD, AbsTol=sonic.Tolerances.CompZero);
            testCase.verifyEqual(act_dec_RAD, exp_dec_RAD, AbsTol=sonic.Tolerances.CompZero);
            
        end
    end
    
end