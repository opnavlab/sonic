classdef LunarLambert < sonic.Reflectance
    % This software is made available under the MIT license. See SONIC/LICENSE
    % for details.

    properties
        AL                  double          % Albedo. nxmxk 3D matrix
        modelType           string          % Type of model as determined by the defined coefficients
        coeffs              double          % polynomial coefficients for polyPhase implementation
    end

    methods

        function obj = LunarLambert(AL, modelType, coeffs)
            %% obj = LunarLambert(AL, modelType, coeffs)
            %   Instantiates a sonic LunarLambert object.
            %
            %   Inputs:
            %       - AL (n,m,k): Albedo matrix. Each n,m index contains a
            %           vector of length k that contains the coefficients of an
            %           albedo polynomial in increasing order.
            %       - modelType (string): Implementation type of LunarLambert
            %           model. Current supported inputs are:
            %               - "GaskellApprox"
            %               - "McEwenCubic"
            %               - "polyPhase"
            %       - coeffs (n,m,k): ONLY REQUIRED FOR "polyPhase"
            %           modelType. Each n,m pair is a point at which the
            %           reflectance will be calculated, while the third
            %           dimension k contains polynomial coefficients for
            %           that point in increasing order.
            %
            %   Outputs:
            %       - obj (sonic.LunarLambert): Lunar Lambert object.
            %
            %   References:
            %       -https://articles.adsabs.harvard.edu/pdf/1996LPI....27..841M
            %           McEwen's cubic model
            %       -https://doi.org/10.1111/j.1945-5100.2008.tb00692.x
            %           Gaskell's approximation
            %
            %   Last revised: 09/30/2024
            %   Last author: Priyal Soni

            arguments
                AL          (:,:,:)             double {mustBeInRange(AL, 0, 1)}
                modelType   (1,1)               string {mustBeMember(modelType, {'GaskellApprox', 'McEwenCubic', 'polyPhase'})}
                coeffs      (:,:,:)             double = []
            end

            % Validate model type specified
            switch modelType

                case "GaskellApprox"
                    assert(isempty(coeffs),'sonic:LunarLambert:unexpectedCoeffsInput', ...
                            'Polynomial coefficients are not specified for the Gaskell approximation approach.')

                case "McEwenCubic"
                    assert(isempty(coeffs),'sonic:LunarLambert:unexpectedCoeffsInput', ...
                            'Polynomial coefficients are not specified for the Gaskell approximation approach.')

                case "polyPhase"
                    % coefficients must be provided to use polyPhase
                    assert(~isempty(coeffs),'sonic:LunarLambert:invalidPolyPhaseCoeffs', ...
                            'When using polyPhase for LunarLambert, a column vector of polynomial coefficients must be provided as a third input.') 

                    % check albedo poly coeffs combo
                    polyALcheck = (all(size(AL,1:2) == [1 1]) || all(size(coeffs,1:2) == [1 1])) || (all(size(AL,1:2) == size(coeffs,1:2)));
                    assert(polyALcheck,'sonic:LunarLambert:invalidAlbedoCoeffsCombo', ...
                            'To be a valid input combination, either the albedo or polynomials must have dimension 1 in n and m, or they must be the same size in n and m')
   
                    obj.coeffs = coeffs;

                otherwise
                    error('sonic:LunarLambert:invalidModelType', ...
                        ['Unknown Lunar Lambert input. ' ...
                        '\nPlease specify Lunar Lambert model as one of the following supported model types:' ...
                        '\n\t1. "GaskellApprox"'...
                        '\n\t2. "McEwenCubic"'...
                        '\n\t3. "polyPhase"']);
            end

            % Assign properties
            obj.AL = AL;
            obj.modelType = modelType;

        end

        function r = refl(obj, inputType, var1, var2, var3)
            %% r = refl(obj, input_type, var1, var2, var3)
            %   calculates the reflectance of a surface as described by the
            %   Lunar Lambert model.
            %
            %   Inputs:
            %       - obj (sonic.LunarLambert): The Lunar Lambert reflectance
            %           model object.
            %       - inputType (string):
            %               - "iep" for incidence/emission/phase angles
            %               - "iea" for incidence/emission/azimuth angles
            %               - "ien" for incidence/emission/normal vectors
            %                   -- Must be column vectors, (3xn)
            %       - var1: Must be a double or PointsS2 object
            %       - var2: Must be a double or PointsS2 object
            %       - var3: Must be a double or PointsS2 object
            %
            %   Outputs:
            %       - r (:,: double): matrix of BRDF values
            %
            %   Last revised: 09/30/2024
            %   Last author: Priyal Soni

            arguments
                obj         (1,1)   sonic.Reflectance
                inputType   (1,1)   string {mustBeMember(inputType, {'iep', 'iea', 'ien'})}
                var1        (:,:)       {sonic.Reflectance.mustBeDoubleOrPointsS2(var1)}
                var2        (:,:)       {sonic.Reflectance.mustBeDoubleOrPointsS2(var2)}
                var3        (:,:)       {sonic.Reflectance.mustBeDoubleOrPointsS2(var3)}
            end

            % Input Validation (in Reflectance.m)
            [var1, var2, var3] = obj.checkInputs(inputType, var1, var2, {var3});

            switch inputType
                case "iep"
                    % Check input sizes of albedo and polyPhase (if
                    % applicable)
                    if obj.modelType == "polyPhase" 
                        iepCheck = all(size(var1) == size(obj.AL,1:2)) || all(size(var1) == size(obj.coeffs,1:2)) || (all(size(obj.AL,1:2) == [1 1]) && all(size(obj.coeffs,1:2) == [1 1]));
                        assert(iepCheck,'sonic:LunarLambert:refl:invalidInput', ...
                                        'Please ensure input sizes are valid.')
                    else
                        iepCheck = all(size(var1) == size(obj.AL,1:2)) || all(size(obj.AL,1:2) == [1 1]);
                        assert(iepCheck,'sonic:LunarLambert:refl:invalidInput', ...
                                        'Please ensure input sizes are valid.')
                    end

                    % Assign properties
                    i = var1;
                    e = var2;
                    g = var3;

                    % Calculate cosine i,e for BRDF calculation
                    cosi = cos(i);
                    cose = cos(e);
                        
                case "iea"
                    % Check input sizes of albedo and polyPhase (if
                    % applicable)
                    if obj.modelType == "polyPhase" 
                        ieaCheck = all(size(var1) == size(obj.AL,1:2)) || all(size(var1) == size(obj.coeffs,1:2)) || (all(size(obj.AL,1:2) == [1 1]) && all(size(obj.coeffs,1:2) == [1 1]));
                        assert(ieaCheck,'sonic:LunarLambert:refl:invalidInput', ...
                                        'Please ensure input sizes are valid.')
                    else
                        ieaCheck = all(size(var1) == size(obj.AL,1:2)) || all(size(obj.AL,1:2) == [1 1]);
                        assert(ieaCheck,'sonic:LunarLambert:refl:invalidInput', ...
                                        'Please ensure input sizes are valid.')
                    end
                    
                    % Assign properties
                    i = var1;
                    e = var2;
                    phi = var3;
                    
                    % Calculate cosine i,e for BRDF calculation
                    cosi = cos(i);
                    cose = cos(e);
                    
                    % Calculate phase angles
                    g = obj.calcPhase(i, e, phi);

                case "ien"
                    % Check input sizes of albedo and polyPhase (if
                    % applicable)
                    if obj.modelType == "polyPhase" 
                        ienCheck = all([1, size(var1,2)] == size(obj.AL,1:2)) || all([1, size(var1,2)] == size(obj.coeffs,1:2)) || (all(size(obj.AL,1:2) == [1 1]) && all(size(obj.coeffs,1:2) == [1 1]));
                        assert(ienCheck,'sonic:LunarLambert:refl:invalidInput', ...
                                        'Please ensure input sizes are valid.')
                    else
                        ienCheck = all([1, size(var1,2)] == size(obj.AL,1:2)) || all(size(obj.AL,1:2) == [1 1]);
                        assert(ienCheck,'sonic:LunarLambert:refl:invalidInput', ...
                                        'Please ensure input sizes are valid.')
                    end
                    
                    % Assign properties
                    u_i = var1;
                    u_e = var2;
                    u_n = var3;

                    % Calculate incidence and emission angles for each facet
                    [cosi, cose] = obj.calcCosIE(u_i, u_e, u_n);
                    
                    
                    % Calculate phase angles for each facet
                    g = obj.calcPhase(u_i, u_e);
            end
            
            % Use refl_LL to calculate BRDFs given input parameters
            r = refl_LL(obj, cosi, cose, g);

        end

        function r = refl_LL(obj, cosi, cose, g)
            %% r = refl_LL(obj, cosi, cose, g)
            %   calculates the reflectance of a surface as described by the
            %   lunar-Lambert model.
            %
            %   Inputs:
            %       - obj (sonic.LunarLambert): The Lunar Lambert reflectance
            %           model object.
            %       - cosi: cosine of incidence angles
            %       - cose: cosine of emission angles
            %       - g: phase angles [rad]
            %
            %   Outputs:
            %       - r (:,: double): matrix of BRDF values
            %
            %   Last revised: 09/30/2024
            %   Last author: Priyal Soni

            % Calculate reflectance based on obj.modelType variable
            switch obj.modelType
                case "GaskellApprox"
                    % Define g0 constant as given by Gaskell
                    g0 = pi/3;

                    % Calculate the phase angle-dependent scalar function (beta)
                    L = exp(-g./g0);

                case "McEwenCubic"
                    A = -0.019;
                    B = 0.242e-3;
                    C = -1.46e-6;

                    % Must convert phase angle to degrees to match the
                    % coefficient units in McEwen's paper
                    % Validated against Figure 3 in Gaskell reference
                    gDeg = sonic.Units.RADtoDEG(g);
                    L = 1 + A*gDeg + B*gDeg.^2 + C*gDeg.^3;

                case "polyPhase"
                    % Calculate reflectance using weighted coefficients
                    % provided

                    % Depth of phase polynomial
                    polyCoeffs = obj.coeffs;
                    d = size(polyCoeffs, 3);

                    % Check if given 3D matrix of polynomial coeffs
                    % or one set of polynomial coeffs for all points
                    if d > 1
                        % Calculate L
                        L = polyCoeffs(:,:,1);
                        for g_power = 1:d-1
                            gexp = g.^(g_power);
                            polyCoeffs_i = polyCoeffs(:,:,g_power+1);
                            L = L + gexp.*polyCoeffs_i;
                        end
                    else
                        % obj.coeffs = nxmx1, not phase dependent
                        L = polyCoeffs;
                    end

                otherwise
                    % This shouldn't be thrown because it should error in
                    % the constructor instantiation, but just in case
                    error('sonic:LunarLambert:refl_LL:invalidModelType', ...
                        ['Unknown Lunar Lambert input. ' ...
                        '\nPlease specify Lunar Lambert model as one of the following supported model types:' ...
                        '\n\t1. "GaskellApprox"'...
                        '\n\t2. "McEwenCubic"'...
                        '\n\t3. "polyPhase"']);
            end

            % Calculate albedo to use for reflectance computation
            % Depth of albedo polynomial
            polyAL = obj.AL;
            d = size(polyAL, 3);

            % Check if given 3D matrix of polynomial coeffs or scalar
            if d > 1
                % Calculate A
                A = polyAL(:,:,1);
                for g_power = 1:d-1
                    gexp = g.^(g_power);
                    polyAL_i = polyAL(:,:,g_power+1);
                    A = A + gexp.*polyAL_i;
                end
            else
                % if the third dimension is 1,only scalar albedos were
                % defined
                A = polyAL;
            end

            % Calculate BRDF using the provided phase
            r = (A/pi).*((1-L) + (L.*2)./(cosi+cose));

        end
    end
end

