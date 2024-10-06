classdef OrenNayar < sonic.Reflectance
    % This software is made available under the MIT license. See SONIC/LICENSE
    % for details.

    properties
        AL                  double      % Albedo (n x m)
        sigma               double      % Roughness (n x m) [rad]
        calcInterrefl       logical     % T/F Calculate interreflectance
    end

    methods
        function obj = OrenNayar(AL, sigma, calcInterrefl)
            %% obj = OrenNayar(AL, sigma, calcInterrefl)
            %   Instantiates a sonic OrenNayar object.
            %
            %   Inputs:
            %       - AL (n,m): Albedo matrix between 0 and 1
            %           Note: can be scalar.
            %       - sigma (n,m): Surface roughness matrix 
            %           Note: can be scalar.
            %       - calcInterrefl: Logical (0 = false (default), 1 = true)
            %           for whether the user wishes to include interreflectance
            %           in their modeling.
            %
            %   Outputs:
            %       - obj (sonic.OrenNayar): Oren-Nayar reflectance model object.
            %
            %   References:
            %       - "Generalization of the Lambertian model..."
            %           https://doi.org/10.1007/BF01679684
            %
            %   Last revised: 10/02/2024
            %   Last author: Priyal Soni

            arguments
                AL              (:,:)     double {mustBeInRange(AL, 0, 1)}
                sigma           (:,:)     double {mustBePositive(sigma)}
                calcInterrefl   (1,1)     {sonic.OrenNayar.checkInterreflInput(calcInterrefl)} = 0
            end

            % Check if the size of roughness parameter and albedo match
            if all(size(sigma) == size(AL)) || all(size(AL) == [1, 1]) || all(size(sigma) == [1, 1])
                % Assign properties
                obj.AL = AL;
                obj.sigma = sigma;
            else
                error('sonic:OrenNayar:invalidInput', ...
                    'Please ensure all model parameter dimensions match');
            end

            % Assign the calcInterrefl option
            obj.calcInterrefl = calcInterrefl;

        end

        function r = refl(obj, inputType, var1, var2, var3)
            %% r = refl(obj, input_type, var1, var2, var3)
            %   calculates the BRDF of a surface as described by the
            %   Oren-Nayar model.
            %
            %   Inputs:
            %       - obj (sonic.OrenNayar): Oren-Nayar reflectance
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
            %   Last revised: 10/02/2024
            %   Last author: Priyal Soni

            arguments
                obj         (1,1)       sonic.Reflectance
                inputType   (1,1)       string {mustBeMember(inputType, {'iep', 'iea', 'ien'})}
                var1        (:,:)       {sonic.Reflectance.mustBeDoubleOrPointsS2(var1)}
                var2        (:,:)       {sonic.Reflectance.mustBeDoubleOrPointsS2(var2)}
                var3        (:,:)       {sonic.Reflectance.mustBeDoubleOrPointsS2(var3)}
            end

            % Input Validation (in Reflectance.m)
            [var1, var2, var3] = obj.checkInputs(inputType, var1, var2, {var3});

            % Albedo and roughness parameter input validation
            switch inputType
                case "iep"
                    % Check sizes of albedo/roughness match var1 (i) size
                    if (all(size(var1) == size(obj.AL)) || size(obj.AL, 2) == 1)
                        % Assign properties
                        i = var1;
                        e = var2;
                        g = var3;
                    else
                        error('sonic:OrenNayar:invalidInput', ...
                            'Please ensure input sizes match reflectance model parameters (e.g. AL)');
                    end

                    % Calculate emission azimuth angles from the
                    % plane of incidence
                    cosPhi = obj.calcCosAzimuth("phaseAng", i, e, g);

                case "iea"
                    % Check sizes of albedo/roughness match var1 (i) size
                    if (all(size(var1) == size(obj.AL)) || size(obj.AL, 2) == 1)
                        % Assign properties
                        i = var1;
                        e = var2;
                        phi = var3;
                    else
                        error('sonic:OrenNayar:invalidInput', ...
                            'Please ensure input sizes match reflectance model parameters (e.g. AL)');
                    end

                    % Calculate cos(phi) for the BRDF calculation
                    cosPhi = cos(phi);

                case "ien"

                    % Check if size of albedo input is either 1xn or 1x1
                    if ~(((size(var1, 2) == size(obj.sigma, 2) && size(obj.sigma, 1) == 1) || all(size(obj.sigma) == [1, 1])))

                        error('sonic:Hapke:invalidInput', ...
                            ['Please ensure input sizes match roughness parameter input, sigma.'...
                            '\nFor 3xn vector inputs, sigma should be 1xn or 1x1.']);
                    end

                    % Check sizes of albedo/roughness match var1 (u_i) size
                    if (size(var1, 2) == size(obj.AL, 2) || size(obj.AL, 2) == 1)
                        % Assign properties
                        u_i = var1;
                        u_e = var2;
                        u_n = var3;
                    else
                        error('sonic:OrenNayar:invalidInput', ...
                            'Please ensure input sizes match reflectance model parameters (e.g. AL)');
                    end

                    % Calculate incidence and emission angles for each facet
                    [i, e] = obj.calcIE(u_i, u_e, u_n);

                    % Calculate emission azimuth angles from the
                    % plane of incidence
                    cosPhi = abs(obj.calcCosAzimuth("unitVec", u_i, u_e, u_n));
            end

            % Call refl_ON to calculate the BRDF given i, e, a angles
            r = obj.refl_ON(i, e, cosPhi);

        end


        function r = refl_ON(obj, i, e, cosPhi)
            %% r =  refl_ON(obj, i, e, cosPhi)
            %   Calculates the reflectance of a surface as described by the
            %   Oren-Nayar model. This function uses the Oren-Nayar
            %   approximation rather than directly evaluating the model.
            %
            %   Inputs:
            %       - obj (sonic.OrenNayar): The Oren-Nayar reflectance
            %           model object.
            %       - i: incidence angles
            %       - e: emission angles
            %       - cosPhi: cosine of emission vector azimuth angle as
            %           measured from the plane of incidence
            %
            %   Outputs:
            %       - r (:,:) (double): matrix of BRDF values
            %
            %   Last revised: 10/01/2024
            %   Last author: Priyal Soni

            arguments
                obj         (1,1)   sonic.Reflectance
                i           (:,:)   double
                e           (:,:)   double
                cosPhi      (:,:)   double
            end

            % Determine alpha and beta from incidence and emission angles
            max_ie = max(i, e);
            min_ie = min(i, e);

            % Calculate the constants for the model fit
            C1 = 1 - 0.5.*obj.sigma.^2./(obj.sigma.^2 + 0.33);

            C3 = 0.125*(obj.sigma.^2./(obj.sigma.^2 + 0.09)).*((4.*max_ie.*min_ie)./pi^2).^2;

            % C2 constant depends on cosPhi
            l = cosPhi >= 0;
            C2 = 0.45*(obj.sigma.^2./(obj.sigma.^2 + 0.09)).*ones(size(i));
            C2(l) = C2(l).*sin(max_ie(l));
            C2(~l) = C2(~l).*(sin(max_ie(~l)) - (2*min_ie(~l)/pi).^3);

            % Calculate reflectance due to direct illumination
            Lr1 = (obj.AL./pi).*(C1 + cosPhi.*C2.*tan(min_ie) + ...
                (1 - abs(cosPhi)).*C3.*tan((max_ie+min_ie)./2));

            % Calculate reflectance due to interreflections
            if obj.calcInterrefl == 1
                Lr2 = 0.17*(obj.AL.^2/pi).*(obj.sigma.^2./(obj.sigma.^2 + 0.13)).* ...
                    (1 - cosPhi.*(2*min_ie./pi).^2);
            else
                Lr2 = 0;
            end
            
            % Calculate total reflectance
            r = Lr1 + Lr2;
        end
    end

    methods (Static)
        function checkInterreflInput(calcInterrefl)
            %% checkInterreflInput(calcInterrefl)
            %   Checks if input is either logical or 0 or 1 for input
            %   validation.
            %
            %   Inputs:
            %       - calcInterrefl: Logical (0 = false (default), 1 = true)
            %           for whether the user wishes to include interreflectance
            %           in their modeling.
            % 
            %   Outputs:
            %       - N/A for input validation functions
            % 
            %   Last revised: 10/02/2024
            %   Last author: Priyal Soni

            if ~(class(calcInterrefl) == "logical" || calcInterrefl == 0 || calcInterrefl == 1)
                error('sonic:OrenNayar:invalidInput', ...
                    'Please ensure calcInterrefl is either a logical or 0 or 1.');
            end
        end
    end
end
