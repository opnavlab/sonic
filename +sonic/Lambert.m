classdef Lambert < sonic.Reflectance

    properties
        AL          double    % Albedo
    end

    methods

        function obj = Lambert(AL)
            %% obj = Lambert(AL)
            %   Instantiates a sonic Lambert object.
            %
            %   Inputs:
            %       - AL (n,m): Albedo matrix between 0 and 1
            %           Note: can be scalar.
            %
            %   Outputs:
            %       - obj (sonic.Lambert): Lambert reflectance model object
            %
            %   References:
            %       - "Photometria sive de mensura et gradibus luminis..."
            %           Jean-Henri Lambert
            %       - Spacecraft Optical Navigation, forthcoming textbook
            %           by John Christian
            %
            %   Last revised: 09/30/2024
            %   Last author: Priyal Soni
            arguments
                AL          (:,:)     double {mustBeInRange(AL, 0, 1)}
            end

            obj.AL = AL;
        end

        function r = refl(obj, inputType, var1, var2, varargin)
            %% r = refl(obj, input_type, var1, var2, varargin)
            %   calculates the BRDF of a surface as described by the
            %   Lambert reflectance model.
            %
            %   Inputs:
            %       - obj (sonic.Lambert): Lambert reflectance
            %           model object.
            %       - inputType (string):
            %               - "iep" for incidence/emission/phase angles
            %               - "iea" for incidence/emission/azimuth angles
            %               - "ien" for incidence/emission/normal vectors
            %                   -- Must be column vectors, (3xn)            
            %       - var1: Must be a double or PointsS2 object
            %       - var2: Must be a double or PointsS2 object
            %       - varargin: - REQUIRED var3 input for "ien"
            %                   - optional var3 input for "iea" and "iep"
            %                   - Must be a double or PointsS2 object
            %
            %   Outputs:
            %       - r (:,:) (double): matrix of BRDF values
            %
            %   Last revised: 10/02/2024
            %   Last author: Priyal Soni

            arguments
                obj         (1,1)   sonic.Reflectance
                inputType   (1,1)   string {mustBeMember(inputType, {'iep', 'iea', 'ien'})}
                var1        (:,:)   {sonic.Reflectance.mustBeDoubleOrPointsS2(var1)}
                var2        (:,:)   {sonic.Reflectance.mustBeDoubleOrPointsS2(var2)}
            end

            arguments (Repeating)
                % var3 (required for inputType = "ien")
                varargin     % Must be Double or PointsS2, but comes packaged in a cell due to varargin           
            end
            
            switch inputType
                case "ien"
                    % Input Validation (in Reflectance.m)
                    [var1, var2, var3] = obj.checkInputs(inputType, var1, var2, varargin);
                    
                    % Assign Inputs
                    u_i = var1;
                    u_e = var2;
                    u_n = var3;
                    
                    % Check if size of albedo input is either 1xn or 1x1
                    if ((size(u_i, 2) == size(obj.AL, 2) && size(obj.AL, 1) == 1) || all(size(obj.AL) == [1, 1]))
                        
                        % Calculate i,e angles needed for BRDF calculation
                        [i, e] = obj.calcIE(u_i, u_e, u_n);

                    else
                        error('sonic:Lambert:invalidInput', ...
                            ['Please ensure input sizes match albedo input, AL.'...
                            '\nFor 3xn vector inputs, AL should be 1xn or 1x1.']);
                    end

                otherwise % "iea" and "iep"
                    % Input Validation (in Reflectance.m)
                    [var1, var2] = obj.checkInputs(inputType, var1, var2, varargin);
                    
                    % Assign Inputs
                    i = var1; 
                    e = var2;
                    
                    % Check if size of albedo input is either nxm or 1x1
                    if ~(all(size(var1) == size(obj.AL)) || all(size(obj.AL) == [1, 1]))
                        error('sonic:Lambert:invalidInput', ...
                            ['Please ensure input sizes match albedo input, AL.'...
                            '\nFor nxm inputs, AL should be nxm or 1x1.']);
                    end
            end
            
            % Calculate the BRDF
            r = (obj.AL./pi).*ones((size(i)));

        end
    end
end

