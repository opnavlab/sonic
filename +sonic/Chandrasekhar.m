classdef Chandrasekhar < sonic.Reflectance
    % This software is made available under the MIT license. See SONIC/LICENSE
    % for details.

    properties
        AL_ss           double      % Albedo - Single Scattering (n x m)    [0, 1]
        HFuncType       string      % Type of H-function approximation: Either "exact", "rational", or "linear"
    end

    methods
        function obj = Chandrasekhar(AL_ss, HFuncType)
            %% obj = Chandrasekhar(AL_ss, HFuncType)
            %   Instantiates a sonic Chandrasekhar object.
            %
            %   Inputs:
            %       - AL_ss (n,m): Single-scattering albedo matrix between
            %           0 and 1. Note: can be scalar.
            %       - HFuncType: (1x1 string): Implementation type of
            %           Chandrasekhar H-Function. Current supported inputs
            %           are:
            %               - "rational"
            %               - "linear"
            %               - "exact"
            %
            %   Outputs:
            %       - obj (sonic.Chandrasekhar): Chandrasekhar reflectance
            %           model object
            %
            %   References:
            %       - Radiative Transfer, S. Chandrasekhar
            %       - Spacecraft Optical Navigation, forthcoming textbook
            %           by John Christian
            %
            %   Last revised: 09/30/2024
            %   Last author: Priyal Soni

            arguments
                AL_ss       (:,:)    double {mustBeInRange(AL_ss,0,1)}
                HFuncType   (1,1)    string {mustBeMember(HFuncType,{'exact','rational','linear'})}
            end
            
            % Assign properties
            obj.AL_ss = AL_ss;
            obj.HFuncType = HFuncType;
        end

        function r = refl(obj, inputType, var1, var2, varargin)
            %% r = refl(obj, input_type, var1, var2, varargin)
            %   calculates the BRDF of a surface as described by the
            %   Chandrasekhar reflectance model.
            %
            %   Inputs:
            %       - obj (sonic.Chandrasekhar): Chandrasekhar reflectance
            %           model object.
            %       - inputType (string):
            %               - "iep" for incidence/emission/phase angles
            %               - "iea" for incidence/emission/azimuth angles
            %               - "ien" for incidence/emission/normal vectors
            %                   -- Must be column vectors, (3xn)            
            %       - var1: Must be a double or PointsS2 object
            %       - var2: Must be a double or PointsS2 object
            %       - varargin: REQUIRED var3 input for "ien"
            %                   optional var3 input for "iea" and "iep"
            %
            %   Outputs:
            %       - r (:,:) (double): matrix of BRDF values
            %
            %   Last revised: 10/02/2024
            %   Last author: Priyal Soni

            arguments
                obj         (1,1)   sonic.Reflectance
                inputType   (1,1)   string {mustBeMember(inputType, {'iep', 'iea', 'ien'})}
                var1                {sonic.Reflectance.mustBeDoubleOrPointsS2(var1)}
                var2                {sonic.Reflectance.mustBeDoubleOrPointsS2(var2)}
            end

            arguments (Repeating)
                % var3 (required for inputType = "ien")
                varargin     % Must be Double or PointsS2, but comes packaged in a cell due to varargin
            end
            
            % Albedo parameter input validation
            switch inputType
                case "ien"
                    % Input Validation (in Reflectance.m)
                    [var1, var2, var3] = obj.checkInputs(inputType, var1, var2, varargin);
                    
                    % Assign properties
                    u_i = var1;
                    u_e = var2;
                    u_n = var3;

                    % Check if size of albedo input is either 1xn or 1x1
                    if ((size(u_i, 2) == size(obj.AL_ss, 2) && size(obj.AL_ss, 1) == 1) || all(size(obj.AL_ss) == [1, 1]))

                        % Calculate i,e angles needed for BRDF calculation
                        [i, e] = obj.calcIE(u_i, u_e, u_n);

                    else
                        error('sonic:Chandrasekhar:invalidInput', ...
                            ['Please ensure input sizes match albedo input, AL_ss.'...
                            '\nFor 3xn vector inputs, AL_ss should be 1xn or 1x1.']);
                    end

                otherwise % "iea" and "iep"
                    % Input Validation (in Reflectance.m)
                    [var1, var2] = obj.checkInputs(inputType, var1, var2, varargin);
                    
                    % Assign properties
                    i = var1;
                    e = var2;

                    % Check if size of albedo input is either nxm or 1x1
                    if ~(all(size(var1) == size(obj.AL_ss)) || all(size(obj.AL_ss) == [1, 1]))
                        error('sonic:Chandrasekhar:invalidInput', ...
                            ['Please ensure input sizes match albedo input, AL_ss.'...
                            '\nFor nxm inputs, AL_ss should be nxm or 1x1.']);
                    end
            end
            
            % Calculate Lommel-Seeliger reflectance
            LS = sonic.LommelSeeliger(obj.AL_ss);
            r_LS = LS.refl("iep", i, e);
            
            % Calculate H(cos(i))*H(cos(e))
            Hcosi=sonic.Math.HFunction(cos(i), obj.AL_ss, obj.HFuncType);
            Hcose=sonic.Math.HFunction(cos(e), obj.AL_ss, obj.HFuncType);
            HcosiHcose = Hcosi.*Hcose;
        
            % Calculate the BRDF
            r = r_LS.*HcosiHcose;
        end
    end
end