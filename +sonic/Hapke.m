classdef Hapke < sonic.Reflectance
    % This software is made available under the MIT license. See SONIC/LICENSE
    % for details.

    properties
        AL_ss           double      % Albedo - Single Scattering (n x m)    [0, 1]
        hapkeModelType  string      % Type of Hapke model: "IMSA" or "MIMSA"
        HFuncType       string      % Type of H-function approximation: "exact", "rational", or "linear"
        SSPF            string      % Type of Single-Scattering Phase Function: "isotropic", "Rayleigh", or "HG1" (Henyey-Greenstein, one-sided)
        xi              double      % Scattering parameter [0, 1]
    end
    properties (Constant)
        n = 20                      % Number of terms used for summation (if MIMSA option is used to compute P and/or M)  
    end
    properties (SetAccess = private)
        xiFlat                      % Scattering parameter (converted to vector for matrix math)
        AL_ssFlat                   % Albedo - Single Scattering (converted to vector for matrix math)
    end

    methods
        function obj = Hapke(AL_ss, hapkeModelType, HFuncType, SSPF, varargin)
            %% obj = Hapke(AL_ss, HFuncType, SSPF, varargin)
            %   Instantiates a sonic Hapke object.
            %
            %   Inputs:
            %       - AL_ss (nxm double): Single-scattering albedo matrix.
            %       - hapkeModelType (1x1 string): Implementation type of
            %         Hapke model
            %       - HFuncType (1x1 string): Implementation type of 
            %           Chandrasekhar H-Function. Current supported inputs
            %           are:
            %               - "rational"
            %               - "linear"
            %               - "exact"
            %       - SSPF (1x1 string): Implementation type of 
            %           Single-Scattering Phase Function. Current supported 
            %           inputs are:
            %               - "isotropic"
            %               - "Rayleigh"
            %               - "HG1"
            %       - varargin (nx1 double): scattering parameter (xi); 
            %         ONLY REQUIRED FOR "HG1" SSPF
            %
            %   Outputs:
            %       - obj (1x1 sonic.Hapke): Hapke object.
            %
            %   Last revised: 10/4/2024
            %   Last author: Jennifer Nolan
            %
            
            arguments
                AL_ss      (:,:) double {mustBeInRange(AL_ss,0,1)}
                hapkeModelType (1,1) string {mustBeMember(hapkeModelType,{'IMSA','MIMSA'})}
                HFuncType  (1,1) string {mustBeMember(HFuncType,{'exact','rational','linear'})}
                SSPF       (1,1) string {mustBeMember(SSPF,{'isotropic','Rayleigh','HG1'})}
            end
            
            arguments (Repeating)
                varargin (:,:) double {mustBeInRange(varargin,-1,1)}
            end
            
            % Assign inputted values to Hapke object
            obj.AL_ss = AL_ss;
            obj.hapkeModelType = hapkeModelType;
            obj.HFuncType = HFuncType;
            obj.SSPF = SSPF;
            obj.AL_ssFlat = reshape(AL_ss,[1,length(AL_ss(:))]);

            % Check if xi was entered when SSPF is HG1
            if SSPF == 'HG1'
                if ~all(size(varargin) == [1,1])
                    error('sonic:Hapke:invalidInput', ...
                    'Please ensure all inputs are included for SSPF = "HG1" (e.g. xi).');
                else
                    % Check if the size of roughness parameter and albedo match
                    if all(size(varargin{1}) == size(AL_ss)) || all(size(obj.AL_ss) == [1, 1]) || all(size(varargin{1}) == [1, 1])

                   % Assign properties
                    obj.xi = varargin{1};
                    obj.xiFlat = reshape(varargin{1},[1,length(obj.xi(:))]);

                    else
                         error('sonic:Hapke:invalidInput', ...
                    'Please ensure the dimensions of xi match other inputs (e.g., AL_ss).');
                    end
                end
            else
                if  ~isempty(varargin)
                        error('sonic:Hapke:invalidInput', ...
                    ['No xi value is used for SSPF = "isotropic" or SSPF = "Rayleigh".\n' ...
                    'Please specify different SSPF or remove xi input.']);
                end
            end
               

        end
        function r = refl(obj, inputType, var1, var2, var3)
            %% r = refl(obj, inputType, var1, var2, var3)
            %   Creates object containing Hapke values
            %
            %   Inputs:
            %       - obj (1x1 sonic.Hapke): object containing Hapke
            %         parameters
            %       - inputType (1x1 string): defines subsequent set of
            %         inputs.  Current supported inputs are:
            %               - "iep" for incidence/emission/phase angles
            %               - "iea" for incidence/emission/azimuth angles
            %               - "ien" for incidence/emission/normal vectors
            %                   -- Must be column vectors, (3xn)
            %       - var1 
            %       - var2 
            %       - var3
            %
            %   Outputs:
            %       - r (:,:) (double): matrix of Hapke BRDF values 
            %
            %   References:
            %       [1] J. A. Christian, "Spacecraft Optical Navigation"
            % 
            %   Last revised: 10/4/2024
            %   Last author: Jennifer Nolan
            
            arguments
                obj         (1,1) sonic.Hapke
                inputType   (1,1) string {mustBeMember(inputType,["iep","iea","ien"])}
                var1        {sonic.Reflectance.mustBeDoubleOrPointsS2(var1)}
                var2        {sonic.Reflectance.mustBeDoubleOrPointsS2(var2)}
                var3        {sonic.Reflectance.mustBeDoubleOrPointsS2(var3)}
            end
            
            % Input Validation (in Reflectance.m)
            [var1, var2, var3] = obj.checkInputs(inputType, var1, var2, {var3});
                    
            % Verify input sizes for each inputType option
            switch inputType
                case "ien"
                    
                    % Assign inputted variables to specified values
                    u_i = var1;
                    u_e = var2;
                    u_n = var3;

                    % Check if size of albedo input is either 1xn or 1x1
                    if obj.SSPF == 'HG1' && ~(((size(u_i, 2) == size(obj.xi, 2) && size(obj.xi, 1) == 1) || all(size(obj.xi) == [1, 1])))

                        error('sonic:Hapke:invalidInput', ...
                            ['Please ensure input sizes match scattering parameter input, xi.'...
                            '\nFor 3xn vector inputs, xi should be 1xn or 1x1.']);
                    end
                    
                    % Check if size of albedo input is either 1xn or 1x1
                    if ((size(u_i, 2) == size(obj.AL_ss, 2) && size(obj.AL_ss, 1) == 1) || all(size(obj.AL_ss) == [1, 1]))

                        % Calculate i,e, and p angles needed for BRDF calculation
                        [i, e] = obj.calcIE(u_i, u_e, u_n);
                        cosg = obj.calcCosPhase(u_i, u_e);

                    else
                        error('sonic:Hapke:invalidInput', ...
                            ['Please ensure input sizes match albedo input, AL_ss.'...
                            '\nFor 3xn vector inputs, AL_ss should be 1xn or 1x1.']);
                    end

                    outputSize = [1,size(u_i,2)];

                case "iea"

                    % Assign inputted variables to specified values
                    i = var1;
                    e = var2;
                    phi = var3;

                    % Check if size of albedo input is either nxm or 1x1
                    if ~(all(size(var1) == size(obj.AL_ss)) || all(size(obj.AL_ss) == [1, 1]))
                        error('sonic:Hapke:invalidInput', ...
                            ['Please ensure input sizes match albedo input, AL_ss.'...
                            '\nFor nxm inputs, AL_ss should be nxm or 1x1.']);
                    end

                    % Calculate phase angle needed for BRDF calculation
                    cosg = obj.calcCosPhase(i, e, phi);

                    outputSize = size(i);

                case "iep"
                    
                    % Assign inputted variables to specified values
                    i = var1;
                    e = var2;
                    g = var3;
                    
                    % Check if size of albedo input is either nxm or 1x1
                    if ~(all(size(var1) == size(obj.AL_ss)) || all(size(obj.AL_ss) == [1, 1]))
                        error('sonic:Hapke:invalidInput', ...
                            ['Please ensure input sizes match albedo input, AL_ss.'...
                            '\nFor nxm inputs, AL_ss should be nxm or 1x1.']);
                    end

                    cosg = cos(g);

                    outputSize = size(i);                    

            end
            cosg = reshape(cosg,[1,length(cosg(:))]);
            cosi = reshape(cos(i),[1,length(i(:))]);
            cose = reshape(cos(e),[1,length(e(:))]);

            % Define a Lommel-Seeliger object and the Lommel-Seeliger
            % reflectance model for the given incidence angle 
            % and emission angle
            LS = sonic.LommelSeeliger(obj.AL_ss);
            r_LS = LS.refl("iep", i, e);
            r_LSFlat = reshape(r_LS,[1,length(r_LS(:))]);
            
            % Combine the SSPF approximation, the
            % multiple-scattering term, and the Lommel-Seeliger 
            % reflectance model to form the Hapke reflectance model
            P = calcP(obj, cosg);
            M = calcM(obj, cosi, cose);
            r = reshape(r_LSFlat.*(P + M), outputSize);

        end

        function [P] = calcP(obj, cosg)
            %% P = calcSSPF(obj, p)
            %   Compute the selected single-scattering phase function for 
            %   the given phase angle
            %
            %   Inputs:
            %       - obj (1x1 sonic.Hapke): object containing Hapke
            %         parameters
            %       - cosg (nxm double): cosine of phase angles (rad)
            %
            %   Outputs:
            %       - P (nxm double): matrix of values for specified SSPF
            %
            %   References:
            %       [1] L. G. Henyey and J. L. Greenstein, 
            %           "Diffuse radiation in the Galaxy," 
            %           doi: 10.1086/144246.
            %       [2] Lord Rayleigh (Strutt, J. W.), "X. On the 
            %           electromagnetic theory of light,"
            %           doi: 10.1080/14786448108627074
            %       [3] J. A. Christian, "Spacecraft Optical Navigation"
            %
            %   Last revised: 10/4/2024
            %   Last author: Jennifer Nolan
            %
            arguments
                obj     (:,:) sonic.Hapke
                cosg    (:,:) double 
            end            
            
            % Specify which type of SSPF should be used
            switch obj.SSPF
                
                % If specified, use  the isotropic SSPF, which by
                % definition is unity
                case "isotropic"
                    P = ones(size(cosg));               
    
                % If specified, compute the Rayleigh SSPF using Eq. 6.60 in
                % [3] (originally referenced in [2])
                case "Rayleigh"
                    P = (3/4).*(1+(cosg).^2);
    
                % If specified, compute the one-sided Henyey-Greenstein
                % SSPF using Eq. 6.61 in [3] (originally referenced
                % in [1])
                case "HG1"
                    P = (1-obj.xiFlat.^2)./(1 + 2*obj.xiFlat.*cosg + obj.xiFlat.^2).^(3/2);      
            end
        end

        function [M] = calcM(obj, cosi, cose)
            %% M = calcM(obj, cosi, cose)
            %   Compute the multiple-scattering component for the given 
            %   SSPF type at the user-provided i, e, and g values
            %
            %   Inputs:
            %       - obj (1x1 sonic.Hapke): object containing Hapke
            %         parameters
            %       - cosi (nxm double): cosine of incidence angles
            %       - cose (nxm double): cosine of emission angles
            %
            %   Outputs:
            %       - P (nxm double): matrix of values for specified SSPF
            %
            %   References:
            %       [1] L. G. Henyey and J. L. Greenstein, 
            %           "Diffuse radiation in the Galaxy," 
            %           doi: 10.1086/144246.
            %       [2] Lord Rayleigh (Strutt, J. W.), "X. On the 
            %           electromagnetic theory of light,"
            %           doi: 10.1080/14786448108627074
            %       [3] J. A. Christian, "Spacecraft Optical Navigation"
            %
            %   Last revised: 10/4/2024
            %   Last author: Jennifer Nolan
            %
            arguments
                obj   (:,:) sonic.Hapke
                cosi  (:,:) double 
                cose  (:,:) double 
            end

            % The multiple-scattering component for IMSA and MIMSA are
            % identical for "isotropic" and "Raleigh" SSPFs
            if obj.SSPF == "isotropic" || obj.SSPF == "Rayleigh" || obj.hapkeModelType == "IMSA"

                % Determine the Chandrasekhar H-Function for the cos(i)
                % and cos(e) and compute M
                Hcosi=sonic.Math.HFunction(cosi,obj.AL_ssFlat,obj.HFuncType);
                Hcose=sonic.Math.HFunction(cose,obj.AL_ssFlat,obj.HFuncType);
                M = Hcosi.*Hcose-1;

            else

                % Specify vector of all summation terms to be included
                nTerms = 1:1:obj.n;

                % Specify which SSPF type should be used
                switch obj.SSPF
                    
                    % Define the bn coefficients for HG1
                    case "HG1"
                        bn = (2.*nTerms-1).*(-obj.xiFlat).^nTerms;
                end

                % Calculate the MIMSA multiple-scattering component
                % using summation 
                M = MIMSAMultipleScattering(obj, bn, cosi, cose);
            end

        
        end

        function [M] = MIMSAMultipleScattering(obj,bn,cosi,cose)
            %% M = MIMSAMultipleScattering(obj,bn,cosi,cose)
            %   Compute the MIMSA multiple-scattering component for the 
            %   given SSPF type using n-term summation
            %
            %   Inputs:
            %       - obj (1x1 sonic.Hapke): object containing Hapke
            %         parameters
            %       - bn (nxm double): 
            %       - cosi (nxm double): cosine of incidence angles
            %       - cose (nxm double): cosine of emission angles
            %
            %   Outputs:
            %       - P (nxm double): matrix of values for specified SSPF
            %
            %   References:
            %       [1] L. G. Henyey and J. L. Greenstein, 
            %           "Diffuse radiation in the Galaxy," 
            %           doi: 10.1086/144246.
            %       [2] Lord Rayleigh (Strutt, J. W.), "X. On the 
            %           electromagnetic theory of light,"
            %           doi: 10.1080/14786448108627074
            %       [3] J. A. Christian, "Spacecraft Optical Navigation"
            %
            %   Last revised: 10/4/2024
            %   Last author: Jennifer Nolan
            %

            % Specify vector of all summation terms to be included
            nTerms = 1:1:obj.n;

            % Determine the coefficient A_n
            An = zeros(size(nTerms));
            nOdd = nTerms(mod(nTerms,2)==1);
            An(mod(nTerms,2)==1) = (((-1).^((nOdd+1)./2))./nOdd).*...
                                (cumprod(nOdd)./cumprod(nOdd+1));
    
            % Determine the legendre polynomial evaluated at cos(g)
            scriptP_n_cosi = zeros(obj.n,length(cosi));
            scriptP_n_cose = zeros(obj.n,length(cose));
    
            for j=1:obj.n
                legendreCosi = legendre(j,cosi);
                legendreCose = legendre(j,cose);
                scriptP_n_cosi(j,:) = legendreCosi(1,:);
                scriptP_n_cose(j,:) = legendreCose(1,:);
            end
               
            % Determine the scriptP values for both cosi and cose
            scriptP_cosi = 1 + sum(An'.*bn'.*scriptP_n_cosi,1);
            scriptP_cose = 1 + sum(An'.*bn'.*scriptP_n_cose,1);
            
            % Determine the scalar scriptP, dependent solely on An and
            % bn
            scriptP = 1 - sum((An.^2).*bn);
    
            % Determine the Chandrasekhar H-Function for the cos(i)
            % and cos(e)
            Hcosi=sonic.Math.HFunction(cosi,obj.AL_ssFlat,obj.HFuncType);
            Hcose=sonic.Math.HFunction(cose,obj.AL_ssFlat,obj.HFuncType);
    
            % Compute the multiple-scattering component, M
            M=scriptP_cosi.*(Hcose-1)+scriptP_cose.*(Hcosi-1)+scriptP.*(Hcose-1).*(Hcosi-1);
        end
            
    end
end