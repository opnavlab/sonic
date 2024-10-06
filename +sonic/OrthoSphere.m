classdef OrthoSphere < sonic.OrthoRender
    % This software is made available under the MIT license. See SONIC/LICENSE
    % for details.

    properties
        R       (1,1)       double   % radius of sphere [pix]
        g       (1,1)       double   % phase angle [rad]
    end

    properties (SetAccess = private)
        imgSize     (1,1)      double   % image size [pix]
    end

    methods
        function obj = OrthoSphere(R, g)
            %% obj = OrthoSphere(R, g)
            %   Instantiates a sonic OrthoSphere object. 
            % 
            %   OrthoSphere assumes observer is nadir-viewing, and the 
            %   light source is always in the horizontal plane. The phase
            %   angle is angle between the observer and the light source.
            %   
            %   OrthoSphere also renders a sphere with a homogeneous
            %   surface. As a result, all reflectance models used in
            %   conjunction with OrthoSphere should contain scalar model
            %   parameters like albedo (AL). 
            %
            %   Inputs:
            %       - R (1,1) (double): Radius of sphere [pix]. R > 0.
            %       - g (1,1) (double): Phase Angle [rad] between 0 and pi
            %
            %   Outputs:
            %       - obj (sonic.OrthoSphere): OrthoSphere object.
            %
            %   Last revised: 09/30/2024
            %   Last author: Priyal Soni

            arguments
                R       (1,1)       double  {mustBeGreaterThan(R, 0)}
                g       (1,1)       double  {sonic.OrthoSphere.mustBe0ToPi(g)}
            end
            
            % Assign radius and phase angle
            obj.g = g;
            obj.R = R;
            
            % Calculate and assign image size based on radius (1x1)
            obj.imgSize = 2*ceil(R);
            
            % Calculate and assign normal vectors (unit) (3xn) 
            obj.u_n = calcNormals(obj);
            
            % Calculate and assign incidence vectors (unit) (3xn)            
            obj.u_i = calcIncidenceVecs(obj);
            
            % Calculate and assign emission vectors (unit) (3xn)            
            obj.u_e = calcEmissionVecs(obj);
  
            % Calculate and assign emission azimuth angles (1xn) [rad]
            obj.phi = abs(sonic.Reflectance.calcAzimuth("unitVec", obj.u_i, obj.u_e, obj.u_n));

            % Calculate and assign incidence and emission angles (1xn) [rad]
            [obj.i, obj.e] = sonic.Reflectance.calcIE(obj.u_i, obj.u_e, obj.u_n);
        end

        function u_n = calcNormals(obj)
            %% u_n = calcNormals(obj)
            %   Calculates normal vectors for all points on the sphere.
            %
            %   Inputs:
            %       - obj (sonic.Orthosphere) - OrthoSphere class object
            %   Outputs:
            %       - u_n (3xn) (double): unit vectors for all the normals on
            %           the sphere
            %
            %   Last revised: 09/30/2024
            %   Last author: Priyal Soni

            arguments
                obj     sonic.OrthoSphere
            end

            % Span the y and z axes from -R to R
            y = linspace(-obj.R, obj.R, obj.imgSize);
            z = linspace(-obj.R, obj.R, obj.imgSize);

            % Create the y-z plane as a grid
            [y, z] = meshgrid(y, z);

            % Reshape y and z to be vectors for easy calculations
            y = reshape(y, (obj.imgSize)^2, 1);
            z = reshape(z, (obj.imgSize)^2, 1);

            % Calculate x component
            x = sqrt(obj.R^2 - z.^2 - y.^2);

            % If x is imaginary, set it to 0
            % x is only imaginary when R^2 - z^2 - y^2 is negative, which
            % occurs when the point is outside of the sphere
            l = imag(x) == 0;
            x(~l) = NaN;

            % Convert points into vectors and normalize
            n_vec = [x'; y'; z'];
            u_n = n_vec./vecnorm(n_vec, 2, 1);
        end

        function u_i = calcIncidenceVecs(obj)
            %% u_i = calcIncidenceVecs(obj)
            %   Calculates incidence vectors for all points on the sphere.
            %
            %   Inputs:
            %       - obj (sonic.Orthosphere) - OrthoSphere class object
            %   Outputs:
            %       - u_i (3xn) (double): incidence unit vectors for all 
            %           points on the sphere.
            %
            %   Last revised: 09/30/2024
            %   Last author: Priyal Soni

            arguments
                obj     sonic.OrthoSphere
            end
                
            % Use the phase angle input
            u_i = [cos(obj.g); sin(obj.g); 0].*ones(3, (2*obj.R)^2);
        end

        function u_e = calcEmissionVecs(obj)
            %% u_e = calcEmissionVecs(obj)
            %   Calculates emission vectors for all points on the sphere.
            %
            %   Inputs:
            %       - obj (sonic.Orthosphere) - OrthoSphere class object
            %   Outputs:
            %       - u_i (3xn) (double): emission unit vectors for all 
            %           points on the sphere.
            %
            %   Last revised: 09/30/2024
            %   Last author: Priyal Soni

            arguments
                obj     sonic.OrthoSphere
            end
                
            % Always nadir pointing for all points on sphere due to
            % orthographic projection assumption
            u_e = [1; 0; 0].*ones(3, (2*obj.R)^2);
        end


        function [r, sphereMask] = calcRefl(obj, reflObj)
            %% [r, sphereMask] = calcRefl(obj, reflObj)
            %   Calculates BRDF values for all points on the sphere given a
            %   particular reflectance model
            %
            %   Inputs:
            %       - obj (sonic.Orthosphere) - OrthoSphere class object
            %       - reflObj - Any class under sonic.Reflectance
            % 
            %   Outputs:
            %       - r (nxn) (double): BRDF values for all points in image
            %           with n = obj.imgSize. 
            %               Note: Points outside of sphere are assigned -1
            %               in output BRDF matrix
            %       - sphereMask (nxn) (logical): Logical mask to identify
            %               points on sphere with n = obj.imgSize =
            %               2*ceil(R)
            %
            %   Last revised: 09/30/2024
            %   Last author: Priyal Soni

            arguments
                obj         sonic.OrthoSphere
                reflObj     sonic.Reflectance
            end

            % Check type of Reflectance Model
            RM = class(reflObj);
            RM = RM(7:end);

            % Check model specific parameters to ensure scalar inputs
            switch RM
                case "LunarLambert"
                    % Check albedo is scalar or a single polynomial
                    if ~all(size(reflObj.AL, 1:2) == [1, 1])
                        error('sonic:OrthoSphere:invalidInput', ...
                            ['OrthoSphere renders homogeneous spheres.', ...
                            '\nOnly one albedo scalar or polynomial is used for all points on sphere.',...
                            '\nPlease ensure polynomial coefficients are (1, 1, n) with n > 0.']);
                    end
                    
                    % Check coeffs is scalar or a single polynomial
                    if reflObj.modelType == "polyPhase" && all(size(reflObj.coeffs, 1:2) ~= [1, 1])
                        error('sonic:OrthoSphere:invalidInput', ...
                            ['OrthoSphere renders homogeneous spheres.', ...
                            '\nOnly one phase function polynomial is used for all points on sphere.',...
                            '\nPlease ensure polynomial coefficients are (1, 1, n) with n > 0.']);
                    end

                case "OrenNayar"
                    % Check sigma is scalar
                    if ~all(size(reflObj.sigma) == [1, 1])
                        error('sonic:OrthoSphere:invalidInput', ...
                            'Please ensure sigma is a scalar value.');
                    end

                    % Check albedo is scalar
                    if ~all(size(reflObj.AL) == [1, 1])
                        error('sonic:OrthoSphere:invalidInput', ...
                            'Please ensure AL is a scalar value.');
                    end

                case "Hapke"
                    % Check sigma is scalar
                    if ~isempty(reflObj.xi) && ~all(size(reflObj.xi) == [1, 1])
                        error('sonic:OrthoSphere:invalidInput', ...
                            'Please ensure xi is a scalar value.');
                    end

                    % Check albedo is scalar
                    if all(size(reflObj.AL_ss) ~= [1, 1])
                        error('sonic:OrthoSphere:invalidInput', ...
                            'Please ensure AL_ss is a scalar value.');
                    end

                otherwise
                    % Check albedo is a scalar input
                    if isprop(reflObj, 'AL') && ~all(size(reflObj.AL) == [1, 1])
                        error('sonic:OrthoSphere:invalidInput', ...
                            'Please ensure AL is a scalar vsalue.');
                    elseif isprop(reflObj, 'AL_ss') && ~all(size(reflObj.AL_ss) == [1, 1])
                        error('sonic:OrthoSphere:invalidInput', ...
                            'Please ensure AL_ss is a scalar value.');
                    end
            end

            % Check if the angle between source and patch is within +/- pi/2 to ensure
            % the illumination aligns with the orientation
            sphereMask = obj.i <= pi/2 & obj.i >= -pi/2;

            % Generate BRDF matrix
            r = -1*ones(obj.imgSize);

            % Call OrthoRender reflectance function
            r(sphereMask) = obj.renderRefl(reflObj, obj.i(sphereMask), obj.e(sphereMask), obj.phi(sphereMask));

            % Reshape sphereMask for output
            sphereMask = reshape(sphereMask, obj.imgSize, obj.imgSize);
        end
    end

    methods (Static)
        function mustBe0ToPi(g)
            %% mustBe0ToPi(a)
            %   Input validation function for phase angle input
            %
            %   Inputs:
            %       - g: phase angle [rad]
            %
            %   Outputs:
            %       - N/A for input validation functions
            %
            %   Last revised: 09/30/2024
            %   Last author: Priyal Soni

            % Check if phase angle is within [0, pi]
            if g < 0 || g > pi
                error('sonic:Reflectance:invalidInput', ...
                    'Please ensure phase angle (g) is between 0 and pi');
            end
        end
    end
end
