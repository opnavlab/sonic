classdef (Abstract) OrthoRender
    % This software is made available under the MIT license. See SONIC/LICENSE
    % for details.
    
    properties
        u_n     (3,:)   double      % Normal vectors (unit)
        u_i     (3,:)   double      % Incidence vectors (unit)
        u_e     (3,:)   double      % Emission vectors (unit)
        i       (:,:)   double      % Incidence angles [rad]
        e       (:,:)   double      % Emission angles [rad]
        phi     (:,:)   double      % Emission azimuth angles [rad]
    end

    methods
        function r = renderRefl(obj, reflObj, i, e, phi)
            %% r = renderRefl(obj, reflObj, i, e, a)
            %   Calculates BRDF values for all points on the render given a
            %   particular reflectance model
            % 
            %   Current Rendering Options Supported: OrthoSphere
            %       - More coming soon.
            %
            %   Inputs:
            %       - obj (sonic.Orthosphere) - OrthoSphere class object
            %       - reflObj - Any class under sonic.Reflectance
            %       - i (nxm) (double) - incidence angles [rad]
            %       - e (nxm) (double) - emission angles [rad]
            %       - phi (nxm) (double) - emission azimuth angles [rad]
            %
            %   Outputs:
            %       - r (nxm) (double): BRDF values for all points given
            %
            %   Last revised: 09/30/2024
            %   Last author: Priyal Soni

            arguments
                obj         sonic.OrthoSphere
                reflObj     sonic.Reflectance
                i           (:,:) double
                e           (:,:) double
                phi         (:,:) double
            end

            % Call the reflectance object's refl function to calculate BRDF
            r = reflObj.refl("iea", i, e, phi);
        end
    end
end