classdef OrthoPatch < sonic.OrthoRender

    properties
        % grayscale image with height encoded as int b/t 0 and 255, inclusive
        terrain (:, :)

        % albedos for each terrain point between 0 and 1
        albedo (:, :) double {mustBeInRange(albedo, 0, 1)}

        % scale of unit height
        htscale (1, 1) double = 1

        % scale of unit length (in plane)
        lgscale (1, 1) double = 1

        l_i

        % Actual properties of OrthoPatch to be used in rendering
        % x is the number of columns, y is the number of rows in terrain
        facets % height values at facet centers ((y-1), (x-1))
        %normals % normal vectors at facet centers (3, (x-1) * (y-1))
        %slopesX % slope in the x-direction at facet center ((y-1), (x-1))
        %slopesY % slope in the y-direction at facet center ((y-1), (x-1))

    end

    methods

        % Local patch plane frame
        % x defined as positive first dimension in DEM
        % y defined as positive second dimension in DEM
        % z defined as x cross y

        function obj = OrthoPatch(terrain, albedo, l_inc_x, l_inc_y, l_inc_z, htscale, lgscale)
            obj.terrain = double(terrain * htscale);
            obj.albedo = albedo;
            obj.l_i = [l_inc_x; l_inc_y; l_inc_z];
            obj.htscale = htscale;
            obj.lgscale = lgscale;

            [obj.facets, obj.u_n] = calcFacets(obj);

            %obj.normals = calcNormals(obj);
            obj.u_e = calcEmissionVectors(obj);
            obj.u_i = calcIncidenceVectors(obj);
            obj.a = sonic.Reflectance.calcAzimuth(obj.u_i, obj.u_e, obj.u_n);
            [obj.i, obj.e] = sonic.Reflectance.calcIE(obj.u_i, obj.u_e, obj.u_n);

            %changeMatrix = [0, 0, 1; 1, 0, 0; 0, 1, 0];
            %obj.u_e = changeMatrix * obj.u_e;
            %obj.u_i = changeMatrix * obj.u_i;
            %obj.normals = changeMatrix * obj.normals;

            %obj.shadows = applyShadows(obj);
        end

        function [facets, norms] = calcFacets(obj)
            rows = size(obj.terrain, 1);
            cols = size(obj.terrain, 2);
            facets = zeros([rows - 1, cols - 1]);
            % norms = zeros(3, numel(facets));

            for j = 1:size(facets, 1)
                for k = 1:size(facets, 2)
                    % vec1 = [((j+1) - j) * obj.lgscale; ((k+1) - k) * obj.lgscale; (obj.terrain(j+1, k+1) - obj.terrain(j, k))];
                    % vec2 = [((j+1) - j) * obj.lgscale; (k - (k+1)) * obj.lgscale; (obj.terrain(j+1, k) - obj.terrain(j, k+1))];
                    % normal = cross(vec2, vec1);
                    % normal = normal/norm(normal);
                    facets(j, k) = mean([obj.terrain(j, k), obj.terrain(j+1, k), obj.terrain(j, k+1), obj.terrain(j+1, k+1)]);
                    % elem = k * (size(facets, 1) - 1) + j;
                    % norms(1, elem) = normal(1);
                    % norms(2, elem) = normal(2);
                    % norms(3, elem) = normal(3);
                end
            end

            norm_sobelY = (1/8) * [1, 0, -1; 2, 0, -2; 1, 0, -1];
            norm_sobelX = transpose(norm_sobelY);

            convX = conv2(facets, norm_sobelX, "same");
            convY = conv2(facets, norm_sobelY, "same");

            slopesX = (1 / obj.lgscale) .* convX;
            slopesY = (1 / obj.lgscale) .* convY;

            norms = zeros(3, numel(facets));

            for m = 1:numel(facets)
                norm_temp = [-slopesX(m), -slopesY(m), 1];
                norm_temp = norm_temp / norm(norm_temp);
                norms(1, m) = norm_temp(1);
                norms(2, m) = norm_temp(2);
                norms(3, m) = norm_temp(3);
                %norms(3, m) = 1;
            end


        end

        function result = calcEmissionVectors(obj)
            result = [0; 0; 1] .* ones(3, size(obj.facets, 1) * size(obj.facets, 2));
        end

        function result = calcIncidenceVectors(obj)
            result = obj.l_i .* ones(3, size(obj.facets, 1) * size(obj.facets, 2));
        end

        function r = calcRefl(obj, reflObj)

            n = length(obj.i);

            % Check type of Reflectance Model
            RM = class(reflObj);
            RM = RM(7:end);
            
            if isprop(reflObj, 'AL') && size(reflObj.AL, 2) > 1
                reflObj.AL = reshape(reflObj.AL, [1, n]);
            elseif isprop(reflObj, 'AL_ss') && size(reflObj.AL_ss, 2) > 1
                reflObj.AL_ss = reshape(reflObj.AL_ss, [1, n]);
            end

            % this is gonna change later lol
            switch RM
                case "LunarLambert"
                    switch reflObj.modelType
                        case "McEwenApprox"
                            if size(reflObj.g0, 2) > 1
                                reflObj.g0 = reshape(reflObj.g0, [1, n]);
                            end
                        case "McEwenCubic"
                            if size(reflObj.A, 2) > 1
                                reflObj.A = reshape(reflObj.A, [1, n]);
                                reflObj.B = reshape(reflObj.B, [1, n]);
                                reflObj.C = reshape(reflObj.C, [1, n]);
                            end
                        case "coeff"
                            if size(reflObj.A, 2) > 1
                                reflObj.A = reshape(reflObj.A, [1, n]);
                                reflObj.B = reshape(reflObj.B, [1, n]);
                            end
                    end
                case "Oren-Nayar"
                    if size(reflObj.sigma, 2) > 1
                        reflObj.sigma = reshape(reflObj.sigma, [1, n]);
                    end
            end

            % Call OrthoRender reflectance function
            r = renderRefl(obj, reflObj);

            % Set reflectance to 0 if incidence angle is greater than 90
            % degrees
            r(abs(obj.i) > pi/2) = 0;

            % reshape reflectance to be square matrix
            r = reshape(r, size(obj.facets, 1), size(obj.facets, 2));

        end

        function showPatch(obj)
            figure
            title('Facets and normals');
            hold on;
            axis equal
            %surf(obj.facets); % Plot the surface
            [ycoords, xcoords] = meshgrid(1:size(obj.facets, 1), 1:size(obj.facets, 2));

            xcoords = xcoords * obj.lgscale;
            ycoords = ycoords * obj.lgscale;

            surf(xcoords, ycoords, obj.facets);

            xnorms = obj.u_n(1, :)';
            ynorms = obj.u_n(2, :)';
            znorms = obj.u_n(3, :)';

            xcoords = reshape(xcoords, [], 1);
            ycoords = reshape(ycoords, [], 1);
            zcoords = reshape(obj.facets, [], 1);

            quiver3(xcoords, ycoords, zcoords, xnorms, ynorms, znorms, 'Color', 'r');
            hold off;

            % figure
            % title('original terrain');
            % surf(obj.terrain); % Plot the surface
            % hold on;
            % axis equal
            % hold off;

            % figure
            % arr = reshape(obj.shadows, size(obj.facets, 1), size(obj.facets, 2));
            % imshow(arr);
        end

    end
    
    methods (Static)
        function showRefl(r, varargin)
            figure();
            imshow(r(2:size(r, 1)-1, 2:size(r, 2)-1))

            if size(varargin, 2) == 1
                title(varargin{1});
            end
                
        end
    end

end