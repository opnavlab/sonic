classdef Planes3 < sonic.GeometryP3
    
    properties     
        n               (1, 1)      uint64   
        inf_planes      (1, :)      logical
        has_inf_planes  (1, 1)      logical

        p3              (4, :)      double
        dist            (1, :)      double
        normal          (3, :)      double
    end

    methods

        function obj = Planes3(varargin)
        %% obj = Planes3(varargin)
        %   Constructs a Plane object in P3/R3. This object can contain a
        %   multitude of planes, including planes at infinity. 
        %   
        %   Inputs:
        %       - coords (4xn double): A set of coordinates specifying 
        %         `n` planes in P3
        %       
        %           OR
        %
        %       - normal (3xn double): Normal vector for each plane,
        %         with the normals as the columns of this matrix
        %       - dist (1xn double): Signed distance from origin for `n`
        %         planes
        %   Outputs:
        %       - obj (sonic.Planes3): Planes3 object, which can contain
        %         multiple planes. 
        %
        %   Last revised: 5/06/24
        %   Last author: Michael Krause

            switch nargin
                case 1

                    % Pull out the raw inputs:
                    raw_coords = varargin{1};

                    % Grab the size of the input.
                    in_size = size(raw_coords);
            
                    % Number of columns = number of planes.
                    obj.n = in_size(2);
                    % Number of rows should be 4 (for P3 representation).
                    num_rows = in_size(1);
                    
                    % Ensure we're given 4 rows:
                    if num_rows ~= 4
                        error('sonic:Planes3:malformedVector', ...
                            ['If parameterizing a plane via four ' ...
                            'coordinates, the input must be a ' ...
                            '4-by-n matrix.']);
                    end

                    % Normalize the P3 representation:
                    obj.p3 = raw_coords./vecnorm(raw_coords);

                    % Consider if we're encountering an infinite plane by
                    % comparing the magnitudes of the distance to the
                    % normal vector:
                    dist_ratios = vecnorm(obj.p3(1:3, :))./abs(obj.p3(4, :));
                    obj.inf_planes = dist_ratios < sonic.Tolerances.InfPlaneDistRatio;
                    obj.has_inf_planes = any(obj.inf_planes);

                    % Save the normals and distances:
                    obj.normal = obj.p3(1:3, :);
                    obj.dist = -obj.p3(4, :);

                    % And then overwrite accordingly for the planes at
                    % infinity:
                    if obj.has_inf_planes
                        obj.normal(:, obj.inf_planes) = [0;0;0];
                        obj.dist(obj.inf_planes) = Inf;
                        obj.p3(:, obj.inf_planes) = [0;0;0;1];
                    end
                    
                case 2
                    
                    norm_vec = varargin{1};
                    dist = varargin{2};

                    norm_vec_size = size(norm_vec);
                    dist_size = size(dist);

                    % Ensure normal is 3xn matrix:
                    if norm_vec_size(1) ~= 3
                        error('sonic:Planes3:malformedNormalVec', ...
                            ['If parameterizing a plane via a normal ' ...
                            'and distance, must specify the normal vector ' ...
                            'then distance as arguments. Normal vector ' ...
                            'must be a 3xn matrix.']);
                    end

                    % Ensure distnace is 1xn vector:
                    if dist_size(1) ~= 1
                        error('sonic:Planes3:malformedDistance', ...
                            ['If parameterizing a plane via a normal ' ...
                            'and distance, must specify the normal vector ' ...
                            'then distance as arguments. Distance ' ...
                            'must be an 1xn vector.']);
                    end

                    % Ensure # of columns agrees:
                    if norm_vec_size(2) ~= dist_size(2)
                        error('sonic:Planes3:mismatchedDims', ...
                            ['When specifying the normal vectors ' ...
                            'and distances as arguments, the numbers '...
                            'columns must agree.']);
                    end
                    
                    % Save off number of planes:
                    obj.n = norm_vec_size(2);

                    % Check for planes at infinity, as specified by the
                    % paralle distance being Inf:
                    obj.inf_planes = isinf(dist);
                    obj.has_inf_planes = any(obj.inf_planes);
                
                    % Save the P3 repr, normals and distances:
                    obj.normal = norm_vec;
                    obj.dist = dist;
                    raw_p3 = [norm_vec; -dist];
                    obj.p3 = raw_p3./vecnorm(raw_p3);

                    % And then overwrite accordingly for the planes at
                    % infinity:
                    if obj.has_inf_planes
                        obj.normal(:, obj.inf_planes) = [0;0;0];
                        obj.p3(:, obj.inf_planes) = [0;0;0;1];
                    end

                otherwise
                    error('sonic:Planes3:IncorrectParameterization', ...
                        ['Must input a 4-element P3 description of a ' ...
                        'plane, or a normal/distance pair.']);
            end
                
        end

    end
end

