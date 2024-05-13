classdef Points3 < sonic.GeometryP3
    
    % This can contain a multitude of points. 
    properties
        n                   (1, 1) uint64
        p3                  (4, :) double
        inf_points          (1, :) logical
        has_inf_points      (1, 1) logical
    end
    
    properties (Dependent)
        r3                  (3, :) double
        r3_finite           (3, :) double
        p3_finite           (4, :) double
        p3_infinite         (4, :) double
    end

    methods

        function obj = Points3(pts)
        %% obj = Points3(pts)
        %   Constructs a Points object in P3/R3. This object can contain a
        %   multitude of points, including points at infinity. 
        %   
        %   Inputs:
        %       - pts: A set of points, either in P3 or R3. The source of
        %         the points is inferred from the dimension of the input
        %         points:
        %           - (3xn double): points are in R3
        %           - (4xn double): points are in P3
        %
        %   Outputs:
        %       - obj (sonic.Points3): Points3 object, a container for
        %         a multitude of 3D points.
        %
        %   Last revised: 2/15/24
        %   Last author: Michael Krause
            

            % Grab the size of the input.
            in_size = size(pts);
    
            % Number of columns = number of points. Good to know.
            obj.n = in_size(2);
            % Number of rows determines if P3 or R3.
            in_dim = in_size(1);

            switch in_dim
                case 3  % Points are in R3
                    % Check for points at infinity:
                    inf_pts = any(pts == Inf, 1);

                    % Save off points and metadata, noting the infinite
                    % points accordingly:
                    obj.p3 = [pts; ones(1, obj.n)];
                    obj.p3(4, inf_pts) = 0;
                    obj.inf_points = inf_pts;
                    obj.has_inf_points = any(inf_pts);
                case 4  % Points are in P3

                    % Check for points at infinity and points with negative
                    % homogenous component:
                    inf_pts = abs(pts(4, :)) < sonic.Tolerances.HomNorm;

                    % Normalize the infinite points accordingly:
                    pts(:, ~inf_pts) = pts(:, ~inf_pts)./pts(4, ~inf_pts);
                    pts(4, inf_pts) = 0;

                    % Save off points and metadata:
                    obj.p3 = pts;
                    obj.inf_points = inf_pts;
                    obj.has_inf_points = any(inf_pts);
                otherwise
                    error('sonic:Points3:IncorrectDimension', ...
                        'Must input a 3 or 4-by-n matrix of 3D points.');
            end
                
        end

        function val = get.r3(obj)
        %% val = get.r3(obj)
        %   Returns the R3 representation of all points. Points at infinity
        %   have the value infinity for all coordinates. 
        %   
        %   Inputs:
        %       - obj (sonic.Points3): A Points3 object with n points
        %
        %   Outputs:
        %       - val (3xn): The R3 realization of these points
        %
        %   Last revised: 2/15/24
        %   Last author: Michael Krause

            val = obj.p3(1:3, :);
            val(:, obj.inf_points) = Inf;
        end

        function val = get.r3_finite(obj)
        %% val = get.r3_finite(obj)
        %   Returns the R3 representation of all finite points. Points at 
        %   infinity are excluded.
        %   
        %   Inputs:
        %       - obj (sonic.Points3): A Points3 object with m finite
        %         points
        %
        %   Outputs:
        %       - val (3xm): The R3 realization of all finite points
        %
        %   Last revised: 2/15/24
        %   Last author: Michael Krause

            val = obj.p3(1:3, ~obj.inf_points);
        end

        function val = get.p3_finite(obj)
        %% val = get.p3_finite(obj)
        %   Returns the P3 representation of all finite points. Points at 
        %   infinity are excluded.
        %   
        %   Inputs:
        %       - obj (sonic.Points3): A Points3 object with m finite
        %         points
        %
        %   Outputs:
        %       - val (4xm): The P3 realization of all finite points
        %
        %   Last revised: 2/15/24
        %   Last author: Michael Krause


            val = obj.p3(:, ~obj.inf_points);
        end

        function val = get.p3_infinite(obj)
        %% val = get.p3_infinite(obj)
        %   Returns the P3 representation of all infinite points. Points at 
        %   infinity have zero as their homogenous coordinate. 
        %   
        %   Inputs:
        %       - obj (sonic.Points3): A Points3 object with k infinite
        %         points
        %
        %   Outputs:
        %       - val (4xk): The P3 realization of all infinite points
        %
        %   Last revised: 2/15/24
        %   Last author: Michael Krause


            val = obj.p3(:, obj.inf_points);
        end

    end
end

