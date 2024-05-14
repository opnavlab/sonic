classdef Points2 < sonic.GeometryP2

    % This can contain a multitude of points. 
    properties
        n                   (1, 1) uint64 
        p2                  (3, :) double
        inf_points          (1, :) logical
        has_inf_points      (1, 1) logical
    end
    
    properties (Dependent)
        r2                  (2, :) double
        r2_finite           (2, :) double
        p2_finite           (3, :) double
        p2_infinite         (3, :) double
    end

    methods

        function obj = Points2(pts)
        %% obj = Points2(pts)
        %   Constructs a Points object in P2/R2. This object can contain a
        %   multitude of points, including points at infinity. 
        %   
        %   Inputs:
        %       - pts: A set of points, either in P2 or R2. The source of
        %         the points is inferred from the dimension of the input
        %         points:
        %           - (2xn double): points are in R2
        %           - (3xn double): points are in P2
        %
        %   Outputs:
        %       - obj (sonic.Points2): Points2 object, a container for
        %         a multitude of 2D points.
        %
        %   Last revised: 2/15/24
        %   Last author: Michael Krause

            % Grab the size of the input.
            in_size = size(pts);

            % Number of columns = number of points. Good to know.
            obj.n = in_size(2);
            % Number of rows determines if P2 or R2.
            in_dim = in_size(1);

            switch in_dim
                case 2  % Points are in R2
                    
                    % Check for points at infinity:
                    inf_pts = any(pts == Inf, 1);

                    % Save off points and metadata, noting the infinite
                    % points accordingly:
                    obj.p2 = [pts; ones(1, obj.n)];
                    obj.p2(3, inf_pts) = 0;
                    obj.inf_points = inf_pts;
                    obj.has_inf_points = any(inf_pts);
                
                case 3  % Points are in P2

                    % Check for points at infinity and points with negative
                    % homogenous component:
                    inf_pts = abs(pts(3, :)) < sonic.Tolerances.HomNorm;

                    % Normalize the infinite points accordingly:
                    pts(:, ~inf_pts) = pts(:, ~inf_pts)./pts(3, ~inf_pts);
                    pts(3, inf_pts) = 0;        

                    % Save off points and metadata:
                    obj.p2 = pts;
                    obj.inf_points = inf_pts;
                    obj.has_inf_points = any(inf_pts);
                otherwise
                    error('sonic:Points2:IncorrectDimension', ...
                        'Must input a 2 or 3-by-n matrix of 2D points.');
            end
        end

        function val = get.r2(obj)
        %% val = get.r2(obj)
        %   Returns the R2 representation of all points. Points at infinity
        %   have the value infinity for both coordinates. 
        %   
        %   Inputs:
        %       - obj (sonic.Points2): A Points2 object with n points
        %
        %   Outputs:
        %       - val (2xn): The R2 realization of these points
        %
        %   Last revised: 2/15/24
        %   Last author: Michael Krause

            val = obj.p2(1:2, :);
            val(:, obj.inf_points) = Inf;
        end

        function val = get.r2_finite(obj)
        %% val = get.r2_finite(obj)
        %   Returns the R2 representation of all finite points. Points at 
        %   infinity are excluded.
        %   
        %   Inputs:
        %       - obj (sonic.Points2): A Points2 object with m finite
        %         points
        %
        %   Outputs:
        %       - val (2xm): The R2 realization of all finite points
        %
        %   Last revised: 2/15/24
        %   Last author: Michael Krause

            val = obj.p2(1:2, ~obj.inf_points);
        end

        function val = get.p2_finite(obj)
        %% val = get.p2_finite(obj)
        %   Returns the P2 representation of all finite points. Points at 
        %   infinity are excluded.
        %   
        %   Inputs:
        %       - obj (sonic.Points2): A Points2 object with m finite
        %         points
        %
        %   Outputs:
        %       - val (3xm): The P2 realization of all finite points
        %
        %   Last revised: 2/15/24
        %   Last author: Michael Krause

            val = obj.p2(:, ~obj.inf_points);
        end

        function val = get.p2_infinite(obj)
        %% val = get.p2_infinite(obj)
        %   Returns the P2 representation of all infinite points. Points at 
        %   infinity have zero as their homogenous coordinate. 
        %   
        %   Inputs:
        %       - obj (sonic.Points2): A Points2 object with k infinite
        %         points
        %
        %   Outputs:
        %       - val (3xk): The P2 realization of all infinite points
        %
        %   Last revised: 2/15/24
        %   Last author: Michael Krause

            val = obj.p2(:, obj.inf_points);
        end

    end
end