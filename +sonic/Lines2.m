% This software is made available under the MIT license. See SONIC/LICENSE
% for details.
classdef Lines2 < sonic.GeometryP2

    properties
        % Number of lines represented in this object:
        n                   (1, 1)      uint64

        % Primarily storing parameterization in P2:
        p2                  (3, :)      double

        % Line angles:
        angle               (1, :)      double

        % Perpendicular distances:
        perp_dist           (1, :)      double

        % Record any lines at infinity:
        inf_lines           (1, :)      logical
        has_inf_lines       (1, 1)      logical
    end

    properties (Dependent)
        p2_finite           (3, :)      double
        p2_infinite         (3, :)      double
        angle_finite        (1, :)      double
        perp_dist_finite    (1, :)      double
    end
    
    methods
        function obj = Lines2(varargin)
        %% obj = Lines2(varargin)
        %
        %   Container for Lines in P2. Can specify the lines either
        %   directly in P2 coordinates, or using angle and perpendicular
        %   distance (see Inputs list for details). If specified with angle
        %   and perpendicular distance, the convention for conversion into
        %   P2 coordinates is:
        %       
        %           Coords = [cos(angle); sin(angle); -distance];
        %
        %   Note that this object can contain an arbitrary number of lines.
        %
        %   Inputs:
        %       - coeffs (3xn double): P2 coordinates of `n` lines, with 
        %         the coordinates of each line as the columns of the matrix.
        %   
        %           OR
        %   
        %       - angle (n double): Vector of angles in radians (measured 
        %         clockwise from the +x axis) for each line, with the n-th 
        %         entry cooresponding to the n-th line. Must have the same 
        %         length as `perp_dist`.
        %       - perp_dist (n double): Vector of perpendicular distances 
        %         from the origin to each line, with the n-th entry
        %         cooresponding to the n-th line. Must have same length as
        %         `angle`.
        %
        %   Outputs:
        %       - obj (1x1 sonic.Lines2): Lines2 object, which can contain
        %         a multitude of lines.
        %
        %   Last revised: 4/05/24
        %   Last author: Michael Krause

            switch nargin
                case 1

                    % Sanity check dimensions of P2 repr.
                    coeffs = varargin{1};
                    in_size = size(coeffs);

                    % P2 coefficients should lie on the columns: 
                    if in_size(1) ~= 3
                        error('sonic:Lines2:P2ReprBadSize', ...
                            ['When creating Lines2 object with coordinates ' ...
                            'in P2, must specify input as a 3-by-n ' ...
                            'matrix, with the P2 coordinates of ' ...
                            'each line forming the columns of ' ...
                            'this matrix.']);
                    end

                    small = abs(coeffs) < sonic.Tolerances.SmallNumber;
                     badScale = any(all(small));
                     if badScale
                        error('sonic:Lines2:invalidInput', ...
                            ['Line [0; 0; 0] not member of P2. If ' ...
                            'this error arises from scaling issues, ' ...
                            'consider normalizing the line coefficients.']);
                     end

                    % Number of lines = number of columns:
                    obj.n = in_size(2);
    
                    throughOrigin = abs(coeffs(3, :)) < sonic.Tolerances.HomNorm;
                   
                    % Check for lines at infinity. When we normalize by the
                    % third component, if the first two components are
                    % close enough to zero, we consider this to be a line
                    % at infinity. 
                    norm_coeffs = coeffs./coeffs(3, :);
                    ang_coeffs_norm = vecnorm(norm_coeffs(1:2, :));
                    obj.inf_lines = ang_coeffs_norm < sonic.Tolerances.LineNorm;
                    obj.has_inf_lines = any(obj.inf_lines);

                    % For lines at infinity, set the P2 repr to be [0;0;1].
                    coeffs(1:2, obj.inf_lines) = 0;
                    coeffs(3, obj.inf_lines) = 1;
                    
                    if obj.n>size(coeffs(1:2, obj.inf_lines),2) + size(coeffs(1:2, throughOrigin),2)
                        % For lines not at infinity, normalize such that the
                        % norm of the first two components is 1.
                        coeffs(:, ~obj.inf_lines & ~throughOrigin) = ...
                            norm_coeffs(:, ~obj.inf_lines  & ~throughOrigin )./ang_coeffs_norm(~obj.inf_lines  & ~throughOrigin);
                    end

                    % account separately for line that goes through the origin
                    coeffs(:,throughOrigin) = coeffs(:,throughOrigin)./vecnorm(coeffs(:,throughOrigin));

                    % Store the P2 representation and calculate the angles
                    % and perpendicular distances. For lines at infinity,
                    % angle = NaN, perp distance = Inf.
                    obj.p2 = coeffs;

                    obj.angle = sonic.Math.atan2c(coeffs(2, :), coeffs(1, :));
                    obj.angle(obj.inf_lines) = NaN;

                    obj.perp_dist = -coeffs(3, :);
                    obj.perp_dist(obj.inf_lines) = Inf;
                case 2
                    
                    % In this case, both should be 1-by-n vectors. 
                    ang = varargin{1};
                    perp_dist = varargin{2};

                    % Make sure these are explicitly 1-D vectors of equal
                    % length. Anything else, and we have reason to be 
                    % suspicious the inputs are incorrect.
                    size_dist = size(perp_dist);
                    size_theta = size(ang);

                    if (length(size_dist) ~= 2 || ~any(size_dist == 1)) || ...
                            (length(size_theta) ~= 2 || ~any(size_theta == 1)) || ...
                            (length(perp_dist) ~= length(ang))
                    
                        error('sonic:Lines2:invalidParamSizes', ...
                            ['When creating Lines2 object with ' ...
                            'angles and perpendicular distances, ' ...
                            'must specify both parameters ' ...
                            'as one-dimensional vectors of equal ' ...
                            'length.']);
                    end

                    % Number of lines = number of columns:
                    obj.n = size_dist(2);

                    % Ensure that these are row vectors:
                    perp_dist = perp_dist(:)';
                    ang = ang(:)';

                    % Now, check for lines at infinity. This will happen
                    % when the perpendicular distance is Inf or if the
                    % angle is NaN. 
                    obj.inf_lines = isnan(ang) | isinf(perp_dist);
                    obj.has_inf_lines = any(obj.inf_lines);

                    % Store the angle and perpendicular distances, ensuring
                    % that infinite cases are noted accordingly:
                    obj.angle = ang;
                    obj.angle(obj.inf_lines) = NaN;
                    obj.perp_dist = perp_dist;
                    obj.perp_dist(obj.inf_lines) = Inf;

                    % Assemble into P2, also noting the infinite lines
                    % accordingly:
                    obj.p2 = [cos(ang); sin(ang); -perp_dist];
                    obj.p2(1:2, obj.inf_lines) = 0;
                    obj.p2(3, obj.inf_lines) = 1;
                    
                otherwise
                    error('sonic:Lines2:invalidParameterization', ...
                        ['Invalid number of parameters specified. ' ...
                        'Provide a single 3-by-n matrix of P2 ' ...
                        'coordinates or two equal length vectors ' ...
                        'of angle and perpendicular distance.']);
            end
        end

        %% Dependent Properties:
            
        function val = get.p2_finite(obj)
        %% val = get.p2_finite(obj)
        %   Returns the P2 representation of all finite lines. Lines at 
        %   infinity are excluded.
        %   
        %   Inputs:
        %       - obj (sonic.Lines2): A Lines2 object with m finite
        %         lines
        %
        %   Outputs:
        %       - val (3xm): The P2 realization of all finite lines
        %
        %   Last revised: 4/18/24
        %   Last author: Michael Krause

            val = obj.p2(:, ~obj.inf_lines);

        end

        function val = get.p2_infinite(obj)
        %% val = get.p2_infinite(obj)
        %   Returns the P2 representation of all infinite lines.
        %   
        %   Inputs:
        %       - obj (sonic.Lines2): A Lines2 object with k lines at
        %         infinity
        %
        %   Outputs:
        %       - val (3xk): The P2 realization of all infinite lines
        %
        %   Last revised: 4/18/24
        %   Last author: Michael Krause

            val = obj.p2(:, obj.inf_lines);

        end

        function val = get.angle_finite(obj)
        %% val = get.angle_finite(obj)
        %   Returns the angle (in radians) of every finite line. Angle is
        %   measured +CCW from the +x axis. 
        %   
        %   Inputs:
        %       - obj (sonic.Lines2): A Lines2 object with m finite
        %         lines
        %
        %   Outputs:
        %       - val (1xm): The angle of each line, in radians.
        %
        %   Last revised: 4/18/24
        %   Last author: Michael Krause

            val = obj.angle(~obj.inf_lines);

        end

        function val = get.perp_dist_finite(obj)
        %% val = get.perp_dist_finite(obj)
        %   Returns the perpendicular distance of every finite line. 
        %   
        %   Inputs:
        %       - obj (sonic.Lines2): A Lines2 object with m finite
        %         lines
        %
        %   Outputs:
        %       - val (1xm): The perpendicular distance of each line
        %
        %   Last revised: 4/18/24
        %   Last author: Michael Krause

            val = obj.perp_dist(~obj.inf_lines);

        end
        
    end

end

