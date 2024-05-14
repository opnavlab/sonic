classdef PointsS2
    % Points lying on the 3D unit sphere (i.e., they only have two degrees
    % of freedom)

    % This can contain a multitude of points. 
    properties
        n       (1, 1) uint64
        u       (3, :) double
        ra_dec  (2, :) double
    end

    methods

        function obj = PointsS2(pts)
        %% obj = PointsS2(pts)
        %   Constructs a Points object containing points lying on a sphere
        %   (i.e., with two degrees of freedom).This object can contain a
        %   multitude of points.
        %   
        %   Inputs:
        %       - pts: A set of points, either as unit vectors or in a 
        %         set of right ascension/declination pairs. The source of
        %         the points is inferred from the dimension of the input
        %         points:
        %           - (2xn double): points are RA/Dec pairs (in radians)
        %           - (3xn double): points are unit vectors
        %
        %   Outputs:
        %       - obj (sonic.PointsS2): PointsS2 object, a container for
        %         a multitude of spherical points.
        %
        %   Last revised: 2/18/24
        %   Last author: Michael Krause

            % Grab the size of the input.
            in_size = size(pts);

            % Number of columns = number of points. Good to know.
            obj.n = in_size(2);
            % Number of rows determines if RA/DEC or unit vector.
            in_dim = in_size(1);

            switch in_dim
                case 2  % Points are RA/DEC pairs:
                   
                    obj.ra_dec = pts;
                    obj.u = sonic.SphereCoords.sphereToCart(pts(1, :)', pts(2, :)');

                case 3  % Points are unit vectors:

                    % Check if these are indeed unit vectors:
                    u_norms = vecnorm(pts);
                    if any(abs(u_norms - 1) > sonic.Tolerances.UnitVecNorm1)
                        error('sonic:PointsS2:notUnitVectors', ...
                            ['Vectors must be unit vectors. ' ...
                            'At least one vector does not have unit norm']);
                    end

                    obj.u = pts;
                    [ra, dec] = sonic.SphereCoords.cartToSphere(pts);
                    obj.ra_dec = [ra'; dec'];

                otherwise
                    error('sonic:PointsS2:IncorrectDimension', ...
                        ['Must input a 2 or 3-by-n matrix of S2 points.' ...
                        'If 2-by-n, must be right ascension and ' ...
                        'declination specified in radians. If 3-by-n, ' ...
                        'must be matrix of unit vectors.']);
            end
        end
    end
end