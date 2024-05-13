classdef Lines3 < sonic.GeometryP3 

    properties
        % Number of lines represented in this object:
        n                       (1, 1)      uint64

        % Primarily encoding these lines in Plucker coords:
        plucker                 (6, :)      double
    end

    properties (Dependent)

        % Dual is just flip of Plucker coords:
        plucker_dual            (6, :)      double

        % Matrix is calculated from the Plucker coords on the fly:
        plucker_mtx             (4, 4, :)   double

        % Matrix dual is calculated from the Plucker dual coords on the
        % fly:
        plucker_mtx_dual        (4, 4, :)   double

        % Can extract moment and direction of a line from the Plucker
        % coordinates as well:
        moment                  (3, :)      double
        direction               (3, :)      double
        %% Perpendicular point
    end
    
    methods
        function obj = Lines3(varargin)
        %% obj = Lines3(varargin)
        %
        %   Container for Lines in P3. Can hold a multitude of lines, and
        %   can be parameterized by Plucker coordinates, Plucker matrices,
        %   or points and directions. See `Inputs` for details. 
        %
        %   Inputs:
        %       - plucker_coords (6xn double): A set of Plucker coordinates
        %       defining `n` lines in P3. Each column corresponds to a
        %       line. 
        %
        %       OR
        %
        %       - plucker_mtx (4x4xn double): A set of Plucker matrices
        %       defining `n` lines in P3. Each 4x4 matrix corresponds to a
        %       line. 
        %
        %       OR
        %
        %       - points (3xn double): A set of `n` points, where the n-th
        %       line passes through the n-th point.
        %       - direc (3xn double): A set of `n` unit vectors, where the
        %       n-th line is parallel to this unit vector. 
        %
        %   Outputs:
        %       - obj (1x1 sonic.Lines3): Lines3 object, containing
        %       information about `n` lines. 
        %
        %   Last revised: 4/08/24
        %   Last author: Michael Krause
        
            switch nargin
                case 1
    
                    % If given a single input, we either have a set of 
                    % Plucker  vectors (6 by n) or we have set of Plucker 
                    % matrices (4 by 4 by n). We thus filter by the number 
                    % of dimensions and length of these dimensions.
    
                    plucker_param = varargin{1};
                    param_size = size(plucker_param);
                    % If we are given Plucker coords directly:
                    if (length(param_size) == 2) && param_size(1) == 6
                   
                            % Enforce the Plucker-Grassman relation:
                            pg_rel = plucker_param(3, :).*plucker_param(4, :) + ...
                                plucker_param(1, :).*plucker_param(6, :) + ...
                                plucker_param(5, :).*plucker_param(2, :);

                            if any(abs(pg_rel) > sonic.Tolerances.PGRelation)
                                error('sonic:Lines3:PGRelViolated', ...
                                    ['At least one line specified ' ...
                                    'violates the Plucker-Grassman ' ...
                                    'constraints, and thus is not ' ...
                                    'properly formed.']);
                            end

                            % Just directly assign in:
                            obj.plucker = plucker_param;
    
                            % Number of lines = number of columns:
                            obj.n = param_size(end);
    
                    % If we're given Plucker matrices:
                    elseif (length(param_size) == 2 || length(param_size) == 3) && ...
                            (param_size(1) == 4 && param_size(2) == 4)

                        % Verify that Plucker matrix is skew-symmetric:
                        if ~all(sonic.Lines3.isSkewSym(plucker_param))
                            error('sonic:Lines3:PluckerMtxNotSkewSym', ...
                                ['When parameterizing a line using ' ...
                                'Plucker matrices, the matrices must ' ...
                                'be skew-symmetric. At least one matrix ' ...
                                'is not skew-symmetric (as determined) ' ...
                                'with MATLAB''s `issymmetric()`.']);
                        end

                        if length(param_size) == 2
                            % Number of lines = 1:
                            obj.n = 1;
                        else    % length = 3
                            % Number of lines = number of items along 3rd dim:
                            obj.n = param_size(end);
                        end
    
                        % Enforce the Plucker-Grassman relation:
                        pg_rel = plucker_param(1, 4, :).*plucker_param(2, 3, :) + ...
                            plucker_param(1, 2, :).*plucker_param(3, 4, :) + ...
                            plucker_param(2, 4, :).*plucker_param(3, 1, :);

                        if any(abs(pg_rel) > sonic.Tolerances.PGRelation)
                            error('sonic:Lines3:PGRelViolated', ...
                                ['At least one line specified ' ...
                                'violates the Plucker-Grassman ' ...
                                'constraints, and thus is not ' ...
                                'properly formed.']);
                        end

                        % Extract the appropriate elements from the matrix:
                        obj.plucker = [plucker_param(1, 2, :);
                                       plucker_param(3, 1, :);
                                       plucker_param(1, 4, :);
                                       plucker_param(2, 3, :);
                                       plucker_param(2, 4, :);
                                       plucker_param(3, 4, :)];


                    else
                        error('sonic:Lines3:ParameterizationSizeMismatch', ...
                            ['When instantiating a Lines3 object with ' ...
                            'a single argument, this argument must be ' ...
                            'either a 6xn matrix of Plucker coordinates ' ...
                            'or a 4x4xn matrix of Plucker matrices. ' ...
                            'See <a href="matlab: help(''sonic.Lines3'')">' ...
                            'the documentation of the Line3 object</a> ' ...
                            'for details on correct usage.']);
                    end
    
                case 2
    
                    % If given two inputs, these correspond to points and 
                    % unit vectors (directions). 

                    points = varargin{1};
                    direc = varargin{2};

                    % Check dimensions:
                    pts_size = size(points);
                    dir_size = size(direc);

                    if (length(pts_size) ~= 2) || (pts_size(1) ~= 3)
                        error('sonic:Lines3:BadPointsSize', ...
                            ['When instantiating a Lines3 object with ' ...
                            'points and unit vectors, the points and ' ...
                            'unit vector matrices must both be 3xn, ' ...
                            'with points and unit vectors lying along ' ...
                            'the columns. '...
                            'See <a href="matlab: help(''sonic.Lines3'')">' ...
                            'the documentation of the Line3 object</a> ' ...
                            'for details on correct usage.']);
                    end

                    if (length(dir_size) ~= 2) || (dir_size(1) ~= 3)
                        error('sonic:Lines3:BadDirectionsSize', ...
                            ['When instantiating a Lines3 object with ' ...
                            'points and unit vectors, the points and ' ...
                            'unit vector matrices must both be 3xn, ' ...
                            'with points and unit vectors lying along ' ...
                            'the columns. '...
                            'See <a href="matlab: help(''sonic.Lines3'')">' ...
                            'the documentation of the Line3 object</a> ' ...
                            'for details on correct usage.']);
                    end

                    if pts_size(2) ~= dir_size(2)
                        error('sonic:Lines3:PointDirectionMismatch', ...
                            ['When instantiating a Lines3 object with ' ...
                            'points and unit vectors, the points and ' ...
                            'unit vector matrices must both be 3xn, ' ...
                            'with points and unit vectors lying along ' ...
                            'the columns. '...
                            'See <a href="matlab: help(''sonic.Lines3'')">' ...
                            'the documentation of the Line3 object</a> ' ...
                            'for details on correct usage.']);
                    end

                    % The resultant moments will lie along the columns:
                    moment = cross(points, direc);

                    % The moments and directions can be used to form a
                    % Plucker matrix. We can short-circuit this a bit and
                    % directly extract the indices of interest:
                    obj.plucker = [-moment(3, :);
                                   -moment(2, :);
                                    direc(1, :);
                                   -moment(1, :);
                                    direc(2, :);
                                    direc(3, :)];

                    obj.n = pts_size(2);

                otherwise
                    error('sonic:Lines3:InvalidParameterization', ...
                            ['Parameterization of Line3 not recognized. ' ...
                            'See <a href="matlab: help(''sonic.Lines3'')">' ...
                            'the documentation of the Line3 object</a> ' ...
                            'for details on correct usage.']);

                    
            end

        end

        function val = get.plucker_dual(obj)
        %% val = get.plucker_dual(obj)
        %
        %   Returns the dual of the Plucker coordinates for the lines
        %   contained within the Lines3 object.
        %
        %   Inputs:
        %       - obj (1x1 sonic.Lines3): Lines3 object, containing `n`
        %       lines.
        %   
        %   Outputs:
        %       - val (6xn double): Dual of the Plucker coordinates
        %       corresponding to each line. Coordinates make up the columns
        %       of this matrix. 
        %
        %   Last revised: 04/08/24
        %   Last author: Michael Krause
        
            c = obj.plucker;
            val = flip(c);
        
        end

        function val = get.plucker_mtx(obj)
        %% val = get.plucker_mtx(obj)
        %
        %   Returns the matrix of the Plucker coordinates for the lines
        %   contained within the Lines3 object.
        %
        %   Inputs:
        %       - obj (1x1 sonic.Lines3): Lines3 object, containing `n`
        %       lines.
        %   
        %   Outputs:
        %       - val (4x4xn double): Plucker coordinates matrices
        %       corresponding to each line. Concatenated along the third
        %       dimension.
        %
        %   Last revised: 04/08/24
        %   Last author: Michael Krause
        
            % Grab the Plucker coords and reshape:
            c = reshape(obj.plucker, [6, 1, obj.n]);
            z = zeros(1, 1, obj.n);
            val = [    z,      c(1,:,:), -c(2,:,:),  c(3,:,:);
                   -c(1,:,:),         z,  c(4,:,:),  c(5,:,:);
                    c(2,:,:), -c(4,:,:),         z,  c(6,:,:);
                   -c(3,:,:), -c(5,:,:), -c(6,:,:),         z];
        
        end

        function val = get.plucker_mtx_dual(obj)
        %% val = get.plucker_mtx_dual(obj)
        %
        %   Returns the matrix of the Plucker coordinates for the duals of
        %   the lines contained within the Lines3 object.
        %
        %   Inputs:
        %       - obj (1x1 sonic.Lines3): Lines3 object, containing `n`
        %       lines.
        %   
        %   Outputs:
        %       - val (4x4xn double): Plucker coordinates matrices
        %       corresponding to the dual of each line. Concatenated along 
        %       the third dimension.
        %
        %   Last revised: 04/08/24
        %   Last author: Michael Krause

            % Grab the Plucker dual and reshape:
            c = reshape(obj.plucker_dual, [6, 1, obj.n]);
            z = zeros(1, 1, obj.n);
            val = [    z,      c(1,:,:), -c(2,:,:),  c(3,:,:);
                   -c(1,:,:),         z,  c(4,:,:),  c(5,:,:);
                    c(2,:,:), -c(4,:,:),         z,  c(6,:,:);
                   -c(3,:,:), -c(5,:,:), -c(6,:,:),         z];
        
        end

        function val = get.moment(obj)
        %% val = get.moment(obj)
        %
        %   Returns the moments of each line contained in the Lines3
        %   object. 
        %
        %   Inputs:
        %       - obj (1x1 sonic.Lines3): Lines3 object, containing `n`
        %       lines.
        %   
        %   Outputs:
        %       - val (3xn double): Moments of each line, with the moments
        %       making up the columns. 
        %
        %   Last revised: 04/08/24
        %   Last author: Michael Krause

            val = [-obj.plucker(4, :);
                   -obj.plucker(2, :);
                   -obj.plucker(1, :)];
            
        end

        function val = get.direction(obj)
        %% val = get.direction(obj)
        %
        %   Returns the direction of each line contained in the Lines3
        %   object. 
        %
        %   Inputs:
        %       - obj (1x1 sonic.Lines3): Lines3 object, containing `n`
        %       lines.
        %   
        %   Outputs:
        %       - val (3xn double): Directions of each line, with the 
        %       directions making up the columns. 
        %
        %   Last revised: 04/08/24
        %   Last author: Michael Krause

            val = [obj.plucker(3, :);
                   obj.plucker(5, :);
                   obj.plucker(6, :)];
            
        end

    end

    methods (Static, Hidden)

        function res = isSkewSym(A)
        %% res = isSkewSym(A)
        %
        %   Checks for skew-symmetry across m-by-m-by-n matrices, returning
        %   an n-by-1 vector of whether each is skew-symmetric. Invokes
        %   MATLAB's builtin `issymmetric()` method internally.
        %
        %   Inputs:
        %       - A (mxmxn double): 3D matrix, square along the first two
        %       dimensions, to check for skew-symmetry.
        %   
        %   Outputs:
        %       - res (1xn logical): Vector of skew-symmetry check results,
        %       one per square matrix.
        %
        %   Last revised: 5/06/24
        %   Last author: Michael Krause

            sz = size(A);
            switch length(sz)
                case 2
                    n = 1;
                case 3
                    n = sz(3);
                otherwise
                    error('sonic:Lines3:isSkewSym:invalidDimensions', ...
                        ['Invalid dimensions for argument. Must be ' ...
                        'm-by-m-by-n (n >= 1).']);
            end

            res = false(1, n);
            for idx = 1:n
                res(idx) = issymmetric(A(:, :, idx), "skew");
            end

        end

    end

end

