% This software is made available under the MIT license. See SONIC/LICENSE
% for details.
classdef (Abstract) GeometryP3 < sonic.MeetJoinable & matlab.mixin.Heterogeneous
    
    % For now, only dealing with flats (points, lines, planes):
    methods (Static)

        function res = join(varargin)
        %% res = join(varargin)
        %
        %   Performs the join of an arbitrary list of GeometryP3 objects.
        %   Specific combinations of possible joins are:
        %       - Point & Point
        %       - Point & Line
        %       - Point & Plane
        %       - Line & Line
        %   and these will be composed together automatically by this
        %   method, depending on the arguments passed in. 
        %
        %   Inputs:
        %       - varargin (n sonic.GeometryP3): A number of subclasses of
        %         GeometryP3. 
        %   Outputs:
        %       - res (1x1 sonic.GeometryP3 or double): The result of the
        %         join operation, which could be another GeometryP3 object or
        %         simply a number. 
        %
        %   Last revised: 5/06/24
        %   Last author: Michael Krause

            % Specify options for joining geometry s.t. there is a 1:1
            % correspondence between the handles used for the joins and the
            % set of inputs types into each handle. 
            join_handles = {@sonic.GeometryP3.joinPoints,     ...
                            @sonic.GeometryP3.joinPointLine,  ...
                            @sonic.GeometryP3.joinPointPlane, ...
                            @sonic.GeometryP3.joinLineLine };

            type_choices = {{'Points3', 'Points3'}, ...
                            {'Points3', 'Lines3' }, ...
                            {'Points3', 'Planes3' }, ...
                            {'Lines3' , 'Lines3' }};

            res = sonic.MeetJoinable.meetJoinWrapper('join', ...
                join_handles, type_choices, varargin{:});

        end

        function res = meet(varargin)    
        %% res = meet(varargin)  
        %
        %   Performs the meet of an arbitrary list of GeometryP3 objects.
        %   Specific combinations of possible meets are:
        %       - Plane & Plane
        %       - Plane & Line
        %       - Plane & Point
        %       - Line & Line
        %   and these will be composed together automatically by this
        %   method, depending on the arguments passed in. 
        %
        %   Inputs:
        %       - varargin (n sonic.GeometryP3): A number of subclasses of
        %         GeometryP3. 
        %   Outputs:
        %       - res (1x1 sonic.GeometryP3 or double): The result of the
        %         meet operation, which could be another GeometryP3 object or
        %         simply a number. 
        %
        %   Last revised: 5/06/24
        %   Last author: Michael Krause 

            % Specify options for meeting geometry s.t. there is a 1:1
            % correspondence between the handles used for the meets and the
            % set of inputs types into each handle. 
            meet_handles = {@sonic.GeometryP3.meetPlanePlane, ...
                            @sonic.GeometryP3.meetPlaneLine , ...
                            @sonic.GeometryP3.meetPlanePoint, ...
                            @sonic.GeometryP3.meetLineLine };

            type_choices = {{'Planes3', 'Planes3' }, ...
                            {'Planes3', 'Lines3' }, ...
                            {'Planes3', 'Points3'}, ...
                            {'Lines3',  'Lines3' }};

            res = sonic.MeetJoinable.meetJoinWrapper('meet', ...
                meet_handles, type_choices, varargin{:});

        end

    end

    % Join Methods:
    methods (Static, Hidden)

        function line = joinPoints(pt1, pt2)
        %% line = joinPoints(pt1, pt2)
        %
        %   Performs the join of two P3 points, the result of which will be
        %   a line (assuming uniqueness of the points). 
        %
        %   Inputs:
        %       - pt1 (1x1 sonic.Points3): A point in P3. If pt2 contains
        %         multiple points, must contain only one point.
        %       - pt2 (1x1 sonic.Points3): A point in P3. If pt1 contains
        %         multiple points, must contain only one point.
        %
        %   Outputs:
        %       - line (1x1 sonic.Lines3): A line or set of lines in P3
        %
        %   Last revised: 5/06/24
        %   Last author: Michael Krause

            arguments
                pt1     (1, 1)      sonic.Points3
                pt2     (1, 1)      sonic.Points3
            end

            [obj_mult, obj_single] = ...
                sonic.MeetJoinable.verifyOneObjWithMultipleN(pt1, pt2);
            
            % Need to iterate over the object with multiple entries.
            % Store each computation as a Plucker mtx:
            lines = zeros(4, 4, obj_mult.n);
            B = obj_single.p3;
            for idx = 1:obj_mult.n
                A = obj_mult.p3(:, idx);
                lines(:, :, idx) = A*B' - B*A';
            end
            % Package back into a single Lines3 object:
            line = sonic.Lines3(lines);
            
        end

        function plane = joinPointLine(pt, line)
        %% plane = joinPointLine(pt, line)
        %
        %   Performs the join of a P3 point and a P3 line, the result of
        %   which (assuming the point does not lie on the line) is a plane
        %   in P3. 
        %
        %   Inputs:
        %       - pt (1x1 sonic.Points3): A point in P3. Must contain only
        %         one point.
        %       - line (1x1 sonic.Lines3): A line in P3. Must contain only
        %         one line.
        %
        %   Outputs:
        %       - plane (1x1 sonic.Plane3): A plane in P3
        %
        %   Last revised: 5/06/24
        %   Last author: Michael Krause

            [obj_mult, obj_single] = ...
                sonic.MeetJoinable.verifyOneObjWithMultipleN(line, pt);

            planes = zeros(4, obj_mult.n);
            switch class(obj_mult)
                case 'sonic.Lines3'
                    B = obj_single.p3;
                    for idx = 1:obj_mult.n
                        A = obj_mult.plucker_mtx_dual(:, :, idx);
                        planes(:, idx) = A*B;
                    end
                case 'sonic.Points3'
                    A = obj_single.plucker_mtx_dual;
                    for idx = 1:obj_mult.n
                        B = obj_mult.p3(:, idx);
                        planes(:, idx) = A*B;
                    end
                otherwise
                    error('sonic:GeometryP3:joinPointLine:incorrectTypes', ...
                        ['Must pass in a Points3 and Lines3 object ' ...
                        'to joinPointLine().']);
            end
            
            % Package back into a single Planes3 object:
            plane = sonic.Planes3(planes);
        end

        function val = joinPointPlane(pt, plane)
        %% val = joinPointPlane(pt, plane)
        %
        %   Performs the join of a P3 point and a P3 plane, the result of
        %   which is a number.
        %
        %   Inputs:
        %       - pt (1x1 sonic.Points3): A point in P3. Must contain only
        %         one point.
        %       - plane (1x1 sonic.Plane3): A plane in P3
        %
        %   Outputs:
        %       - val (1x1 double): Join of point and plane in P3
        %
        %   Last revised: 3/28/24
        %   Last author: Michael Krause


            [obj_mult, obj_single] = ...
                sonic.MeetJoinable.verifyOneObjWithMultipleN(pt, plane);
            
            % Need to iterate over the object with multiple entries.
            % Store each computation as a Plucker mtx:
            val = zeros(1, obj_mult.n);
            B = obj_single.p3;
            for idx = 1:obj_mult.n
                A = obj_mult.p3(:, idx);
                val(idx) = A'*B;
            end

        end

        function val = joinLineLine(line1, line2)
        %% val = joinLineLine(line1, line2)
        %
        %   Performs the join of two P3 lines, the result of
        %   which is a number. 
        %
        %   Inputs:
        %       - line (1x1 sonic.Lines3): A line in P3. If line2 contains
        %         multiple lines, must contain only one line.
        %       - line (1x1 sonic.Lines3): A line in P3. If line1 contains
        %         multiple lines, must contain only one line.
        %
        %   Outputs:
        %       - val (1x1 double): Join of line and line in P3
        %
        %   Last revised: 4/18/24
        %   Last author: Michael Krause

            [obj_mult, obj_single] = ...
                sonic.MeetJoinable.verifyOneObjWithMultipleN(line1, line2);
    
            % Need to iterate over the object with multiple entries.
            val = zeros(1, obj_mult.n);
            B = obj_single.plucker;
            for idx = 1:obj_mult.n
                A = obj_mult.plucker(:, idx);

                % Boils down to u_A'*m_B + u_B'*m_A. That is, a mutual
                % evaluation of the Plucker-Grassman constraint. If they're
                % the same time, this will be zero. 
                val(idx) = A(3)*B(4) + A(4)*B(3) + A(5)*B(2) + A(2)*B(5) + A(6)*B(1) + A(1)*B(6);
            end
        end

    end

    % Meet Methods:
    methods (Static, Hidden)

        function line = meetPlanePlane(plane1, plane2)
        %% line = meetPlanePlane(plane1, plane2)
        %
        %   Performs the meet of two P3 planes, the result of
        %   which is a P3 line (provided the two planes are not coplanar).
        %
        %   Inputs:
        %       - plane1 (1x1 sonic.Plane3): A plane in P3
        %       - plane2 (1x1 sonic.Plane3): A plane in P3
        %
        %   Outputs:
        %       - line (1x1 sonic.Lines3): Meet of plane and plane in P3
        %
        %   Last revised: 5/06/24
        %   Last author: Michael Krause

            [obj_mult, obj_single] = ...
                sonic.MeetJoinable.verifyOneObjWithMultipleN(plane1, plane2);
            
            % Need to iterate over the object with multiple entries.
            % Store each computation as a Plucker mtx:
            dual_lines = zeros(4, 4, obj_mult.n);
            B = obj_single.p3;
            for idx = 1:obj_mult.n
                A = obj_mult.p3(:, idx);
                dual_lines(:, :, idx) = A*B' - B*A';
            end
            % Package back into a single Lines3 object, while inverting the
            % dual that we found:
            line = sonic.Lines3(sonic.Lines3(dual_lines).plucker_dual);
        
        end

        function point = meetPlaneLine(plane, line)
        %% point = meetPlaneLine(plane, line)
        %
        %   Performs the meet of a P3 plane and a P3 line, the result of
        %   which is a P3 point (provided the line and plane are not
        %   coplanar). 
        %
        %   Inputs:
        %       - plane (1x1 sonic.Plane3): A plane in P3
        %       - line (1x1 sonic.Lines3): A line in P3. Must contain only
        %         one line.
        %
        %   Outputs:
        %       - point (1x1 sonic.Points3): Meet of plane and line in P3
        %
        %   Last revised: 5/06/24
        %   Last author: Michael Krause

            [obj_mult, obj_single] = ...
                sonic.MeetJoinable.verifyOneObjWithMultipleN(line, plane);

            points = zeros(4, obj_mult.n);
            switch class(obj_mult)
                case 'sonic.Lines3'
                    B = obj_single.p3;
                    for idx = 1:obj_mult.n
                        A = obj_mult.plucker_mtx(:, :, idx);
                        points(:, idx) = A*B;
                    end
                case 'sonic.Planes3'
                    A = obj_single.plucker_mtx;
                    for idx = 1:obj_mult.n
                        B = obj_mult.p3(:, idx);
                        points(:, idx) = A*B;
                    end
                otherwise
                    error('sonic:GeometryP3:meetPlaneLine:incorrectTypes', ...
                        ['Must pass in a Planes3 and Lines3 object ' ...
                        'to meetPlaneLine().']);
            end
            
            % Package back into a single Points3 object:
            point = sonic.Planes3(points);
        

        end

        function val = meetPlanePoint(plane, pt)
        %% val = meetPlanePoint(plane, pt)
        %
        %   Performs the meet of a P3 plane and a P3 point, the result of
        %   which is a number (signifying the signed distance from the
        %   point to the plane). 
        %
        %   Inputs:
        %       - plane (1x1 sonic.Plane3): A plane in P3
        %       - pt (1x1 sonic.Points3): A point in P3. Must contain only
        %         one point.
        %
        %   Outputs:
        %       - line (1x1 sonic.Lines3): Meet of plane and point in P3
        %
        %   Last revised: 5/06/24
        %   Last author: Michael Krause

            [obj_mult, obj_single] = ...
                sonic.MeetJoinable.verifyOneObjWithMultipleN(plane, pt);
            
            % Need to iterate over the object with multiple entries.
            % Store each computation as a Plucker mtx:
            val = zeros(1, obj_mult.n);
            B = obj_single.p3;
            for idx = 1:obj_mult.n
                A = obj_mult.p3(:, idx);
                val(idx) = A'*B;
            end

        end

        function val = meetLineLine(line1, line2)
        %% val = meetLineLine(line1, line2)
        %
        %   Performs the meet of two P3 lines, the result of
        %   which is a number.
        %
        %   Inputs:
        %       - line1 (1x1 sonic.Lines3): A line in P3. If line2 contains
        %         multiple lines, must contain only one line.
        %       - line2 (1x1 sonic.Lines3): A line in P3. If line1 contains
        %         multiple lines, must contain only one line.
        %
        %   Outputs:
        %       - val (1x1 double): Meet of two lines in P3
        %
        %   Last revised: 5/06/24
        %   Last author: Michael Krause
        
            [obj_mult, obj_single] = ...
                sonic.MeetJoinable.verifyOneObjWithMultipleN(line1, line2);
    
            % Need to iterate over the object with multiple entries.
            val = zeros(1, obj_mult.n);
            B = obj_single.plucker;
            for idx = 1:obj_mult.n
                A = obj_mult.plucker(:, idx);

                % Boils down to u_A'*m_B + u_B'*m_A. That is, a mutual
                % evaluation of the Plucker-Grassman constraint. If they're
                % the same time, this will be zero. 
                val(idx) = A(3)*B(4) + A(4)*B(3) + A(5)*B(2) + A(2)*B(5) + A(6)*B(1) + A(1)*B(6);
            end
        end

    end
end

