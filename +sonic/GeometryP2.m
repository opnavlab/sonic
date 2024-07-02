% This software is made available under the MIT license. See SONIC/LICENSE
% for details.
classdef (Abstract) GeometryP2 < sonic.MeetJoinable & matlab.mixin.Heterogeneous

    % For now, only dealing with flats (points, lines, planes):
    methods (Static)

        function res = join(varargin)
            %% res = join(varargin)
            %
            %   Performs the join of an arbitrary list of GeometryP2 objects.
            %   Specific combinations of possible joins are:
            %       - Point & Point
            %       - Point & Line
            %   and these will be composed together automatically by this
            %   method, depending on the arguments passed in.
            %
            %   Inputs:
            %       - varargin (n sonic.GeometryP2): A number of subclasses of
            %         GeometryP2.
            %   Outputs:
            %       - res (1x1 sonic.GeometryP2 or double): The result of the
            %         join operation, which could be another GeometryP2 object or
            %         simply a number.
            %
            %   Last revised: 3/29/24
            %   Last author: Michael Krause

            % Specify options for joining geometry s.t. there is a 1:1
            % correspondence between the handles used for the joins and the
            % set of inputs types into each handle.
            join_handles = {@sonic.GeometryP2.joinPoints, ...
                @sonic.GeometryP2.joinPointLine};

            type_choices = {{'Points2', 'Points2'}
                {'Points2', 'Lines2'}};

            res = sonic.MeetJoinable.meetJoinWrapper('join', ...
                join_handles, type_choices, varargin{:});

        end

        function res = meet(varargin)
            %% res = meet(varargin)
            %
            %   Performs the meet of an arbitrary list of GeometryP2 objects.
            %   Specific combinations of possible meets are:
            %       - Line & Line
            %       - Point & Line
            %   and these will be composed together automatically by this
            %   method, depending on the arguments passed in.
            %
            %   Inputs:
            %       - varargin (n sonic.GeometryP2): A number of subclasses of
            %         GeometryP2.
            %   Outputs:
            %       - res (1x1 sonic.GeometryP2 or double): The result of the
            %         meet operation, which could be another GeometryP2 object or
            %         simply a number.
            %
            %   Last revised: 3/29/24
            %   Last author: Michael Krause

            % Specify options for meeting geometry s.t. there is a 1:1
            % correspondence between the handles used for the meets and the
            % set of inputs types into each handle.
            meet_handles = {@sonic.GeometryP2.meetLineLine, ...
                @sonic.GeometryP2.meetPointLine, ...
                @sonic.GeometryP2.meetLineConic, ...
                @sonic.GeometryP2.meetConicConic};

            type_choices = {{'Lines2', 'Lines2'}, ...
                {'Points2', 'Lines2'},...
                {'Lines2','Conic'},...
                {'Conic', 'Conic'}};

            res = sonic.MeetJoinable.meetJoinWrapper('meet', ...
                meet_handles, type_choices, varargin{:});

        end

    end

    % Join Methods:
    methods (Static, Hidden)

        function line = joinPoints(pt1, pt2)
        %% line = joinPoints(pt1, pt2)
        %
        %   Performs the join of two P2 points, the result of which will be
        %   a line (assuming uniqueness of the points). 
        %
        %   Inputs:
        %       - pt1 (1x1 sonic.Points2): A point in P2. If pt2 contains
        %         multiple points, must contain only one point.
        %       - pt2 (1x1 sonic.Points2): A point in P2. If pt1 contains
        %         multiple points, must contain only one point.
        %
        %   Outputs:
        %       - line (1x1 sonic.Lines2): A line or set of lines in P2
        %
        %   Last revised: 5/06/24
        %   Last author: Michael Krause

            [obj_mult, obj_single] = ...
                sonic.MeetJoinable.verifyOneObjWithMultipleN(pt1, pt2);
      
            % Need to iterate over the object with multiple entries:
            lines = zeros(3, obj_mult.n);
            for idx = 1:obj_mult.n
                lines(:, idx) = cross(obj_mult.p2(:, idx), obj_single.p2);
            end
            % Package back into a single Lines2 object:
            line = sonic.Lines2(lines);

        end

        function val = joinPointLine(pt, line)
        %% val = joinPointLine(pt, line)
        %
        %   Performs the join of a P2 point and a P2 line, the result of
        %   which is a number, which is zero if the point lies along the
        %   line. 
        %
        %   Inputs:
        %       - pt (1x1 sonic.Points2): A point in P2. If line contains
        %         multiple lines, must contain only one point.
        %       - line (1x1 sonic.Lines2): A line in P2. If pt contains
        %         multiple points, must contain only one line.
        %
        %   Outputs:
        %       - val (1x1 double): Result of the join operation. Zero if
        %         the point lies along the line.
        %
        %   Last revised: 5/06/24
        %   Last author: Michael Krause

            [obj_mult, obj_single] = ...
                sonic.MeetJoinable.verifyOneObjWithMultipleN(pt, line);
       
            % Need to iterate over the object with multiple entries:
            val = zeros(1, obj_mult.n);
            for idx = 1:obj_mult.n
                val(idx) = obj_mult.p2(:, idx)'*obj_single.p2;
            end

        end

    end

    % Meet Methods:
    methods (Static, Hidden)

        function pt = meetLineLine(line1, line2)
        %% pt = meetLineLine(line1, line2)
        %
        %   Performs the meet of two P2 lines, the result of
        %   which is a point in P2.
        %
        %   Inputs:
        %       - line1 (1x1 sonic.Lines2): A line in P2. If line2 contains
        %         multiple lines, must contain only one line.
        %       - line2 (1x1 sonic.Lines2): A line in P2. If line1 contains
        %         multiple lines, must contain only one line.
        %
        %   Outputs:
        %       - pt (1x1 sonic.Points2): Meet of two lines in P2
        %
        %   Last revised: 5/06/24
        %   Last author: Michael Krause

            [obj_mult, obj_single] = ...
                sonic.MeetJoinable.verifyOneObjWithMultipleN(line1, line2);
        
            % Need to iterate over the object with multiple entries:
            pts = zeros(3, obj_mult.n);
            for idx = 1:obj_mult.n
                pts(:, idx) = cross(obj_mult.p2(:, idx), obj_single.p2);
            end
            % Package back into a single Points2 object:
            pt = sonic.Points2(pts);

        end

        function val = meetPointLine(pt, line)
        %% val = meetPointLine(pt, line)
        %
        %   Performs the meet of a P2 point and a P2 line, the result of
        %   which is a number, which is zero if the point lies along the
        %   line. 
        %
        %   Inputs:
        %       - pt (1x1 sonic.Points2): A point in P2. If line contains
        %         multiple lines, must contain only one point.
        %       - line (1x1 sonic.Lines2): A line in P2. If pt contains
        %         multiple points, must contain only one line.
        %
        %   Outputs:
        %       - val (1x1 double): Result of the meet operation. Zero if
        %         the point lies along the line.
        %
        %   Last revised: 5/06/24
        %   Last author: Michael Krause

            [obj_mult, obj_single] = ...
                sonic.MeetJoinable.verifyOneObjWithMultipleN(pt, line);

            % Need to iterate over the object with multiple entries:
            val = zeros(1, obj_mult.n);
            for idx = 1:obj_mult.n
                val(idx) = obj_mult.p2(:, idx)'*obj_single.p2;
            end

        end

        function pts = meetLineConic(line,conic)
        %% pts = meetLineConic(line,conic)
        %
        %   Performs the meet of a line and a conic, the result of
        %   which are two points in P2
        %
        %   Inputs:
        %       - line (a sonic.Line2 object): A line in P2.
        %       - conic (a sonic.Conic object): A non-degenerate conic.
        %   Outputs:
        %       - pts (2x1 sonic.Points2): two points of intersection of
        %         line and conic
        %
        %   Last revised: 04/22/2024
        %   Last author: Michela Mancini

            locus = conic.locus;

            pts = sonic.Conic.meetLineConic(line,locus);

        end

        function pts = meetConicConic(conic1,conic2)
        %% pts = meetConicConic(conic1,conic2)
        %
        %   Performs the meet of two conics, the result of which are
        %   four points in P2.
        %
        %   Inputs:
        %       - conic1 (a sonic.Conic object): A non-degenerate conic.
        %       - conic2 (a sonic.Conic object): A non-degenerate conic.
        %   Outputs:
        %       - pts (4x1 sonic.Points2): four points of intersection of
        %         the two conics
        %
        %   Last revised: 04/22/2024
        %   Last author: Michela Mancini

        locus1 = conic1.locus;
        locus2 = conic2.locus;

        pts = sonic.Conic.meetConicConic(locus1,locus2);

        end

    end

end

