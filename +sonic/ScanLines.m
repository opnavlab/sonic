classdef ScanLines
    properties
        lines (:, 1) sonic.Lines2 = sonic.Lines2.empty()
        ispos (:, 1) logical 
    end

    methods
        function obj = ScanLines(lines, ispos)
        %% obj = ScanLines(lines, ispos)
        %   Constructs a sonic.ScanLines object which contains all the
        %   specified lines and their properties as well as a list of
        %   direction flags for each line which contains 
        %
        %   Inputs
        %       - lines (1,1 sonic.Lines2): Lines2 object containing all
        %         the lines for a scan 
        %       - ispos (n,1 logical): list of direction flagss for each
        %         line indicating the desired direction of scan. true for a
        %         positive direction (in direction of increasing u or top to
        %         bottom for verticle lines) or false for a negative 
        %         direction (in direction of decreasing u or bottom to top    
        %         for verticle lines).
        %
        %   Outputs
        %       - obj (1,1 sonic.ScanLines): ScanLines object containing
        %         the lines for a scan and their directionality
        %
        %   Last revised: 04/26/24
        %   Last author: Ava Thrasher
            arguments
                lines (1,1) sonic.Lines2
                ispos (:,1) logical
            end

            obj.lines = lines;
            obj.ispos = ispos;

        end
    end

    methods (Static)

        function obj = getParallelLines(img_obj, drho, u_illum_obj)
        %% obj = getParallelLines(img_obj, drho, u_illum_obj)
        %   Constructs a sonic.ScanLines object which contains all the
        %   lines which are parallel, intersect with the image, and are
        %   spaced by drho
        %
        %   Inputs
        %       - img_obj (1,1 sonic.Image): Image for which to generate 
        %         parallel scan lines 
        %       - drho (1,1 double): Perpendicular distance between 
        %         parallel lines
        %       - u_illum_obj (1,1 sonic.Points2): Illumination direction
        %         in the image specified from the image top left corner
        %
        %   Outputs
        %       - obj (1,1 sonic.ScanLines): ScanLines object containing
        %         the lines for a scan and their directionality, which 
        %         are parallel and intersect with the specified image
        %
        %   Last revised: 05/01/24
        %   Last author: Ava Thrasher
            arguments
                img_obj (1,1) sonic.Image
                drho (1,1) double
                u_illum_obj (1,1) sonic.Points2
            end
            
            % Get size of image 
            rows = img_obj.rows;
            cols = img_obj.cols;

            % Get theta based on illumination direction
            u_illum = u_illum_obj.r2;
            theta = atan2(u_illum(1),-u_illum(2));

            % Find dirFalg based on direction of illumination          
            if theta <= pi && theta > 0
                dirFlag = 1;
            else
                dirFlag = 0;
            end

            % Find max and min rho for given theta
            toler = sonic.Tolerances.SmallAngle;
            if abs(theta) < toler
                rho_min = 0;
                rho_max = cols;
            elseif abs(theta - pi/2) < toler
                rho_min = 0;
                rho_max = rows;
            elseif abs(theta + pi/2) < toler
                rho_min = -rows;
                rho_max = 0;
            elseif abs(abs(theta) - pi) < toler
                rho_min = -cols;
                rho_max = 0;
            elseif theta < pi/2 && theta > 0
                rho_min = 0;
                rho_max = (rows*tan(theta)+cols)*cos(theta);
            elseif theta < 0 && theta > -pi/2
                rho_min = rows*sin(theta);
                rho_max = cols*cos(theta);
            elseif theta < pi && theta > pi/2
                rho_min = -cols*sin(theta);
                rho_max = -rows*cos(theta);
            elseif theta < -pi/2 && theta >-pi
                rho_min = (rows*tan(theta)+cols)*cos(theta);
                rho_max = 0;
            end

            % Generate Lines2 object with scanlines in image
            perp_dists = rho_min:drho:rho_max;
            numLines = size(perp_dists,2);
            lines = sonic.Lines2(theta.*ones(numLines,1), perp_dists);

            % Create ScanLines object
            dirFlags = logical(dirFlag.*ones(numLines,1));
            obj = sonic.ScanLines(lines,dirFlags);

        end
    end

end