classdef Project 

    methods (Static)
    
        function [proj_obj, did_proj] = pinhole(to_proj, proj_behind_cam)
        %% [proj_obj, did_proj] = pinhole(to_proj, proj_behind_cam)
        %   Performs a rectilinear projection of 3D points/lines/quadrics 
        %   into a 2D image plane. Assumes a camera origin at z = 0, with
        %   an image plane lying at z = 1.
        %   
        %   NOTE: Currently only supports points and Ellipsoids.
        %
        %   Inputs:
        %       - to_proj (1x1 sonic.Points3 OR sonic.PointsS2 OR 
        %         sonic.Ellipsoid): Object to project. Must be expressed in 
        %         the camera frame.
        %       - proj_behind_cam (1x1 logical): Flag indicating whether
        %         to project points that lie behind the plane of the camera.
        %         OPTIONAL: defaults to FALSE. 
        %
        %   Outputs:
        %       - proj_obj (1x1): Set of projected object(s),
        %         lying on the image plane at z = 1.
        %       - did_proj (1xn logical): Mapping of which points were
        %         successfully projected. If proj_behind_cam was set to false
        %         (as it is by default), only points in front of the camera
        %         will be projected. Thus, it will be a 1xn vector (where n
        %         is the number of points contained in to_proj), with m true
        %         entries, where m points lie in front of the camera. 
        %
        %   Last revised: 03/06/24
        %   Last author: Michael Krause

            arguments
                to_proj             (1, 1)      
                proj_behind_cam     (1, 1)      logical     = false
            end

            if isa(to_proj, 'sonic.Points3') || isa(to_proj, 'sonic.PointsS2')

                % Depending on input type, need to grab the R3 repr 
                % differently. Rest of algorithm is the same though.
                if isa(to_proj, 'sonic.Points3')
                    % Pull out all finite points.
                    pts_raw = to_proj.r3_finite;
                else % Else is PointsS2
                    % Pull out unit vec repr of points:
                    pts_raw = to_proj.u;
                end

                if proj_behind_cam
                    % Allow all points to project:
                    pts_to_proj = pts_raw;
                    did_proj = true(1, size(pts_raw, 2));
                else
                    % Filter to those with z > 0:
                    did_proj = pts_raw(3, :) > 0;
                    pts_to_proj = pts_raw(:, did_proj);
                end

                % Then just divide by the z coordinate.
                pts_plane = pts_to_proj./pts_to_proj(3, :);

                % Then insert back into a Points2 object.
                proj_obj = sonic.Points2(pts_plane(1:2, :));     

            elseif isa(to_proj, 'sonic.Ellipsoid')
                Ap = to_proj.shapeMat;
                RC2P = to_proj.att;
                Ac = RC2P.dcm'*Ap*RC2P.dcm;
                rc = to_proj.pos;
                rc = rc.r3;
                C = Ac*(rc*rc')*Ac-(rc'*Ac*rc-1)*Ac;

                if ~proj_behind_cam && rc(3) < 0
                    did_proj = false;
                    proj_obj = sonic.Conic.empty;
                else
                    did_proj = true;
                    proj_obj = sonic.Conic(C,'locus');
                end

            else
                error('sonic:PerspectiveProject:pinhole:invalidType', ...
                    ['Invalid type encountered while trying to '...
                     'project. Must be a GeometryP3 object or Ellipsoid.']);
            end
        
        end

        function [proj_obj, did_proj] = stereographic(to_proj, proj_behind_cam)
        %% [proj_obj, did_proj] = stereographic(to_proj, proj_behind_cam)
        %   Performs a stereographic projection of 3D points/lines/quadrics 
        %   into a 2D image plane. Assumes a camera origin at z = -1, with
        %   an image plane lying at z = 0.
        %   
        %   NOTE: Currently only supports points.
        %
        %   Inputs:
        %       - to_proj (1x1 sonic.Points3 OR sonic.PointsS2): Object to
        %         project. Currently only supports projection of points,
        %         which must either be 3D points or points lying on the
        %         sphere.
        %       - proj_behind_cam (1x1 logical): Flag indicating whether
        %         to project points that lie behind the plane of the camera.
        %         OPTIONAL: defaults to FALSE. 
        %
        %   Outputs:
        %       - proj_obj (1x1 sonic.Points2): Set of projected points,
        %         lying on the image plane at z = 0.
        %       - did_proj (1xn logical): Mapping of which points were
        %         successfully projected. If proj_behind_cam was set to false
        %         (as it is by default), only points in front of the camera
        %         will be projected. Thus, it will be a 1xn vector (where n
        %         is the number of points contained in to_proj), with m true
        %         entries, where m points lie in front of the camera. 
        %
        %   Last revised: 03/06/24
        %   Last author: Michael Krause

            arguments
                to_proj             (1, 1)      % Eventually make this sonic.GeometryP3
                proj_behind_cam     (1, 1)      logical     = false
            end

            if isa(to_proj, 'sonic.Points3') || isa(to_proj, 'sonic.PointsS2')

                % Depending on input type, need to grab the R3 repr 
                % differently. Rest of algorithm is the same though.
                if isa(to_proj, 'sonic.Points3')
                    % Pull out all finite points.
                    pts_raw = to_proj.r3_finite;
                else % Else is PointsS2
                    % Pull out unit vec repr of points:
                    pts_raw = to_proj.u;
                end

                if proj_behind_cam
                    % Allow all points to project:
                    pts_to_proj = pts_raw;
                    did_proj = true(1, size(pts_raw, 2));
                else
                    % Filter to those with z > -1:
                    did_proj = pts_raw(3, :) > -1;
                    pts_to_proj = pts_raw(:, did_proj);
                end

                % Divide by 1 + z:
                pts_plane = pts_to_proj./(1 + pts_to_proj(3, :));

                % Then insert back into a Points2 object.
                proj_obj = sonic.Points2(pts_plane(1:2, :));     

            else
                error('sonic:PerspectiveProject:stereographic:invalidType', ...
                    ['Invalid type encountered while trying to '...
                     'project. Must be Points3 or PointsS2 object.']);
            end
        
        end

        function [rot_obj] = rotate(to_rot, rot_att)
        %% [rot_obj] = rotate(to_rot, rot_att)
        %
        %   Rotates an object relative to the observer. Currently only
        %   supports rotating points in P3, i.e., sonic.Points3 and
        %   sonic.PointsS2.
        %
        %   Inputs:
        %       - to_rot (1x1 sonic.Points3 or sonic.PointsS2): Object to
        %         rotate. Currently only supports rotating 3D points.
        %       - rot_att (1x1 sonic.Attitude): Object to Camera attitude. 
        %
        %   Outputs:
        %       - rot_obj (1x1 sonic.Points3 or sonic.PointsS2): Rotated
        %         object. Return type will match the input type.
        %
        %   Last revised: 03/06/24
        %   Last author: Michael Krause
            
            arguments
                to_rot      (1, 1)  % Make this sonic.GeometryP3 eventually
                rot_att     (1, 1)      sonic.Attitude            
            end

            % For now just support points.
            if isa(to_rot, 'sonic.PointsS2') || isa(to_rot, 'sonic.Points3')
                    
                % Pull out the points depending on the object:
                if isa(to_rot, 'sonic.Points3')
                    raw_pts = to_rot.r3;
                else % is PointsS2
                    raw_pts = to_rot.u;
                end

                % Do the actual rotation:
                rot_pts = rot_att.dcm*raw_pts;

                % Package up result commensurate with input type:
                if isa(to_rot, 'sonic.Points3')
                    rot_obj = sonic.Points3(rot_pts);
                else % is PointsS2
                    rot_obj = sonic.PointsS2(rot_pts);
                end

            else
                error('sonic:PerspectiveProject:rotate:invalidType', ...
                        ['Attempted to rotate invalid type. Valid ' ...
                        'types are sonic.Points3 or sonic.PointsS2.']);
            end
        
        end

        function [cropped_obj, crop_map] = crop(to_crop, hfov_RAD, vfov_RAD, dist_model)
        %% [cropped_obj, crop_map] = crop(to_crop, hfov_RAD, vfov_RAD, dist_model)
        %
        %   Crops a set of points lying in P2/R2 based on the specified
        %   horizontal and vertical FOV values. Optionally considers the
        %   effects of a distortion model when calculating the bounds to
        %   crop. 
        %
        %   Inputs:
        %       - to_crop (1x1 sonic.Points2): 2D scene to crop. Currently
        %         only supports cropping a set of 2D points. 
        %       - hfov_RAD (1x1 double): Horizontal field-of-view, in
        %         radians. This is a full-angle FOV.
        %       - vfov_RAD (1x1 double): Vertical field-of-view, in
        %         radians. This is a full-angle FOV.
        %       - dist_model (1x1 sonic.DistortionModel): OPTIONAL, if
        %         specified will consider the effect of the distortion model
        %         when calculating the cropping bounds. 
        %
        %   Outputs:
        %       - cropped_obj (1x1 sonic.Points2): Set of points which were
        %         contained within the calculated cropping bounds. 
        %       - crop_map (1xn logical): Mapping of which points were
        %         successfully projected. It will be a 1xn vector (where n
        %         is the number of points contained in to_crop), with m true
        %         entries, where m points lie within the calculated cropping
        %         bounds. 
        %
        %   Last revised: 03/15/24
        %   Last author: Michael Krause
            
            arguments
                to_crop     (1, 1)  sonic.Points2
                hfov_RAD    (1, 1)  double
                vfov_RAD    (1, 1)  double
                dist_model  (1, 1)  sonic.DistortionModel = sonic.Pinhole()
            end

            switch class(dist_model)
                case 'sonic.BrownConrady'

                    % With the BC model, we know that we'll observe
                    % either pincushion or barrel distortion, and
                    % that'll essentially be decided by the sign of the
                    % first nonzero k coefficient. This will affect
                    % how we calculate our cropping bounds.

                    % If leading k is negative, take the midpoints of 
                    % each side, distort them, and use those for the 
                    % crop bounds. (This is barrel distortion).

                    % If leading k is positive, take the corners,
                    % distort them, and use those for the crop bounds.
                    % (This is pincushion distortion).

                    % Grab the signs of each k term:
                    k_signs = sign([dist_model.k1, dist_model.k2, dist_model.k3]);
                    
                    % Find the sign of the first non-zero k term:
                    major_sign = 0;
                    for idx = 1:length(k_signs)
                        if k_signs(idx) ~= 0
                            major_sign = k_signs(idx);
                            break;
                        end
                    end

                    switch major_sign
                        case -1
                            % Distort midpoints of edges to get bounds
                            xmax_raw = sin(hfov_RAD/2);
                            ymax_raw = sin(vfov_RAD/2);
                            
                            dist_bnds = dist_model.undistort(sonic.Points2(...
                                [xmax_raw,        0, -xmax_raw,         0; ...
                                 0       , ymax_raw,         0, -ymax_raw] ...
                            ));

                            xmax = dist_bnds.r2(1, 1);
                            ymax = dist_bnds.r2(2, 2);
                            xmin = dist_bnds.r2(1, 3);
                            ymin = dist_bnds.r2(2, 4);

                        case 0
                            % Just normal bounds
                            xmax = sin(hfov_RAD/2);
                            ymax = sin(vfov_RAD/2);
                            xmin = -xmax;
                            ymin = -ymax;
                        case 1
                            % Distort corners to get bounds
                            xmax_raw = sin(hfov_RAD/2);
                            ymax_raw = sin(vfov_RAD/2);
                            
                            dist_bnds = dist_model.undistort(sonic.Points2(...
                                [xmax_raw,  xmax_raw, -xmax_raw, -xmax_raw; ...
                                 ymax_raw, -ymax_raw, -ymax_raw,  ymax_raw] ...
                            ));

                            xmax = max(dist_bnds.r2(1, 1:2));
                            xmin = min(dist_bnds.r2(1, 3:4));
                            ymax = max(dist_bnds.r2(2, [1, 4]));
                            ymin = min(dist_bnds.r2(2, [2, 3]));

                        otherwise
                            error('sonic:PerspectiveProject:crop:invalidSign', ...
                                ['Invalid sign detected while attempting ' ...
                                 'to calculate cropping bounds with ' ...
                                 'a distortion model. This should not ' ...
                                 'happen, make sure there are no nans.']);
                    end
                    
                case 'sonic.Pinhole'
                    % There's no distortion, no need to adjust.
                    xmax = sin(hfov_RAD/2);
                    ymax = sin(vfov_RAD/2);
                    xmin = -xmax;
                    ymin = -ymax;
                otherwise
                    error('sonic:PerspectiveProject:crop:badDistortionModel', ...
                        ['Invalid distortion model used to calculate ' ...
                        'crop bounds. Must be a subclass of ' ...
                        'sonic.DistortionModel.']);
            end

            % Check x and y points, then combine into a single `valid`
            % vector:
            pts_x_valid = to_crop.r2(1, :) > xmin & to_crop.r2(1, :) < xmax;
            pts_y_valid = to_crop.r2(2, :) > ymin & to_crop.r2(2, :) < ymax;
            crop_map = pts_x_valid & pts_y_valid;

            % Do the cropping:
            cropped_obj = sonic.Points2(to_crop.r2(:, crop_map));
            
        end

    end

end