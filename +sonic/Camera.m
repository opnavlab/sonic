classdef Camera
    
    properties
        d_x             (1, 1)      double
        d_y             (1, 1)      double
        alpha           (1, 1)      double
        u_p             (1, 1)      double
        v_p             (1, 1)      double
        dist_model      (1, 1)      sonic.DistortionModel = sonic.Pinhole
    end

    properties (Dependent)
        K               (3, 3)      double
        Kinv            (3, 3)      double
        hfov_RAD        (1, 1)      double
        vfov_RAD        (1, 1)      double
        ifov_h_RAD      (1, 1)      double
        ifov_v_RAD      (1, 1)      double
    end
    
    % Constructor:
    methods
        function obj = Camera(varargin)
        %% obj = Camera(varargin)
        %   
        %   Instantiates a model for a framing camera. Depending on the
        %   arguments, either invokes `constructFromFOV()` or
        %   `constructFromK()` to complete the model instantiation. See
        %   those methods below for details.
        %   
        %   Inputs:
        %       - varargin: see `constructFromFOV()` and `constructFromK()`.
        %
        %   Outputs:
        %       - obj (sonic.Camera): Camera object, containing camera
        %       intrinsics and a distortion model.
        %
        %   Last revised: 03/15/24
        %   Last author: Michael Krause

            switch nargin
                case 5
                    obj = obj.constructFromFOV(varargin{:});
                case 6
                    obj = obj.constructFromK(varargin{:});
                otherwise
                    error('sonic:Camera:invalidArgs', ...
                        ['Must specify arguments concerning camera ' ...
                        'FOV or camera intrinsics. See the inputs for ' ...
                        '<a href="matlab: help(''sonic.Camera.' ...
                        'constructFromFOV'')">instantiation with ' ...
                        'FOV info</a> or the inputs for <a href="matlab: ' ...
                        'help(''sonic.Camera.constructFromK'')">' ...
                        'instantiation from camera intrinsics</a>.']);
            end

        end
    end

    % Constructor helper methods:
    methods (Access=private)
        
        function obj = constructFromFOV(obj, fov, fov_type, direction, res, dist_model)
        %% obj = constructFromFOV(obj, fov, fov_type, direction, res, dist_model)
        %   
        %   Creates a model for a framing camera using information about 
        %   the camera's FOV or instantaneous FOV, as well as sensor 
        %   resolution. Also consumes a distortion model. Note that x/y 
        %   directions here coorespond with the canonical camera
        %   model, i.e.: +z along boresight, +x to the right looking out of
        %   the camera, +y down looking out of the camera. As well, u/v
        %   cooresponds to x/y directions, but is centered in the top-left
        %   of the image. 
        %
        %   NOTE: When constructing a camera model with this method, square
        %   pixels and no shear are assumed. If this is not the case,
        %   construct the camera model directly from the camera intrinsics.
        %   
        %   Inputs:
        %       - fov (1x1 double): Field-of-view or IFOV, in radians. In
        %       either case, these are full-angle values.
        %       - fov_type (1x1 string): String specifying whether the
        %       previous argument was an FOV ('fov') or IFOV ('ifov').
        %       - direction (1x1 string): String indicating the direction
        %       of the aforementioned FOV/IFOV: 'h' for horizontal
        %       FOV/IFOV, and 'v' for vertical FOV/IFOV.
        %       - res (2 double): Resolution of the sensor, specified as
        %       [# of rows, # of cols]
        %       - dist_model (1x1 sonic.DistortionModel): Distortion model
        %       used to transform projected points prior to being passed
        %       through the camera intrinsics.
        %
        %   Outputs:
        %       - obj (sonic.Camera): Camera object, containing camera
        %       intrinsics and a distortion model.
        %
        %   Last revised: 03/15/24
        %   Last author: Michael Krause

            arguments
                obj             (1, 1)      sonic.Camera
                fov             (1, 1)      double
                fov_type        (1, 1)      string
                direction       (1, 1)      string
                res             (2, 1)      double      % [rows, cols]
                dist_model      (1, 1)      sonic.DistortionModel
            end

            switch lower(direction)
                case 'h'
                    major_res = res(2);
                case 'v'
                    major_res = res(1);
                otherwise
                    error('sonic:Camera:constructFromFOV:invalidFOVDir', ...
                        ['Must specify FOV direction as one of `h` ' ...
                        '(for horizontal FOV/IFOV) or `v` (for ' ...
                        'vertical FOV/IFOV).']);
            end
            
            switch lower(fov_type)
                case 'fov'
                    obj.d_x = major_res/(2*tan(fov/2));
                    obj.d_y = obj.d_x;
                case 'ifov'
                    obj.d_x = 1/fov;
                    obj.d_y = obj.d_x;
                otherwise
                    error('sonic:Camera:constructFromFOV:invalidFOVFlag', ...
                        ['Must specify FOV/IFOV type with an ' ...
                        'appropriate flag: `fov` for FOV, ' ...
                        '`ifov` for IFOV. ']);
            end

            obj.u_p = res(2)/2;
            obj.v_p = res(1)/2;
            obj.alpha = 0;

            obj.dist_model = dist_model;

        end

        function obj = constructFromK(obj, d_x, d_y, alpha, u_p, v_p, dist_model)
        %% obj = constructFromK(obj, d_x, d_y, alpha, u_p, v_p, dist_model)
        %   
        %   Creates a model for a framing camera using elements of the
        %   camera calibration matrix. Also consumes a distortion model. 
        %   Note that x/y directions here coorespond with the canonical 
        %   camera model, i.e.: +z along boresight, +x to the right looking
        %   out of the camera, +y down looking out of the camera. As well, 
        %   u/v cooresponds to x/y directions, but is centered in the 
        %   top-left of the image.
        %   
        %   Inputs:
        %       - d_x (1x1 double): Focal length/pixel pitch, x direction
        %       - d_y (1x1 double): Focal length/pixel pitch, y direction
        %       - alpha (1x1 double): Pixel shear 
        %       - u_p (1x1 double): Center point, u coordinates (pixels)
        %       - v_p (1x1 double): Center point, v coordinates (pixels)
        %       - dist_model (1x1 sonic.DistortionModel): Distortion model
        %       used to transform projected points prior to being passed
        %       through the camera intrinsics.
        %
        %   Outputs:
        %       - obj (sonic.Camera): Camera object, containing camera
        %       intrinsics and a distortion model.
        %
        %   Last revised: 03/15/24
        %   Last author: Michael Krause
            
            arguments
                obj             (1, 1)      sonic.Camera
                d_x             (1, 1)      double
                d_y             (1, 1)      double
                alpha           (1, 1)      double
                u_p             (1, 1)      double
                v_p             (1, 1)      double
                dist_model      (1, 1)      sonic.DistortionModel
            end

            % Just assign the inputs to the object:
            obj.d_x = d_x;
            obj.d_y = d_y;
            obj.alpha = alpha;
            obj.u_p = u_p;
            obj.v_p = v_p;

            obj.dist_model = dist_model;
        
        end

    end

    % Dependent property methods:
    methods

        function val = get.K(obj)
        %% val = get.K(obj)
        %   
        %   Returns the camera intrinsics as a 3x3 matrix.
        %   
        %   Inputs:
        %       - obj (1x1 sonic.Camera): Camera object, containing camera
        %       intrinsics and a distortion model.
        %
        %   Outputs:
        %       - val (3x3 double): Camera intrinsics matrix, such that:
        %                                       [ d_x, alpha, u_p;
        %                                   K =     0,   d_y, v_p;
        %                                           0,     0,   1];
        %
        %   Last revised: 03/06/24
        %   Last author: Michael Krause

            val = [obj.d_x, obj.alpha, obj.u_p
                        0,    obj.d_y, obj.v_p
                        0,         0,       1];
        
        end

        function val = get.Kinv(obj)
        %% val = get.K(obj)
        %   
        %   Returns the inverse of the camera intrinsics as a 3x3 matrix.
        %   Note that this inverse is calculated analytically, not
        %   numerically. 
        %   
        %   Inputs:
        %       - obj (1x1 sonic.Camera): Camera object, containing camera
        %       intrinsics and a distortion model.
        %
        %   Outputs:
        %       - val (3x3 double): Inverse of the camera intrinsics 
        %       matrix (see `K` for details).
        %
        %   Last revised: 03/13/24
        %   Last author: Michael Krause
        
            dxdy = obj.d_x*obj.d_y;
        
            val = [1/obj.d_x, -obj.alpha/dxdy, (obj.alpha*obj.v_p - obj.d_y*obj.u_p)/dxdy
                   0,               1/obj.d_y, -obj.v_p/obj.d_y
                   0,                       0, 1];
            
        end

        function val = get.hfov_RAD(obj)
        %% val = get.hfov_RAD(obj)
        %   
        %   Returns the camera's horizontal field-of-view (FOV) in radians.
        %   This is the full-angle FOV.
        %   
        %   Inputs:
        %       - obj (1x1 sonic.Camera): Camera object, containing camera
        %       intrinsics and a distortion model.
        %
        %   Outputs:
        %       - val (1x1 double): Horizontal field-of-view of the camera,
        %       in radians.
        %
        %   Last revised: 03/06/24
        %   Last author: Michael Krause

            val = 2*atan(obj.u_p/obj.d_x);
        end

        function val = get.vfov_RAD(obj)
        %% val = get.vfov_RAD(obj)
        %   
        %   Returns the camera's vertical field-of-view (FOV) in radians.
        %   This is the full-angle FOV. 
        %   
        %   Inputs:
        %       - obj (1x1 sonic.Camera): Camera object, containing camera
        %       intrinsics and a distortion model.
        %
        %   Outputs:
        %       - val (1x1 double): Vertical field-of-view of the camera,
        %       in radians.
        %
        %   Last revised: 03/06/24
        %   Last author: Michael Krause

            val = 2*atan(obj.v_p/obj.d_y);
        end

        function val = get.ifov_h_RAD(obj)
        %% val = get.ifov_h_RAD(obj)
        %   
        %   Returns the camera's horizontal instantaneous field-of-view 
        %   (IFOV) in radians.
        %   
        %   Inputs:
        %       - obj (1x1 sonic.Camera): Camera object, containing camera
        %       intrinsics and a distortion model.
        %
        %   Outputs:
        %       - val (1x1 double): Horizontal instantaneous field-of-view 
        %       (IFOV) of the camera, in radians.
        %
        %   Last revised: 03/13/24
        %   Last author: Michael Krause

            val = 1/obj.d_x;     % IFOV_x
        
        end

        function val = get.ifov_v_RAD(obj)
        %% val = get.ifov_v_RAD(obj)
        %   
        %   Returns the camera's vertical instantaneous field-of-view 
        %   (IFOV) in radians.
        %   
        %   Inputs:
        %       - obj (1x1 sonic.Camera): Camera object, containing camera
        %       intrinsics and a distortion model.
        %
        %   Outputs:
        %       - val (1x1 double): Vertical instantaneous field-of-view 
        %       (IFOV) of the camera, in radians.
        %
        %   Last revised: 03/13/24
        %   Last author: Michael Krause

            val = 1/obj.d_y;     % IFOV_y
        
        end
    end

    % Public methods:
    methods 

        function [proj_px, proj_map] = synthImage(obj, to_proj, attitude, velocity, opts)
        %% [proj_px, proj_map] = synthImage(obj, to_proj, attitude, velocity, opts)
        %   
        %   Given 3D geometry (currently only supports 3D points, i.e.,
        %   Points3 or PointsS2), projects this geometry into a synthetic
        %   image per the camera model.
        %   
        %   Inputs:
        %       - obj (1x1 sonic.Camera): Camera object, containing camera
        %       intrinsics and a distortion model.
        %       - to_proj (1x1 sonic.Points3 or sonic.PointsS2): Points to
        %       project into an image.
        %       - attitude (1x1 sonic.Attitude): Attitude of the object
        %       relative to the camera
        %       - velocity (3x1 double): OPTIONAL: Velocity of the camera 
        %       in the inertial frame. Defaults to zero if not supplied.
        %       - opts.projectBehindCamera (1x1 logical): OPTIONAL: 
        %       defaults to false. If true, also projects points that lie 
        %       behind the camera. See sonic.Project for more details. 
        %       - opts.cropFOV (2 double or []): OPTIONAL: If
        %       specified, allows the user to override the FOV of
        %       the camera when cropping the projected image. Can specify
        %       as a 2 element vector, [hfov_RAD, vfov_RAD], or [] to not
        %       crop at all. Defaults to the HFOV/VFOV of the camera. 
        %
        %   Outputs:
        %       - proj_px (1x1 sonic.Points2): 2D points lying in the image
        %       plane, in pixel coordinates. Represents the projection of
        %       the provided points per the camera model.
        %       - proj_map (1xn logical): Given that `to_proj` contains n
        %       points, this is 1-by-n logical vector containing m true
        %       values (where m points were successfully projected into the
        %       resultant image), corresponding to each index of a 
        %       successful projection. 
        %
        %   Last revised: 03/15/24
        %   Last author: Michael Krause

            arguments
                obj                         (1, 1)  sonic.Camera
                to_proj                     (1, 1)  % Eventually, sonic.GeometryP3
                attitude                    (1, 1)  sonic.Attitude
                velocity                    (3, 1)  double      = [0;0;0];
                opts.projectBehindCamera    (1, 1)  logical     = false
                opts.cropFOV                (:, 1)  double      = [obj.hfov_RAD; obj.vfov_RAD]
            end

            % Rotate the points in the scene per our camera attitude:
            to_proj_rot = sonic.Project.rotate(to_proj, attitude);

            % Do a pinhole projection:
            [proj_raw, did_proj] = ...
                sonic.Project.pinhole(to_proj_rot, opts.projectBehindCamera);

            % Crop to just capture the relevant points:
            if ~(isempty(opts.cropFOV) || length(opts.cropFOV) ~= 2)
                [proj_crop, did_crop] = sonic.Project.crop(proj_raw, ...
                    opts.cropFOV(1), opts.cropFOV(2), obj.dist_model);
            else
                proj_crop = proj_raw;
                did_crop = true(1, proj_crop.n);
            end

            % Aberrate
            proj_ab = sonic.Aberration.aberrate(proj_crop, velocity);

            % Distort:
            proj_dist = obj.dist_model.distort(proj_ab);

            % Put it through the camera intrinsics:
            proj_px = sonic.Points2(obj.K*proj_dist.p2);

            % If we're not doing any custom cropping, let's do one final
            % pixel-space crop to ensure we're in the image:
            if ~(isempty(opts.cropFOV) || length(opts.cropFOV) ~= 2 || ...
                any(opts.cropFOV ~= [obj.hfov_RAD; obj.vfov_RAD]))
            
                % Pixel coords are at the center of the pixel, so
                % technically (since we have subpixel points), we need to
                % evaluate bounds at the 1/2 pix offset point.
                x_valid = proj_px.r2(1, :) > -0.5 & ...
                    proj_px.r2(1, :) < (2*obj.u_p)+0.5;

                y_valid = proj_px.r2(2, :) > -0.5 & ...
                    proj_px.r2(2, :) < (2*obj.u_p)+0.5;

                did_px_crop = x_valid & y_valid;
                proj_px = sonic.Points2(proj_px.r2(:, did_px_crop));

                % Calculate the projection map:
                crop_map = did_crop;
                crop_map(did_crop) = did_px_crop;
                proj_map = did_proj;
                proj_map(did_proj) = crop_map;

            else
                % Calculate the projection map:
                proj_map = did_proj;
                proj_map(did_proj) = did_crop;
            end            

        end
        
    end
end

