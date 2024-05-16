classdef Image

    properties (SetAccess = protected)
        DNmat       double 
        rows        double
        cols        double
    end 

    methods

        function obj = Image(DNmat)
        %% obj = Image(DNmat)
        %   Constructs an image object that contains the digital number
        %   values and image size
        %   
        %   Inputs:
        %       - DNmat  (nxm double): matrix containing the digital number 
        %         values of each pixel in an image 
        %
        %   Outputs:
        %       - obj (1x1 sonic.Image): Image object, encoding an image
        %
        %   Last revised: 2/29/24
        %   Last author: Ava Thrasher
            arguments
                DNmat   double
            end

            obj.DNmat = DNmat;
            
            % get and store size of image
            [obj.rows, obj.cols] = size(DNmat); 
        end

        function flatImgObj = flattenImage(obj,med_area)
        %% flatImgObj = flattenImage(obj,med_area)
        %   Flattens an image by subtracting off the MATLAB median filter 
        %   of the image 
        %   
        %   Inputs:
        %       - obj (1x1 sonic.Image): matrix containing the digital 
        %         number values of each pixel in an image 
        %       - med_area (1x1 double): optional input neighborhood size
        %         to calculate median
        %
        %   Outputs:
        %       - flatImgObj (1x1 sonic.Image): Image object of the 
        %         flattened image
        %
        %   Last revised: 2/29/24
        %   Last author: Ava Thrasher
            arguments
                obj         (1,1) sonic.Image
                med_area    (1,1) double        = 25
            end

            % get digital number matrix
            img = obj.DNmat;

            % median filtering value subtraction
            med_filter = medfilt2(img,[med_area med_area]);
            flat_img = img - med_filter;
            
            % store in img object
            flatImgObj = sonic.Image(flat_img);
        end

        function noiseThresh = estNoiseRand(obj,sigMultiple,randNum)
        %% noiseThresh = estNoiseRand(obj,sigMultiple,randNum)
        % Estimates the background noise of an image by sampling random
        % pixels
        %   
        %   Inputs:
        %       - obj (1x1 sonic.Image): Sonic image object of which the 
        %         background noise will be estimated
        %       - sigMultiple (1x1 double): optional input that will adjust
        %         how many sigma multiples are used    
        %       - randNum (1x1 double): optional total number of random
        %         pixel pairs selected for estimation
        %
        %   Outputs:
        %       - noiseThresh (1x1 double): digital number cut off for
        %         background noise
        %
        %   Last revised: 2/29/24
        %   Last author: Ava Thrasher
            arguments
                obj         (1,1) sonic.Image
                sigMultiple (1,1) double = 3
                randNum     (1,1) double = 100
            end
            % get DN of image
            img = obj.DNmat;
            % get size of image
            n = obj.rows;
            m = obj.cols;
            pixels = n*m;

            % get random set of pixels
            rand_1 = randi(pixels,randNum,1);
            rand_2 = randi(pixels,randNum,1);
            rand_px_1 = img(rand_1);
            rand_px_2 = img(rand_2);

            % for each pair, take difference of brightness then take stdd dev/sqrt2
            diff = rand_px_1-rand_px_2;
            sig_nf = std(diff)/sqrt(2);
            noiseThresh = sigMultiple*sig_nf;
        end

        function noiseThresh = estNoiseSort(obj,sigMultiple,cutOff)
        %% noiseThresh = estNoiseSort(obj,sigMultiple,cutOff)
        % Estimates the background noise of an image by sorting the digital
        % numbers and determining a percentage DN cutoff
        %   
        %   Inputs:
        %       - obj (1x1 sonic.Image): Sonic image object of which the 
        %         background noise will be estimated
        %       - sigMultiple (1x1 double): Optional sigma multiplier 
        %       - cutOff (1x1 double): Optional cutoff percentage ranging
        %         between 0 and 1 that determines the upper bound for the
        %         sorted pixels to be considered
        %
        %   Outputs:
        %       - noiseThresh (1x1 double): digital number cut off for
        %         background noise
        %
        %   Last revised: 2/29/24
        %   Last author: Ava Thrasher
            arguments
                obj             (1,1)   sonic.Image
                sigMultiple     (1,1)   double = 3
                cutOff          (1,1)   double = 0.8
            end
            % get DN of image
            img = obj.DNmat;

            % create vector and sort all pixels based on brightness
            sortedDNs = sort(img(:));

            % sort ascending cut off at percentage
            noise = sortedDNs(1:round(length(sortedDNs)*cutOff));

            % determine standard deviation
            noiseThresh = std(noise)*sigMultiple;  
    
        end

        function DN = bilinInterp(obj,pts_obj)
        %% DN = bilinInterp(obj,pts_obj)
        %   Inputs
        %       - obj (1,1 sonic.Image): Image to sample
        %       - pts_obj (1,1 sonic.Points2): Contains the u,v coordinates
        %         to sample in the image (u's in row 1, v's in row 2)
        %         where the first pixel is centered at 1,1
        %   Outputs
        %       - DN (n,1 double): Sampled pixel values
        %
        %   Last revised: 04/10/24
        %   Last author: Ava Thrasher

            % Get image matrix
            img = obj.DNmat;
            
            % Get all the u,v coordinates
            pts = pts_obj.r2;
            
            % Initialize space for sampled pixel values
            DN = zeros(size(pts,2),1);

            % For each point in pt object
            for i = 1:size(pts,2)
                % Get current point
                u = pts(1,i);
                v = pts(2,i);

                % Get nearest neighbors
                u1 = floor(u);
                u2 = ceil(u);
                v1 = floor(v);
                v2 = ceil(v);

                % Get pixel values
                f11 = img(v1,u1);
                f12 = img(v1,u2);
                f21 = img(v2,u1);
                f22 = img(v2,u2);

                % Thresholding to see if linear interpolation is fine
                thresh = sonic.Tolerances.bilinInterpThresh;
                if abs(u - round(u)) <= thresh && abs(v - round(v)) <= thresh
                    DN_i = img(round(v), round(u));

                elseif abs(u - round(u)) <= thresh
                    % Calculate weight between v values
                    w1 = (v2 - v);
                    w2 = (v - v1);

                    % Calculate sample
                    DN_i = w1*f11 + w2*f21;

                elseif abs(v - round(v)) <= thresh
                    % Calculate weight between u values
                    w1 = (u2 - u);
                    w2 = (u - u1);

                    % Calculate sample
                    DN_i = w1*f11 + w2*f12;

                else 
                    % Calculate weights
                    w11 = (u2 - u)*(v2 - v);
                    w12 = (u - u1)*(v2 - v);
                    w21 = (u2 - u)*(v - v1);
                    w22 = (u - u1)*(v - v1);
    
                    % Calculate sample
                    DN_i = w11*f11 + w12*f12 + w21*f21 + w22*f22;
                end

                % Store sample
                DN(i) = DN_i;
            end

        end

        function [pxVals, pxLocs, lineLengths] = getLine(obj, scan_lines)
        %% [pxVals, pxLocs, lineLengths] = getLine(obj, scan_lines)
        % Given a line object in uv space, getLine will find pixel 
        % values along that line. This assumes 1,1 is centered in the 
        % image's top left pixel. 
        %
        % Inputs:
        %   - obj (1,1 sonic.Image): Image object to find pixel values in
        %   - scanlines (1,1 sonic.ScanLines): ScanLines object containing
        %     the lines and direction for the lines to get
        %
        % Outputs:
        %   - pxVals (1,n cell): cell array containing the pixel
        %     values along the n lines specified. Order of cell array
        %     will be the same as the input line order.
        %   - pxLocs (1,n cell): cell array containing Points2 objects
        %     of each [u,v] coordinate sampled for each line. The order
        %     of vector will be the same as the input line order.
        %   - lineLengths (n,1 double): contains length of each line in the
        %     image in the same order that lines were provided
            arguments
                obj (1,1) sonic.Image
                scan_lines (1,1) sonic.ScanLines
            end

            % Get size of image
            umax = obj.cols;
            vmax = obj.rows;

            % Get lines and directions
            lines = scan_lines.lines;
            ispos = scan_lines.ispos;

            % Define image bounds as lines
            ustart = sonic.Lines2([1;0;-1]);
            uend = sonic.Lines2([1;0;-umax]);
            vstart = sonic.Lines2([0;1;-1]);
            vend = sonic.Lines2([0;1;-vmax]);

            % Storage space for sampled values and locations
            lineLengths = zeros(lines.n,1);
            pxVals = cell(1,lines.n);
            pxLocs = cell(1,lines.n);      

            % For each line, get the DN values
            for i = 1:lines.n
                % Get angle and distance
                ang = lines.angle(i);
                rho = lines.perp_dist(i);

                % Get line slope from angle 
                u_slope = [sin(ang); -cos(ang)];

                % Find line intersection with image borders
                line_i = sonic.Lines2(ang,rho);
                uMeet = sonic.GeometryP2.meet(ustart,line_i);
                umaxMeet = sonic.GeometryP2.meet(uend,line_i);
                vMeet = sonic.GeometryP2.meet(vstart,line_i);
                vmaxMeet = sonic.GeometryP2.meet(vend,line_i);

                % Check which are within image bounds
                active_check = [uMeet.r2(2) <= vmax && uMeet.r2(2) >= 1;
                                umaxMeet.r2(2) <= vmax && umaxMeet.r2(2) >= 1;
                                vMeet.r2(1) <= umax && vMeet.r2(1) >= 1;
                                vmaxMeet.r2(1) <= umax && vmaxMeet.r2(1) >= 1;];

                % If no valid intersections are found, skip to next line
                if sum(active_check) < 2
                    pxVals{i} = [];
                    pxLocs{i} = sonic.Points2.empty(0,1);
                    continue
                end

                % Get edge points
                intersects = [uMeet; umaxMeet; vMeet; vmaxMeet];
                edgePts = intersects(active_check);
                
                % Get initial point
                u0 = edgePts(1);

                % Get distance between active edge points
                ustop = edgePts(2);
                d = sqrt((ustop.r2(1) - u0.r2(1))^2 + (ustop.r2(2) - u0.r2(2))^2);
                % Step by length of one pixel
                t = 0:1:d;

                % Order sample points based on direction flag
                dirFlag = ispos(i);
                if dirFlag
                    % If positive flag, ensure ascending u values
                    startStop = [u0.r2, ustop.r2];
                    [~, ascend] = sort([u0.r2(1), ustop.r2(1)]);
                    startStop = startStop(:,ascend);
                    u = startStop(:,1) + t.*u_slope;
                    % Edge case of vertical lines, sort top to bottom
                    if u0.r2(1) == ustop.r2(1)
                        [~, ascend] = sort([u0.r2(2), ustop.r2(2)],'ascend');
                        startStop = startStop(:,ascend);
                        u = startStop(:,1) + t.*u_slope;
                    end
                elseif ~dirFlag
                    % If negative flag, ensure descending u values
                    startStop = [u0.r2, ustop.r2];
                    [~, descend] = sort([u0.r2(1), ustop.r2(1)],'descend');
                    startStop = startStop(:,descend);
                    u = startStop(:,1) + t.*u_slope;
                    % Edge case of vertical lines, sort bottom to top
                    if u0.r2(1) == ustop.r2(1)
                        [~, descend] = sort([u0.r2(2), ustop.r2(2)],'descend');
                        startStop = startStop(:,descend);
                        u = startStop(:,1) + t.*u_slope;
                    end
                else
                    eid = 'sonic:Image:InvalidDirFlag';
                    msg = 'Direction flags must indicate negative (false) or positive (true) directionality for the line output.';
                    error(eid,msg)
                end
               
                % Store length of line
                lineLengths(i) = size(t,2);

                % Points object 
                samplePts = sonic.Points2(u);

                % Get bilinearly interpolated pixel values
                sampleDN = obj.bilinInterp(samplePts);

                % Save to pxVals and pxLocs
                pxVals{i} = sampleDN;
                pxLocs{i} = samplePts;
            end
        end

        function cleaned_img = removeSmearCCD(obj, t_exp, t_inout)
        %% cleaned_img = removeSmearCCD(obj, t_exp, t_inout)
        % Given the clock cycling time of a CCD focal array and an image 
        % containing smear, this function will remove the smear and output
        % the corrected image.
        %   
        %   Inputs:
        %       - obj (1x1 sonic.Image): Sonic image object which contains
        %         smear to be removed
        %       - t_exp (1x1 double): Image exposure time  
        %       - t_inout (1x1 OR 1x2 double): CCD readin/readout time 
        %         per row. If input as a 1,1 the read in and out time are
        %         assumed to be the same. If input as a 1,2 the first index
        %         is the read in time and the second index is the read out
        %         time.
        %
        %   Outputs:
        %       - cleaned_img (1x1 sonic.Image): Image object with smear
        %         removed
        %
        %   Last revised: 5/14/24
        %   Last author: Ava Thrasher
        arguments
            obj     (1,1) sonic.Image
            t_exp   (1,1) double
            t_inout double 
        end

            % Ensure read in and out times are input correctly
            one_check = all(size(t_inout) == [1,1]);
            two_check = all(size(t_inout) == [1,2]);
            if ~one_check && ~two_check
                eid = 'sonic:Image:removeSmearCCD:invalidReadTime';
                msg = 't_inout must be either a scalar (where read in and out time are assumed to be the same) or a 1x2 vector where the first index is read in time and the second index is read out time.';
                error(eid,msg)
            end

            % Ensure that the image is square
            n = obj.rows;
            m = obj.cols;
            if n~=m
                eid = 'sonic:Image:removeSmearCCD:invalidImageSize';
                msg = 'Image must be square.';
                error(eid,msg)
            end
            
            % Check if in and out times are close enough together to
            % consider them the same
            if two_check
                t_diff = abs(t_inout(2) - t_inout(1));
                if t_diff < eps
                    t_inout = t_inout(1);
                    one_check = true;
                end
            end

            im_raw = obj.DNmat;
            if one_check
                eps_inout = t_inout/t_exp;
                % analytic inversion of M
                Minv = (1/(1-eps_inout)).*(eye(n) + eps_inout/(eps_inout - n*eps_inout - 1).*ones(size(n)));
                im_clean = Minv*im_raw;
                cleaned_img = sonic.Image(im_clean);
            else
                eps_in = t_inout(1)/t_exp;
                eps_out = t_inout(2)/t_exp;
                L = tril(ones(n),-1);
                D = eps_in*L + eps_out*L';
                M = eye(n) + D;
                im_clean = M\im_raw;
                cleaned_img = sonic.Image(im_clean);
            end

        end

    end
end