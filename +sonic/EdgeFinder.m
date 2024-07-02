% This software is made available under the MIT license. See SONIC/LICENSE
% for details.
classdef EdgeFinder

    methods (Static)
        
        function [subPixelEdge, subPixelNormals] = getSubPixel(im_obj,pixelEdge,BlurWidth,kernel_type,M)
        %% [subPixelEdge, subPixelNormals] = getSubPixel(im_obj,pixelEdge,BlurWidth,kernel_type,M)
        % The getSubPixel function computes the subpixel edge 
        % location and orientation from an image, pixel-level edge map, and 
        % information about the image blur. This function is intended for 
        % use on ISOLATED EDGES (that is, edges that are not within about 
        % 5 pixels of another edge). The algorithm is based on the method 
        % published in:
        %
        %    Christian, J.A., "Accurate Planetary Limb Localization for Image-Based 
        %    Spacecraft Navigation," Journal of Spacecraft and Rockets, Vol. 54, 
        %    No. 3, 2017, pp 708-730. doi: 10.2514/1.A33692.
        %    Online: https://arc.aiaa.org/doi/abs/10.2514/1.A33692 
        %
        %   Inputs:
        %       - im_obj (1x1 sonic.Image): grayscale image containing the
        %         body on which to find limb points
        %       - pixelEdge (1,1 sonic.Points2): pixel level points in the
        %         image where edges are found. Must be rounded to nearest 
        %         pixel. 
        %       - BlurWidth (1x1 double) width of Gaussian blur in image
        %       - kernel_type (1x1 double) flag for type of Zernike moment 
        %         kernels. Input 5 for 5x5 moment or 7 for 7x7 moment.
        %       - M (nxn logical): optional forground mask to only include
        %         a masked part of the image in the search
        %
        %   Outputs:
        %       - subPixelEdge (1x1 sonic.Points2): edge points 
        %       - subPixelNormals (1x1 sonic.Points2): normal of edge
        %         points, corresponding to edge points by index
        %
        % The edge normals are defined positive going from background (darker
        % region) to foreground (brighter region).
        % 
        %   Last revised: 3/28/2024
        %   Author: Dr. John A. Christian
        %   Edited: Ava Thrasher
        
            arguments
                im_obj (1,1) sonic.Image
                pixelEdge (1,1) sonic.Points2
                BlurWidth (1,1) double
                kernel_type (1,1) double {sonic.EdgeFinder.mustBe5or7} = 5
                M = []
            end
                  
            % Get image matrix
            img = im_obj.DNmat;

            % Get image size
            nrows = im_obj.rows;
            ncols = im_obj.cols;
            
            % Create mask where points are in the image
            E = false(nrows, ncols);
            for idx = 1:pixelEdge.n
                 E(pixelEdge.r2(2, idx), pixelEdge.r2(1, idx)) = true;
            end

            % Check input image size
            if nrows<5 || ncols<5
                error('sonic:EdgeFinder:getSubPixel:invalidImageInput', ...
                            'Input image (img) must be at least 5x5.')
            end
                                             
            if BlurWidth > 1.5 && kernel_type == 5
                error('sonic:EdgeFinder:getSubPixel:BlurWidthTooLarge',...
                    'Blur width is too large for the 5x5 kernels used here. Need to implement with larger kernels (e.g. 7x7).')
            end
            
            % Check for optional foreground mask 
            if ~isempty(M)
                
                % Check type and size of the mask
                if ~isa(M,'logical')
                    error('sonic:EdgeFinder:getSubPixel:invalidMType',...
                        'Mask input (M) must be a logical array.')
                end
                
                if ~all( [nrows,ncols] == size(M) )
                    error('sonic:EdgeFinder:getSubPixel:invalidMSize',...
                        'Input image (img) and mask (M) must be arrays of the same size.')
                end
            
                % Zero out a two-pixel border around the image periphery
                M_NoBorder = false(nrows,ncols);
                if kernel_type == 5
                    M_NoBorder(3:nrows-2,3:ncols-2) = M(3:nrows-2,3:ncols-2);
                elseif kernel_type == 7
                    M_NoBorder(4:nrows-3,4:ncols-3) = M(4:nrows-3,4:ncols-3);
                end
                
                % Apply mask to edge locations
                E = E & M_NoBorder;
                
            else
                
                % Zero out a pixel border around the image periphery
                % depending on Zernike kernel size
                E_temp = false(nrows,ncols);
                if kernel_type == 5
                    E_temp(3:nrows-2,3:ncols-2) = E(3:nrows-2,3:ncols-2);
                elseif kernel_type == 7
                    E_temp(4:nrows-3,4:ncols-3) = E(4:nrows-3,4:ncols-3);
                end
                E = E_temp;
                
            end
            
            % 5x5 Zernike moment kernels (see reference paper for details)
            ReM11_55 = sonic.Constants.ReM11_55;
            ImM11_55 = ReM11_55';
            M20_55 = sonic.Constants.M20_55;
            
            % 7x7 Zernike moment kernels
            ReM11_77 = sonic.Constants.ReM11_77;
            ImM11_77 = ReM11_77';
            M20_77 = sonic.Constants.M20_77;
            
            % Perform subpixel edge localization   
            
            % Convert blur width to ramp width, then normalize into unit 
            % circle (divide by 5/2=2.5 since we have a 5x5 kernel).
            if BlurWidth < 1/sqrt(12)
                warning('sonic:EdgeFinder:getSubPixel:BlurWidthTooSmall',...
                    'BlurWidth is less than expected due to pixel quantization. Consider increasing to at least 0.289. Having a blur width that is too low may result in a small bias in edge location (although still usually better than the originl pixel-level result).')
            end
            w = BlurWidth * 1.66 / 2.5;
            
            % Find indices of edges in E
            E_idx = find(E);
            [E_v,E_u]  = ind2sub([nrows,ncols], E_idx); 
            NumEdges = length(E_idx);
            
            % Initialize loop parameters for speed
            EdgesFound = 0;
            u_temp  = zeros(NumEdges,1);
            v_temp  = zeros(NumEdges,1);
            nu_temp = zeros(NumEdges,1);
            nv_temp = zeros(NumEdges,1);
            
            % Loop through all pixel-level edges
            for ii = 1:NumEdges
                
                % Get current row and cloumn
                row_ii = E_v(ii);  %v-direction is vertical (row number)
                col_ii = E_u(ii);  %u-direction is horizontal (column number)
                EdgeGuess = [col_ii; row_ii];  %pixel-level edge guess
                
                if kernel_type == 5
                    % Get 5x5 window centered around current pixel location
                    img_window = img(row_ii-2:row_ii+2,col_ii-2:col_ii+2);
                elseif kernel_type == 7
                    % 7x7 window:
                    img_window = img(row_ii-3:row_ii+3,col_ii-3:col_ii+3);
                else
                    error('sonic:EdgeFinder:getSubPixel:invalidInput',...
                        'Must specify `kernel_type` as 5 or 7.');
                end
            
                % Compute image correlations at current edge location
                if kernel_type == 5
                    CorrImA11 = ImM11_55 .* img_window; 
                    CorrReA11 = ReM11_55 .* img_window;
                    CorrA20   =   M20_55 .* img_window; 
                elseif kernel_type == 7
                    CorrImA11 = ImM11_77 .* img_window; 
                    CorrReA11 = ReM11_77 .* img_window;
                    CorrA20   =   M20_77 .* img_window; 
                end
                ImA11 = sum(CorrImA11(:));
                ReA11 = sum(CorrReA11(:));
                A20   = sum(  CorrA20(:));
                
                % Compute edge normal
                phi = atan2( ImA11, ReA11 );
                n = [cos(phi); sin(phi)];
                
                % Get ratio of moments
                A11_prime = ReA11*cos(phi) + ImA11*sin(phi);
                R = A20 / A11_prime; 
                
                % Compute edge location in unit circle
                if w > 0.001
                    w2 = w^2;
                    l = ( 1 - w2 - sqrt( (w2-1)^2 - 2*w2*R ) ) / w2;
                else
                    l = R;
                end
                
                % Prevent corrections that take you too close to the edge 
                % of the unit circle (bad points are removed)
                if abs(l) < 0.9
                    
                    % Compute subpixel edge location
                    e = EdgeGuess + 2.5*l*n;
                    
                    % Increment edge counter
                    EdgesFound = EdgesFound + 1;
                    
                    % Store edge and normal direction
                    u_temp(EdgesFound,1)  = e(1);
                    v_temp(EdgesFound,1)  = e(2);
                    nu_temp(EdgesFound,1) = n(1);
                    nv_temp(EdgesFound,1) = n(2);
                    
                end 
            end
            
            % Output results to sonic.Points2 object
            if EdgesFound == 0
                subPixelEdge = sonic.Points2.empty();
                subPixelNormals = sonic.Points2.empty();
            else
                subPixelEdge = sonic.Points2([u_temp(1:EdgesFound), v_temp(1:EdgesFound)]');
                subPixelNormals = sonic.Points2([nu_temp(1:EdgesFound), nv_temp(1:EdgesFound)]');
            end
        end

        function E = sobel(im_obj,gradThresh)
        %% E = sobel(im_obj,gradThresh)
        %
        % Computes edge locations based on local extrema in Gu or Gv gradients
        %
        %   Inputs:
        %       - im_obj (1,1 sonic.Image): image object to find edges in
        %       - gradThresh (1,1 double): Threshhold above which a gradient i
        %         considered as a potential edge
        %   Outputs:
        %       - E (nxm logical): mask which is true where edges were detecte
        
            % Get digital number matrix
            im = im_obj.DNmat;
            
            % Find image rows and columns
            nrows = im_obj.rows;
            ncols = im_obj.cols;
            
            % Initalize gradient matrices for u and v directions
            Gu = zeros([nrows,ncols]);
            Gv = Gu;
            
            % Compute gradient with normalized Sobel operator 
            % Accomplished by shifting image instead of convolution 
            % Ignores all pixels around the image border
            Gu(3:end-2,3:end-2) = -(1/8)*im(2:end-3,2:end-3) + ...
                                  -(2/8)*im(3:end-2,2:end-3) + ...
                                  -(1/8)*im(4:end-1,2:end-3) + ...
                                   (1/8)*im(2:end-3,4:end-1) + ...
                                   (2/8)*im(3:end-2,4:end-1) + ...
                                   (1/8)*im(4:end-1,4:end-1);
            
            Gv(3:end-2,3:end-2) = -(1/8)*im(2:end-3,2:end-3) + ...
                                  -(2/8)*im(2:end-3,3:end-2) + ...
                                  -(1/8)*im(2:end-3,4:end-1) + ...
                                   (1/8)*im(4:end-1,2:end-3) + ...
                                   (2/8)*im(4:end-1,3:end-2) + ...
                                   (1/8)*im(4:end-1,4:end-1);                   
            
                               
            % Compute the absolute value of Gu and Gv                   
            Gu = abs(Gu);
            Gv = abs(Gv);
            
            % Compute magnitude of gradient
            % RSS is slow, so we approximate using the absoulte value
            G = Gu + Gv;  
            
            % Find all pixel coordinates where gradient is above the threshold
            check1 = G( 4:end-3, 4:end-3 ) > gradThresh;
            
            % Find graident extrema that are ``mostly'' vertical
            % A vertical extrema has Gu bigger than Gv
            % To be a local exterema, query Gu must be bigger than the Gu to the left and right
            check2 = Gu( 4:end-3, 4:end-3 ) > Gv( 4:end-3, 4:end-3 ) & ...
                     Gu( 4:end-3, 4:end-3 ) > Gu( 4:end-3, 3:end-4 ) & ...
                     Gu( 4:end-3, 4:end-3 ) > Gu( 4:end-3, 5:end-2 ) ;
            
            % Find graident extrema that are ``mostly'' horizontal   
            % A horizontal extrema has Gv bigger than Gu
            % To be a local exterema, query Gv it is also bigger than the Gv above and below
            check3 = Gv( 4:end-3, 4:end-3 ) > Gu( 4:end-3, 4:end-3 ) & ...
                     Gv( 4:end-3, 4:end-3 ) > Gv( 3:end-4, 4:end-3 ) & ...
                     Gv( 4:end-3, 4:end-3 ) > Gv( 5:end-2, 4:end-3 ) ;
            
            % Initalize output as a logical     
            E = false( [nrows,ncols] );
            
            % Edges are those locations where the gradient is strong enough AND there
            % is a local extrema in Gu or Gv
            E( 4:end-3, 4:end-3 ) = check1 & (check2 | check3);
        
        end

        function limb_pts = runScan(im_obj, scan_lines, min_consec, run_start_thresh)
        %% limb_pts = runScan(im_obj, min_consec, run_start_thresh)
        %  Performs a scan across the provided image in the direction of
        %  sun illumination to find candidate edge points.
        % 
        %   Inputs: 
        %       - im_obj (1x1 sonic.Image): SONIC Image object containing 
        %         the limb to find points on.
        %       - scans_lines (1x1 sonic.ScanLines): Lines in the image to
        %         perform scan along.
        %       - min_consec (1x1 double): Minimum number of consecutive
        %         bright pixels found on a scan for the first bright
        %         pixel to be considered a limb point. Used to keep hot
        %         pixels, stars, or debries from being falsely identified
        %         as an edge
        %       - run_start_thresh (1x1 double): Minimum digital number to 
        %         flag as a limb point 
        %
        %   Outputs:
        %       - limb_pts (1x1 sonic.Points2): SONIC Points2 object 
        %         containing the determined u v coordinate limb points. U
        %         in first row V in second row.
        %
        %   Last revised: 3/28/24
        %   Last author: Ava Thrasher
            arguments
                im_obj (1,1) sonic.Image
                scan_lines (1,1) sonic.ScanLines
                min_consec (1,1) double
                run_start_thresh (1,1) double
            end  

            % Get pixel values along lines from image
            [pxVals, pxLocs, lineLengths] = im_obj.getLine(scan_lines);

            % preallocate space for limb points 
            limb_pts = nan(2, length(pxVals));

            %  Scan each line for a potential edge point
            for i = 1:length(pxVals)
                % Get current scan line and pixel coordinate
                pxVals_i = pxVals{i};
                pxLocs_i = pxLocs{i};

                % Iterate over current scan 
                run = 0;
                for j = 1:lineLengths(i)
                    % Get current pixel value
                    pxVal_j = pxVals_i(j);

                    % Check if it is above the threshold
                    if pxVal_j >= run_start_thresh
                        % If yes, increase run
                        run = run + 1;
                    else
                        run = 0;
                    end

                    % Check if minimum run size has been met
                    if run >= min_consec
                        % If yes, save head point as edge and move to next
                        % line
                        limb_pts(:,i) = pxLocs_i.r2(:,j-run);
                        break
                    end

                end
            end

            % Remove nans from the limb points
            limb_pts = limb_pts(:,~isnan(limb_pts(1,:)));
            
            % Store candidate points to sonic.Points2 object
            limb_pts = sonic.Points2(limb_pts);
        end

        function activeBorders = getActiveBorder(u_illum_obj)
        %% activeBorders = getActiveBorder(u_illum_obj)
        %   Determines which image borders are active, meaning which image
        %   borders have light coming in
        %
        %   Inputs:
        %       - u_illum_obj (1x1 sonic.Points2): Sonic Points2 object
        %         containing the uv direction of illumination from the 
        %         sun to the body in the image
        %   Outputs:
        %       - active_mask (1x4 logical): Logical array containing true
        %         for and active edge and false for an inactive edge in 
        %         the clockwise order of [top, right, bottom, left]
            
            % Edge normals defined as positive pointing inwards (MAKE ONE BIG MATRI AND DO TRANPOE MULTIPLE)
            edges = [0 1; -1 0; 0 -1; 1 0];
                     
            % Get illumination direction
            u_illum = u_illum_obj.r2;

            % Determine active edges where n*u_illum is positive
            activeBorders = ((edges*u_illum)'>0)';

        end

    end

    methods (Static,Hidden=true,Access=protected)

        function mustBe5or7(kernel_type)
            % Test for input to be either 5 or 7
            if kernel_type ~= 5 && kernel_type ~= 7
                eid = 'sonic:EdgeFinder:InvalidKernelType';
                msg = 'Zernike kernel type must be either 5 or 7.';
                error(eid,msg)
            end
        end

    end

end