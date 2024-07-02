% This software is made available under the MIT license. See SONIC/LICENSE
% for details.
classdef Centroider

    methods (Static)

        function [centroids,comps] = COB(imgObj, thresh, options) 
        %% centroids = COB(imgObj, thresh, options) 
        % Estimate centers of brightness within an image by removing 
        % background noise, and using the center of brightness equation
        %   
        %   Inputs:
        %       - imgObj (1x1 sonic.Image): Sonic image to detect
        %         centroids in
        %       - thresh (1x1 double): threshold value for noise usually
        %         determined by running img.estNoiseRand or 
        %         img.estNoiseSort
        %       - minStarSize (1x1 double): Optional cutoff for minimum
        %         cluster size. Good for removing hot pixels.
        %       - maxStarSize (1x1 double): Optional cutoff for maximum
        %         cluster size.
        %       - nonMaxSupDist (1x1 double): Optional minimum distance 
        %         between cluster centroids. Will always save centroids
        %         with largest DN within this area. Used for non maximum
        %         suppression.
        %       - buffer (1x1 double): Optional no centroids found to be
        %         within this many pixels from the edge will be 
        %         considered. Often used when a median filter has been 
        %         applied.
        %
        %   Outputs: 
        %       - centroids (1x1 sonic.Points2 or empty 0x0 sonic.Points2): 
        %         Points2 object containing the centroids found in the 
        %         image. If no centroids are detected, an empty Points2 
        %         will be output.
        %       - comps (nx1 cell array): Contains the linear indices of
        %         the connected components corresponding to the centroids
        %
        %   Last revised: 3/7/24
        %   Last author: Ava Thrasher
            arguments %% ENFORCE POSITIVE INPUT %%
                imgObj                  (1,1)   sonic.Image
                thresh                  (1,1)   double
                options.minStarSize     (1,1)   double {mustBePositive}      = 9
                options.maxStarSize     (1,1)   double {mustBePositive}      = 50
                options.nonMaxSupDist   (1,1)   double {mustBePositive}      = 5
                options.buffer          (1,1)   double {mustBeNonnegative}   = 0
            end
            
            % get parameters             
            minStarSize = options.minStarSize; 
            maxStarSize = options.maxStarSize; 
            nonMaxSupDist = options.nonMaxSupDist;
            buffer = options.buffer;

            % get size of image
            n = imgObj.rows;
            m = imgObj.cols;

            % get image matrix
            img = imgObj.DNmat;
    
            % mask where image intensity is above the threshold to a binary image
            % DEPENDENCY: IMAGE PROCESSING TOOLBOX
            BW = imbinarize(img, thresh);
    
            % find connected components
            % DEPENDENCY: IMAGE PROCESSING TOOLBOX
            CC = bwconncomp(BW);
    
            % preallocate space for centers
            centers = nan(2,length(CC.PixelIdxList));
            
            % preallocate space for saving CCs
            comps = cell(2,length(CC.PixelIdxList));

            % only save connected components within the size range and not
            % beyond the image edge boundary
            for i = 1:length(CC.PixelIdxList)
                group = CC.PixelIdxList{i};
                [row, col] = ind2sub([n,m], group);
                dn = img(group); 
        
                % if this group is too small, or in the buffer, skip
                if length(group) < minStarSize || length(group) > maxStarSize || any(row<buffer) || ...
                        any(row>n-buffer) || any(col<buffer) || any(col>m-buffer)
                    continue;
                end

                % moment method with the flattened image
                xc = sum(col.*dn)/sum(dn);
                yc = sum(row.*dn)/sum(dn); 

                % distance between this point and all others for non max
                % suppression
                dists = sqrt((xc - centers(1,:)).^2 + (yc - centers(2,:)).^2);

                % create mask where distance is less than search radius
                distsMask = dists < nonMaxSupDist;

                % if there are centers too close, check DNs and only take
                % brightest
                if any(distsMask)
                    xcClose = centers(1,distsMask);
                    ycClose = centers(2,distsMask);
    
                    % get brightness from flat image of both centers
                    cDN = img(round(yc), round(xc));
                    cCloseDN = img(round(ycClose), round(xcClose));
    
                    % if any of the previous pixels are brighter, skip this
                    % pixel
                    if all(cCloseDN >= cDN)
                        continue;
                    % otherwise remove all of the previous ones and then 
                    % proceed with adding the current
                    else
                        centers(:,distsMask) = nan(2,sum(distsMask));
                    end
                end

                centers(:,i) =[xc; yc];
                comps{i} = group;

            end
    
            % remove any nans from centers
            nanMask = isnan(centers(1, :));
            centers(:, nanMask) = [];
            
            % remove all empty cells from the connected components
            comps = comps(~cellfun('isempty',comps));

            if isempty(centers)
                % if no centroids were found, output an empty points object
                centroids = sonic.Points2.empty;
                comps = {};
            else
                % otherwise, save centers as P2 points
                centroids = sonic.Points2(centers);
            end
        end

    end
end