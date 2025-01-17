% This software is made available under the MIT license. See SONIC/LICENSE
% for details.
classdef Kvector
    % Last Revised: 5/1/2024
    % Last author: Sebastien Henry

    % Kvector object, using the representation in References [1] and [2]

    % [1] Mortari, Daniele, et al. "The pyramid star identification technique."
    % Navigation 51.3 (2004): 171-183.
    %  https://doi.org/10.1002/j.2161-4296.2004.tb00349.x

    % [2] Christian, John A., and Joy Arulraj. "Review of the k-Vector and 
    % its Relation to Classical Data Structures." 
    % Journal of Guidance, Control, and Dynamics 45.11 (2022): 2121-2127.
    % https://doi.org/10.2514/1.G006688

    % I, J are catalog entries index of pair with interstar angle
    % note: do not put size on interstar_angle, Is, and Js, as this will
    % significantly slow down the init function!

    % no size restriction on cos_interstar_angle, Is, and Js because it
    % will drastically slow down their creation after
    properties
        n                   (1, 1)  int64 % number of pairs
        cos_interstar_angle         double % 1xn
        Is                          int64  % 1xn, always smaller than Js
        Js                          int64  % 1xn
        k                   (1, :)  double % the "k" vector
        bin_width           (1, 1)  double % constant bin width in cos(angle)
        n_bins              (1, 1)  int64  % number of bins in histogram
        ymin                (1, 1)  double % min cos interstar angle
        ymax                (1, 1)  double % max cos interstar angle
        m                   (1, 1)  double % slope
        q                   (1, 1)  double % intercept
        et                  (1, 1)  double 
        r_obs_AU            (3, 1)  double
    end
    
    methods
        function obj = Kvector(cat, min_angle, max_angle, max_mag, ...
                bin_width_cos_angle, et, r_obs_AU)
            %% obj = Kvector(varargin)
            %   Instantiate the Kvector table
            %   no preallocation but fast append in MATLAB
            %
            %   Inputs:
            %       - hip_cat (sonic.StarCatalog): Star Catalog, or
            %         subset thereof that we wish to create a k-ver
            %       - min_angle (1x1 double): the minimum interstar angle
            %         that will be stored in the k-vector (radians)
            %       - max_angle (1x1 double): the maximum interstar angle
            %         that will be stored in the k-vector (radians)
            %       - max_Vmag (1x1 double): maximum visual magnitude
            %       - bin_width_cos_angle (1x1): OPTIONAL
            %         the width of a bin of the histogram reprentation of
            %         the k-vector. Default to 0, which means there will be
            %         as many bins as interstar angle pairs.
            %       - et (1x1 double): OPTIONAL Epoch at which to evaluate
            %         the catalog, in TDB (consistent with SPICE ephemeris
            %         time)
            %       - r_obs_AU (3x1 double): OPTIONAL The position of the
            %         observer in ICRF, in astronomical units (AU).
            %   Outputs:
            %       - obj (sonic.Kvector): Kvector object
            %
            %   Last revised: 4/7/24
            %   Last author: Sebastien Henry
            arguments
            cat                 (1, 1)      sonic.StarCatalog
            min_angle           (1, 1)      double
            max_angle           (1, 1)      double
            max_mag            (1, 1)      double
            bin_width_cos_angle (1, 1)      double = 0
            et                  (1, 1)      double = 0 % default at time 0
            r_obs_AU            (3, 1)      double = [0;0;0] % default at SSB
            end

            % Check that the full star catalog has been provided:
            if ~isempty(cat.filter_map)
                error('sonic:StarId:nonBaseCatalogEntered', ...
                    'A full, unfiltered Star Catalog must be provided.')
            end

            % filter stars that are too faint
            if isa(cat,"sonic.Hipparcos")
                cat = cat.filter(cat.Vmag < max_mag);
            elseif isa(cat,"sonic.USNOGNC")
                cat = cat.filter(cat.gmag < max_mag);
            end

            % register the time and central position
            obj.et = et;
            obj.r_obs_AU = r_obs_AU; 

            eval = cat.eval(obj.et, obj.r_obs_AU);
            
            for i = 1:eval.n-1
                % first star direction
                ui = eval.u(:,i);
                
                % compute all possible interstar angle at once
                cos_ij_remaining = ...
                    ui'*eval.u(:,(i+1):end);
                
                % only keep the interstar angles between max_angle and
                % min_angle
                locs = cos_ij_remaining >= cos(max_angle) ...
                        & cos_ij_remaining <= cos(min_angle);
                cos_ij_valid = cos_ij_remaining(locs);
                
                % transform from boolean to indices
                locs = find(locs);
                
                % record interstar angles and catalog indices
                % these vectors are not preallocated because
                % 1) we don't know how large they will grow and they may
                % grow very large
                % 2) append elements in this direction is pretty fast
                % compared to other methods to append
                obj.cos_interstar_angle(end+1:end+length(cos_ij_valid)) = ...
                    cos_ij_valid;
                obj.Is(end+1:end+length(locs)) = ... 
                    cat.filter_map(i*uint64(ones(1,length(locs))));
                obj.Js(end+1:end+length(locs)) = ...
                    cat.filter_map(i+uint64(locs));

            end

            % sort by interstar angle
            [obj.cos_interstar_angle, indices] = sort(obj.cos_interstar_angle);
            obj.Is = obj.Is(indices);
            obj.Js = obj.Js(indices);
            
            % number of interstar angle pairs
            obj.n = length(obj.cos_interstar_angle);

            % fetch min and max cos of interstar angle possible in kvector
            obj.ymin = obj.cos_interstar_angle(1);
            obj.ymax = obj.cos_interstar_angle(end);
            span = obj.ymax - obj.ymin;
            
            % compute number of bins given the user specified parameter
            % too many bins, defaults to obj.n
            if bin_width_cos_angle > 0
                n_bins = ceil(span / bin_width_cos_angle);
                obj.n_bins = min(n_bins, obj.n);
            else
                obj.n_bins = obj.n;
            end

            obj.bin_width = span / double(obj.n_bins);

            % tolerance 
            Xsi = sonic.Tolerances.SlopeAdjust * ...
                max(abs(obj.ymin), abs(obj.ymax));

            % compute slope of the kvector
            obj.m = (span + 2 * Xsi) / (double(obj.n_bins));
            
            % compute the intercept of the kvector
            obj.q = obj.ymin - obj.m - Xsi;
            
            % count the number of elements per bin
            binedges = linspace(obj.ymin, obj.ymax, obj.n_bins+1);
            [N, ~] = histcounts(obj.cos_interstar_angle, binedges);

            % create k-vector
            obj.k = [0, cumsum(N)];
        end
        
        function [min_idx, max_idx] = query(obj, min_angle, max_angle)
            %% [min_idx, max_idx] = query(obj, min_angle, max_angle)
            %  Query the Kvector table and return the indices that
            %  define the limits for the interstar angles
            %
            %   Inputs:
            %       - obj (1x1 sonic.Kvector): Kvector object
            %       - min_angle (1x1 double): the minimum interstar angle 
            %         (radians)
            %       - max_angle (1x1 double): the maximum interstar angle 
            %         (radians)
            %   Outputs:
            %       - min_idx (1x1, double): the approximate index
            %         corresponding to min_angle
            %       - max_idx (1x1, double): the approximate index
            %         corresponding to max_angle
            %
            %   Last revised: 4/12/24
            %   Last author: Sebastien Henry

            arguments
                obj         (1, 1)  sonic.Kvector
                min_angle   (1, 1)  double
                max_angle   (1, 1)  double
            end

            % input check
            if min_angle > max_angle
                error('sonic:KvectorError:invalidInput', ...
                    'min_angle must be smaller than max_angle');
            end
            
            % compute cosine of angles
            cos_min_angle = cos(min_angle);
            cos_max_angle = cos(max_angle);
            
            % min_angle > max_angle of kvec
            if cos_min_angle < obj.cos_interstar_angle(1)
                min_idx = [];
                max_idx = [];
                return
            end
            
            % max_angle < min_angle of kvec
            if cos_max_angle > obj.cos_interstar_angle(end)
                min_idx = [];
                max_idx = [];
                return
            end

            % approximate indices that bracket max angle
            j_bottom = floor((cos_max_angle - obj.q)/ obj.m);
            j_bottom = max(j_bottom, 1);

            min_idx = obj.k(j_bottom);
            if min_idx == 0
                min_idx = 1;
            end

            % approximate indices that bracket min angle
            j_top = ceil((cos_min_angle-obj.q)/obj.m);
            j_top = min(j_top, length(obj.k));
            max_idx = obj.k(j_top);

        end
    end
end



