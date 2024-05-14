classdef (Abstract) StarCatalog

    properties
        const_data      (:, 1)      struct
    end

    methods

        function obj = StarCatalog()
        %% obj = StarCatalog()
        %
        %   Constructor for StarCatalog object. Note that StarCatalog is an
        %   abstract class, so it is not possible to directly instantiate
        %   this object. Instead, this constructor is intended for implicit
        %   use by any subclasses, such that they'll all be able to have
        %   access to constellation data.
        %
        %   Inputs:
        %       - None.
        %
        %   Outputs:
        %       - obj (1x1 concrete subclass of sonic.StarCatalog): As this
        %         will be implicitly called by any subclass of StarCatalog,
        %         this will be an instance of that subclass type, containing
        %         constellation data.
        %
        %   Last revised: 3/29/24
        %   Last author: Michael Krause
            
            % Have the StarCatalog object take care of the constellation
            % data. This will let any subclass of it have this available. 
            DB_LOC = '+sonic/+data/constdata.mat';
            const_data = load(DB_LOC);
            obj.const_data = const_data.const_struct;

        end

        function const = raDecToConstellation(obj, ra_RAD, dec_RAD)
        %% const = raDecToConstellation(obj, ra_RAD, dec_RAD)
        %
        %   Given a particular right ascension and declination on the
        %   celestial sphere, this will return a char array indicating the
        %   constellation that this point lies within. If, for whatever
        %   reason, no match is found, an empty char array is returned.
        %
        %   Inputs:
        %       - obj (1x1 concrete subclass of sonic.StarCatalog): object
        %         containing constellation data.
        %       - ra_RAD (1x1 double): Right ascension to query, in
        %         radians. 
        %       - dec_RAD (1x1 double): Declination to query, in radians. 
        %
        %   Outputs:
        %       - const (1x1 char array): Name of resultant constellation,
        %         where the point specified by (ra_RAD, dec_RAD) lies. 
        %
        %   Last revised: 3/29/24
        %   Last author: Michael Krause

            % Get this into a SONIC Points3 object. Note that this point
            % isn't *really* in P3 (it's in S2), but this will make it
            % easier to work with when calculating a signed distance to the
            % planes of interest later on. 
            targ_pt = sonic.Points3(sonic.SphereCoords.sphereToCart(ra_RAD, dec_RAD));

            % Let's do a prefiltering step first to check how close we are
            % to a given constellation. This helps prioritize which
            % constellations we should check first. 
            cand_idxs = obj.findNearestConstCands(targ_pt);

            % Want to find the index of the polygon we're lying in (where
            % this corresponds to the index of our constellation data)
            in_poly = -1;

            % Now we call the actual checking algorithm, and break when we
            % find a match:
            for const_idx = 1:length(cand_idxs)
                isInPolygon = obj.checkInSphericalPolygon(...
                    targ_pt, cand_idxs(const_idx));

                if isInPolygon
                    in_poly = const_idx;
                    break;
                end
            end

            if in_poly > 0
                % If found, return the constellation name:
                const = obj.const_data(cand_idxs(in_poly)).name;      
            else
                % If, for some reason, we could not find a match, we return
                % an empty string. This is to keep commonality with the
                % return type.
                const = '';
            end 
                    
        end

        function isInPolygon = checkInSphericalPolygon(obj, targ_pt, const_idx)
        %% isInPolygon = checkInSphericalPolygon(obj, targ_pt, const_idx)
        %
        %   Given a particular point and a spherical polygon (encoded as a
        %   series of unit vectors specifying vertices of the polygon),
        %   this checks if the point is in the spherical polygon. Note that
        %   this method implicitly assumes that each spherical polygon
        %   subtends less than 180 degrees on the celestial sphere. 
        %
        %   Inputs:
        %       - obj (1x1 concrete subclass of sonic.StarCatalog): object
        %         containing constellation data.
        %       - targ_pt (1x1 sonic.Points3): Unit vector specifying the
        %         query point.
        %       - const_idx (1x1 double): Index into the constellation data
        %         (contained within `obj`), yielding information on the
        %         spherical polygon we're attempting to test. 
        %
        %   Outputs:
        %       - isInPolygon (1x1 logical): True if point lies within
        %         spherical polygon, false otherwise. 
        %
        %   Last revised: 3/29/24
        %   Last author: Michael Krause

            % If the # of crossings is even, we're in the polygon, so we
            % need to keep track of these. 
            crossings = 0;

            % Convert the interior point and zero to be points in SONIC:
            int_pt = sonic.Points3(obj.const_data(const_idx).refstar);
            zero_pt = sonic.Points3([0;0;0]);

            % Short-circuit the generic join() logic since we know the
            % types here. It's a little faster. 
            int_pt_plane = sonic.GeometryP3.joinPointLine(targ_pt, ...
                    sonic.GeometryP3.joinPoints(int_pt, zero_pt));
            
            % Pull out the polygon boundaries (points and planes) from the
            % constellation data:
            poly_bounds = obj.const_data(const_idx).bounds;
            poly_bound_planes = obj.const_data(const_idx).bound_planes;

            % Iterate over each side of the polygon:
            for side_idx_1 = 1:length(poly_bounds)

                % Grab the (precomputed) plane corresponding to this side
                % (which passes through the origin of the sphere):
                plane = poly_bound_planes(side_idx_1);

                % We also need the points at the ends of the side.
                side_idx_2 = side_idx_1 + 1;
                if side_idx_2 > length(poly_bounds)
                    side_idx_2 = side_idx_2 - length(poly_bounds);
                end
                side_bnd_1 = poly_bounds(side_idx_1);
                side_bnd_2 = poly_bounds(side_idx_2);

                % Get signed distance for the known interior point and the
                % point of interest relative to the boundary plane. Given
                % that the points here don't strictly lie in P3 (they're in
                % S2!), we just do the dot product rather than a meet, to
                % reinforce this logic.
                dist_pt = plane.p3'*targ_pt.p3;
                dist_int = plane.p3'*int_pt.p3;

                % Get signed distance for each boundary point relative to
                % the plane between the known interior point and the point
                % of interest. Given that the points here don't strictly 
                % lie in P3 (they're in S2!), we just do the dot product 
                % rather than a meet, to reinforce this logic.
                dist_side_pt1 = int_pt_plane.p3'*side_bnd_1.p3;
                dist_side_pt2 = int_pt_plane.p3'*side_bnd_2.p3;

                % This checks if a given plane is crossed, and if so,
                % whether it was crossed within the bounds of the side:
                if (sign(dist_pt) ~= sign(dist_int)) && ...
                    (sign(dist_side_pt1) ~= sign(dist_side_pt2))
                
                    % Iterate the number of crossings
                    crossings = crossings + 1;
                
                end

            end

            % If we have an even number of crossings (including zero as
            % even), then we're in the spherical polygon!
            isInPolygon = (mod(crossings, 2) == 0);

        end

        function sort_mapping = findNearestConstCands(obj, targ_pt)
        %% sort_mapping = findNearestConstCands(obj, targ_pt)
        %
        %   Given a set of spherical polygons, each with a known point
        %   lying inside of the polygon, this method takes the dot product
        %   between our query point and each known point. These dot
        %   products are sorted, and polygons with the largest dot products
        %   (i.e., those whose known points lie "close" to the query
        %   points, in the angular sense), are checked first. This
        %   significantly accelerates the constellation matching algorithm.
        %
        %   Inputs:
        %       - obj (1x1 concrete subclass of sonic.StarCatalog): object
        %         containing constellation data.
        %       - targ_pt (1x1 sonic.Points3): Unit vector specifying the
        %         query point.
        %
        %   Outputs:
        %       - sort_mapping (1xm double): Given that `obj` contains
        %         information on `n` constellations, this returns the indices
        %         of the closest `m` constellations ordered such that
        %         the dot products corresponding to each index are in
        %         descending order. Note that m < n, where m is the number of
        %         constellations where the dot product with the query point
        %         is positive.
        %
        %   Last revised: 3/29/24
        %   Last author: Michael Krause   

            % Grab the dot-products between a point known to lie in each
            % polygon, and the point itself:
            const_pts = [obj.const_data.refstar];
            dists = const_pts'*targ_pt.r3;

            % Sort them (we'll check in this order), and only retain those
            % that are pointing in the same direction. This places a
            % natural limitation on the angular extent of a given spherical
            % polygon (less than 180 deg), but for all constellations, this
            % is fine. 
            [sorted_dists, sort_mapping] = sort(dists, 'descend');
            sort_mapping = sort_mapping(sorted_dists > 0);
            
        end

    end

    methods (Static)

        

        function parseConstellationFiles(const_dir)
        %% const_struct = parseConstellationFiles(const_dir)
        %   The IAU provides files delineating constellation boundaries.
        %   This function grabs those files, parses them, and outputs a
        %   struct with each constellation's boundaries.
        %
        %   The saved file contains a struct array containing 89
        %   entries, each with a name and 2-by-n list of RA/DEC
        %   vertices defining the boundary of that constellation. Note
        %   that while there are 88 constellations, Serpens is split
        %   into two distinct parts.
        %   
        %   Inputs:
        %       - const_dir (char): directory where IAU constellation files
        %         live
        %
        %   Last revised: 2/14/24
        %   Last author: Michael Krause

            % The IAU provides constellation boundaries in a series of text
            % files. This function consumes those files and processes them
            % into a single struct of constellation boundaries.

            const_files = dir(fullfile(const_dir, '*.txt'));

            % Loop backwards to preallocate the resultant boundary
            for idx = length(const_files):-1:1
                const_bounds = sonic.StarCatalog.parseConstFile(fullfile(...
                    const_files(idx).folder, const_files(idx).name));

                const_name = split(const_files(idx).name, '.');
                const_struct(idx).name = const_name{1};
                const_struct(idx).bounds = const_bounds;
            end

            % Save this as a .mat file rather than returning the struct
            save('+sonic/+data/constdata.mat');
            
        end

        function const_bounds = parseConstFile(const_file)
        %% const_bounds = parseConstFile(const_file)
        %   Parses an individual IAU constellation boundary file to yield
        %   the boundary list.
        %   
        %   Inputs:
        %       - const_file (char): constellation boundary filename
        %   Outputs:
        %       - const_bounds (2xn double): List of `n` RA/DEC (in
        %         degrees) boundary points
        %
        %   Last revised: 2/14/24
        %   Last author: Michael Krause
            

            fh = fopen(const_file, 'r');

            dec_DEG = [];
            ra_HMS = zeros(3, 0);

            while ~feof(fh)
                raw_line = fgetl(fh);
                comps = split(raw_line, '|');
                ra_str_HMS = strip(comps{1});
                dec_DEG(end+1) = str2double(strip(comps{2}));
                ra_split = split(ra_str_HMS, ' ');
                ra_HMS(:, end+1) = str2double(ra_split);
            end
            fclose(fh);
            ra_DEG = sonic.Units.HMStoDEG(ra_HMS);

            const_bounds = [ra_DEG; dec_DEG];
            
        end
        
    end

end