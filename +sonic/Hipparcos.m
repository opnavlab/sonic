% This software is made available under the MIT license. See SONIC/LICENSE
% for details.
classdef Hipparcos < sonic.StarCatalog

    properties
        n                   (1, 1) uint64
        id                  (:, 1) uint64
        Vmag                (:, 1) double
        ra_RAD              (:, 1) double
        dec_RAD             (:, 1) double
        pm_ra_RPS           (:, 1) double
        pm_dec_RPS          (:, 1) double
        plx_RAD             (:, 1) double
        bayer               (:, 1) cell
        const_bayer         (:, 1) cell
        const_radec         (:, 1) cell        
    end

    properties (SetAccess=private)
        filter_map          (:, 1) uint64
    end

    methods
    
        function obj = Hipparcos(varargin)
        %% obj = Hipparcos(varargin)
        %   Instantiate the Hipparcos catalog from a SONIC-generated
        %   Hipparcos datafile.
        %   
        %   Inputs:
        %       - [No Input]: uses the internal cache location to find the 
        %         SONIC-generated Hipparcos datafile
        %         OR
        %       - hip_db (table): table from the SONIC-generated Hipparcos
        %         datafile
        %   Outputs:
        %       - obj (sonic.Hipparcos): Hipparcos catalog object
        %
        %   Last revised: 3/29/24
        %   Last author: Michael Krause
        
            % Add in a nice error handler here
            if nargin < 1
                % Load in the parsed data file:
                DB_LOC = '+sonic/+data/hipdata.mat';
                hip_db = load(DB_LOC);
                hip_db = hip_db.hip_db;
            else
                hip_db = varargin{1};
            end

            % Split out the raw database into the constituent parts:
            obj.n = size(hip_db, 1);
            obj.id = hip_db.HIP;
            obj.Vmag = hip_db.Vmag;
            obj.ra_RAD = sonic.Units.DEGtoRAD(hip_db.RA_DEG);
            obj.dec_RAD = sonic.Units.DEGtoRAD(hip_db.DEC_DEG);
            obj.pm_ra_RPS = sonic.Units.MASPYtoRPS(hip_db.pmRA_MASPY);
            obj.pm_dec_RPS = sonic.Units.MASPYtoRPS(hip_db.pmDEC_MASPY);
            obj.plx_RAD = sonic.Units.MAStoRAD(hip_db.Plx_MAS);
            obj.bayer = hip_db.Bayer;
            obj.const_bayer = hip_db.Const;
            obj.const_radec = hip_db.ConstRaDec;

            obj.filter_map = [];
        end

        function filtered_cat = filter(obj, filt_vec)
        %% filtered_cat = filter(obj, filt_vec)
        %   Filter the catalog based on the supplied binary vector, which
        %   can be arbitrarily crafted.
        %   
        %   Inputs:
        %       - obj (sonic.Hipparcos): Hipparcos catalog object
        %       - filt_vec (1xm double): Vector of indices to filter to 
        %         (where `m` will be the size of the resultant catalog), OR
        %         a logical vector of size `n`, with `m` true values, where
        %         only those marked `true` will be retained by the filter. 
        %   Outputs:
        %       - filtered_cat (sonic.Hipparcos): Filtered Hipparcos 
        %         catalog object, with `n` corresponding to the number of 
        %         `true` found in bin_vec.
        %
        %   Last revised: 2/14/24
        %   Last author: Michael Krause
        
        arguments
            obj         (1, 1)      sonic.Hipparcos
            filt_vec    (1, :)      double
        end

            bin_vec = logical(filt_vec);
        
            if all(filt_vec == bin_vec) && (length(filt_vec) ~= obj.n)
                error('sonic:Hipparcos:filter:invalidFilterLength', ...
                    ['Invalid filter length (=%d, must =%d). Filter ' ...
                    'must be a 1-by-n logical vector, where n is the ' ...
                    'number of elements in the catalog being ' ...
                    'filtered, OR a 1-by-m vector of indices to ' ...
                    'retain while filtering.'], ...
                    length(filt_vec), obj.n);
            end

            % The rationale here is that users can concoct their own filter
            % with standard Matlab syntax, and that'll just form a binary
            % vector or vector of indices. 

            % This will then copy the catalog to a new object, grab each
            % property, and filter accordingly. A couple of properties are
            % special cases.
            filtered_cat = obj;
            cat_props = properties(filtered_cat);

            for idx = 1:length(cat_props)
                switch cat_props{idx}
                    case 'n'
                        % Save the new number of stars directly:
                        filtered_cat.n = sum(bin_vec);
                    case 'filter_map'
                        % Save the filter map in case it's needed for
                        % traceability
                        if ~all(filt_vec == bin_vec)
                            filtered_cat.filter_map = filt_vec;
                        else
                            filtered_cat.filter_map = find(filt_vec);
                        end
                    case 'const_data'
                        % Do nothing here, just preserve the constellation
                        % data. 
                    otherwise
                        % Filter either with binary vec or the vector of
                        % indices, whichever applies. 
                        if all(filt_vec == bin_vec)
                            filtered_cat.(cat_props{idx}) = ...
                                filtered_cat.(cat_props{idx})(bin_vec);
                        else
                            filtered_cat.(cat_props{idx}) = ...
                                filtered_cat.(cat_props{idx})(filt_vec);
                        end
                        
                end
            end
        
        end

        function u = eval(obj, et, r_obs_AU)
        %% u = eval(obj, et, r_obs_AU)
        %   Evaluate the Hipparcos catalog per the 5-parameter standard
        %   model.
        %   
        %   Inputs:
        %       - obj (sonic.Hipparcos): Hipparcos catalog object
        %       - et (1x1 double): Epoch at which to evaluate the catalog,
        %         in TDB (consistent with SPICE ephemeris time)
        %       - r_obs_AU (3x1 double): The position of the observer in 
        %         ICRF, in astronomical units (AU).
        %   Outputs:
        %       - u (1x1 sonic.Points3): The resultant star unit vectors,
        %         packed into a Points3 object. 
        %
        %   Last revised: 2/14/24
        %   Last author: Michael Krause

        arguments
            obj         (1, 1)      sonic.Hipparcos
            et          (1, 1)      double
            r_obs_AU    (3, 1)      double
        end

            % Evaluate the catalog here according to the 5-parameter
            % Hipparcos standard model. 
            t_epoch = -2.761289418143435e+08; % Hipparcos epoch in SPICE TDB
            dt = et - t_epoch;
            
            % Define the LVLH frame on the celestial sphere where proper motion
            % occurs. These equations are carried out in a vectorized form, and the
            % horizontal concatenion + transpose is a function of the input vectors
            % being columns. 
            p = [-sin(obj.ra_RAD), cos(obj.ra_RAD), zeros(obj.n, 1)]';
            q = [-sin(obj.dec_RAD).*cos(obj.ra_RAD), -sin(obj.dec_RAD).*sin(obj.ra_RAD), cos(obj.dec_RAD)]';
            e0 = sonic.SphereCoords.sphereToCart(obj.ra_RAD, obj.dec_RAD);   % Obtain our initial e0 from RA and DEC.
        
            % Basic standard model equation. Note r is in AU, so no need for any
            % conversion term there. Normalization ensures these are unit vecs.
            u = e0 + dt.*(obj.pm_ra_RPS'.*p + obj.pm_dec_RPS'.*q) - (obj.plx_RAD'.*r_obs_AU);
            u = sonic.PointsS2(u./vecnorm(u));
        
        end
    
    end

    % The following methods are used to generate the Hipparcos datafile,
    % which is loaded upon instantiation of a Hipparcos object. This should
    % only need to be run once per installation. 
    methods (Static)
    
        function parseData(hip_db_path, bayer_db_path, save_loc)
        %% parseData(hip_db_path, bayer_db_path, save_loc)
        %   Load the raw Hipparcos catalog and the associated Bayer lookup
        %   tables, parse and combine them, and save it to the specified
        %   location.
        %   
        %   Inputs:
        %       - hip_db_path (char): Path to Hipparcos catalog
        %       - bayer_db_path (char): Path to Bayer catalog
        %       - save_loc (char): Path to save processed datafile
        %
        %   Last revised: 2/14/24
        %   Last author: Michael Krause

            % Process these two data sources into a single unified .mat
            % file. 

            hip_data_raw = readtable(hip_db_path);
            hip_parse_targs = sonic.Hipparcos.getHipparcosParseTargs();
            hip_data = sonic.Hipparcos.parseTable(hip_data_raw, hip_parse_targs);

            bayer_data_raw = readtable(bayer_db_path);
            bayer_parse_targs = sonic.Hipparcos.getBayerParseTargs();
            bayer_data = sonic.Hipparcos.parseTable(bayer_data_raw, bayer_parse_targs);

            hip_db = sonic.Hipparcos.matchHipBayer(hip_data, bayer_data);
            save(save_loc, 'hip_db');

        end
    
        function parsed_tab = parseTable(raw_tab, parse_targs)
        %% parsed_tab = parseTable(raw_tab, parse_targs)
        %   Given a table, extracts only the columns of interest and
        %   renames them as specified.
        %   
        %   Inputs:
        %       - raw_tab (table): an arbitrary table
        %       - parse_targs (1xn struct): a struct array, where each
        %         entry indicates a column to select and a desired name for
        %         each column
        %   Outputs:
        %       - parsed_tab (table): retained columns from the input
        %         table, renamed accordingly
        %
        %   Last revised: 2/14/24
        %   Last author: Michael Krause
        
            % Pull out the info we'll need to trim down the raw tables:
            targ_cols = [parse_targs.col];
            targ_names = {parse_targs.varname};

            % Trim the tables and rename vars accordingly:
            parsed_tab = raw_tab(:, targ_cols);
            parsed_tab.Properties.VariableNames = targ_names;

            % Scan for rows with NaN HIP IDs and remove them:
            bad_rows = isnan(table2array(parsed_tab(:, 1)));

            parsed_tab = parsed_tab(~bad_rows, :);
        
        end

        function comb_data = matchHipBayer(hip_data, bayer_data)
        %% comb_data = matchHipBayer(hip_data, bayer_data)
        %   Given parsed tables of Hipparcos and Bayer catalog data, 
        %   matches the two and concatenates the tables appropriately.
        %   
        %   Inputs:
        %       - hip_data (table): table containing (parsed) Hipparcos 
        %         catalog
        %       - bayer_data (table): table containing (parsed) Bayer 
        %         catalog
        %   Outputs:
        %       - comb_data (table): combined table, with Hipparcos and
        %         Bayer data correlated
        %
        %   Last revised: 2/14/24
        %   Last author: Michael Krause

            [~, hip_matches] = ismember(hip_data.HIP, bayer_data.HIP);
            hip_match_idxs = find(hip_matches);

            % Need to create an "auxiliary" (sparse) Bayer designation 
            % table that we can just append to the Hipparcos table. This 
            % concatenates the matched data.

            % MATLAB is annoying about char arrays in tables. We know this,
            % what we're doing isn't breaking any rules, so temporarily
            % shut off warnings for this.
            warning('off', 'MATLAB:table:PreallocateCharWarning');
                
            % Preallocate the auxiliary table:
            bayer_aux = table('Size', [height(hip_data), width(bayer_data)-1], ...
                'VariableTypes', varfun(@(x) class(x{1}), bayer_data(1, 2:end), 'OutputFormat', 'cell'), ...
                'VariableNames', bayer_data.Properties.VariableNames(2:end));

            % Warnings back on.
            warning('on', 'MATLAB:table:PreallocateCharWarning');

            % Fill in the auxiliary table
            for idx = 1:length(hip_match_idxs)
                hip_idx = hip_match_idxs(idx);
                bayer_idx = hip_matches(hip_match_idxs(idx));

                bayer_aux(hip_idx, :) = bayer_data(bayer_idx, 2:end);
            end

            % Concatenate the resultant tables:
            comb_data = [hip_data, bayer_aux];

        end

        %% Parsing Parameterization for Hipparcos/Bayer Catalogs:
        function hip_parse_targs = getHipparcosParseTargs()
        %% hip_parse_targs = getHipparcosParseTargs()
        %   Yields the struct for appropriately parsing the raw Hipparcos
        %   table, used by sonic.Hipparcos.parseTable()
        %   
        %   Outputs:
        %       - hip_parse_targs (nx1 struct): a struct array, where each
        %         entry indicates a column to select and a desired name for
        %         each column
        %
        %   Last revised: 2/14/24
        %   Last author: Michael Krause

            % Hip ID
            hip_parse_targs(1).col =        2;
            hip_parse_targs(1).varname =    'HIP';            

            % Vmag
            hip_parse_targs(2).col =        6;
            hip_parse_targs(2).varname =    'Vmag';

            % Right Ascension (deg)
            hip_parse_targs(3).col =        9;
            hip_parse_targs(3).varname =    'RA_DEG';

            % Declination (deg)
            hip_parse_targs(4).col =        10;
            hip_parse_targs(4).varname =    'DEC_DEG';

            % Parallax
            hip_parse_targs(5).col =        12;
            hip_parse_targs(5).varname =    'Plx_MAS';

            % Proper Motion, Right Ascension
            hip_parse_targs(6).col =        13;
            hip_parse_targs(6).varname =    'pmRA_MASPY';

            % Proper Motion, Declination
            hip_parse_targs(7).col =        14;
            hip_parse_targs(7).varname =    'pmDEC_MASPY';

        end

        function bayer_parse_targs = getBayerParseTargs()
        %% bayer_parse_targs = getBayerParseTargs()
        %   Yields the struct for appropriately parsing the raw Bayer
        %   table, used by sonic.Hipparcos.parseTable()
        %   
        %   Outputs:
        %       - bayer_parse_targs (nx1 struct): a struct array, where 
        %         each entry indicates a column to select and a desired name 
        %         for each column
        %
        %   Last revised: 2/14/24
        %   Last author: Michael Krause

            % Hip ID
            bayer_parse_targs(1).col =      5;
            bayer_parse_targs(1).varname =  'HIP';

            % Bayer Des
            bayer_parse_targs(2).col =      9;
            bayer_parse_targs(2).varname =  'Bayer';

            % Constellation
            bayer_parse_targs(3).col =      10;
            bayer_parse_targs(3).varname =  'Const';

        end

        
    end

end