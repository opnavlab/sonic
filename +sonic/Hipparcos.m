% This software is made available under the MIT license. See SONIC/LICENSE
% for details.
classdef Hipparcos < sonic.StarCatalog

    properties (SetAccess = {?sonic.StarCatalog})
        id                  (:, 1) uint64
        Vmag                (:, 1) double
        bayer               (:, 1) cell
        const_bayer         (:, 1) cell
        const_radec         (:, 1) cell   
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
                DB_LOC = '+sonic/+data/hipparcos.mat';
                hip_db = load(DB_LOC);
                hip_db = hip_db.hip_db;
            else
                hip_db = varargin{1};
            end

            % Split out the raw database into the constituent parts:
            obj.n = size(hip_db, 1);
            obj.id = hip_db.HIP;
            obj.Vmag = hip_db.Vmag;
            obj.t_epoch = -2.761289418143435e+08; % Hipparcos epoch J1991.25 in SPICE TDB
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