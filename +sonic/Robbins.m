% This software is made available under the MIT license. See SONIC/LICENSE
% for details.
classdef Robbins
    % This class instantiates a Robbins database object, which contains 
    % information about ~1.3 million lunar crater positions and sizes, as
    % documented in the Robbins crater database. The database can be found
    % at this link: https://astrogeology.usgs.gov/search/map/moon_crater_database_v1_robbins
    % The current version of the database in sonic was obtained from the
    % link above, as downloaded on Oct 21, 2024
    %
    % More information about the Robbins crater database, including details
    % about how it was generated can be found at this reference: 
    %    Robbins, S. J. (2019). A new global database of lunar impact 
    %    craters >1-2 km: 1. Crater locations and sizes, comparisons with 
    %    published databases, and global analysis. Journal of Geophysical 
    %    Research: Planets, 124(4), pp. 871-892. 
    %    https://doi.org/10.1029/2018JE005592
    %
    % TODO: Add property descriptions
    properties
        n                       (1, 1) uint64     
        id                      (:, 1) string
        latCir_RAD              (:, 1) double
        lonCir_RAD              (:, 1) double
        latElp_RAD              (:, 1) double
        lonElp_RAD              (:, 1) double
        latElpDev_RAD           (:, 1) double
        lonElpDev_RAD           (:, 1) double
        diamCir_KM              (:, 1) double
        diamCirDev_KM           (:, 1) double
        majAxElp_KM             (:, 1) double
        minAxElp_KM             (:, 1) double
        angElp_RAD              (:, 1) double
        majAxElpDev_KM          (:, 1) double
        minAxElpDev_KM          (:, 1) double
        angElpDev_RAD           (:, 1) double
        eccElp                  (:, 1) double
        elptElp                 (:, 1) double
        eccElpDev               (:, 1) double
        elptElpDev              (:, 1) double
        arcRim                  (:, 1) double
        numPts                  (:, 1) uint16   % uint16 (0-65535) max is 8088
        centerPtM_KM            (3, :) double
        surfNorm                (3, :) double
        diskQuadEnv_KM          (4, 4, :) double 
        conicsENU_KM            (5, :) double
        attQuatENUtoM_vs        (4, :) double
        dqUsesCircInfo          (:, 1) logical
    end

    properties (SetAccess=private)
        filterMap           (:, 1) uint64
    end

    methods
    
        function obj = Robbins(varargin)
        %% obj = Robbins(varargin)
        %   Instantiate the Robbins crater database from a SONIC-generated
        %   Robbins datafile.
        %   
        %   Inputs:
        %       - [No Input]: uses the internal cache location to find the 
        %         sonic-generated Robbins datafile
        %         OR
        %       - robTab (table): table from the sonic-generated Robbins
        %         datafile
        %   Outputs:
        %       - obj (sonic.Robbins): Robbins crater database object
        %
        %   Last revised: Jan 11, 2024
        %   Last author: Tara Mina
        
            loadParsedRobbins=true;
            if nargin == 1
                rbn = varargin{1};
                obj.checkValidRobbinsDatafile(rbn);
                loadParsedRobbins=false;
            elseif nargin > 1
                error('sonic:Robbins:invalidInput', ...
                    ['Too many input arguments. Either construct a ', ...
                    'Robbins() object without inputs or provide the ', ...
                    'Robbins datafile as a table object.']);
            end

            if loadParsedRobbins
                % Load in the parsed data file
                dbLoc = '+sonic/+data/robbinsData.mat';
                rbnMat = load(dbLoc);
                rbn = rbnMat.robTab;

                % Load in center points and surface normals
                obj.centerPtM_KM = rbnMat.centerPtM_KM;
                obj.surfNorm = rbnMat.surfNorm;

                % save disk quadric, attitude, conics info
                obj.diskQuadEnv_KM = rbnMat.diskQuadEnv_KM;
                obj.dqUsesCircInfo = rbnMat.dqUsesCircInfo';
                obj.attQuatENUtoM_vs = rbnMat.attQuatENUtoM_vs;
                obj.conicsENU_KM = rbnMat.conicsExplENU_KM;
            end

            % Extract key terms from raw database
            
            % number of craters
            obj.n = size(rbn, 1);

            % crater IDs
            obj.id = rbn.CRATER_ID;

            % lat/lon of circle and ellipse fits + standard deviation
            obj.latCir_RAD = sonic.Units.DEGtoRAD(rbn.LAT_CIRC_IMG);
            obj.lonCir_RAD = sonic.Units.DEGtoRAD(rbn.LON_CIRC_IMG);
            obj.latElp_RAD = sonic.Units.DEGtoRAD(rbn.LAT_ELLI_IMG);
            obj.lonElp_RAD = sonic.Units.DEGtoRAD(rbn.LON_ELLI_IMG);
            obj.latElpDev_RAD = sonic.Units.DEGtoRAD(rbn.LAT_ELLI_SD_IMG);
            obj.lonElpDev_RAD = sonic.Units.DEGtoRAD(rbn.LON_ELLI_SD_IMG);

            % diameter of circle fit + standard deviation
            obj.diamCir_KM = rbn.DIAM_CIRC_IMG;
            obj.diamCirDev_KM = rbn.DIAM_CIRC_SD_IMG;

            % ellipse fit major/minor axes and clocking angle + std dev
            obj.majAxElp_KM = rbn.DIAM_ELLI_MAJOR_IMG;
            obj.minAxElp_KM = rbn.DIAM_ELLI_MINOR_IMG;
            obj.angElp_RAD = sonic.Units.DEGtoRAD(rbn.DIAM_ELLI_ANGLE_IMG);
            obj.majAxElpDev_KM = rbn.DIAM_ELLI_MAJOR_SD_IMG;
            obj.minAxElpDev_KM = rbn.DIAM_ELLI_MINOR_SD_IMG;
            obj.angElpDev_RAD = sonic.Units.DEGtoRAD(rbn.DIAM_ELLI_ANGLE_SD_IMG); 

            % eccentricity/ellipticity of ellipse fits + standard deviation
            obj.eccElp = rbn.DIAM_ELLI_ECCEN_IMG;
            obj.elptElp = rbn.DIAM_ELLI_ELLIP_IMG;
            obj.eccElpDev = rbn.DIAM_ELLI_ECCEN_SD_IMG; 
            obj.elptElpDev = rbn.DIAM_ELLI_ELLIP_SD_IMG; 
            
            % arc around crater rim, number of points sampled along rim
            obj.arcRim = rbn.ARC_IMG;
            obj.numPts = rbn.PTS_RIM_IMG;

            % define an empty filter map array
            obj.filterMap = [];
        end

        function filteredCat = filter(obj, filtVec)
            %% filtered_cat = filter(obj, filtVec)
            %   Filter the catalog based on the supplied binary vector, 
            %   which can be arbitrarily crafted.
            %   
            %   Inputs:
            %       - obj (sonic.Robbins): Robbins crater catalog object
            %       - filtVec (1D double, can be 1xm or mx1): Vector of  
            %         indices to filter to (`m` to be size of resultant 
            %         catalog), OR a logical vector of size `n`, with `m` 
            %         true values, where only those marked `true` will be 
            %         retained by the filter. 
            %   Outputs:
            %       - filteredCat (sonic.Robbins): Filtered Robbins 
            %         catalog object, with `n` corresponding to number of 
            %         `true` found in bin_vec.
            %
            %   Last revised: Jan 13, 2024
            %   Last author: Tara Mina
            
            arguments
                obj        (1, 1)      sonic.Robbins
                filtVec    (:, :)      double
            end

            assert( size(filtVec,1)==1 || size(filtVec,2)==1, ...
                'filtVec must be 1D (either a column or row vector)' );
            % make filtVec a column vector
            filtVec = reshape(filtVec, [length(filtVec),1]);
            binVec = logical(filtVec);
            isLogicalArr = all(filtVec == binVec);
        
            if isLogicalArr 
                if (length(filtVec) ~= obj.n)
                    error('sonic:Robbins:filter:invalidFilterLength', ...
                        ['Invalid filter length (=%d, must = %d). Filter ' ...
                        'must be a 1-by-n logical vector, where n is the ' ...
                        'number of elements in the catalog being ' ...
                        'filtered, OR a 1-by-m vector of indices to ' ...
                        'retain while filtering.'], ...
                        length(filtVec), obj.n);
                end
                filtVec = binVec;
            end

            % NOTE: added these checks to validate that if filtVec is an
            % array of indices, they are integers and in the bounds of the
            % catalog
            if ~isLogicalArr
                assert(all(round(filtVec)==filtVec), ...
                    'sonic:Robbins:filter:invalidFilterIndices', ...
                    ['Invalid values for non-logical filter vector. ', ...
                    'Filter vector must have all indices be integers ']);
                assert(min(filtVec)>0 && max(filtVec)<=obj.n, ...
                    'sonic:Robbins:filter:invalidFilterIndices', ...
                    ['Invalid values for non-logical filter vector. ', ...
                    'Filter vector must have all indices > 0 and ', ...
                    '<= to the number of elements in the Robbins ', ...
                    'catalog.']);
            end

            % The rationale here is that users can concoct their own filter
            % with standard Matlab syntax, and that'll just form a binary
            % vector or vector of indices. 

            % This will then copy the catalog to a new object, grab each
            % property, and filter accordingly. A couple of properties are
            % special cases.
            filteredCat = obj;
            catProps = properties(filteredCat);

            vectorVals = ["centerPtM_KM", "surfNorm", ...
                "attQuatENUtoM_vs", "conicsENU_KM"];

            for idx = 1:length(catProps)
                switch catProps{idx}
                    case 'n'
                        % Save the new number of craters directly:
                        filteredCat.n = sum(binVec);
                    case 'filterMap'
                        % Save the filter map in case it's needed for
                        % traceability
                        if isLogicalArr
                            filteredCat.filterMap = find(filtVec);
                        else
                            filteredCat.filterMap = filtVec;
                        end
                    otherwise
                        propStr = catProps{idx};
                        currVals = filteredCat.(propStr);
                        
                        if strcmp(propStr, 'diskQuadEnv_KM')
                            filteredCat.(propStr) = ...
                                currVals(:,:,filtVec);
                        elseif ismember(propStr, vectorVals)
                            filteredCat.(propStr) = ...
                                currVals(:,filtVec);
                        else
                            filteredCat.(propStr) = currVals(filtVec);
                        end
                end
            end
        end

        function result = subsref(obj, S)
            % Overload the subsref method to support parentheses ()
            switch S(1).type
                case '()'
                    % Call the filter method when parentheses are used
                    result = obj.filter(S(1).subs{:});
                otherwise
                    % Fallback to built-in subsref for other types
                    result = builtin('subsref', obj, S);
            end
        end
    end

    methods(Static, Access=public)
        function [dq, conicENU, attPtoM, centerPt_M, uVec] = ...
                getMCMFDiskQuadric(lat_rad, lon_rad, majAx_km, ...
                minAx_km, ang_rad)

            % Get conic envelope
            explConic = [0; 0; majAx_km/2; minAx_km/2; ang_rad];
            conicENU = sonic.Conic(explConic);
            conicEnv = conicENU.envelope;

            % get normal vector ("up") based on latitude longitude 
            % (assuming sphere model now -- could change to ellipsoid)
            uVec = sonic.SphereCoords.raDecToCart(lon_rad, lat_rad);

            % get direction of the major axis (p in perifocal frame)
            eVec = [-sin(lon_rad); cos(lon_rad); 0];
            nVec = [-sin(lat_rad)*cos(lon_rad); -sin(lat_rad)*sin(lon_rad); cos(lat_rad)];

            % Moon is approximately a sphere
            moonRadius_KM = sonic.Constants.MOON_radii_KM(1); 
            centerPt_M = moonRadius_KM*uVec;

            % transformation matrix from MCMF to local ENU frame, def as P
            T_PtoM = [eVec, nVec, uVec];
            attPtoM = sonic.Attitude(T_PtoM);

            % Get disk quadric envelope
            H_M = [eVec, nVec, centerPt_M];
            H_Mk = [H_M; 0,0,1];
            dqEnv_M2 = H_Mk * conicEnv * H_Mk';
            dq = sonic.Quadric(dqEnv_M2, 'envelope');
        end
    end

    methods(Static, Access=private)
        function checkValidRobbinsDatafile(rbn)
            arguments
                rbn         (:, 21)      table
            end
        end
        function throwInvalidInputDatafileError()
            error('sonic:Robbins:invalidInputDatafile', ...
                    ['Please provide a valid Robbins ', ...
                    'crater database file as a table object.']);
        end
    end
end
