classdef USNOGNC < sonic.StarCatalog

    properties (SetAccess = {?sonic.StarCatalog})
        gnc_id              (:, 1) uint64
        hip_id              (:, 1) uint64
        gaia_id             (:, 1) uint64
        e_ra_MAS            (:, 1) uint64
        e_de_MAS            (:, 1) uint64
        e_pm_ra_RPS         (:, 1) double
        e_pm_dec_RPS        (:, 1) double
        e_plx_RAD           (:, 1) double
        cat_2016            (:, 1) string
        gmag                (:, 1) double
        e_gmag              (:, 1) double    
        bpmag               (:, 1) double
        e_bpmag             (:, 1) double
        rpmag               (:, 1) double
        e_rpmag             (:, 1) double
        jmag                (:, 1) double
        e_jmag              (:, 1) double
        hmag                (:, 1) double
        e_hmag              (:, 1) double
        kmag                (:, 1) double
        e_kmag              (:, 1) double  
        neighbor_id         (:, 1) double
        neighbor_d_AS       (:, 1) double
    end

    methods

        function obj = USNOGNC(varargin)
        %% obj = USNOGNC(varargin)
        %   Instantiate the U.S. Naval Observatory GNC catalog from the 
        %   internally stored GNC catalog datafile.
        %   
        %   Inputs:
        %       - [No Input]: uses the internal cache location to find the 
        %         SONIC-generated USNO GNC datafile
        %         OR
        %       - gnc_db (table): table from the SONIC-generated USNO GNC
        %         datafile
        %   Outputs:
        %       - obj (sonic.USNOGNC): USNOGNC catalog object
        %
        %   Last revised: 10/30/24
        %   Last author: Ava Thrasher
        
            % Add in a nice error handler here
            if nargin < 1
                % Load in the parsed data file:
                DB_LOC = '+sonic/+data/usnognc.mat';
                gnc_db = load(DB_LOC);
                gnc_db = gnc_db.usnognc;
            else
                gnc_db = varargin{1};
            end

            % Split out the raw database into the constituent parts:
            obj.n               = size(gnc_db, 1);
            obj.gnc_id          = gnc_db.gnc_id;
            obj.hip_id          = gnc_db.hip_id;
            obj.gaia_id         = gnc_db.gaiadr3_id;
            obj.t_epoch         = 5.049216681839228e+08; % GNC epoch J2016 in SPICE TDB

            % Star right ascension and declinations stored in radians
            obj.ra_RAD          = sonic.Units.DEGtoRAD(gnc_db.ra_2016);
            obj.dec_RAD         = sonic.Units.DEGtoRAD(gnc_db.de_2016);
            obj.e_ra_MAS        = gnc_db.e_ra_2016;
            obj.e_de_MAS        = gnc_db.e_de_2016;

            % proper motion
            obj.pm_ra_RPS       = sonic.Units.MASPYtoRPS(gnc_db.pmra);
            obj.pm_dec_RPS      = sonic.Units.MASPYtoRPS(gnc_db.pmde);
            obj.e_pm_ra_RPS     = gnc_db.e_pmra;
            obj.e_pm_dec_RPS    = gnc_db.e_pmde;

            % parallax
            obj.plx_RAD         = sonic.Units.MAStoRAD(gnc_db.plx);
            obj.e_plx_RAD       = sonic.Units.MAStoRAD(gnc_db.e_plx);

            % catalog designation
            obj.cat_2016        = gnc_db.cat_2016;

            % Magnitudes by different pass band filters
            obj.gmag            = gnc_db.gmag;
            obj.e_gmag          = gnc_db.e_gmag;
            obj.bpmag           = gnc_db.bpmag;
            obj.e_bpmag         = gnc_db.e_bpmag;
            obj.rpmag           = gnc_db.rpmag;
            obj.e_rpmag         = gnc_db.e_rpmag;
            obj.jmag            = gnc_db.jmag;
            obj.e_jmag          = gnc_db.e_jmag;
            obj.hmag            = gnc_db.hmag;
            obj.e_hmag          = gnc_db.e_hmag;
            obj.kmag            = gnc_db.kmag;
            obj.e_kmag          = gnc_db.e_kmag;

            % Nearest neighbor information
            obj.neighbor_id     = gnc_db.neighbor_id;
            obj.neighbor_d_AS   = gnc_db.neighbor_distance_arcsec;

            obj.filter_map = [];
        end

    end

end