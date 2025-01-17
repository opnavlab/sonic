% This software is made available under the MIT license. See SONIC/LICENSE
% for details.
classdef StarId

    methods (Static)

        function [matches, att_ICRF2C] = interstarAngle( ...
                kvec, cat, pointsS2_meas, tol, max_angle, min_matches)
            %% [matches, att_ICRF2C] = interstarAngle(kvec, cat, pointsS2_meas, tol, max_angle, min_matches)
            %
            %   Takes in multiple line of sight measurements, in the camera
            %   frame, and matches them to stars in the given star catalog.
            %   Also returns an attitude estimate from ICRF to camera.
            %
            %   Inputs:
            %       - kvec (sonic.Kvector): a k-vector object
            %       - hip_cat (sonic.StarCatalog): star catalog
            %         NOTE: it MUST be an unfiltered catalog
            %       - points2_meas (sonic.Points2): unit-normalized vectors
            %         measured in the camera frame
            %       - tol (1x1 double): angular tolerance, a match is
            %         detected if the angle between measurement and catalog
            %         is less than tolerance
            %       - max_angle (1x1, double): the maximum angle between
            %         stars. Usually constrained by the camera FOV.
            %       - min_matches (1x1, uint64): the minimum number of
            %         pairs in matches to deem the solution correct. Minimum
            %         is 4. Default is 5.
            %   Outputs:
            %       - matches (2xn int): pairs of measurement-hipparcos index.
            %         In the i-th pair, element matches(1,i) is the number of
            %         the measurement and matches(2,i) is the corresponding
            %         index in the hipparcos catalog (!!! not the id!!!)
            %       - att_ICRF2C (sonic.Attitude): attitude of the camera
            %
            %   Last revised: 05/01/24
            %   Last author: Sebastien Henry
            arguments
                kvec            (1, 1)  sonic.Kvector
                cat             (1, 1)  sonic.StarCatalog
                pointsS2_meas   (1, 1)  sonic.PointsS2
                tol             (1, 1)  double
                max_angle       (1, 1)  double
                min_matches     (1, 1)  uint64 = 5
            end

            % check that user wants at least 4 matches
            if min_matches < 4
                error('sonic:StarId:minMatchesTooFew', ...
                    'Must require a minimum of 4 matches.');
            end

            % check that the full catalog has been provided:
            if ~isempty(cat.filter_map)
                error('sonic:StarId:nonBaseCatalogEntered', ...
                    'A full, unfiltered Star Catalog must be provided.')
            end

            n = pointsS2_meas.n;
            pointsS2_cat = cat.eval(kvec.et, kvec.r_obs_AU);

            % needs to be a fresh catalog!!!
            cat_idx = 1:cat.n;

            % use smart pattern of for loops
            % see EPS from
            % "Fast and robust kernel generators for star trackers"
            % https://doi.org/10.1016/j.actaastro.2017.02.016
            count = 0;
            for dj = 1:(n-2)
                for dk = 1:(n-dj-1)
                    for ii = 1:3
                        for i = ii:3:(n-dj-dk)
                            j = i + dj;
                            k = j + dk;

                            % for each triad, verify
                            ijk = [i,j,k];
                            [matches] = ...
                                sonic.StarId.checkTriad( ...
                                kvec, cat_idx, pointsS2_meas, pointsS2_cat, ...
                                ijk, tol, max_angle, min_matches);
                            count = count+1;

                            % if triad works out, then return
                            if size(matches, 2) >= min_matches
                                % recompute att_ICRF2C and return
                                u_cam   = pointsS2_meas.u(:,matches(1,:));
                                u_world = pointsS2_cat.u(:, matches(2,:));
                                att_ICRF2C = ...
                                    sonic.Attitude.solveWahbasProblem( ...
                                    sonic.PointsS2(u_cam), ...
                                    sonic.PointsS2(u_world));

                                return
                            end
                        end
                    end
                end
            end
        end

        function [matches, att_ICRF2C, valid] = checkTriad( ...
                kvec, cat_idx, pointsS2_meas, pointsS2_cat, ...
                ijk, tol, max_angle, min_matches)
            %% [matches, att_ICRF2C] = check_triad(kvec, cat_idx, pointsS2_meas, pointsS2_cat, ijk, tol, max_angle, min_matches)
            %
            %  Check if a triad is a valid triad, given a k-vector table.
            %
            %   Inputs:
            %       - kvec (1x1sonic.Kvector): a k-vector table
            %       - hip_cat_idx (:x1 uint64): indices of the catalog
            %       - points2_meas (sonic.Points2): unit-normalized vectors
            %         measured in the camera frame
            %       - pointsS2_cat (sonic.PointsS2): catalog line of sights
            %       - ijk (3x1 uint64): the indexes of the measurement
            %         triad we are checking
            %       - tol (1x1 double): angular tolerance
            %       - max_angle (1x1, double): the maximum angle between
            %         stars. Usually constrained by the camera FOV.
            %       - min_matches (1x1 uint64): minimum matches to stop
            %         and return
            %   Outputs:
            %       - matches (2x: int): pairs of measurement and id.
            %         In the i-th pair, element matches(1,i) is the number of
            %         the measurement and matches(2,i) is the corresponding
            %         catalog id
            %       - att_ICRF2C (sonic.Attitude): attitude of the camera
            %       - valid (1x1 boolean): whether sufficient matches were
            %       found
            %
            %   Last revised: 04/15/24
            %   Last author: Sebastien Henry

            arguments
                kvec            (1, 1)  sonic.Kvector
                cat_idx     (:, 1)  uint64
                pointsS2_meas   (1, 1)  sonic.PointsS2
                pointsS2_cat    (1, 1)  sonic.PointsS2
                ijk             (3, 1)  uint64
                tol             (1, 1)  double
                max_angle       (1, 1)  double
                min_matches     (1, 1)  uint64 = 5
            end
            valid = false;
            cosmaxangle = cos(max_angle);

            u_cat_all = pointsS2_cat.u;

            i = ijk(1);
            j = ijk(2);
            k = ijk(3);

            cosij = pointsS2_meas.u(:,i)'* pointsS2_meas.u(:,j);
            cosik = pointsS2_meas.u(:,i)'* pointsS2_meas.u(:,k);
            cosjk = pointsS2_meas.u(:,j)'* pointsS2_meas.u(:,k);

            angij = acos(cosij);
            angik = acos(cosik);
            angjk = acos(cosjk);

            angij_min = angij - tol;
            angij_max = angij + tol;
            [idxij_min, idxij_max] = kvec.query(angij_min, angij_max);
            idxij = idxij_min:idxij_max;

            angik_min = angik - tol;
            angik_max = angik + tol;
            [idxik_min, idxik_max] = kvec.query(angik_min, angik_max);
            idxik = idxik_min:idxik_max;

            angjk_min = angjk - tol;
            angjk_max = angjk + tol;
            [idxjk_min, idxjk_max] = kvec.query(angjk_min, angjk_max);
            idxjk = idxjk_min:idxjk_max;

            % create candidate Is and Js
            Is = [kvec.Is(idxij), kvec.Js(idxij)];
            Js = [kvec.Js(idxij), kvec.Is(idxij)];

            % create candidate Ks and verification Is
            Is_prime = [kvec.Is(idxik), kvec.Js(idxik)];
            Ks       = [kvec.Js(idxik), kvec.Is(idxik)];

            % create verification Js and Ks
            Js_prime = [kvec.Is(idxjk), kvec.Js(idxjk)];
            Ks_prime = [kvec.Js(idxjk), kvec.Is(idxjk)];

            % first pair check
            check_Is = ismember(Is, Is_prime);
            check_Js = ismember(Js, Js_prime);
            Is_reduced = Is(check_Is & check_Js);
            Js_reduced = Js(check_Is & check_Js);

            % second pair check
            check_Is_reverse = ismember(Is_prime, Is_reduced);
            check_Ks         = ismember(Ks, Ks_prime);
            Is_prime_reduced = Is_prime(check_Is_reverse & check_Ks);
            Ks_reduced       = Ks(check_Is_reverse & check_Ks);

            % third pair check
            check_Js_reverse = ismember(Js_prime, Js_reduced);
            check_Ks_reverse = ismember(Ks_prime, Ks_reduced);
            Js_prime_reduced = Js_prime(check_Js_reverse & check_Ks_reverse);
            Ks_prime_reduced = Ks_prime(check_Js_reverse & check_Ks_reverse);

            IJKs = [];

            % Three loops to determine whether any of the candidates check
            % out. For a valid triangle, we need that
            % I checks with Iprime, J with Jprime, and K with Kprime.
            for idx_I = 1:length(Is_reduced)
                I = Is_reduced(idx_I);
                J = Js_reduced(idx_I);

                for idx_Iprime = 1:length(Is_prime_reduced)
                    I_prime = Is_prime_reduced(idx_Iprime);
                    K       = Ks_reduced(idx_Iprime);

                    if (I == I_prime)
                        for idx_Jprime = 1:length(Js_prime_reduced)
                            J_prime = Js_prime_reduced(idx_Jprime);
                            K_prime = Ks_prime_reduced(idx_Jprime);

                            if J_prime == J && K_prime == K
                                IJKs(:,end+1) = [I;J;K];
                            end
                        end
                    end
                end
            end

            IJKs = uint64(unique(IJKs', 'rows'))';

            % returns empty
            matches = [];
            att_ICRF2C = sonic.Attitude.empty;
            if isempty(IJKs)
                return
            end

            % if there are valid triads
            u_cam   = pointsS2_meas.u(:,ijk);

            % compute the center of the star triad
            u_triad_center = mean(u_cam, 2);
            u_triad_center = u_triad_center / norm(u_triad_center);

            % [~, locs] = ismember(IJKs, cat_idx');

            for ii = 1:size(IJKs, 2)
                % use kvec to evaluate
                u_world = pointsS2_cat.u(:, IJKs(:,ii));

                att_ICRF2C = sonic.Attitude.solveWahbasProblem( ...
                    sonic.PointsS2(u_cam), sonic.PointsS2(u_world));

                % reproject triad center to ICRF
                u_triad_center_world = att_ICRF2C.dcm'*u_triad_center;

                % reproject all measurement to world
                u_meas_world = att_ICRF2C.dcm'*pointsS2_meas.u;

                % only keep stars near the center of
                % the triad
                idx = (u_triad_center_world'*u_cat_all >= cosmaxangle);
                u_cat_reduced = u_cat_all(:,idx);
                cat_idx_reduced = cat_idx(idx);
                costol = cos(tol);

                % compute the interstar angles between
                % all measurement u and catalog u
                cos_all = u_meas_world'*u_cat_reduced;

                % if any angle are smaller than tolerance store matches
                indices = find(cos_all >=  costol);
                [meas_indices,cat_indices] = ...
                    ind2sub(size(cos_all),indices);

                % quickly prune indices that have double assignment
                % does not look at smallest reprojection error for
                % speed
                [~, prune] = unique(meas_indices);
                meas_indices = meas_indices(prune);
                cat_indices = cat_indices(prune);

                [~, prune2] = unique(cat_indices);
                meas_indices = meas_indices(prune2);
                cat_indices = cat_indices(prune2);

                matches = [meas_indices';
                    cat_idx_reduced(cat_indices)'];


                if size(matches, 2) >= min_matches
                    valid = true;
                    return
                end
            end

        end
    end
end