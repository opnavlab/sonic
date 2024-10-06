classdef (Abstract) Reflectance
    % This software is made available under the MIT license. See SONIC/LICENSE
    % for details.

    properties

    end

    methods (Static)

        function [var1, var2, varargout] = checkInputs(inputType, var1, var2, varargin)
            %% [var1, var2, varargout] = checkInputs(inputType, var1, var2, varargin)
            %   Validates all reflectance model inputs to ensure size
            %   compatibility and value
            %
            %   Inputs:
            %       - inputType (string):
            %               - "iep" for incidence/emission/phase angles
            %               - "iea" for incidence/emission/azimuth angles
            %               - "ien" for incidence/emission/normal vectors
            %                   -- Must be column vectors, (3xn)
            %       - var1
            %       - var2
            %       - varargin: REQUIRED for "ien" as var3
            %                   OPTIONAL for "iea" and "iep" depending on
            %                   the reflectance model
            %
            %   Outputs:
            %       - var1: same as input
            %       - var2: same as input
            %       - varargout: same as varargin
            %                   will be output if varargin is used for var3
            %
            %   Last revised: 09/30/2024
            %   Last author: Priyal Soni

            arguments
                inputType   (1,1)   string {mustBeMember(inputType, {'iep', 'iea', 'ien'})}
                var1        (:,:)   {sonic.Reflectance.mustBeDoubleOrPointsS2(var1)}
                var2        (:,:)   {sonic.Reflectance.mustBeDoubleOrPointsS2(var2)}
            end

            arguments (Repeating)
                varargin            % Must be Double or PointsS2, but comes packaged in a cell due to varargin
            end

            % Ensure varargin is only 0 or 1 input
            if size(varargin{1}, 2) > 1
                error('sonic:Reflectance:invalidInput', ...
                    'Please check your inputs. Maximum of 4 inputs allowed');

                % Check if the cell array is empty (var3 is not needed)
            elseif size(varargin{1}, 2) ~= 0
                var3 = varargin{1}{1};
            end

            % Check var1, var2 and optional var3 parameters based on
            % inputType
            switch inputType
                case "iep"
                    % Check input sizes of var1 and var2 and ensure they
                    % are both within [0, pi/2]
                    if ~(all(size(var1) == size(var2)) && ...
                            all(var1 <= pi/2, 'all') && all(var1 >= 0, 'all') && all(~isnan(var1), 'all') && ...
                            all(var2 <= pi/2, 'all') && all(var2 >= 0, 'all') && all(~isnan(var2), 'all'))
                        error('sonic:Reflectance:invalidInput', ...
                            'Please ensure input sizes are consistent and i, e are within [0, pi/2].');
                    end

                    % If var3 was entered, check that var3 is positive,
                    % satisfies the spherical triangle inequality, matches
                    % the size of the other inputs, and does not contain
                    % NaNs
                    if exist("var3", "var")
                        if ~(all(size(var1) == size(var3)) && ...
                                all(var3 >= 0, 'all') && ...
                                all(var3 <= var1 + var2 + sonic.Tolerances.PhaseAng, 'all') && ...
                                all(~isnan(var3), 'all'))

                            error('sonic:Reflectance:invalidInput', ...
                                'Please ensure size of phase angle input matches i and e and is a valid spherical angle pair.');
                        else
                            varargout = {var3};
                        end
                    end

                case "iea"
                    % Check input sizes of var1 and var2 and ensure they
                    % are both within [0, pi/2]
                    if ~(all(size(var1) == size(var2)) && ...
                            all(var1 <= pi/2, 'all') && all(var1 >= 0, 'all') && all(~isnan(var1), 'all') && ...
                            all(var2 <= pi/2, 'all') && all(var2 >= 0, 'all') && all(~isnan(var2), 'all'))
                        error('sonic:Reflectance:invalidInput', ...
                            'Please ensure input sizes are consistent and i, e are within [0, pi/2].');
                    end

                    % If var3 was entered, check that var3 matches
                    % the size of the other inputs and does not contain
                    % NaNs
                    if exist("var3", "var")
                        if ~all(size(var1) == size(var3)) || all(isnan(var3), 'all')
                            error('sonic:Reflectance:invalidInput', ...
                                'Please ensure size of azimuth angle input matches i and e.');
                        else
                            varargout = {var3};
                        end
                    end

                case "ien"
                    % Ensure normals (u_n) were input
                    if ~exist("var3", "var")
                        error('sonic:Reflectance:invalidInput', ...
                            'Please ensure all inputs are included for input_type = "ien".)');
                    end

                    % Check if all inputs are in PointsS2 class
                    if ~(isa(var1,"sonic.PointsS2") && ...
                            isa(var2,"sonic.PointsS2") && ...
                            isa(var3,"sonic.PointsS2"))

                        error('sonic:Reflectance:invalidInput', ...
                            'Please ensure inputs are in the PointsS2 class.');
                    end

                    % Check if all inputs are same size
                    if ~(var1.n == var2.n && var1.n == var3.n)
                        error('sonic:Reflectance:invalidInput', ...
                            'Please ensure input sizes are consistent.');
                    end

                    % Check is u_n has any NaN
                    if all(isnan(var1.n), 'all') || all(isnan(var2.n), 'all') || all(isnan(var3.n), 'all')
                        error('sonic:Reflectance:invalidInput', ...
                            'Please ensure u_i, u_e, u_n do not contain NaN.');
                    end

                    % Ouput only the unit vectors rather than PointsS2 obj
                    var1 = var1.u;
                    var2 = var2.u;
                    varargout = {var3.u};

                otherwise
                    error('sonic:Reflectance:invalidInput', ...
                        ['Please check your inputs.', ...
                        '\n\t For inputType = "iep" and "iea",'...
                        '\n\t\t var1 = incidence angles', ...
                        '\n\t\t var2 = emission angles', ...
                        '\n\t\t var3 = phase angle (iep) or azimuth angle (iea) (optional)', ...
                        '\n\t For inputType = "ien",'...
                        '\n\t\t var1 = incidence vectors (3xn)', ...
                        '\n\t\t var2 = emission vectors (3xn)', ...
                        '\n\t\t var3 = normal vectors (3xn)']);
            end
        end

        function [i, e] = calcIE(u_i, u_e, u_n)
            %% [i, e] = calcIE(u_i, u_e, u_n)
            %   Calculates incidence and emission angles (measured from the
            %   normal vector) given the incidence, emission and normal
            %   unit vectors.
            %
            %   Inputs:
            %       - u_i (3xn) (double): incidence unit vectors
            %       - u_e (3xn) (double): emission unit vectors
            %       - u_n (3xn) (double): normal unit vectors
            %
            %   Outputs:
            %       - i (1xn) (double): incidence angles [rad]
            %       - e (1xn) (double): emission angles [rad]
            %
            %   Last revised: 09/30/2024
            %   Last author: Priyal Soni

            % Use calcCosIE to calculate cos(i) and cos(e)
            [cosi, cose] = sonic.Reflectance.calcCosIE(u_i, u_e, u_n);

            % Take the inverse cosine for i, e, output
            i = acos(cosi);
            e = acos(cose);
        end

        function [cosi, cose] = calcCosIE(u_i, u_e, u_n)
            %% [cosI, cosE] = calcCosIE(u_i, u_e, u_n)
            %   Calculates cosine of the incidence and emission angles
            %   (measured from the normal vector) given the incidence,
            %   emission and normal unit vectors.
            %
            %   Inputs:
            %       - u_i (3xn) (double): incidence unit vectors
            %       - u_e (3xn) (double): emission unit vectors
            %       - u_n (3xn) (double): normal unit vectors
            %
            %   Outputs:
            %       - cosi (1xn) (double): cosine of incidence angles
            %       - cose (1xn) (double): cosine of emission angles
            %
            %   Last revised: 09/30/2024
            %   Last author: Priyal Soni

            % Calculate the cosine of incidence and emission angles based
            % on the dot product between the vectors and the normals
            cosi = dot(u_i, u_n);
            cose = dot(u_e, u_n);

            % Perform argument check to ensure that any values
            % close to +1 or -1 are set as such to avoid any
            % complex outputs when taking cosine
            cosi = sonic.Reflectance.cosineArgCheck(cosi);
            cose = sonic.Reflectance.cosineArgCheck(cose);
        end

        function g = calcPhase(varargin)
            %% g = calcPhase(varargin)
            %   Calculates the phase angles given either
            %   incidence, emission and azimuth angles or incidence and
            %   emission unit vectors
            %
            %   Inputs:
            %       - varargin:
            %           1. Incidence, emission and azimuth angles (nxm)
            %               (double) [rad]
            %           2. Incidence and emission unit vectors (3xn)
            %               (double)
            %
            %   Outputs:
            %       - g (1xn) (double): phase angles [rad]
            %
            %   Last revised: 09/30/2024
            %   Last author: Priyal Soni

            % Check size of inputs to determine which method to calculate
            % phase.
            if size(varargin, 2) == 3
                % If there are 3 inputs, option (1) from above.
                i = varargin{1};
                e = varargin{2};
                phi = varargin{3};

                % Use calcCosPhase to calcualte cosg
                g = acos(sonic.Reflectance.calcCosPhase(i, e, phi));

            else
                % If there are 2 inputs, option (2) from above.
                u_i = varargin{1};
                u_e = varargin{2};

                % Use calcCosPhase to calcualte cosg
                g = acos(sonic.Reflectance.calcCosPhase(u_i, u_e));
            end
        end

        function cosg = calcCosPhase(varargin)
            %% cosg = calcCosPhase(varargin)
            %   Calculates cosine of the phase angles given either
            %   incidence, emission and azimuth angles or incidence and
            %   emission unit vectors
            %
            %   Inputs:
            %       - varargin:
            %           1. Incidence, emission and azimuth angles (nxm)
            %               (double) [rad]
            %           2. Incidence and emission unit vectors (3xn)
            %               (double)
            %
            %   Outputs:
            %       - cosg (1xn) (double): phase angles [rad]
            %
            %   Last revised: 09/30/2024
            %   Last author: Priyal Soni

            % Check case (1) from above
            if size(varargin, 2) == 3
                i = varargin{1};
                e = varargin{2};
                phi = varargin{3};

                % Calculate cosPhase baseg on i, e, phi
                cosg = cos(i).*cos(e) + sin(i).*sin(e).*cos(phi);

                % Perform argument check to ensure that any values
                % close to +1 or -1 are set as such to avoid any
                % complex outputs when taking cosine
                cosg = sonic.Reflectance.cosineArgCheck(cosg);
            else
                u_i = varargin{1};
                u_e = varargin{2};

                % Calculate cosPhase based on angle between vectors
                cosg = dot(u_i, u_e);

                % Perform argument check to ensure that any values
                % close to +1 or -1 are set as such to avoid any
                % complex outputs when taking cosine
                cosg = sonic.Reflectance.cosineArgCheck(cosg);
            end
        end

        function phi = calcAzimuth(method, var1, var2, var3)
            %% calcAzimuth(method, var1, var2, var3)
            %   Calculates the emission azimuth angles given
            %   either incidence, emission and normal unit vectors or
            %   incidence, emission, and phase angles.
            %
            %   Inputs:
            %       - method (string): Type of azimuth calculation inputs
            %           1. 'unitVec': Unit Vector Inputs
            %           2. 'phaseAng': Phase angle with incidence and
            %               emission angle inputs [rad]
            %       - var1/var2/var3:
            %           1. Incidence, emission and normal unit vectors
            %               (3xn) (double)
            %           2. Incidence, emission, phase angles (nxm)
            %               (double) [rad]
            %
            %   Outputs:
            %       - cosPhi (1xn) (double): cosine of the emission azimuth
            %           angles
            %
            %   Last revised: 09/30/2024
            %   Last author: Priyal Soni

            arguments
                method  (1,1)   string {mustBeMember(method, {'unitVec', 'phaseAng'})}
                var1    (:,:)   double
                var2    (:,:)   double
                var3    (:,:)   double
            end

            % Use calcCosAzimuth to calculate cos(phi)
            switch method
                case 'unitVec'
                    [cosPhi, quadCheck] = sonic.Reflectance.calcCosAzimuth(method, var1, var2, var3);
                    phi = acos(cosPhi);
                    phi(quadCheck) = -phi(quadCheck);
                otherwise
                    phi = acos(sonic.Reflectance.calcCosAzimuth(method, var1, var2, var3));
            end
        end

        function [cosPhi, varargout] = calcCosAzimuth(method, var1, var2, var3)
            %% cosPhi = calcCosAzimuth(method, var1, var2, var3)
            %   Calculates cosine of the emission azimuth angles given
            %   either incidence, emission and normal unit vectors or
            %   incidence, emission, and phase angles.
            %
            %   Inputs:
            %       - method (string): Type of azimuth calculation inputs
            %           1. 'unitVec': Unit Vector Inputs
            %           2. 'phaseAng': Phase angle with incidence and
            %               emission angle inputs
            %       - var1/var2/var3:
            %           1. Incidence, emission and normal unit vectors
            %               (3xn) (double)
            %           2. Incidence, emission, phase angles (nxm)
            %               (double) [rad]
            %
            %   Outputs:
            %       - cosPhi (1xn) (double): cosine of the emission azimuth
            %           angles
            %
            %   Last revised: 09/30/2024
            %   Last author: Priyal Soni


            arguments
                method  (1,1)   string {mustBeMember(method, {'unitVec', 'phaseAng'})}
                var1    (:,:)   double
                var2    (:,:)   double
                var3    (:,:)   double
            end

            switch method
                case 'unitVec'
                    u_i = var1;
                    u_e = var2;
                    u_n = var3;

                    % Define local xyz frame for each patch (n = z, x = u_i in plane)
                    dot_uin = dot(u_i, u_n);
                    u_i_plane = u_i - u_n.*dot_uin;
                    x = u_i_plane./vecnorm(u_i_plane, 2, 1);

                    % Project the incidence and emission vectors onto the local xy plane
                    dot_uen = dot(u_e, u_n);
                    u_e_plane = u_e - u_n.*dot_uen;
                    u_e_plane = u_e_plane./vecnorm(u_e_plane, 2, 1);

                    % Use projection to calculate cosine of orientation angle phi_e
                    cosPhi = dot(u_e_plane, x);

                    % Determine if any u_i == u_n
                    % If so, the plane of incidence does not exist
                    % and azimuth is measured from the plane of emission
                    zeroIncidence_mask = all(abs(u_i-u_n) < sonic.Tolerances.SmallAngle);
                    cosPhi(zeroIncidence_mask) = 1;

                    % Check if any u_e == u_i. Then the azimuth angle is
                    % automatically zero.
                    zeroEmission_mask = all(abs(u_e-u_i) < sonic.Tolerances.SmallAngle);
                    cosPhi(zeroEmission_mask) = 1;

                    % Perform argument check to ensure that any values
                    % close to +1 or -1 are set as such to avoid any
                    % complex outputs when taking cosine
                    cosPhi = sonic.Reflectance.cosineArgCheck(cosPhi);

                    % quadrant check for azimuth angle
                    quadCheck = dot(cross(x, u_e_plane), u_n) < 0;
                    varargout = {quadCheck};

                case 'phaseAng'
                    i = var1;
                    e = var2;
                    g = var3;

                    % Calculate cosPhi based on spherical trigonometry
                    cosPhi = (cos(g) - cos(i).*cos(e))./(sin(i).*sin(e));

                    % Perform argument check to ensure that any values
                    % close to +1 or -1 are set as such to avoid any
                    % complex outputs when taking cosine
                    cosPhi = sonic.Reflectance.cosineArgCheck(cosPhi);

                    % Determine if any u_i == u_n
                    % If so, the plane of incidence does not exist
                    % and azimuth is measured from the plane of emission
                    zeroIncidence_mask = all(abs(i) < sonic.Tolerances.SmallAngle);
                    cosPhi(zeroIncidence_mask) = 1;

                    % Check if any e == i. Then the azimuth angle is
                    % automatically zero.
                    zeroEmission_mask = all(abs(i-e) < sonic.Tolerances.SmallAngle);
                    cosPhi(zeroEmission_mask) = 1;

                otherwise
                    error('sonic:Reflectance:invalidInput', ...
                        ['Please choose either...',...
                        '\n\t1. ''unitVec'' to calculate azimuth angles given unit vectors for incidence (var1), emission (var2) and normal (var3) vectors', ...
                        '\n\t2. ''phaseAng'' to calculate azimuth angles given incidence (var1), emission (var2) and phase (var3) angles']);
            end
        end

        function x = cosineArgCheck(x)
            %% x = cosineArgCheck(x)
            %   Checks x to ensure any value sufficiently close to +1 or -1
            %   are set to +1 and -1 (respectively) to avoid complex
            %   evaluation in cos(x)
            %
            %   Inputs:
            %       - x (n,m) (double): matrix of values that will be
            %       evaluated using cos(x)
            %
            %   Outputs:
            %       - x (n,m) (double): corrected matrix with any values
            %       within sonic.Tolerances.CosArgCheck and greater than 1
            %       are set to +1, and similar with -1.
            %
            %   Last revised: 09/30/2024
            %   Last author: Priyal Soni

            % Perform cosine argument check before acos()
            l1 = abs(x-1.0) < sonic.Tolerances.CosArgCheck & x > 1; % Checks if it's close enough to +1
            l2 = abs(x+1.0) < sonic.Tolerances.CosArgCheck & x < -1; % Checks if it's close enough to -1

            x(l1) = 1.0;
            x(l2) = -1.0;
        end

        function mustBeDoubleOrPointsS2(a)
            %% mustBeDoubleOrPointsS2(a)
            %   Input validation function for reflectance model inputs.
            %
            %   Inputs:
            %       - a: typically var1, var2, or var3 in a reflctance
            %           object refl() function.
            %
            %   Outputs:
            %       - N/A for input validation functions
            %
            %   Last revised: 09/30/2024
            %   Last author: Priyal Soni

            if ~(isa(a,'double') || isa(a,'sonic.PointsS2'))
                error('sonic:Reflectance:invalidInput', ...
                    'Please ensure input class is double or sonic.PointsS2');
            end
        end
    end

end