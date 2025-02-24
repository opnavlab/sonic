% This software is made available under the MIT license. See SONIC/LICENSE
% for details.
classdef Attitude

    properties (SetAccess = protected)
        % Quaternion representation of attitude (this is what we will store
        % all attitudes as, fundamentally)
        quat_v          (3, 1)  double
        quat_s          (1, 1)  double
    end

    properties (SetAccess = private)
        % These are the "derivative" forms of attitude that we want users
        % to be able to access. Rather than using dependent properties,
        % which would recompute these every time on call, we want to
        % precompute these ahead of time and have them ready for access.

        dcm             (3, 3)  double
        rotation_vec    (3, 1)  double
    end

    methods

        function obj = Attitude(raw_att, quat_flag)
            %% obj = Attitude(raw_att, quat_flag)
            %   Constructs an Attitude object, which parameterizes different
            %   representations of the same attitude. Note that this class
            %   adopts the passive representation of Attitude.
            %
            %   Inputs:
            %       - raw_att: A valid parameterization of attitude. Valid
            %         input representations are:
            %           - (3x1 double): Rotation vector
            %           - (4x1 double): Quaternion
            %           - (3x3 double): Direction Cosine Matrix (DCM)
            %         The type of attitude parameterization is inferred from the
            %         number of elements in the supplied parameterization.
            %       - quat_flag (char/string): If inputting a quaternion, this
            %         is REQUIRED. Specify 'vs' or 'sv' to indicate the
            %         convention of the quaternion as scalar-last or
            %         scalar-first, respectively.
            %
            %   Outputs:
            %       - obj (sonic.Attitude): Attitude object, encoding a
            %         particular attitude which can be viewed in multiple
            %         parameterizations.
            %
            %   Last revised: 2/15/24
            %   Last author: Michael Krause

            arguments
                raw_att
                quat_flag = []
            end

            % For now, let's accept inputs of a DCM, rotation vector, or
            % quaternion. In each case, we should sanitize the user input,
            % and calculate the other representations of attitude.
            [m,n] = size(raw_att);

            if (m==3 && n==1)
                % This is a rotation vector.
                % Check to make sure it is a valid rotation vector.
                sonic.Attitude.verifyValidRotationVec(raw_att);
                obj.rotation_vec = raw_att;

                % Then convert to quaternion/DCM:
                [qv, qs] = sonic.Attitude.rotationVecToQuat(raw_att);
                obj.quat_v = qv;
                obj.quat_s = qs;

                obj.dcm = sonic.Attitude.rotationVecToDCM(raw_att);

            elseif (m==4 && n==1)
                % This is a quaternion.
                % Check to make sure it is a valid quaternion.
                sonic.Attitude.verifyValidQuat(raw_att);

                if ~isempty(quat_flag)
                    switch quat_flag
                        case 'sv'
                            % Scalar-first:
                            qv = raw_att(2:4);
                            qs = raw_att(1);
                        case 'vs'
                            % Scalar-last:
                            qv = raw_att(1:3);
                            qs = raw_att(4);
                        otherwise
                            error('sonic:Attitude:incorrectQuatConvention', ...
                                ['Must specify quaternion convention as one' ...
                                'of ''vs'' (vector/scalar, i.e., vector-first) ' ...
                                'or ''sv'' (scalar/vector, i.e., scalar-first).']);
                    end
                else
                    error('sonic:Attitude:enterQuatConvention', ...
                        ['No quaternion convention supplied. When ' ...
                        'entering a quaternion, must specify ' ...
                        'quaternion convention as one ' ...
                        'of ''vs'' (vector/scalar, i.e., vector-first) ' ...
                        'or ''sv'' (scalar/vector, i.e., scalar-first).']);
                end

                % Force vector to be a column:
                obj.quat_v = qv(:);
                obj.quat_s = qs;

                % Then convert to rotation vector/DCM:
                obj.dcm = sonic.Attitude.quatToDCM(obj.quat_v, obj.quat_s);
                obj.rotation_vec = sonic.Attitude.quatToRotationVec(obj.quat_v, obj.quat_s);

            elseif (m==3 && n==3)
                % This is a DCM.
                % Check to make sure it is a valid DCM.
                sonic.Attitude.verifyValidDCM(raw_att);

                obj.dcm = raw_att;

                % Then convert to quaternion/rotation vector:
                [qv, qs] = sonic.Attitude.DCMtoQuat(raw_att);
                obj.quat_v = qv;
                obj.quat_s = qs;

                obj.rotation_vec = sonic.Attitude.DCMtoRotationVec(raw_att);

            else
                % Toss an error here.
                error('sonic:Attitude:invalidInput', ...
                    ['Attitude representation not recognized. ' ...
                    'Must enter a 3x1-element rotation vector, '...
                    '4x1-element quaternion, or '...
                    '3x3 DCM.']);

            end

        end

        function new_att = comp(obj, obj2)
            %% new_att = comp(obj, obj2)
            %
            %   Composes attitudes together. This may be thought of in terms of
            %   DCMs as:
            %
            %               new_att.dcm = obj.dcm*t_obj.dcm;
            %
            %   Inputs:
            %       - obj (1x1 sonic.Attitude): Attitude object
            %       - obj2 (1x1 sonic.Attitude): Another attitude object to
            %         compose with the calling object
            %
            %   Outputs:
            %       - new_att (1x1 sonic.Attitude): The resultant composed
            %         attitude object
            %
            %   Last revised: 3/15/24
            %   Last author: Michael Krause


            % Carry out composition in terms of quaternions.
            % Note that q_mult(...) automatically renormalizes the
            % resultant quaternion, so we don't need to worry about that
            % here.
            [qv_comp, qs_comp] = sonic.Attitude.q_mult(...
                obj.quat_v, obj.quat_s, obj2.quat_v, obj2.quat_s);
            new_att = sonic.Attitude([qv_comp; qs_comp], 'vs');
        end

        function new_att = inv(obj)
            %% new_att = inv(obj)
            %
            %   Inverses attitude. This may be thought of in terms of
            %   DCMs as:
            %
            %               new_att.dcm = obj.dcm';
            %
            %   Inputs:
            %       - obj (1x1 sonic.Attitude): Attitude object
            %
            %   Outputs:
            %       - new_att (1x1 sonic.Attitude): The resultant inversed
            %         attitude object
            %
            %   Last revised: 10/29/2024
            %   Last author: Sebastien Henry

            new_att = sonic.Attitude(obj.dcm');
        end

        function [rot_obj] = rotate(obj, to_rot)
            %% [rot_obj] = rotate(to_rot, rot_att)
            %
            %   Rotates an object relative to the observer. Currently only
            %   supports rotating points in P3, i.e., sonic.Points3 and
            %   sonic.PointsS2.
            %
            %   Inputs:
            %       - to_rot (1x1 sonic.Points3 or sonic.PointsS2): Object to
            %         rotate. Currently only supports rotating 3D points.
            %       - rot_att (1x1 sonic.Attitude): Object to Camera attitude.
            %
            %   Outputs:
            %       - rot_obj (1x1 sonic.Points3 or sonic.PointsS2): Rotated
            %         object. Return type will match the input type.
            %
            %   Last revised: 11/08/24
            %   Last author: Sebastien Henry

            arguments
                obj         (1, 1)      sonic.Attitude
                to_rot      (1, 1)  % Make this sonic.GeometryP3 eventually
            end

            % For now just support points.
            if isa(to_rot, 'sonic.PointsS2') || isa(to_rot, 'sonic.Points3')

                % Pull out the points depending on the object:
                if isa(to_rot, 'sonic.Points3')
                    if to_rot.has_inf_points
                        % If infinite points, need to use 4x1 p3
                        raw_pts = to_rot.p3;
                        rot_mat = [obj.dcm, zeros(3,1); zeros(1,3), 1];
                    else
                        % Otherwise 3x1 r3 suffices
                        raw_pts = to_rot.r3;
                        rot_mat = obj.dcm;
                    end
                else % is PointsS2
                    raw_pts = to_rot.u;
                    rot_mat = obj.dcm;
                end

                % Do the actual rotation:
                rot_pts = rot_mat*raw_pts;

                % Package up result commensurate with input type:
                if isa(to_rot, 'sonic.Points3')
                    rot_obj = sonic.Points3(rot_pts);
                else % is PointsS2
                    rot_obj = sonic.PointsS2(rot_pts);
                end

            else
                error('sonic:Attitude:rotate:invalidType', ...
                    ['Attempted to rotate invalid type. Valid ' ...
                    'types are sonic.Points3 or sonic.PointsS2.']);
            end
        end


        function multiplied_obj = mtimes(obj, to_mult)
            %% multiplied_obj = mtimes(obj, to_mult)
            %   Pre-multiply the att to an object result depending on the
            %   type of to_mult
            %
            %   Inputs:
            %       - to_mult: a SONIC object, either
            %           - (1x1 sonic.PointsS2): In this case, the output is
            %           the attitude matrix times the points s2
            %           - (1x1 sonic.Points3): In this case, the output is
            %           the attitude matrix times the points p3
            %           - (1x1 sonic.Attitude): In this case, the output is
            %           the compose of the attitude matrices
            %   Last revised: 11/07/2024
            %   Last author: Jennifer Nolan

            switch class(to_mult)
                case 'sonic.Points3'
                    multiplied_obj = obj.rotate(to_mult);
                case 'sonic.PointsS2'
                    multiplied_obj = obj.rotate(to_mult);
                case 'sonic.Attitude'
                    multiplied_obj = obj.comp(to_mult);
                otherwise
                    error('sonic:Attitude:mtimes:invalidInputType', ...
                        ['Only Points3, PointsS2, and Attitude can be ' ...
                        'multiplied by an Attitude'])
            end

        end

    end

    %% Conversion Functions
    methods (Static)

        %% Conversion functions to quaternion:

        function [qv, qs] = DCMtoQuat(T)
            %% [qv, qs] = DCMtoQuat(T)
            %   Converts a DCM to a quaternion.
            %
            %   Inputs:
            %       - T (3x3 double): Valid DCM
            %   Outputs:
            %       - qv (3x1 double): Vector component of the quaternion.
            %       - qs (1x1 double): Scalar component of the quaternion.
            %
            %   Last revised: 2/15/24
            %   Last author: Michael Krause

            % DO SMALL ANGLE APPROX HERE!

            [ax, ang_RAD] = sonic.Attitude.DCMtoAxAng(T);
            [qv, qs] = sonic.Attitude.axAngToQuat(ax, ang_RAD);
        end

        function [qv, qs] = rotationVecToQuat(rotation_vec)
            %% [qv, qs] = rotationVecToQuat(rotation_vec)
            %   Converts an Euler rotation vector to a quaternion.
            %
            %   Inputs:
            %       - rotation_vec: (3x1 double): Valid rotation vector
            %   Outputs:
            %       - qv (3x1 double): Vector component of the quaternion.
            %       - qs (1x1 double): Scalar component of the quaternion.
            %
            %   Last revised: 2/29/24
            %   Last author: John Christian

            ang_RAD = norm(rotation_vec);

            if ang_RAD < sonic.Tolerances.SmallAngle
                qs = sqrt( 1 - ang_RAD^2/4 );
                qv = rotation_vec / 2;
            else
                ax = rotation_vec./ang_RAD;
                [qv, qs] = sonic.Attitude.axAngToQuat(ax, ang_RAD);
            end

        end

        function [qv, qs] = axAngToQuat(ax, ang_RAD)
            %% [qv, qs] = axAngToQuat(ax, ang_RAD)
            %   Converts an Euler axis/angle to a quaternion.
            %
            %   Inputs:
            %       - ax: (3x1 double): Valid Euler rotation axis
            %       - ang_RAD: (1x1 double): Valid Euler rotation angle, in
            %         radians.
            %   Outputs:
            %       - qv (3x1 double): Vector component of the quaternion.
            %       - qs (1x1 double): Scalar component of the quaternion.
            %
            %   Last revised: 2/15/24
            %   Last author: Michael Krause

            qv = sin(ang_RAD/2).*ax;
            qs = cos(ang_RAD/2);
        end

        %% Conversion functions to DCM:

        function T = quatToDCM(qv, qs)
            %% T = quatToDCM(qv, qs)
            %   Converts a quaternion (expressed as separate vector/scalar
            %   parts) to a direction cosine matrix (DCM).
            %
            %   Inputs:
            %       - qv (3x1 double): Vector component of the quaternion.
            %       - qs (1x1 double): Scalar component of the quaternion.
            %   Outputs:
            %       - T (3x3 double): DCM representation of attitude.
            %
            %   Last revised: 2/15/24
            %   Last author: Michael Krause

            % Rodrigues' formula for quaternions:
            T = (qs^2 - qv'*qv)*eye(3) - 2*qs*sonic.Math.crossmat(qv) + 2*(qv*qv');

        end

        function T = rotationVecToDCM(rotation_vec)
            %% T = rotationVecToDCM(rotation_vec)
            %   Converts an Euler rotation vector to a direction cosine
            %   matrix (DCM).
            %
            %   Inputs:
            %       - rotation_vec (3x1 double): Valid rotation vector
            %   Outputs:
            %       - T (3x3 double): DCM representation of attitude.
            %
            %   Last revised: 2/15/24
            %   Last author: Michael Krause

            if norm(rotation_vec) < sonic.Tolerances.SmallAngle
                T = eye(3) - sonic.Math.crossmat(rotation_vec);
            else
                % Decompose into Euler ax/ang and then use Rodrigues' to
                % assemble T.
                [ax, ang_RAD] = sonic.Attitude.rotationVecToAxAng(rotation_vec);

                % Rodrigues' formula for ax/ang:
                T = sonic.Attitude.axAngToDCM(ax, ang_RAD);
            end

        end

        function T = axAngToDCM(ax, ang_RAD)
            %% T = axAngToDCM(ax, ang_RAD)
            %   Converts an Euler axis/angle pair to a direction cosine
            %   matrix (DCM).
            %
            %   Inputs:
            %       - ax (3x1 double): Valid Euler rotation axis
            %       - ang_RAD (1x1 double): Valid Euler rotation angle, in
            %           radians.
            %   Outputs:
            %       - T (3x3 double): DCM representation of attitude.
            %
            %   Last revised: 2/15/24
            %   Last author: Michael Krause

            % Euler finite rotation formula for ax/ang:
            T = cos(ang_RAD)*eye(3) - sin(ang_RAD)*sonic.Math.crossmat(ax) + (1 - cos(ang_RAD))*(ax*ax');

        end

        %% Conversion functions to rotation vectors:

        function rotation_vec = quatToRotationVec(qv, qs)
            %% rotation_vec = quatToRotationVec(qv, qs)
            %   Converts a quaternion (expressed as separate vector/scalar
            %   parts) to an Euler rotation vector.
            %
            %   Inputs:
            %       - qv (3x1 double): Vector component of the quaternion.
            %       - qs (1x1 double): Scalar component of the quaternion.
            %   Outputs:
            %       - rotation_vec: (3x1 double): Euler rotation vector.
            %
            %   Last revised: 2/15/24
            %   Last author: Michael Krause

            % SMALL ANGLE HERE???

            [ax, ang_RAD] = sonic.Attitude.quatToAxAng(qv, qs);
            rotation_vec = sonic.Attitude.axAngToRotationVec(ax, ang_RAD);
        end

        function rotation_vec = DCMtoRotationVec(T)
            %% rotation_vec = DCMtoRotationVec(T)
            %   Converts a direction cosine matrix (DCM) to an Euler rotation
            %   vector.
            %
            %   Inputs:
            %       - T (3x3 double): Valid DCM
            %   Outputs:
            %       - rotation_vec (3x1 double): Euler rotation vector.
            %
            %   Last revised: 2/15/24
            %   Last author: Michael Krause

            % DEFINITELY NEED SMALL ANGLE HERE!!!! (This should be an easy
            % one too).

            % Use the inverse of Rodrigues' much like we do in the quat
            % formulation.
            [ax, ang_RAD] = sonic.Attitude.DCMtoAxAng(T);
            rotation_vec = sonic.Attitude.axAngToRotationVec(ax, ang_RAD);
        end

        function rotation_vec = axAngToRotationVec(ax, ang)
            %% rotation_vec = axAngToRotationVec(ax, ang)
            %   Converts an Euler axis/angle pair (angle in radians) to an
            %   Euler rotation vector.
            %
            %   Inputs:
            %       - ax (3x1 double): Valid Euler rotation axis
            %       - ang_RAD (1x1 double): Valid Euler rotation angle, in
            %         radians.
            %   Outputs:
            %       - rotation_vec: (3x1 double): Euler rotation vector.
            %
            %   Last revised: 2/15/24
            %   Last author: Michael Krause

            rotation_vec = ax*ang;
        end

        %% Conversion functions to axis/angle:

        function [ax, ang_RAD] = DCMtoAxAng(T)
            %% [ax, ang_RAD] = DCMtoAxAng(T)
            %   Converts a direction cosine matrix (DCM) to an Euler axis/angle
            %   pair (angle in radians).
            %
            %   Inputs:
            %       - T (3x3 double): Valid DCM
            %   Outputs:
            %       - ax (3x1 double): Euler rotation axis
            %       - ang_RAD (1x1 double): Euler rotation angle, in radians
            %
            %   Last revised: 2/15/24
            %   Last author: Michael Krause

            ang_RAD = acos((trace(T) - 1)/2);
            ax = (1/(2*sin(ang_RAD)))*[ T(2, 3) - T(3, 2);
                T(3, 1) - T(1, 3);
                T(1, 2) - T(2, 1)];

            if abs(ang_RAD) < sonic.Tolerances.SmallAngle
                ax = [0;0;0];
            end

        end

        function [ax, ang_RAD] = rotationVecToAxAng(rotation_vec)
            %% [ax, ang_RAD] = rotationVecToAxAng(rotation_vec)
            %   Converts an Euler rotation vector to an Euler axis/angle
            %   pair (angle in radians).
            %
            %   Inputs:
            %       - rotation_vec (3x1 double): Euler rotation vector.
            %   Outputs:
            %       - ax (3x1 double): Euler rotation axis
            %       - ang_RAD (1x1 double): Euler rotation angle, in radians
            %
            %   Last revised: 2/15/24
            %   Last author: Michael Krause

            ang_RAD = norm(rotation_vec);

            if ang_RAD < sonic.Tolerances.SmallAngle
                ax = [0;0;0];
            else
                ax = rotation_vec./ang_RAD;
            end

        end

        function [ax, ang_RAD] = quatToAxAng(qv, qs)
            %% [ax, ang_RAD] = quatToAxAng(qv, qs)
            %   Converts a quaternion (expressed as separate vector/scalar
            %   parts) to an Euler axis/angle pair (angle in radians).
            %
            %   Inputs:
            %       - qv (3x1 double): Vector component of the quaternion.
            %       - qs (1x1 double): Scalar component of the quaternion.
            %   Outputs:
            %       - ax (3x1 double): Euler rotation axis
            %       - ang_RAD (1x1 double): Euler rotation angle, in radians
            %
            %   Last revised: 2/15/24
            %   Last author: Michael Krause

            ang_RAD = 2*atan2(norm(qv), qs);

            % Then pull out the axis:
            if abs(ang_RAD) < sonic.Tolerances.SmallAngle
                ax = [0;0;0];       % Update this.
            else
                ax = qv./norm(qv);
            end

        end

    end

    %% Attitude util functions
    methods (Static)
        function ang_RAD = angleBetween(att1, att2)
            %% angle = angleBetween(att1, att2)
            % Computes the angle between two attitude objects
            %
            %   Inputs:
            %       - att1 (1x1 sonic.Attitude): attitude object 1
            %       - att2 (1x1 sonic.Attitude): attitude object 2
            %   Outputs:
            %       - ang_RAD (1x1 double): angle between the two attitudes
            %
            %   Last revised: 4/15/24
            %   Last author: Sebastien Henry

            % Compose one attitude with the inverse of the other attitude
            att_comp = att1 * att2.inv();

            % Compute the angle from the quaternion
            ang_RAD = 2*acos(att_comp.quat_s);
        end

        function att = solveWahbasProblem(pointsS2_cam, pointsS2_world)
            %% att = solveWahbasProblem(pointsS2_cam, pointsS2_world)
            %   Compute attitude from world to camera using directions in
            %   the camera frame and corresponding directions in the world
            %   frame. Requires at least 2 points.
            %
            %   Inputs:
            %       - pointsS2_cam (1x1 sonic.PointsS2): directions in the
            %         camera frame
            %       - pointsS2_world (1x1 sonic.PointsS2): directions in the
            %         world frame
            %   Outputs:
            %       - att (1x1 sonic.Attitude): attitude object that goes
            %         from world to camera
            %   Last revised: 4/15/24
            %   Last author: Sebastien Henry

            arguments
                pointsS2_cam    (1, 1) sonic.PointsS2
                pointsS2_world  (1, 1) sonic.PointsS2
            end

            if pointsS2_cam.n < 2 || pointsS2_world.n < 2 || ...
                    (pointsS2_world.n ~= pointsS2_cam.n)

                error('sonic:Attitude:solveWahbasProblem:badSizes', ...
                    ['Must supply at least 2 matching points in the ' ...
                    'camera and world frames. Number of points ' ...
                    'in each frame must correspond.']);
            end

            [U, ~, V] = svd(pointsS2_cam.u*pointsS2_world.u');
            M = diag([1, 1, det(U)*det(V)]);
            T = U*M*V';

            att = sonic.Attitude(T);
        end

        function dvecDCM_drotVec = dDCMdRotationVec(rotVec)
            % Compute the partial of a vectorized DCM with respect to a rotation
            % vector.

            % Decompose rotationVec into axis / angle
            [n, theta] = sonic.Attitude.rotationVecToAxAng(rotVec);

            % Compute dvec(crossmat(n)) / dn
            dveccrossn_dn = [
                0,  0,  0;
                0,  0,  1;
                0, -1,  0;
                0,  0, -1;
                0,  0,  0;
                1,  0,  0;
                0,  1,  0;
                -1,  0,  0;
                0,  0,  0];

            % If small angle, return the simplified version
            if theta < sonic.Tolerances.SmallAngle
                % approx R = I - crossmat(rotVec)
                dvecDCM_drotVec = -dveccrossn_dn;
                return
            end

            % some pre-computations
            ctheta = cos(theta);
            stheta = sin(theta);
            crossn = sonic.Math.crossmat(n);
            outern = n*n';

            % Partial of outer product
            n1 = n(1);
            n2 = n(2);
            n3 = n(3);
            dvecoutern_dn = [
                2*n1,    0,    0;
                n2,     n1,    0;
                n3,      0,   n1;
                n2,     n1,    0;
                0,    2*n2,    0;
                0,      n3,   n2;
                n3,      0,   n1;
                0,      n3,   n2;
                0,       0, 2*n3];

            % Partial of DCM with respect to axis
            dvecDCM_dn = -stheta*dveccrossn_dn + (1-ctheta)*dvecoutern_dn;

            % Partial of axis with respect to rotation vector
            dn_drotVec = 1/theta;

            % Partial of axis with respect to angle
            dn_dtheta = -rotVec/theta^2;

            % Partial of DCM with respect to angle
            I3 = eye(3);
            dvecDCM_dtheta = -stheta*I3(:) - ctheta*crossn(:) + stheta*outern(:);

            % Partial of angle with respect to rotation vector
            dtheta_drotVec = n';

            % Partial of DCM with respect to rotation vector
            % Application of the chain rule
            dvecDCM_drotVec = dvecDCM_dtheta*dtheta_drotVec + ...
                dvecDCM_dn*dn_drotVec + dvecDCM_dn*(dn_dtheta*dtheta_drotVec);
        end


    end

    %%  Verification functions:

    methods (Static)

        function verifyValidDCM(cand_dcm)
            %% verifyValidDCM(cand_dcm)
            %   Verifies that a DCM is valid. Errors if invalid, runs to
            %   completion if valid. Currently, the validity checks are:
            %       - Ensure determinant is +1.
            %
            %   Inputs:
            %       - cand_dcm (3x3 double): Possible DCM
            %
            %   Last revised: 2/15/24
            %   Last author: Michael Krause

            % Want to implement some (speedy) checks to ensure the user has
            % provided sane inputs.

            % Verify dcm*dcm' is identity.

            % For now, just do a quick check that we have a +1 det.
            % Depending on how this works out, may want to adjust this to
            % be within a certain tolerance of 1.
            if abs(det(cand_dcm)-1) > sonic.Tolerances.DCMDetOne
                error('sonic:Attitude:verifyValidDCM:DetNot1', ...
                    'A valid right-handed DCM should have determinant = +1.');
            end

            if abs(trace((cand_dcm*cand_dcm')-eye(3))) > sonic.Tolerances.DCMOrthogCheck
                error('sonic:Attitude:verifyValidDCM:NotOrthog', ...
                    'Columns of a DCM must be orthogonal (i.e., T*T'' = I)');
            end

        end

        function verifyValidQuat(cand_quat)
            %% verifyValidQuat(cand_quat)
            %   Verifies that a quaternion is valid. Errors if invalid, runs to
            %   completion if valid. Currently, the validity checks are:
            %       - Ensure norm of quaternion is 1.
            %
            %   Inputs:
            %       - cand_quat (4x1 double): Possible quaternion
            %
            %   Last revised: 2/15/24
            %   Last author: Michael Krause

            % For now, just check for a norm of 1. Again, likely want to
            % change this to be a tolerance-based check.
            if abs(norm(cand_quat)-1) > sonic.Tolerances.QuatNormOne
                error('sonic:Attitude:verifyValidQuat:NormNot1', ...
                    'A valid unit quaternion should have norm = 1.');
            end

        end

        function verifyValidRotationVec(cand_rotation_vec)
            %% verifyValidRotationVec(cand_rotation_vec)
            %   Verifies that a rotation vector is valid. Errors if invalid,
            %   runs to completion if valid. Currently, there are no validity
            %   checks, and this serves as a placeholder.
            %
            %   Inputs:
            %       - cand_rotation_vec (3x1 double): Possible rotation vector
            %
            %   Last revised: 2/15/24
            %   Last author: Michael Krause

            % No screening for now.
            if ~isreal(cand_rotation_vec)
                error('sonic:Attitude:verifyValidRotationVec:NotReal', ...
                    ['Imaginary entries found. A valid rotation vector' ...
                    ' should be real.']);
            end
            if ~allfinite(cand_rotation_vec)
                error('sonic:Attitude:verifyValidRotationVec:NotFinite', ...
                    ['NaN or Inf entries found. A valid rotation vector' ...
                    ' should be finite.']);
            end
        end

    end

    %% Spherical Linear Interpolation (SLERP):
    methods (Static)

        function interp_att = interpolate(att1, att2, t)
            %% interp_att = interpolate(att1, att2, t)
            %
            %   Given two Attitude objects and a parameter for the fraction of
            %   arclength between the attitudes, interpolates the desired
            %   intermediate attitude. Internally, uses a quaternion-based
            %   spherical linear interpolation (SLERP) routine.
            %
            %   Inputs:
            %       - att1 (1x1 sonic.Attitude): an Attitude object
            %       - att2 (1x1 sonic.Attitude): an Attitude object
            %       - t (1x1 double): a fraction of the spherical arclength
            %         between `att1` and `att2` at which to interpolate the
            %         attitude.
            %
            %   Outputs:
            %       - interp_att (1x1 sonic.Attitude): a new Attitude object
            %         representing the desired interpolated attitude.
            %
            %   Last revised: 4/18/24
            %   Last author: Michael Krause

            % Extract quaternion components:
            qv1 = att1.quat_v;
            qs1 = att1.quat_s;
            qv2 = att2.quat_v;
            qs2 = att2.quat_s;

            % Call slerp
            [qvi, qsi] = sonic.Attitude.q_slerp(qv1, qs1, qv2, qs2, t);

            % Stuff back into a new attitude object:
            interp_att = sonic.Attitude([qvi; qsi], 'vs');

        end

    end

    % Hidden support methods (keep quaternion math abstracted away)
    methods (Static, Hidden)

        function [qvr, qsr] = q_slerp(qv1, qs1, qv2, qs2, t)
            %% [qvr, qsr] = q_slerp(qv1, qs1, qv2, qs2, t)
            %
            %   Performs a quaternion spherical linear interpolation (SLERP)
            %   routine.
            %
            %   See Pg. 248 of ACM SIGGRAPH Computer Graphics, Vol 19, Issue 3
            %   "Animating rotation with quaternion curves" by Ken Shoemake,
            %   https://doi.org/10.1145/325165.325242
            %
            %   Inputs:
            %       - qv1 (3x1 double): Vector part of the first quaternion
            %       - qs1 (1x1 double): Scalar part of the first quaternion
            %       - qv2 (3x1 double): Vector part of the second quaternion
            %       - qs2 (1x1 double): Scalar part of the second quaternion
            %       - t (1x1 double): a fraction of the spherical arclength
            %         between the two quaternions at which to interpolate the
            %         resultant quaternion.
            %
            %   Outputs:
            %       - qvr (3x1 double): Vector part of the interpolated
            %         quaternion.
            %       - qsr (1x1 double): Scalar part of the interpolated
            %         quaternion.
            %
            %   Last revised: 4/18/24
            %   Last author: Michael Krause

            q1 = [qv1; qs1];
            q2 = [qv2; qs2];

            theta = acos(q1'*q2);
            sin_th = sin(theta);

            q_int = (sin((1 - t)*theta)/sin_th)*q1 + (sin(t*theta)/sin_th)*q2;

            qvr = q_int(1:3);
            qsr = q_int(4);

        end

        function [qvr, qsr] = q_mult(qv1, qs1, qv2, qs2)
            %% [qvr, qsr] = q_mult(qv1, qs1, qv2, qs2)
            %
            %   Multiplies two quaternions, in the non-Hamiltonian sense.
            %
            %   Inputs:
            %       - qv1 (3x1 double): Vector part of the first quaternion
            %       - qs1 (1x1 double): Scalar part of the first quaternion
            %       - qv2 (3x1 double): Vector part of the second quaternion
            %       - qs2 (1x1 double): Scalar part of the second quaternion
            %
            %   Outputs:
            %       - qvr (3x1 double): Vector part of the multiplied
            %         quaternion.
            %       - qsr (1x1 double): Scalar part of the multiplied
            %         quaternion.
            %
            %   Last revised: 4/18/24
            %   Last author: Michael Krause

            % Ensure the vectors are columns:
            qv1 = qv1(:);
            qv2 = qv2(:);

            qvr = qs1*qv2 + qs2*qv1 - cross(qv2, qv1);
            qsr = qs2*qs1 - qv2'*qv1;

            % Renormalize here:
            q_full = [qvr; qsr];
            q_full = q_full./norm(q_full);

            % And re-extract the components:
            qvr = q_full(1:3);
            qsr = q_full(4);

        end

    end

end