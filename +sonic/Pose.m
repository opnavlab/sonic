% This software is made available under the MIT license. See SONIC/LICENSE
% for details.
classdef Pose

    properties
        att                 sonic.Attitude
        t                   sonic.Points3
        T                   (4, 4) double % transformation matrix
    end

    methods
        function obj = Pose(arg1, arg2)
            %% obj = pose(arg1, arg2)
            %   Constructs a Pose object.
            %   A pose contains both an Attitude and a Points3
            %   The pose from A to B is such that:
            %   x_B = T * x_A = [pose.att.dcm, pose.t] * x_A
            %
            %   Inputs:
            %       - arg1
            %         If only one argument
            %           - T (4x4 double): A matrix containing the DCM and
            %           translation such that x_B.p3 = T * x_A.p3
            %         If two arguments
            %           - att (1x1 sonic.Attitude): the attitude of the pose
            %       - arg2 (sonic.Points3): the translation of the point in B
            %
            %   Outputs:
            %       - obj (sonic.Pose): Pose object.
            %
            %   Last revised: 10/31/24
            %   Last author: Sebastien Henry


            if nargin == 1
                [m, n] = size(arg1);
                if ~(m==4 && n==4)
                    error('sonic:Pose:invalidSize', ['the transformation matrix T', ...
                        'should be 4x4']);
                end

                if abs(arg1(4,4) - 1) > sonic.Tolerances.CompZero
                    error('sonic:Pose:invalidSE3', ['the transformation matrix T' ...
                        'should have its last entry as 1']);
                end

                if any( abs(arg1(4,1:3)) > sonic.Tolerances.CompZero)
                    error('sonic:Pose:invalidSE3', ['the transformation matrix T' ...
                        'should have its last row as [0 0 0 1]']);
                end


                obj.T = arg1;

                obj.att = sonic.Attitude(obj.T(1:3, 1:3));
                obj.t = sonic.Points3(obj.T(:, 4));

            elseif nargin == 2
                if ~isa(arg1,'sonic.Attitude')
                    error('sonic:Pose:invalidInputType', ['the attitude ' ...
                        'should be a sonic.Attitude']);
                end
                obj.att = arg1;

                % make sure that the translation is a Points3
                if ~isa(arg2,'sonic.Points3')
                    error('sonic:Pose:invalidInputType', ['the translation ' ...
                        'should be a sonic.Points3']);
                end
                if ~(arg2.n ==1)
                    error('sonic:Pose:invalidPoints3', ['the translation ' ...
                        'Points3 should contain exactly one point']);
                elseif (arg2.has_inf_points)
                    error('sonic:Pose:invalidPoints3', ['currently ' ...
                        'does not support points at infinity']);
                end

                obj.t = arg2;
                obj.T = [obj.att.dcm, obj.t.r3;
                    0 0 0 1];
            end

        end

        function newPose = comp(obj, obj2)
            %% [rot_obj] = transform(obj, to_transform)
            %
            %   Transforms a point such that
            %   x_B = P * x_A = [pose.att.dcm, pose.t] * x_A
            %
            %   Inputs:
            %       - to_transform (1x1 sonic.Points3): Object to
            %         transform. Can have more than n=1.
            %
            %   Outputs:
            %       - transformed_obj (1x1 sonic.Points3 or sonic.PointsS2):
            %           Rotated object. Return type will match the input type.
            %
            %   Last revised: 10/31/24
            %   Last author: Sebastien Henry

            newPose = sonic.Pose(obj.T*obj2.T);

        end

        function newPose = inv(obj)
            %% [rot_obj] = inv(obj)
            %
            %   Inverse a pose
            %
            %   Outputs:
            %       - newPose (1x1 sonic.Pose): Inversed Pose
            %
            %   Last revised: 10/31/24
            %   Last author: Sebastien Henry

            att_inv = obj.att.inv();
            newPose = sonic.Pose(att_inv, ...
                sonic.Points3(-att_inv.rotate(obj.t).r3));

        end


        function transformed_obj = transform(obj, to_transform)
            %% transformed_obj = transform(obj, to_transform)
            %
            %   Transforms a point such that
            %   x_B = T * x_A = [pose.att.dcm, pose.t] * x_A
            %
            %   Inputs:
            %       - to_transform (1x1 sonic.Points3): Object to
            %         transform. Can have more than n=1.
            %
            %   Outputs:
            %       - transformed_obj (1x1 sonic.Points3 or sonic.PointsS2): Rotated
            %         object. Return type will match the input type.
            %
            %   Last revised: 10/31/24
            %   Last author: Sebastien Henry

            transformed_obj = sonic.Points3(obj.T * to_transform.p3);
        end

        function multiplied_obj = mtimes(obj, to_mult)
            %% multiplied_obj = mtimes(obj, to_mult)
            % Pre-multiply the pose to an object
            % result depends on the type of to_mult
            %
            %   Inputs:
            %       - to_mult: a SONIC object, either
            %           - (1x1 sonic.Points3): In this case
            %           the output the pose matrix times the points p3.
            %           - (1x1 sonic.Pose): In this case, the output is the
            %           compose of the poses.
            %
            %   Outputs:
            %       - multiplied (1x1 sonic.Points3 or sonic.Pose): object
            %       pre-multiplied by the pose.
            %   Last revised: 1/07/24
            %   Last author: Sebastien Henry

            switch class(to_mult)
                case 'sonic.Points3'
                    multiplied_obj = obj.transform(to_mult);
                case 'sonic.Pose'
                    multiplied_obj = obj.comp(to_mult);
                otherwise
                    error('sonic:Pose:mtimes:invalidInputType', ['Only' ...
                        'Points3 and Pose can be multiplied by a Pose'])
            end
        end

    end
end