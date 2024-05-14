classdef Quadric

    properties (SetAccess = protected)
        % Quadric internally represented as the quadric locus
        locus         (4,4) double

        envelope      (4,4) double

        type                string

        proper              logical
    end

    methods

        function obj = Quadric(raw_quadric, inputType)
            %% obj = Quadric(raw_quadric, type)
            %
            %   Instantiates a Quadric object.
            %
            %   Inputs:
            %       - raw_quadric (4x4 double): a 4x4 matrix representing the
            %         quadric
            %       - type (string or char): specifying what type of input.
            %         "locus", "envelope"
            %
            %   Outputs:
            %       - obj (1,1 sonic.Quadric): Quadric object containing locus
            %         and envelope
            %
            %   Last revised: 05/07/24
            %   Last author:Michela Mancini
            arguments
                raw_quadric double
                inputType string
            end

            if inputType == "locus"
                % save locus and envelope
                obj.locus = raw_quadric;
                obj.proper = sonic.Quadric.isProperQuadric(obj.locus);
                obj.type = sonic.Quadric.getType(obj.locus);
                if obj.proper
                    obj.envelope = inv(obj.locus);
                end


            elseif inputType == "envelope"
                % save locus and envelope
                obj.envelope = raw_quadric;
                obj.proper = sonic.Quadric.isProperQuadric(obj.envelope);
                if obj.proper
                    obj.locus = inv(obj.envelope);
                    obj.type = sonic.Quadric.getType(obj.locus);
                else
                    obj.type = "diskQuadric";
                end

            else
                error('sonic:Quadric:invalidInput', ...
                    ['Quadric representation must be specified as a string or char in the second input. ' ...
                    'Please specify either "locus" or "envelope". ']);
            end

        end
    end


    methods (Static)

        function [proper] = isProperQuadric(matrix)
            %% proper = isProperQuadric(matrix)
            %   Determines whether the quadric locus represents a proper
            %   quadric
            %
            %   Inputs:
            %       - matrix (4x4 double): quadric locus or quadric envelope
            %   Outputs:
            %       - proper (1x1 logical): true if the quadric is proper,
            %         false otherwise
            %
            %   Last revised: 4/24/24
            %   Last author: Michela Mancini

            % checking condition number instead of determinant to overcome
            % scale issues
            tolCond = sonic.Tolerances.CondTol;

            rcondNumber = rcond(matrix);

            if rcondNumber >= tolCond
                proper = 1;
            else
                proper = 0;
            end
        end

        function[type] = getType(locus)
            %% type = getType(locus)
            %   Determines type of quadric
            %
            %   Inputs:
            %       - locus (4x4 double): quadric locus matrix
            %   Outputs:
            %       - type (1xn string)
            %
            %   Last revised: 4/24/24
            %   Last author: Michela Mancini

            % get tolerances
            SphereTol = sonic.Tolerances.SphereTol;

            % calculate determinant
            deter = det(locus);

            % check if the quadric is proper
            proper = sonic.Quadric.isProperQuadric(locus);

            % take upper left 3x3 minor
            locus44 = locus(1:3,1:3);
            r44 = rank(locus44);
            eig44 = eig(locus44);
            eig44 = sort(eig44);

            if proper
                if r44==3
                    if sign(eig44(1))==sign(eig44(2)) && sign(eig44(1))==sign(eig44(3))
                        % eigenvalues are of the same sign
                        
                        diffEig = abs([diff(eig44);eig44(3)-eig44(1)]);
                        
                        if diffEig(3)<SphereTol 
                            % if difference between the biggest and
                            % smallest is small, the three eigenvalues are
                            % all the same
                            type = "sphere";
                        elseif diffEig(1)<SphereTol
                            % first and second eigenvalue are the same
                            if abs(eig44(3))<abs(eig44(1))
                                type = "prolateSpheroid";
                            else
                                type = "oblateSpheroid";
                            end
                        elseif diffEig(2)<SphereTol
                            % second and third eigenvalue are the same
                            if abs(eig44(1))<abs(eig44(2))
                                type = "prolateSpheroid";
                            else
                                type = "oblateSpheroid";
                            end
                        else
                            % all eigenvalues are different
                            type = "ellipsoid";
                        end
                    else
                        % eigenvalues have different sign
                        type = "hyperboloid";
                    end
                elseif r44==2
                    if deter>0
                        type = "hyperbolicParaboloid";
                    else
                        type ="ellipticParaboloid";
                    end
                end
            else
                r = rank(locus);
                if r==3
                    if r44==3
                        type = "cone";
                    elseif r44 == 2 || r44 ==1
                        type = "cylinder";
                    end
                end
            end

        end


    end

end