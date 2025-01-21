% This software is made available under the MIT license. See SONIC/LICENSE
% for details.
classdef Quadric

    properties (SetAccess = protected)
        % Quadric internally represented as the quadric locus
        locus         (4,4) double

        envelope      (4,4) double

        type                string

        proper              logical

        stable              logical
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
            %   Last revised: 11/14/24
            %   Last author:Michela Mancini
            arguments
                raw_quadric double
                inputType string
            end

            if inputType == "locus"
                % save locus and envelope
                obj.locus = sonic.Math.A_toDet1(raw_quadric);
                obj.proper = sonic.Quadric.isProperQuadric(obj.locus,"locus");
                obj.stable = sonic.Quadric.isNumericallyStable(obj.locus,"locus");
                obj.type = sonic.Quadric.getType(obj.locus,"locus");
                obj.envelope = sonic.Quadric.locusToEnvelope(obj.locus);

            elseif inputType == "envelope"
                % save locus and envelope
                obj.envelope =sonic.Math.A_toDet1(raw_quadric);
                obj.proper = sonic.Quadric.isProperQuadric(obj.envelope,"envelope");
                obj.stable = sonic.Quadric.isNumericallyStable(obj.envelope,"envelope");
                obj.locus = sonic.Quadric.envelopeToLocus(obj.envelope);
                obj.type = sonic.Quadric.getType(obj.envelope,"envelope");
            else
                error('sonic:Quadric:invalidInput', ...
                    ['Quadric representation must be specified as a string ' ...
                    'or char in the second input. ' ...
                    'Please specify either "locus" or "envelope". ']);
            end

        end

        function [centeredQuadric,center] = center(quadric)
            %% centeredQuadric = center(quadric)
            %
            %   Takes a quadric object and translates it by its center.
            %
            %   Inputs:
            %       - quadric (sonic.Quadric object): the original quadric
            %   Outputs:
            %       - centeredQuadric (sonic.Quadric object): the quadric
            %       translated by its center
            %       - center (sonic.Points3 object): new origin of the
            %       coordinate system (center of the original quadric)
            %
            %   Last revised: 11/14/2024
            %   Last author: Michela Mancini
            arguments
                quadric (1,1) sonic.Quadric
            end
            oldEnvelope = quadric.envelope;

            % find the center as the pole of the plane at infinity
            center = sonic.Points3(oldEnvelope*[0;0;0;1]);
            
            centeredQuadric = quadric.translate(center);
        end

        function [translatedQuadric] = translate(quadric,translation)
            %% translatedQuadric = translate(quadric,translation)
            %
            %   Takes a quadric object and a translation vector and outputs
            %   the quadric object shifted by the translation vector
            %
            %   Inputs:
            %       - quadric (sonic.Quadric object): the original quadric
            %       - translation (sonic.Points3 object): the translation
            %       vector (origin of the new coordinate system)
            %   Outputs:
            %       - translatedQuadric (sonic.Quadric object): the quadric
            %       translated by the translation vector
            %
            %   Last revised: 11/14/2024
            %   Last author: Michela Mancini
            arguments
                quadric (1,1) sonic.Quadric
                translation (1,1) sonic.Points3
            end
            quadricType = quadric.type;

            if strcmp(quadricType,"diskQuadric")
                oldEnvelope = quadric.envelope;
                translatedEnvelope = sonic.Quadric.translateQuadricMatrix(oldEnvelope,"envelope",translation);
                translatedQuadric = sonic.Quadric(translatedEnvelope,"envelope");
            else
                oldLocus = quadric.locus;
                translatedLocus = sonic.Quadric.translateQuadricMatrix(oldLocus,"locus",translation);
                translatedQuadric = sonic.Quadric(translatedLocus,"locus");
            end
        end   


        function [rotatedQuadric] = rotate(quadric,DCM)
            %% rotatedQuadric = rotate(quadric,translation)
            %
            %   Takes a quadric object and a translation vector and outputs
            %   the quadric object shifted by the translation vector
            %
            %   Inputs:
            %       - quadric (sonic.Quadric object): the original quadric
            %       - DCM (3x3 double): the rotation matrix
            %   Outputs:
            %       - rotatedQuadric (sonic.Quadric object): the quadric
            %       rotated by the DCM
            %
            %   Last revised: 11/14/2024
            %   Last author: Michela Mancini
            arguments
                quadric (1,1) sonic.Quadric
                DCM (3,3) double
            end
            quadricType = quadric.type;

            R = [DCM,[0;0;0];
                     0,0,0,1];
                
            if strcmp(quadricType,"diskQuadric")
                oldEnvelope = quadric.envelope;
                rotatedEnvelope = R*oldEnvelope*R';
                rotatedQuadric = sonic.Quadric(rotatedEnvelope,"envelope");
            else
                oldLocus = quadric.locus;
                RStar = sonic.Math.adjoint4x4(R);
                rotatedLocus = RStar'*oldLocus*RStar;
                rotatedQuadric = sonic.Quadric(rotatedLocus,"locus");
            end
        end   

    end


    methods (Static, Access=private)

        function [newMatrix, origin] = translateQuadricMatrix(matrix, matrixType,varargin)
            %% [newMatrix, origin] = translateQuadricMatrix(matrix, matrixType)
            %   Translates the quadric locus or quadric envelope by its
            %   center. If the center is not defined (two parallel planes)
            %   or if the center is at infinity, the function returns the
            %   original matrix
            %
            %   Inputs:
            %       - matrix (4x4 double): quadric locus or quadric envelope
            %       - matrixType (string or char): either quadric or envelope
            %   Outputs:
            %       - newMatrix (4x4 double): translated quadric locus or
            %       envelope
            %       - origin (1x1 sonic.Points3): origin of the new
            %       coordinate system
            %
            %   Last revised: 11/01/24
            %   Last author: Michela Mancini

            matrix = sonic.Math.A_toDet1(matrix);

            if nargin == 2
                if strcmp(matrixType,"locus")
                    A = matrix(1,1);
                    B = matrix(1,2)*2;
                    C = matrix(1,3)*2;
                    D = matrix(1,4)*2;
                    E = matrix(2,2);
                    F = matrix(2,3)*2;
                    G = matrix(2,4)*2;
                    H = matrix(3,3);
                    I = matrix(3,4)*2;

                    % set the center as new origin (center is the pole of the plane
                    % at infinity)
                    center = [D*F^2 - C*F*G + 2*B*G*H - B*F*I + 2*C*E*I - 4*D*E*H
                        C^2*G - C*D*F - B*C*I + 2*B*D*H - 4*A*G*H + 2*A*F*I
                        B^2*I - B*C*G - B*D*F + 2*C*D*E + 2*A*F*G - 4*A*E*I
                        - 2*H*B^2 + 2*B*C*F - 2*E*C^2 - 2*A*F^2 + 8*A*H*E];

                elseif strcmp(matrixType,"envelope")

                    % find the center as third column of the envelope
                    center = matrix*[0;0;0;1];
                end

                centerScale = log10(norm(center)/sonic.Tolerances.SmallNumber);

                if centerScale>-10 && centerScale<=0
                    % the point is not [0;0;0,0] so we can divide by biggest
                    % element
                    center = center/max(abs(center));
                elseif isinf(centerScale)
                    % if the point is [0;0;0;0] then the quadric is two parallel
                    % planes and we output the original matrix
                    newMatrix = matrix;
                    origin = sonic.Points3([0;0;0;1]);
                    return
                end

                origin = sonic.Points3(center);

                if origin.has_inf_points
                    % for the moment, this implementation is not translating
                    % objects with center at infinity
                    newMatrix = matrix;
                    origin = sonic.Points3([0;0;0;1]);
                    return
                end
            elseif nargin == 3
                origin = varargin{1};
            end


            if strcmp(matrixType,"locus")
                % adjoint of the transformation matrix that translates the
                % origin
                P_star = [origin.p3(4), 0,0, origin.p3(1);
                    0, origin.p3(4),0, origin.p3(2);
                    0, 0,origin.p3(4), origin.p3(3);
                    0,0,0,origin.p3(4)];

                newMatrix = P_star'*matrix*P_star;
            elseif strcmp(matrixType,"envelope")

                P = [origin.p3(4), 0,0, -origin.p3(1);
                    0, origin.p3(4),0, -origin.p3(2);
                    0, 0,origin.p3(4), -origin.p3(3);
                    0,0,0,origin.p3(4)];
                newMatrix = P*matrix*P';
            end
            newMatrix = sonic.Math.A_toDet1(newMatrix);


        end

        function [envelopeDetOne] = locusToEnvelope(locus)
            %% [envelope] = locusToEnvelope(locus)
            %   Converts locus representation to envelope
            %
            %   Inputs:
            %       - locus (4x4 double): quadric locus
            %   Outputs:
            %       - envelope (4x4 double): quadric envelope
            %
            %   Last revised: 7/11/24
            %   Last author: Michela Mancini

            % use adjoint to avoid failures with degenerate quadrics

            envelope = sonic.Math.adjoint4x4(locus);

            envelopeDetOne = sonic.Math.A_toDet1(envelope);

        end


        function [locusDetOne] = envelopeToLocus(envelope)
            %% [envelope] = envelopeToLocus(locus)
            %   Converts envelope representation to locus
            %
            %   Inputs:
            %       - envelope (4x4 double): quadric envelope
            %   Outputs:
            %       - locus (4x4 double): quadric locus
            %
            %   Last revised: 7/11/24
            %   Last author: Michela Mancini

            locus = sonic.Math.adjoint4x4(envelope);

            locusDetOne = sonic.Math.A_toDet1(locus);

        end


    end

    methods (Static)

        function stable = isNumericallyStable(matrix, matrixType)
            %% stable = isNumericallyStable(matrix)
            %   Determines whether the quadric matrix representation is
            %   numerically stable or not considering if an origin shift
            %   makes the previously degenerate matrix a non-degenerate one
            %
            %   Inputs:
            %       - matrix (4x4 double): quadric locus or envelope
            %   Outputs:
            %       - stable (1x1 logical): true if the quadric is proper, or
            %       degenerate and remains degenerate after origin shift,
            %       false if origin shift transforms degenerate quadric to
            %       non-degenerate quadric
            %
            %   Last revised: 10/25/24
            %   Last author: Michela Mancini

            stable = true;

            rcondNumber = rcond(matrix);

            proper = rcondNumber >= sonic.Tolerances.CondTol;

            if ~proper
                % perform change of coordinates to try improving numerical
                % behavior
                if strcmp(matrixType,"locus") || strcmp(matrixType,"envelope")
                    [new_matrix, ~] = sonic.Quadric.translateQuadricMatrix(matrix,matrixType);
                else
                    error('sonic:Quadric:isNumericallyStable:invalidInput', ...
                        ['Second input should be a string or char ' ...
                        'specifying if matrix is locus or envelope.']);
                end

                rcondNumber = rcond(new_matrix);

                if rcondNumber >= sonic.Tolerances.CondTol
                    stable = false;
                end
            end
        end

        function [proper] = isProperQuadric(matrix, matrixType)
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
            %   Last revised: 10/30/24
            %   Last author: Michela Mancini

            % checking condition number instead of determinant to overcome
            % scale issues
            rcondNumber = rcond(matrix);
            proper = rcondNumber >= sonic.Tolerances.CondTol;

            if ~proper
                [translatedMatrix,~] = sonic.Quadric.translateQuadricMatrix(matrix,matrixType);
                rcondNumber = rcond(translatedMatrix);
                proper = rcondNumber >= (10*sonic.Tolerances.CondTol);
            end

        end

        function[type] = getType(matrix, matrixType)
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

            if strcmp(matrixType,"locus")
                % get tolerances
                SphereTol = sonic.Tolerances.SphereTol;

                % check if the quadric is proper
                proper = sonic.Quadric.isProperQuadric(matrix, "locus");

                % stable = sonic.Quadric.isNumericallyStable(matrix,"locus");
                %
                % if ~stable
                [matrix,~] = sonic.Quadric.translateQuadricMatrix(matrix,"locus");
                % end

                % calculate determinant
                deter = det(matrix);

                % take upper left 3x3 minor
                locus44 = matrix(1:3,1:3);
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
                    r = rank(matrix);
                    if r==3
                        if r44==3
                            type = "cone";
                        elseif r44 == 2 || r44 ==1
                            type = "cylinder";
                        end
                    end
                end

            elseif strcmp(matrixType,"envelope")


                % check if the quadric is proper
                proper = sonic.Quadric.isProperQuadric(matrix, "envelope");

                % stable = sonic.Quadric.isNumericallyStable(matrix,"envelope");

                % if ~stable
                [matrix,~] = sonic.Quadric.translateQuadricMatrix(matrix,"envelope");
                % end

                if proper
                    locus = sonic.Quadric.envelopeToLocus(matrix);
                    type = sonic.Quadric.getType(locus,"locus");
                else
                    r = rank(matrix);
                    if r<3
                        error('sonic:Quadric:getType:invalidInput', ...
                            ['Only degenerate envelopes of rank 3 ' ...
                            '(disk quadrics) can be handled.'])
                    end
                    type = "diskQuadric";
                end
            else
                error('sonic:Quadric:getType:invalidInput', ...
                    ['Second input should be locus or ' ...
                    'envelope, as a char or a string variable.']);
            end
        end

    end

end