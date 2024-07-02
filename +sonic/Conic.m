% This software is made available under the MIT license. See SONIC/LICENSE
% for details.
classdef Conic

    properties (SetAccess = protected)
        % normalized implicit coefficients [A;B;C;D;E;F] satisfying the
        % equation Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
        implicit         (6,1) double
    end

    properties (SetAccess = private)
        % conic locus of the form [A, B/2, D/2; B/2, C, E/2; D/2, E/2, F]
        % satisfies all points lying on the conic
        locus            (3,3) double

        % inverse of conic locus satisfying all lines tangent to the conic
        envelope         (3,3) double

        % values for the explicit representation of a conic
        % [xc; yc; a; b; psi]
        explicit         (5,1) double = zeros(5,1);

        % conic classification (ellipse, circle, parabola, etc...)
        type            string

        % degeneracy of conic
        proper          logical

        % representation of degenerate conics as two lines
        degenerateLocus sonic.Lines2

    end

    methods

        function obj = Conic(raw_conic, varargin)
        %% obj = Conic(raw_conic, varargin)
        %   Creates a Conic object.
        %   
        %   Inputs:
        %       - raw_conic (1,6 double): Implicit representation of the
        %         conic. If inputs with 6 parameters, this is assumed to 
        %         be the coefficients of the implicit equation 
        %         [A;B;C;D;E;F] satisfying the equation 
        %         Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0.
        %
        %       OR
        %
        %       - raw_conic (1,5 double): Explicit representation of the
        %         conic. If inputs with 5 parameters, this is assumed to 
        %         be [xc; yc; a; b; psi] where xc and yc are the center, 
        %         a and b are the semi major and minor axis, and psi is 
        %         the conic rotation.
        %
        %       OR
        %
        %       - raw_conic (3,3 double): Matrix containing either the
        %         locus or envelope representation of the conic. Must 
        %         specify using varargin whether this is a locus or 
        %         envelope as shown below.
        %       - varargin (string): "locus" for the locus
        %         representation and "envelope" for the envelope
        %         representation.
        %
        %       OR
        %
        %       - raw_conic (1,1 sonic.Lines2): Lines2 object containing
        %         the degenerate conic representation of 2 lines. Must
        %         contain exactly two lines.
        %   Outputs:
        %       - obj (1,1 sonic.Conic): Conic object containing properties
        %         of the input conic.
        %
        %   Last revised: 4/19/24
        %   Last author: Michela Mancini
            if strcmp(class(raw_conic),"sonic.Lines2")
                
                if raw_conic.n ~= 2
                    error('sonic:Conic:invalidInput', ...
                                    ['Degenerate conic as ' ...
                                    'sonic.Lines2 object must contain',...
                                    'two lines']);
                 end
                
                obj.degenerateLocus = raw_conic;

                % convert to conic locus
                obj.locus = sonic.Conic.degenerateToLocus(obj.degenerateLocus);

                % calculate normalized implicit coefficients 
                obj.implicit = sonic.Conic.locusToImplicit(obj.locus);

                % calculate envelope
                obj.envelope = sonic.Conic.locusToEnvelope(obj.locus);

                % the input is a degenerate conic
                obj.proper = false;

                % get conic type
                obj.type = sonic.Conic.getType(obj.locus);
           
            else
                switch numel(raw_conic)

                    case 6
                        % convert to conic locus
                        obj.locus = sonic.Conic.implicitToLocus(raw_conic);

                        % calculate normalized implicit coefficients to match
                        % locus sign
                        obj.implicit = sonic.Conic.locusToImplicit(obj.locus);

                        % calculate envelope
                        obj.envelope = sonic.Conic.locusToEnvelope(obj.locus);

                        % check if conic is proper
                        obj.proper = sonic.Conic.isProperConic(obj.locus);

                        if obj.proper
                            % convert to explicit
                            obj.explicit = sonic.Conic.implicitToExplicit(obj.implicit);
                        else
                            % convert to two lines
                            obj.degenerateLocus = sonic.Conic.implicitToDegenerate(obj.implicit);
                        end

                        % get conic type
                        obj.type = sonic.Conic.getType(obj.locus);


                    case 5

                        % explicit [xc; yc; a; b; psi]
                        obj.explicit = raw_conic;

                        % convert to conic locus
                        obj.locus = sonic.Conic.explicitToLocus(obj.explicit);

                        % convert to implicit
                        obj.implicit = sonic.Conic.locusToImplicit(obj.locus);

                        % calculate envelope
                        obj.envelope = sonic.Conic.locusToEnvelope(obj.locus);

                        % the input is a non-degenerate conic
                        obj.proper = true;

                        % get conic type
                        obj.type = sonic.Conic.getType(obj.locus);

                    case 9
                        if isempty(varargin)
                            % if the locus vs envelope flag is not input, error
                            error('sonic:Conic:invalidInput', ...
                                ['Type of 3x3 matrix not defined. ' ...
                                'Must specify locus or envelope as a string or char '...
                                'variable for the second input.'])
                        else
                            if strcmp(varargin{1}, "locus")

                                % normalize locus to 1
                                obj.locus = sonic.Math.A_toDet1(raw_conic);

                                % calculate implicit coeffs
                                obj.implicit = sonic.Conic.locusToImplicit(obj.locus);

                                % calculate envelope
                                obj.envelope = sonic.Conic.locusToEnvelope(obj.locus);

                                % check conic degeneracy
                                obj.proper = sonic.Conic.isProperConic(obj.locus);

                                if obj.proper
                                    % convert to explicit
                                    obj.explicit = sonic.Conic.locusToExplicit(obj.locus);
                                else
                                    % convert to two lines
                                    obj.degenerateLocus = sonic.Conic.locusToDegenerate(obj.locus);
                                end

                                % get conic type
                                obj.type = sonic.Conic.getType(obj.locus);

                            elseif strcmp(varargin{1}, "envelope")
                                % normalize envelope
                                obj.envelope = sonic.Math.A_toDet1(raw_conic);

                                % check conic degeneracy
                                obj.proper = sonic.Conic.isProperConic(obj.envelope);

                                if ~obj.proper
                                    % if degenerate envelope, error
                                    error('sonic:Conic:invalidInput',...
                                        ['Degenerate conics are not accepted '...
                                        'in terms of the conic envelope.']);
                                else

                                    % calculate the locus
                                    obj.locus = sonic.Conic.envelopeToLocus(obj.envelope);

                                    % get the normalized implicit coefficients
                                    obj.implicit = sonic.Conic.locusToImplicit(obj.locus);

                                    % convert to explicit
                                    obj.explicit = sonic.Conic.implicitToExplicit(obj.implicit);

                                    % get conic type
                                    obj.type = sonic.Conic.getType(obj.locus);
                                end

                            else
                                % if the flag was input incorrectly, error
                                error('sonic:Conic:invalidInput', ...
                                    ['Type of 3x3 matrix uncertain. ' ...
                                    'Must specify locus or envelope as a string or char '...
                                    'variable for the second input.'])
                            end
                        end
                    otherwise
                        % if the input doesn't fit any of these types, error
                        error('sonic:Conic:invalidInput', ...
                            'Conic representation not recognized.');
                end
            end
        end
    end

    methods (Static)
    
        function [proper] = isProperConic(matrix)
        %% proper = isProperConic(matrix)
        %   Determines whether the conic locus represents a proper conic
        %   
        %   Inputs:
        %       - matrix (3x3 double): conic locus or conic envelope
        %   Outputs:
        %       - proper (1x1 logical): true if the conic is proper, false
        %         otherwise
        %
        %   Last revised: 4/19/24
        %   Last author: Michela Mancini

            % Checking condition number instead of determinant to overcome
            % scale issues
            tolCond = sonic.Tolerances.CondTol;
    
            rcondNumber = rcond(matrix);
    
            if rcondNumber >= tolCond
                proper = 1;
            else
                proper = 0;
            end
        end

        % given implicit --> convert to other forms
        function [locus_det1] = implicitToLocus(implicit)
        %% locus_det1 = implicitToLocus(implicit)
        %   Converts implicit representation to locus representation of
        %   conic
        %   
        %   Inputs:
        %       - locus_det1 (6x1 double): Vector containing implicit 
        %         coefficients [A;B;C;D;E;F]
        %   Outputs:
        %       - locus_det1 (3x3 double): Conic locus matrix
        %
        %   Last revised: 2/15/24
        %   Last author: Ava Thrasher
            A = implicit(1);
            B = implicit(2);
            C = implicit(3);
            D = implicit(4);
            E = implicit(5);
            F = implicit(6);

            locus = [A, B/2, D/2; B/2, C, E/2; D/2, E/2, F];
            locus_det1 = sonic.Math.A_toDet1(locus);
        end
        
        function [explicit] = implicitToExplicit(implicit)
        %% [explicit] = implicitToExplicit(implicit_normalized)
        %   Converts implicit representation to explicit
        %   
        %   Inputs:
        %       - implicit_normalized (6x1 double): Vector containing implicit 
        %         coefficients [A;B;C;D;E;F]
        %   Outputs:
        %       - explicit (5x1 double): Explicit representation of conic 
        %         [xc; yc; a; b; psi]
        %
        %   Last revised: 2/15/24
        %   Last author: Ava Thrasher
           
            locus = sonic.Conic.implicitToLocus(implicit);
            explicit = sonic.Conic.locusToExplicit(locus);
        end

        function [degenerateConic] = implicitToDegenerate(implicit)
        %% [degenerateConic] = implicitToDegenerate(implicit)
        %   Converts implicit representation to degenerate conic
        %   
        %   Inputs:
        %       - implicit_normalized (6x1 double): Vector containing implicit 
        %         coefficients [A;B;C;D;E;F]
        %   Outputs:
        %       - degenerateConic (1x1 Lines2 object): representation of
        %         the degenerate conic in terms of two lines
        %
        %   Last revised: 4/15/24
        %   Last author: Michela Mancini
            locus = sonic.Conic.implicitToLocus(implicit);
            degenerateConic = sonic.Conic.locusToDegenerate(locus);
        end

        % given conic locus --> convert to other forms
        function [implicit_normalized] = locusToImplicit(locus)
        %% [implicit_normalized] = locusToImplicit(locus)
        %   Converts locus representation to implicit
        %   
        %   Inputs:
        %       - locus (3x3 double): conic locus
        %   Outputs:
        %       - implicit_normalized (6x1 double): normalized implicit
        %         coefficients
        %
        %   Last revised: 2/15/24
        %   Last author: Ava Thrasher
            A = locus(1);
            B = locus(2)*2;
            C = locus(5);
            D = locus(3)*2;
            E = locus(6)*2;
            F = locus(9);
            implicit = [A; B; C; D; E; F];
            implicit_normalized = implicit./norm(implicit);
        end

        function [explicit] = locusToExplicit(locus)
        %% [explicit] = locusToExplicit(locus)
        %   Converts locus representation to explicit
        %   
        %   Inputs:
        %       - locus (3x3 double): conic locus
        %   Outputs:
        %       - explicit (5x1 double): explicit representation
        %         of the conic
        %
        %   Last revised: 4/19/24
        %   Last author: Michela Mancini

            psiConicTol = sonic.Tolerances.psiConicTol;

            % calculate determinant
            deter = det(locus);

            % using eig because we need to account for different signs for
            % hyperbola and parabola
            %[U,S,~] = svd(locus(1:2,1:2)); % for a symmetric matrix, singular value = |eigenvalue|            
            [U,S] = eig(locus(1:2,1:2));

            [S_sort,in_sort] = sort(diag(S));
            U_sort = U(:,in_sort);
           
            l1 = abs(S_sort(1));
            l2 = abs(S_sort(2));
            v2 = U_sort(:,2);

            A = locus(1);
            B = locus(2);
            C = locus(5);
            D = locus(3);
            E = locus(6);
            F = locus(9);
            disc = (B^2 - A*C);

            type = sonic.Conic.getType(locus);

            if ismember(type,["ellipse","circle","hyperbola"])
                % set semi-major and semi-minor axis radii (make sure that a is
                % always the larger one)
                a = sqrt(deter/(l1*l2^2)); % IF NORMALIZED TO 1, THIS IS COMPLEX AS WRITTEN IN THE PAPER (IF BOTH EIGVALS ARE POSITIVE)
                b = sqrt(deter/(l1^2*l2));

                xc = (C*D - B*E)/disc;
                yc = (A*E - B*D)/disc;

                % find psi (only calculate if a and b are NOT within tol of each other, otherwise, psi=0)
                if abs(a-b) < psiConicTol
                    psi = 0;
                else
                    psi = atan2(v2(2),v2(1));
                end

                if strcmp(type,"hyperbola")
                    a = -a;
                    b = -b;
                end

            elseif strcmp(type,"parabola")

                envelope = inv(locus);

                % determine the finite focus of the parabola as pole on the
                % line on the right
                focus = sonic.Points2(locus\[2*envelope(1,3);2*envelope(2,3);-envelope(1,1)-envelope(2,2)]);

                % determine the center as pole of the line at infinity
                center = sonic.Points2(locus \ [0;0;1]);

                % find the axis of the parabola joining the center and the
                % focus
                axisParabola = cross(focus.p2,center.p2); %sonic.GeometryP2.joinPoints(focus,center);

                y_q = -axisParabola(1)/axisParabola(2);
                y_p = -axisParabola(3)/axisParabola(2);

                coeff2 = (y_p*y_q*C+B*y_p+D+E*y_q);
                coeff3 = (C*y_p^2+2*E*y_p+F);

                % coordinates of the vertex of the parabola 
                xc = -coeff3/(2*coeff2);
                yc = y_p + xc*y_q;

                a = Inf;
                
                b = Inf;
              
                % we need to extract the eigenvector corresponding to the
                % zero eigenvalue, using this approach is shorter
                psi = atan2(center.p2(2),center.p2(1));

            else
                error('sonic:Conic:locusToExplicit:invalidInput',...
                    ['explicit parameterization only allowed for proper ',...
                    'conics.'])
            end

            % store
            explicit = [xc; yc; a; b; psi];
        end

        function [envelope] = locusToEnvelope(locus)
        %% [envelope] = locusToEnvelope(locus)
        %   Converts locus representation to envelope
        %   
        %   Inputs:
        %       - locus (3x3 double): conic locus
        %   Outputs:
        %       - envelope (3x3 double): conic envelope
        %
        %   Last revised: 4/15/24
        %   Last author: Michela Mancini

            % use adjoint to avoid failures with degenerate conics
            locus_det1 = sonic.Math.A_toDet1(locus);
    
            envelope = sonic.Math.adjoint3x3(locus_det1);
        
        end


        function [degenerateConic] = locusToDegenerate(locus)
        %% [degenerateConic] = locusToDegenerate(locus)
        %   Converts degenerate conic from locus to two lines
        %   
        %   Inputs:
        %       - locus (3x3 double): conic locus
        %   Outputs:
        %       - degenerateConic (1x1 Lines2 object): two lines
        %         representing the degenerate conic
        %
        %   Last revised: 4/15/24
        %   Last author: Michela Mancini            
         
            tolNum = sonic.Tolerances.SmallNumber;
    
            A = locus(1,1);
            B = locus(1,2);
            C = locus(2,2);
    
            % sample four points on the conic 
            line_coord = [0,1,1,5;
                          1,1,2,1;
                          1,1,1,1];
    
            y_q = -line_coord(1,:)./line_coord(2,:);
            coeff1 = (A+2*B*y_q+C*y_q.^2);
    
            [~,ind] = sort(abs(coeff1),'descend');
    
            line1 = sonic.Lines2(line_coord(:,ind(1)));
            line2 = sonic.Lines2(line_coord(:,ind(2)));
    
            points1 = sonic.Conic.meetLineConic(line1,locus);
            points2 = sonic.Conic.meetLineConic(line2,locus);
    
            pts = [points1.p2,points2.p2];
    
            r = rank(locus);
            if r>1
                % construct homography
                X = pts(:,1);
                Y = pts(:,2);
                Z = pts(:,3);
                E = pts(:,4);
    
                lambdas = [X Y Z]\E;
    
                H=[lambdas(1)*X,lambdas(2)*Y,lambdas(3)*Z];
    
                % transform conic matrix
                Cp = transpose(H)*locus*H;
    
                % extract lines in the new coordinates
                if abs(Cp(2,3))<tolNum
                    l1 = [1;0;0];
                    l2 = [0;Cp(1,2);Cp(1,3)];
                elseif abs(Cp(1,3))<tolNum
                    l1 = [0;1;0];
                    l2 = [Cp(1,2);0;Cp(2,3)];
                end
    
                % go back to original coordinates
                L1 = inv(transpose(H))*l1;
                L2 = inv(transpose(H))*l2;
                degenerateConic = sonic.Lines2([L1,L2]);
    
            else
    
                X = pts(:,1);
                Y = rand(3,1);
                Z = pts(:,3);
                E = rand(3,1);
    
                lambdas = [X Y Z]\E;
    
                H=[lambdas(1)*X,lambdas(2)*Y,lambdas(3)*Z];
    
                l12 = [0;1;0];
                L12 = inv(transpose(H))*l12;
                degenerateConic = sonic.Lines2([L12,L12]);
            end

        end

        % given explicit values --> convert to other forms
        function [locus] = explicitToLocus(explicit)
        %% locus = explicitToLocus(explicit)
        %   Converts explicit representation to conic locus
        %   
        %   Inputs:
        %       - explicit (5x1 double): explicit representation
        %         of the conic
        %   Outputs:
        %       - locus (3x3 double): conic locus
        %
        %   Last revised: 4/18/24
        %   Last author: Michela Mancini
            xc = explicit(1);
            yc = explicit(2);
            a = explicit(3);
            b = explicit(4);
            psi = explicit(5);

            if isinf(a)
                error('sonic:Conic:explicitToLocus:invalidInput', ...
                    ['The explicit parameterization of a parabola '....
                    'cannot be converted to other forms.']);     
            else
                if a>0
                    % for ellipse
                    A = a^2*sin(psi)^2 + b^2*cos(psi)^2;
                    B = 2*(b^2 - a^2)*cos(psi)*sin(psi);
                    C = a^2*cos(psi)^2 + b^2*sin(psi)^2;
                else
                    % for hyperbola
                    A = b^2*cos(psi)^2 - a^2*sin(psi)^2;
                    B= 2*(a^2+b^2)*cos(psi)*sin(psi);
                    C = b^2*sin(psi)^2 - a^2*cos(psi)^2;               
                end
                D = -2*A*xc - B*yc;
                E = -B*xc - 2*C*yc;
                F = A*xc^2 + B*xc*yc + C*yc^2 - a^2*b^2;
            end

            locus_unnorm = [A, B/2, D/2; B/2, C, E/2; D/2, E/2, F];
            locus = sonic.Math.A_toDet1(locus_unnorm);
        end

        function [implicit_normalized] = explicitToImplicit(explicit)
        %% implicit_normalized = explicitToImplicit(explicit)
        %   Converts explicit representation to normalized implicit
        %   coefficients
        %   
        %   Inputs:
        %       - explicit (5x1 double): explicit representation
        %         of the conic
        %   Outputs:
        %       - implicit_normalized (6x1 double): normalized implicit
        %         coefficients
        %
        %   Last revised: 2/15/24
        %   Last author: Ava Thrasher
            xc = explicit(1);
            yc = explicit(2);
            a = explicit(3);
            b = explicit(4);
            psi = explicit(5);

            A = a^2*sin(psi)^2 + b^2*cos(psi)^2;
            B = 2*(b^2 - a^2)*cos(psi)*sin(psi);
            C = a^2*cos(psi)^2 + b^2*sin(psi)^2;
            D = -2*A*xc - B*yc;
            E = -B*xc - 2*C*yc;
            F = A*xc^2 + B*xc*yc + C*yc^2 - a^2*b^2;

            implicit = [A; B; C; D; E; F];
            implicit_normalized = implicit./norm(implicit);
        end
        
        function [implicit_normalized]  = degenerateToImplicit(degenerateLocus)
        %% implicit_normalized = degenerateToImplicit(degenerateLocus)
        %   Converts the degenerate conic representation in terms of two
        %   lines to the implicit representation
        %   
        %   Inputs:
        %       - degenerateLocus (1x1 Lines2 object): lines representing
        %         the degenerate conic
        %    Outputs:
        %       - implicit_normalized (6x1 double): normalized implicit
        %         coefficients
        %
        %   Last revised: 4/15/24
        %   Last author: Michela Mancini
            locus = sonic.Conic.degenerateToLocus(degenerateLocus);
    
            A = locus(1,1); B = 2*locus(1,2); C = locus(2,2);
            D = 2*locus(1,3); E = 2*locus(2,3); F = locus(3,3);
    
            implicit = [A; B; C; D; E; F];
            implicit_normalized = implicit./norm(implicit);
        end

        function [locus] = degenerateToLocus(degenerateLocus)
        %% locus = degenerateToLocus(degenerateLocus)
        %   Converts the degenerate conic representation in terms of two
        %   lines to the conic locus
        %   
        %   Inputs:
        %       - degenerateLocus (1x1 Lines2 object): lines representing
        %         the degenerate conic
        %    Outputs:
        %       - locus (3x3 double): conic locus
        %
        %   Last revised: 4/15/24
        %   Last author: Michela Mancini

            [twoLines] = degenerateLocus.p2;
    
            L1 = twoLines(:,1);
            L2 = twoLines(:,2);
    
            locus = L2 * transpose(L1) + L1 * transpose(L2);
        end

        function [locus] = envelopeToLocus(envelope)
        %% [envelope] = envelopeToLocus(locus)
        %   Converts envelope representation to locus
        %   
        %   Inputs:
        %       - envelope (3x3 double): conic envelope
        %   Outputs:
        %       - locus (3x3 double): conic locus
        %
        %   Last revised: 4/16/24
        %   Last author: Michela Mancini

            envelope_det1 = sonic.Math.A_toDet1(envelope);
            
            locus = sonic.Math.adjoint3x3(envelope_det1);

        end

        function[type] = getType(locus)
        %% type = getType(locus)
        %   Determines type of conic
        %   
        %   Inputs:
        %       - locus (3x3 double): conic locus matrix
        %   Outputs:
        %       - type (1xn string): string signifying the type of the
        %         conic. Will be one of: "circle", "ellipse", "hyperbola",
        %         "parabola", "point", "intersecting_lines", "parallel_lines"
        %        
        %
        %   Last revised: 4/18/24
        %   Last author: Michela Mancini

            % get tolerances for checks
            discTol = sonic.Tolerances.ConicDisc; 
            circTol = sonic.Tolerances.ConicCircle;

            A = locus(1,1);
            B = locus(1,2)*2;
            C = locus(2,2);

            % calculate discriminant
            disc = B^2 - 4*A*C;

            if sonic.Conic.isProperConic(locus)
                if abs((A - C)^2 + B^2) < circTol
                    type = "circle";
                elseif disc <= -discTol  
                    type = "ellipse";
                elseif disc >= discTol 
                    type = "hyperbola";
                else % elseif abs(disc) < discTol
                    type = "parabola";
                end
                
            % added check for degeneracy
            else
                if disc <= -discTol
                    type = "point";
                elseif disc >= discTol
                    type = "intersecting_lines";
                else
                    type = "parallel_lines";
                end
            end            
        end


        function [p] = locusToParameter(locus)
        %% p = locusToParameter(locus)
        %   Determine the semi-latus rectum for a proper conic
        %   
        %   Inputs:
        %       - locus (3x3 double): conic locus matrix
        %   Outputs:
        %       - p (1x1 double): semi-latus rectum
        %
        %   Last revised: 4/19/24
        %   Last author: Michela Mancini
        
            if ~sonic.Conic.isProperConic(locus)
                error('sonic:Conic:locusToParameter:invalidInput', ...
                            'Conic should be proper.');
            end
            explicit = sonic.Conic.locusToExplicit(locus);
    
            a = explicit(3);
            b = explicit(4);
            if isinf(a)
                envelope = sonic.Conic.locusToEnvelope(locus);
    
                % determine focus coordinates
                focus = sonic.Points2(locus\[2*envelope(1,3);2*envelope(2,3);-envelope(1,1)-envelope(2,2)]);
    
                % determine vertex coordinates
                xc = explicit(1); yc = explicit(2);
    
                % distance between the focus and the vertex
                rmin = norm([xc;yc]-focus.p2(1:2)/focus.p2(3));
    
                % semi-latus rectum is twice the periapsis distance
                p = rmin*2;
    
            else
                if a>0 
                    %ellipse
                    p = b^2/a;
                else 
                    %hyperbola
                    p = -b^2/a;
                end
            end
        end

        function [locus] = parabolaToLocus(explicit_parabola)
        %% locus = parabolaToLocus(explicit_parabola)
        %   Determine the conic locus associated with a parabola with a
        %   specified vertex, tilt angle and semi-latus rectum
        %   
        %   Inputs:
        %       - explicit_parabola (4x1 double): contains a
        %         parameterization of the parabola in terms of vertex
        %         coordinates (xc,yc), tilt angle (psi) and the semi-latus
        %         rectum (p). The expected ordering is the following: 
        %         explicit_parabola = [xc; yc; psi; p]. 
        %   Outputs:
        %       - locus (3x3 double): conic locus
        %
        %   Last revised: 4/22/24
        %   Last author: Michela Mancini

            xc = explicit_parabola(1);
            yc = explicit_parabola(2);
            psi = explicit_parabola(3);
            p = explicit_parabola(4);
    
            A = sin(psi)^2;
            C = cos(psi)^2;
            B = -cos(psi)*sin(psi);
            D = -A*xc-B*yc-p*cos(psi);
            E = -C*yc-B*xc-p*sin(psi);
            F = A*xc^2+C*yc^2+2*B*xc*yc+2*p*(xc*cos(psi)+yc*sin(psi));
    
            locus = [A,B,D;
                     B,C,E;
                     D,E,F];
        end

        function pts = meetLineConic(line,locus)
        %% pts = meetLineConic(line,locus)
        %
        %   Performs the meet of a line and a conic, the result of
        %   which are two points in P2
        %
        %   Inputs:
        %       - line (a sonic.Line2 object): A line in P2.
        %       - locus (a 3x3 double array): A conic locus
        %   Outputs:
        %       - pts (1x1 sonic.Points2): two points of intersection of
        %         line and conic
        %
        %   Last revised: 04/22/2024
        %   Last author: Michela Mancini

            tolNum = sonic.Tolerances.SmallNumber;

            A = locus(1,1);
            B = locus(1,2);
            C = locus(2,2);
            D = locus(1,3);
            E = locus(2,3);
            F = locus(3,3);
 
            line_coord = line.p2;

            % choosing whether writing x as a function of y or viceversa
            % considering which line coordinate is the biggest (and can be
            % at the denominator)
            [m,ind] = max(abs(line_coord(1:2)));
 
            if m>=tolNum
                switch ind
                    case 1
                        % write x as a function of y from the line equation
                        x_q = -line_coord(2)/line_coord(1);
                        x_p = -line_coord(3)/line_coord(1);
    
                        % substitute the expression of x inside the conic equation:
                        % the remaining equation is quadratic in y
                        coeff1 = (C+2*B*x_q+A*x_q^2);
                        coeff2 = (x_p*x_q*A+B*x_p+E+D*x_q);
                        coeff3 = (A*x_p^2+2*D*x_p+F);
    
                        delta = coeff2^2-coeff1*coeff3;
    
                        y1 = (-coeff2+sqrt(delta))/coeff1;
                        x1 = x_p + y1 * x_q;
                        y2 = (-coeff2-sqrt(delta))/coeff1;
                        x2 = x_p + y2 * x_q;
    
                    case 2
                        % write y as a function of x from the line equation
                        y_q = -line_coord(1)/line_coord(2);
                        y_p = -line_coord(3)/line_coord(2);
    
                        % substitute the expression of x inside the conic equation:
                        % the remaining equation is quadratic in y
                        coeff1 = (A+2*B*y_q+C*y_q^2);
                        coeff2 = (y_p*y_q*C+B*y_p+D+E*y_q);
                        coeff3 = (C*y_p^2+2*E*y_p+F);
    
                        delta = coeff2^2-coeff1*coeff3;
    
                        x1 = (-coeff2+sqrt(delta))/coeff1;
                        y1 = y_p + x1 * y_q;
                        x2 = (-coeff2-sqrt(delta))/coeff1;
                        y2 = y_p + x2 * y_q;
                end
            else
                % the line is the line at infinity
                delta = B^2-A*C;
                if abs(A)>abs(C)                    
                    x1 = (-B+sqrt(delta))/A;
                    x2 = (-B-sqrt(delta))/A;
                    y1 = 1;
                    y2 = 1;
                else
                    y1 = (-B+sqrt(delta))/C;
                    y2 = (-B-sqrt(delta))/C;
                    x1 = 1;
                    x2 = 1;
                end
            end

            points = [x1,x2;y1,y2;1,1];

            pts = sonic.Points2(points);
        end

        function pts = meetConicConic(locus1,locus2)
        %% pts = meetConicConic(locus1,locus2)
        %
        %   Performs the meet of two conics, the result of which are
        %   four points in P2.
        %
        %   Inputs:
        %       - locus1 (a 3x3 double array): A non-degenerate conic locus
        %       - locus2 (a 3x3 double array): A non-degenerate conic locus
        %   Outputs:
        %       - pts (1x1 sonic.Points2): four points of intersection of
        %         the two conics
        %
        %   Last revised: 04/19/2024
        %   Last author: Michela Mancini

            tolConic = sonic.Tolerances.IsOnConic;
            
            if ~(sonic.Conic.isProperConic(locus1) && sonic.Conic.isProperConic(locus2))
                error('sonic:Conic:meetConicConic:degenerateConics', ...
                    ['The intersection of two conics requires the conics '...
                    'to be non-degenerate.']);
            end
    
            % determine vertices of the common self-polar trangle
            [vertices,~]=eig(locus2\locus1);
    
            X=vertices(:,1);
            Y=vertices(:,2);
            Z=vertices(:,3);
    
            dist_X_C1C2 = abs(transpose(X)*locus1*X)+abs(transpose(X)*locus2*X);
            dist_Y_C1C2 = abs(transpose(Y)*locus1*Y)+abs(transpose(Y)*locus2*Y);
            dist_Z_C1C2 = abs(transpose(Z)*locus1*Z)+abs(transpose(Z)*locus2*Z);
    
            justThreePoints=false;
    
            % checks if there is only one point of tangency
            if (dist_X_C1C2<tolConic && dist_Z_C1C2<tolConic) || (dist_X_C1C2<tolConic && dist_Y_C1C2<tolConic)
                Y=X;
                justThreePoints=true;
            elseif (dist_Y_C1C2<tolConic && dist_Z_C1C2<tolConic)
                Y=Z;
                justThreePoints=true;
            end
    
            if justThreePoints
                % select two points on the conic by intersecting the conic
                % with one line
    
                A = locus1(1,1);
                B = locus1(1,2);
                C = locus1(2,2);
    
                % choosing between three possibilities to avoid division by 
                % zero when intersecting the line and the conic
                lines = [0,1,1;
                         1,3,1;
                         1,1,1];
                y_q = -lines(1,:)./lines(2,:);
                coeff1 = (A+2*B*y_q+C*y_q.^2);
    
                [~,ind]=max(abs(coeff1));
    
                % intersect conics and best line
                twoPointsOnConic = sonic.Conic.meetLineConic(sonic.Lines2(lines(:,ind)),locus1);
    
                % define points for homography
                Z = twoPointsOnConic.p2(:,1);
                E = twoPointsOnConic.p2(:,2);
    
                % compute polar lines
                line1 = sonic.Lines2(locus1*Y);
                line2 = sonic.Lines2(locus1*Z);
    
                X = sonic.GeometryP2.meetLineLine(line1,line2).p2;
            else
    
                % choosing the unit point so that it is not aligned with two of
                % the vertices of the self-polar triangle.
                Evec = [1,1,1,1;
                        1,2,5,7;
                        1,1,1,1];
    
                xdir = X/norm(X); 
                ydir = Y/norm(Y);
                zdir = Z/norm(Z);
    
                check = zeros(4,1);
                for i = 1 : 4
                    edir = Evec(:,i)/norm(Evec(:,i));
                    check1 = rcond([xdir,ydir,edir]);
                    check2 = rcond([xdir,zdir,edir]);
                    check3 = rcond([ydir,zdir,edir]);
                    check(i) = min([check1,check2,check3]);
                end
                [~,bestEind]=max(check);
                E = Evec(:,bestEind);
            end
           
            % determine scaling elements for homography
            lambdas = [X Y Z]\E;
    
            % 3x3 array, homography matrix
            H=[lambdas(1)*X,lambdas(2)*Y,lambdas(3)*Z];
    
            % conic locus in the new frame
            locus2p = transpose(H)*locus2*H;
            locus1p = transpose(H)*locus1*H;
    
            if justThreePoints
                % solve the intersection between a parabola and the
                % transformed conic
                a2 = locus2p(1,1); d2 = locus2p(1,3); e2 = locus2p(2,3); f2=locus2p(3,3);
    
                x1p = (-d2+sqrt(d2^2-f2*(a2+2*e2)))/(a2+2*e2);
                x2p = (-d2-sqrt(d2^2-f2*(a2+2*e2)))/(a2+2*e2);
                y1p = x1p^2;
                y2p = x2p^2;
    
                % intersection points in the new frame
                pointsNewFrame = [x1p,x2p;y1p,y2p;1,1];
    
                pointsOldFrame=[H*pointsNewFrame,Y,Y];
    
            else
                % solve intersection of conics in the new frame, which both
                % have a diagonal locus matrix
    
                a1 = locus1p(1,1); b1 = locus1p(2,2); c1 = locus1p(3,3);
                a2 = locus2p(1,1); b2 = locus2p(2,2); c2 = locus2p(3,3);
    
                % points of intersection in the new coordinates
                Ysq = (c1*a2/a1-c2)/(-a2/a1*b1+b2);
                Xsq = (-b1*Ysq-c1)/a1;
                x1p = sqrt(Xsq);
                y1p =sqrt(Ysq);
    
                pointsNewFrame = [x1p, x1p, -x1p,-x1p;
                                  y1p, -y1p, y1p,-y1p;
                                  ones(1,4)];
    
                pointsOldFrame = H*pointsNewFrame;
    
            end
    
            pts = sonic.Points2(pointsOldFrame);

        end

    end
end