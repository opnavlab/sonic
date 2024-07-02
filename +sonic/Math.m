% This software is made available under the MIT license. See SONIC/LICENSE
% for details.
classdef Math

    methods (Static)

        function adj = adjoint3x3(matrix)
        %% adj = adjoint3x3(matrix)
        %   Determines the adjoint of the given 3x3 matrix
        %
        %   Inputs:
        %       - matrix (3x3 double): any 3x3 matrix
        %   Outputs:
        %       - adj (3x3 double): adjoint of the given matrix
        %
        %   Last revised: 4/19/24
        %   Last author: Michela Mancini

        A = matrix(1,1);
        B = matrix(1,2);
        C = matrix(2,2);
        D = matrix(1,3);
        E = matrix(2,3);
        F = matrix(3,3);

        adj = [C*F - E^2, D*E - B*F, B*E - C*D;
            D*E - B*F, A*F - D^2, B*D - A*E;
            B*E - C*D, B*D - A*E, A*C - B^2];
        end
    
        function T = crossmat(v)
        %% T = crossmat(v)
        %   3x3 skew-symmetric cross-product matrix
        %
        %   Inputs:
        %       - v (3 double): vector to convert to cross-product matrix
        %   Outputs:
        %       - T (3x3 double): skew-symmetric matrix, such that multiplying this
        %         matrix with another 3x1 vector is equivalent to the 
        %         cross-product of those two vectors.
        %
        %   Last revised: 2/14/24
        %   Last author: Michael Krause
        
            if numel(v) ~= 3
                error('sonic:crossmat:invalidSize', ...
                    'Must input a vector of length 3.');
            end
        
            T = [0, -v(3), v(2);
                 v(3), 0, -v(1);
                 -v(2), v(1), 0];
        
        end

        function val = atan2c(y, x)
        %% val = atan2c(y, x)
        %
        %   Returns the 4-quadrant arctangent for complex values of y and
        %   x, in radians.
        %
        %   Inputs:
        %       - y (1xn double): Numerator of the tangent (Y/X). May be
        %         complex.
        %       - x (1xn double): Denominator of the tangent (Y/X). May be
        %         complex.
        %
        %   Outputs:
        %       - val (1xn double): The 4-quadrant arctangent of Y/X, in 
        %         radians. If Y and X are real, this is equivalent to 
        %         atan2(). 
        %
        %   Last revised: 5/07/24
        %   Last author: Michael Krause  

            arguments
                y   (1, :)  double
                x   (1, :)  double
            end

            % Ensure y and x are the same length:
            if length(y) ~= length(x)
                error('sonic:Math:atan2c:dimensionMismatch', ...
                    'Dimensions of `y` and `x` do not match.');
            end

            % Note signs for later:
            sy = sign(real(y));
            sx = sign(real(x));

            % Calculate the "raw" atan value using complex logs:
            z = y./x;
            raw_atan = (1/(2i)).*log((1i - z)./(1i + z));

            % Clamp real and imaginary parts if close to machine precision:
            realpart = real(raw_atan);
            imagpart = imag(raw_atan);

            zero_real = abs(realpart) < sonic.Tolerances.CompZero;
            realpart(zero_real) = 0;

            zero_imag = abs(imagpart) < sonic.Tolerances.CompZero;
            imagpart(zero_imag) = 0;
            
            % Assemble the final value:
            val = realpart + 1i.*imagpart;

            for idx = 1:length(val)
                % Handle the four-quadrant cases:
                switch sx(idx)
                    case 1
                        % Do nothing.
                    case 0
                        switch sy(idx)
                            case 1
                                val(idx) = pi/2;
                            case 0
                                % Technically undefined, but set to 0 to match 
                                % Matlab's conventions. 
                                val(idx) = 0;        
                            case -1
                                val(idx) = -pi/2;
                        end
                    case -1
                        switch sy(idx)
                            case 1
                                val(idx) = val(idx) + pi;
                            case 0
                                val(idx) = val(idx) + pi;
                            case -1
                                val(idx) = val(idx) - pi;
                        end
                end
            end

        end

        function A_det1 = A_toDet1(A)
        %% A_det1 = A_toDet1(A)
        %   Normalizes matrix to determinant of 1 if the matrix is
        %   non-singular, otherwise returns the original function
        %   
        %   Inputs:
        %       - A (nxn double): matrix
        %   Outputs:
        %       - A_det1 (nxn string): matrix with determinant of 1
        %
        %   Last revised: 4/18/24
        %   Last author: Michela Mancini
            Adim = size(A);
            if Adim(1) ~= Adim(2)
                error('A_toDet1:NonSquareMtx', 'Matrix A must be square.');
            end
            
            N = Adim(1);
            detA = det(A);
            % difference tolerance not with conic
            detTol = sonic.Tolerances.ConicLocusDetOne;
            % added check for degeneracy
            if abs(detA) >= detTol
                A_det1 = A./nthroot(det(A), N);
            else
                A_det1 = A;
            end
        end

    end

end