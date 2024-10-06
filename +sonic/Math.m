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

        function [minSingVal, minRSingVect] = getSmallestRSingVector(A)
        %% [minSingVal, minRSingVect] = getSmallestRSingVector(A)
        %
        %   Get right singular vector corresponding to smallest magnitude
        %   magnitude value
        %
        %   Inputs:
        %       - A (mxn double): matrix for which we want to extract 
        %         smallest singular value 
        %
        %   Outputs:
        %       - minSingVal (1x1 double) -- smallest magnitude singular
        %            value
        %       - minRSingVect (nx1 double) -- right singular vector 
        %         corresponding to minSingVal
        %
        %   Last revised: Sep 24, 2024
        %   Last author: Tara Mina
            arguments
                A                 (:, :)     double     {mustBeNonempty}
            end

            % get matrix size for analysis
            [nRow,nCol] = size(A);

            % get singular values and right singular vectors
            if nRow >= nCol
                % if nRow > nCol, no need to compute all nRow singular 
                % values of U (since we're not using it)
                [~,S,V] = svd(A, "econ");
            else
                error('This function requires number of rows >= number of columns.')
            end
            singVals = diag(S);

            % if at least 2 columns in V, check no repeated smallest 
            % singular values -- if so, throw an error
            if (nCol>=2 && abs(singVals(end-1) - singVals(end)) < ...
                sonic.Tolerances.CompZero)
                error(['Repeated smallest singular values in ', ...
                    'the input matrix, A.'])
            end

            % get singular value corresponding to column index
            minRSingVect = V(:,end);
            minSingVal = singVals(end);

        end
        
        function minEigvec = solveGenEigvalSmallestLambda(M,N,boolNorm)
        %% minEigvec = solveGenEigvalSmallestLambda(M,N,boolNorm)
        %
        %   Helper function to solve the generalized eigenvalue problem,
        %   returning the eigenvector corresponding to the smallest
        %   magnitude eigenvalue \lambda
        %
        %   Generalized eigenvalue problem finds x and \lambda such that: 
        %                   M x = \lambda N x
        %   And this function returns the eigenvector x corresponding to 
        %   the eigen value \lambda with the smallest magnitude
        %
        %   Inputs:
        %       - M (nxn double): square matrix which corresponds to the 
        %         generalized eigenvalue problem, defined above
        %       - N (nxn double): square matrix normalization factor which
        %         corresp to generalized eigenvalue problem, defined above
        %       - boolNorm (boolean): if True, L2-normalizes the
        %         eigenvector
        %
        %   Outputs:
        %       - minEigvec (nx1 double) -- solution to the above-defined
        %         generalized eigenvalue problem, corresponding to the
        %         eigenvalue \lambda with the smallest magnitude
        %
        %   Last revised: Sep 24, 2024
        %   Last author: Tara Mina

            arguments
                M                    (:, :)      double
                N                    (:, :)      double
                boolNorm            (1, 1)      logical = true
            end
            assert( size(M,1)==size(M,2), 'Matrix M must be square');
            assert( size(N,1)==size(N,2), 'Matrix N must be square');
            assert( size(M,1)==size(N,1), ...
                'Matrix M and N must be the same dimensions');

            % get generalized eigenvalues and eigenvectors
            [Vg, Dg] = eig(M,N); 

            % get absolute eigenvalues, find smallest magnitude via sorting
            absEigvals = diag(abs(Dg));
            [absEigvalSorted, eigvalIdcSorted] = sort(absEigvals, 'descend');

            % from sorted vector, check if more than 2 columns, and if so,
            % check if minimum generalized eigenvalue is repeated
            if length(absEigvalSorted)>=2 && ...
                    abs(absEigvalSorted(end-1)-absEigvalSorted(end)) ...
                    < sonic.Tolerances.CompZero
                error(['Repeated minimum generalized eigenvalues ', ...
                    'in the input matrix pair (M,N)']);
            end

            % get minimum generalized eigval idx & corresponding eigvector
            minEigvalColi = eigvalIdcSorted(end);
            minEigvec = Vg(:,minEigvalColi);
        
            % normalize eigenvector in L2, if flag is true
            if boolNorm
                minEigvec = minEigvec / norm(minEigvec);
            end
        end

    function [val] = HFunction(mu,AL_ss,evalMethod)
        %% Hval = HFunction(mu,AL_ss,evalMethod)
        %   Computes the Chandrasekhar's H-function for a given set of
        %   directional parameters (mu) an single-scattering albedo 
        %   (AL_ss) using a user-defined evaluation method
        %
        %   Inputs:
        %       - mu (nxm double): matrix of values on the domain 
        %         [0,1] (can be a scalar as well)
        %       - AL_ss (nxm double): matrix of values on the domain 
        %         [0,1] (can be a scalar as well)
        %       - evalMethod (string): defines evaluation method of 
        %         H-Function
        %       (must be "exact", "numeric", or "linear")
        %       - Note: mu and AL_ss must have the same dimensions unless mu and/or
        %       AL_ss is a scalar value
        %   Outputs:
        %       - Hval (nxm double): matrix of H-function values for either the
        %       pairwise (if both mu and AL_ss are matrices) or for each mu/AL_ss given
        %       a scalar AL_ss/mu
        %
        %   References:
        %       [1] D. W. N. Stibbs, and R. E. Weir, 
        %           "On the H-Functions for Isotropic Scattering", 
        %           https://doi.org/10.1093/mnras/119.5.512
        %       [2] B. Hapke, "Bidirectional Reflectance Spectroscopy: 
        %           1. Theory", https://doi.org/10.1029/JB086iB04p03039
        %       [3] J. A. Christian, "Spacecraft Optical Navigation"
        %
        %   Last revised: 9/24/24
        %   Last author: Jennifer Nolan
        %
        arguments
            mu         (:,:) double {mustBeInRange(mu,0,1)}
            AL_ss      (:,:) double {mustBeInRange(AL_ss,0,1)}
            evalMethod (1,1) string {mustBeMember(evalMethod,{'exact','rational','linear'})}
        end


        % Verify Input Dimensions
        % Check if mu and AL_ss have different dimensions
        if any(size(mu)~=size(AL_ss))

            % Repeat mu input for all AL_ss values if mu is scalar
            if size(mu,1)*size(mu,2)==1
                mu = mu.*ones(size(AL_ss));

            % Repeat AL_ss input for all mu values if AL_ss is scalar
            elseif size(AL_ss,1)*size(AL_ss,2)==1
                AL_ss = AL_ss.*ones(size(mu));

            % If neither mu nor AL_ss are scalar, error on grounds of
            % non-compatible dimensions
            else
                error('sonic:HFunction:IncompatibleMatrices', 'Either the dimensions of the mu and AL_ss must be identical or at least one input must be a scalar.');
            end
        end

        % Specify NaN values and only continue with computation on
        % non-NaN doubles
        nanInds = isnan(mu) | isnan(AL_ss);
        mu = mu(~nanInds);
        AL_ss = AL_ss(~nanInds);

        % Call Specified H-Function Calculation Method:
        % Execute the "exact" evaluation of H-Function, following the
        % Stibbs-Weir method described in [1]
        if evalMethod == "exact"
            Hval = sonic.Math.exactHFunction(mu, AL_ss);

        % Execute the "rational" evaluation of H-Function, following the
        % rational apprximation described in [2]
        elseif evalMethod == "rational"
            Hval = sonic.Math.rationalHFunction(mu,AL_ss);
        
        % Execute the "linear" evaluation of H-Function, following the
        % linear approximation described in [3]
        elseif evalMethod == "linear"
            Hval = sonic.Math.linearHFunction(mu,AL_ss);

        end

        % Construct output where computed values are placed at their
        % original indices and NaNs avoid computation and are
        % reinserted here
        val = zeros(size(nanInds));
        val(nanInds)=NaN;
        val(~nanInds)=Hval;
    end

    function [Hval] = linearHFunction(mu,AL_ss)
        %% Hval=linearHFunction(mu,AL_ss)
        %   Uses the linearized method [1] to compute the 
        %   Chandrasekhar's H-function for a given set of directional 
        %   parameters (mu) and single-scattering albedos (AL_ss).
        %
        %   Inputs:
        %       - mu (nxm double): matrix of values on the domain 
        %         [0,1] (can be a scalar as well)
        %       - AL_ss (nxm double): matrix of values on the domain 
        %         [0,1] (can be a scalar as well)
        %       - Note: mu and AL_ss mustst have the same dimensions 
        %         unless mu and/or AL_ss is a scalar value
        %   Outputs:
        %       - val (nxm double): matrix of H-function values for 
        %         either the pairwise (if both mu and AL_ss are 
        %         matrices) or for each mu/AL_ss given a scalar 
        %         AL_ss/mu
        %
        %   References:
        %       [1] J. A. Christian, "Spacecraft Optical Navigation"
        %
        %   Last revised: 9/30/24
        %   Last author: Jennifer Nolan
        %
        arguments
            mu       (:,:) double {mustBeInRange(mu,0,1)}
            AL_ss      (:,:) double {mustBeInRange(AL_ss,0,1)}
        end

        % Define the values of gamma and epsilon, parameters dependent 
        % only on single-scattering albedo
        gamma = sqrt(1-AL_ss);
        epsilon = (1-gamma)./(1+gamma);
        
        % Evaluate the linear H-function for the given set of 
        % (i,e,AL_ss) inputs using Eq. 6.129 from [1]
        Hval=(1+epsilon).*(1-(epsilon./2)+epsilon.*mu);

    end

    function [H] = rationalHFunction(mu, AL_ss)
        %% Hval=rationalHFunction(mu,AL_ss)
        %   Uses the rational approximation method [1] to compute the 
        %   Chandrasekhar's H-function for a given set of directional 
        %   parameters (mu) and single-scattering albedos (AL_ss).
        %
        %   Inputs:
        %       - mu (nxm double): matrix of values on the domain 
        %         [0,1] (can be a scalar as well)
        %       - AL_ss (nxm double): matrix of values on the domain 
        %         [0,1] (can be a scalar as well)
        %       - Note: mu and AL_ss must have the same dimensions 
        %         unless mu and/or AL_ss is a scalar value
        %   Outputs:
        %       - val (nxm double): matrix of H-function values for 
        %         either the pairwise (if both mu and AL_ss are 
        %         matrices) or for each mu/AL_ss given a scalar 
        %         AL_ss/mu
        %
        %   References:
        %       [1] B. Hapke, "Bidirectional Reflectance Spectroscopy: 
        %           1. Theory", https://doi.org/10.1029/JB086iB04p03039
        %       [2] J. A. Christian, "Spacecraft Optical Navigation"
        %
        %   Last revised: 9/30/24
        %   Last author: Jennifer Nolan
        %
        arguments
            mu         (:,:) double {mustBeInRange(mu,0,1)}
            AL_ss      (:,:) double {mustBeInRange(AL_ss,0,1)}
        end
    
        % Define the value gamma, a parameter dependent only of
        % single-scattering albedo using Eq. 8 in [1]
        gamma = sqrt(1-AL_ss);

        % Compute the H-function value directly for the given pairs of 
        % (mu,AL_ss) inputs using Eq. 10 in [1] (Eq. 6.128 in [2])
        H = (1+2.*mu)./(1+2.*gamma.*mu);
    end

    function [H] = exactHFunction(mu, AL_ss)
        %% Hval=exactFunction(mu,AL_ss)
        %   Uses the Stibbs-Weird method [1] to compute the 
        %   Chandrasekhar's H-function for a given set of directional 
        %   parameters (mu) and single-scattering albedos (AL_ss).
        %
        %   Inputs:
        %       - mu (nxm double): matrix of values on the domain 
        %         [0,1] (can be a scalar as well)
        %       - AL_ss (nxm double): matrix of values on the domain 
        %         [0,1] (can be a scalar as well)
        %       - Note: mu and AL_ss mustst have the same dimensions 
        %         unless mu and/or AL_ss is a scalar value
        %   Outputs:
        %       - val (nxm double): matrix of H-function values for 
        %         either the pairwise (if both mu and AL_ss are 
        %         matrices) or for each mu/AL_ss given a scalar 
        %         AL_ss/mu
        %
        %   References:
        %       [1] D. W. N. Stibbs, and R. E. Weir, 
        %           "On the H-Functions for Isotropic Scattering", 
        %           https://doi.org/10.1093/mnras/119.5.512
        %
        %   Last revised: 9/24/24
        %   Last author: Jennifer Nolan
        %
        arguments
            mu       (:,:) double {mustBeInRange(mu,0,1)}
            AL_ss      (:,:) double {mustBeInRange(AL_ss,0,1)}
        end
    
        % Resahpe inputs to vectors for simplified matrix math and
        % store the original dimensions to reshape the output
        inputShape = size(mu);
        mu = mu(:);
        AL_ss = AL_ss(:);
    
        % The Stibbs-Weird [1] approach defines the H-function as:
        % H(mu,AL_ss) = exp[I(mu,AL_ss)], where I is defined by the 
        % integral of a function, f(theta), from 0 to pi/2. Through the
        % Stibbs-Weir approach, additionally, I=I2+I3+I4.These values 
        % will be denoted as I2, I3, and I4 here
        
        % 1) COMPUTE I2
        % The first derivative of f(theta) exhibits problematic
        % behavior at the upper bound of the integral (theta=pi/2) when
        % mu approaches 0. As a result, the integral of f(theta)
        % cannot be evaluated for all mu and AL_ss in its current
        % form. To resolve this issue, a function (g(theta) is defined
        % (Eq. 12 in [1]). The g(theta) function is removed from the 
        % original I integral and added back as a separate integral, 
        % which defines I2 (second integral in Eq. 13 in [1]).

        % The convergence of the series described by Eq. 14.1 and 
        % Eq. 14.2 [1] vary with the value of mu. The point at which 
        % their convergence is equal is when 
        % mu=(1-mu)/(1+mu)=sqrt(2)-1, so this will be the cutoff 
        % value between using Eq. 14.1 and Eq. 14.2 for the inputs
        muConvergenceCutOff = sqrt(2)-1;
        belowConvergenceCutOff = mu<=muConvergenceCutOff & mu>0;
        
        % Define the values of the inputted mu and AL_ss that lie
        % below the convergence cutoff value
        muBelowCutOff = mu(belowConvergenceCutOff);
        AL_ssBelowCutOff = AL_ss(belowConvergenceCutOff);
        
        % Define the values of the inputted mu and AL_ss that lie
        % above the convergence cutoff value
        muAboveCutOff = mu(~belowConvergenceCutOff);
        AL_ssAboveCutOff = AL_ss(~belowConvergenceCutOff);
    
        % The calculation of I2 is contains a summation of an infinite
        % series, as defined in Eq. 14.1 and Eq. 14.2. The number of
        % terms of this infinite series used to compute I2 in this
        % function must be defined. (TODO pick a value).
        nSeriesTerms = 6;
        constantSeriesTerms = (2.*(0:nSeriesTerms)+1);
    
        % Define a matrix to populate the computed values of 
        % I2, where I2 = (1/pi)*integral(g(theta))dtheta from 
        % theta=0 to theta=pi/2 (second integral in Eq. 13 in [1])
        I2 = zeros(size(mu));
    
        % The summation term for (var,AL_ss) pairs less than the mu
        % cutoff is determined (summation in Eq. 14.1 in [1])
        if any(belowConvergenceCutOff)
            belowConvergenceSummation = sum(...
                ((repmat(muBelowCutOff,1,nSeriesTerms+1).^constantSeriesTerms)...
                ./((constantSeriesTerms).^2)),2);
        
        % The full expression defined in Eq. 14.1 in [1] is evaluated for all 
        % (var,AL_ss) pairs below the convergence cutoff
            I2(belowConvergenceCutOff) = ...
                0.5.*AL_ssBelowCutOff.*(...
                0.5.*log(muBelowCutOff).*log((1-muBelowCutOff)./(1+muBelowCutOff))+...
                belowConvergenceSummation);
        end
        
        % The summation term for (var,AL_ss) pairs greater than the mu
        % cutoff is determined (summation in Eq. 14.2 in [2])
        if any(~belowConvergenceCutOff)
            aboveConvergenceSummation = sum(((1./(constantSeriesTerms.^2)).*...
                ((1-(repmat(muAboveCutOff,1,nSeriesTerms+1)))./...
                (1+(repmat(muAboveCutOff,1,nSeriesTerms+1)))).^constantSeriesTerms),2);

        % The full expression defined in Eq. 14.2 in [1] is evaluated for all 
        % (var,AL_ss) pairs above the convergence cutoff
            I2(~belowConvergenceCutOff) = ...
                0.5.*AL_ssAboveCutOff.*((pi^2/8)-aboveConvergenceSummation);
        end

        % When mu=0, the I2 expression for values below the
        % convergence cutoff (Eq. 14/1 in [1]) is undefined. However,
        % the limit of I2 as mu approaches 0 is 0. Therefore, at
        % mu=0, I2=0.
        I2(mu==0)=0;

        % 2) COMPUTE I4
        % The behavior of f(theta) is discontinuous for small theta
        % when the single-scattering albedo, AL_ss, approaches unity. To
        % resolve this discontinuity, an h(theta) is defined as the
        % difference between f(theta) when AL_ss=1 and AL_ss~=1. The
        % h(theta) function is removed from the original I integral and
        % added back as a separate integral, which defines I3 (second
        % integral in Eq. 17 in [1]).
        
        % Define a matrix to populate the computed values of 
        % I4, where I4 = (1/pi)*integral(h(theta))dtheta from 
        % theta=0 to theta=pi/2 (second integral in Eq. 17 in [1])
        I4 = zeros(size(mu));

        % Examining the equation for I4 (Eq. 18 in [1]), it is evident
        % that the value of I4 is undefined when AL_ss=0. The limit
        % demonstrates that I4 approaches zero as AL_ss approaches zero.
        % Therefore, at these undefined points, I4=0.
        undefinedI4 = AL_ss == 0;
        I4(undefinedI4) = 0;
        muDefinedI4 = mu(~undefinedI4);
        AL_ssDefinedI4 = AL_ss(~undefinedI4);
        
        % The full expression defined in Eq. 18 in [1] is evaluated for 
        % all (var,AL_ss) pairs that yield defined values of I4
        I4(~undefinedI4) = (-2.*muDefinedI4./pi).*...
            sqrt(3.*(1-AL_ssDefinedI4)./AL_ssDefinedI4).*...
            atan((pi./2).*sqrt(AL_ssDefinedI4./(3.*(1-AL_ssDefinedI4))));
    
        % 3) COMPUTE I3
        % With the more problematic behavior of f(theta) being managed
        % by the manipulation of g(theta) and h(theta) functions, as
        % seen in in the two previous sections, the value of I3 can be
        % calculated directly by an integral, whose integrand is
        % f(theta)-g(theta)-h(theta) (Eq. 17 in [1]).
        
        % Define a matrix to populate the computed values of 
        % I3, where 
        % I3 = (1/pi)*integral(f(theta)-g(theta)-h(theta))dtheta from 
        % theta=0 to theta=pi/2 (first integral in Eq. 17 in [1])
        I3 = zeros(size(mu));
        
        % The integral of the I3 integrand is evaluated using a 
        % numerical quadrature approach (as recommended by [1]) from 
        % theta=0 to theta=pi/2 for each (mu,AL_ss) pair
        for i=1:size(mu,1)
    
            % Define the computeI3Integrand function, representing 
            % f(theta)-g(theta)-h(theta), as the integrand and specify
            % the evaluation is for the current (mu,AL_ss) pair.
            % Note: see "computeI3Integrand" function for details
            I3Integrand = @(theta)sonic.Math.computeI3Integrand(mu(i),AL_ss(i),theta);
    
            % The quad (integral) functions in MATLAB utilize default
            % tolerance values for achieving sufficiently accurate
            % integral evaluations. This tolerance is expressedly
            % defined here for the numeric quadrature apporach. The
            % default value is TBD because TBD. The quad (integral) is
            % used to evaluate the defined integrand (first integrand
            % in Eq. 17 in [1]) from theta=0 to theta=pi and scalled by
            % (1/pi) to define I3 for the given (var,AL_ss) pair
            tolIntegral = 1E-8;
            I3(i) = (1./pi).*quad(I3Integrand,0,pi./2,tolIntegral);
            % I3(i) = (1./pi).*integral(I3Integrand,0,pi./2,'AbsTol',tol);
        end

        % The Stibbs-Weird [1] approach defines the H-function as:
        % H(mu,AL_ss) = exp[I(mu,AL_ss)], where I is defined by the 
        % integral of a function, f(theta), from 0 to pi/2. Through the
        % Stibbs-Weir approach, additionally, I=I2+I3+I4.These values 
        % will be denoted as I2, I3, and I4 here

        I3(AL_ss==0)=0;

        % 4) COMPUTE I
        % The value of I(mu,AL_ss)=I1+I2+I3 (Eq. 21 in [1])
        I = I2+I3+I4;
        
        % 5) COMPUTE H
        % The value of H(mu,AL_ss) is the exponential of I(mu,AL_ss).
        H = exp(I);

        % Finally, the output H-function values must be returned to the
        % original shape of the inputted of the inputs. This value was
        % preserved at the start of the function and ensure H returns
        % the appropriate nxm matrix expected from the nxm inputs
        H = reshape(H,inputShape);
    end
    
    function [I3Integrand] = computeI3Integrand(mu,AL_ss,theta)
        %% I3Integrand = computeI3Integrand(mu,AL_ss,theta)
        %   Determine the value of the I3 integrand component from the
        %   Stibbs-Weird method [1] to compute the Chandrasekhar's 
        %   H-function. The functions f (Eq. 11.2 in [1]), g (Eq. 12 in
        %   [1]), and h (Eq. 16 in [1]) form the ingtegrand according
        %   to Eq. 17 in [1].
        %
        %   Inputs:
        %       - mu (double): scalar value on the domain [0,1]
        %       - AL_ss (double): scalar value on the domain [0,1]
        %       - theta (1xn double): 
        %   Outputs:
        %       - I3Integrand (1xn matrix): vector of values 
        %         representing the ingtegrand 
        %         f(theta)-g(theta)-h(theta) for each inputted theta
        %
        %   Last revised: 9/26/24
        %   Last author: Jennifer Nolan

        %  f(theta) denotes the integrand of the simplified expression
        % for I(var,AL_ss) (Eq. 10 in [1]). The derivative of the
        % H-function with respect to mu is infinite when mu=1.
        % Similarly, the derivative of the the H-function with respect
        % to AL_ss is infinite when AL_ss=1. The introduction of the
        % g(theta) and h(theta) functions address these issues to allow
        % for the I3 integrand to be evaluated on its expected domain.

        % 1) COMPUTE f(theta)
        % f(theta) denotes the integrand of the simplified expression
        % for I(var,AL_ss) (Eq. 10 in [1]). However, this function  is 
        % discontinuous as AL_ss approaches unity and theta approaches 
        % 0, so additional manipulation is required before integrating

        % Initialize an "f" vector to store the values of f(theta)
        f = zeros(size(theta));
    
        % When theta=0, the function f(theta) is undefined. However,
        % when AL_ss is unity, its limit approaches 2.*mu, such that 
        % the value of f(0)=2.*mu when AL_ss=1 (see note below Eq. 15 
        % in [1])
        f(theta==0 & AL_ss==1) = 2.*mu;
    
        % Even when AL_ss is not unity, the function f(theta) is 
        % undefined when theta=0. In this case, the limit approaches 0, 
        % such that the value of f(0)=0 when AL_ss~=1 (see note below 
        % Eq. 15 in [1])
        f(theta==0 & AL_ss~=1) = 0;
    
        % When theta is not 0, the function f(theta) is defined, and 
        % the value of f(theta) can be directly computed using Eq. 11.2
        % in [1]
        if any(theta~=0)
            thetaNonzero = theta(theta~=0);
            f(theta~=0) = AL_ss.*...
                        (thetaNonzero.*csc(thetaNonzero).^2-cot(thetaNonzero))./...
                        (1-AL_ss.*thetaNonzero.*cot(thetaNonzero)).*...
                        atan(mu.*tan(thetaNonzero));
        end
    
        % 2) COMPUTE g(theta)
        % To address the problematic behavior of f(theta) when 
        % theta=pi/2 and mu approaches 0, a function g(theta) is
        % established which allows a numerical quadrature approach to
        % be used when mu is small and AL_ss is not unity (the first
        % edge case discussed in [1]).

        % Since g(theta) is continuous for all mu, AL_ss, and theta
        % on their respective domains, g(theta) can be computed 
        % directly for all proper inputs using Eq. 12 in [1]
        g = (pi/2).*AL_ss.*atan(mu.*tan(theta));   
    
        % 3) COMPUTE h(theta)
        % To address the problematic behavior of f(theta) as AL_ss 
        % approaches unity, a function h(theta) (not to be confused
        % with the overall H-function) is defined as the difference 
        % between f(theta) when AL_ss=1 and AL_ss~=1. This function
        % allows a numerical quadrature approach to be used when AL_ss
        % is near unity (the second edge case discussed in [1])

        % Initialize an "h" vector to store the values of h(theta)
        h = zeros(size(theta));

        % When AL_ss is unity and theta=0, the function is undefined but
        % its limit approaches unity, such that the value of h(0)=1
        % when AL_ss=1
        h(theta==0 & AL_ss==1) = 1;
    
        % Besides the "theta=0, AL_ss=1" case, the function h(theta) is
        % defined, and the value of h(theta) can be computed using 
        % Eq. 16 in [1]
        h(~(theta==0 & AL_ss==1)) = -2.*mu.*(1-AL_ss)./...
            (1-AL_ss+(1/3).*AL_ss.*(theta(~(theta==0 & AL_ss==1)).^2));
        
        % COMPUTE I3 Integrand
        % Combining the f, g, and h function values allows the I3 to be
        % evaluated using a numerical quadrature approach on the whole
        % domain of mu, AL_ss, and theta with no discontinuity. The
        % integrand in the expression of I3 is
        % f(theta)-g(theta)-h(theta) as shown by Eq. 17 in [1]
        I3Integrand = f-g-h;
    
    end  

    end

end