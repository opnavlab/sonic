classdef EllipseFitter
    % This software is made available under the MIT license. 
    % See SONIC/LICENSE for details.
    
    % TODO: change variable names to camel case (in Math too)
    methods(Static)
        function [conicObj, implSoln, explSoln] = ...
                fitEllipse(ellipsePts, method_str)
        %% [conicObj implSoln, explSoln] = fitEllipse(ellipsePts, method_str)
        %
        %   Performs an ellipse fitting to (x,y) points, given a particular
        %   fitting method. Currently can perform least-squares, semi-hyper
        %   least squares and hyper least squares fits.
        %
        %   NOTE: Returns both implicit and explicit solutions, since at 
        %   times the conic's implicit solution won't be stable enough to
        %   convert back to explicit solution without recentering
        %
        %   NOTE: The code in this work is largely based off of the
        %   following references:
        %     1. Kanatani, K., & Rangarajan, P. (2011). "Hyper least 
        %        squares fitting of circles and ellipses." Computational 
        %        Statistics & Data Analysis, 55(6), 2197-2208.
        %     2. Krause, M., Price, J., & Christian, J. A., "Analytical 
        %        Methods in Crater Rim Fitting and Pattern Recognition." 
        %        AAS/AIAA Astrodynamics Specialist Conference, 
        %        Paper AAS 23-350, Big Sky, MT, 12-17 Aug 2023.
        %   
        %   Inputs:
        %       - ellipsePts (sonic.Points2): X,Y vals of 2D points to fit
        %       - method_str (1x1 string): String specifying which type of
        %         ellipse fitting operation to perform. Can be 'ls' for
        %         least-squares, 'shls' for semi-hyper least squares, or
        %         'hls' for hyper least squares 
        %
        %   Outputs:
        %       - conicObj (1,1) sonic.Conic: conic object for solution
        %       - implSoln (6x1): values for the implicit representation
        %         of the ellipse 
        %       - explSoln (5x1): values for the explicit representation 
        %         of the ellipse [xc; yc; a; b; psi]
        %
        %   Last revised: Sep 24, 2024
        %   Last author: Tara Mina
        
            arguments
                ellipsePts         (1, 1)      sonic.Points2
                method_str         (1, 1)      string
            end

            % extract x and y points
            xyvals_2Darr = ellipsePts.r2;
            xvals = xyvals_2Darr(1,:)';
            yvals = xyvals_2Darr(2,:)';

            % compute center of mass (helps with numerical precision issues
            % if conic is too far from the origin) and get centered points
            xcom = mean(xvals); ycom = mean(yvals);
            xcent = xvals-xcom;
            ycent = yvals-ycom;
            
            % get stacked 6D xi vector representation for each 2D point, 
            % which is used for the different least-squares problem
            % formulations
            xiVects = sonic.EllipseFitter.getXiVects(xcent,ycent);

            % fit crater rim based on user's method of choice
            switch method_str
                case 'ls'
                    implSoln = sonic.EllipseFitter.solveTLS(xiVects);
                case 'shls'
                    implSoln = sonic.EllipseFitter.solveSHLS(xiVects);
                case 'hls'
                    implSoln = sonic.EllipseFitter.solveHLS(xiVects);
                otherwise
                    error(['sonic:EllipseFitter:invalidFittingMethod ', ...
                        'Specified crater rim fitting method of ''', ...
                        method_str, ''' is not implemented. Please ', ...
                        'specify one of the impelmented methods ', ...
                        '(''ls'' for least-squares, ''shls'' for ', ...
                        'semi hyper-least squares, or ''hls'' for ', ...
                        'hyper-least squares.']);
            end
            
            % get corresponding explicit solution
            explSoln = sonic.Conic.implicitToExplicit(implSoln);
            
            % add back center of mass
            explSoln = explSoln + [xcom; ycom; 0; 0; 0];
            implSoln = sonic.Conic.explicitToImplicit(explSoln);

            % get and return sonic.Conic object
            conicObj = sonic.Conic(explSoln);
            
            % if conic is not proper, throw warning
            if ~conicObj.proper
                if ( abs(explSoln(3))>sonic.Tolerances.CompZero ...
                        && abs(explSoln(4))>sonic.Tolerances.CompZero )
                    warning(['sonic:EllipseFitter:improperConicButStableExplicit ', ...
                        'Current output Conic object is not proper. '...
                        'Consider using explicit solution vector, '...
                        'which is more numerically stable.']);
                else
                    warning(['sonic:EllipseFitter:improperConic ', ...
                        'Current output Conic object is not proper.'])
                end
            % otherwise, if conic fit is not an ellipse, throw warning
            elseif ~strcmp(conicObj.type,"ellipse")
                warning(['sonic:EllipseFitter:fitNotAnEllipse ', ...
                    'Current fit is not an ellipse, but of type: ', ...
                    conicObj.type])
            end

        end

        function Pa = getEllipseCovariance(implSoln, ellipsePts, sigsqrXY)
        %% Pa = getEllipseCovariance(implSoln, ellipsePts, sigsqrXY)
        %
        %   Compute analytical covariance of implicit solution for
        %   ellipse fitting.
        %
        %   Derivation of HLS covariance originally from Kanatani 
        %   & Rangarajan (2011) and is also readily accessible in 
        %   Krause, Price, & Christian (2023). This function uses the same 
        %   notation as in Krause, et al. (2023).
        % 
        %   NOTE: The choice of normalization factor, N, of the different
        %   ellipse fitting methods affects the bias of the implicit
        %   ellipse solution, but not the standard deviation. Thus, this
        %   same function can be used to compute the covariance for all
        %   ellipse fitting methods.
        %
        %   NOTE: This function implements Equation (45h) in 
        %   Krause, et al. (2023). But, there are 2 typos:
        %      1. last term M_pseudoinverse should be transposed (though
        %         because M is symmetric, its pseudoinverse should also be
        %         symmetric)
        %      2. P_a should have a (1/n^2) factor  multiplied to it 
        %         (should also be in Eqs. 45 b-h), where n is number of  
        %         point samples
        %
        %   Inputs:
        %       - implSoln (6x1 double): implicit solution, as computed by
        %         a particular ellipse fitting method
        %       - ellipsePts (sonic.Points2): X,Y vals of 2D points on 
        %         the edge of the ellipse, to fit
        %       - sigsqrXY -- variance (squared std dev) of x and y
        %         coordinates (assumed that noise in 2D points are 
        %         isometric)
        %
        %   Outputs:
        %       - Pa -- covariance (6x6) of the implicit solution 
        %
        %   Last revised: Sep 24, 2024
        %   Last author: Tara Mina

            arguments
                implSoln            (6, 1)      double
                ellipsePts          (1, 1)      sonic.Points2
                sigsqrXY            (1, 1)      double
            end

            % number of coefficients for implicit ellipse representation
            ncoeff = 6;

            % get xiVects
            xyvals_2Darr = ellipsePts.r2;
            xvals = xyvals_2Darr(1,:)';
            yvals = xyvals_2Darr(2,:)';
            xiVects = sonic.EllipseFitter.getXiVects(xvals,yvals);

            % Compute matrix used for total least squares, called M
            [M, xiOuterAll] = sonic.EllipseFitter.computeLSMatrix(xiVects);
            Mpinv = pinv(M);

            % Get number of sampled points around ellipse
            nsamp = size(xiVects,1);

            % get R0xi (computational effort for this is most of the work
            % to compute the Taubin normalization factor, so we do that)
            xvals = xiVects(:,4);   yvals = xiVects(:,5);
            [~, R0xiAll] = ...
                sonic.EllipseFitter.computeTaubinNormFactor(xvals,yvals);
            RxiAll = sigsqrXY*R0xiAll;

            % Loop through all observations to compute inner matrix
            nobs = length(xvals);
            sumInnerMat = zeros(ncoeff,ncoeff);
            for i=1:nobs
                intermTerm = implSoln' * RxiAll(:,:,i) * implSoln;
                innerMat = (1/nsamp)^2 * intermTerm * xiOuterAll(:,:,i);
                sumInnerMat = sumInnerMat + innerMat;
            end
            
            % Final covariance matrix
            Pa = Mpinv * sumInnerMat * Mpinv';

        end
    end

    methods(Static,Access=private)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% helper functions for computing LS, SHLS, and HLS solutions
        
        function [M, xiOuterAll] = computeLSMatrix(xiVects)
        %% [M, xiOuterAll] = computeLSMatrix(xiVects)
        %
        %   Helper function for ellipse fitting. Computes the least-squares
        %   matrix (also defined in Kanatani & Rangarajan, 2021, as well as
        %   Krause et al., 2023)
        %
        %   Least-squares problem is analogous to finding a which
        %   minimizes a^\top M a (where M is the LS matrix). This matrix is
        %   also utilized for the other fitting methods (e.g., shls, hls).
        %   LS solution defined in EllipseFitter.solveTLS
        %
        %   Also return xiOuterAll, which is the outer-product of the xi
        %   vector, since this term (used to get M), can be reused in HLS
        %
        %   Inputs:
        %       - xiVects (nx6 double): stacked set of 6D vector 
        %         representations of each 2D ellipse point. Definition is
        %         provided in sonic.EllipseFitter.getXiVects function.
        %         Using same notation as in Krause et al., 2023.
        %
        %   Outputs:
        %       - M (6x6): least-squares matrix 
        %       - xiOuterAll (6x6xn): array of 6D outer products (6x6
        %         matrices) for each of the n total 6D xi vectors in
        %         xiVects. This is used to compute M, and can be reused 
        %         downstream for HLS fitting, so we output this term too.
        %
        %   Last revised: Sep 24, 2024
        %   Last author: Tara Mina
        %   Last author for speed-up: Sebastien Henry

            arguments
                xiVects            (:, 6)      double   {mustBeNonempty}
            end

            % number of coefficients for implicit ellipse representation
            ncoeff = 6;

            % get number of samples, or 2D points around the ellipse
            nsamples = size(xiVects, 1);
            
            % Do outer product for each 6D vector in xi (each row of xi),
            % in a vectorized manner

            % convert shape for xi vectors
            xiVectsDeepCols = reshape(xiVects', [ncoeff, 1, nsamples]);
            xiVectsDeepRows = reshape(xiVects', [1, ncoeff, nsamples]);

            % store all outer products individually too, since used in
            % future HLS normalization factor computation
            xiOuterAll = xiVectsDeepCols .* xiVectsDeepRows;

            % get LS matrix from the outer products and normalize
            M = sum(xiOuterAll,3);
            M = (1/nsamples) * M;
        end
        
        function implSoln = solveTLS(xiVects)
        %% implSoln = solveTLS(xiVects)
        %
        %   Get total-least-squares solution, analogous to finding a which
        %   minimizes a^\top M a (where M is the LS matrix). This matrix is
        %   also utilized for the other fitting methods (e.g., shls, hls).
        %
        %   Definition of M in EllipseFitter.computeLSMatrix
        %
        %   Inputs:
        %       - xiVects (nx6 double): stacked set of 6D vector 
        %         representations of each 2D ellipse point. Definition is
        %         provided in sonic.EllipseFitter.getXiVects function.
        %         Using same notation as in Krause et al., 2023.
        %
        %   Outputs:
        %       - implSoln (6x1): values for the implicit representation
        %         of the ellipse
        %
        %   Last revised: Sep 24, 2024
        %   Last author: Tara Mina

            arguments
                xiVects            (:, 6)      double   {mustBeNonempty}
            end

            % get R singular vector corresp to smallest singular value
            [~,implSoln] = sonic.Math.getSmallestRSingVector(xiVects);
            
        end

        function [NTaub, R0xiAll] = computeTaubinNormFactor(xvals,yvals)
        %% [NTaub, R0xiAll] = computeTaubinNormFactor(xvals,yvals)
        %
        %   Get normalization factor for Taubin's ellipse fitting method.
        %   This normalization factor NTaub is also necessary for 
        %   semi-hyper least squares and hyper-least squares methods, since
        %   NTaub is an additive term for the normalization factors of these 
        %   other ellipse fitting methods as well.
        %
        %   NOTE: Taubin's method, Semi-hyper LS, and hyper LS methods aim
        %   to solve the LS problem while constraining the solution vector
        %   such that a^\top N a = constant, where N is a particular 6x6
        %   normalization matrix, which is specific to the ellipse fitting
        %   method. Careful choice of N ensures that the solution is
        %   unbiased. Taubin's method with NTaub as the normalization
        %   factor is still statistically baised, but tends to outperform
        %   unconstrained LS.
        %
        %   The definition of Taubin's normalization factor, NTaub, is
        %   also available in Kanatani & Rangarajan (2011) and Krause, 
        %   Price, & Christian (2023). This function uses the same notation
        %   as in Krause, et al. (2023).
        %
        %   The original Taubin's method corresponds to this reference:
        %   Taubin, G., 1991. "Estimation of planar curves, surfaces, and 
        %   non-planar space curves defined by implicit equations with 
        %   applications to edge and range image segmentation."
        %   IEEE Trans. Pattern Anal. Mach. Intell. 13, 1115–1138.
        %
        %   NOTE: Also returning R0xiAll, since this is computed for the
        %   normalization factor, but can be re-used for the future HLS
        %   ellipse fitting computation as well
        %
        %   Inputs:
        %       - xvals (nx1 double): X vals of 2D points to fit
        %       - yvals (nx1 double): Y vals of 2D points to fit 
        %
        %   Outputs:
        %       - NTaub (6x6): Taubin's normalization matrix 
        %       - R0xiAll (6x6xn): outer product of first differential of
        %         all n 6D xi vectors (def in EllipseFitter.getXiVects), 
        %         which corresponds to the covariance matrix of each 6D xi 
        %         vector, when the original (x,y) pixel variance is 1.
        %         Used for future HLS normalization factor computation, so 
        %         this term is also returned.
        %
        %   Last revised: Sep 24, 2024
        %   Last author: Tara Mina

            arguments
                xvals            (:, 1)      double   {mustBeNonempty}
                yvals            (:, 1)      double   {mustBeNonempty}
            end

            % get Taubin's normalization matrix
            % initialize the normalization factor
            ncoeff = 6;
            NTaub = zeros(ncoeff,ncoeff);
            
            % also return R0xiAll since this can be reused in HLS
            % store R0xiAll (can be used for future HLS comp as well)
            nsamples = length(xvals);
            R0xiAll = zeros(ncoeff,ncoeff,nsamples);
            for i = 1:nsamples
                dxidx = zeros(6,2);
                dxidx(1,1) = 2*xvals(i);
                dxidx(2,1) = yvals(i);
                dxidx(2,2) = xvals(i);
                dxidx(3,2) = 2*yvals(i);
                dxidx(4,1) = 1; dxidx(5,2) = 1;
                R0xi = dxidx * dxidx' ;     R0xiAll(:,:,i) = R0xi;
                NTaub = NTaub + R0xi;
            end
            NTaub = (1/nsamples)*NTaub;
        end
        
        function [Nshls, R0xiAll] = computeSHLSNormFactor(xiVects)
        %% [Nshls, R0xiAll] = computeSHLSNormFactor(xiVects)
        %
        %   Get normalization factor for semi hyper-least squares (SHLS) 
        %   solution for ellipse fitting.
        %
        %   Derivation of SHLS originally from Kanatani & Rangarajan (2011)
        %   and is also readily accessible in Krause, Price, & Christian 
        %   (2023). This function uses the same notation as in 
        %   Krause, et al. (2023).
        %
        %   NOTE: This normalization matrix Nshls is also one of the
        %   factors for the normalization matrix used in HLS (defined in
        %   EllipseFitter.computeHLSNormFactor)
        %
        %   NOTE: Also returning R0xiAll, since this is computed for the
        %   normalization factor, but can be re-used for the future HLS
        %   ellipse fitting computation, as well (which would also
        %   call this function to get Nshls when computing its  
        %   normalization factor).
        %
        %   Inputs:
        %       - xiVects (nx6 double): stacked set of 6D vector 
        %         representations of each 2D ellipse point. Definition is
        %         provided in sonic.EllipseFitter.getXiVects function.
        %         Using same notation as in Krause et al., 2023.
        %
        %   Outputs:
        %       - Nshls (6x6): semi-hyper LS normalization matrix
        %       - R0xiAll (6x6xn): outer product of first differential of
        %         all n 6D xi vectors (def in EllipseFitter.getXiVects), 
        %         which corresponds to the covariance matrix of each 6D xi 
        %         vector, when the original (x,y) pixel variance is 1.
        %         Used for future HLS normalization factor computation, so 
        %         this term is also returned. 
        %
        %   Last revised: Sep 24, 2024
        %   Last author: Tara Mina  

            arguments
                xiVects            (:, 6)      double   {mustBeNonempty}
            end
        
            % extract x and y points
            xvals = xiVects(:,4); 
            yvals = xiVects(:,5);
        
            % compute Taubin's normalization factor, NTaub
            [NTaub, R0xiAll] = ...
                sonic.EllipseFitter.computeTaubinNormFactor(xvals,yvals);
        
            % compute final normalization factor Nshls, using notation 
            % from Krause, et al., (2023)
            xiCent = mean(xiVects,1)';
            eVec = [1; 0; 1; 0; 0; 0];
            B = xiCent * eVec';
            SB = sonic.EllipseFitter.computeSFcn(B);
            Nshls = NTaub + 2*SB;
        end
        
        function implSoln = solveSHLS(xiVects)
        %% implSoln = solveSHLS(xiVects)
        %
        %   Get the semi hyper-least squares (SHLS) solution for ellipse
        %   fitting (return the implicit solution).
        %
        %   Derivation of SHLS originally from Kanatani & Rangarajan (2011)
        %   and is also readily accessible in Krause, Price, & Christian 
        %   (2023). This function uses the same notation as in 
        %   Krause, et al. (2023).
        %
        %   Inputs:
        %       - xiVects (Nx6 double): stacked set of 6D vector 
        %         representations of each 2D ellipse point. Definition is
        %         provided in sonic.EllipseFitter.getXiVects function.
        %         Using same notation as in Krause et al., 2023.
        %
        %   Outputs:
        %       - implSoln (6x1): values for the implicit representation
        %         of the ellipse 
        %
        %   Last revised: Sep 24, 2024
        %   Last author: Tara Mina   

            arguments
                xiVects            (:, 6)      double   {mustBeNonempty}
            end

            % get M, least-squares matrix 
            [M, ~] = sonic.EllipseFitter.computeLSMatrix(xiVects);

            % get normalization factor for SHLS
            [Nshls,~] = sonic.EllipseFitter.computeSHLSNormFactor(xiVects);
        
            % solve generalized eigenvalue problem (for the smallest 
            % magnitude lambda) to get the SHLS solution
            implSoln = ...
                sonic.Math.solveGenEigvalSmallestLambda(M, Nshls);
        
        end

        function [Nhls,M] = computeHLSNormFactor(xiVects)
        %% [Nhls,M] = computeHLSNormFactor(xiVects)
        %
        %   Get normalization factor for hyper-least squares (HLS) 
        %   solution for ellipse fitting.
        %
        %   Derivation of HLS originally from Kanatani & Rangarajan (2011)
        %   and is also readily accessible in Krause, Price, & Christian 
        %   (2023). This function uses the same notation as in 
        %   Krause, et al. (2023). 
        %
        %   NOTE: In Krause, et al., the computation of HLS normalization
        %   factor corresponds to Eq. (44)
        %
        %   NOTE: Also returning M, the least-squares matrix, since this
        %   is necessary to compute the HLS normalization factor, but will 
        %   also be needed again for HLS when solving the generalized
        %   eigenvalue problem
        %
        %   Inputs:
        %       - xiVects (nx6 double): stacked set of 6D vector 
        %         representations of each 2D ellipse point. Definition is
        %         provided in sonic.EllipseFitter.getXiVects function.
        %         Using same notation as in Krause et al., 2023.
        %
        %   Outputs:
        %       - Nhls (6x6): hyper least squares normalization matrix
        %       - M (6x6): least-squares matrix  
        %
        %   Last revised: Sep 24, 2024
        %   Last author: Tara Mina  

            arguments
                xiVects            (:, 6)      double   {mustBeNonempty}
            end
            
            % Get key values necessary to compute the HLS normalization
            % factor, Nhls
            nsamples = size(xiVects, 1);
            [M, xiOuterAll] = ...
                sonic.EllipseFitter.computeLSMatrix(xiVects);
            [Nshls, R0xiAll] = ...
                sonic.EllipseFitter.computeSHLSNormFactor(xiVects);
        
            % compute pseudoinverse of LS matrix M
            Mpinv = pinv(M);

            % Number of coefficients for implicit ellipse representation
            ncoeff = 6;

            % compute summation terms for evaluating Nhls
            sumTerm = zeros(ncoeff, ncoeff);
            for i = 1:nsamples
                % extract pre-computed values we need
                R0xi = R0xiAll(:,:,i);
                xi = xiVects(i,:)';
                xiOuter_i = xiOuterAll(:,:,i);
            
                % compute individual terms of the ith iteration in the sum
                term1 = trace(Mpinv * R0xi)*xiOuter_i;
                term2 = xi'*Mpinv*xi*R0xi;
                term3 = 2*sonic.EllipseFitter.computeSFcn(...
                    R0xi*Mpinv*xiOuter_i);
                
                % add up terms for this iteration of the sum
                sumTerm = sumTerm + term1 + term2 + term3;
            end
            
            % finally, get normalization term by adding this summation term
            Nhls = Nshls - (1/(nsamples*nsamples))*sumTerm;
        end
        
        function implSoln = solveHLS(xiVects)
        %% implSoln = solveHLS(xiVects)
        %
        %   Get the hyper-least squares (HLS) solution for ellipse
        %   fitting (return the implicit solution).
        %
        %   Derivation of HLS originally from Kanatani & Rangarajan (2011)
        %   and is also readily accessible in Krause, Price, & Christian 
        %   (2023). This function uses the same notation as in 
        %   Krause, et al. (2023).
        %
        %   Inputs:
        %       - xiVects (nx6 double): stacked set of 6D vector 
        %         representations of each 2D ellipse point. Definition is
        %         provided in sonic.EllipseFitter.getXiVects function.
        %         Using same notation as in Krause et al., 2023.
        %
        %   Outputs:
        %       - implSoln (6x1): values for the implicit representation
        %         of the ellipse 
        %
        %   Last revised: Sep 24, 2024
        %   Last author: Tara Mina

            arguments
                xiVects            (:, 6)      double   {mustBeNonempty}
            end
            
            % get HLS normalization factor, Nhls
            [Nhls, M] = ...
                sonic.EllipseFitter.computeHLSNormFactor(xiVects);
        
            % solve generalized eigenvalue problem (for the smallest
            % magnitude lambda)
            implSoln = ...
                sonic.Math.solveGenEigvalSmallestLambda(M, Nhls);
        
        end

        function SB = computeSFcn(B)
        %% SB = computeSFcn(B)
        %
        %   Helper function S(.) for matrix symmetrization, used to compute
        %   SHLS & HLS normalization factors (also defined in  
        %   Kanatani & Rangarajan, 2011 and Krause, et al., 2023)
        %
        %   Inputs:
        %       - B (nxn double): square matrix to be symmetrized
        %
        %   Outputs:
        %       - SB (nxn double) -- symmetrized version of B matrix, 
        %         defined by S(B) in Kanatani & Rangarajan (2011), between 
        %         Eqs. (39) and (40)
        %
        %   Last revised: Sep 24, 2024
        %   Last author: Tara Mina

            arguments
                B                    (:, :)      double
            end
            assert( size(B,1)==size(B,2), 'Matrix B must be square');

            SB = 0.5*(B + B');

        end
   
        function xiVects = getXiVects(xvals,yvals)
        %% xiVects = getXiVects(xvals,yvals)
        %
        %   Helper function to get the 6D vector representations of each 2D
        %   ellipse point, which is used to perform all ellipse fitting
        %   methods.
        %
        %   Corresponds to the \xi vector defined in 
        %   Kanatani & Rangarajan (2011) in Eq. (12) as well as in
        %   Krause, et al. (2023) in Eq. (3)
        %
        %   Inputs:
        %       - xvals (nx1 double): X vals of 2D points to fit
        %       - yvals (nx1 double): Y vals of 2D points to fit
        %
        %   Outputs:
        %       - xiVects (nx6 double) -- stacked set of 6D vectors of 
        %         each of the n 2D ellipse points. In terms of a given 2D 
        %         point (x,y), the corresp 6D \xi vector is defined as:
        %            [  x^2  xy  y^2  x  y  1  ]^T
        %
        %   Last revised: Sep 24, 2024
        %   Last author: Tara Mina

        % stacked set of 6D vector 
        %         representations of each 2D ellipse point. Definition is
        %         provided in sonic.EllipseFitter.getXiVects function.
        %         Using same notation as in Krause et al., 2023.

            arguments
                xvals               (:, 1)      double
                yvals               (:, 1)      double
            end
            assert( size(xvals,1)==size(yvals,1), ...
                'Vectors xvals and yvals must be the same length');

            nsamples = size(xvals,1);
            xiVects = [ xvals.*xvals , xvals.*yvals , yvals.*yvals , ...
                xvals , yvals , ones(nsamples,1) ];
        end

    end  % methods
end  % class