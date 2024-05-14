classdef BrownConrady < sonic.DistortionModel

    properties
        p1  (1,1) double
        p2  (1,1) double
        k1  (1,1) double
        k2  (1,1) double
        k3  (1,1) double
    end

    methods

        function obj = BrownConrady(p1,p2,k1,k2,k3)
        %% obj = BrownConrady(p1,p2,k1,k2,k3)
        %   Constructs a BrownConrady object, which contains the parameters 
        %   of the Brown Conrady camera distortion model 
        %   
        %   Inputs:
        %       - p1  (1x1 double): Vertical decentering 
        %       - p2  (1x1 double): Horizontal decentering
        %       - k1  (1x1 double): Order 2 radial distortion
        %       - k2  (1x1 double): Order 4 radial distortion
        %       - k3  (1x1 double): Order 6 radial distortion 
        %
        %   Outputs:
        %       - obj (sonic.BrownConrady): BrownConrady object, encoding
        %         the parameters of the distortion model 
        %
        %   Last revised: 2/21/24
        %   Last author: Ava Thrasher
            arguments
                p1  (1,1) double
                p2  (1,1) double
                k1  (1,1) double
                k2  (1,1) double
                k3  (1,1) double
            end
            obj.p1 = p1; 
            obj.p2 = p2;
            obj.k1 = k1;
            obj.k2 = k2;
            obj.k3 = k3;
        end

        function [pointsd] = distort(obj,points)
        %% [pointsd] = distort(obj,points)
        %  Takes image plane coordinates and distorts them based on the 
        %  Brown Conrady model
        %   
        %   Inputs:
        %       - obj (1x1 sonic.BrownConrady)
        %       - points (1x1 sonic.Points2): image plane coordinates of 
        %         points to distort
        %
        %   Outputs:
        %       - pointsd (1x1 sonic.Points2): image plane coordinates of 
        %         distorted points
        %
        %   Last revised: 2/28/24
        %   Last author: Ava Thrasher
            arguments
                obj (1,1) sonic.BrownConrady
                points (1,1) sonic.Points2
            end
            % get x and y coordinates
            x = points.r2(1,:);
            y = points.r2(2,:);

            % precompute x^2, y^2, and x*y
            x2 = x.^2;
            y2 = y.^2;
            xy = x.*y;

            % calculate r squared
            r2 = x2 + y2;

            % find coefficient CHANGE TO HAVE R2 TERM IN, WITH NEW
            % EXPONENTS
            kcoeff = (1 + obj.k1*r2 + obj.k2*r2.^2 + obj.k3*r2.^3);

            % distort points CHANGE TO HAVE R2 TERM
            xd = kcoeff.*x + 2*obj.p1*xy + obj.p2*(r2 + 2*x2);
            yd = kcoeff.*y + obj.p1*(r2 + 2*y2) + 2*obj.p2*xy;

            % package distorted points into a Points2 object
            pointsd = sonic.Points2([xd;yd]);
        end

        function [points, convergeMask] = undistort(obj, pointsd)
        %% [points] = undistort(obj, pointsd)
        %  Takes distorted image plane coordinates and undistorts them 
        %  using a fixed point iteration and provided Brown Conrady model
        %   
        %   Inputs:
        %       - obj (1x1 sonic.BrownConrady)
        %       - pointsd (1x1 sonic.Points2): image plane coordinates of 
        %         distorted points
        %
        %   Outputs:
        %       - points (1x1 sonic.Points2): image plane coordinates of 
        %         points to distort
        %       - convergeMask (1xn boolean): array containing true if the 
        %         i-th point converged, or false if it did not
        %
        %   Last revised: 2/28/24
        %   Last author: Ava Thrasher
            arguments
                obj (1,1) sonic.BrownConrady
                pointsd (1,1) sonic.Points2
            end

            % get x and y distorted coordinates
            xd = pointsd.r2(1,:);
            yd = pointsd.r2(2,:);

            % fixed point iteration to determine undistorted points
            checkPtx = xd;
            checkPty = yd;
            
            % initialize iteration counter
            count = 0;
            % initialize boolean for terminating iteration
            continueIter = true;

            % iterate while errors are above tolerance and iteration count
            % is below the maximum
            while continueIter == true

                % compute distorted guess points 
                outPt = obj.distort(sonic.Points2([checkPtx;checkPty]));
                outPtx = outPt.r2(1,:);
                outPty = outPt.r2(2,:);

                % error between actual distorted points and computed
                errx = outPtx - xd;
                erry = outPty - yd;
                
                % increase iteration counter
                count = count + 1;

                % get convergance mask
                convergeMask = abs(errx) < sonic.Tolerances.UndistortTolBC & abs(erry) < sonic.Tolerances.UndistortTolBC;

                % terminating conditions, assign final estimate if met
                if count == sonic.Tolerances.MaxIters || all(convergeMask)
                    continueIter = false;
                    finalPtx = checkPtx;
                    finalPty = checkPty;
                end
                
                % otherwie apply error to the check point and continue
                checkPtx = checkPtx - errx;
                checkPty = checkPty - erry;

            end

            % once converged, package distorted points into a Points2 object
            points = sonic.Points2([finalPtx;finalPty]);
        end
    end
end