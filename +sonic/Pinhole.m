classdef Pinhole < sonic.DistortionModel

    methods

        function obj = Pinhole()
        %% obj = Pinhole()
        %   Constructs a Pinhole distortion model, which is simply a
        %   "pass-through" distortion model.
        %   
        %   Inputs:
        %       None.
        %
        %   Outputs:
        %       - obj (1x1 sonic.Pinhole): Pinhole distortion model object
        %
        %   Last revised: 03/06/24
        %   Last author: Michael Krause
        end

        function [pointsd] = distort(obj,points)
        %% [pointsd] = distort(obj,points)
        %  Takes image plane coordinates and distorts them based on the 
        %  Pinhole ("passthrough") model.
        %   
        %   Inputs:
        %       - obj (1x1 sonic.Pinhole): Pinhole distortion model object
        %       - points (1x1 sonic.Points2): image plane coordinates of 
        %           points to distort
        %
        %   Outputs:
        %       - pointsd (1x1 sonic.Points2): image plane coordinates of 
        %           distorted points
        %
        %   Last revised: 03/06/24
        %   Last author: Michael Krause
            arguments
                obj         (1,1)       sonic.Pinhole
                points      (1,1)       sonic.Points2
            end

            pointsd = points;
        end

        function [points, convergeMask] = undistort(obj, pointsd)
        %% [points, convergeMask] = undistort(obj, pointsd)
        %  Takes "distorted" image plane coordinates and undistorts them 
        %  (notably, no real undistortion occurs since this is a 
        %  passthrough model).
        %   
        %   Inputs:
        %       - obj (1x1 sonic.Pinhole): Pinhole distortion model object
        %       - pointsd (1x1 sonic.Points2): image plane coordinates of 
        %           distorted points
        %
        %   Outputs:
        %       - points (1x1 sonic.Points2): image plane coordinates of 
        %           undistorted points
        %       - convergeMask (1xn logical): array containing true if the 
        %           ith point converged, or false if it did not. This is to
        %           maintain syntax compatibility with other distortion
        %           models, and is always true for this model. 
        %
        %   Last revised: 03/06/24
        %   Last author: Michael Krause
            arguments
                obj         (1,1)       sonic.Pinhole
                pointsd     (1,1)       sonic.Points2
            end

            points = pointsd;
            convergeMask = true([1, points.n]);

        end
    end
end