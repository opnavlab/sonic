classdef (Abstract) DistortionModel
    methods (Abstract)
        pointsd = distort(obj,points)

        [points, convergeMask] = undistort(obj, pointsd)
        
    end
end