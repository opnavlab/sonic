% This software is made available under the MIT license. See SONIC/LICENSE
% for details.
classdef (Abstract) DistortionModel
    methods (Abstract)
        pointsd = distort(obj,points)

        [points, convergeMask] = undistort(obj, pointsd)
        
    end
end