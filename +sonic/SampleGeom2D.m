% This software is made available under the MIT license. See SONIC/LICENSE
% for details.
classdef SampleGeom2D

    methods (Static)

        function [conicPoints] = conicPts(obj, points, thetaRange)
        %% [conicPoints] = conicPts(obj, points, thetaRange)
        %   Creates a Points2 object of points around the conic. The
        %   points are spaced out evenly in terms of angular separation.
        %   
        %   Inputs:
        %       - obj (1x1 sonic.Conic object): Conic object to generate
        %         points on
        %       - points (1,1 double): OPTIONAL. Number of points to 
        %         generate. Defaults to 360 points.
        %       - thetaRange (1,2 double): OPTIONAL. Angular position (in 
        %         radians) to sample the conic from. Defaults to [-pi, pi].
        %   Outputs:
        %       - conicPoints (1x1 sonic.Points2 object): Points2 object
        %         containing the generated points
        %
        %   Last revised: 10/30/24
        %   Last author: Michela Mancini

            arguments
                obj
                points (1,1) double = 360
                thetaRange (1,2) double = [-pi, pi]
            end
            % get start and end theta
            theta1 = thetaRange(1);
            theta2 = thetaRange(2);
            
            % generate angles evenly spaced out
            if theta1 < theta2
                t = linspace(theta1, theta2, points);
            else
                
                t = linspace(mod(theta1,pi)-pi, theta2, points);
            end
            
            % get a and b from the conic
            if ~obj.proper
                error('sonic:SampleGeom2D:invalidInput', ...
                                ['Conic provided is not proper']);
            end

            expVals = obj.explicit;
            x0 = expVals(1);
            y0 = expVals(2);
            a = expVals(3);
            b = expVals(4);
            psi = expVals(5);


            if isinf(a)
                % parabola
                p = sonic.Conic.locusToParameter(obj.locus);
                x = sinh(t).^2;
                y = sqrt(2*p)*sinh(t);
            else
                if a>0
                    % ellipse
                    x = a*cos(t);
                    y = b*sin(t);
                else
                    % hyperbola
                    x = [a*cosh(t),-a*cosh(t)];
                    y = [b*sinh(t),-b*sinh(t)];
                end
            end

            % pair values
            points = [x; y]; %2xn

            % rotate points by psi
            rotPoints = [cos(psi) -sin(psi); sin(psi) cos(psi)]*points;

            % translate points by x0 and y0 and account for origin position
            transformedPoints = rotPoints + [x0;y0];

            % store as a point object
            conicPoints = sonic.Points2(transformedPoints);
        end
    end    
end