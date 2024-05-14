classdef Ellipsoid < sonic.Quadric

    properties (SetAccess = protected)
        shapeMat (3,3) double
        pos (1,1) sonic.Points3 = sonic.Points3([0;0;0])
        att (1,1) sonic.Attitude = sonic.Attitude(eye(3))
    end

    methods
        function obj = Ellipsoid(principleAxes, att_C2P, rc) 
        %% obj = Ellipsoid(principleAxes, att_C2P, rc)
        %   Constructs an ellipsoid object from the principle axis radii.
        %   If attitude and position of ellipsoid relative to observer are
        %   provided, the ellipsoid will be represented relative to the
        %   observer. Otherwisse, the ellipsoid will be represented
        %   relative to the principle axis frame.
        %
        %   Inputs
        %       - principleAxis (3x1 double): The principle axis radii of
        %         the ellipsoid in x,y,z order
        %       - att_O2P (1,1 sonic.Attitude): OPTIONAL attitude of the 
        %         ellipsoid relative to the camera frame
        %       - rc (1,1 sonic.Points3): OPTIONAL position of the
        %         ellipsoid center relative to the camera
        %
        %   Outputs
        %       - obj (1,1 sonic.Ellipsoid): Object containing the
        %         ellipsoid
            arguments
                principleAxes (3,1) double
                att_C2P       (1,1) sonic.Attitude = sonic.Attitude(eye(3))
                rc           (1,1) sonic.Points3 = sonic.Points3([0;0;0])
            end

            % Construct the shape matrix in the principle axis frame
            shapeMat_P = diag(1./(principleAxes.^2));

            % Rotate shape matrix to camera frame (C)
            shapeMat_C = att_C2P.dcm'*shapeMat_P*att_C2P.dcm; 

            % Transform ellipsoid relative to the camera
            ellipsoid = [shapeMat_C, -shapeMat_C*rc.r3; -rc.r3'*shapeMat_C', rc.r3'*shapeMat_C*rc.r3 - 1];

            % create Quadric from the locus
            obj@sonic.Quadric(ellipsoid,'locus');

            % Set orientation, position, and shape matrix of ellipsoid with
            % respect to the camera
            obj.pos         = rc;
            obj.att         = att_C2P;
            obj.shapeMat    = shapeMat_P;
        end
    end
end