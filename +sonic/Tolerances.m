% This software is made available under the MIT license. See SONIC/LICENSE
% for details.
classdef Tolerances
    properties (Constant)
        % HomNorm (Points2.m, Points3.m) 
        % In both cases, used for testing if the norm of the homogenous
        % coordinate of a point is "zero". 
        HomNorm             double = 1e-10

        % SmallAngle (Attitude.m, ScanLines.m)
        % Used for determining if an angle is sufficiently small to be
        % treated as a small angle, and approximated accordingly.
        SmallAngle          double = sqrt(eps)

        % DCMDetOne (Attitude.m)
        % Used as a threshold to check if the determinant of a potential
        % DCM is one.
        DCMDetOne           double = 1e-7

        % DCMOrthogCheck (Attitude.m)
        % Used as a threshold to check if a potential DCM is orthogonal via
        % comparing the trace of (T*T' - I) to zero.
        DCMOrthogCheck      double = sqrt(eps)

        % ConicLocusDetOne (Conic.m)
        % Used as a threshold for checking if the locus of a conic has a
        % determinant of one.
        ConicLocusDetOne    double = 1e-14

        % ConicDisc (Conic.m)
        % Used as a threshold for various checks involving the discriminant
        % of a conic. 
        ConicDisc           double = 1e-14

        % ConicCircle (Conic.m)
        % Used as a threshold for checking if a given conic is a circle.
        ConicCircle         double = 1e-14

        % QuatNormOne (Attitude.m)
        % Used as a threshold for checking if the norm of a potential
        % quaternion is one. 
        QuatNormOne         double = 1e-9

        % InfPlaneDistRatio (Plane3.m) 
        % Used for testing the ratio between the normal vector of a plane
        % and the distance. If sufficiently small, then the plane is at
        % infinity.
        InfPlaneDistRatio   double = 1e-10 
                                  
        % psiConicTol (Conic.m)
        % Used for testing how close psi is to zero.
        psiConicTol         double = 1e-10                     

        % undistortTolBC (BrownConrady.m)
        % used to test convergence when undistorting points
        UndistortTolBC double = 1e-12

        % MaxIters (BrownConrady.m)
        % used to cap the amount of iteration when undistorting points
        MaxIters double = 1000      

        % UnitVecNorm1 (PointsS2.m)
        % Used for testing if unit vectors have magnitude 1.
        UnitVecNorm1        double = 1e-7

        % SlopeAdjust (Kvector.m)
        % adjustment of the kvector slope due to floating point precision
        SlopeAdjust         double = 2.22e-16

        % bilinInterpThresh (Image.m)
        % used in bilinear interpolation to see if interploation is linear
        bilinInterpThresh       double = 1e-6

        % CondTol (Conic.m and Quadric.m)
        % used to check if the conic/quadric locus is degenerate
        CondTol                 double = 1e-10;

        % SphereTol (Quadric.m)
        % used to check if the ellipsoid is a sphere
        SphereTol               double = 1e-10;

        % IsOnConic (GeometryP2.m)
        % Used for testing if a point lies on the conic
        IsOnConic double = 1e-7

        % SmallNumber (GeometryP2.m and Conic.m)
        % Used for testing if division by zero occurs
        SmallNumber double = 1e-7
        
        % LineNorm (Lines2.m) 
        % Used for testing if the norm of the first two homogenous
        % coordinates of a line are "zero".
        LineNorm             double = 1e-10

        % CompZero (Math.m) 
        % Used for testing if a given component (real or imaginary) of a 
        % complex number is zero.
        CompZero             double = 1e-12

        % PGRelation (Lines3.m)
        % Used for testing if the Plucker-Grassman relation is satisfied
        % for a given line.
        PGRelation          double = 1e-12

    end
end

