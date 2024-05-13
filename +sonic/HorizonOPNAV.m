classdef HorizonOPNAV

    methods (Static)
        function r_C = CRA(limbPts, radii_P, att_P2C)
        %% r_C = CRA(limbPts, radii_P, att_P2C, cam)
        %   Performs horizon based optical navigation using the Christian-Robinson
        %   algorithm as detailed in "A Tutorial on Horizon-Based Optical 
        %   Navigation and Attitude Determination With Space Imaging Systems"
        %
        %   Reference DOI: 10.1109/ACCESS.2021.3051914
        %
        %   Inputs:
        %       - limbPts (1x1 sonic.Points2) Image plane coordinates of 
        %           the limb
        %       - radii_P (3x1 double) Celestial body's principle ai radii 
        %           as described in the planet's principle axis frame
        %       - att_P2C (1x1 sonic.Attitude) Attitude transformation from the 
        %           celestial body's principle axis frame to the camera frame
        % 
        %   Outputs:
        %       - r_C (1x1 sonic.Points3) vector from the camera to the center of the celestial
        %           body expressed in the camera frame
        %
        %   Last revised: 3/27/24
        %   Last author: Ava Thrasher
        
            arguments
                limbPts (1,1) sonic.Points2
                radii_P (3,1) double
                att_P2C (1,1) sonic.Attitude
            end

            % Get rotation matrix
            T_P2C = att_P2C.dcm;
            T_C2P = T_P2C';
        
            % Compute D eq 98
            Dinv = diag(radii_P);
            D = diag(1./radii_P);
        
            % Compute R = DT eq 102
            R = D*T_C2P;
                    
            % Get homogeneous pixel coordinates (p2) 
            % %%
            xbar = limbPts.p2;
        
            % Initialize storage for H
            H = zeros(length(xbar),3);
            % for i 1 to n
            for i = 1:length(xbar)
                % xbar_i' = R ubari eq 95
                xbari = R*xbar(:,i);
        
                % s_i' = xbar_i'/norm(xbar_i') eq 105
                si = xbari./norm(xbari);
                H(i,:) = si'; % eq 109
            end
            
            % Compute LS solution for n eq 108
            n = H\ones(length(H),1);
        
            % Compute rc eq 119
            r_C = (1/sqrt(n'*n - 1))*T_P2C*Dinv*n;
        
            % Output to Points3 object
            r_C = sonic.Points3(r_C);
        end
    end
end