% This software is made available under the MIT license. See SONIC/LICENSE
% for details.
classdef Aberration

    methods (Static)

        function u_adj = aberrate(u_in, v_KMS)
        %% u_adj = aberrate(u_in, v_KMS)
        %  Implements stellar aberration given the line of sights in frame
        %  a and velocity of b with respect to a, computes line of sights
        %  in frame b
        % 
        %   Inputs:
        %       - u_in (1x1 sonic.Points2 or sonic.PointsS2): either a 
        %         SONIC Points2 object or PointsS2 object containing 
        %         vectors to correct via stellar aberration.
        %       - v_KMS (3x1 double): velocity vector of b 
        %         with respect to a.
        % 
        %   Outputs:
        %       - u_adj (1x1 sonic.Points2 or sonic.PointsS2): SONIC 
        %         PointsS2 object containing the LOS vectors from the 
        %         observer to the stars corrected for stellar aberration  
        %
        %   Last revised: 3/5/24
        %   Last author: Ava Thrasher
    
            % extract points
            switch class(u_in)
                case 'sonic.Points2'
                    % if it's input as p2 points, convert to unit vector
                    u_obs = u_in.p2;
                    u_obs = u_obs./vecnorm(u_obs);
                case 'sonic.PointsS2'
                    % if it's input as spherical grab unit vectors
                    u_obs = u_in.u;
                otherwise
                    error('sonic:Aberration:aberrate:IncorrectInput', ...
                        'First input must be of class sonic.Points2 or sonic.PointsS2.');
            end

            c_KMS = sonic.Constants.c_KMS;
            beta = v_KMS/c_KMS;
            [~,num] = size(u_obs);
            
            % Calculate the Lorentz factor
            gamma = 1/sqrt(1 - (v_KMS'*v_KMS)/c_KMS^2);

            % initialize storage for adjusted observations
            u_adj = zeros(3,num);
    
            for i = 1:num 
                % get observation
                u_i = u_obs(:,i);
    
                % Calculate u adjusted
                coeff = 1/(1 + beta'*u_i);
                BcU = cross(beta, u_i);
                BcBcU = cross(beta, BcU);
                u_adj_i = coeff*(u_i + beta - (1 - gamma)/gamma*BcBcU);
    
                % store adjusted observation
                u_adj_i = u_adj_i./norm(u_adj_i);
                u_adj(:,i) = u_adj_i;
    
            end

            % output aberrated points
            switch class(u_in)
                case 'sonic.Points2'
                    % if it's input as p2 points, output as p2
                    u_adj = sonic.Points2(u_adj);
                case 'sonic.PointsS2'
                    % if it's input as spherical, output as spherical
                    u_adj = sonic.PointsS2(u_adj);
            end
            
        end
        
    end

end