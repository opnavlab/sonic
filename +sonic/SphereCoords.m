classdef SphereCoords
    
    methods (Static)
        
        function v = sphereToCart(ra_RAD, dec_RAD, r)
        %% v = sphereToCart(ra_RAD, dec_RAD, r)
        %   Converts equally sized vectors of right ascension and declination 
        %   (given in radians) to vectors. If no `r` argument is supplied,
        %   it will convert them to unit vectors.
        %   
        %   Inputs:
        %       - ra_RAD (nx1 double): Column vector of right ascension 
        %         values, in radians.
        %       - dec_RAD (nx1 double): Column vector of declination 
        %         values, in radians.
        %       - r (nx1 double): OPTIONAL: Column vector of vector
        %         magnitude, defaults to unit magnitude if not specified.
        %   Outputs:
        %       - v (3xn double): Matrix of resultant vectors, with the vectors 
        %         making up the columns. The n-th vector corresponds to the
        %         n-th (RA, DEC) pair. If r is not specified, these are
        %         all unit vectors. 
        %
        %   Last revised: 03/01/24
        %   Last author: Michael Krause

            arguments
                ra_RAD      (:, 1)      double
                dec_RAD     (:, 1)      double
                r           (:, 1)      double      = 1
            end

            % Check to make sure we have equal lengths of our vectors. 
            veclens = [length(ra_RAD), length(dec_RAD), length(r)];
            if ~(all(veclens == veclens(1)) || ...
                    (veclens(1) == veclens(2) && veclens(3) == 1))
                error('sonic:SphereCoords:sphereToCart:mismatchedValues', ...
                    ['Must specify an equal number of right ascension ' ...
                    'and declination values, as column vectors. If ' ...
                    'specifying an `r` value as well, that must ' ...
                    'also have the same length and be a column vector.']);
            end
        
            v = (r').*([cos(dec_RAD).*cos(ra_RAD), cos(dec_RAD).*sin(ra_RAD), sin(dec_RAD)]');
        end

        function [ra_RAD,dec_RAD,r] = cartToSphere(v)
        %% [ra_RAD,dec_RAD,r] = cartToSphere(v)
        %   Converts vectors (v) in cartesian coordinates to spherical 
        %   coordinates given by right ascension (RA) and declination 
        %   (DEC), both in radians. Also returns the magnitude of each
        %   vector.
        % 
        %   Inputs:
        %       - v (3xn double): Matrix of vectors, with the vectors making
        %         up the columns. The n-th unit vector corresponds to the
        %         n-th (RA, DEC) pair.
        % 
        %   Outputs:
        %       - ra_RAD (nx1 double): Column vector of right ascension 
        %         values, in radians.
        %       - dec_RAD (nx1 double): Column vector of declination 
        %         values, in radians.
        %       - r (nx1 double): Magnitude of each vector.
        %
        %   Last revised: 03/04/24
        %   Last author: Michael Krause

            arguments
                v           (3, :)      double
            end

            % Grab the magnitudes and output them as columns:  
            r = vecnorm(v);
            r = r';     % Flip this to be a column vector
            
            % Calculate RA and DEC using components of u. Since we're using
            % atan2, don't need to normalize the vectors, they'll be fine.
            dec_RAD = atan2(v(3, :), sqrt(v(1, :).^2 + v(2, :).^2));
            ra_RAD = atan2(v(2, :), v(1, :));
            
            % Transpose to convert back to column vectors:
            dec_RAD = dec_RAD';
            ra_RAD = ra_RAD';
        
        end

    end
end

