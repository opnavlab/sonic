% This software is made available under the MIT license. See SONIC/LICENSE
% for details.
classdef SphereCoords
    
    methods (Static)        
        function v = raDecToCart(ra_RAD, dec_RAD, r)
        %% v = raDecToCart(ra_RAD, dec_RAD, r)
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
                error('sonic:SphereCoords:raDecToCart:mismatchedValues', ...
                    ['Must specify an equal number of right ascension ' ...
                    'and declination values, as column vectors. If ' ...
                    'specifying an `r` value as well, that must ' ...
                    'also have the same length and be a column vector.']);
            end
        
            v = (r').*([cos(dec_RAD).*cos(ra_RAD), cos(dec_RAD).*sin(ra_RAD), sin(dec_RAD)]');
        end

        function [ra_RAD,dec_RAD,r] = cartToRaDec(v)
        %% [ra_RAD,dec_RAD,r] = cartToRaDec(v)
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

        function cartVec = azZenToCart(az, zen, radi)
        %% cartVec = azZenToCart(az, zen, radius)
        %   Converts azimuth and zenith points on a sphere of radius radi
        %   to 3D vectors in Cartesian coordinates.
        % 
        %   Inputs:
        %       - az (1xn double): vector of azimuth values, or angles wrt
        %         +x axis, in radians
        %       - zen (1xn double): vector of zenith values, or angles wrt
        %         +z axis, in radians (must be between 0 and pi)
        %       - radi (1x1 double): radius of sphere, defaults to 1
        % 
        %   Outputs:
        %       - cartVec (3xn double): corresponding 3D vector in
        %         Cartesian coords (x,y,z)
        %
        %   Last revised: 10/01/24
        %   Last author: Tara Mina
            arguments
                az           (1, :)      double
                zen          (1, :)      double
                radi         (1, 1)      double = 1
            end
            assert(all(zen >= 0 & zen <= pi), ...
                ["sonic:SphereCoords:zenithOutOfBoundsError ", ...
                "all zenith values in zen must be between 0 and pi."])
            cartVec = radi * ...
                [sin(zen).*cos(az); sin(zen).*sin(az); cos(zen)];
        end

        function [az, zen, radi] = cartToAzZen(cartVec)
        %% [az, zen, radi] = cartToAzZen(cartVec)
        %   Converts 3D vectors in Cartesian coordinates to azimuth and 
        %   zenith points on a sphere of radius radi.
        % 
        %   Inputs:
        %       - cartVec (3xn double): corresponding 3D vector in
        %         Cartesian coords (x,y,z)
        % 
        %   Outputs:
        %       - az (1xn double): vector of azimuth values, or angles wrt
        %         +x axis, in radians (between 0 and 2*pi)
        %       - zen (1xn double): vector of zenith values, or angles wrt
        %         +z axis, in radians (between 0 and pi)
        %       - radi (1x1 double): radius of sphere (positive)
        %
        %   Last revised: 10/01/24
        %   Last author: Tara Mina
            arguments
                cartVec      (3, :)      double
            end
            radi = vecnorm(cartVec,2,1);
            radi_xy = vecnorm(cartVec(1:2,:),2,1);
            az = wrapTo2Pi(atan2(cartVec(2,:), cartVec(1,:)));
            zen = pi/2 - atan2(cartVec(3,:), radi_xy);
        end
    end
end

