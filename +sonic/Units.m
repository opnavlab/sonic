% This software is made available under the MIT license. See SONIC/LICENSE
% for details.
classdef Units
       
    % If we want to be really slick, we could come up with a more standard
    % unit conversion convention (define conversions as edges along a
    % unit graph), but for now we'll keep it simple.

    methods (Static)        

        function val_DEG = RADtoDEG(val_RAD)
        %% val_DEG = RADtoDEG(val_RAD)
        %   Conversion function from radians to degrees.
        %   Supports vectorized input.
        %   
        %   Inputs:
        %       - val_RAD (1xn double): vector of angles in radians
        %   Outputs:
        %       - val_DEG (1xn double): vector of angles in degrees
        %
        %   Last revised: 2/14/24
        %   Last author: Michael Krause

            val_DEG = val_RAD.*(180/pi);
            
        end

        function val_RAD = DEGtoRAD(val_DEG)
        %% val_RAD = DEGtoRAD(val_DEG)
        %   Conversion function from degrees to radians.
        %   Supports vectorized input.
        %   
        %   Inputs:
        %       - val_DEG (1xn double): vector of angles in degrees
        %   Outputs:
        %       - val_RAD (1xn double): vector of angles in radians
        %
        %   Last revised: 2/14/24
        %   Last author: Michael Krause
        
            val_RAD = val_DEG.*(pi/180);
            
        end

        function val_HMS = DEGtoHMS(val_DEG)
        %% val_HMS = DEGtoHMS(val_DEG)
        %   Conversion function from degrees to Hours/Mins/Seconds.
        %   Supports vectorized input.
        %   
        %   Inputs:
        %       - val_DEG (1xn double): vector of angles in degrees
        %   Outputs:
        %       - val_HMS (3xn double): vector of angles, first row =
        %         hours, second row = minutes, third row = seconds
        %
        %   Last revised: 2/14/24
        %   Last author: Michael Krause

            % Grab the hours:
            hour = floor(val_DEG./15);
    
            % Then the minutes:
            min_raw = 4.*mod(val_DEG, 15);
            min = floor(min_raw);
    
            % Then the seconds:
            sec = 60.*mod(min_raw, 1);
    
            % And combine:
            val_HMS = [hour; min; sec];

        end

        function val_DEG = HMStoDEG(val_HMS)
        %% val_DEG = HMStoDEG(val_HMS)
        %   Conversion function from Hours/Mins/Seconds to degrees.
        %   Supports vectorized input.
        %   
        %   Inputs:
        %       - val_HMS (3xn double): vector of angles, first row =
        %         hours, second row = minutes, third row = seconds
        %   Outputs:
        %       - val_DEG (1xn double): vector of angles in degrees
        %
        %   Last revised: 2/14/24
        %   Last author: Michael Krause

            val_DEG = (val_HMS(1, :).*15) + (val_HMS(2, :)./4) + (val_HMS(3, :)./240);
        end

        function val_ARCSEC = uRADtoARCSEC(val_uRAD)
        %% val_ARCSEC = uRADtoARCSEC(val_uRAD)
        %   Conversion function from microradians to arcseconds.
        %   Supports vectorized input.
        %   
        %   Inputs:
        %       - val_uRAD (1xn double): vector of angles in microradians
        %   Outputs:
        %       - val_ARCSEC (1xn double): vector of angles in arcseconds
        %
        %   Last revised: 2/18/24
        %   Last author: Michael Krause

            val_ARCSEC = val_uRAD.*((3.6e-3)*(180/pi));
            
        end

        function val_uRAD = ARCSECtouRAD(val_ARCSEC)
        %% val_uRAD = ARCSECtouRAD(val_uRAD)
        %   Conversion function from arcseconds to microradians.
        %   Supports vectorized input.
        %   
        %   Inputs:
        %       - val_ARCSEC (1xn double): vector of angles in arcseconds
        %   Outputs:
        %       - val_uRAD (1xn double): vector of angles in microradians
        %
        %   Last revised: 2/18/24
        %   Last author: Michael Krause

            val_uRAD = val_ARCSEC./((3.6e-3)*(180/pi));
            
        end

        function val_DEG = MAStoDEG(val_MAS)
        %% val_DEG = MAStoDEG(val_MAS)
        %   Conversion function from milliarcseconds to degrees.
        %   Supports vectorized input.
        %   
        %   Inputs:
        %       - val_MAS (1xn double): vector of angles in milliarcseconds
        %   Outputs:
        %       - val_DEG (1xn double): vector of angles in degrees
        %
        %   Last revised: 2/14/24
        %   Last author: Michael Krause

            val_DEG = val_MAS./(1000*3600);
            
        end

        function val_RAD = MAStoRAD(val_MAS)
        %% val_RAD = MAStoRAD(val_MAS)
        %   Conversion function from milliarcseconds to radians.
        %   Supports vectorized input.
        %   
        %   Inputs:
        %       - val_MAS (1xn double): vector of angles in milliarcseconds
        %   Outputs:
        %       - val_RAD (1xn double): vector of angles in radians
        %
        %   Last revised: 2/14/24
        %   Last author: Michael Krause

            val_RAD = sonic.Units.DEGtoRAD(...
                sonic.Units.MAStoDEG(val_MAS) ...
            );
            
        end

        function val_PS = PYtoPS(val_PY)
        %% val_PS = PYtoPS(val_PY)
        %   Conversion function from [unit] per year to [unit] per second.
        %   Supports vectorized input.
        %   
        %   Inputs:
        %       - val_PY (1xn double): vector of values per year
        %   Outputs:
        %       - val_PS (1xn double): vector of values per second
        %
        %   Last revised: 2/14/24
        %   Last author: Michael Krause

            val_PS = val_PY./(365.25*24*3600);
            
        end

        function val_RPS = MASPYtoRPS(val_MASPY)
        %% val_RPS = MASPYtoRPS(val_MASPY)
        %   Conversion function from milliarcsecond per year to radians 
        %   per second.
        %   Supports vectorized input.
        %   
        %   Inputs:
        %       - val_MASPY (1xn double): vector of angular rates in
        %         milliarcseconds per year
        %   Outputs:
        %       - val_RPS (1xn double): vector of angular rates in radians
        %         per second
        %
        %   Last revised: 2/14/24
        %   Last author: Michael Krause
        

            val_RPS = sonic.Units.DEGtoRAD(...
                sonic.Units.PYtoPS(...
                    sonic.Units.MAStoDEG(val_MASPY) ...
                ) ...
            );
        
        
        end

        function val_KM = AUtoKM(val_AU)
        %% val_KM = AUtoKM(val_AU)
        %   Conversion function from astronomical units to kilometers, as
        %   defined in the SI Brochure, 9th edition:
        %       https://www.bipm.org/documents/20126/41483022/SI-Brochure-9-EN.pdf
        %   Supports vectorized input.
        %   
        %   Inputs:
        %       - val_AU (1xn double): vector of distances in astronomical
        %         units (AU)
        %   Outputs:
        %       - val_KM (1xn double): vector of distances in kilometers
        %
        %   Last revised: 10/4/24
        %   Last author: Tara Mina

            val_KM = val_AU.*149597870.7;
            
        end

        function val_AU = KMtoAU(val_KM)
        %% val_AU = KMtoAU(val_KM)
        %   Conversion function from kilometers to astronomical units, as
        %   defined in the SI Brochure, 9th edition:
        %       https://www.bipm.org/documents/20126/41483022/SI-Brochure-9-EN.pdf
        %   Supports vectorized input.
        %   
        %   Inputs:
        %       - val_KM (1xn double): vector of distances in kilometers
        %   Outputs:
        %       - val_AU (1xn double): vector of distances in astronomical
        %         units (AU)
        %
        %   Last revised: 10/4/24
        %   Last author: Tara Mina

            val_AU = val_KM./149597870.7;
            
        end


    end
end

