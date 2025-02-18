classdef Haslam < GlobalSkyModelBase
    
    properties
        
    end
    
    properties (SetAccess = private)
        spectral_index(1,1) double {mustBeReal,mustBeFinite,mustBeNegative} = -2.6
        data(:,1) double
        
        Tcmb(1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative} = 2.73
    end
    
    methods 
        function obj = Haslam(freq_unit,spectral_index)
            % HASLAM class constructor method
            % Reproduces the haslam.py class in pygdsm
            %
            % Notes
            % -----
            % This is a crude model; the sky's spectral index changes with sidereal time
            % and a singular spectral index is not realistic. Here are some measured values
            % at different frequencies:
            %
            %     Frequency        Spectral Index        Reference     Notes
            %     50--100 MHz     si = -2.56 +/- 0.03    Mozdzen+19    LST 0--12 hr
            %                           2.46 +/- 0.01                  LST 18.2 hr
            %     90--190 MHz     si = -2.61 +/- 0.01    Mozdzen+16    LST 0--12 hr
            %                          −2.50 +/- 0.02                  LST 17.7 hr
            %     0.4--4.7 GHz 	si = −2.85 +/- 0.05    Dickinson+19
            %     0.4--22.8 GHz 	si = −2.88 +/- 0.03    Dickinson+19
            %     4.7--22.8 GHz 	si = −2.91 +/- 0.04    Dickinson+19
            %     22.8--44 GHz 	si = −2.85 +/- 0.14    Dickinson+19
            %
            % References
            % ----------
            % Mozdzen et al (2016), EDGES high-band, https://doi.org/10.1093/mnras/stw2696
            % Mozdzen et al (2019), EDGES low-band, https://doi.org/10.1093/mnras/sty3410
            % Dickinson et al (2019), C-BASS experiment, https://doi.org/10.1093/mnras/stz522
            %
            % Inputs
            % - freq_unit: {'Hz',('MHz'),'GHz'}
            % - spectral_index: negative real number (-2.6)
            %
            % Outputs
            % - obj:    Haslam object
            %
            % Dependencies
            % -
            %
            % Created: 2021-09-06, Dirk de Villiers
            % Updated: 2021-09-06, Dirk de Villiers
            %
            % Tested : Matlab R2021a
            %  Level : 1
            %   File : 
            %
            % Example
            %   Hmap = Haslam;
            %   Hmap = Hmap.generate(408);
            %   Hamp.view(1,true)
        
            if nargin > 0 && ~isempty(freq_unit), obj.freq_unit = freq_unit; end
            if nargin > 1 && ~isempty(spectral_index), obj.spectral_index = spectral_index; end
            
            % Set the path
            P = fileparts(mfilename('fullpath'));
            obj.dataPath = [P,'\data\haslam408_dsds_Remazeilles2014.fits'];
            
            % Load the data
            d_ =  fitsread(obj.dataPath,'BinaryTable');
            d_ = transpose(d_{1});
            obj.data = d_(:);
        end
        
        function [obj, map_out] = generate(obj,freqs)
            %     [obj, map_out] = generate(obj,freqs) 
            %     Generate a global sky model at a given frequency or frequencies
            %     
            %     Parameters
            %     ----------
            %     freq: scalar or array
            %     Frequency for which to return GSM model (unit as in the object)
            %     
            %     Returns
            %     -------
            %     map_out
            %     Global sky model in healpix format, with NSIDE=512. Output map
            %     is in galactic coordinates, and in antenna temperature units (K).
            
            assert(min(size(freqs))  == 1, 'freqs must be vector')
            freqs_mhz = freqs.*obj.freqScale;
            
            map_out = (obj.data - obj.Tcmb).*(freqs_mhz./408.0).^(obj.spectral_index) + obj.Tcmb;
            
            obj.generated_map_data = map_out;
            obj.generated_map_freqs = freqs;
        end
        
        function obj = setTcmb(Tcmb)
            % SETTCMB sets the CMB temperature in the object
            % obj = setTcmb(Tcmb)
            % If Tcmb = 0 the full temperature is scled by the power law.
            % If Tcmb > 0, it is removed from the scaling and added back.
            % Default is 2.73 K
            
            obj.Tcmb = Tcmb;
            if ~isempty(obj.generated_map_data)
                obj = obj.generate(obj.generated_map_freqs);
            end
        end
    end
    
end