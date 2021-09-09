classdef GlobalSkyModelBase
    % Shared base superclass for GlobalSkyModel, GlobalSkyModel2016 and
    % Haslam classes
    properties (SetAccess = protected)
        freq_unit(1,:) char {mustBeMember(freq_unit,{'Hz','kHz','MHz','GHz'})} = 'MHz'
        generated_map_data(:,:) double = []     % Always in galactic coordinates, on a HEALpix ring sampling grid
        generated_map_freqs(1,:) double = []
    end
    
    properties (SetAccess = private)
        coorSys = 'GalLongLat'        % Can be set to {'Horiz','RAdec','GalLongLat'} setCoorSys
        location(1,3) double = [(-30.721745)*pi/180, (21.411701)*pi/180,  300.0000]  % Earth location in [rad rad mASL]
        UTCtime(1,1) datetime = datetime(2018,1,1,0,0,0)
    end
    
    properties (SetAccess = protected, Hidden = true)
        dataPath
    end
         
    properties (SetAccess = private, Hidden = true)
        xy(:,2) double  % The current local grid
    end
    properties (Dependent = true, Hidden = true)
       freqScale 
       longlat  % Longitude and latitude of the original galactic grid (N x 2)
    end
    
    properties (Dependent = true)
        Nside 
        Nf
        xyHorizon
        xySun
        julDate
    end
    
    properties (Constant = true, Hidden = true)
        signPhi = -1;
        astroGrids = {'Horiz','RAdec','GalLongLat'}
    end
    
    methods
        function freqScale = get.freqScale(obj)
            switch obj.freq_unit
                case 'Hz'
                    freqScale = 1e-6;
                case 'kHz'
                    freqScale = 1e-3;
                case 'MHz'
                    freqScale = 1;
                case 'GHz'
                    freqScale = 1e3;
            end
        end
        
        function Nside = get.Nside(obj)
            Nside = sqrt(size(obj.generated_map_data,1)./12);
        end
        
        function Nf = get.Nf(obj)
            Nf = size(obj.generated_map_data,2);
        end
        
        function longlat = get.longlat(obj)
            tp = pix2ang(obj.Nside);
            tp = [tp{:}];
            th = tp(1,:);
            ph = obj.signPhi.*tp(2,:);
            longlat = [wrap2pi(ph(:)),pi/2 - th(:)];
        end
        
        function xyHorizon = get.xyHorizon(obj)
            az = linspace(-pi,pi,1001).';
            alt = zeros(size(az));
            if strcmp(obj.coorSys,'Horiz')
                xyHorizon = [az,alt];
            else
                xyHorizon = wrap2pi(celestial.coo.horiz_coo([az,alt],obj.julDate,obj.location(1:2),'e'));
                if strcmp(obj.coorSys,'GalLongLat')
                    xyHorizon = celestial.coo.coco(xyHorizon,'j2000.0','g','r','r');
                end
            end
        end
        
        function xySun = get.xySun(obj)
            [RA,dec,~] = celestial.SolarSys.suncoo1(obj.julDate);
            xySun = [RA,dec];
            switch obj.coorSys
                case 'Horiz'
                    xySun = wrap2pi(celestial.coo.horiz_coo(xySun,obj.julDate,obj.location(1:2),'h'));
                case 'GalLongLat'
                    xySun = celestial.coo.coco(xySun,'j2000.0','g','r','r');
            end
            
        end
        
        function julDate = get.julDate(obj)
            julDate = convert.date2jd([obj.UTCtime.Day,obj.UTCtime.Month,obj.UTCtime.Year,obj.UTCtime.Hour,obj.UTCtime.Minute,obj.UTCtime.Second]);
        end
        
        function obj = setTime(obj,UTCtime)
            % SETTIME sets the time
            % obj = setTime(obj,UTCtime)
            
            obj.UTCtime = UTCtime;
            obj = obj.setCoorSys(obj.coorSys);  % This always goes from gal -> whereever, so update will happen automatically
        end
        
        function obj = setLocation(obj,location)
            % SETLOCATION sets the location of the observer
            % obj = = setLocation(obj,location)
            
            obj.location  = location;
            obj = obj.setCoorSys(obj.coorSys);  % This always goes from gal -> whereever, so update will happen automatically
        end
        
        function obj = setCoorSys(obj,coorSys)
            % SETCOORSYS sets the coordinate system
            % obj = setCoorSys(obj,coorSys)
            % Input: coorSys can be anything in obj.astroGrids
            
            assert(ismember(coorSys,obj.astroGrids),'Unkown coorSys. See obj.astroGrids for allowable names')
            
            if strcmp(coorSys,'GalLongLat')
                obj.xy = obj.longlat;
            else %if strcmp(coorSys,'RAdec') || strcmp(coorSys,'Horiz') 
                % Always calculate this - needed for both transforms
                equCoords = celestial.coo.coco([obj.longlat],'g','j2000.0','r','r');
                obj.xy = [wrap2pi(equCoords(:,1)),wrap2pi(equCoords(:,2))];
                if strcmp(coorSys,'Horiz')  % Update if needed
                    horzCoords = wrap2pi(celestial.coo.horiz_coo([obj.xy],obj.julDate,obj.location(1:2),'h'));
                    obj.xy = horzCoords(:,1:2);
                end
            end
            obj.coorSys = coorSys;
        end
        
        function obj = underSample(obj,sampleFactor)
           % UNDERSAMPLE returns the objects on a coarses grid
           % obj = underSample(obj,sampleFactor)
           % SampleFactor is a positive integer indicating 2^sampleFactor reduction
           % in nSide with default (1)
           
           if nargin < 2 || isempty(sampleFactor), sampleFactor  = 1; end
           assert(sampleFactor >= 0 && ~mod(sampleFactor,1),'Expected integer value for sampleFactor') 
           
           % Iteratively make the thing smaller
           for ii = 1:sampleFactor
               % Get nested indexes
               L = size(obj.generated_map_data,1);
               iNest = nest2ring(obj.Nside,(1:L).');
               dNest = zeros(L/4,obj.Nf);
               for ff = 1:obj.Nf
                   % Taking the raw average here - reshape to 4 rows, add them
                   % and normalise by 4
                   dNest(:,ff) = (sum(reshape(obj.generated_map_data(iNest,ff),4,L/4),1)./4).';
               end
               iRing = ring2nest(obj.Nside./2,(1:nSide2nPix(obj.Nside./2)).');
               obj.generated_map_data = dNest(iRing,:);
           end
        end
        
%         function [RAdec] = coor2RAdec(obj)
%             % COOR2RADEC returns the [RA,Dec] coordinates of the object
%             equCoords = celestial.coo.coco([obj.longlat],'g','j2000.0','r','r');
%             RAdec = [wrap22pi(equCoords(:,1)),wrap2pi(equCoords(:,2))]; 
%         end
        
%         function [azAlt,UTCtime,location] = coor2hor(obj,UTCtime,location)
%             % COOR2HOR returns the horizontal (az,alt) coordinates of the
%             % object
%             % [azAlt,locTime,location] = coor2hor(obj,UTCtime,location)
%             %
%             % Inputs:
%             % UTCtime in datetime format (default = 2018/01/01 , 00:00)
%             % location as [long (rad) lat (rad) alt (m)] with default = [0.3292   -0.5922  300.0000] 
%             
%             if nargin < 2 || isempty(UTCtime), UTCtime = datetime(2018,1,1,0,0,0); end
%             if nargin < 3 || isempty(location), location = [0.3292   -0.5922  300.0000]; end
%             
%             julDate = convert.date2jd([UTCtime.Day,UTCtime.Month,UTCtime.Year,UTCtime.Hour,UTCtime.Minute,UTCtime.Second]);
%             horzCoords = wrap2pi(celestial.coo.horiz_coo([obj.coor2RAdec],julDate,location(1:2),'h'));
%             azAlt = horzCoords(:,1:2);
%         end
        
        function plotProj(obj,idx,logged,projectionType,R,par1,par2)
            % PLOTPROJ plots a projection of the map
            % plotProj(obj,idx,longlat,projectionType,R,par1,par2)
            % idx is the frequency index to plot
            % See GlobalSkyModelBase.project for details on the input
            % structure
            
            if nargin < 2 || isempty(idx), idx = 1; end
            if nargin < 3 || isempty(logged), logged = false; end
            if nargin < 4 || isempty(projectionType), projectionType = 'm'; end
            if nargin < 5 || isempty(R), R = 1; end
            if nargin < 6, par1 = []; end
            if nargin < 7, par2 = []; end
            
            gmap = obj.generated_map_data(:,idx);
            if logged, gmap = log2(gmap); end
            
            if isempty(obj.xy), obj.xy = obj.longlat; end
            
            [xyProj] = GlobalSkyModelBase.project(obj.xy,projectionType,R,par1,par2);
            
%             plot3(xy(:,1),xy(:,2),gmap,'.')
            gridDelaunay = delaunay(xyProj(:,1),xyProj(:,2));
            h = trisurf(gridDelaunay,xyProj(:,1),xyProj(:,2),gmap.*0,gmap);
            set(h, 'LineStyle', 'none')
            axis equal
            % axis off
            view([0,90])
            hold on
        end
        
        function plotHorizon(obj,projectionType,R,par1,par2,style)
            % PLOTHORIZON plots the horizon on the current axis
            % plotHorizon(obj,projectionType,R,par1,par2,style)
            % See plotProj for input details.
            % style is the linestyle of the plot ('w-')
            
            % TODO: Fix this
            
            if nargin < 2 || isempty(projectionType), projectionType = 'm'; end
            if nargin < 3 || isempty(R), R = 1; end
            if nargin < 4, par1 = []; end
            if nargin < 5, par2 = []; end
            if nargin < 6 || isempty(style), style = 'w-'; end
            
            xyProj = GlobalSkyModelBase.project(obj.xyHorizon,projectionType,R,par1,par2);
            plot(xyProj(:,1),xyProj(:,2),style)
        end
        
        function plotSun(obj,projectionType,R,par1,par2,style)
            % PLOTSUN plots the approximate sun position
            
            % TODO: Fix this
            
            if nargin < 2 || isempty(projectionType), projectionType = 'm'; end
            if nargin < 3 || isempty(R), R = 1; end
            if nargin < 4, par1 = []; end
            if nargin < 5, par2 = []; end
            if nargin < 6 || isempty(style), style = 'w*'; end
            
            xyProj = GlobalSkyModelBase.project(obj.xySun,projectionType,R,par1,par2);
            plot(xyProj(:,1),xyProj(:,2),style)

        end
        
        function view(obj, idx, logged)
            %     View generated map using mollweide projection.
            % 
            %     Parameters
            %     ----------
            %     idx: int (1)
            %         index of map to view. Only required if you generated maps at
            %         multiple frequencies.
            %     logged: logical (false)
            %         Take the log of the data before plotting. Defaults to
            %         False..
            
            
           assert(~isempty(obj.generated_map_data),'No GSM map has been generated yet. Run generate() first.')
           if nargin > 1 && ~isempty(idx)
               assert(numel(idx) == 1,'Scalar idx expected')
               gmap = obj.generated_map_data(:,idx);
               freq = obj.generated_map_freqs(idx);
           else
               gmap = obj.generated_map_data(:,1);
               freq = obj.generated_map_freqs(1);
           end
           
           if nargin < 3, logged = false; end
           if logged, gmap = log2(gmap); end
           
           healpixPlotMollweide(gmap)
%            title(['Global Sky Model at ', num2str(freq), ' MHz from the ', obj.basemap, ' map'])
        end
        
        function write_fits(obj,filename)
           fitswrite(obj.generated_map_data,filename); 
        end
    end
    
    methods (Abstract = true)
        [obj,map_out] = generate(obj,freqs)
        
    end
    
    methods (Static = true)
        function [xy] = project(longlat,projectionType,R,par1,par2)
            % PROJECT projects the spherical map to 2D plane
            % [xy] = project(obj,projectionType,R,par1,par2)
            % Uses MAAT toolbox celestial.proj.projectcoo
            % Input  : - Matrix of [Longitude,Latitude], in radians.
            %          - Projection Type
            %            'a' - Aitoff (default)
            %            'm' - Mollweide equal-area
            %            'h' - Hammer
            %            'p' - Parabolic
            %            's' - Sinusoidal
            %            'l' - Lambert equal-area cylindrical
            %            'b' - Behrmann equal-area cylindrical
            %            't' - Tristan Edwards cylindrical
            %            'P' - Peters cylindrical
            %            'G' - Gall Orthographic cylindrical
            %            'B' - Balthasart cylindrical
            %            'c' - General cylindrical, opt par 1 is [Stand_Long, Stand_Lat] (radians)
            %            'C' - Cassini
            %            'x' - XY projection, no transformation.
            %            'r' - polar projection (from north pole).
            %            'A' - Albers equal-area, Par 1 is CenCoo and 2 is for ParLat (radians)
            %            'g' - Gnomonic nonconformal projection. Par1 is [Long_cen, Lat_cen]
            %            'M' - Mercator projection. Par1 is long_cen.
            %            'o' - Bonne projection. Par1 is [central_long, standard_parallel]
            %            'S' - Stereographic projection. Par1 is [central_long, standard_parallel]
            %            #.# - Conic projection. where #.# is height of apex.
            %          - Radius scale parameters, default is 1.
            %          - Optional parameters 1.
            %          - Optional parameters 2.
            % Output : - Matrix of [x,y] positions
            
            assert(size(longlat,2)==2,'Expect 2 column matrix for longlat')
            if nargin < 4, par1 = []; end
            if nargin < 5, par2 = []; end
            
            if contains('amhpslbtPGBCxr',projectionType) || isempty(par1)
                Nin = 4;
            elseif contains('cgMoS',projectionType) || isempty(par2)
                Nin = 5;
            else
                Nin = 6;
            end
            switch Nin
                case 4
                    [xy(:,1),xy(:,2)] = celestial.proj.projectcoo(longlat(:,1),longlat(:,2),projectionType,R);
                case 5
                    [xy(:,1),xy(:,2)] = celestial.proj.projectcoo(longlat(:,1),longlat(:,2),projectionType,R,par1);
                case 6
                    [xy(:,1),xy(:,2)] = celestial.proj.projectcoo(longlat(:,1),longlat(:,2),projectionType,R,par1,par2);
            end
            
        end
        
    end
end