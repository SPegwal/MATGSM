classdef GlobalSkyModelBase
    % Shared base superclass for GlobalSkyModel, GlobalSkyModel2016 and
    % Haslam classes
    properties
        projectionType = 'm'    % See MAAT toolbox celestial.proj.projectcoo for details 
        projectionR = 1         % See MAAT toolbox celestial.proj.projectcoo for details
        projectionPar1 = []     % See MAAT toolbox celestial.proj.projectcoo for details
        projectionPar2 = []     % See MAAT toolbox celestial.proj.projectcoo for details
        
        gridPlotStepDeg(1,1) double = 10;   % Density of grid lines for plots
        FoV = 235;              % Field of view for projections with singularities (deg) 
    end
    
    properties (SetAccess = protected)
        freq_unit(1,:) char {mustBeMember(freq_unit,{'Hz','kHz','MHz','GHz'})} = 'MHz'
        generated_map_data(:,:) double = []     % Always in galactic coordinates, on a HEALpix ring sampling grid
        generated_map_freqs(1,:) double = []
    end
    
    properties (SetAccess = private)
        gridType = 'GalLongLat'        % Can be set to {'Horiz','RAdec','GalLongLat'} setCoorSys
        location(1,3) double = [(-30.721745), (21.411701),  300.0000]  % Earth location in [Lat(deg) Long(deg) mASL]
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
       xyGridAz 
       xyNESW
       xyVerifyMarkers
    end
    
    properties (Dependent = true)
        Nside 
        Nf
        xyHorizon
        xySun
        xyMoon
        julDate
    end
    
    properties (Constant = true, Hidden = true)
        signPhi = 1;
        astroGrids = {'Horiz','RAdec','GalLongLat'}
        verifyMarkers = {'GC','VelaSNR','Cygnus','Cas-A','Cen-A','Tau-A','Orion-A','LMC','SMC'}
        verifyMarkerGalCoors = [0,0; -96,-3.3; 77,2; 111.75,-2.11; -50.5,19.42; -175.42,-5.79; -151,-19.36; -79.5,-32.85; -57.2,-44.3];
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
            
            xyHorizon = [az,alt];
            if ~strcmp(obj.gridType,'Horiz')
                xyHorizon = wrap2pi(celestial.coo.horiz_coo(xyHorizon,obj.julDate,deg2rad(fliplr(obj.location(1:2))),'e'));
                if strcmp(obj.gridType,'GalLongLat')
                    xyHorizon = celestial.coo.coco(xyHorizon,'j2000.0','g','r','r');
                end
            end
            xyHorizon = wrap2pi(xyHorizon);
        end
        
        function xyGridAz = get.xyGridAz(obj)
            % Make vectors of required grid points
            xyGridAz.latVect = deg2rad(-90+obj.gridPlotStepDeg:obj.gridPlotStepDeg:90-obj.gridPlotStepDeg);
            xyGridAz.longVect = deg2rad(-180+obj.gridPlotStepDeg:obj.gridPlotStepDeg:180);
            xyGridAz.Nlat = length(xyGridAz.latVect);
            xyGridAz.Nlong = length(xyGridAz.longVect);
            xyGridAz.Nlines = xyGridAz.Nlat + xyGridAz.Nlong;
            
            % Pack the lines as increasing parallels, and then increasing
            % meridians
            Nplot = 201;
            gridlines = zeros(Nplot,2,xyGridAz.Nlines);
            for nnLat = 1:xyGridAz.Nlat
                az = linspace(-pi,pi,Nplot).';
                alt = ones(size(az)).*xyGridAz.latVect(nnLat);
                gridLines(:,:,nnLat) = [az,alt];
            end
            for nnLong = 1:xyGridAz.Nlong
                alt = linspace(-pi/2,pi/2,Nplot).';
                az = ones(size(alt)).*xyGridAz.longVect(nnLong);
                gridLines(:,:,xyGridAz.Nlat + nnLong) = [az,alt];
            end
           
            gridLines = reshape(gridLines,Nplot*xyGridAz.Nlines,2);
            if ~strcmp(obj.gridType,'Horiz')
                gridLines = wrap2pi(celestial.coo.horiz_coo(gridLines,obj.julDate,deg2rad(fliplr(obj.location(1:2))),'e'));
                if strcmp(obj.gridType,'GalLongLat')
                    gridLines = celestial.coo.coco(gridLines,'j2000.0','g','r','r');
                end
            end
            gridLines = reshape(gridLines,Nplot,2,xyGridAz.Nlines);
            xyGridAz.gridLines = wrap2pi(gridLines);
        end
        
        function xyNESW = get.xyNESW(obj)
            az = deg2rad(0:90:270).';
            alt = zeros(size(az));
            xyNESW = [az,alt];
            if ~strcmp(obj.gridType,'Horiz')
                xyNESW = wrap2pi(celestial.coo.horiz_coo(xyNESW,obj.julDate,deg2rad(fliplr(obj.location(1:2))),'e'));
                if strcmp(obj.gridType,'GalLongLat')
                    xyNESW = celestial.coo.coco(xyNESW,'j2000.0','g','r','r');
                end
            end
            xyNESW = wrap2pi(xyNESW);
        end
        
        function xySun = get.xySun(obj)
            sunStruct = celestial.SolarSys.get_sun(obj.julDate,deg2rad(fliplr(obj.location(1:2))));
            xySun = [sunStruct.RA,sunStruct.Dec];
            switch obj.gridType
                case 'Horiz'
                    xySun = [sunStruct.Az,sunStruct.Alt];
                case 'GalLongLat'
                    xySun = celestial.coo.coco(xySun,'j2000.0','g','r','r');
            end
            xySun = wrap2pi(xySun);
        end
        
        function xyMoon = get.xyMoon(obj)
            moonStruct = celestial.SolarSys.get_moon(obj.julDate,deg2rad(fliplr(obj.location(1:2))));
            xyMoon = [moonStruct.RA,moonStruct.Dec];
            switch obj.gridType
                case 'Horiz'
                    xyMoon = [moonStruct.Az,moonStruct.Alt];
                case 'GalLongLat'
                    xyMoon = celestial.coo.coco(xyMoon,'j2000.0','g','r','r');
            end
            xyMoon = wrap2pi(xyMoon);
        end
        
        function xyVerifyMarkers = get.xyVerifyMarkers(obj)
            
            xyVerifyMarkers = wrap2pi(deg2rad([obj.signPhi,1].*obj.verifyMarkerGalCoors));
            
            if ~strcmp(obj.gridType,'GalLongLat')
                xyVerifyMarkers = celestial.coo.coco(xyVerifyMarkers,'g','j2000.0','r','r');
                xyVerifyMarkers = [wrap2pi(xyVerifyMarkers(:,1)),wrap2pi(xyVerifyMarkers(:,2))];
                if any(strcmp(obj.gridType,{'Horiz'}))  % Update if needed
                    xyVerifyMarkers = wrap2pi(celestial.coo.horiz_coo(xyVerifyMarkers,obj.julDate,deg2rad(fliplr(obj.location(1:2))),'h'));
                    xyVerifyMarkers = xyVerifyMarkers(:,1:2);
                end
            end
        end
        
        function julDate = get.julDate(obj)
            julDate = convert.date2jd([obj.UTCtime.Day,obj.UTCtime.Month,obj.UTCtime.Year,obj.UTCtime.Hour,obj.UTCtime.Minute,obj.UTCtime.Second]);
        end
        
        function obj = setTime(obj,UTCtime)
            % SETTIME sets the time
            % obj = setTime(obj,UTCtime)
            
            obj.UTCtime = UTCtime;
            obj = obj.changeGrid(obj.gridType);  % This always goes from gal -> whereever, so update will happen automatically
        end
        
        function obj = setLocation(obj,location)
            % SETLOCATION sets the location of the observer
            % obj = = setLocation(obj,location)
            
            obj.location  = location;
            obj = obj.changeGrid(obj.gridType);  % This always goes from gal -> whereever, so update will happen automatically
        end
        
        function obj = changeGrid(obj,gridType)
            % CHANGEGRID sets the coordinate system grid
            % obj = changeGrid(obj,gridType)
            % Input: changeGrid can be anything in obj.astroGrids
            
            assert(ismember(gridType,obj.astroGrids),'Unkown coorSys. See obj.astroGrids for allowable names')
            
            if strcmp(gridType,'GalLongLat')
                obj.xy = obj.longlat;
            else %if any(strcmp(coorSys,{'RAdec','Horiz'}))
                % Always calculate this - needed for all three transforms
                equCoords = celestial.coo.coco([obj.longlat],'g','j2000.0','r','r');
                obj.xy = [wrap2pi(equCoords(:,1)),wrap2pi(equCoords(:,2))];
                if any(strcmp(gridType,{'Horiz'}))  % Update if needed
                    horzCoords = wrap2pi(celestial.coo.horiz_coo([obj.xy],obj.julDate,deg2rad(fliplr(obj.location(1:2))),'h'));
                    obj.xy = horzCoords(:,1:2);
                end
            end
            obj.gridType = gridType;
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
        
        function [xy,z] = project(obj,longlat)
            % PROJECT projects the spherical map to 2D plane
            % [xy] = project(obj,longlat)
            % Uses MAAT toolbox celestial.proj.projectcoo
            % Input  : - Matrix of [Longitude,Latitude], in radians.
            % Output : - Matrix of [x,y] positions
            
            assert(size(longlat,2)==2,'Expect 2 column matrix for longlat')
            
            if contains('amhpslbtPGBCxr',obj.projectionType) || isempty(obj.projectionPar1)
                Nin = 4;
            elseif contains('cgMoS',obj.projectionType) || isempty(obj.projectionPar2)
                Nin = 5;
            else
                Nin = 6;
            end
            switch Nin
                case 4
                    [xy(:,1),xy(:,2)] = celestial.proj.projectcoo(longlat(:,1),longlat(:,2),obj.projectionType,obj.projectionR);
                case 5
                    [xy(:,1),xy(:,2)] = celestial.proj.projectcoo(longlat(:,1),longlat(:,2),obj.projectionType,obj.projectionR,obj.projectionPar1);
                case 6
                    [xy(:,1),xy(:,2)] = celestial.proj.projectcoo(longlat(:,1),longlat(:,2),obj.projectionType,obj.projectionR,obj.projectionPar1,obj.projectionPar2);
            end
            z = [];
        end
        
        function plotProj(obj,idx,logged)
            % PLOTPROJ plots a projection of the map
            % plotProj(obj,idx,logged)
            % idx is the frequency index to plot
            % logged is logical to plot in log scale (true)
            
            if nargin < 2 || isempty(idx), idx = 1; end
            if nargin < 3 || isempty(logged), logged = true; end
            
            
            gmap = obj.generated_map_data(:,idx);
            if logged, gmap = log2(gmap); end
            
            if isempty(obj.xy), obj.xy = obj.longlat; end
            
            [xyProj,zProj] = obj.project(obj.xy);
            if strcmp(obj.projectionType,'top')
                xyProj = xyProj(zProj >= 0,:);
                gmap = gmap(zProj >= 0);
            elseif strcmp(obj.projectionType,'bot')
                xyProj = xyProj(zProj <= 0,:);
                gmap = gmap(zProj <= 0);
            end
            gridDelaunay = delaunay(xyProj(:,1),xyProj(:,2));
            h = trisurf(gridDelaunay,xyProj(:,1),xyProj(:,2),gmap.*0,gmap);
            set(h, 'LineStyle', 'none')
            axis equal
            axis off
            view([0,90])
            hold on
            colorbar('location','SouthOutside')
        end
        
        function plotSkyView(obj,idx,logged,details)
            % PLOTSKYVIEW plots the sky view at the current time and place
            % plotSkyView(obj,idx,logged,paramsPlot)
            % idx is the frequency index to plot
            % logged is logical to plot in log scale (true)
            % details is logicals requesting plotting of 
            %   [Directions, Sun, Moon, verifyMarkers]
            
            if nargin < 2 || isempty(idx), idx = 1; end
            if nargin < 3 || isempty(logged), logged = true; end
            if nargin < 4 || isempty(details), details = [0,0,0,0]; end
            
            if numel(details) < 4
                details_ = zeros(1,4);
                details_(1:numel(details)) = details;
                details = details_;
            end
            
            gmap = obj.generated_map_data(:,idx);
            if logged, gmap = log2(gmap); end
            
            if isempty(obj.xy), obj.xy = obj.longlat; end
            
            obj.projectionType = 'S';
            obj.projectionPar1 = deg2rad([0,90]);
            obj = obj.changeGrid('Horiz');
            valInd = find(obj.xy(:,2) >= 0);
            [xyProj,~] = obj.project(obj.xy(valInd,:));
            
            gridDelaunay = delaunay(xyProj(:,1),xyProj(:,2));
            h = trisurf(gridDelaunay,xyProj(:,1),xyProj(:,2),gmap(valInd).*0,gmap(valInd));
            set(h, 'LineStyle', 'none')
            axis equal
            axis off
            view([0,90])
            hold on
            colorbar('location','SouthOutside')
            if details(1), obj.plotDirections('k'); end
            if details(2) && obj.xySun(2) >= 0, obj.plotSun; end
            if details(3) && obj.xyMoon(2) >= 0, obj.plotMoon; end
            if details(4)
                idxVerify = find(obj.xyVerifyMarkers(:,2) >= 0);
                obj.plotVerifyMarkers([],idxVerify);
            end
        end
        
        function plotHorizon(obj,style)
            % PLOTHORIZON plots the horizon on the current axis
            % plotHorizon(obj,style)
            % style is the linestyle of the plot ('w-')
            
            if nargin < 2 || isempty(style), style = 'w-'; end
            
            xyProj = obj.project(obj.xyHorizon);
            plot(xyProj(:,1),xyProj(:,2),style)
        end
        
        function plotGridAz(obj,style)
            % PLOTGRIDAZ plots the azimuthal grid on the current axis
            % plotGridAz(obj,style)
            % style is the linestyle of the plot ('w-')
            
            if nargin < 2 || isempty(style), style = 'w-'; end
            
            for ll = 1:obj.xyGridAz.Nlines
                xyProj = obj.project(obj.xyGridAz.gridLines(:,:,ll));
                plot(xyProj(:,1),xyProj(:,2),style,'linewidth',0.5), hold on
            end
        end
        
        function plotDirections(obj,color)
            % PLOTDIRECTIONS prints the cardinal directions on the current axis
            % plotDirections(obj,color)
            % color is the color of the text ('w')
            
            if nargin < 2 || isempty(color), color = 'w'; end
            
            textStr = 'NESW';
            xyProj = obj.project(obj.xyNESW);
            for tt = 1:4
                text(xyProj(tt,1),xyProj(tt,2),textStr(tt),'Color',color);
            end
        end
        
        function plotSun(obj,style)
            % PLOTSUN plots the approximate sun position
            % plotSun(obj,style)
            % style is the markerstyle of the plot ('w*')
            
            if nargin < 2 || isempty(style), style = 'w*'; end
            
            [xyProj,~] = obj.project(obj.xySun);
            plot(xyProj(:,1),xyProj(:,2),style)
        end

        function plotMoon(obj,style)
            % PLOTMOON plots the approximate moon position
            % plotMoon(obj,style)
            % style is the markerstyle of the plot ('w*')
            
            if nargin < 2 || isempty(style), style = 'wo'; end
            
            [xyProj,~] = obj.project(obj.xyMoon);
            plot(xyProj(:,1),xyProj(:,2),style)
        end
        
        
        function plotVerifyMarkers(obj,color,idxPlot)
            % PLOTVERIFYMARKERS plots the verification markers
            % plotVerifyMarkers(obj,color,idxPlot)
            % color is the text color ('w')
            % idxPlot is a vector of which indexes to plot (ones)
            
            if nargin < 2 || isempty(color), color = 'w'; end
            if nargin < 3 || isempty(idxPlot), idxPlot = ones(size(obj.verifyMarkers)); end
            
            if numel(idxPlot) < length(obj.verifyMarkers)
                idxPlot_ = ones(size(obj.verifyMarkers));
                idxPlot_(1:numel(idxPlot)) = idxPlot;
                idxPlot = idxPlot_;
            end
            
            xyProj = obj.project(obj.xyVerifyMarkers);
            for tt = 1:length(obj.verifyMarkers)
                if idxPlot(tt)
                    text(xyProj(tt,1),xyProj(tt,2),obj.verifyMarkers(tt),'Color',color);
                end
            end
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
%         
    end
end