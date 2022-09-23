classdef GlobalSkyModelBase_II
    % Shared base superclass for GlobalSkyModel, GlobalSkyModel2016 and
    % Haslam classes
    properties
        projectionType = 'm'    % See MAAT toolbox celestial.proj.projectcoo for details (top/bot is direction cosine hemispheres)
        projectionR = 1         % See MAAT toolbox celestial.proj.projectcoo for details
        projectionPar1 = []     % See MAAT toolbox celestial.proj.projectcoo for details
        projectionPar2 = []     % See MAAT toolbox celestial.proj.projectcoo for details
    end
    
    properties (SetAccess = protected)
        freq_unit(1,:) char {mustBeMember(freq_unit,{'Hz','kHz','MHz','GHz'})} = 'MHz'
        generated_map_data(:,:) double = []     % Always in galactic coordinates, on a HEALpix ring sampling grid
        generated_map_freqs(1,:) double = []
    end
    
    properties (SetAccess = public)
        coorSys = 'GalLongLat'        % Can be set to {'Horiz','RAdec','GalLongLat'} setCoorSys
        location(1,3) double %= [(-30.721745), (21.411701),  300.0000]  % Earth location in [Lat(deg) Long(deg) mASL]
        UTCtime(1,1) datetime = datetime(2018,1,1,0,0,0)
    end
    
    properties (SetAccess = protected, Hidden = true)
        dataPath
    end
         
    properties (SetAccess = private, Hidden = true)
        xy(:,2) double  % The current local grid
        markers_xy (:,2) double %current markers location
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
        xyMoon
        julDate
        signPhi
    end
    
    properties (Constant = true, Hidden = true)
        
        astroGrids = {'Horiz','RAdec','GalLongLat'}
        verifyMarkers = {'GC','VelaSNR','Cygnus','Cas-A','Cen-A','Tau-A','Orion-A','LMC','SMC','SUN','MOON'};
        verifyMarkerGalCoors = [0,0; -96,-3.3; 77,2; 111.75,-2.11; -50.5,19.42; -175.42,-5.79; -151,-19.36; -79.5,-32.85; -57.2,-44.3].*(pi/180);

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
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@        
        function signPhi = get.signPhi(obj)
            switch obj.coorSys
                case {'GalLongLat', 'RAdec'}
                    signPhi = -1;
                case 'Horiz'
                    signPhi = 1;
            end
        end

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@        
        function Nside = get.Nside(obj)
            Nside = sqrt(size(obj.generated_map_data,1)./12);
        end
        
        function Nf = get.Nf(obj)
            Nf = size(obj.generated_map_data,2);
        end
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@        
        function longlat = get.longlat(obj)
            tp = pix2ang(obj.Nside);
            tp = [tp{:}];
            th = tp(1,:);
            ph = tp(2,:);%obj.signPhi.*
            longlat = [wrapToPi(ph(:)),pi/2 - th(:)];
        end
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        function xyHorizon = get.xyHorizon(obj)
            az = linspace(min(obj.longlat(:,1)),max(obj.longlat(:,1)),1001).'; 
            alt = zeros(size(az));
            if strcmp(obj.coorSys,'Horiz')
                xyHorizon = [az,alt];%equator/horizon...horizontal coordinate
            else
                %changing position of equator/horizon from hoirzontal to equatorial coordinates 
                xyHorizon = wrapToPi(horiz_coo([az,alt],obj.julDate,deg2rad(fliplr(obj.location(1:2))),'e'));
                if strcmp(obj.coorSys,'GalLongLat')
                    %changing position of equator/horizon from equatorial to galactic coordinates 
                    xyHorizon = celestial.coo.coco(xyHorizon,'j2000.0','g','r','r');
                end
            end
            xyHorizon = wrapToPi([obj.signPhi,1].*xyHorizon);
        end
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@        
        function xySun = get.xySun(obj)
            sunStruct = celestial.SolarSys.get_sun(obj.julDate,deg2rad(fliplr(obj.location(1:2))));%This gives sun's location in equitorial coordinates and horizontal coordinates
            xySun = [sunStruct.RA,sunStruct.Dec];
            switch obj.coorSys
                case 'Horiz'
                    xySun = [sunStruct.Az,sunStruct.Alt];
                case 'GalLongLat'
                    xySun = celestial.coo.coco(xySun,'j2000.0','g','r','r');%converting sun's position into galactic coordinates
            end
            xySun = wrapToPi([obj.signPhi,1].*xySun);
            
        end
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@              
        function xyMoon = get.xyMoon(obj)
            moonStruct = celestial.SolarSys.get_moon(obj.julDate,deg2rad(fliplr(obj.location(1:2))));
            xyMoon = [moonStruct.RA,moonStruct.Dec];
            switch obj.coorSys
                case 'Horiz'
                    xyMoon = [moonStruct.Az,moonStruct.Alt];
                case 'GalLongLat'
                    xyMoon = celestial.coo.coco(xyMoon,'j2000.0','g','r','r');
            end
            xyMoon = wrapToPi([obj.signPhi,1].*xyMoon);
            
        end
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@              
        function julDate = get.julDate(obj)
            julDate = convert.date2jd([obj.UTCtime.Day,obj.UTCtime.Month,obj.UTCtime.Year,obj.UTCtime.Hour,obj.UTCtime.Minute,obj.UTCtime.Second]);
        end
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@        
        function obj = setTime(obj,UTCtime)
            % SETTIME sets the time
            % obj = setTime(obj,UTCtime)
            
            obj.UTCtime = UTCtime;
            obj = obj.setCoorSys(obj.coorSys);  % This always goes from gal -> whereever, so update will happen automatically
        end
   
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@        
        function obj = setCoorSys(obj,coorSys)
            % SETCOORSYS sets the coordinate system
            % obj = setCoorSys(obj,coorSys)
            % Input: coorSys can be anything in obj.astroGrids
            
            assert(ismember(coorSys,obj.astroGrids),'Unkown coorSys. See obj.astroGrids for allowable names')
            
            if strcmp(coorSys,'GalLongLat')
                obj.xy = [obj.signPhi,1].*obj.longlat;
                obj.markers_xy = wrapToPi([obj.signPhi,1].*obj.verifyMarkerGalCoors);
                obj.markers_xy = [obj.markers_xy;obj.xySun;obj.xyMoon];
            else
                %if strcmp(coorSys,'RAdec') || strcmp(coorSys,'Horiz')
                % Always calculate this - needed for both transforms
                
                equCoords = celestial.coo.coco([obj.longlat],'g','j2000.0','r','r');
                obj.xy = wrapToPi([obj.signPhi,1].*equCoords);%[wrapToPi(equCoords(:,1)),wrapToPi(equCoords(:,2))]; obj.xy= [obj.signPhi,1].*obj.xy;   
                
                equCoords_markers = celestial.coo.coco([obj.verifyMarkerGalCoors],'g','j2000.0','r','r');
                obj.markers_xy = wrapToPi([obj.signPhi,1].*equCoords_markers);%[wrapToPi(equCoords_markers(:,1)),wrapToPi(equCoords_markers(:,2))];

                if strcmp(coorSys,'Horiz')  % Update if needed
                    horzCoords = wrapToPi(horiz_coo([obj.xy],obj.julDate,deg2rad(fliplr(obj.location(1:2))),'h'));
                    obj.xy = [obj.signPhi,1].*horzCoords(:,1:2);

                    
                    horzCoords_markers = wrapToPi(horiz_coo([obj.markers_xy],obj.julDate,deg2rad(fliplr(obj.location(1:2))),'h'));
                    obj.markers_xy = ([obj.signPhi,1].*horzCoords_markers(:,1:2));

                end
                obj.markers_xy = [obj.markers_xy;obj.xySun;obj.xyMoon];
            end
            obj.coorSys = coorSys;
        end
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@        
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
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@       
        function plotProj(obj,idx,logged)
            % plotProj plots a projection of the map
            % plotProj(obj,idx,longlat)
            % idx is the frequency index to plot
            % See GlobalSkyModelBase.project for details on the input
            % structure
            %nargin
            if nargin < 2 || isempty(idx), idx = 1;
            end
            if nargin < 3 || isempty(logged), logged = false;
            end
            
            
            gmap = obj.generated_map_data(:,idx);
            if logged, gmap = log2(gmap);
            end
            
            if isempty(obj.xy)
                obj.xy = [obj.signPhi,1].*obj.longlat;   
            
            end


            if strcmp(obj.projectionType,'top')
                [xy_top,xy_top_m]=obj.rotate(0,-90,0);
                axesm stereo                 
                scatterm(xy_top(:,2)*(180/pi),xy_top(:,1)*(180/pi),12,gmap,'filled');
                colorbar ('southoutside');tightmap;axis off;
                scatterm(xy_top_m(:,2)*180/pi,xy_top_m(:,1)*(180)/pi,30,'k','filled');
                textm(xy_top_m(:,2)*180/pi,xy_top_m(:,1)*(180)/pi,obj.verifyMarkers,...
                    'Color','k','FontSize',10,'EdgeColor','b','VerticalAlignment','top','LineStyle','none');
                title(['Sky View of Global Sky Model 2016 in ',obj.coorSys,' Coordinate System']);
            else

                axesm mollweid                     
                scatterm(obj.xy(:,2)*180/pi,obj.xy(:,1)*180/pi,12,gmap,'filled');
                %colorbar ('south'); colormap jet;caxis([4 8]);
                
                gridm;plabel;mlabel;                
                title(['Global Sky Model 2016 in ',obj.coorSys,' Coordinate System']);
                
                setm(gca,'GLineWidth',0.01,'Gcolor',[0 0 0],'GLineStyle','--',...
                    'MLineLocation',15,'MLabelParallel','equator',...
                    'FontColor',[0.4660 0.6740 0.1880],'FontSize',18,...
                    'LabelFormat','signed','PLabelMeridian','prime','Frame','off')
                colorbar ('southoutside');tightmap;axis off;
                if strcmp(obj.coorSys, 'GalLongLat')
                    mlabelzero22pi;
                end
            end
                
        end
  
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
        function [longlat_p,longlat_m]= rotate(obj,rot_x,rot_y,rot_z)
            gmap = obj.generated_map_data(:,1);
            gmap = log2(gmap);
            eq=obj.xyHorizon;
            eq_ph=eq(:,1);eq_th=eq(:,2);
        
            [x_cart,y_cart,z_cart] = sph2cart(obj.xy(:,1).',obj.xy(:,2).',1);% converting spherical to cartesian
            [x_cart_m,y_cart_m,z_cart_m] = sph2cart(obj.markers_xy(:,1).',obj.markers_xy(:,2).',1);% converting spherical to cartesian
            [x_cart_eq,y_cart_eq,z_cart_eq] = sph2cart(eq_ph,eq_th,1);
            Rx=rotx(-rot_x);Ry=roty(-rot_y);Rz=rotz(rot_z);
        
            %Rotation around X i.e. around axis out of the screen +ve
            %angles for clockwise and -ve angles for counter-clockwise
            %rotation
            new_top=Rx*[x_cart;y_cart;z_cart];%rotating in cartesian
            x1_cart=new_top(1,:);y1_cart=new_top(2,:);z1_cart=new_top(3,:);
            [xy_rx(:,1),xy_rx(:,2),~] = cart2sph(x1_cart,y1_cart,z1_cart);%converting back to spehrical for plotting
        
            new_top_m=Rx*[x_cart_m;y_cart_m;z_cart_m];%rotating in cartesian
            x1_cart_m=new_top_m(1,:);y1_cart_m=new_top_m(2,:);z1_cart_m=new_top_m(3,:);
            [markers_xy_rx(:,1),markers_xy_rx(:,2),~] = cart2sph(x1_cart_m,y1_cart_m,z1_cart_m);
        
            new_top_eq=Rx*[x_cart_eq,y_cart_eq,z_cart_eq].';
            x1_cart_eq=new_top_eq(1,:);y1_cart_eq=new_top_eq(2,:);z1_cart_eq=new_top_eq(3,:);
            [eq_ph1,eq_th1,~]=cart2sph(x1_cart_eq,y1_cart_eq,z1_cart_eq);
        
            %************************************************************************************************************************
        
            %Rotation around Y i.e. movement in N-S direction +ve angles
            %moves in North direction and -ve angles moves in
            %South direction
            if (rot_x==0)
                new_NS=Ry*[x_cart;y_cart;z_cart];
                new_NS_m=Ry*[x_cart_m;y_cart_m;z_cart_m];
                new_NS_eq=Ry*[x_cart_eq,y_cart_eq,z_cart_eq].';
            else
                new_NS=Ry*new_top;
                new_NS_m=Ry*new_top_m;
                new_NS_eq=Ry*new_top_eq;
            end
            x2_cart=new_NS(1,:);y2_cart=new_NS(2,:);z2_cart=new_NS(3,:);
            [xy_ry(:,1),xy_ry(:,2),~] = cart2sph(x2_cart,y2_cart,z2_cart);
        
            x2_cart_m=new_NS_m(1,:);y2_cart_m=new_NS_m(2,:);z2_cart_m=new_NS_m(3,:);
            [markers_xy_ry(:,1),markers_xy_ry(:,2),~] = cart2sph(x2_cart_m,y2_cart_m,z2_cart_m);
        
            x2_cart_eq=new_NS_eq(1,:);y2_cart_eq=new_NS_eq(2,:);z2_cart_eq=new_NS_eq(3,:);
            [eq_ph2,eq_th2,~]=cart2sph(x2_cart_eq,y2_cart_eq,z2_cart_eq);
            %************************************************************************************************************************
        
            %Rotation around Z i.e. movement in E-W direction +ve angles
            %moves in East direction and -ve angles moves in West
            %direction
            if(rot_x~=0 && rot_y==0)
                new_EW=Rz*new_top;
                new_EW_m=Rz*new_top_m;
                new_EW_eq=Rz*new_top_eq;
            elseif (rot_y~=0)
                new_EW=Rz*new_NS;
                new_EW_m=Rz*new_NS_m;
                new_EW_eq=Rz*new_NS_eq;
        
            elseif (rot_x==0 && rot_y==0)
                new_EW=Rz*[x_cart;y_cart;z_cart];
                new_EW_m=Rz*[x_cart_m;y_cart_m;z_cart_m];
                new_EW_eq=Rz*[x_cart_eq,y_cart_eq,z_cart_eq].';
            end
            x3_cart=new_EW(1,:);y3_cart=new_EW(2,:);z3_cart=new_EW(3,:);
            [xy_rz(:,1),xy_rz(:,2),~] = cart2sph(x3_cart,y3_cart,z3_cart);
        
            x3_cart_m=new_EW_m(1,:);y3_cart_m=new_EW_m(2,:);z3_cart_m=new_EW_m(3,:);
            [markers_xy_rz(:,1),markers_xy_rz(:,2),~] = cart2sph(x3_cart_m,y3_cart_m,z3_cart_m);
        
            x3_cart_eq=new_EW_eq(1,:);y3_cart_eq=new_EW_eq(2,:);z3_cart_eq=new_EW_eq(3,:);
            [eq_ph3,eq_th3,~]=cart2sph(x3_cart_eq,y3_cart_eq,z3_cart_eq);
            %************************************************************************************************************************
        
            if(rot_x~=0 && rot_y==0 && rot_z==0)
                longlat_p = [wrapToPi(xy_rx(:,1)),xy_rx(:,2)];
                longlat_m = [wrapToPi(markers_xy_rx(:,1)),markers_xy_rx(:,2)];
                eq_p=[eq_ph1.',eq_th1.'];
            elseif (rot_x==0 && rot_y~=0 && rot_z==0)||(rot_x~=0 && rot_y~=0 && rot_z==0)
                longlat_p = [wrapToPi(xy_ry(:,1)),xy_ry(:,2)];
                longlat_m = [wrapToPi(markers_xy_ry(:,1)),markers_xy_ry(:,2)];
                eq_p=[eq_ph2.',eq_th2.'];
            elseif(rot_x==0 && rot_y==0 && rot_z~=0)||...
                    (rot_x==0 && rot_y~=0 && rot_z~=0)||...
                    (rot_x~=0 && rot_y==0 && rot_z~=0)||...
                    (rot_x~=0 && rot_y~=0 && rot_z~=0)
                longlat_p = [wrapToPi(xy_rz(:,1)),xy_rz(:,2)];
                longlat_m = [wrapToPi(markers_xy_rz(:,1)),markers_xy_rz(:,2)];
                eq_p=[eq_ph3.',eq_th3.'];
            elseif (rot_x==0 && rot_y==0 && rot_z==0)
                longlat_p = [wrapToPi(obj.xy(:,1)),obj.xy(:,2)];
                longlat_m = [wrapToPi(obj.markers_xy(:,1)),obj.markers_xy(:,2)];
                eq_p=eq;
            end
        
        
            if strcmp(obj.projectionType,'m')
                figure

                axesm mollweid
                scatterm(longlat_p(:,2)*180/pi,longlat_p(:,1)*180/pi,10,gmap,'filled');%colormap jet;colorbar;caxis([4 8]);
                gridm;plabel;mlabel;
                setm(gca,'GLineWidth',0.01,'Gcolor',[0 0 0],'GLineStyle','--',...
                    'MLineLocation',15,'MLabelParallel','equator',...
                    'FontColor',[0.4660 0.6740 0.1880],'FontSize',18,...
                    'LabelFormat','signed','PLabelMeridian','prime','Frame','off');
                colorbar ('southoutside');tightmap;axis off;
                
                title({['Rotated Projection in ',obj.coorSys,' Coordinate System'];...
                    ['Rotaion East-West =',num2str(rot_z),...
                    ', North-South =',num2str(rot_y),', Top =',num2str(rot_x)]});
                if strcmp(obj.coorSys, 'GalLongLat')
                    mlabelzero22pi;
                end
        
                scatterm(eq_p(:,2)*180/pi,eq_p(:,1)*180/pi,2,'w','filled');
        
                scatterm(longlat_m(:,2)*180/pi,longlat_m(:,1)*180/pi,30,'k','filled');
                textm(longlat_m(:,2)*180/pi,longlat_m(:,1)*180/pi,obj.verifyMarkers,...
                    'Color','k','FontSize',10,'EdgeColor','b','VerticalAlignment','top','LineStyle','none');
            end
        end

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        function plotmarkers(obj)
            % plotmarkers plots the horizon 
            % plot approximate sun and moon positions 
            % on the current axis

            if strcmp(obj.projectionType,'m')
                scatterm(obj.xyHorizon(:,2)*180/pi,obj.xyHorizon(:,1)*180/pi,2,'w','filled');
                scatterm(obj.markers_xy(:,2)*180/pi,obj.markers_xy(:,1)*180/pi,30,'k','filled');
                textm(obj.markers_xy(:,2)*180/pi,obj.markers_xy(:,1)*180/pi,obj.verifyMarkers,...
                    'Color','k','FontSize',10,'EdgeColor','b','VerticalAlignment','top','LineStyle','none');  
            end
        end        
        
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@         
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
%%*************************************************************************************************%%
%wrap2pi>>wrapToPi
%PhTh2Dircos>>coo2cosined
%%*************************************************************************************************%%