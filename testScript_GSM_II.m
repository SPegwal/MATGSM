close all
clc

logged = true;

fPlot = 408;

G0 = GlobalSkyModel2016;
G0 = G0.generate(fPlot);
G1 = G0.underSample(3);

G1 = G1.setTime(datetime(2021,09,20,12,0,0));

% labels allowed for setLocation( 'Null Island', 
% 'REACH', 'SKA', 'HERA','PAPER' , 'GBT' , 'GMRT' , 
% 'FAST', 'ASKAP', 'ATCA', 'MWA', 'PaST', 'LOFAR', 'CHIME');

[a,b,c]=setLocation('REACH');
G1.location=[a,b,c];


figure,
G1.plotProj(1,true)
G1.plotmarkers

G2 = G1.setCoorSys('RAdec');
figure,
G2.plotProj(1,true)
G2.plotmarkers


G3 = G1.setCoorSys('Horiz');
figure,
G3.plotProj(1,true)
G3.plotmarkers


G4 = G3;
G4.projectionType = 'top';
figure,G4.plotProj(1,true)
G4.plotmarkers

% % for rotation
% % example: G1.rotate(x,y,z);
% % x= rotation along the axis coming out of the screen
% (+ve clockwise/-ve anti-clockwise)
% % y= rotation along N-S axis (+ve towards North/-ve South)
% % z= rotation along E-W axis (+ve towards East/-ve West)
% % Euler roation implemented so avoid using more than 2 axis rotation
% possibility of GIMBAL lock situation