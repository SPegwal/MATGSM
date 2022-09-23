function [lat,long,alt] = setLocation(label)
%setLocation sets the location of the observer
%input: label from the list 
%output: location coordinates of that label
% Tested : matlab 9.10
%     By : Saurabh Pegwal                 Jan 2022
% Example: [a,b,c]=setLocation(label)
Location_label = { 'Null Island', 'REACH', 'SKA', 'HERA','PAPER' , 'GBT' , 'GMRT' , 'FAST', 'ASKAP', 'ATCA', 'MWA', 'PaST', 'LOFAR', 'CHIME'};       

L_coord=[0,0,0; -30.838750, 21.374920, 1150.0000; -30.72113, 21.411128,1080; -30.721339829017378, 21.4282424879626, 1080;...
    -30.722400, 21.427800, 1080;  38.43302873568111, -79.83982982125139,807.43;...
    19.091220708385247, 74.0505174019285, 680; 25.652889, 106.856778, 5080;...
    -26.69699102820478, 116.63105557519636, 380; -30.312987696087397, 149.56440104951383, 180;...
    -26.702994018015414, 116.67038222390741, 1060.00; 42.924200, 86.716000, 2500;...
    52.90156227803543, 6.849166859086357, 140; 49.320833, -119.623611, 545;];

assert(ismember(label,Location_label),'Unkown location. Allowable Location label are: Null Island, REACH, SKA, HERA, PAPER , GBT, GMRT, FAST, ASKAP, ATCA, MWA, PaST, LOFAR, CHIME');

for i=1:length(L_coord(:,1))
    if strcmp(Location_label{i}, label)
        lat=L_coord(i,1);long=L_coord(i,2);alt=L_coord(i,3);
    end
end
end