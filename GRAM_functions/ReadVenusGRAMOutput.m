% Script to read VenusGRAM output files for processing in MATLAB.
% DELETE ALL LIST.txt FILES before running this
clear
clc

%% INPUTS -----
% Folder where output files reside
OutputDir = 'D:\Documents\GRAMsuite\Venus\results';

%% ------------------------------
outFiles = dir(OutputDir);
time = [];
h_km = [];
lat = [];
lon = [];
rho = [];
T = [];
u = [];
v = [];
Ls = [];
% CO2 & concentrations..

for i = 3:length(outFiles)
    
    outFileAdd = fullfile(OutputDir,outFiles(i).name);
    % using readtable() for convenience, textscan() can also be used for a
    % lower level handling of data (like in GPS-TEC project)
    fprintf('Reading file: %s \n',outFiles(i).name)

    output = readtable(outFileAdd);
    
    % process data here..
    time = [time;output.Time];
    h_km = [h_km;output.Height];
    lat = [lat;output.Lat];
    lon = [lon;output.LonE];
    rho = [rho;output.Denkgm3];
    T = [T;output.Temp];
    u = [u;output.EWind];
    v = [v;output.NWind];
    Ls = [Ls;output.Ls];
    
    % putting all together in one array for sorting tasks
    allOuts = [time h_km lat lon rho T u v Ls];
end

% useful sorting
% by height: h_10 = allOuts(allOuts(:,2)==10,:);

% interpolating data onto a finer grid/ onto 2d grid for contour plotting
% [X,Y] = ndgrid(-90:10:90,0:10:360);
% U = griddata(h0(:,3),h0(:,4),h0(:,7),X,Y);
% contourf(Y,X,U,'EdgeColor','none')

% % lat-height plot for u, at a given lon
 lon0 = allOuts(allOuts(:,4)== 180,:);
 [X,Y] = ndgrid(-90:1:90,0:1:250);
 U = griddata(lon0(:,3),lon0(:,2),lon0(:,7),X,Y);
 contourf(X,Y,U,20,'EdgeColor','none')
 colorbar
 
 % plot vertical temperature for some lats
 lat0 = allOuts(allOuts(:,3)==0,:);
 plot(lat0(:,6),lat0(:,2),'.')