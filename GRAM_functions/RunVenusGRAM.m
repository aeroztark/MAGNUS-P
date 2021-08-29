% This script runs VenusGRAM program for a large number of inputs files provided
% in a directory. The input files are themselves created using
% CreateVenusGRAMInput.m script from an iteration schedule listed in a csv file. 

% Just need to provide the VenusGRAM exe location and input files' folder

clear
clc
%% INPUTS -----
% Folder where VenusGRAM executable ('a.exe') resides
VenusGRAMdir = 'D:\Documents\GRAMsuite\Venus';

% Folder where input files reside
InputsDir = 'D:\Documents\GRAMsuite\Venus\inputs';
% Desired folders for output files are already specified in the input files. 

%%----------
inFiles = dir(InputsDir);

for i = 3:length(inFiles) % skip '.' and '..' files
    
    inputFileAdd = fullfile(InputsDir,inFiles(i).name);  % full address  
    
    % Running VenusGRAM
    cd(VenusGRAMdir); % enter the dir
    system(sprintf('echo %s | a.exe',inputFileAdd)); % execute the command with input file
    fprintf('\nsuccess!\n')
end

fprintf('\nAll VenusGRAM outputs successfully generated!\n')