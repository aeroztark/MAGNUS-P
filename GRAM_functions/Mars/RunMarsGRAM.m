% This script runs MarsGRAM program for a large number of inputs files provided
% in a directory. The input files are themselves created using
% CreateMarsGRAMInput.m script from an iteration schedule listed in a csv file. 

% Just need to provide the MarsGRAM exe location and input files' folder

clear
clc
%% INPUTS -----
% Folder where MarsGRAM executable ('marsgram_M10.exe') resides
MarsGRAMdir = 'D:\Documents\GRAMsuite\Mars\Release1.0_Nov10\Executables';

% Folder where input files reside
InputsDir = 'D:\Documents\GRAMsuite\Mars\Release1.0_Nov10\RunResults\MCS_demo\Inputs';
% Desired folders for output files are already specified in the input files. 

%%----------
inFiles = dir(InputsDir);

for i = 3:length(inFiles) % skip '.' and '..' files
    
    inputFileAdd = fullfile(InputsDir,inFiles(i).name);  % full address  
    
    % Running MarsGRAM
    cd(MarsGRAMdir); % enter the dir
    system(sprintf('echo %s | marsgram_M10.exe',inputFileAdd)); % execute the command with input file
    fprintf('\nsuccess!\n')
end

fprintf('\nAll MarsGRAM outputs successfully generated!\n')