%% Script to set up the parameters object for the figure 1 experiment (d=4, c=0.5, various p, L).
% You can not only assign values here, but have differential processing
% (e.g. to do different things on your desktop or cluster).
% Required members are described as they appear below.
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

clear parameters; % In case parameters previously held the name of the parameters file

% You can use this boolean to flag different folders for local and cluster runs:
if exist('/home/jlia0386/') > 0
    isCluster = false;
else
    isCluster = true;
end

parameters.EZ = [@P1];


% - folder - directory where all of the files are to be stored.
if isCluster
    parameters.folder = './results/Stability/[@P1]';
else
    parameters.folder = './results/';
end

