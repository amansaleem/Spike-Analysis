%
% EXPTPICK (formerly SPIKES) defines a few environment variables and launches the Graphical User Interfaces 
%

% 2009-03-06 AZ
% Renamed spikes and split into function and script versions. This SCRIPT 
% version, called simply by its name:
% >> ExptPick
% is just an alias for:
global DEMO DIRS PICK

ExptPickFcn;

% This allows all the variables which later programs depend on to be left
% in the active workspace.