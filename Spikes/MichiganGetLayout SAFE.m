function [arrayLayout, arrayNo, shanks, tet] = MichiganGetLayout(animal, iseries, nchan, arrayStyleNumber)
% MichiganGetLayout - gives the physical arrangement of Michigan array channels in an NxN matrix
%
% [arrayLayout, shanks, tet] = MichiganGetLayout operates on
% current choice of animal and series (global PICK)
%
% [arrayLayout, shanks, tet] = MichiganGetLayout(animal, iseries) 
% lets you specify animal and series
%
% [arrayLayout, shanks, tet] = MichiganGetLayout(animal, iseries, arrayStyleNumber)
%
% It returns a matrix with the layout of the channels, and the number of
% shanks, whether or not in tetrode configuration and whether they are 
% tetrodes if tet = 1
%
% If the argument arraySerialNumber is not provided, the function might
% prompt the user to specify the serial number of the array that was used in a
% particular series.
%
% WE SHOULD REPLACE THIS ASAP WITH THE APPROPRIATE CALLS TO OBJECT Probe
%
% 2010-07    AS created (adapted from UtahGetLayout)
% 2010-12-08 AS added  option 'NOT A MICHIGAN' and allowed for extra channels.
% 2011-05-18 MC added option to load from PICK, removed tetLayout as output
% 2011-06-30 ND added options for 64 channel Buzsaki and Polytrode1A
% 2011-07-06 ND added a check for nchans if not given as input argument
% then assumes nchan should just be # of electrode --- not yet implemented

global DIRS PICK

nos = {'***NOT A MICHIGAN***',...
    '16-channel: A1x16','16-channel: A4x4','16-channel: A4x1tet','16-channel: A2x2tet',...
    '32-channel: A4x2tet','32-channel: A4x8','32-channel: A1x32','32-channel: A8x4','32-channel: A2x16',...
    '64-channel: A8x8','64-channel: Buzsaki 8x8','64-channel: Polytrode1A', '64-channel: A4x16','64-channel: A4x16 flipped banks',...
    };

if nargin < 2
    animal = PICK.animal;
    iseries = PICK.iseries;
end

if nargin <= 3 
    arrayNo = [];
end

if nargin == 4 && ~ismember(arrayStyleNumber, nos)
    error('Invalid Array Style Number');
elseif nargin == 4 && ismember(arrayStyleNumber, nos)
    arrayNo = arrayStyleNumber;
end

dataDir = strcat(DIRS.Cerebus, '\', animal, '\');
matFileName = strcat(dataDir, sprintf('layoutSeries%03d.mat', iseries));

% check if we have access to the machine
tmp = ls(DIRS.Cerebus);
if isempty(tmp)
    error('Cannot access %s', DIRS.Cerebus);
end

% check if a file with the layout already exists
if exist(matFileName, 'file')
    load(matFileName);
    if ~exist('arrayNo', 'var')
        arrayNo = [];
    else
        return;
    end
end

% if the *.mat file is not found in the data directory and the serial number is not provided as an argument,
% get the serial number from user

if isempty(arrayNo)
    [index, ok] = listdlg('PromptString', 'Select a style number for the Michigan array:',...
        'ListSize', [250 180], 'SelectionMode','single', 'ListString', nos);
    if ~ok
        error('No Array Serial Number specified.');
    else
        arrayNo = nos{index};
    end
end

% determine arrangement of banks
switch arrayNo
    
    case '***NOT A MICHIGAN***'
        A = (1:nchan)';
        
        shanks = 1;
        tet = 0;
        arrayLayout = A;
        
    case '16-channel: A1x16'
        
        A = [9 8 10 7 13 4 12 5 15 2 16 1 14 3 11 6]';
        
        shanks = 1;
        tet = 0;
        arrayLayout = A;
        if ~exist('nchan','var'), nchan=16; end
        
    case '16-channel: A4x4'
        
        A = [1 3 2 6]';
        B = [7 4 8 5]';
        C = [13 10 12 9]';
        D = [14 16 11 15]';
        
        shanks = 4;
        tet = 0;
        arrayLayout = [A B C D];
        if ~exist('nchan','var'), nchan=16; end
        
    case '16-channel: A4x1tet'
        
        %         Order is top-down and then left-right (in each shank)
        A = [3 6 2 1]';
        B = [4 5 8 7]';
        C = [10 9 12 13]';
        D = [16 15 11 14]';
        
        shanks = 4;
        tet = 1;
        arrayLayout = [A B C D];
        if ~exist('nchan','var'), nchan=16; end
        tetLayout = [1 2 3 4; 5 6 7 8;9 10 11 12; 13 14 15 16]'; %[1 2 3 4];
        
    case '16-channel: A2x2tet'
        
        %         Order is top-down and then left-right (in each shank)
        A = [2 3 7 5]';
        B = [1 6 8 4]';
        C = [12 10 14 15]';
        D = [13 9 11 16]';
        
        shanks = 2;
        tet = 1;
        arrayLayout = [A B C D];
        if ~exist('nchan','var'), nchan=16; end
        tetLayout = [1 2 3 4; 5 6 7 8;9 10 11 12; 13 14 15 16]'; %[[1 2]' [3 4]'];
        
    case '32-channel: A4x2tet'
        
        %         Order is top-down and then left-right (in each shank)
        A = [4 2 7 5]';
        B = [3 1 8 6]';
        C = [12 10 15 13]';
        D = [11 9 16 14]';
        E = [20 18 23 21]';
        F = [19 17 24 22]';
        G = [28 26 31 29]';
        H = [27 25 32 30]';
        
        shanks = 4;
        tet = 1;
        arrayLayout = [A B C D E F G H];
        if ~exist('nchan','var'), nchan=32; end
        tetLayout = [1:4; 5:8; 9:12; 13:16; 17:20; 21:24; 25:28; 29:32]'; %[[1 2]' [3 4]' [5 6]' [7 8]'];
        
    case '32-channel: A4x8'
        
        %         Order is top-down and then left-right (in each shank)
        A = [5 4 6 3 7 2 8 1]';
        B = [13 12 14 11 15 10 16 9]';
        C = [21 20 22 19 23 18 24 17]';
        D = [29 28 30 27 31 26 32 25]';
        
        shanks = 4;
        tet = 0;
        arrayLayout = [A B C D];
        if ~exist('nchan','var'), nchan=32; end
        
    case '32-channel: A2x16'
        
        %         Order is top-down and then left-right (in each shank)
        A = [9 8 10 7 11 6 12 5 13 4 14 3 15 2 16 1]';
        B = [25 24 26 23 27 22 28 21 29 20 30 19 31 18 32 17]';
           
        shanks = 2;
        tet = 0;
        arrayLayout = [A B];
        if ~exist('nchan','var'), nchan=32; end
        
    case '32-channel: A1x32'
        
        %         Order is top-down and then left-right (in each shank)
        A = [17 16 18 15 19 14 20 13 21 12 22 11 23 10 24 9 25 8 26 7 27 6 28 5 ...
            29 4 30 3 31 2 32 1]';
        
        shanks = 1;
        tet = 0;
        arrayLayout = A;
        if ~exist('nchan','var'), nchan=32; end
        
    case '32-channel: A8x4'
        
        %         Order is top-down and then left-right (in each shank)
        A = [3 2 4 1]';
        B = [7 6 8 5]';
        C = [11 10 12 9]';
        D = [15 14 16 13]';
        E = [19 18 20 17]';
        F = [23 22 24 21]';
        G = [27 26 28 25]';
        H = [31 30 32 29]';
        
        shanks = 8;
        tet = 0;
        arrayLayout = [A B C D E F G H];
        if ~exist('nchan','var'), nchan=32; end
        %         tetLayout = [[1 2]' [3 4]' [5 6]' [7 8]'];
        
    case '64-channel: A8x8'
        
        %         Order is top-down, shanks go left-right
        A = [28 23 21 29 26 25 19 27]';
        B = [20 15 13 22 18 17 11 24]';
        C = [12 7 3 14 10 9 16 1]';
        D = [5 4 30 2 31 6 32 8]';
        E = [61 60 63 35 59 34 57 33]';
        F = [58 53 51 62 56 55 49 64]';
        G = [50 45 43 52 48 47 41 54]';
        H = [42 37 36 44 40 39 38 46]';
        
        shanks = 8;
        tet = 0;
        arrayLayout = [A B C D E F G H];
        if ~exist('nchan','var'), nchan=64; end

    case '64-channel: Buzsaki 8x8'
        
        %         Order is top-down the left edge then back up the right
        %         edge, shanks go left-right
        A = [27 25 29 23 28 21 26 19]'; % (1:8)';
        B = [24 17 22 15 20 13 18 11]'; % (9:16)';
        C = [16 9 14 7 12 3 10 1]'; % (17:24)';
        D = [8 6 2 4 5 30 31 32]'; % (25:32)';
        E = [33 34 35 60 61 63 59 57]'; % (33:40)';
        F = [64 55 62 53 58 51 56 49]'; % (41:48)';
        G = [54 47 52 45 50 43 48 41]'; % (49:56)';
        H = [46 39 44 37 42 36 40 38]'; % (57:64)';
        
        shanks = 8;
        tet = 0;
        arrayLayout = [A B C D E F G H];
        if ~exist('nchan','var'), nchan=64; end

    case '64-channel: Polytrode1A'
        
        %         Order is top-down and then left-right Right now calling
        %         it three shank even though it is one shank with three
        %         columns of recording sites.
        A = [11 12 10 13 9 14 8 15 7 16 17 18 6 5 4 3 2 1]';
        B = [27 26 28 25 29 24 30 23 31 22 32 21 33 20 34 19 35 36]';
        C = [43 42 44 41 45 40 46 39 47 38 48 37 49 50 51 52 53 54]';
        
        shanks = 3;
        tet = 0;
        arrayLayout = [A B C];
        if ~exist('nchan','var'), nchan=54; end
        
    case '64-channel: A4x16'
        %         Order is top-down and then left-right (in each shank)
        A = [24 19 17 26 22 21 15 28 20 23 13 29 18 25 11 27]';
        B = [8 1 6 10 2 3 4 12 5 7 30 14 31 9 32 16]';
        C = [64 57 55 59 62 63 53 61 58 60 51 35 56 34 49 33]';
        D = [46 41 39 48 44 43 37 50 42 45 36 52 40 47 38 54]';
        
        shanks = 4;
        tet = 0;
        arrayLayout = [A B C D];
        if ~exist('nchan','var'), nchan=64; end
        
    case '64-channel: A4x16 flipped banks'
       %         Order is top-down and then left-right (in each shank)
        A = [9 3 1 13 6 5 27 8 4 7 26 12 2 11 25 15]';
        B = [28 17 22 32 18 19 20 29 21 23 14 30 16 24 10 31]';
        C = [48 37 33 43 46 47 36 45 42 44 35 51 41 49 34 55]';
        D = [62 56 52 64 60 59 57 38 58 61 53 39 54 63 50 40]';
        
        shanks = 4;
        tet = 0;
        arrayLayout = [A B C D];
        if ~exist('nchan','var'), nchan=64; end        
        
    otherwise
        error('Unknown Array Style: %s.', arrayNo);
end

arrayLayout = arrayLayout(:);

if length(arrayLayout) <= nchan
    arrayLayout = [arrayLayout' (length(arrayLayout)+1):nchan]';
end
% make sure target directory exists (this is helpful for online analyses,
% in which case the data might not have been copied to zserver yet)
if ~exist(dataDir, 'dir')
    mkdir(dataDir);
end

if tet
    save(matFileName, 'arrayLayout', 'arrayNo', 'shanks', 'tet', 'tetLayout');
else
    save(matFileName, 'arrayLayout', 'arrayNo', 'shanks', 'tet');
end