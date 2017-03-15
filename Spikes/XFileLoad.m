function x = XFileLoad ( xfilename, quiet )
% XFileLoad loads an x file
%
% x = XFileLoad ( xfilename )
% 
% XFileLoad ( xfilename,'quiet') suppresses any text
% output (important e.g. when called from a Visual Basic program).
%
% needs to be extended so it reads more properties of the xfile
%
% 2006-08 Matteo Carandini
% 2007-xx MC changed the scanning because Matlab 7.4 was choking
% 2010-10 MC changed so it reads the parameter defaults
% 2011-02 MC added '.x' if needed, introduced SetDefaultDirs

if nargin<2
    quiet = 'loud';
end

if nargin<1
    error('Must specify the x file to load');
end

global DIRS

if ~strcmp(xfilename(end-1:end),'.x')
    xfilename = [xfilename '.x'];
end

x.name      = '';
x.parnames  = [];
x.pardefs   = [];
x.npars     = [];

if ~isfield(DIRS,'xfiles')
    SetDefaultDirs;
    % error('Need to specify DIRS.xfiles');
end

if isempty(DIRS.xfiles)
    xfiledir = fullfile(fileparts(DIRS.data),'xfiles');
else
    xfiledir = DIRS.xfiles;
end

if ~isdir(xfiledir)
    error('Problem with x file directory : it is not a directory');
end

fullxfilename = fullfile(xfiledir,xfilename);

if ~exist(fullxfilename,'file')
    errordlg(['Could not find xfile ' fullxfilename ],'Spikes', 'modal');
    return;
end

xfileptr = fopen(fullxfilename,'r');

if xfileptr < 0
    errordlg(['xfile ' fullxfilename ' exists but I could not open it'],'Spikes');
    return;
end

if ~strcmp(quiet,'quiet')
    disp([ 'Reading ' xfilename ]);
end

ss = fscanf(xfileptr,'%c');

controlMs  = findstr(ss,10);	% the control-M characters
endoflines = findstr(ss,13);	% the end-of-line characters

begoflines = 1+[0, endoflines(1:end-1), controlMs(1:end-1)];
% the good lines start with a number
goodlines = find( ss(begoflines)>='0' & ss(begoflines)<='9' );
if isempty(goodlines)
    error('Trouble reading the xfile');
end

starthere = begoflines(goodlines(1));
thisstring = ss(starthere:end);
[foo1 x.parnames x.pardefs x.pardefaults foo2 foo3] = strread(thisstring,'%d %s %q %d %s %s');

successflag = 1;

if ~successflag
    disp('something wrong reading the xfile');
    return;
end

fclose(xfileptr);

x.name  = xfilename;
x.npars = length(x.pardefs);
