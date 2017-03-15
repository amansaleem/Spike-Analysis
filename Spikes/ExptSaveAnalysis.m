function success = ExptSaveAnalysis(expt, chans, spikedir,datadir,overwriteflag)
% ExptSaveAnalysis saves the Chan info and the Unit info (calling UnitSave)
%
% success = ExptSaveAnalysis(expt,chans,spikedir,datadir)
%
% success = ExptSaveAnalysis(expt,chans,spikedir,datadir,'overwrite')
% does not ask if you want to overwrite existing files.
%

% 2000 Matteo Carandini created
% 2007-11 MC removed warning about NaN channel names
% 2007-11 MC forced deletion of previous files for relevant channels
% 2009-01 SK print information about the working directory in the command window before saving units
% 2009-02 AZ Allow archiving of previous spike files

if nargin < 5
   overwriteflag = 'ask';
end

success = 0;

nchans = length(chans);

% see if the animal directory exists
if ~exist(fullfile(spikedir,expt.animal),'dir')
   [status, msg] = mkdir(spikedir,expt.animal);
   if ~isempty(msg), errordlg(msg,'Spike Sorter','modal'); return; end
   % if a log file exists, copy it into the new dir
   logfilerelpath = fullfile(expt.animal,[expt.animal '.txt']);
   if exist(fullfile(datadir,logfilerelpath),'file')
      [status, msg] = copyfile(...
         fullfile(datadir,logfilerelpath),...
         fullfile(spikedir,logfilerelpath));
      if ~isempty(msg), errordlg(msg,'Spike Sorter','modal'); return; end
   end
end

% see if the series directory exists
if ~exist(fullfile(spikedir,expt.animal,num2str(expt.iseries)),'dir')
   [status, msg] = mkdir(fullfile(spikedir,expt.animal),num2str(expt.iseries));
   if ~isempty(msg), errordlg(msg,'Spike Sorter','modal'); return; end
end

% see if all the spikes have been given a proper name 
for ichan = 1:nchans
   cellids = chans(ichan).cellids(1:chans(ichan).nprots);
%    if any(isnan(cellids))
%       warndlg(['Channel ' num2str(ichan) ' has one or more prototypes that are not assigned a number and will not be saved.'],...
%          'Spike Sorter', 'modal');
%    end
   if length(cellids)~=length(unique(cellids))
      errordlg(['Channel ' num2str(ichan) ' has two or more prototypes that are assigned the same number! Not saving.'],...
         'Spike Sorter', 'modal');
      return
   end
end

% we add a suffix because Spike Sorter could be run on various types of
% data in a given experiment
switch chans(1).channame
    case 101
        DataType =  'Michigan';
    case 501
        DataType = 'CerebusTraces';
    otherwise
        DataType = 'Multispike'; % this is also for traces
end

fullfilename = ChanGetFileName(spikedir,expt.animal,expt.iseries,expt.iexp,DataType);
fprintf(['\nWorking directory: %s%c%s%c%d\n'], spikedir, filesep, expt.animal, filesep, expt.iseries);

if exist(fullfilename,'file')
    % this section is run only if discrim_pars_Michigan_expt_x.mat
    % exists...    
   if ~strcmp(overwriteflag,'overwrite')
      ButtonName=questdlg(['Spike files already exist for this experiment. ',...
          'Archive them to ''..' filesep 'archive',filesep,''', or overwrite them?'], ...
         'Spike Sorter', ...
         'Archive','Overwrite','Cancel', 'Archive');
   else
      ButtonName = 'Overwrite';
   end
   % AZ 2009-02: Allow archiving of previous spike files
   if strcmp(ButtonName, 'Archive')
       if ~isdir([spikedir filesep expt.animal filesep 'archive'])
           mkdir([spikedir filesep expt.animal filesep 'archive']);
       end
       old.fileinfo = dir(fullfilename);
       [path,name,ext] = fileparts(fullfilename);
       old.filename = [spikedir filesep expt.animal filesep 'archive' filesep name ...
           '-' datestr(old.fileinfo.date,'yyyymmddHHMM') ext];
       movefile(fullfilename,old.filename);
       fprintf(['   Previously sorted spike file\n   ''%s''\n   moved to\n   ''%s''\n'],...
           fullfilename,old.filename);
   elseif strcmp(ButtonName, 'Cancel')
       return;
   end
end

%% now go ahead and delete all previous files for those channels:

for ichan = 1:nchans
    % A hack to find the value of "pos":
    filename = UnitGetFilename( expt.animal, expt.iseries, expt.iexp, chans(ichan).channame, 1 );
    pos = findstr('_u001',filename);
    d = dir(fullfile(spikedir, [ filename(1:pos) '*' ]));
    if ~isempty(d)
        fprintf(1, 'Deleting files:\n');
        for ifile = 1:length(d),
            fprintf(1, '%s\n', d(ifile).name);
        end
    end
    delete(fullfile(spikedir, [ filename(1:pos) '*' ]));
end
   
%%

fprintf(1, 'Saving spike discrimination information for all the channels...');
save(fullfilename,'chans');
fprintf(1, 'done\n');

% copy the pfile to spikedir
pfilename = sprintf('%s_%d_%d.p',expt.animal, expt.iseries, expt.iexp);
pfileoldpath = fullfile( datadir,expt.animal,int2str(expt.iseries),int2str(expt.iexp),pfilename);
pfilenewpath = fullfile( spikedir,expt.animal,int2str(expt.iseries),pfilename);
[status, msg] = copyfile(pfileoldpath,pfilenewpath,'f'); % MC 2010-03-02 added 'f'
if ~isempty(msg), errordlg(msg,'ExptSaveAnalysis','modal'); return; end

% make a list of which cells were present
neighborhood = {};
for ichan = 1:nchans
   for iprot = 1:chans(ichan).nprots
      icell = chans(ichan).cellids(iprot); % the cell id associated with the prototype
      if ~isnan(icell)
         strcell = num2str(1000+icell); strcell = strcell(2:4);
         neighborhood{end+1} = [ expt.animal '_' num2str(expt.channames(ichan)) '_' strcell];
      end
   end
end

%%  save the data Unit by Unit
for ichan = 1:nchans
   for iprot = 1:chans(ichan).nprots
      icell = chans(ichan).cellids(iprot); % the cell id associated with the prototype
      if ~isnan(icell)
         unit = UnitMake(expt,chans(ichan),iprot,neighborhood);
         UnitSave(unit,spikedir);
      end
   end
end

success = 1;
