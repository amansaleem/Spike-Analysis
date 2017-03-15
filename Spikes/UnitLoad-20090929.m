function units = UnitLoad( spikedir, animal, iseries, iexp, ichan, icell, batchflag, fixtimingflag )
% UnitLoad loads a Unit data structure
%
% units = UnitLoad( spikedir, animal, iseries, iexp ) returns an array of units
%
% units = UnitLoad( ExptTag ) is a shortcut for the same call (ExptTag
% should have fields animal, iseries, iexp)
%
% unit = UnitLoad( spikedir, animal, iseries, iexp, ichan, icell )
% lets you specify which ichan and which icell you are interested in
% (multiple ichan values accepted as a vector).
%
% unit = UnitLoad( spikedir, animal, iseries, iexp, [], icell )
% returns all units with any ichan, for a specified icell
% 
% units = UnitLoad( spikedir, animal, iseries, iexp, 'batch' ) and
% unit = UnitLoad( spikedir, animal, iseries, iexp, ichan, icell, 'batch' )
% does not bug you with error messages in modal windows.
%
% unit = UnitLoad( spikedir, animal, iseries, iexp, ichan, icell, 'batch', 'fix' )
% fixes timing for buggy xfiles only when it actually finds a problem
% [DEFAULT: does not fix]
% unit = UnitLoad( spikedir, animal, iseries, iexp, ichan, icell, 'batch', 'fixalways' )
% fixes timing for all buggy xfiles, independent on whether it finds a problem or not
%
% part of Spikes
%
% 2000-12 MC created
% 2003-02 MC fixed 'nobatch' bug, fixed removal of old stim duration fields
% 2003-03 VM made it fit for traces. When there is no field 'unit.spiketimes', all the
%            UnitGet and UnitPlot funtions look for a field 'unit.traces'.
% 2003-03 VB batchflag was sometimes not initialized. corrected.
% 2003-04 VM added fixtimingflag
% 2003-05 MC VM added a field "nrepeats"
% 2003-08 VB corrected small bug by adding remove flag.
% 2008-03 LB MC introduced call to UnitCreate to correct for missing fields
%   and added a few checks to avoid crashes for corrupted data converted from *.nev files
% 2009-05 MC changed the error message when the unit does not exist
% 2009-09 AZ can now set ichan to [] to get all units for any ichan, given
%   a specified icell.  can also set ichan to multiple values

global DIRS

if nargin == 1
    ExptTag = spikedir;
    units = UnitLoad(DIRS.spikes, ExptTag.animal, ExptTag.iseries, ExptTag.iexp);
    return
end

%%

if nargin < 8
   fixtimingflag = 'nofix';
end

if nargin<7 || isempty(batchflag)
   batchflag = 'nobatch';
end

if nargin==5
   batchflag = ichan;
   clear ichan;
end

if nargin >= 6
   if isempty(ichan)
      clear ichan;
   elseif isempty(icell)
      clear icell;
   end
end

if nargin<4
   error('You must at least specify the spikedir, the animal, the series, and the experiment');
end

if exist('ichan','var') && length(ichan)==1
   
   filename = UnitGetFilename( animal, iseries, iexp, ichan, icell );
   
   fullfilename = fullfile(spikedir,filename);
   
   if ~exist(fullfilename,'file')
      errmsg = ['Cannot find a unit with filename ' fullfilename ')'];
      if strcmp(batchflag,'nobatch')
         errordlg(errmsg,'UnitLoad', 'modal');
      else
         disp(errmsg);
      end
      units = [];
      return;
   end
   
   load(fullfilename);
   units = unit;
else
   
   % find all the units saved for this experiment
   if exist('icell','var') && ~isempty(icell)
      dd = dir(fullfile(spikedir,UnitGetFilename( animal, iseries, iexp, '*', icell )));
      if exist('ichan','var') && length(ichan) > 1
         % found_ichans = ichans of files found from given ichan vector, icell
         found_ichans = vertcat(dd.name);
         found_ichans = found_ichans(:,18:end-9);
         found_ichans = mat2cell(  found_ichans,ones(size(found_ichans,1),1));
         found_ichans = str2double(found_ichans);
         [junk,found_ichansIX] = intersect(found_ichans,ichan);
         dd = dd(found_ichansIX);
      end
   else
      dd = dir(fullfile(spikedir,UnitGetFilename( animal, iseries, iexp, '*', '*' )));
   end
   nunits = length(dd);
   if nunits<1
      errmsg = ['Could not find any saved units for animal ' animal ' series ' num2str(iseries) ' exp ' num2str(iexp) ];
      if strcmp(batchflag,'nobatch')
         errordlg(errmsg,'Spike Sorter', 'modal');
      else
         disp(errmsg);
      end
      units = [];
      return;
   end
   
   units = repmat(UnitCreate,[1,nunits]);
   for iunit = 1:nunits
      filename = dd(iunit).name;
      [yeah1, yeah2, chan, unit ] = sscanf(filename,[animal '_s%2d_e%2d_c%1d_u%3d.mat'],1);
      pp = findstr(filename,'_');
      chan = sscanf(filename(pp(3)+2:pp(4)-1),'%d');
      unit = sscanf(filename(pp(4)+2:end),'%d');
      % the following line can give an error if file is not there...
      load(fullfile(spikedir,UnitGetFilename( animal, iseries, iexp, chan, unit )));
      unit = UnitCreate(unit);
      units(iunit) = unit;
   end
   
end

if isfield(unit(1),'spiketimes') && ~isempty(unit(1).spiketimes)
    nstim = size(units(1).spiketimes,1);
    nrpts = size(units(1).spiketimes,2);
elseif isfield(unit(1),'traces') && ~isempty(unit(1).traces)
    nstim = size(units(1).traces,1);
    nrpts = size(units(1).traces,2);
else 
    expname = sprintf('%s-%d-%d',animal,iseries,iexp);
    errmsg = ['units in ' expname ' contain neither spiketimes nor traces, couldn`t determine nstim and nrpts'];
    if strcmp(batchflag,'nobatch')
        errordlg(errmsg,'Spike Sorter', 'modal');
    else
        disp(errmsg);
    end
end

nunits = length(units);

%--------------------------------------------------------------------------------

% Let's be optimistic at first
for iunit = 1:nunits
   units(iunit).trustedtiming = 1;
   units(iunit).correctedspiketimes = 0;
end

% TO ENSURE BACKWARD COMPATIBILITY
% check if you are loading one of the old structs
flagremovestimdur = 0;
if isfield(units(iunit),'stimdurmatrix')
	removeflags = zeros(1,length(units));  
   for iunit = 1:length(units) 
      if all(size(units(iunit).stimdurmatrix)==[nstim nrpts] )  
         units(iunit).stimdurs = units(iunit).stimdurmatrix;
         removeflags(iunit)=1;
      end
   end
   if all(removeflags);      units = rmfield(units,'stimdurmatrix');end;
   % or perhaps an even older one
elseif isfield(units,'stimdur') % it is one of the really old structs
   disp('Warning: You are loading an old data set');
   disp('Stimulus duration will be assumed to be the same across repeats');
   for iunit = 1:length(units) 
      units(iunit).stimdurs = units(iunit).stimdur*ones(1,nrpts);
      units(iunit).trustedtiming = 0;
   end
   units = rmfield(units,'stimdur');
end

% HACK for old files in which stimdurs is not saved
if (~isfield(units,'stimdurs') || isempty(units(1).stimdurs)) && isempty(units(1).source)
   disp('Warning: You are loading an old data set');
   disp('Stimulus duration was not saved.');
   disp('---> As a hack I will assume it is the time of the last spike.');
   
   stimdurs = ones(nstim,nrpts); % hack: first set them all to 1 second
   for istim = 1:nstim
      allspiketimes = [];
      for iunit = 1:nunits
         allspiketimes = [allspiketimes [units(iunit).spiketimes{istim,:}] ];
      end
      if ~isempty(allspiketimes)
         stimdurs(istim,:) = max(allspiketimes); % hack, set them to the last spike
      end
   end
   
   stimdurs(:) = max(stimdurs(:)); % SUPER HACK!

   for iunit = 1:length(units)
      units(iunit).stimdurs = stimdurs;
   end
end

%-----------------------------------------------------------------------------------
%						experiments for which we should not trust the timing
%-----------------------------------------------------------------------------------


% FOR THE FIRST 3 CAT EXPERIMENTS, DATA COLLECTION CONTINUED WELL AFTER STIMULUS END
buggyanimals = {'CAT001', 'CAT002', 'CAT003',};
if any(strmatch(upper(unit.animal),buggyanimals,'exact'))
   disp('---------- WARNING. For cats 1, 2 and 3, data acquisition went on beyond stimulus end.');
   % load the protocol
   protocol = ProtocolLoad(unit.animal,unit.iseries,unit.iexp);
   if ~isempty(protocol) & ~isempty(protocol.pfiledurs)
      disp('---------- (truncating to hypothetical stimulus duration)');
      for istim = 1:nstim
         pfiledur = protocol.pfiledurs(istim);
         for irpt = 1:nrpts
            for iunit = 1:nunits
               units(iunit).trustedtiming = 0;
               units(iunit).correctedspiketimes = 1;
               units(iunit).stimdurs(istim,irpt) = pfiledur;
               spiketimes = units(iunit).spiketimes{istim,irpt};
               if ~isempty(spiketimes)
                  spiketimes = spiketimes(find(spiketimes<pfiledur));
                  units(iunit).spiketimes{istim,irpt} = spiketimes;
               end
            end
         end
      end
   end
end

% FOR A LARGE PART OF CATZ027, PHOTODIODE WAS ALWAYS ON
% AND DATA COLLECTION STARTED WAY BEFORE STIMULUS (AND WENT ON A LITTLE TOO LONG AS WELL)
% (actually, we should rescue those expts in which timing was ok)
if strmatch(upper(unit.animal),'CATZ027','exact')
   disp('---------- WARNING. For exp CATZ027, data acquisition was badly timed. Do not trust timing data.');
   % load the protocol
   protocol = ProtocolLoad(unit.animal,unit.iseries,unit.iexp);
   if ~isempty(protocol) & ~isempty(protocol.pfiledurs)
      disp('---------- (truncating to hypothetical stimulus duration)');
      for istim = 1:nstim
         pfiledur = protocol.pfiledurs(istim);
         for irpt = 1:nrpts
            for iunit = 1:nunits
               units(iunit).trustedtiming = 0;
               units(iunit).correctedspiketimes = 1;
               recdur = units(iunit).stimdurs(istim,irpt);
               extradur = recdur - pfiledur;
               spiketimes = units(iunit).spiketimes{istim,irpt};
               if ~isempty(spiketimes)
                  spiketimes = spiketimes - extradur;
                  goodspikes = find(spiketimes>0 & spiketimes<pfiledur);
                  units(iunit).spiketimes{istim,irpt} = spiketimes(goodspikes);
               end
               units(iunit).stimdurs(istim,irpt) = pfiledur;
            end
         end
      end
   end
end

% SOME XFILES CAUSE A TIMING PROBLEM BECAUSE OF A BUG IN THE SIZE OF THE FRAMES
% IF fixtimingflag == 'fix' THE NEXT LINES CORRECT THE SPIKETIMES WHEN POSSIBLE
if strcmp(fixtimingflag(1:3),'fix') & isfield(unit,'spiketimes') & ~isempty(unit.spiketimes) 
   % Load the protocol and myscreen
   protocol = ProtocolLoad(units(1).animal,units(1).iseries,units(1).iexp);
   myscreen = ScreenLogLoad(units(1).animal,units(1).iseries,units(1).iexp);
   % Check whether the one can trust the timing
   framesequence = [];
   if ~isempty(protocol)
      [status,framesequence] = ProtocolTimingIsBuggy(protocol,myscreen);
      % Check whether one can trust the timing
      if status
         for iunit = 1:nunits
            units(iunit).trustedtiming = 0;
            units(iunit).framesequence = framesequence;
            %disp('----> WARNING !!!');
            disp(sprintf('----> Do not trust timing for animal %s series %i exp %i',...
               units(iunit).animal,units(iunit).iseries,units(iunit).iexp));
         end
      end
   end
   % Correct the spiketimes if possible and if needed
   if status == 1 & ~isempty(framesequence) & ~isempty(myscreen)
      for iunit = 1:nunits
         protdur = [];
         extradur = [];
         for istim = 1:nstim
            nframes = length(framesequence{istim});
            protdur(istim) = nframes/myscreen.FrameRate;
            meandur = mean(units(iunit).stimdurs(istim,:));
            extradur(istim) = meandur - protdur(istim);
         end
         
         if any(extradur > 1/myscreen.FrameRate) | strcmp(fixtimingflag,'fixalways')
            units(iunit).correctedspiketimes = 1;
            %buggystimuli = find(extradur > 1/myscreen.FrameRate);
            %disp('Buggy stimuli:');
            %if isempty(buggystimuli)
            %   disp('none');
            %else
            %   disp(sprintf('%i ',buggystimuli));
            %end
            disp('----> Adjusted spiketimes for all stimuli');
            %disp('----> Might have eliminated early spikes');
            for istim = 1:nstim
               % Shift all repeats by the same amount
               for irpt = 1:nrpts
                  recdur = units(iunit).stimdurs(istim,irpt);
                  spiketimes = units(iunit).spiketimes{istim,irpt};
                  if ~isempty(spiketimes)
                     spiketimes = spiketimes - extradur(istim);
                     goodspikes = find(spiketimes>0 & spiketimes<protdur(istim));
                     units(iunit).spiketimes{istim,irpt} = spiketimes(goodspikes);
                  end
                  units(iunit).stimdurs(istim,irpt) = protdur(istim);
               end
            end
         end
      end
   end
end



%-----------------------------------------------------------------------------------
%								add an "id" field to uniquely identify units
%-----------------------------------------------------------------------------------

for iunit = 1:nunits
   units(iunit).id = sprintf('%s.%d.%d',units(iunit).animal,units(iunit).ichan,units(iunit).icell);
end

%-----------------------------------------------------------------------------------
%								add an "nrepeats" field, pretty useful
%-----------------------------------------------------------------------------------
if exist('nrpts', 'var') % some corrupted experiments converted from *.nev don't have nrpts
    for iunit = 1:nunits
        units(iunit).nrepeats = nrpts;
    end
end

