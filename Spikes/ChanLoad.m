function chans = ChanLoad(spikedir,animal,iseries,iexp,DataType)
% ChanLoad loads a channel data structure
%
%	chans = ChanLoad(spikedir,animal,iseries,iexp)
%
% chans = ChanLoad(spikedir,animal,iseries,iexp,DataType)
% lets you specify the type of data, which can be 'Michigan', or
% 'Multispike' (DEFAULT).
%
% 2000 Matteo Carandini
% 2002-01 MC added extension of tt
% 2007-09 MC added 5th argument

if nargin<5
    DataType = 'Multispike';
end

chanfilename = ChanGetFileName(spikedir,animal,iseries,iexp,DataType);
if ~ exist(chanfilename,'file')
    warndlg('No file exists -- initializing channels');
    chans = ChanInitialize;
    return
end

load(chanfilename);

if exist('CHANS','var') == 1 % you are loading one of the old files
   chans = CHANS;
   for ichan = 1:length(CHANS)
      % BRUTALLY MAKE UP THE FIELDS THAT WERE NOT DEFINED:
      chans(ichan).iexptchan = ichan;
      chans(ichan).threshsign = 1; 	
      chans(ichan).mindur = 2; 		
   end
elseif exist('chans','var') == 1 % you are loading one of the new files
   % do nothing
else
   disp(['Could not load ' chanfilename]);
end

% if the channame field is not set: 
if ~isfield(chans,'channame')
   for ichan = 1:length(chans)
      chans(ichan).channame = []; % it really should be obtained from the channames field of Expt.
   end
end

if ~isfield(chans,'filtercoeffs')
   for ichan = 1:length(chans)
      chans(ichan).filtercoeffs = []; 
   end
end

% THIS IS UNBELIEVABLY SLOW, SO I COMMENTED IT OUT
% remove the field vv from spikes
% for ichan = 1:length(chans)
%   chans(ichan).spikes = rmfield(chans(ichan).spikes, 'vv');
% end

% ----- extend the duration of tt and of prototypes -----

tlim = [-1.50 4.00];	

for ichan = 1:length(chans)
   if chans(ichan).tt(1) > -1
      tt = chans(ichan).tt;
      dt = mean(diff(tt));
      pad_l = sort((tt(1)-dt):-dt:tlim(1));	% padding on the left
      pad_r = (tt(end)+dt):dt:tlim(2);			% padding on the right
      newtt = sort([ pad_l tt pad_r  ]);
      newprots = [ NaN*ones(chans(ichan).nprots,length(pad_l)), chans(ichan).prots, NaN*ones(chans(ichan).nprots,length(pad_r)) ];
      chans(ichan).tt = newtt;
      chans(ichan).prots = newprots;
   end
end
