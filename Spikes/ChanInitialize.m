function chans = ChanInitialize(expt,chanlist)
% ChanInitialize initializes a channel structure
%
%	chans = ChanInitialize(expt)
%
%	chans = ChanInitialize(expt,chanlist)
%
%  chans = ChanInitialize;
%
% 2000 Matteo Carandini
% 2007-09 MC corrected a couple of bugs

% MC 2007-09-10 changed defaults:
% nprots from 2 to 1
% mindur from 1 to 2

if nargin == 0
   chans = struct(...
      'done',0,...
      'spikes',[],...
      'thresh',[],...
      'threshsign',[],...
      'mindur',1.2,...
      'tt',[],...
      'nprots',1,...
      'frozenprots',[1 0 0 0],...
      'prots',[],...
      'mincorr',0.1,...
      'cellids',[999 NaN NaN NaN],...
      'cc',[],...
      'iexptchan',[],...
      'channame',[],...
      'filtercoeffs',[]);
   return;
end

if nargin<2
   chanlist = 1:expt.nchans;
else
   if any(chanlist>expt.nchans | chanlist<1)
      error('argument chanlist has problems');
   end
end

% preallocation:
chans = repmat( ChanInitialize, length(chanlist), 1 );

for iichan = 1:length(chanlist)
   iexptchan = chanlist(iichan);
   
   chans(iichan) = ChanInitialize;
   
   chans(iichan).iexptchan = iexptchan;
   chans(iichan).channame = expt.channames(iexptchan);
   
   if ~any(isnan(expt.threshsigns)) && length(expt.threshsigns)==expt.nchans
      % use the threshold sign that was set during the exp
      chans(iichan).threshsign = expt.threshsigns(iexptchan);
   else
      chans(iichan).threshsign = -1; 
   end
   
   if ~any(isnan(expt.threshvalues)) && length(expt.threshvalues)==expt.nchans
      % use the threshold that was set during the exp
      chans(iichan).thresh = expt.threshvalues(iexptchan);
   else
      % place the thresh in a smart place, at 2.5 stds
      v = [expt.data{iexptchan}{:,:}];
      meanv = mean(v);
      minv = min(v);
      maxv = max(v);
      if maxv-meanv > meanv-minv, sgn = 1; else sgn = -1; end
      chans(iichan).thresh = meanv + sgn*2.5*std(single(v)); 
   end 
   
   dt = 1000/expt.samplerate;
   % this was the standard: tlim = [-0.75 1.50];	
   tlim = [-1.50 4.00];	
   chans(iichan).tt = sort([ 0:-dt:tlim(1) , dt:dt:tlim(2) ]);
   
   chans(iichan).filtercoeffs = expt.filtercoeffs;
   
end



