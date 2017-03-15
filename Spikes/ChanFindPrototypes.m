function chan = ChanFindPrototypes(vv, chan, irepeat)
% ChanFindPrototypes finds prototype spikes using unsupervised learning.
%
% This function uses unsupervised learning to find chan.prots, a matrix of
% prototype spikes. 
%
% chan = ChanFindPrototypes(vv, chan,irepeat) works only on spikes obtained
% during repeat irepeat 
%
% chan.nprots determines how many prototypes should be found
%
% 1999-10 Matteo Carandini 
% 2000-01 MC
% 2000-03 MC
% 2000-09 MC
% 2010-03-09 MC candidates considered only in relevant time intervals

if nargin<2
   error('Must give at least 2 arguments.');
end

if ~isfield(chan,'prots')
   chan.prots = [];
end

if isempty(chan.nprots), chan.nprots = 0; end

if chan.nprots == 0, return; end

nprots = size(chan.prots,1);
nsamples = length(chan.tt);

% added by MC 2010-03-09 so distance is only computed in tight time
% interval
ii = (chan.tt > -0.5 & chan.tt < chan.mindur);

if chan.nprots < nprots, error('nprots is less than existing number of prots'); end

% initialize the missing prots
if chan.nprots>nprots
   chan.prots((nprots+1):chan.nprots,:) = zeros(chan.nprots-nprots,nsamples);
end

if isempty(chan.cc) || length(chan.cc)~=length(chan.spikes)
    unassigned = ones(size(chan.spikes));
else
    unassigned = (chan.cc==0);
end

if nargin == 3
    p = find([chan.spikes(:).irpt] == irepeat);
    vv = vv(p,:); % a matrix whose rows are candidate spikes
    unassigned = unassigned(p);
end
ncands = size(vv,1);

if ncands<1
   % return to the calling function. New prots (if any) are all zero
   return
end

if chan.nprots>nprots
   % VERY OLD: add random prototypes
   % chan.prots(nprots+1:chan.nprots,:) = rand(chan.nprots-nprots,nsamples)-0.5;
   % OLD: add prototypes made of the mean v
   % chan.prots((nprots+1):chan.nprots,:) = ones(chan.nprots-nprots,1)*mean(vv,1);
   % NEW: add prototypes made of the mean v for the unassigned candidates
   % (if they exist)
   if any(unassigned)
       chan.prots((nprots+1):chan.nprots,:) = ones(chan.nprots-nprots,1)*mean(vv(unassigned,:),1); % was find(unassigned)
   else
       chan.prots((nprots+1):chan.nprots,:) = ones(chan.nprots-nprots,1)*mean(vv,1);
   end
end

if size(chan.prots,1) ~= chan.nprots || size(chan.prots,2) ~= nsamples
   error('Something wrong with the size of the spike prototype matrix. Maybe one of your repeats has no candidate spikes?');
end


ws = [ 0.5 0.5 0.5 0.4 0.2 0.1 0.05 0.02 0.01 0.005 0.002 0.001];
ws = ws(1:2:end);

if all(chan.frozenprots(1:chan.nprots))
   return
end

notfrozen = ~chan.frozenprots;

niters = length(ws);
dd = zeros(chan.nprots,1); % the distance from each prototype to a given candidate
% classif = NaN*zeros(ncands,niters); % this is useful to track the performance
fprintf('Finding prototypes');
for iiter = 1:niters
    fprintf('.'); % disp(['Iteration ' num2str(iiter)]);
    % schedule of weights is important
    for icand = 1:ncands
        if ~all(isnan(vv(icand,:)))
            % find the closest prototype
            for iprot = 1:chan.nprots
                % replaced var with more basic sum, after suggestion by Mark Staeheli
                dd(iprot) = nansum((vv(icand,ii)-chan.prots(iprot,ii)).^2)/nsamples;
                % MC added ii 2010-03-09
            end
            mindd = nanmin(dd);
            closestprot= find(dd==mindd(1)); closestprot= closestprot(1);
            
            % Bring the prototype closer to the candidate
            if notfrozen(closestprot)
                chan.prots(closestprot,:) = (1-ws(iiter))*chan.prots(closestprot,:)+ws(iiter)*vv(icand,:);
            end
            
            % TO AVOID DEAD PROTOTYPES move all the prototypes in that direction, with a smaller step
            for iprot = 1:chan.nprots
                if notfrozen(iprot)
                    nans = find(isnan(chan.prots(iprot,:)));
                    chan.prots(iprot,nans) = vv(icand,nans);
                    chan.prots(iprot,:) = (1-ws(iiter)/1000)*chan.prots(iprot,:)+ws(iiter)/1000*vv(icand,:);
                end
            end
        end
    end
end
fprintf('done\n');

