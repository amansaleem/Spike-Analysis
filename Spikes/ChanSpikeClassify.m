function chan = ChanSpikeClassify(vv,chan)
% ChanSpikeClassify classifies spikes
%
% The syntax is chan = ChanSpikeClassify(vv,chan), where vv are candidate spikes, and
% chan.prots is an array of prototype spikes (nprots x length)
% chan.mindur specifies how many ms after time 0 to include (default if undefined: 2)
% chan.tt
%
% The results are: 
% chan.cc, a vector that assigns each candidate to a class (0 is unassigned) 
% 
% chan.mincorr specifies the minimal required correlation with a prototype for 
% the final assignment of a spike to a class (default if undefined: 0.85). 
%
% 1999-10 Matteo Carandini 
% 2000-01 MC
% 2000-09 MC
% 2004-02 MC based discrimination on distance rather than correlation

if nargin<2
   error('Call with one argument is obsolete. Must now call with two arguments');
end

fprintf('Classifying candidate spikes...');

prots = chan.prots;

ncands = size(vv,1);

nprots = size(prots,1);
c = zeros(ncands,1); % all unclassified in the beginning
if nprots == 0
   chan.cc = c;
   return
end

% drop data points after mindur, so they don't affect the correlation
if isfield(chan,'mindur') && isfield(chan,'tt') && ~isempty(chan.mindur) && ~isempty(chan.tt)
   samplelist = find(chan.tt<=chan.mindur & chan.tt>=-0.5 & ~isnan(prots(1,:)));
   vv = vv(:,samplelist);
   prots = prots(:,samplelist);
   % tt = chan.tt(samplelist); % not sure this is needed?
end
nsamples = size(vv,2);

mincorr = 0.85;
if isfield(chan,'mincorr') 
   if ~isempty(chan.mincorr)
      mincorr	= chan.mincorr;
   end
end

mvv = vv - mean(vv,2)*ones(1,nsamples);
mprots = prots - mean(prots,2)*ones(1,nsamples);

% ncands by nprots, correlation bet cand and prot

cc = ( mvv * mprots')/(nsamples-1)./sqrt(var(mvv,[],2)* var(mprots,[],2)'); 

dd = zeros(ncands,nprots); % matrix of distances
for iprot = 1:nprots
    dd(:,iprot) = std(mvv - ones(ncands,1)*mprots(iprot,:),0,2);
end

[foo, c] = min(dd,[],2);
for icand = 1:ncands
   if all(cc(icand,:)<mincorr)
      c(icand) = 0;
  end
end

%% Graphics

persistent ScatterFig

if isempty(ScatterFig), ScatterFig = figure; end
figure(ScatterFig);
if nprots > 1
    set(ScatterFig,'visible','on');
    % hack to avoid NaNs, added by MC 2009-11-02:
    vv(~isfinite(vv)) = 0;
    [coeff, score] = princomp( vv, 'econ' );
    [h,ax] = gplotmatrix( [[0,0,0]; score(:,1:3)], [], [0; c], 'wrgbc', '.', 6, 'off', 'hist', {'PC1','PC2','PC3'} );
    % [h,ax,bigax] = gplotmatrix( [[0,0]; score(:,1:2)], [], [0; c], 'wrgbc', '.', 6, 'off', 'hist', {'PC1','PC2'} );
    set(ax,'Color',0.3*[1 1 1]);
else
    set(ScatterFig,'visible','off');
end

%% write the output

chan.cc = c;

fprintf('...done\n');


