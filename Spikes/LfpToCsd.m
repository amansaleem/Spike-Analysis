function myCsd = LfpToCsd( lfp, icSpacing, n, upscale, sigma, options )
% LfpToCsd computes the current source density from the local field
% potential
% 
% myCsd = LfpToCsd( lfp, icSpacing)
% gets the current source density for lfp tensor of size nx X nz X nt.
% icSpacing is the inter-contact spacing in mm
%
% myCsd = LfpToCsd( lfp, icSpacing, n)
% n*icSpacing is the differentiation grid. n can be either 1 or 2.
% DEFAULT: 1
%
% myCsd = LfpToCsd( lfp, icSpacing, n, upscale)
% upscale is the amount by which the data are upsampled (data will
% automatically be smoothed) upsampling applied in depth dimension
% DEFAULT: 8
%
% myCsd = LfpToCsd( lfp, icSpacing, n, upscale, sigma)
% sigma is the width of the Gaussian used for filtering. If 0 no additional spatial
% filtering will be applied (except for that done by upsampling) filtering
% applied in depth dimension
% DEFAULT: 0.8
%
% myCsd = LfpToCsd( lfp, icSpacing, n, upscale, sigma, options) options is
% a structure that has four fields, plotFlag, smooth_win_z, smooth_win_t
% and interp_factor. 
% If plotFlag is true then displays the computed CSD as an image with the
% LFP traces superimposed (DEFAULT : 0). 
% smooth_win_z is the standard deviation of a guassian that is used to
% smooth the CSD in the depth dimension. If 0 no additional spatial
% filtering will be applied except for that done by upsampling (DEFAULT :
% 0.8). Units are in z-samples.
% smooth_win_t is the standard deviation of a gaussian that is used to
% smooth the CSD in the time dimension. If 0 no temporal filtering will be
% applied (DEFAULT : 0). Units are in t-samples.
% interp_factor is the multiple by which the CSD is upsampled (DEFAULT :
% 8).

% For reasons of backward compatibility, parameters upscale and sigma are
% reassigned to options.interp_factor and options.smooth_win_z,
% respectively. if options is just a flag, then assumed to be plotFlag.
% 
% 
% 2009      Neel Dhruv created it
% 2009-09   LB edited it
% 2010-12   ND edited the defaults and required set sigma==0 if want no extra filtering
% 2010-12   ND edited the smoothing to use filtfilt rather than filter
% 2010-12   ND added option of doing smoothing in both z-dim and time. have
% to select it (line 99) to make it work
% 2011-01   ND made options the input parameter to supercede plotFlag and
% to incorporate options for the smoothing and the interpolation. if both
% smooth_win_t and smooth_win_t are defined automatically does 2-D
% smoothing. should be backwards compatible.


if exist('options','var')
    if ~isstruct(options), options = struct('plotFlag',options); end
else
    options = struct();
end

if nargin < 5 || isempty(sigma), sigma = 0.8; end
if ~isfield(options,'smooth_win_z'), options.smooth_win_z = sigma; end

if nargin < 4 || isempty(upscale), upscale = 8; end
if ~isfield(options,'interp_factor'), options.interp_factor = upscale; end

if ~isfield(options,'plotFlag'), options.plotFlag = false; end
if ~isfield(options,'smooth_win_z'), options.smooth_win_z = 0.8; end
if ~isfield(options,'smooth_win_t'), options.smooth_win_t = 0; end
if ~isfield(options,'interp_factor'), options.interp_factor = 8; end

if nargin < 3 || isempty(n)
    n = 1;
end

if ~ismember(n, [1 2])
    error('<LfpToCsd> input n can only be 1 or 2 (is %d)', n);
end

if n ~= 1
    error('<LfpToCsd> n other than 1 not implemented yet');
end

[nr, nz, nt] = size(lfp);

%% quick hack added by Matteo on 21 Dec 2011 to deal with nr>1

if nr>1
    % warning('Because nr > 1, we are going into a quick hack');
    myCsd = [];
    for ir = 1:nr
        myCsd(ir,:,:) = LfpToCsd( lfp(ir,:,:), icSpacing );
    end
    return
end

%% go on as usual


% set up differential filter for computing 1-D partial second derivative in depth 
% dd = [1,-2,1];
% myCsd = filter(dd, 1, lfp, [], 2) / (icSpacing^2); % differentiating filter and scale
% myCsd = myCsd(1,3:end,:); % corresponds directly to the iterative solution below
myContact = 0;
myCsd = nan(1, nz-2, nt);
for isite = 2:size(lfp,2)-1
    myContact = myContact + 1;
    myCsd(1,myContact,:) = (lfp(1,isite-1,:) - 2*lfp(1,isite,:) + lfp(1,isite+1,:))/icSpacing^2;
end
myCsd = myCsd/1000; % convert from uV/mm^2 to mV/mm^2

% upsample
yy = shiftdim(myCsd,1);
yy = resample(double(yy), options.interp_factor, 1); % careful, edge effects can occur here if first and last line of yy are not close to zero
myCsd = permute(yy,[3 1 2]);

% smooth
vector_z = ((-min(round(3*options.smooth_win_z),2)*options.interp_factor): ...
    min(round(3*options.smooth_win_z),2)*options.interp_factor); % this looks complicated because need to keep the length of this down in order to do filter
gg_z = gaussian( options.smooth_win_z*options.interp_factor, vector_z );
vector_t = -round(3*options.smooth_win_t):round(3*options.smooth_win_t);
gg_t = gaussian( options.smooth_win_t, vector_t );
gg_zt = gg_z'*gg_t; % 2-D gaussian filter

%%% please note, using filter results in a shift by 0.5*length(gg) so
%%% filtering in depth is going to result in depths being too low or
%%% filtering in time is going to result in a spurious extra delay
%     myCsd = filter(gg, 1, myCsd, [], 3); % smoothing in time
%     myCsd = filter(gg, 1, myCsd, [], 2); % smoothing in depth

% filter using filtfilt
if ~isempty(options.smooth_win_z) && options.smooth_win_z~=0
    filteredCsdz = filtfilt(gg_z, 1, squeeze(myCsd)); % smoothing in depth
    myCsd = permute(filteredCsdz, [3 1 2]);
elseif ~isempty(options.smooth_win_t) && options.smooth_win_t~=0
    filteredCsdt = filtfilt(gg_t, 1, squeeze(myCsd)'); % smoothing in time
    myCsd = permute(filteredCsdt, [3 2 1]);
end
% 2-D filter using filter2
if ~isempty(options.smooth_win_z) && ~isempty(options.smooth_win_t) && ...
        options.smooth_win_z~=0 && options.smooth_win_t~=0
    filteredCsd2d = filter2(gg_zt, squeeze(myCsd)); % smoothing in depth and time
    myCsd = permute(filteredCsd2d, [3 1 2]);
end

% timing information could be useful for this plot
if options.plotFlag
    
    yAx = 1 + options.interp_factor : nz*options.interp_factor - options.interp_factor;
    
    figure('Name', 'CSD with LFPs');
    imagesc(1:size(myCsd,3), yAx, squeeze(myCsd)); hold on; % plot the csd
    for iz = 1 : nz
        plot(1:size(lfp,3), -1*0.1*squeeze(lfp(1,iz,:)) + (options.interp_factor - (options.interp_factor -1)*0.5) + (iz-1)*options.interp_factor, 'k'); % we need the -1 here because the y-axis goes from small numbers on top to large numbers on the bottom
    end
    set(gca, 'YLim', [0 nz*options.interp_factor + options.interp_factor], 'YTick', (options.interp_factor - (options.interp_factor -1)*0.5):options.interp_factor:nz*options.interp_factor, 'YTickLabel', 1:nz, ...
        'Box', 'off', 'TickDir', 'out');
    xlabel('Time (samples)');
    ch = colorbar;
    set(get(ch,'XLabel'),'String','mV/mm^2')
end
