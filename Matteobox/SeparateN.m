function [v,w,s] = SeparateN(rrr,ShowGraphics,ShowGraphicsFull)
% SeparateN takes a N-d matrix and separates it into N vectors
%
% v = SeparateN(rrr) takes a N-dim matrix with dimensions n1 x n2 x ... x nN and
% returns N vectors v{i}, i = 1:N, each with dimension ni x 1
%
% [v,w,s] = SeparateN(rrr) also returns w, the n-1 dim projections, and,
% s, the scaling factors.
%
% Note that each vector will have norm 1, so the amplitudes are
% meaningless.
%
% If N is 2, you are probably better off calling MakeSeparable so you also
% get a scaling factor, residuals, etc.
%
% If some of these dimensions are separable, then the relevant vectors are
% meaningful.
%
% For example, if the data are nstim x nt x nsites, and the different
% sites have different tuning, it does not make sense to look at v{1} and
% v{3}. But if all the sites and all the stimuli give responses with the
% same time course, then v{2} is meaningful (it is THE time course).
%
% NEW FROM FEB 2009: code removes the grand mean before running.
% NEW FROM MAY 2010: no it doesn't!
%
% v = SeparateN(rrr,ShowGraphics) lets you specify whether you want to see
% the individual graphics that are output by MakeSeparable [DEFAULT,
% ShowGraphics = 1] or not (ShowGraphics = 0). Note that this could be
% impractical for huge data sets.
%
% v = SeparateN(rrr,ShowGraphics,ShowGraphicsFull) lets you specify whether
% you want to see all the graphics that show predictions of the separable
% model with one dimension pulled out at a time. For data that has many
% dimensions, this will produce a lot of figures. [DEFAULT, ShowGraphicsFull = 0]
%
% See also: MakeSeparable

% 2007-08 Matteo Carandini
% 2009-02 MC added option to show graphics
% 2009-02 MC imposed removal of grand mean
% 2009-02 ND added option to show graphics of the full optimal separable model
% 2009-04 MC added output of ws (MUST PUT IN HELP TOO)
% 2009-04 ND added output of s (MUST PUT IN HELP TOO)

%% 
if nargin<3
    ShowGraphicsFull = 0;
end
if nargin<2
    ShowGraphics = 0;
end

%% Set up

nd = ndims(rrr);
nn = size(rrr);
n = numel(rrr);

if nd<2
    error('This function requires at least 2-dim input');
end

% rrr = rrr - mean(rrr(:));

%% Do the work

v = {};
w = {};
s = {};
for id = 1:nd
    sss = shiftdim(rrr,id-1); % now the first dimension is the id
    rr = reshape(sss,nn(id),n/nn(id));
    [w{id},v{id},s{id}] = MakeSeparable(rr,ShowGraphics);
    newdims = nn; newdims(id) = [];
    w{id} = reshape(w{id},newdims);
end
clear foo

%% Graphics

if ShowGraphicsFull
    figure; ax = [];
    for id = 1:nd
        ax(id) = subplot(nd,1,id);
        plot(v{id});
    end
end

%% Done

return







%% code to test the function (2 dims)

n1 = 10;
n2 = 20;

r = {};
r{1} = gaussian(2,1:n1);
r{2} = gaussian(2,(1:n2)-n2/2);

nd = length(r);

figure; ax = [];
for id = 1:nd
    ax(id) = subplot(nd,1,id);
    plot(r{id});
end

rrr = zeros(n1,n2);
for i1 = 1:n1
    for i2 = 1:n2
        rrr(i1,i2) = r{1}(i1)*r{2}(i2);
    end
end

v = SeparateN(rrr);

%% code to test the function (3 dims)

n1 = 10;
n2 = 20;
n3 = 30;

r = {};
r{1} = gaussian(2,1:n1);
r{2} = gaussian(2,(1:n2)-n2/2);
r{3} = sawtooth(1:n3);

nd = length(r);

figure; ax = [];
for id = 1:nd
    ax(id) = subplot(4,1,id);
    plot(r{id});
end

rrr = zeros(n1,n2,n3);
for i1 = 1:n1
    for i2 = 1:n2
        for i3 = 1:n3
            rrr(i1,i2,i3) = r{1}(i1)*r{2}(i2)*r{3}(i3);
        end
    end
end

rrr = rrr + + normrnd(0,.2,[n1 n2 n3]);

v = SeparateN(rrr);

%% code to test the function (4 dims)

n1 = 10;
n2 = 20;
n3 = 30;
n4 = 40;

r = {};
r{1} = gaussian(2,1:n1);
r{2} = gaussian(2,(1:n2)-n2/2);
r{3} = sawtooth(1:n3);
r{4} = square(1:n4);

nd = length(r);

figure; ax = [];
for id = 1:nd
    ax(id) = subplot(4,1,id);
    plot(r{id});
end

rrr = zeros(n1,n2,n3,n4);
for i1 = 1:n1
    for i2 = 1:n2
        for i3 = 1:n3
            for i4 = 1:n4
                rrr(i1,i2,i3,i4) = r{1}(i1)*r{2}(i2)*r{3}(i3)*r{4}(i4);
            end
        end
    end
end

rrr = rrr + + normrnd(0,.2,[n1 n2 n3 n4]);

v = SeparateN(rrr);

