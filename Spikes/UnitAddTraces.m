function u = UnitAddTraces( u, gain, offset)
% Plots the tuning for traces in an expt
%USAGE: u is a unit structure
%       gain is gain of multiclamp software
%       offset is voltage outside of cell
% 2010-06 Matteo and Bilal -- from ExptShowTraces


SetDefaultDirs;

if nargin<3
    offset = 0;
end

if nargin<2
    gain = 1;
end

if nargin < 1
    error('We do need an input');
end

expt = ExptLoad(u);

if isempty(expt)
    fprintf('No data?\n');
    return
end

p = ProtocolLoad(expt);

u.traces = cell(expt.nstim,expt.nrepeats);

for istim = 1:expt.nstim
    for irpt = 1:expt.nrepeats
        if ~isempty(expt.data{1}{istim, irpt})
            vv = double(expt.data{1}{istim, irpt});
            vv = (vv*3276.8)/1000; % takes care of unitspervolt in ExptLoad line 214; then converted to raw mV per V
            vv = ((vv*gain) + offset);
            u.traces{istim,irpt} = vv;
        end
    end
end

