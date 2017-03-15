function yy = expfunc(pars,tt)
% EXPFUNC exponential decay function
%
%	yy = expfunc([y0,dy,tau],tt) gives yy = (y0+dy)- dy* exp(-tt/tau);
%
% In other words, y0 is the value at time 0, y0+dy is the value at time
% infinity.
%
%	yy = expfunc([y0,dy,tau, t0],tt)
% gives yy = (y0+dy)- dy* exp(-(tt-t0)/tau);
% 
% 1996 Matteo Carandini
% 2013-11 MC added 4th parameter

% part of the Matteobox toolbox

y0  = pars(1);
dy  = pars(2);
tau = pars(3);
t0 = 0;

if length(pars)==4, t0 = pars(4); end

yy = (y0+dy)- dy* exp(-(tt-t0)/tau);
