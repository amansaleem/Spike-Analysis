function p = circle( ctr, rad, opts, opts2 )
% CIRCLE draws a circle
%	
%	p = circle( ctr, rad, opts )
%
%	circle( ctr, rad ) defaults opts to '-'
%
% 1995 Matteo Carandini
% 2002 VM added the handle as an output
% part of the Matteobox toolbox

if nargin == 2; opts = '-'; end
if nargin <4; opts2= 0.75; end

theta = linspace(-pi,pi);

p = plot( ctr(1)+rad*cos(theta), ctr(2)+rad*sin(theta), opts , 'linewidth', opts2);
 
