function f = supertitle(s,pos,interprtr)
% SUPERTITLE 	makes a big title over all subplots
%
%	supertitle(S) writes the string S as a title
%	It returns the handle if you want it.
%
%	supertitle(S,pos) lets you specify a position be 0 and 1. Default is .95
%
%   ------ added by ND on 2012_06_20-------
%   supertitle(S,pos,interprtr) lets you specify an interpreter. Default is
%   'tex' if you want literal characters use 'none'

% 1995,1996 Matteo Carandini
% part of the Matteobox toolbox

if nargin < 2, pos = .95; end
if nargin < 3, interprtr = 'tex'; end

currax = gca;
axes('position',[.5 pos 0.1 0.1 ],'xcolor',[0 0 0],'visible','off');

if ~isempty(s)
	f = text(0 , 0, s, 'units', 'normalized',...
		'verticalalignment','top','horizontalalignment','center','interpreter',interprtr);
end
axes(currax);
