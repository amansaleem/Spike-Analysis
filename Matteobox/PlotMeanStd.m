function PlotMeanStd( xx, yyy, c )
% Plots the mean and standard deviation
%
% PlotMeanStd( xx, yyy ). Number of rows of yyy should be same as number of
% points of xx.
%
% PlotMeanStd( xx, yyy, c ) lets you specify the color c. DEFAULT: 'k'.
%
% 2012-03 MC

if nargin<3
    c = 'k';
end

npoints = length(xx(:));
if size(yyy,1)~=npoints
    error('Number of rows of yyy should be same as number of points of xx');
end

mm = mean( yyy, 2 ); 
% ee = std( yyy, [], 2 );
ee = sem( yyy');
errorbar( xx, mm, ee, [c '.']); hold on
plot( xx, mm, [c 'o'], 'markerfacecolor',c); 

