function ww = robustassgn_mtx( vv, ii )
% ROBUSTASSGN robust subscript assignment (NaNs for subscripts out of range)
%
% ww = robustassgn( vv, ii )
% is like ww = vv(ii) but returns NaN where the ii are out of range.
%
% 2006-03 Matteo Carandini
%
% 2012_06_20 ND made it work for matrices in that it uses all rows and then
% robustly assigns the columns

[nrows,ncols] = size(vv);

goodcols = ( ii >= 1 & ii <= ncols );
ww = NaN(nrows,length(ii));
ww(:,goodcols ) = vv(:,ii(goodcols));
