function [ pref, circvar, amp ] = circstats( yy, oo )
%CIRCSTATS circular statistics (circular variance, preferred angle)
%
% [ pref, circvar, amp ] = circstats( yy, oo ) returns the preferred angle and
% the circular variance of data yy with angles oo. yy can be a tensor with
% dimensions [no, nx, ny], where no is the length of oo. oo is assumed in
% radians, between 0 and 2*pi.
% 
% Careful, if you are in the orientation domain (0 to 180 deg) you have to
% map that to [0 2*pi].
%
% See also CIRCSTATS360, CIRCCORR.

% 2004-11 Matteo Carandini
% 2011-05 MC added amp output
    
% nstim   = size(yy,1);
nx      = size(yy,2);

if ndims(yy) == 3
    ny      = size(yy,3);
else
    ny = 1;
end

if length(unique(oo))<3
    disp('Cannot do circular statistics with less than 3 angles');
    pref 	    = zeros(nx, ny);
    circvar 	= zeros(nx, ny);
    return
end

if ny==1
    pp = exp(1i*oo(:)); % in orientation code: *2
    % mm = sum( yy(:).*pp )./ sum( yy(:) );
    f1 = sum( yy(:).*pp );
    f0 = sum( yy(:) );
    pow = sum( abs(pp).^2 );
    % now orientations in [0 pi] have been stretched onto [0 2*pi]
else
    pp = repmat( exp(1i*oo(:)), [ 1, nx, ny ] ); % in orientation code: *2
    % mm = squeeze(sum( yy.*pp, 1 ))./ squeeze(sum( yy, 1 ));
    f1 = squeeze(sum( yy.*pp, 1 ));
    f0 = squeeze(sum( yy, 1 ));
    pow = squeeze(sum( abs(pp).^2, 1 ));

end


% mm(isnan(mm)) = 0;
f1(isnan(f1)) = 0;
f0(isnan(f0)) = 0;

% pref    = angle(mm); 
% circvar = 1 - abs(mm);  

pref = angle(f1);
circvar = 1 - abs(f1)./abs(f0);
amp  = abs(f1)./abs(pow);

% % in optical imaging code:
% pref = angle(mm)+pi; % now from 0 to 2*pi
% pref = pref/2 - pi/2; % finally, bet -pi/2 and pi/2

 
return 

%-------------------------------------------------------
%           code to test the function
%-------------------------------------------------------

oo = 2* pi* [0 1/4 1/2 3/4 ]';

yy = [ 0 10 20 10 ];
yy = [ 10 10 10 10 ];
yy = [ 0 10 0 0 ];
[ pref, circvar ] = circstats( yy, oo );

% figure; plot(oo,yy);

nx = 3;
ny = 4;
yy = zeros(4,nx,ny);
for ix = 1:nx
    for iy = 1:ny
        yy(:,ix,iy) = [ 0 10 0 0 ];
    end
end
[ pref, circvar ] = circstats( yy, oo );

