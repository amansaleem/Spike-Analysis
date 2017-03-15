function f = oritune180(pars,oris) 
% ORITUNE180 	gaussians living on a circle, for orientation tuning
% 
%		oritune([Op, Rp, Ro, sigma], oris), where
%		Op is the preferred direction (bet 0 and 180)
%		Rp is the response to the preferred orientation;
%		Ro is the background response (useful only in some cases)
%		sigma is the tuning width;
% 		oris are the orientations, bet 0 and 180.
%
%		oritune([Op, Rp, Ro, sigma], 'spiel'), gives you a little description of the tuning.
%
% See also: ORITUNE, FITORI, CIRCSTATS, CIRCSTATS360

% 2013-05 MC from Oritune

% 1997-1998 Matteo Carandini
% part of the Matteobox toolbox

if length(pars(:))~= 4
   error('You must specify all 4 parameters');
end

Op = pars(1); Rp = pars(2); Ro = pars(3); sigma = pars(4);

if ischar(oris) && strcmp(oris,'spiel')
   disp(['Width at half height is ' num2str(sqrt(log(2)*2)*sigma,3) ' deg']);
else
   anglesp = 180/(2*pi)*angle(exp(1i*(oris-Op)    *(2*pi)/180));   
   f = Ro +  Rp*	exp( -anglesp.^2 ./ (2 * sigma^2));
end


