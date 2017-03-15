function filename = UnitGetFilename( animal, iseries, iexp, ichan, icell )
% UnitGetFilename the file name of a file for a Unit data structure
%
% filename = UnitGetFilename( animal, iseries, iexp, ichan, icell )
% gives you something like CAT123/1/CAT123_s01_e03_c2_u001.mat
%
% filename = UnitGetFilename( animal, iseries )
% gives you something like CAT123/1/
%
% filename = UnitGetFilename( animal, iseries, iexp, '*', '*' )
% is useful for searching
%
% 2000 Matteo Carandini
% 2009-09-29 AZ slightly optimized performance
% 2010-03-08 MC added support for two inputs

% if iseries < 100
% 	strSeries = num2str(100+iseries); strSeries = strSeries(2:3);
% else
%    strSeries = num2str(iseries);
% end

% if iexp < 100
%    strExp = num2str(100+iexp); strExp = strExp(2:3);
% else
%    strExp = num2str(iexp);
% end

if nargin == 2
    filename = sprintf('%s%i%s',[animal filesep],iseries,filesep);
    return
end

    
if ischar(icell)
   strCell = icell;
else
%    if icell < 1000
%       strCell = num2str(1000+icell); strCell = strCell(2:4);
%    else
%       strCell = num2str(icell);
      strCell = sprintf('%03i',icell); %AZ 2009-09-29
%    end
end

if ischar(ichan)
   strChan = ichan;
else
%    strChan = num2str(ichan);
   strChan = sprintf('%i',ichan); % AZ 2009-09-29
end

% filename = [ animal filesep num2str(iseries) filesep ...
%    animal '_s' strSeries '_e' strExp '_c' strChan '_u' strCell '.mat'];

% AZ 2009-09-29
filename = sprintf('%s%i%s_s%02i_e%02i_c%s_u%s.mat',...
   [animal filesep],iseries,[filesep animal],iseries,iexp,strChan,strCell);