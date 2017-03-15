function filename = UnitGetFilename( animal, iseries, iexp, ichan, icell )
% UnitGetFilename the file name of a file for a Unit data structure
%
% filename = UnitGetFilename( animal, iseries, iexp, ichan, icell )
% gives you something like CAT123/1/CAT123_s01_e03_c2_u001.mat
%
% filename = UnitGetFilename( animal, iseries, iexp, '*', '*' )
%
% 2000 Matteo Carandini

if iseries < 100
	strSeries = num2str(100+iseries); strSeries = strSeries(2:3);
else
   strSeries = num2str(iseries);
end

if iexp < 100
   strExp = num2str(100+iexp); strExp = strExp(2:3);
else
   strExp = num2str(iexp);
end

if ischar(icell)
   strCell = icell;
else
   if icell < 1000
      strCell = num2str(1000+icell); strCell = strCell(2:4);
   else
      strCell = num2str(icell); 
   end
end

if ischar(ichan)
   strChan = ichan;
else
   strChan = num2str(ichan); 
end

filename = [ animal filesep num2str(iseries) filesep ...
   animal '_s' strSeries '_e' strExp '_c' strChan '_u' strCell '.mat'];
