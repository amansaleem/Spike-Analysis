function exptinfos = ExptReadLogFile(datadir,animal)
% ExptReadLogFile reads the log file for an experiment
%
% exptinfos = ExptReadLogFile(datadir,animal)
%
% Example: exptinfos = ExptReadLogFile('Z:\cat','CATZ005')
%
% 2001-09 Matteo Carandini
% 2002-03 MC fixed little bug with the position of the spaces
% 2002-06 VB crashed when comment line was empty. corrected.
% 
% datadir = 'Z:\cat';
% animal = 'CATZ005';

logfilename = fullfile(datadir,animal,[animal '.txt']);

if exist(logfilename)~=2
   error('cannot find log file');
end

logfilehandle = fopen(logfilename,'r');

keepongoing = 1;

exptinfos = struct('iseries',[],'iexp',[],'StartTime',[],'StartComment',[],'datafile',[],'EndTime',[],'EndComment',[]);

mystring = fgets( logfilehandle );
spaces = findstr(mystring, ' '); s1 = spaces(1);
% in some log files, the date is 01-Dec-01, in others it is 12/01/01...

while keepongoing
   mystring = fgets( logfilehandle );
   spaces = findstr(mystring, ' ');    
   if (length(mystring)==1&&mystring == -1) || isempty(spaces)  % MC added this last ||
      keepongoing = 0;
   else
     
      s1 = spaces(1);
      day = mystring(1:(s1-1));
      time = mystring(s1+[1:5]);
      foo = mystring((s1+11):end);
      spaces = findstr(foo,' ');
      if length(spaces)>1 % comment line not empty.
          firstword = foo(1:(spaces(1)-1));
          switch firstword
          case 'Starting'
             
             [iseriesexp, count, errmsg, nextindex ] = sscanf(foo,'Starting Series %d Exp %d ');
             periods = findstr(foo, '.');
             comment = foo((periods(1)+1):end);
             
             thisexpt = struct('iseries',[],'iexp',[],'StartTime',[],'StartComment',[],'datafile',[],'EndTime',[],'EndComment',[]);
             thisexpt.iseries	= iseriesexp(1);
             thisexpt.iexp		= iseriesexp(2);
             thisexpt.StartTime = [day ' ' time];
             thisexpt.StartComment	= mydeblank(sscanf(comment,'%c'));
             thisexpt.datafile = sprintf('%s_%d_%d', animal, thisexpt.iseries, thisexpt.iexp);
             
             exptinfos(end+1) = thisexpt;
          case 'Interrupted'
             
             periods = findstr(foo, '.');
             comment = foo((periods(1)+1):end);
             
             exptinfos(end).EndTime = [day ' ' time];
             exptinfos(end).EndComment	= mydeblank(sscanf(comment,'%c'));
             
             
          case 'Completed'
             
             periods = findstr(foo, '.');
             comment = foo((periods(1)+1):end);
             
             exptinfos(end).EndTime = [day ' ' time];
             exptinfos(end).EndComment	= mydeblank(sscanf(comment,'%c'));
          end % end switch firstword
      end % end length(spaces)>1
   end % end if mystring == -1
end

fclose(logfilehandle);

exptinfos = exptinfos(2:end);

return

%-----------------------------------------------------------------------

function str = mydeblank(str)
% A HACK TO REMOVE CRAPPY END CHARACTERS
while ~isempty(str) & any(char([32 13 10]) == str(end))
   str(end) = '';
end
