function [status,framesequence] = ProtocolTimingIsBuggy(protocol,myscreen)
% ProtocolTimingIsBuggy checks whether a protocol might have timing problems
%
% [status,framesequence] = ProtocolTimingIsBuggy(protocol,myscreen)
%
% status is: 0 if the xfiles that are not buggy
%            1 if the xfile is buggy, but the timing can be corrected
%            2 if the xfile is buggy, and timing is hard to correct
% framesequence is a structure containing the frame sequence as it should
% have been shown (can be used to correct the bug when UnitLoad doesn't)
%
% 2003-04 VM made it
% 2005-12 VM replaced FrameRate with StimFrameRate

if nargin < 2
   myscreen = [];
end

status = 0;
framesequence = [];

switch protocol.xfile(1:end-2)
   
% Cases where one can easily fix the timing   
case {'visdriftsin','visdriftsin100'}
   % The first frame is buggy
   status = 1;
   if ~isempty(myscreen) & ~isempty(myscreen.StimFrameRate)
      stimdurs = protocol.pars(1,:)/10;
      nframes = round(myscreen.StimFrameRate .* stimdurs);
      for istim = 1:protocol.nstim
         framesequence{istim} = ones(1,nframes(istim));
      end
   end
   
% Cases where the timing is not that easy to fix   
case 'visdriftsin2'
   % The first two frames are buggy
   status = 2;
   if ~isempty(myscreen) & ~isempty(myscreen.StimFrameRate)
      stimdurs = protocol.pars(1,:)/10;
      nframes = round(myscreen.StimFrameRate .* stimdurs);
      for istim = 1:protocol.nstim
         framesequence{istim} = ones(1,nframes(istim)*2);
         framesequence{istim}((1:nframes(istim))*2) = 2;
      end
   end
   
case 'vissweep'
   % There is only one frame, might be fixable easily
   status = 2;
   if ~isempty(myscreen) & ~isempty(myscreen.StimFrameRate)
      tadapt = protocol.pars(1,:)/10;
      tsweep = protocol.pars(2,:)/10;
      twait	= protocol.pars(3,:)/10;
      nupsweep = ceil(tsweep/2*myscreen.StimFrameRate);
      nsweep = max(0,2*nupsweep - 1);
      nadapt = ceil(tadapt*myscreen.StimFrameRate);
      nwait = ceil(twait*myscreen.StimFrameRate);
      for istim = 1:protocol.nstim
         framesequence{istim} = ones(1,nadapt(istim)+nsweep(istim)+nwait(istim));
      end
   end

case 'vissweep2grat'
   % There are two frames, harder to fix
   status = 2;
   if ~isempty(myscreen) & ~isempty(myscreen.StimFrameRate)
      tadapt = protocol.pars(1,:)/10;
      tsweep = protocol.pars(2,:)/10;
      twait	= protocol.pars(3,:)/10;
      nupsweep = ceil(tsweep/2*myscreen.StimFrameRate);
      nsweep = max(0,2*nupsweep - 1);
      nadapt = ceil(tadapt*myscreen.StimFrameRate);
      nwait = ceil(twait*myscreen.StimFrameRate);
      for istim = 1:protocol.nstim
         framesequence{istim} = ones(1,nadapt(istim)+nsweep(istim)+nwait(istim))*2;
         framesequence{istim}(nadapt(istim)+1:nadapt(istim)+nsweep(istim)) = 1;
      end
   end

% Buggy, but unlikey that someone will ever try to fix them
case {'vismradapt';'visringlog'}
   status = 2;
   
end

