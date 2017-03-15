
function PFileInitialize(xfilename)
% Creates a rudimentary pfile given an xfile
%
% PFileInitialize('foo.x') creates a pfile called fooDemo.p in DIRS.xfiles,
% containing only one stimulus (the default one)
%
% 2010-10 Matteo Carandini

global DIRS

% xfilename = 'vmovie3sequentialGrating.x'

x = XFileLoad( xfilename );

pfilename = [ xfilename(1:end-2) 'Demo.p' ];

pfile = fopen(fullfile(DIRS.xfiles,pfilename),'w');

fprintf(pfile,'%s\r\n',xfilename);

fprintf(pfile,'1 %d 0\r\n',x.npars);
fprintf(pfile,'0 1 1\r\n');


for ipar = 1:x.npars
    fprintf(pfile,'%d ',x.pardefaults(ipar));
end
fprintf(pfile,'\r\n');

fclose(pfile);

fprintf(1,'Wrote %s\n',pfilename);
