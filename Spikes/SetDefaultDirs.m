% SetDefaultDirs sets the default directories
% 
% 2013-09-24 Matteo Carandini
% 2013-10-10 AS: line 44 & 47, corrected the variable name

% Did I screw up anything? If so, please see the older version at
% \\zserver.ioo.ucl.ac.uk\Code\Archive\Spikes\Spikes 2013-09-24
% and let me know... Thanks 
% Matteo

% At some point we may want to use IP addresses to determine if we are at
% Rockefeller or at Bath Street. If the former, we should probably use
% zcloneX instead of zserverX. A method to determine IP address is:
% address = java.net.InetAddress.getLocalHost 
% IPaddress = char(address.getHostAddress)

global DIRS serverName server2Name server3Name

if isunix 
    
    serverName    = '/mnt/zserver';
    server2Name   = '/mnt/zserver2';
    server3Name   = '/mnt/zserver3';
    
    DIRS.Temp      = '/tmp';
else   
    
    serverName     = '\\zserver.ioo.ucl.ac.uk';
    server2Name    = '\\zserver2.ioo.ucl.ac.uk';  
    server3Name    = '\\zserver3.ioo.ucl.ac.uk'; 
    
    if isdir('D:\Temp')
        DIRS.Temp       = 'D:\Temp'; 
    else
        DIRS.Temp       = 'C:\Windows\Temp'; 
    end
    
end

if ~isdir(fullfile(serverName,'Data'))
    fprintf('Make sure directory %s is accessible!\n',fullfile(serverName,'Data'));
end
if ~isdir(fullfile(server2Name,'Data'))
    fprintf('Make sure directory %s is accessible!\n',fullfile(server2Name,'Data'));
end
if ~isdir(fullfile(server3Name,'Data'))
    fprintf('Make sure directory %s is accessible!\n',fullfile(server3Name,'Data'));
end


DIRS.data           = fullfile(serverName,'Data','trodes');
DIRS.spikes         = fullfile(serverName,'Data','Spikes');
DIRS.camera         = fullfile(server2Name,'Data','Camera');
DIRS.EyeCamera      = fullfile(server2Name,'Data','EyeCamera'); 
DIRS.EyeTrack       = fullfile(serverName,'Data','EyeTrack');
DIRS.xfiles         = fullfile(serverName,'Data','xfiles');
DIRS.michigan       = fullfile(serverName,'Data','michigan');
DIRS.Cerebus        = fullfile(serverName,'Data','Cerebus');
DIRS.stimInfo       = fullfile(serverName,'Data','stimInfo');
DIRS.behavior       = fullfile(serverName,'Data','behavior');
DIRS.mouselogs      = fullfile(serverName,'Data','logs','mouse','behavior');
DIRS.multichanspikes= fullfile(serverName,'Data','multichanspikes');
DIRS.ball           = fullfile(serverName,'Data','ball');
DIRS.Stacks         = fullfile(server2Name,'Data','Stacks');
DIRS.expInfo        = fullfile(serverName,'Data','expInfo');
