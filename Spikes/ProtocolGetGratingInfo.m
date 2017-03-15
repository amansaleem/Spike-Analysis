function gg = ProtocolGetGratingInfo(prot, graphflag)
% ProtocolGetGratingInfo gives information about gratings in a stimulus
%
% ProtocolGetGratingInfo with no arguments makes a figure with a cartoon of
% the stimulus (getting arguments from ExptPick).
%
% gg = ProtocolGetGratingInfo(prot)
% returns a structure array gg with fields
%
%     'c' 			% contrast, in percent
%     'sf' 		 	% spatial frequency, in cyc/deg
%     'tf' 			% temporal frequency, in Hz
%     'ori' 	 	% orientation, in deg
%     'idiam' 	 	% inner diameter, in deg
%     'odiam' 	 	% outer diameter, in deg
%     'x' 		 	% x center, in deg
%     'y'
%     'phase'
%
% each gg(istim) contains fields that have ngrats rows, where
% ngrats is 0 for blanks, 1 for gratings, 2 for plaids
%
% warning: duration is ignored
%
% gg = ProtocolGetGratingInfo(prot,'graphics') makes a cartoon of the
% stimuli
%
% EXAMPLE:
%
% global DIRS
% DIRS.data = 'Z:\cat'
% DIRS.spikes = 'S:\'
%
% prot = ProtocolLoad('catz005', 4, 5);
% gg = ProtocolGetGratingInfo(prot,'graphics')

% prot = ProtocolLoad('catz036', 18, 14);
% prot = ProtocolLoad('catz036', 18, 13);
% prot = ProtocolLoad('catz036', 20, 8);
% gg = ProtocolGetGratingInfo(prot,'graphics')

% 2002-05 Matteo Carandini
% 2003-10 MC added graphics
% 2003-12 VM added mean luminance
% 2004-05 AB added v2lut2gratdel and vlutgrat
% 2005-12 MC added v2lut2gratdelbar BUT CODE DOES NOT UNDERSTAND STIMULI THAT
% ARE NOT DISKS OR ANNULI
% 2008-04 LB added oglMovie2Grat.x and oglMovie2GratLin.x
% 2009-08 MC added case with no arguments, vdriftsin100.x, msgbox if unknown 
% 2013-05 MC changed the behavior for blank stimuli, now c = 0 rest is NaN
% 2013-05 MC replaced all strmatch with strncmp

%% sort out parameters

global PICK

if nargin<2
    graphflag = '';
end

if nargin<1
    prot = PICK.protocol;
    graphflag = 'graphics';
end

%% go for it

switch prot.xfile
    
    case 'visdriftsin.x'
        
        ngrats = 1;
        cs 		= prot.pars(strncmp('c',	prot.parnames,1),1:prot.nstim);
        sfs 	= prot.pars(strncmp('sf',	prot.parnames,2),1:prot.nstim)/10;
        tfs 	= prot.pars(strncmp('tf',	prot.parnames,2),1:prot.nstim)/10;
        oris 	= prot.pars(strncmp('ori',	prot.parnames,3),1:prot.nstim);
        odiams 	= prot.pars(strncmp('diam',	prot.parnames,4),1:prot.nstim)/10;
        xs 		= prot.pars(strncmp('x',	prot.parnames,1),1:prot.nstim)/10;
        ys 		= prot.pars(strncmp('y',	prot.parnames,1),1:prot.nstim)/10;
        idiams 	= zeros(ngrats, prot.nstim);
        phases	= zeros(ngrats, prot.nstim);   								% phase, in deg
        lmeans 	= zeros(ngrats, prot.nstim) + 0.5;
        
    case 'visdriftsin100.x'
        
        ngrats = 1;
        cs 		= prot.pars(strncmp('c',	prot.parnames,1),1:prot.nstim);
        sfs 	= prot.pars(strncmp('sf100',prot.parnames,5),1:prot.nstim)/100;
        tfs 	= prot.pars(strncmp('tf',	prot.parnames,2),1:prot.nstim)/10;
        oris 	= prot.pars(strncmp('ori',	prot.parnames,3),1:prot.nstim);
        odiams 	= prot.pars(strncmp('diam',	prot.parnames,4),1:prot.nstim)/10;
        xs 		= prot.pars(strncmp('x',	prot.parnames,1),1:prot.nstim)/10;
        ys 		= prot.pars(strncmp('y',	prot.parnames,1),1:prot.nstim)/10;	% y center, in deg
        idiams 	= zeros(ngrats, prot.nstim);
        phases	= zeros(ngrats, prot.nstim);   								% phase, in deg
        lmeans 	= zeros(ngrats, prot.nstim) + 0.5;
        
    case 'vdriftsin100.x'
        
        ngrats = 1;
        cs 		= prot.pars(strncmp('c',	prot.parnames,1),1:prot.nstim);
        sfs 	= prot.pars(strncmp('sf100',prot.parnames,5),1:prot.nstim)/100;
        tfs 	= prot.pars(strncmp('tf',	prot.parnames,2),1:prot.nstim)/10;
        oris 	= prot.pars(strncmp('ori',	prot.parnames,3),1:prot.nstim);
        odiams 	= prot.pars(strncmp('diam',	prot.parnames,4),1:prot.nstim)/10;
        xs 		= prot.pars(strncmp('x',	prot.parnames,1),1:prot.nstim)/10;
        ys 		= prot.pars(strncmp('y',	prot.parnames,1),1:prot.nstim)/10;	% y center, in deg
        idiams 	= zeros(ngrats, prot.nstim);
        phases	= zeros(ngrats, prot.nstim);   								% phase, in deg
        lmeans 	= zeros(ngrats, prot.nstim) + 0.5;
        
    case 'visplaid.x'
        
        ngrats = 2;
        cs 		= prot.pars(strncmp('c',            prot.parnames,1),1:prot.nstim);
        sfs 	= prot.pars(strncmp('sf',           prot.parnames,2),1:prot.nstim)/10;
        tfs 	= [1; 1]*prot.pars(strncmp('tf',    prot.parnames,2),1:prot.nstim)/10;
        oris 	= prot.pars(strncmp('ori',          prot.parnames,3),1:prot.nstim);
        odiams 	= [1; 1]*prot.pars(strncmp('diam',	prot.parnames,4),1:prot.nstim)/10;
        xs 		= [1; 1]*prot.pars(strncmp('x',		prot.parnames,1),1:prot.nstim)/10;
        ys 		= [1; 1]*prot.pars(strncmp('y',		prot.parnames,1),1:prot.nstim)/10;	% y center, in deg
        phases	= zeros(ngrats, prot.nstim);   	% phase, in deg
        idiams 	= zeros(ngrats, prot.nstim);
        lmeans 	= zeros(ngrats, prot.nstim) + 0.5;
        
    case 'vis2luts2grats.x'
        
        ngrats = 2;
        cs 		= prot.pars(strncmp('c',        prot.parnames,1),1:prot.nstim);
        sfs 		= prot.pars(strncmp('sf',   prot.parnames,2),1:prot.nstim)/10;
        tfs 		= prot.pars(strncmp('tf',   prot.parnames,2),1:prot.nstim)/10;
        oris     = prot.pars(strncmp('ori',     prot.parnames,3),1:prot.nstim);
        idiams 	= prot.pars(strncmp('idiam',	prot.parnames,5),1:prot.nstim)/10;
        odiams 	= prot.pars(strncmp('odiam',	prot.parnames,5),1:prot.nstim)/10;
        xs 		= prot.pars(strncmp('x',		prot.parnames,1),1:prot.nstim)/10;
        ys 		= prot.pars(strncmp('y',		prot.parnames,1),1:prot.nstim)/10;	% y center, in deg
        phases	= prot.pars(strncmp('ph',		prot.parnames,2),1:prot.nstim);   	% phase, in deg
        lmeans 	= zeros(ngrats, prot.nstim) + 0.5;
        
    case 'vis2luts2grats100.x'
        
        ngrats = 2;
        cs 		= prot.pars(strncmp('c',        prot.parnames,1),1:prot.nstim);
        sfs 		= prot.pars(strncmp('sf',   prot.parnames,2),1:prot.nstim)/100;
        tfs 		= prot.pars(strncmp('tf',   prot.parnames,2),1:prot.nstim)/10;
        oris 		= prot.pars(strncmp('ori',  prot.parnames,3),1:prot.nstim);
        idiams 	= prot.pars(strncmp('idiam',	prot.parnames,5),1:prot.nstim)/10;
        odiams 	= prot.pars(strncmp('odiam',	prot.parnames,5),1:prot.nstim)/10;
        xs 		= prot.pars(strncmp('x',		prot.parnames,1),1:prot.nstim)/10;
        ys 		= prot.pars(strncmp('y',		prot.parnames,1),1:prot.nstim)/10;	% y center, in deg
        phases	= prot.pars(strncmp('ph',		prot.parnames,2),1:prot.nstim);   	% phase, in deg
        lmeans 	= zeros(ngrats, prot.nstim) + 0.5;
        
    case 'vis2grat2win.x'
        
        ngrats = 2;
        cs 		= prot.pars(strncmp('c',        prot.parnames,1),1:prot.nstim);
        sfs 		= prot.pars(strncmp('sf',   prot.parnames,2),1:prot.nstim)/10;
        tfs 		= prot.pars(strncmp('tf',   prot.parnames,2),1:prot.nstim)/10;
        oris 		= prot.pars(strncmp('ori',  prot.parnames,3),1:prot.nstim);
        idiams 	= prot.pars(strncmp('idiam',	prot.parnames,5),1:prot.nstim)/10;
        odiams 	= prot.pars(strncmp('odiam',	prot.parnames,5),1:prot.nstim)/10;
        xs 		= [1; 1]*prot.pars(strncmp('x',	prot.parnames,1),1:prot.nstim)/10;
        ys 		= [1; 1]*prot.pars(strncmp('y',	prot.parnames,1),1:prot.nstim)/10;	% y center, in deg
        phases	= zeros(ngrats, prot.nstim);   								% phase, in deg
        lmeans 	= zeros(ngrats, prot.nstim) + 0.5;
        
    case 'vismovie2grat.x'
        
        ngrats = 2;
        cs 		= prot.pars(strncmp('c',        prot.parnames,1),1:prot.nstim);
        sfs 		= prot.pars(strncmp('sf',   prot.parnames,2),1:prot.nstim)/100;
        tfs 		= prot.pars(strncmp('tf',   prot.parnames,2),1:prot.nstim)/10;
        oris 		= prot.pars(strncmp('ori',  prot.parnames,3),1:prot.nstim);
        idiams 	= prot.pars(strncmp('idiam',	prot.parnames,5),1:prot.nstim)/10;
        odiams 	= prot.pars(strncmp('odiam',	prot.parnames,5),1:prot.nstim)/10;
        xs 		= [1; 1]*prot.pars(strncmp('x',	prot.parnames,1),1:prot.nstim)/10;
        ys 		= [1; 1]*prot.pars(strncmp('y',	prot.parnames,1),1:prot.nstim)/10;	% y center, in deg
        phases	= prot.pars(strncmp('ph',		prot.parnames,2),1:prot.nstim);		    	% phase, in deg
        lmeans 	= zeros(ngrats, prot.nstim) + 0.5;
        
    case {'vissweep2grat.x','vsweep2grat.x'}
        
        ngrats = 2;
        cs 		= prot.pars(strncmp('c',        prot.parnames,1),1:prot.nstim);
        sfs 	= prot.pars(strncmp('sf',       prot.parnames,2),1:prot.nstim)/100;
        tfs 		= prot.pars(strncmp('tfmin',prot.parnames,5),1:prot.nstim)/10;
        oris 		= prot.pars(strncmp('ori',  prot.parnames,3),1:prot.nstim);
        idiams 	= [0; 1]*prot.pars(strncmp('idiam',	prot.parnames,5),1:prot.nstim)/10;
        odiams 	= [prot.pars(strncmp('idiam',	prot.parnames,5),1:prot.nstim)/10; ...
            prot.pars(strncmp('odiam',          prot.parnames,5),1:prot.nstim)/10];
        xs 		= [1; 1]*prot.pars(strncmp('x',	prot.parnames,1),1:prot.nstim)/10;
        ys 		= [1; 1]*prot.pars(strncmp('y',	prot.parnames,1),1:prot.nstim)/10;	% y center, in deg
        phases	= prot.pars(strncmp('ph',		prot.parnames,2),1:prot.nstim);		    	% phase, in deg
        lmeans 	= prot.pars(strncmp('l',		prot.parnames,1),1:prot.nstim)/100;
        
    case 'v2lut2gratdelbar.x'
        
        ngrats = 2;
        cs 		= prot.pars(strncmp('c',    prot.parnames,1),1:prot.nstim);
        sfs 	= prot.pars(strncmp('sf100',prot.parnames,5),1:prot.nstim)/100;
        tfs 	= prot.pars(strncmp('tf',   prot.parnames,2),1:prot.nstim)/10;
        oris 	= prot.pars(strncmp('ori',  prot.parnames,3),1:prot.nstim);
        idiams 	= zeros(2,prot.nstim);
        odiams 	= zeros(2,prot.nstim)*NaN;
        xs 		= prot.pars(strncmp('x',	prot.parnames,1),1:prot.nstim)/10;
        ys 		= prot.pars(strncmp('y',	prot.parnames,1),1:prot.nstim)/10;	% y center, in deg
        phases	= prot.pars(strncmp('ph',	prot.parnames,2),1:prot.nstim);   	% phase, in deg
        lmeans 	= zeros(ngrats, prot.nstim) + 0.5;
        
    case 'v2lut2gratdel.x'
        
        ngrats = 2;
        cs 		= prot.pars(strncmp('c',    prot.parnames,1),1:prot.nstim);
        sfs 	= prot.pars(strncmp('sf100',prot.parnames,5),1:prot.nstim)/100;
        tfs 	= prot.pars(strncmp('tf',   prot.parnames,2),1:prot.nstim)/10;
        oris 		= prot.pars(strncmp('ori',prot.parnames,3),1:prot.nstim);
        idiams 	= prot.pars(strncmp('idiam',prot.parnames,5),1:prot.nstim)/10;
        odiams 	= prot.pars(strncmp('odiam',prot.parnames,5),1:prot.nstim)/10;
        xs 		= prot.pars(strncmp('x',	prot.parnames,1),1:prot.nstim)/10;
        ys 		= prot.pars(strncmp('y',	prot.parnames,1),1:prot.nstim)/10;	% y center, in deg
        phases	= prot.pars(strncmp('ph',	prot.parnames,2),1:prot.nstim);   	% phase, in deg
        lmeans 	= zeros(ngrats, prot.nstim) + 0.5;
        
    case 'vlutgrat.x'
        
        ngrats = 1;
        cs 		= prot.pars(strncmp('c',    prot.parnames,1),1:prot.nstim);
        sfs 		= prot.pars(strncmp('sf100',prot.parnames,5),1:prot.nstim)/100;
        tfs 		= prot.pars(strncmp('tf',prot.parnames,2),1:prot.nstim)/10;
        oris 		= prot.pars(strncmp('ori',prot.parnames,3),1:prot.nstim);
        odiams 	= prot.pars(strncmp('diam',	prot.parnames,4),1:prot.nstim)/10;
        xs 		= prot.pars(strncmp('x',	prot.parnames,1),1:prot.nstim)/10;
        ys 		= prot.pars(strncmp('y',	prot.parnames,1),1:prot.nstim)/10;	% y center, in deg
        phases	= prot.pars(strncmp('phase',prot.parnames,5),1:prot.nstim);   	% phase, in deg
        idiams 	= zeros(ngrats, prot.nstim);
        lmeans 	= zeros(ngrats, prot.nstim) + 0.5;
        
    case 'vmovie2gratinter.x'
        
        ngrats = 2;
        cs 		= prot.pars(strncmp('c',    prot.parnames,1),1:prot.nstim);
        sfs 	= prot.pars(strncmp('sf',   prot.parnames,2),1:prot.nstim)/100;
        tfs 	= prot.pars(strncmp('tf',   prot.parnames,2),1:prot.nstim)/10;
        oris 	= prot.pars(strncmp('ori',	prot.parnames,3),1:prot.nstim);
        xs 		= prot.pars(strncmp('x',    prot.parnames,1),1:prot.nstim)/10;
        ys 		= prot.pars(strncmp('y',	prot.parnames,1),1:prot.nstim)/10;	% y center, in deg
        phases	= prot.pars(strncmp('tph',	prot.parnames,3),1:prot.nstim);   	% phase, in deg  TEMPORAL
        lmeans 	= zeros(ngrats, prot.nstim) + 0.5;
        
        shapes = prot.pars(strncmp('shape',		prot.parnames,5),1:prot.nstim);
        idiams 	= zeros(ngrats, prot.nstim)*NaN;
        odiams 	= zeros(ngrats, prot.nstim)*NaN;
        for istim = 1:prot.nstim
            if all(shapes(:,istim)==0) % it is a disk
                idiams(:,istim) 	= prot.pars(strncmp('dima',	prot.parnames,4),istim)/10;
                odiams(:,istim) 	= prot.pars(strncmp('dimb',	prot.parnames,4),istim)/10;
            end
        end
        
    case {'vmovie2gratinterok.x', 'oglMovie2Grat.x', 'oglMovie2GratLin.x'}
        
        ngrats = 2;
        cs 		= prot.pars(strncmp('c',		prot.parnames,1),1:prot.nstim);
        sfs 	= prot.pars(strncmp('sf',		prot.parnames,2),1:prot.nstim)/100;
        tfs 	= prot.pars(strncmp('tf',		prot.parnames,2),1:prot.nstim)/10;
        oris 	= prot.pars(strncmp('ori',	    prot.parnames,3),1:prot.nstim);
        xs 		= prot.pars(strncmp('x',		prot.parnames,1),1:prot.nstim)/10;
        ys 		= prot.pars(strncmp('y',		prot.parnames,1),1:prot.nstim)/10;	% y center, in deg
        phases	= prot.pars(strncmp('tph',		prot.parnames,3),1:prot.nstim);   	% phase, in deg  TEMPORAL
        lmeans 	= zeros(ngrats, prot.nstim) + 0.5;
        
        shapes = prot.pars(strncmp('shape',		prot.parnames,5),1:prot.nstim);
        idiams 	= zeros(ngrats, prot.nstim)*NaN;
        odiams 	= zeros(ngrats, prot.nstim)*NaN;
        for istim = 1:prot.nstim
            if all(shapes(:,istim)==0) % it is a disk
                idiams(:,istim) 	= prot.pars(strncmp('dima',	prot.parnames,4),istim)/10;
                odiams(:,istim) 	= prot.pars(strncmp('dimb',	prot.parnames,4),istim)/10; 
            end
        end
        
        if all(isnan(idiams))
            error('Someone please write the code to deal with rectangles!!!');
        end
        
    otherwise
        
        ErrorMessage = sprintf('Do not know the x file %s -- might not be gratings\n',prot.xfile);
        disp(ErrorMessage);
        msgbox(ErrorMessage,'ProtocolGetGratingInfo');
        gg = [];
        return
        
end

% ---------------- in plaid experiments, figure out which stims are single gratings ------------

if ngrats == 2
    
    %-------- these are single gratings because they are identical
    
    gratstims = find( ....
        sfs(1,:)==sfs(2,:) & ...
        tfs(1,:)==tfs(2,:) & ...
        oris(1,:)==oris(2,:) & ...
        idiams(1,:)==idiams(2,:) & ...
        odiams(1,:)==odiams(2,:) & ...
        xs(1,:)==xs(2,:) & ...
        ys(1,:)==ys(2,:) & ...
        phases(1,:)==phases(2,:) & ...
        lmeans(1,:)==lmeans(2,:));
    
    cs(1,gratstims) = cs(1,gratstims)+cs(2,gratstims);
    cs(2,gratstims) = 0;
    
    %---these are single gratings because they abut (grating 1 inside)
    
    gratstims = find( ....
        cs(1,:)==cs(2,:) & ...
        sfs(1,:)==sfs(2,:) & ...
        tfs(1,:)==tfs(2,:) & ...
        oris(1,:)==oris(2,:) & ...
        idiams(2,:)==odiams(1,:) & ...
        xs(1,:)==xs(2,:) & ...
        ys(1,:)==ys(2,:) & ...
        phases(1,:)==phases(2,:) & ...
        lmeans(1,:)==lmeans(2,:));
    
    odiams(1,gratstims) = odiams(2,gratstims);
    cs(2,gratstims) = 0;
    
    %---these are single gratings because they abut (grating 2 inside)
    
    gratstims = find( ....
        cs(1,:)==cs(2,:) & ...
        sfs(1,:)==sfs(2,:) & ...
        tfs(1,:)==tfs(2,:) & ...
        oris(1,:)==oris(2,:) & ...
        idiams(1,:)==odiams(2,:) & ...
        xs(1,:)==xs(2,:) & ...
        ys(1,:)==ys(2,:) & ...
        phases(1,:)==phases(2,:) & ...
        lmeans(1,:)==lmeans(2,:));
    
    odiams(2,gratstims) = odiams(1,gratstims);
    cs(1,gratstims) = 0;
    
end

gg = struct( ...
    'c', 		NaN, ...
    'sf',		NaN, ...
    'tf',		NaN, ...
    'ori',		NaN, ...
    'idiam',	NaN, ...
    'odiam',	NaN, ...
    'x',		NaN, ...
    'y',		NaN, ...
    'phase',	NaN, ...
    'lmean',	NaN, ...
    'ngrats',   0);
    
gg = repmat(gg,[1 prot.nstim]);

for istim = 1:prot.nstim
    
    gratlist = find( cs(:, istim)>0 );
    
    if isempty(gratlist)
        gg(istim).c = 0; % it is a blank, other parameters stay NaN
    else
        
        gg(istim).c  	=  cs(gratlist, istim);
        gg(istim).sf 	= sfs(gratlist, istim);
        gg(istim).tf 	= tfs(gratlist, istim);
        gg(istim).ori 	= oris(gratlist, istim);
        gg(istim).idiam	= idiams(gratlist, istim);
        gg(istim).odiam	= odiams(gratlist, istim);
        gg(istim).x 	= xs(gratlist, istim);
        gg(istim).y 	= ys(gratlist, istim);
        gg(istim).phase  = phases(gratlist, istim);
        gg(istim).lmean  = lmeans(gratlist, istim);
        gg(istim).ngrats = length(gratlist);
    end
end

%% graphics 

if strcmp(graphflag,'graphics')
    
    nstim = length(gg);
    
    % --------------------- figure out bounds ------------------
    
    xlims = [ zeros(nstim,1)+Inf, zeros(nstim,1)-Inf];
    ylims = [ zeros(nstim,1)+Inf, zeros(nstim,1)-Inf];
    for istim = 1:nstim
        g = gg(istim);
        ngrats = length(g.c);
        if ngrats>0
            xlims(istim, :) = [min( g.x - g.odiam/2 ), max( g.x + g.odiam/2 )];
            ylims(istim, :) = [min( g.y - g.odiam/2 ), max( g.y + g.odiam/2 )];
        end
    end
    minx = min(xlims(:,1));
    maxx = max(xlims(:,2));
    miny = min(ylims(:,1));
    maxy = max(ylims(:,2));
    
    % --------------------- make figure ------------------
    
    figure; clf; ax = zeros(nstim,1);
    
    nrows = round(sqrt(nstim));
    ncols = ceil(nstim/nrows);
    
    npix = 128;
    
    xx = linspace(minx,maxx,npix);
    yy = linspace(miny,maxy,npix);
    [xxx, yyy] = meshgrid(xx, yy);
    
    for istim = 1:nstim
        g = gg(istim);
        ngrats = length(g.c);
        ll = zeros(npix,npix) + 0.5;
        cc = zeros(npix,npix);
        window = cell(ngrats,1);
        % The mean luminance
        for igrat = 1:ngrats
            rr2 = (xxx-g.x(igrat)).^2+(yyy-g.y(igrat)).^2;
            window{igrat} = ( rr2<=(g.odiam(igrat)/2)^2 & rr2>=(g.idiam(igrat)/2)^2 );
            ll(window{igrat}) = g.lmean(igrat);
        end
        % The modulation
        for igrat = 1:ngrats
            [wx, wy] = pol2cart( g.ori(igrat)/180*pi, g.sf(igrat) );
            cc = cc +  g.lmean(igrat)*g.c(igrat)/100*window{igrat}.* ...
                sin(2*pi*wx*xxx+2*pi*wy*yyy+pi*g.phase(igrat)/180);
        end
        % The full stimulus
        ii = ll + cc;
        
        ii(isnan(ii)) = 0.5; % to deal with blanks
        
        ax(istim) = subplot(nrows,ncols,istim);
        imagesc(ii,[0 1]);
        title(num2str(istim));
        
    end
    
    colormap gray
    set(ax,'plotboxaspectratio',[maxx-minx maxy-miny 1],'box','off');
    set(ax,'xtick',[],'xcolor','w'); % set(gca,'xtick',[1 npix],'xticklabel',[minx maxx]);
    set(ax,'ytick',[],'ycolor','w'); % set(gca,'ytick',[1 npix],'yticklabel',[miny maxy]);
    
end

