function graphinfo = ProtocolGetGraphInfo(protocol)
% ProtocolGetGraphInfo get info on what to plot for a given protocol
%
% graphinfo = ProtocolGetGraphInfo(protocol)
%
% part of Spikes
%
% 2000-09 MC wrote it
% 2012-02 MC added option to pool when only 2 active pars. cleaned up a bit.


graphinfo = [];
if isempty(protocol.activepars)
   msgbox('Cannot show the tuning for an experiment with no active parameters');
   return
end

nonblanks = setdiff(1:protocol.nstim,protocol.blankstims);

nvals = zeros(length(protocol.activepars),1); % was []
for ipar = 1:length(protocol.activepars)
   iactivepar = protocol.activepars{ipar}(1);
   nvals(ipar) = length(unique(protocol.pars(iactivepar,nonblanks)));
end

% see if the parameters describe a cross 
% if they do, the intersection is a point which which the others share either
% the abscissa or the ordinate.

if length(find(nvals>2))==2
   
   activepars = [protocol.activepars{nvals>2}];
   
   pp = [ protocol.pars(activepars(1),nonblanks)', ...
         protocol.pars(activepars(2),nonblanks)' ];
   
   % find the intersection
   crossvalues = [];
   for ipp = 1:size(pp,1)
      if all( pp(:,1) == pp(ipp,1) | pp(:,2) == pp(ipp,2) )
         crossvalues = pp(ipp,:);
         break
      end
   end 
   
   if ~isempty(crossvalues)
      
      % ax = [];
      for icol = 1:2
         
         othercol = rem(icol,2)+1;
         otherparvalue = crossvalues(othercol);
         
         graphinfo(icol).pos 	= find( pp(:,othercol) == crossvalues(othercol) );
         graphinfo(icol).xx 		= protocol.pars(activepars(icol),nonblanks);
         graphinfo(icol).xlabel	= [protocol.pardefs{activepars(icol)} ' [' protocol.parnames{activepars(icol)} ']' ];
         graphinfo(icol).title 	= [ protocol.parnames{activepars(othercol)} ' = ' num2str(otherparvalue) ];
         
      end
      
      return;
   end
end

if length(protocol.activepars) == 1
   yesdims = 1;
elseif length(protocol.activepars) > 1
   % pick which parameter to put on the x axis
   
   if length(find(nvals>2))==1
      % if only one parameter assumes more than 2 values
      yesdims(1) = find(nvals>2);
   else
      [ pick, ok ] = listdlg('ListString', makeparlist(protocol,nvals),...
         'SelectionMode','single',...
         'PromptString','Parameter for the abscissa:');
      if ~ok, graphinfo = []; return; end
      yesdims(1) = pick;
   end
   remdims = setdiff(1:length(protocol.activepars),yesdims);
end

% MC commented this out and inserted the next if on 2012-02-15
% if length(protocol.activepars) == 2
%    yesdims(2) = remdims;
% elseif length(protocol.activepars) > 2

if length(protocol.activepars) >= 2
        
        % pick which parameter to vary across row
   [ pick, ok ] = listdlg('ListString', makeparlist(protocol,nvals,remdims),...
      'Selectionmode','single',...
      'PromptString','to be varied across rows:',...
      'CancelString','pool them');
   if ok
      yesdims(2) = remdims(pick);
      remdims = setdiff(remdims,yesdims);
   end
end

if length(yesdims)>1 && length(protocol.activepars) == 3
   yesdims(3) = remdims(1);
elseif length(yesdims)>1 && length(protocol.activepars) > 3
   
   % pick which parameter to vary across columns
   [ pick, ok ] = listdlg('ListString', makeparlist(protocol,nvals,remdims),...
      'Selectionmode','single',...
      'PromptString','to be varied across columns:',...
      'CancelString','pool them');
   if ok
      yesdims(3) = remdims(pick);
      remdims = setdiff(remdims,yesdims); %#ok<NASGU>
   end
end

% g_dim is graphics dim, 1 to 3;
% d_dim is data dim

for g_dim = 1:length(yesdims)
   d_dim = yesdims(g_dim);
   if length(protocol.activepars{d_dim})==1
      ipar(g_dim) = protocol.activepars{d_dim};
      
   elseif length(protocol.activepars{d_dim})==2 && ...
         max(protocol.pars(protocol.activepars{d_dim}(1),nonblanks) == ...
         protocol.pars(protocol.activepars{d_dim}(2),nonblanks))
      % if two parameters, with identical values, pick the first
      ipar(g_dim) = protocol.activepars{d_dim}(1);
   else
      stringarray = cell(length(protocol.activepars{d_dim}),1);
      for ii = 1:length(protocol.activepars{d_dim})
         stringarray{ii} = protocol.parnames{protocol.activepars{d_dim}(ii)};
      end
      [ ii, ok ] = listdlg('ListString', stringarray,...
         'Selectionmode','single',...
         'PromptString','These vary together. Pick one:');
      if ~ok, graphinfo = []; return; end
      ipar(g_dim) = protocol.activepars{d_dim}(ii);
   end
end



switch length(ipar) 
   
case 1,	
   
   nrows = 1;
   ncols = 1;
   graphinfo(1,1).xx = protocol.pars(ipar(1),nonblanks);
   graphinfo(1,1).pos = 1:length(protocol.pars(ipar(1),nonblanks)); % MC 9.11.2000
   
case 2,
   ncols = 1;
   
   diffxpar2s = unique(protocol.pars(ipar(2),nonblanks));
   nrows = length(diffxpar2s);
   for irow = 1:nrows
      xpar2 = diffxpar2s(irow);      
      graphinfo(irow,1).xx = protocol.pars(ipar(1),nonblanks);
      graphinfo(irow,1).pos = find(protocol.pars(ipar(2),nonblanks) == xpar2);
   end
   
otherwise
   diffxpar2s = unique(protocol.pars(ipar(2),nonblanks));
   diffxpar3s = unique(protocol.pars(ipar(3),nonblanks));
   nrows = length(diffxpar2s);
   ncols = length(diffxpar3s);
   for irow = 1:nrows
      for icol = 1:ncols
         xpar2 = diffxpar2s(irow); 
         xpar3 = diffxpar3s(icol); 
         graphinfo(irow,icol).xx = protocol.pars(ipar(1),nonblanks);
         graphinfo(irow,icol).pos = ...
            find(	protocol.pars(ipar(2),nonblanks) == xpar2 & ...
            protocol.pars(ipar(3),nonblanks) == xpar3);
      end
   end			
end

if nrows>1
   graphinfo(1,ncols).comment = [ protocol.parnames{ipar(2)} ' =' ];
   
   for irow = 1:nrows
      graphinfo(irow,ncols).midrightcomment = num2str(diffxpar2s(irow));
   end
   
end

if ncols>1
   for icol = 1:ncols
      graphinfo(1,icol).title = [ protocol.parnames{ipar(3)} ' = ' num2str(diffxpar3s(icol)) ];
   end
end


graphinfo(nrows,1).xlabel = [protocol.pardefs{ipar(1)} ' [' protocol.parnames{ipar(1)} ']' ];



%--------------------------------------------------------------------------------

function parlist = makeparlist(protocol, nvals, dimlist)

if nargin<3, dimlist=1:length(protocol.activepars); end

parlist = cell(length(dimlist),1);
for dd = 1:length(dimlist)
   ii = dimlist(dd);
   parlist{dd} = [ num2str(nvals(ii)) ' ' protocol.parnames{protocol.activepars{ii}(1)} ];
   for jj = 2:length(protocol.activepars{ii})
      parlist{dd} = [ parlist{dd} ', ' protocol.parnames{protocol.activepars{ii}(jj)} ];
   end
end
