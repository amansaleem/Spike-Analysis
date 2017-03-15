function FigTuning_callbacks
% FigTuning_callbacks the switchyard for the FigTuning figure
%
% MC 2000-10

switch get(gcbo,'Tag')
   
case 'txtMaxResp'
   
   ax = findobj(gcbf,'userdata','scaleme!');                                                                 
   ytop =str2num( get(gcbo,'string'));                                                                      
   if isempty(ytop) |  ytop<=0 | isempty(ax), set(gcbo,'string',''); else set(ax,'ylim',[0 ytop]); end 
   
otherwise
   
   disp('I do not know this Tag');
   
end
