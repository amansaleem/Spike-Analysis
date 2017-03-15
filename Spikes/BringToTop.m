function BringToTop()
% Brings to the top the objects with the same UserData as the callback object

% 2007-10 Matteo Carandini

tag = get(gcbo,'UserData');
ax = get(gcbo,'Parent');

brothers = get(ax,'children');
twins = findobj(ax,'UserData',tag);
others = setdiff(brothers,twins);

set(ax,'children',[twins; others(end:-1:1)]);

return


%% Example of usage

figure; clf

pp(1) = plot( 1:10, (1:10).^2, 'r','linewidth',4 ); hold on
pp(2) = plot( 1:10, 100-(1:10).^2, 'g','linewidth',4 ); hold on

set(pp(1),'UserData',1);
set(pp(2),'UserData',2);

set(pp,'ButtonDownFcn','BringToTop')
