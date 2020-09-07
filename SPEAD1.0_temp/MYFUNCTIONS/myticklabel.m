function ax2 = ticklabel(ax1)
% TICKLABEL  Shift the tick labels in the X axis
% TICKLABEL(AX) shifts current ticklabels to the right between
% tick marks in axis AX. It only works for the X axis.
% Returns the handler to the hidden axis with the centered ticklabels.
%
% Author M. Vichi (CMCC), from ideas posted on the Mathworks forum
%<http://groups.google.com/group/comp.soft-sys.matlab/browse_thread/thread/cf3b800149f23168/1f0cd24258f45c5d>
%<http://www.mathworks.com/matlabcentral/newsreader/view_thread/128731>

if nargin == 0 
    ax1 = gca
end

ax2 = axes('position',get(ax1,'position'));

%invert the order of the axes
c = get(gcf,'children');
set(gcf,'children',flipud(c))

xlim = get(ax1,'XLim');
xtick1 = get(ax1,'XTick');
delt = diff(xtick1);
xtick2 = [0.5*[delt delt(end)]+xtick1(1:end)];
xticklabels = get(ax1,'XTickLabel');
set(ax2,'Xlim',xlim,'XTick',xtick2,'XTicklabel',xticklabels);
ytick = get(ax1,'YTick');
yticklabels = get(ax1,'YTickLabel');
set(ax2,'Ylim',ylim,'YTick',ytick,'YTicklabel',yticklabels);
set(ax1,'XTickLabel','','Visible', 'on') 


