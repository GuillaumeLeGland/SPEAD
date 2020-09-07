function [mymapSeaWifs]=mycolormapSeawifs(nbins)

%..........................
mapHSV = flipud(hsv);
mapHOT = hot;
mapJET = jet;
%..........................
[msize,nsize]=size(mapJET);
%..........................
%==========================
%..........................
% $$$ jRowHSV = find(mapHSV(:,1)==0 & mapHSV(:,2)==1 & mapHSV(:,3)==1);
% $$$ jRowJET = find(mapJET(:,1)==0 & mapJET(:,2)==1 & mapJET(:,3)==1);
% $$$ Imaphsv = [16:jRowHSV];
% $$$ Imapjet = [jRowJET:msize];
%..........................
jRowHSV = find(mapHSV(:,1)==0 & mapHSV(:,2)==0.0625 & mapHSV(:,3)==1);
jRowJET = 1;
%..........................
Imaphsv = [16:jRowHSV-1];
Imapjet = [jRowJET:msize];
%..........................
MAP1 = mapHSV(Imaphsv,:);
MAP2 = mapJET(Imapjet,:);
%..........................
%==========================
%..........................
MAP = [MAP1;MAP2];
%..........................
MAP([1:6],3) = 0.5;
%..........................
[msize,nsize]=size(MAP);
%..........................
x = [1:msize];
xmin=min(x);
xmax=max(x);
deltax=(xmax-xmin)/(nbins-1);
xi = [xmin:deltax:xmax];
[MAPI]=interp1(x,MAP,xi);
%..........................
mymapSeaWifs = MAPI;
%..........................
