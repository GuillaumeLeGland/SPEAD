function [hp]=myPlot(x,y,myMarker)
    
%...............
[m,n]=size(x);
%...............
if m>1 & n>1 %(m,n) arrays.
    hp=[];
    %...............
    clear myMarker
    %...............
    myMarker{1}='o';
    myMarker{2}='s';
    myMarker{3}='p';
    myMarker{4}='x';
    myMarker{5}='*';
    myMarker{6}='+';
    %...............
    myMarkerSize{1}=[6];
    myMarkerSize{2}=[6];
    myMarkerSize{3}=[8];
    myMarkerSize{4}=[10];
    myMarkerSize{5}=[10];
    myMarkerSize{6}=[10];
    %...............
    for i=1:m
	%..................................................
	hpi=plot(repmat(x(i,:),[2 1]),repmat(y(i,:),[2 1]))
	%..................................................
	set(hpi,'Marker',myMarker{i})
	set(hpi,'MarkerSize',myMarkerSize{i})
	set(hpi,'MarkerEdgeColor',[0 0 0])
	%..................................................
	for j=1:n
	    myMarkerColorj=get(hpi(j),'Color')
	    set(hpi(j),'MarkerFaceColor',myMarkerColorj')
	end
	%..................................................
	hp=[hp,hpi];
	%..................................................
	hold on
    end
    hold off
else %vectors.
    hp=plot(x,y,myMarker)
end
%...............
