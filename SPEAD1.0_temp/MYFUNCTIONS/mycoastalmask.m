function [GLOBE,Landandlitoral]=mycoastalmask(varargin)
%**********************************************
%Use: [GLOBE,Landandlitoral]=mycoastalmask(varargin)
%**********************************************

if length(varargin)==2
    m=varargin{1};
    n=varargin{2};
else
    m=180;
    n=360;
end

load topo
W=topo;
W=rot90(topo');
W=[W(:,181:360),W(:,1:180)]; %centro en el Atlantico.
if m==162 %el mapa va de 80N a 80S.
    W=W(9:170,:);
end
Landfromatlab=find(W>=0);
Landandlitoral=find(W>-100);
GLOBE=zeros(m,n);
GLOBE(Landandlitoral)=0.5;
GLOBE(Landfromatlab)=1;
%.................
figure(1)
imagesc(GLOBE)
colorbar('horiz')
pause(1)
%.................
