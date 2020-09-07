function [GLOBE,Landfromatlab]=myglobeland(varargin)

%[GLOBE,Landfromatlab]=myglobeland(varargin)

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
GLOBE=zeros(m,n);
GLOBE(Landfromatlab)=1;
%.................
% $$$ figure(1)
% $$$ imagesc(GLOBE)
% $$$ colorbar('horiz')
% $$$ pause(1)
%.................
return
