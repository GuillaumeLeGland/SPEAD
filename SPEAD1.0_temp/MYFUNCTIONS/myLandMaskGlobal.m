function [GLOBE,Land] = myLandMaskGlobal(varargin)
%**********************************************
%Use: [GLOBE,Land]=myLandMaskGlobal(varargin)
%**********************************************

mlat = 180;
nlon = 360;
%%keyMapOrientation = 'Default';
keyMapOrientation = 'Traditional'; %USE THIS ONE.

if length(varargin) == 2

    mlat = varargin{1};
    nlon = varargin{2};

elseif length(varargin) == 3 

    mlat = varargin{1};
    nlon = varargin{2};
    keyMapOrientation = varargin{3};

end

if mlat > 180 | nlon > 360
    error('mlat and nlon *must* be [180,360] maximum!!!')
end

load topo 

if strcmp(keyMapOrientation,'Default')
    
    MAP = topo;

elseif strcmp(keyMapOrientation,'Traditional')

    MAP = rot90(topo');
    MAP1 = MAP(:,1:180);
    MAP2 = MAP(:,181:360);
    MAP = [MAP2,MAP1]; %centro en el Atlantico.

end

if mlat==162 %el mapa va de 80N a 80S.
    MAP = MAP(9:170,:);
end

%.................
depthz0 = 0; %[m]
% $$$ depthz0 = -1; %[m]
%.................
Landfromatlab = find(MAP >= depthz0);
%.................
% $$$ GLOBE = zeros(mlat,nlon); %Zero in Ocean.
% $$$ GLOBE(Landfromatlab) = 1.0; %Ones in Land.
%.................
% $$$ GLOBE = ones(mlat,nlon); %Ones in Ocean.
% $$$ GLOBE(Landfromatlab) = 0.0; %Zero in Land.
%.................
GLOBE = ones(mlat,nlon); %Ones in Ocean.
GLOBE(Landfromatlab) = nan; %NANs in Land.
%.................
% $$$ cmap = jet;
% $$$ %%cmap = flipud(bone);
% $$$ figure(1)
% $$$ subplot(2,2,1)
% $$$ imagesc(MAP)
% $$$ colormap(cmap)
% $$$ mycolorbar('vertic')
% $$$ colormap(cmap)
% $$$ subplot(2,2,2)
% $$$ imagesc(GLOBE)
% $$$ mycolorbar('vertic')
% $$$ colormap(cmap)
% $$$ return
%.................

%%%%%%%%
%OUTPUT:
%%%%%%%%
Land = Landfromatlab;
