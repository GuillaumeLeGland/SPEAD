function [ZONAS,TERRA,MAPterra]=myoceanregions() %equivalent to [ZONAS,TERRA,MAP]=corrGLOBAL_byregion_regiones(Land).
%************************************************************
%Use: [ZONAS,TERRA,MAPterra]=myoceanregions()
%----------------------------------------------------------------------
%NOTA: Los valores siguientes (cada uno corresponde a un color) definen
%los 5 Oceanos:
%Artico=5;
%Atlantico=2;
%Pacifico=1;
%Indico=3;
%Southern=4;
%Mediterraneo=6;
%----------------------------------------------------------------------
%************************************************************
%..............................
% $$$ [GLOBE,Land]=myglobeland;
%..............................
load /home/svallina/SERVAL/SER24/PROGRAMMING/MATLAB/PROGRAMAS/CORRELACION2/GLOBAL/Land.txt
%..............................
file='/home/svallina/SERVAL/SER24/PROGRAMMING/MATLAB/PROGRAMAS/VARIOS/oceans.bmp';
[Ain,MAPin]=imread(file,'bmp');
Ain=double(Ain);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CORRIJO UN POCO LA "Ain" QUE ME DEFINE LOS OCEANOS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=Ain;
A(1:20,:)=5;
A(131:end,:)=4;
A(Land)=nan;
%..........................
% $$$ figure(1)
% $$$ subplot(1,2,1)
% $$$ imagesc(Ain)
% $$$ colorbar('horiz')
% $$$ subplot(1,2,2)
% $$$ imagesc(A)
% $$$ colorbar('horiz')
% $$$ pause(1)
% $$$ close all
%..........................

%%%%%%%%%%%%%%%%%%%%%%%
%PONGO NaN EN MI LAND:
%%%%%%%%%%%%%%%%%%%%%%%
A(Land)=nan;

%%%%%%%%%%%%%%%%%%%%%%
%OBTENGO LAS REGIONES:
%%%%%%%%%%%%%%%%%%%%%%
%===========
%Ecuatorial:
%===========
B=ones(180,360)*nan;
B(81:100,:)=999; %10N-10S.
%??????????????????????
% $$$ B=ones(180,360)*nan;
% $$$ B(86:95,:)=999;
% $$$ %??????????????????????
% $$$ B=ones(180,360)*nan; %10N-0N.
% $$$ B(81:90,:)=999;
%??????????????????????

EqAtl=find(B>0 & A==2); %Atlantico Equatorial (10S-10N).
EqInd=find(B>0 & A==3); %Indico Equatorial.
EqPac=find(B>0 & A==1); %Pacifico Equatorial.

%=================
%Subtropical H.N.:
%=================
B=ones(180,360)*nan;
B(51:80,:)=999;

SubAtlN=find(B>0 & A==2); %Atlantico Subtropical (10N-40N)
SubIndN=find(B>0 & A==3); %Indico Subtropical
SubPacN=find(B>0 & A==1); %Pacifico Subtropical

%=================
%Subtropical H.S.:
%=================
B=ones(180,360)*nan;
B(101:130,:)=999;

SubAtlS=find(B>0 & A==2); %Atlantico Subtropical (10S-40S)
SubIndS=find(B>0 & A==3); %Indico Subtropical
SubPacS=find(B>0 & A==1); %Pacifico Subtropical

%==============
%Subpolar H.N.:
%==============
B=ones(180,360)*nan;
B(21:50,:)=999;

NAtl=find(B>0 & A==2); %Atlantico Norte (40N-70N)
NPac=find(B>0 & A==1); %Pacifico Norte

%=========================
%Southern Ocean (40S-65S):
%=========================
SO=find(A==4); %Southern Ocean (40S-90S)
%ææææææææææææææææææææ
%Reduzco a (40S-60S):
%ææææææææææææææææææææ
B=ones(180,360)*nan;
SOlimit=150; %ultima lat del South.Ocean (60S)
B(SO)=1;
B(SOlimit:end,:)=nan;
SO=find(B==1);

%==========
%Antartico: (60S-90S)
%==========
Antarc=find(A==4);
B=ones(180,360)*nan;
B(Antarc)=1;
B(1:SOlimit-1,:)=nan;
Antarc=find(B==1);

%=============
%Mediterraneo:
%=============
Medit=find(A==6); %Mediterranean Sea.

%=======
%Artico:
%=======
B=ones(180,360)*nan;
B(1:20,:)=999;
Artic=find(B>0 & A==5); %Artico (70N-90N)

%%%%%%%%%%%%%%%%%%
%STOCKO LAS ZONAS:
%%%%%%%%%%%%%%%%%%
zona1=NAtl;
zona2=SubAtlN;
zona3=EqAtl;
zona4=SubAtlS;
zona5=NPac;
zona6=SubPacN;
zona7=EqPac;
zona8=SubPacS;
zona9=Medit;
zona10=SubIndN;
zona11=EqInd;
zona12=SubIndS;
zona13=Artic;
zona14=SO;
zona15=Antarc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOTEO LOS OCEANOS Y SUS REGIONES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ PAL=hsv;
% $$$ %PAL=colorcube;
% $$$ %MAP=PAL(1:64-8,:);
% $$$ MAP=PAL(1:64-26,:);
% $$$ m=size(MAP,1);
% $$$ dm=64/m;
% $$$ MAP1=interp1([1:dm:64],MAP(:,1),[1:64]);
% $$$ MAP2=interp1([1:dm:64],MAP(:,2),[1:64]);
% $$$ MAP3=interp1([1:dm:64],MAP(:,3),[1:64]);
% $$$ MAP=[MAP1',MAP2',MAP3'];
% $$$ MAP(end,:)=[1 1 1];

MAPterra=ones(14,3)*nan;
MAPterra(1,:)=[0.40,0.58,0.93]; %cornflower blue (N.Atl)
MAPterra(2,:)=[0,1,1]; %cyan (Sub.Atl.N)
MAPterra(3,:)=[0,0,0.8]; %medim blue (Eq.Atl)
MAPterra(4,:)=[0,0.54,0.54]; %cyan4 (Sub.Atl.S)
MAPterra(5,:)=[1,0.65,0]; %orange (N.Pac)
MAPterra(6,:)=[0.5,0.7,0.2]; %ligth green (Sub.Pacif.N)
MAPterra(7,:)=[0 1 0]; %green (Eq.Pac)
MAPterra(8,:)=[0.13,0.55,0.13]; %forest green (Sub.Pac.S)
MAPterra(9,:)=[1,0,0]; %red (Medit)
MAPterra(10,:)=[0.75,0,0.5]; %anil (Sub.Ind.N)
MAPterra(11,:)=[1,0.41,0.71]; %pink (Eq.Ind)
MAPterra(12,:)=[0.54,0.17,0.88]; %blue violet (Sub.Ind.S)
MAPterra(13,:)=[1 1 1]; %white (Artic)
MAPterra(14,:)=[1,1,0]; %yellow (SO)
MAPterra(15,:)=[0.5 0.5 0.5]; %gray (Antarctic)
MAPterra=[[0 0 0];MAPterra];

TERRA=ones(180,360)*(-1);
TERRA(zona1)=1-0.5;
TERRA(zona2)=2-0.5;
TERRA(zona3)=3-0.5;
TERRA(zona4)=4-0.5;
TERRA(zona5)=5-0.5;
TERRA(zona6)=6-0.5;
TERRA(zona7)=7-0.5;
TERRA(zona8)=8-0.5;
TERRA(zona9)=9-0.5;
TERRA(zona10)=10-0.5;
TERRA(zona11)=11-0.5;
TERRA(zona12)=12-0.5;
TERRA(zona13)=13-0.5;
TERRA(zona14)=14-0.5;
TERRA(zona15)=15-0.5;
TERRA(1,1)=15;

%.........................
% $$$ figure(1)
% $$$ MMAP=hsv;
% $$$ MMAP(1,:)=[0 0 0];
% $$$ AA=A;
% $$$ I=find(A==0);
% $$$ AA(I)=nan;
% $$$ imagesc(AA)
% $$$ colormap(MMAP)
% $$$ figure(2)
% $$$ colormap(MAPterra)
% $$$ imagesc(TERRA)
% $$$ colorbar('horiz')
% $$$ colormap(MAPterra)
% $$$ pause(1)
% $$$ close(1)
% $$$ %close(2)
%.........................

%.............................................
% $$$ [I1,J]=find(TERRA(:,1)==1-0.5);
% $$$ [I2,J]=find(TERRA(:,1)==2-0.5);
% $$$ [I3,J]=find(TERRA(:,1)==3-0.5);
% $$$ [I4,J]=find(TERRA(:,1)==4-0.5);
% $$$ [I5,J]=find(TERRA(:,1)==5-0.5);
% $$$ [I6,J]=find(TERRA(:,1)==6-0.5);
% $$$ [I7,J]=find(TERRA(:,1)==7-0.5);
% $$$ [I8,J]=find(TERRA(:,1)==8-0.5);
% $$$ [I9,J]=find(TERRA(:,1)==9-0.5);
% $$$ [I10,J]=find(TERRA(:,1)==10-0.5);
% $$$ [I11,J]=find(TERRA(:,1)==11-0.5);
% $$$ [I12,J]=find(TERRA(:,1)==12-0.5);
% $$$ [I13,J]=find(TERRA(:,1)==13-0.5);
% $$$ [I14,J]=find(TERRA(:,1)==14-0.5);
% $$$ [I15,J]=find(TERRA(:,1)==15-0.5);
% $$$ 
% $$$ I1
% $$$ I2
% $$$ I3
% $$$ I4
% $$$ I5
% $$$ I6
% $$$ I7
% $$$ I8
% $$$ I9
% $$$ I10
% $$$ I11
% $$$ I12
% $$$ I13
% $$$ I14
% $$$ I15
% $$$ 
% $$$ length(I1)
% $$$ length(I2)
% $$$ length(I3)
% $$$ length(I4)
% $$$ length(I5)
% $$$ length(I6)
% $$$ length(I7)
% $$$ length(I8)
% $$$ length(I9)
% $$$ length(I10)
% $$$ length(I11)
% $$$ length(I12)
% $$$ length(I13)
% $$$ length(I14)
% $$$ length(I15)
%.............................................
ZONAS=ones(15,180*360)*nan;
for i=1:15
    zonai=eval(['zona',num2str(i)]);
    n=length(zonai);
    ZONAS(i,1:n)=zonai';
end
