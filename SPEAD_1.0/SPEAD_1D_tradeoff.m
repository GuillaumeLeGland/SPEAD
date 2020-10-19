% This function plots two trade-offs: gleaner-opportunist and optimal
% temperature for growth (Le Gland, 28/10/2019)
function [ ] = SPEAD_1D_tradeoff(mup0,amup,xaxis,temp0,Q10,yaxis,fignum)

% x-axis value of fx
%nstep = 301; % Number of DIN values considered
nstep = 300;
nx = length(xaxis); % Number of mu vs DIN curves
fx = zeros(nstep,nstep);
for iknp=1:nstep
    % logknp = -2 + (iknp-1)/100;
    % knp = 10^(logknp);
    knp = 0.01*iknp;
    mup = mup0 * knp^amup;
    for iconc=1:nstep
        % logconc = -2 + (iconc-1)/100;
        % conc = 10^(logconc);
        conc = 0.01*iconc;
        fx(iknp,iconc) = mup * conc / (knp + conc);
    end
end
%xcurv = 1 + 100*xaxis; % Indices of conc (or knp) in the fx vs knp (or conc) plots
%nutrcurv    = floor(100*(log10(xaxis)+2)+1);
%nutraxislog = -2:0.01:1;
%nutraxis    = 10.^(nutraxislog);
nutrcurv = floor(100*xaxis);
nutraxis = 0.01:0.01:3;
%xaxis = 0:0.01:2;
        

% y-axis values of fy from 0 to 35 degrees C
% There are several possible expressions for fy
ny  = length(yaxis);
fy  = zeros(191,191);
q10 = zeros(191,1);
for itopt=1:191
    topt = (itopt-1)/5;
    q10(itopt) = Q10.^((topt - temp0)/10);
    for itemp=1:191
        temp = (itemp-1)/5;
        if temp < topt + 5
            fy(itopt,itemp) = exp(0.2*(temp - topt)) .* (topt + 5 - temp)/5 .* q10(itopt);
        else
            fy(itopt,itemp) = 0;
        end
    end
end 
tcurv = floor(1 + 5*yaxis); % Indices of temp (or topt) in the fy vs topt (or temp) plots
taxis = 0:0.2:38;

% I-axis is added, but it is not related to a trait (Le Gland, 27/01/2020)
% Response to PAR is the same for all species
fi = zeros(1,126);
for ipar=1:126;
    fi(ipar) = (13/12)*exp(-(1/12)*log(1/13))*(1-exp(-log(13)*(ipar-1)/25))*exp(-log(13)*(ipar-1)/(25*12));
end

color = ['g','r','k','b','m'];
%color = ['g','r','k','b'];

figure(fignum)
    
subplot(4,2,1)
for i = 1:min(10,nx)
    plot(nutraxis,fx(nutrcurv(i),:),[color(1+mod(i-1,5)),'.-']);
    %semilogx(nutraxis,fx(nutrcurv(i),:),[color(i),'.-']);
    hold on
end
hold off  
grid on
set(gca,'ytick',[0:0.2:1.0])
set(gca,'Ylim',[0 1.0])
% Only valid if xaxis has 4 elements
legend({[repmat('K_N = ',[nx 1]),num2str(xaxis(:),3)]},'FontSize',8,'Location','NorthWest');
%legend({[['K_N = ';'K_N = ';'K_N = ';'K_N = '],num2str(xaxis(:)),[' mmol m^{-3}';' mmol m^{-3}';' mmol m^{-3}';' mmol m^{-3}';]]},'FontSize',15,'Location','NorthWest');    
xlabel('Nutrient (DIN) concentration [mmolN.m^{-3}]')
ylabel('Growth factor  \gamma_n (n.d.)')
title('a) Nutrient-dependent growth factor of different phenotypes')

subplot(4,2,3)
for i=1:min(10,nx)
    semilogx(nutraxis,fx(:,nutrcurv(i)),[color(1+mod(i-1,5)),'.-']);
    %semilogx(nutraxis,fx(:,nutrcurv(i)),[color(i),'.-']);
    hold on
end
hold off
grid on
xlim([0.01 3])
set(gca,'ytick',[0:0.2:1.0])
set(gca,'Ylim',[0 1.0])
% Only valid if xaxis has 4 elements
legend({[repmat('DIN = ',[nx 1]),num2str(xaxis(:))]},'FontSize',8,'Location','NorthWest');
%legend({[['DIN = ';'DIN = ';'DIN = ';'DIN = '],num2str(xaxis(:)),[' mmol m^{-3}';' mmol m^{-3}';' mmol m^{-3}';' mmol m^{-3}';]]},'FontSize',15,'Location','NorthWest'); 
xlabel('Nutrient half-saturation [mmolN.m^{-3}]')
ylabel('Growth factor  \gamma_n (n.d.)')
title('c) Nutrient-dependent growth factor at different nutrient concentrations')

subplot(4,2,5)
for i=1:min(10,nx)
    semilogx(nutraxis,fx(:,nutrcurv(i))/max(fx(:,nutrcurv(i))),[color(1+mod(i-1,5)),'.-']);
    %semilogx(nutraxis,fx(:,nutrcurv(i)),[color(i),'.-']);
    hold on
end
hold off
grid on
xlim([0.01 3])
set(gca,'ytick',[0:0.2:1.2])
set(gca,'Ylim',[0 1.2])
% Only valid if xaxis has 4 elements
legend({[repmat('DIN = ',[nx 1]),num2str(xaxis(:))]},'FontSize',8,'Location','NorthWest');
%legend({[['DIN = ';'DIN = ';'DIN = ';'DIN = '],num2str(xaxis(:)),[' mmol m^{-3}';' mmol m^{-3}';' mmol m^{-3}';' mmol m^{-3}';]]},'FontSize',15,'Location','NorthWest'); 
xlabel('Nutrient half-saturation [mmolN.m^{-3}]')
ylabel('Relative competitive ability  (n.d.)')
title('e) Relative competitive abilities at different nutrient concentrations')

subplot(4,2,2)
for i=1:min(13,ny)
    plot(taxis(51:end),fy(tcurv(i),51:end),[color(1+mod(i-1,5)),'.-']);
    hold on
end
hold off
grid on
set(gca,'xtick',[10:2:38])
xlim([10 38])
set(gca,'Ylim',[0 2.0])
legend({[repmat('Topt = ',[ny 1]),num2str(yaxis(:))]},'FontSize',8,'Location','NorthWest');
%legend({[['T_{opt} = ';'T_{opt} = ';'T_{opt} = ';'T_{opt} = ';'T_{opt} = '],num2str(yaxis(:))]},'FontSize',15,'Location','NorthWest');
xlabel('Temperature (Celsius)')
ylabel('Growth factor  \gamma_T (n.d.)')
title('b) Temperature-dependent growth factor of different phenotypes')

subplot(4,2,4)
for i=1:min(13,ny)
    plot(taxis(51:end),fy(51:end,tcurv(i)),[color(1+mod(i-1,5)),'.-']);
    hold on
end
hold off
grid on
set(gca,'xtick',[10:2:38])
xlim([10 38])
set(gca,'Ylim',[0 2.0])
legend({[repmat('T = ',[ny 1]),num2str(yaxis(:))]},'FontSize',8,'Location','NorthWest');
%legend({[['T = ';'T = ';'T = ';'T = ';'T = '],num2str(yaxis(:))]},'FontSize',15,'Location','NorthWest');
xlabel('Topt (Celsius)')
ylabel('Growth factor  \gamma_T (n.d.)')
title('d) Temperature-dependent growth factor at different temperatures')

subplot(4,2,6)
for i=1:min(13,ny)
    plot(taxis(51:end),fy(51:end,tcurv(i))/max(fy(:,tcurv(i))),[color(1+mod(i-1,5)),'.-']);
    hold on
end
hold off
grid on
set(gca,'xtick',[10:2:38])
xlim([10 38])
set(gca,'ytick',[0:0.2:1.2])
set(gca,'Ylim',[0 1.2])
legend({[repmat('T = ',[ny 1]),num2str(yaxis(:))]},'FontSize',8,'Location','NorthWest');
%legend({[['T = ';'T = ';'T = ';'T = ';'T = '],num2str(yaxis(:))]},'FontSize',15,'Location','NorthWest');
xlabel('Topt (Celsius)')
ylabel('Relative competitive ability (n.d.)')
title('f) Relative competitive abilities at different temperatures')

subplot(4,2,7)
plot([0:1:125],fi,'.-')
hold off
grid on
set(gca,'xtick',[0:25:125])
xlim([0 125])
set(gca,'ytick',[0:0.2:1.2])
set(gca,'Ylim',[0 1.2])
xlabel('PAR (W.m^{-2})')
ylabel('Growth factor (n.d.)')
title('g) PAR-dependent growth factor')


