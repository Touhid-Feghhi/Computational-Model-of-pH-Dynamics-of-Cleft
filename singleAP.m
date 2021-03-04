function [] = difrec()
clear all
clc;
%%%%%%%%%%%%%%%%%%%%%spatial & time resolution %%%%%%%%%%%%%%%%%%%%%%
dimension = 12*1*57; %
SR = dimension/12; %spatial resolution
tfn=1*10^6; %%%%End of simulation
tspan=[0 tfn]; % set time interval in micro seconds
%%%%%%%%%%%%%%%% reaction rates %%%%%%%%%%%%%%%%
%%%%%% bicarbonate and CA %%%%%%
EnAcc =7*10^2 ; % The effect of CA in acceleration of reaction ( Hydration of CO2)
k12= 4.1530*10^(10); % forward reaction rate of comsumption of protons [1/mM*us]
k21=10^(10); % bacward reaction rate (decomposing [H2CO3] to H+ and [HCO3-] [1/us]
k23=EnAcc*10.96*10^(-6); %forward rate of converting bicarbonic acid to carbon dioxide [1/us]
k32=EnAcc*0.0302*10^(-6);%0.0302*10^(-6); %backward rate of bicarbonic acid reaction (converting  carbon dioxide water to bicarbonic acid  [1/us])
%%%%%% phosphate %%%%%%
%k45= 1.4111*10^(8);
k45=1.2589*10^(8);%3.1623*10^(8) %forward rate of phosphate buffer consuming protons   [1/mM*us]
k54=10^(4); %backward rate of phosphate buffer ( decomposing HA to protons and A-) [1/us]
%%%%% Glutamate %%%%%
kgama=5.62*10^(-2); %Glutamte reaction ratios kgama = [Glu-][H+]/[Glu] 
kb = 2.14*10^(-7);  % kb = [Glu2-][H+]/[Glu-] 
%units for 10^3 1/s needs to be revised
k76= 10^(4); %compositiion of Glu
k67= k76*kgama; % decomposition of Glu in reaction Glu >< Glu- +H+
k98= 10^(4); %compositiion of Glu-
k89= k98*kb; % decomposition of Glu- in reaction Glu- >< Glu2- +H+
%%%%% ATP %%%%%
k1110 =10^(4); % production of HATP in raction ATP + H >< HATP
k1011=k1110*10^(3.467);% rate of consumption of ATP in reaction ATP + H >< HATP
k1312=10^(4); % rate of consumption of H2ATP in reaction HATP + H >< H2ATP
k1213=k1312*10^(1.004); %rate of consumption of HATP in reaction HATP + H >< H2ATP
%%%%%%%%%%%%%%%%%%%%%% Initial Concentrations %%%%%%%%%%%%%%%%%%%%%
initialph=7.2; % The initial pH which the cleft is just before the release of protons
svradius= 0.01342; %[um]
cleftwidth=0.013; %[um]
r=sqrt(((4/3)*(svradius)^3)/(cleftwidth));
vesicleinitialpH =5.5; % The initial pH which the synaptic vesicle is filled with glutamates.
releasedproton =10^(-vesicleinitialpH +3); %%%%diffusion or diffusion-reaction%%%% the density of protons
%%%%%%% Realsed Glutamate %%%%%%%
avogadro = 6.022*10^23;
Nglut=8000;
Totalglut=Nglut/(avogadro* 3.14*(4/3)*(svradius*10^(-6))^3)
releasedGlu =0.0533*Totalglut; %%%% concentration of Glu
releasedGluN =0.9467*Totalglut; %%%%concentration of GluN
releasedGlu2N =0*Totalglut; %%%%concentration of Glu2N
%%%%%% Released ATP %%%%%%
TotalATP=0;%100;%150;%2384;
releasedH2ATP =0.0280*TotalATP; %%%%concentration of H2ATP
releasedHATP =0.8773*TotalATP; %%%%concentration of HATP
releasedATP =0.0947*TotalATP; %%%%concentration of ATP
%%%%%%Bicarbonate & phosphate %%%%%%%
%protoninitial= 10^(-initialph +3 ) +releasedproton/2300; %%%reaction only%%%% % initial concentration of protons in cleft [mM] 
protoninitial= 10^(-initialph +3 );   % initial concentration of protons in cleft [mM]
bircbonateinitial=19; %initial value of bicarbonate buffer [mM]
carbonicacidinitial = (k12/k21)*(10^(-initialph +3 ))*bircbonateinitial; % initial value of carbon carbonicacid, here we assume it is in equilbrium [mM]
carbondioxideinitial = (k23/k32)*carbonicacidinitial;  % initial value of carbon dioxide [mM]
phosphateioninitial=0; %initial value of phosphate  ion  buffer  [A]_0 [mM]
phosphateinitial = phosphateioninitial *(10^(-initialph +3 ))*(k45/k54);  %initial value of phosphate buffer  [HA]_0  [mM]
%%%%% glutamate initial concentration in the cleft at pH 7.2%%%%%%
glut=0;%1.72;
Gluinitial =  glut*0.00113128;
Glui2Nnitial =glut*0.00339;
GluNinitial = glut-Gluinitial-Glui2Nnitial;
%%%%% ATP initial concentration in the cleft at pH 7.2%%%%%%
atp=0;%1.4;
HATPinitial = atp*0.156955;
H2ATPinitial =atp*1.2468*10^(-4);
ATPinitial =atp-H2ATPinitial-HATPinitial;
%%%%%%%%%%%%%%%% Steady state concentration in the cleft at pH 7.2%%%%%%%%%%%%%%%%
steadyph=7.2;
%steadyph= initialph; 
steadyGlu =Gluinitial; 
steadyGluN =GluNinitial; 
steadyGlu2N =Glui2Nnitial;
steadyATP =ATPinitial; 
steadyHATP =HATPinitial; 
steadyH2ATP =H2ATPinitial;
steadyproton= 10^(-steadyph +3 ) ;  
steadybircbonate=bircbonateinitial; 
steadycarbonicacid = (k12/k21)*(10^(-steadyph +3 ))*steadybircbonate; 
steadycarbondioxide = (k23/k32)*steadycarbonicacid; % 
steadyphosphateion=phosphateioninitial; 
steadyphosphate = steadyphosphateion *(10^(-steadyph +3 ))*(k45/k54); 
%%%%%%%%%%%%%%%% diffusion coefficients %%%%%%%%%%%%%%%%
gammma = (1/1.6)^2;
DH=gammma*8.69*10^(-3);%diffusion coefficient of protons [um^2/us]
DHCO3 = gammma*1.11*10^(-3); %diffusion coefficient of bicabonate [um^2/us]
DH2CO3 = gammma*1.11*10^(-3);%diffusion coefficient of bicarbonic acid [um^2/us]
DCO2 = gammma*1.71*10^(-3);%diffusion coefficient of carbon dioxide  [um^2/us]
DAH= gammma*1.56*10^(-3);%diffusion coefficient of phosphate  [um^2/us]
DA =  gammma*1.56*10^(3);%diffusion coefficient of phosphate ion  [um^2/us]%
DGlu =gammma*0.76*10^(-3); %diffusion coefficient of Glu[um^2/us]%
DATP = gammma*0.71*10^(-3); %diffusion coefficient of ATP[um^2/us]%
c0(:,:)= zeros(1,12*SR); % creating a zero matrix for initial condition
%%%%%%%%%%%%%%%% setting the initial values for the matrix of all the species %%%%%%%%%%%%%%%%
    for i= 1:SR
    c0(1,i) =protoninitial;
    c0(1,i+SR) = bircbonateinitial;
    c0(1,i+2*SR) = carbonicacidinitial;
    c0(1,i+3*SR) = carbondioxideinitial;
    c0(1,i+4*SR) = phosphateinitial;
    c0(1,i+5*SR) = phosphateioninitial;
    c0(1,i+6*SR) = Gluinitial;
    c0(1,i+7*SR) = GluNinitial;
    c0(1,i+8*SR) = Glui2Nnitial;
    c0(1,i+9*SR) = H2ATPinitial;
    c0(1,i+10*SR) = HATPinitial;
    c0(1,i+11*SR) = ATPinitial;
    end
%%%%%%%%%%%%%%%% ODE15s function %%%%%%%%%%%%%%%%
options = odeset('RelTol',1e-12,'AbsTol',1e-12,'Events',@myEvent); %%% For diffusion-reaction- or reaction only
%options = odeset('RelTol',1e-8,'AbsTol',1e-8,'Stats','on'); %%% For diffusion-reaction- or reaction only
tic,[t,c]=ode15s(@gener,tspan,c0,options);toc %%% For diffusion-reaction- or reaction only
%[t,c]=ode45(@gener,tspan,c0);
function [value, isterminal, direction] = myEvent(t, c)
%%%%Since PMCA can not work at pH higher than 8.8, we end simulation as soon as it hits 8.8%%%%
value      = (c(1) < 10^(-8.8+3));
isterminal = 1;   % Stop the integration
direction  = 0;
end
%%%%%%%%%%%%%%%% making the initila pH of release point= pH of vesicle, if smaller, then we know it's artifact %%%%%%%%%%%%%%%%
L= length(c(:,1));% Number of steps in time corresponding to protons
PH(:,:)= -log10(c(:,:))+3;
%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%
t;
row_header={'Time[us]','PH','[H+]','[HCO3-]','[H2CO3]','[CO2]','[Phosphate]','[Phosphate-ion]','[GLU]','[GLU-]','[GLU2-]','H2ATP','HATP','ATP'};
Datamiddle=[t,PH(:,SR-28),c(:,SR-28),c(:,2*SR-28),c(:,3*SR-28),c(:,4*SR-28),c(:,5*SR-28),c(:,6*SR-28),c(:,7*SR-28),c(:,8*SR-28),c(:,9*SR-28),c(:,10*SR-28),c(:,11*SR-28),c(:,12*SR-28     )];
xlswrite('datamid.xls', row_header,'Sheet1','A1')
xlswrite('datamid.xls', Datamiddle,'Sheet1','A2')
Datasource=[t,PH(:,2),c(:,2),c(:,1*SR+2),c(:,2*SR+2),c(:,3*SR+2),c(:,4*SR+2),c(:,5*SR+2),c(:,6*SR+2),c(:,7*SR+2),c(:,8*SR+2),c(:,9*SR+2),c(:,10*SR+2),c(:,11*SR+2)];
xlswrite('datasource.xls', row_header,'Sheet2','A1')
xlswrite('datasource.xls', Datasource,'Sheet2','A2')
%caldyn =[t,1000*pre(:),1000*post(:),1000*f(:),1000*JJ(:)]; %[uM]
%row_header={'Time[us]','pre','post','total','PMCA'};
%xlswrite('caldyn.xls', row_header,'Sheet3','A1')
%xlswrite('caldyn.xls', caldyn,'Sheet3','A2')
plot(t,PH(:,2))
hold
plot(t,PH(:,27))
%%%%%%%%%%%%%%%% defining the differential equation %%%%%%%%%%%%%%%%
function dcdt = gener(t,c)
h=0.45/(SR) ; % set constants %0.45 and 38 for SR
ratio=r/h;
L(:,:) = zeros(12*SR,12*SR);
for i =1:(SR )
  ep(i) = 1+(1/(2*i));
  en(i) = 1-(1/(2*i));
end
%%%%%%%%%%%%%%%% Diffusion %%%%%%%%%%%%%%%%
%This part is about the boundary condition for diffusion equation
for i =0:11
    L(1 + i*SR, 1+i*SR)=-4;%
    L(1 + i*SR, 2+i*SR)= 4;%
end
%%%%%%%%%%%%%%%%%%%% open-close boundary conditions%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:12
    L((i*SR),(i*SR)-1)=en(SR);
    %L((i*SR),(i*SR)) =-2+ep(SR); %close
    L((i*SR),(i*SR)) =-2;   % open
end
% This part is about polar coordinates!
for i = 2:SR -1
    L(i,i) = -2;
    L(i,i-1) =en(i);
    L(i,i+1) = ep(i);
end
for i = (SR)+2:(2*SR)-1
    L(i,i) = -2;
    L(i,i-1) =en(i-SR);
    L(i,i+1) = ep(i-SR);
end
for i = (2*SR)+2:(3*SR)-1
    L(i,i) = -2;
    L(i,i-1) =en(i-2*SR);
    L(i,i+1) = ep(i-2*SR);
end

for i = (3*SR)+2:4*SR-1
    L(i,i) = -2;
    L(i,i-1) =en(i-3*SR);
    L(i,i+1) = ep(i-3*SR);
end
for i = (4*SR)+2:5*SR-1
    L(i,i) = -2;
    L(i,i-1) =en(i-4*SR);
    L(i,i+1) = ep(i-4*SR);
end
for i = (5*SR)+2:6*SR-1
    L(i,i) = -2;
    L(i,i-1) =en(i-5*SR);
    L(i,i+1) = ep(i-5*SR);
end 
for i = (6*SR)+2:7*SR-1
    L(i,i) = -2;
    L(i,i-1) =en(i-6*SR);
    L(i,i+1) = ep(i-6*SR);
end 
for i = (7*SR)+2:8*SR-1
    L(i,i) = -2;
    L(i,i-1) =en(i-7*SR);
    L(i,i+1) = ep(i-7*SR);
end 
for i = (8*SR)+2:9*SR-1
    L(i,i) = -2;
    L(i,i-1) =en(i-8*SR);
    L(i,i+1) = ep(i-8*SR);
end 
for i = (9*SR)+2:10*SR-1
    L(i,i) = -2;
    L(i,i-1) =en(i-9*SR);
    L(i,i+1) = ep(i-9*SR);
end 
for i = (10*SR)+2:11*SR-1
    L(i,i) = -2;
    L(i,i-1) =en(i-10*SR);
    L(i,i+1) = ep(i-10*SR);
end 
for i = (11*SR)+2:12*SR-1
    L(i,i) = -2;
    L(i,i-1) =en(i-11*SR);
    L(i,i+1) = ep(i-11*SR);
end 
% Diffusion terms due to each element

for i=1:SR
    A(i,1) = DH*c(i);
    B(i,:) = 0;
end
B(SR,:) = ep((SR))*DH*steadyproton;
for i=SR+ 1:SR*(2)
    A(i,1) = DHCO3*c(i);
    B(i,:) = 0;
end
B(2*SR,:) = ep((SR))*DHCO3*steadybircbonate;
for i=SR*(2)+1:SR*(3)
    A(i,1) = DH2CO3*c(i);
    B(i,:) = 0;
end
B(3*SR,:) = ep((SR))*DH2CO3*steadycarbonicacid;
for i= SR*(3)+1:4*SR
    A(i,1) = DCO2*c(i);
    B(i,:) =0;
end
B(4*SR,:) = ep((SR))*DCO2*steadycarbondioxide;
for i= SR*(4)+1:5*SR
    A(i,1) = DAH*c(i);
    B(i,:) = 0;
end
B(5*SR,:) = ep((SR))*DAH*steadyphosphate;
for i= SR*(5)+1:6*SR
    A(i,1) = DA*c(i);
    B(i,:) = 0;
end
B(6*SR,:) = ep((SR))*DA*steadyphosphateion;
for i= SR*(6)+1:7*SR
    A(i,1) = DGlu*c(i);
    B(i,:) = 0;
end
B(7*SR,:) = ep((SR))*DGlu*steadyGlu;
for i= SR*(7)+1:8*SR
    A(i,1) = DGlu*c(i);
    B(i,:) = 0;
end
B(8*SR,:) = ep((SR))*DGlu*steadyGluN;
for i= SR*(8)+1:9*SR
    A(i,1) = DGlu*c(i);
    B(i,:) = 0;
end
B(9*SR,:) = ep((SR))*DGlu*steadyGlu2N;
for i= SR*(9)+1:10*SR
    A(i,1) = DATP*c(i);
    B(i,:) = 0;
end
B(10*SR,:) = ep((SR))*DATP*steadyH2ATP;
for i= SR*(10)+1:11*SR
    A(i,1) = DATP*c(i);
    B(i,:) = 0;
end
B(11*SR,:) = ep((SR))*DATP*steadyHATP;
for i= SR*(11)+1:12*SR
    A(i,1) = DATP*c(i);
    B(i,:) = 0;
end
B(12*SR,:) = ep((SR))*DATP*steadyATP;
%%%%%%%%%%%%%%%% Reaction terms %%%%%%%%%%%%%%%%
  % Here we create the reaction matrix by from the buffer reactions. 
  % Whenevr the index starts at non i=1, we have to shift them!
  % [H+] =c(i),,,,,[HCO3-]= c(i+SR),,,,,[H2CO3] = c(i+2SR),,,,,
  %[CO2]= c(i+3SR),,,,,[HA] = c(i+4SR),,,,,[A-]= c(i+5SR),,,,,
  % [Glu] = c(i+6SR),,,,,[Glu-] = c(i+7SR),,,,,,,[Glu2-] = c(i+8*SR)
  %[H2ATP]= c(i+9SR),,,,,[HATP] = c(i+10SR),,,,,,,[ATP] = c(i+11*SR)
  %d[H+]/dt = -k12[H+][HCO3-]+k21[H2CO3] -k45[H+][A-]+k54[HA]+ k67[Glu]
  %-k76[Glu-][H+] + k89[Glu-] - k98[Glu2-][H+] +
  %k1110[HATP]-k1011[H+][ATP]+ k1312[H2ATP]-k1213[H+][HATP]
  %d[HCO3-]/dt= -k12[H+][HCO3-]+k21[H2CO3]
  %d[H2CO3]/dt = k12[H+][HCo3-] - k21[H2CO3]-k23[H2CO3]+k32[CO2]
  %d[CO2]/dt = k23[H2CO3]- k32[CO2]
  %d[HA]/dt = k45[H+][A-]-k54[HA]
  %d[A-]/dt = -k45[H+][A-]+k54[HA]
  %d[Glu]/dt = -k67[Glu] + k76[H+][Glu-]
  %d[Glu-]/dt = k67[Glu] - k76[Glu-][H+] - k89[Glu-] + k98[Glu2-][H+]
  % d[Glu2-]/dt = k89[Glu-] - k98[Glu2-][H+]
  %d[ATP]/dt=k1110[HATP]-k1011[H+][ATP]
  %d[HATP]/dt=-k1110[HATP]+k1011[H+][ATP] - k1213[H2ATP]+k1312[H+][HATP]
  %d[H2ATP]/dt= k1213[H2ATP]-k1312[H+][HATP]
for i= 1:SR
   reactionmatrix(i,1)  = -k12*c(i)*c(i+SR)+ k21*c(i+2*SR)-k45*c(i)*c(i+5*SR)+ k54*c(i+4*SR)+ k67*c(i+6*SR) -k76*c(i+7*SR)*c(i) + k89*c(i+7*SR) - k98*c(i+8*SR)*c(i) + k1110*c(i+10*SR)-k1011*c(i)*c(i+11*SR)+ k1312*c(i+9*SR)- k1213*c(i)*c(i+10*SR);% d[H+]/dt
end
for i= SR +1 :(2)*SR 
    reactionmatrix(i,1)  = -k12*c(i-SR)*c(i+SR-SR)+ k21*c(i+2*SR-SR);% d[HCO3-]/dt 
end
for i= (2)*SR +1:(3)*SR
    reactionmatrix(i,1)  = k12*c(i-2*SR)*c(i+SR-2*SR) -  k21*c(i+2*SR-2*SR) - k23*c(i+2*SR-2*SR) + k32*c(i+3*SR-2*SR);  % d[H2CO3]/dt 
end

for i= (3)*SR +1:(4)*SR
    reactionmatrix(i,1)  = k23*c(i+2*SR-3*SR) - k32*c(i+3*SR-3*SR); % d[CO2]/dt
end
for i= (4)*SR +1:(5)*SR
    reactionmatrix(i,1)  =  k45*c(i-4*SR)*c(i+5*SR-4*SR) - k54*c(i+4*SR-4*SR) ;   % d[HA]/dt
end
for i= (5)*SR +1:(6)*SR
    reactionmatrix(i,1)  = -k45*c(i-5*SR)*c(i+5*SR-5*SR) + k54*c(i+4*SR-5*SR) ;   % d[A-]/dt
end
for i= (6)*SR +1:(7)*SR
    reactionmatrix(i,1)  = -k67*c(i+6*SR-6*SR) + k76*c(i-6*SR)*c(i+7*SR-6*SR) ;   % d[Glu]/dt
end
for i= (7)*SR +1:(8)*SR
    reactionmatrix(i,1)  = k67*c(i+6*SR -7*SR) - k76*c(i+7*SR- 7*SR)*c(i-7*SR) - k89*c(i+7*SR-7*SR) + k98*c(i+8*SR-7*SR)*c(i-7*SR) ;   % d[Glu-]/dt
end
for i= (8)*SR +1:(9)*SR
    reactionmatrix(i,1)  = k89*c(i+7*SR-8*SR) - k98*c(i+8*SR-8*SR)*c(i-8*SR) ;   % d[Glu2-]/dt
end
for i= (9)*SR +1:(10)*SR
    reactionmatrix(i,1)  = -k1312*c(i+9*SR-9*SR)+k1213*c(i-9*SR)*c(i+10*SR-9*SR); % d[H2ATP]/dt
end
for i= (10)*SR +1:(11)*SR
    reactionmatrix(i,1)  = k1011*c(i-10*SR)*c(i+11*SR-10*SR)-k1110*c(i+10*SR-10*SR)- k1213*c(i+10*SR-10*SR)*c(i-10*SR)+k1312*c(i+9*SR-10*SR); % d[HATP]/dt
end
 for i= (11)*SR +1:(12)*SR
    reactionmatrix(i,1)  =k1110*c(i+10*SR-11*SR)-k1011*c(i-11*SR)*c(i+11*SR-11*SR);% d[ATP]/dt
end
%%%%%%%%%%%%%%%% Release %%%%%%%%%%%%%%%%
s(:,:)=zeros(12*SR,1);
tau_sv =1;% [us]
ts= 10^(-4);
s(2,1)=heaviside(t-ts)*1*releasedproton*(1/tau_sv)*exp(-(t-ts)/tau_sv);
%%%% Glutamate Release %%%%%
s(2+6*SR,1)= 1*releasedGlu*heaviside(t-ts)*(1/tau_sv)*exp(-(t-ts)/tau_sv);
s(2+7*SR,1)= 1*releasedGluN*heaviside(t-ts)*(1/tau_sv)*exp(-(t-ts)/tau_sv);
s(2+8*SR,1)= 1*releasedGlu2N*heaviside(t-ts)*(1/tau_sv)*exp(-(t-ts)/tau_sv);
%%%%% ATP Release%%%%%
s(2+9*SR,1)= 1*releasedH2ATP*heaviside(t-ts)*(1/tau_sv)*exp(-(t-ts)/tau_sv);
s(2+10*SR,1)= 1*releasedHATP*heaviside(t-ts)*(1/tau_sv)*exp(-(t-ts)/tau_sv);
s(2+11*SR,1)= 1*releasedATP*heaviside(t-ts)*(1/tau_sv)*exp(-(t-ts)/tau_sv);
%%%%%%%%%%%%%%%% PMCA %%%%%%%%%%%%%%%%
avogadro = 6.022*10^23;
cvol = 3.14*(0.45)^2*13*10^(-6-6-9+3); %[L] cleft volume in liters
areaa=3.14*(0.45)^2*10^(-6-6+3); %[L] cleft area in liters
denpre =400; % # pmca/um^2
denpost =1700; % # pmca/um^2
npre= 3.14*(0.45)^2*denpre; % # pmca concentration
npost= 3.14*(0.45)^2*denpost; % # pmca concentration
ETpre = npre/(avogadro*cvol); %[M]
ETpost = npost/(avogadro*cvol); %[M]
tau1=50*10^3; %us
tau2=50*10^3; %us
base1 =0.05*10^(-6);% % [M] Ca microdomain initila peak
base2 = 0.05*10^(-6);%%*(heaviside(t - td)); % [M] Ca microdomain initila peak
peak1 =0.4*10^(-6);%+base1;%% [M] Ca microdomain initila peak
peak2 = 1.5*10^(-6);%+base2;% % [M] Ca microdomain initila peak
k1=1.5*10^(2); % 1/M.us
k2= 15*10^(-6); % 1/us
%k3 = 12*10^(-6); %1/us
k3 = 100*10^(-6); %1/us
%km =(k2+k3)/k1; %M
km=0.7*10^(-6);% uM
%vmax=6*2*10^(-6);
td1=0.5*10^3;%us
td2=1*10^3;%us
tp1=1.5*10^3;%us
tp2=39.3*10^3;%us
trise1=0.2*10^3;
trise2=8*10^3;
%pre= (peak1/(tp1-td1))*(heaviside(t-td1)-heaviside(t-tp1))*(t-td1)+ peak1*exp(-(t-tp1)/tau1)*heaviside(t-tp1);%[Ca2+] vs time
%post=(peak2/(tp2-td2))*(heaviside(t-td2)-heaviside(t-tp2))*(t-td2)+ peak2*exp(-(t-tp2)/tau2)*heaviside(t-tp2);%[Ca2+] vs time
%pre=(peak1/(exp((tp1-td1)/trise1)-1))*(heaviside(t-td1)-heaviside(t-tp1))*(exp((t-td1)/trise1)-1)+ peak1*exp(-(t-tp1)/tau1)*heaviside(t-tp1);%[Ca2+] vs time
%post=(peak2/(exp((tp2-td2)/trise2)-1))*(heaviside(t-td2)-heaviside(t-tp2))*(exp((t-td2)/trise2)-1)+ peak2*exp(-(t-tp2)/tau2)*heaviside(t-tp2);%[Ca2+] vs time
%pre =(peak1)*(tanh((t-td1)/(trise1)))*(heaviside(t-td1)-heaviside(t-tp1))+peak1*exp(-(t-tp1)/tau1)*heaviside(t-tp1);%[Ca2+] vs time
%post=(peak2)*(tanh((t-td2)/(trise2)))*(heaviside(t-td2)-heaviside(t-tp2))+peak2*exp(-(t-tp2)/tau2)*heaviside(t-tp2);%[Ca2+] vs time
%post=2*(peak2)*(1/(1+exp(-(t-td2)/trise2))-0.5)*(heaviside(t-td2)-heaviside(t-tp2))+peak2*exp(-(t-tp2)/tau2)*heaviside(t-tp2);
%pre=1.01*(peak1)*(1 - exp(-(t-td1)/trise1))*(heaviside(t-td1)-heaviside(t-tp1))+peak1*exp(-(t-tp1)/tau1)*heaviside(t-tp1);%[Ca2+]
%post=1.01*(peak2)*(1 - exp(-(t-td2)/trise2))*(heaviside(t-td2)-heaviside(t-tp2))+peak2*exp(-(t-tp2)/tau2)*heaviside(t-tp2);%[Ca2+] vs time
b1=5*10^2;
b2=10*10^3;
s1=0.5+0.5*tanh((t-tp1)/b1);%smoothing function
s2=0.5+0.5*tanh((t-tp2)/b2);%smoothing function
pre=(1-s1)*(peak1)*(1 - exp(-(t-td1)/trise1))*(heaviside(t-td1))+s1*peak1*exp(-(t-tp1)/tau1);%[Ca2+] vs time
%pre=2*peak1;%*(heaviside(t-td1))*(1 - exp(-(t-td1)/trise1)); %constant calcium for train
post=(0.988*(1-s2)*(peak2)*(1 - exp(-(t-td2)/trise2))*(heaviside(t-td2))+s2*peak2*exp(-(t-tp2)/tau2));%[Ca2+] vs time
Jbasepre=1000*2*(k3*ETpre*(base1))/(km+(base1)); 
Jbasepost=1000*2*(k3*ETpre*(base2))/(km+(base2)); 
Jpre =1000*2*(k3*ETpre*(pre))/(km+(pre));%- Jbasepre; %[mM/us] I muliplied it to 1000 to convert M to mM and *2 to get 2 protons for each Ca2+
Jpost=1000*2*(k3*ETpost*(post))/(km+(post));%-Jbasepost;
J = Jpre+Jpost; %total Ca2+ extruded by PMCA
PMCA(:,:)=zeros(12*SR,1);
for i=1:SR
    PMCA(i,1)= J;%*((i^2-(i-1)^2)/(SR-1)^2);
end

%diffusionmatrix(:,:) =(1/h^2)*L*A(:,:)+s(:,:);               % Diffusion term --- closed boundary
diffusionmatrix(:,:) =(1/h^2)*L*A(:,:)+ (1/h^2)*B(:,:)+s(:,:); %Diffusion term --- open boundary

%%%%%%%%%%%%%%%% Last defining equation  relating diffusion and reaction - extrusion terms %%%%%%%%%%%%%%%%
%dcdt =    reactionmatrix(:,:); % this option only solves if we had a uniform distribution of elements through reaction. you need to change the initial conditions at line #  
%dcdt = diffusionmatrix(:,:);
%dcdt = reactionmatrix(:,:) + diffusionmatrix(:,:);
%dcdt = diffusionmatrix(:,:) - PMCA(:,:);
dcdt = reactionmatrix(:,:) + diffusionmatrix(:,:)- PMCA(:,:);
%%dcdt = PMCA(:,:);
end
end   