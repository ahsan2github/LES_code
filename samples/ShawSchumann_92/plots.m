clear,clc

SA='../../';

%add utilities to the path temorarily
addpath ../../utilities

%read inputs from LESinputs.txt
readinputs(SA)

z_u=linspace(dz/2,l_z+dz/2,Nz)/h_canopy;
z_w=linspace(0,l_z,Nz)/h_canopy;

SS_U=[0.2081 0.2003 0.1743 0.1455 0.1436 0.1657 0.215 0.2972 0.4454 0.7467 ...
    0.9639 1.0582 1.1524 1.2106 1.2626 1.3028 1.34 1.3742 1.4113 1.4364 1.4705...
    1.4926 1.5238 1.5457 1.562 1.5721 1.5793 1.5805 1.5875 1.6006];
    
au=loadbin('../../output/au.bin',Nz,'l')*u_star;
Mau=mean(au(size(au,1)/2:end,:));
figure;box on;hold on;
plot(Mau/mean(Mau),z_u)
plot(SS_U,z_u,'--b')
xax=get(gca,'Xlim');
plot(xax,[1 1],'--k')
xlabel('U')
ylabel('z/h')
set(gca,'Ylim',[0 3])
