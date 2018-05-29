clear,clc

SA='../../';

%add utilities to the path temorarily
addpath ../../utilities

%read inputs from LESinputs.txt
readinputs(SA)

z_u=linspace(dz/2,l_z+dz/2,Nz);
z_w=linspace(0,l_z,Nz);


    
awt=loadbin('../../output/temperature/awt.bin',Nz,'l')*u_star*scalarScales;
Mwt=mean(awt(size(awt,1)/2:end,:));

figure;
plot(Mwt/surfaceFluxes,z_w/1600)
xlabel('$\langle w^\prime T^\prime \rangle/Q_s$','Interpreter', ...
       'Latex')
ylabel('$z/z_{i0}$','Interpreter','Latex')
title('Nieuwstadt et al (1991) Fig. 2')

aw2=loadbin('../../output/w2.bin',Nz,'l')*u_star^2;
Mw2=mean(aw2(size(aw2,1)/2:end,:));
figure;
ws0=(9.81/theta_0*surfaceFluxes*1600)^(1/3);
plot(Mw2/ws0^2,z_w/1600)
xlabel('$\langle w^{\prime2} \rangle/w_{*0}^2$','Interpreter', ...
       'Latex')
ylabel('$z/z_{i0}$','Interpreter','Latex')
title('Nieuwstadt et al (1991) Fig. 3')

au2=loadbin('../../output/u2.bin',Nz,'l')*u_star^2;
Mu2=mean(au2(size(au2,1)/2:end,:));
figure;
ws0=(9.81/theta_0*surfaceFluxes*1600)^(1/3);
plot(Mu2/ws0^2,z_w/1600)
xlabel('$\langle u^{\prime2} \rangle/w_{*0}^2$','Interpreter', ...
       'Latex')
ylabel('$z/z_{i0}$','Interpreter','Latex')
title('Nieuwstadt et al (1991) Fig. 4')

aT2=loadbin('../../output/temperature/t2.bin',Nz,'l')*u_star^2;
MT2=mean(aT2(size(aT2,1)/2:end,:));
figure;
Ts0=surfaceFluxes/(9.81/theta_0*surfaceFluxes*1600)^(1/3);
plot(MT2/Ts0^2,z_w/1600)
xlabel('$\langle T^{\prime2} \rangle/T_{*}^2$','Interpreter', ...
       'Latex')
ylabel('$z/z_{i0}$','Interpreter','Latex')
title('Nieuwstadt et al (1991) Fig. 5')
