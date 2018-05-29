clear,clc

SA='../../';

%add utilities to the path temorarily
addpath ../../utilities

%read inputs from LESinputs.txt
readinputs(SA)

%averaging period
t_avg=[300000 nsteps]/p_count;

z_u=linspace(dz/2,l_z+dz/2,Nz)/l_z;
z_w=linspace(0,l_z,Nz)/l_z;

%mean velocity profile
au=loadbin('../../output/au.bin',Nz,'l');
Mau=mean(au(t_avg,:));
Ulog=1/0.4*log(z_u*l_z/0.1);
figure;box on;hold on;
plot(z_u,Mau,'ok')
plot(z_u,Ulog,'--k')
set(gca,'Xscale','log','Xlim',[10^-2 1])
title('Fig. 1 Porte-Agel etal (2000)')

%variances
u2=loadbin('../../output/u2.bin',Nz,'l');
Mu2=mean(u2(t_avg,:));
v2=loadbin('../../output/v2.bin',Nz,'l');
Mv2=mean(v2(t_avg,:));
w2=loadbin('../../output/w2.bin',Nz,'l');
Mw2=mean(w2(t_avg,:));

figure;box on;hold on;
p=get(gcf,'Position');
p(3)=p(3)*2;
set(gcf,'Position',p)
subplot(1,3,1);box on;
plot(Mu2,z_u,'-k')
xlabel('$\langle \tilde{u}_1^2 \rangle / u_*^2$','Interpreter','Latex')
ylabel('$z/H$','Interpreter','Latex')
set(gca,'Ylim',[0 1])
subplot(1,3,2);box on;
plot(Mv2,z_u,'-k')
xlabel('$\langle \tilde{u}_2^2 \rangle / u_*^2$','Interpreter','Latex')
set(gca,'Ylim',[0 1])
title('Fig. 3 Porte-Agel etal (2000)')
subplot(1,3,3);box on;
plot(Mw2,z_w,'-k')
xlabel('$\langle \tilde{u}_3^2 \rangle / u_*^2$','Interpreter', 'Latex')

%stress
auw=loadbin('../../output/auw.bin',Nz,'l');
Muw=mean(auw(t_avg,:));
atxz=loadbin('../../output/atxz.bin',Nz,'l');
Mtxz=mean(atxz(t_avg,:));
ustar=sqrt(abs(Mtxz(1)+Muw(1)))*u_star

figure;box on;hold on;
plot(Muw,z_w,'-.k')
plot(Mtxz,z_w,':k')
plot(Muw+Mtxz,z_w,'-k')
xlabel('Normalized stress')
ylabel('$z/H$','Interpreter','Latex')
title('Fig. 4 Porte-Agel etal (2000)')
legend('resolved','SGS','total','Location','NorthWest')

%dynamic coefficient
ScCs=loadbin('../../output/scalar1/cs2pr.bin',Nz,'l');
MScCs=mean(ScCs(t_avg,:));

figure;box on;hold on;
plot(MScCs,z_w,'-k')
xlabel('$Sc_{sgs}^{-1}C_s^2$','Interpreter','Latex')
ylabel('$z/H$','Interpreter','Latex')
set(gca,'Ylim',[0 0.5])
title('Fig. 1 Porte-Agel (2004)')

%scalar variance
theta_star=surfaceFluxes/u_star;
t2=loadbin('../../output/scalar1/t2.bin',Nz,'l');
Mt2=mean(t2(t_avg,:));

figure;box on;hold on;
plot(Mt2/theta_star^2,z_w,'-k')
xlabel('$\langle \theta^{\prime 2} \rangle \theta_*^{-2}$','Interpreter','Latex')
ylabel('$z/H$','Interpreter','Latex')
title('Fig. 5 Porte-Agel (2004)')