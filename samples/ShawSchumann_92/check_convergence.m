addpath ../../utilities

readinputs('../../')

auw=loadbin('../../output/auw.bin',Nz,'l');
atxz=loadbin('../../output/atxz.bin',Nz,'l');
N=size(auw,1);
mom_flux=squeeze(mean(auw+atxz,2));

figure;hold on;box on
plot(linspace(p_count,p_count*N,N),mom_flux)

