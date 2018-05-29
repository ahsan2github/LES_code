addpath ../../utilities

readinputs('../../')

au=loadbin('../../output/au.bin',Nz,'l');
auw=loadbin('../../output/auw.bin',Nz,'l')*u_star^2;
atxz=loadbin('../../output/atxz.bin',Nz,'l')*u_star^2;
%ustar=sqrt(abs(atxz(:,1)+auw(:,1)));
ustar=-sqrt(abs(mean(auw+atxz,2)));

figure;hold on
plot(linspace(p_count,p_count*length(ustar),length(ustar)),ustar)
%xlabel('time steps')
%ylabel('u_*^2')
%plot(auw+atxz,1:Nz)
%plot(au,1:Nz)
xax=get(gca,'Xlim');
plot(xax,[-u_star -u_star],'--k')

% $$$ for t=1:size(auw,1)
% $$$     title(['t=',num2str(t*p_count)])
% $$$ %    plot(au(t,:),1:Nz)
% $$$     plot(auw(t,:)+atxz(t,:),1:Nz)
% $$$     pause(0.1)
% $$$ end
