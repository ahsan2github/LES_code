clear,clc

BENCHMARK=0;

SA='../../';

addpath ../../utilities;

readinputs(SA)

WC{1}=[7.2596 5.3477 1.7532 0.2901 0.1714;
    0.0252 0.0871 0.1559 0.2179 0.2951];
WC{2}=[6.0005 3.9656 2.215 1.4111 0.8437 0.3237 0.1353 0.0887;
       0.0382 0.0912 0.1662 0.2347 0.2961 0.3575 0.4301 0.4955];
WC{3}=[3.8098 2.9151 2.6159 1.8025 1.449 1.0416 0.6341 0.4162 0.1714 0.0617 0.0334;
       0.0344 0.0985 0.1675 0.2281 0.3043 0.3657 0.4344 0.4999 0.558 0.6274 0.686];
WC{4}=[2.4943 2.3334 2.0919 1.9579 1.5011 1.3134 1.0181 0.6689 0.4812 0.2396 0.1593 0.0791 0.0526;
       0.0294 0.0952 0.1681 0.2376 0.3067 0.3724 0.438 0.4998 0.5729 0.6349 0.7006 0.7629 0.8288];
WC{5}=[0.5274 0.5410 0.6605 0.6608 0.8865 0.8734 0.9666 1.0199 1.1262 1.2326 1.1930 1.5116 1.3526 1.1539 0.9418 0.2654 0.1065;
       0.0268 0.0891 0.1623 0.2283 0.3017 0.3639 0.4263 0.4996 0.5655 0.6352 0.701 0.7635 0.822 0.8879 0.9537 1.0194 1.0889];
WC{6}=[0.9584 1.1055 1.0389 0.919 1.0795 0.9727 1.0933 1.0668 1.1204 1.0806 1.0808 1.0009 0.8809 0.8811 0.8146 0.4141 0.0671;
      0.0324 0.0877 0.1645 0.2301 0.3001 0.362 0.432 0.5015 0.5712 0.6371 0.6994 0.7614 0.8271 0.893 0.9588 1.0238 1.089];


       
z_u=linspace(dz/2,l_z+dz/2,Nz);
z_w=linspace(0,l_z,Nz);

zi=1000; %TODO: actually calculate this
ws=(9.81*surfaceFluxes*zi/theta_0)^(1/3); %TODO: actually calculate this

%mean velocity
au=loadbin('../../output/au.bin',Nz,'l')*u_star;
Mau=mean(au(round(size(au,1)/2):end,:));
U=mean(Mau);

figure;box on
plot(Mau,z_u/zi,'-k')
axis([0 4 0 1.25])
xlabel('U (m/s)')
ylabel('z/z_i')
title('Weil et al (2004) Fig. 1a')

%velocity variances
u2=loadbin('../../output/u2.bin',Nz,'l')*u_star^2;
Mu2=mean(u2(round(size(u2,1)/2):end,:));
v2=loadbin('../../output/v2.bin',Nz,'l')*u_star^2;
Mv2=mean(v2(round(size(v2,1)/2):end,:));
w2=loadbin('../../output/w2.bin',Nz,'l')*u_star^2;
Mw2=mean(w2(round(size(w2,1)/2):end,:));

figure;box on;hold on
plot(Mu2/ws^2,z_u/zi,'--k')
plot(Mv2/ws^2,z_u/zi,':k')
plot(Mw2/ws^2,z_w/zi,'-k')
axis([0 0.5 0 1.25])
xlabel('$\sigma_u^2/w_*^2,\,\,\sigma_v^2/w_*^2,\,\,\sigma_w^2/w_*^2$','Interpreter','latex')
ylabel('$z/z_i$','interpreter','latex')
title('Weil et al (2004) Fig. 1b')
h=legend('$\sigma_u^2$','$\sigma_v^2$','$\sigma_w^2$');
set(h,'Interpreter','latex')

%Concentrations
X=[0.12 0.25 0.38 0.50 1.55 2.95];
DX=2*dx;
DZ=2*dz;
C=zeros(length(X),Nz*dz/DZ,2);
for model=0:1
    
    %read release positions
    fid=fopen([SA,'input/release_pos.ini'],'r','l');
    release_pos=fread(fid,inf,'double');
    release_pos=reshape(release_pos,length(release_pos)/3,3);
    fclose(fid);

    %read partilce release indices
    if(model==0)
        fid=fopen([SA,'output/part_frame/attribs_nomodel.bin'],'r','l');
    elseif(model==1)
        fid=fopen([SA,'output/part_frame/attribs_weiliso.bin'],'r','l');
    end
    release_inds=fread(fid,inf,'int');
    release_inds=release_inds(1:2:length(release_inds)-1);
    fclose(fid);
    
for t=1:nframe
    
    if(model==0)
        if( ~exist(sprintf('%s%s%04d%s',SA,'output/part_frame/part_frame_nomodel',t,'.bin'),'file') )
            error('%s%s%04d%s',SA,'output/part_frame/part_frame_nomodel',t,'.bin does not exist')
        end
        fidpart=fopen(sprintf('%s%s%04d%s',SA,'output/part_frame/part_frame_nomodel',t,'.bin'),'r','l');
    elseif(model==1)
        if( ~exist(sprintf('%s%s%04d%s',SA,'output/part_frame/part_frame_weiliso',t,'.bin'),'file') )
            error('%s%s%04d%s',SA,'output/part_frame/part_frame_weiliso',t,'.bin does not exist')
        end
        fidpart=fopen(sprintf('%s%s%04d%s',SA,'output/part_frame/part_frame_weiliso',t,'.bin'),'r','l');
    end
    
    %read particles
    ipart=fread(fidpart,1,'int')/8/3;
    particle=fread(fidpart,ipart*3,'double');
    particle=reshape(particle,ipart,3);
    particle(:,1:3)=particle(:,1:3)*z_i;
    fclose(fidpart);
    
    %subtract off initial x-position
    particle(:,1)=particle(:,1)-release_pos(release_inds(1:ipart),1);
    
    %calculate concentration
    for i=1:length(X)
        x=X(i)*U*zi/ws;
        xi=floor(x/DX)+1;
        
        for p=1:ipart
        
            x_cell=floor(particle(p,1)/DX)+1;
            z_cell=floor(particle(p,3)/DZ)+1;
        
            if(x_cell==xi)
                C(i,z_cell,model+1)=C(i,z_cell,model+1)+1;
            end
            
        end
    end
    
end
end

%normalize C
C=C/DX/DZ/nframe;
Q=npart/((nsteps-start_release)*dt);
C=C*U*zi/Q;

figure;hold on;box on;
set(gcf,'Units','inches','position',[2 2 10 5]);

z_c=DZ/2:DZ:l_z+DZ/2;
for i=1:length(X)
    subplot(2,length(X),i);hold on;box on;
    plot(squeeze(C(i,:,1)),z_c/zi,'-k')
    set(gca,'Ylim',[0 1.5])
    plot(WC{i}(1,:),WC{i}(2,:),'.k','MarkerSize',14)
    title(['X=',num2str(X(i))])
end
for i=1:length(X)
    subplot(2,length(X),length(X)+i);hold on;box on;
    plot(squeeze(C(i,:,2)),z_c/zi,'-k')
    set(gca,'Ylim',[0 1.5])
    plot(WC{i}(1,:),WC{i}(2,:),'.k','MarkerSize',14)
end
        
if(BENCHMARK)
    save benchmark.mat C Mau Mu2 Mv2 Mw2
end
