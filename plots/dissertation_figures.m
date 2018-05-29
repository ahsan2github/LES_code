%Eckman Layer Figures
clear all; clc;
%whitebg('w');
% close all;
%***********************************************************
L_z     = 1500.0;
%  z_i     = 2*4*4*4000/(2*pi);
z_i = 636.626;
p_count = 25; 
ustar   = 0.45;
Ugal    = 0.0;
T_scale = 1;
Ug=0;
Vg      = 0;
f       = 1e-4;                         %Coriolis Force
vonk    = 0.4;
z0      = 0.1;                          %Surface roughness (m)
g       = 9.81;
g_hat   = g*z_i./(ustar.^2);
fgr     = 1;
%***********************************************************

cd '/home/ahsan/LEScode_ahsan_freeForm/output/'
SS ='';
SSs =''; %prefix for scalars empty means not present!
Sprint = '-depsc';
Tend = 220000; % NeutralChannel
Tstart    = Tend*3/10; %Andren's case averages over last 3/10th
dtr       = 0.44019; %physical time in seconds

T         = (p_count:p_count:Tend)*dtr;     %Physical Time (sec)
Nx        = 64; Ny = Nx; Nz = 64;
dx        = 2*pi/Nx; dy = 2*pi/Ny; dz = (L_z/z_i)/(Nz-1);
delta     = fgr*(dx*dy*dz)^(1/3);

opt_t = 1; %plot time statistics
opt_m = 1; %plot mean stats
opt_cs = 1;%plot sgs stats
opt_sp = 1;%plot spectra
lntyp={'-r';'--g';':b';'-.c';'-m';'--y';':k>';...
       '-.r^';':g*';':bs';'-.cs';'--m';':y>';':ks';...
       ':r>';':gs';':b*';':cs';':mo';'-.yo';'-k';'-r'};
%***********************************************************
zH_1 = ((0:1:Nz-1)')/(Nz-1); zH_1(1) = 0.5/(Nz-1);  %For dudz,dvdz,dtdz
zH_2 = ((0:1:(Nz-1))')/(Nz-1);                      %For w,uw,vw,wt
zH_3 = ((0.5:1:Nz)')/(Nz-1);                        %For u,v,T,Pr,Cs
z_1  = zH_1*L_z;
z_2  = zH_2*L_z;
z_3  = zH_3*L_z;
%***********************************************************
S2 = sprintf('%s%s',SS,'au.bin');            au  = loadbin(S2,Nz,'l');
S2 = sprintf('%s%s',SS,'av.bin');            av  = loadbin(S2,Nz,'l');
S2 = sprintf('%s%s',SS,'aw.bin');            aw  = loadbin(S2,Nz,'l');
S2 = sprintf('%s%s',SS,'auw.bin');           auw = loadbin(S2,Nz,'l');
S2 = sprintf('%s%s',SS,'avw.bin');           avw = loadbin(S2,Nz,'l');
S2 = sprintf('%s%s',SS,'atxz.bin');          tuw = loadbin(S2,Nz,'l');
S2 = sprintf('%s%s',SS,'atyz.bin');          tvw = loadbin(S2,Nz,'l');
S2 = sprintf('%s%s',SS,'dudz.bin');          dudz= loadbin(S2,Nz,'l');
S2 = sprintf('%s%s',SS,'dvdz.bin');          dvdz= loadbin(S2,Nz,'l');
S2 = sprintf('%s%s',SS,'u2.bin');             u2 = loadbin(S2,Nz,'l');
S2 = sprintf('%s%s',SS,'v2.bin');             v2 = loadbin(S2,Nz,'l');
S2 = sprintf('%s%s',SS,'w2.bin');             w2 = loadbin(S2,Nz,'l');
S2 = sprintf('%s%s',SS,'Cs_ALL.bin');         Cs = loadbin(S2,Nz,'l');
S2 = sprintf('%s%s',SS,'Cs2_ALL.bin');        Cs2= loadbin(S2,Nz,'l');
S2 = sprintf('%s%s',SS,'beta1.bin');          b1 = loadbin(S2,Nz,'l');
S2 = sprintf('%s%s',SS,'ESGS.bin');         ESGS = loadbin(S2,Nz,'l');

if(isempty(SSs) ~= 1)
    S2 = sprintf('%s%s%s',SS,SSs,'scalarMean.bin');  at  = loadbin(S2,Nz,'l');
    S2 = sprintf('%s%s%s',SS,SSs,'awt.bin');     awt = loadbin(S2,Nz,'l');
    S2 = sprintf('%s%s%s',SS,SSs,'aut.bin');     aut = loadbin(S2,Nz,'l');
    S2 = sprintf('%s%s%s',SS,SSs,'flux_t3.bin'); twt = loadbin(S2,Nz,'l');
    S2 = sprintf('%s%s%s',SS,SSs,'dsdz.bin');    dtdz= loadbin(S2,Nz,'l');
    S2 = sprintf('%s%s%s',SS,SSs,'t2.bin');       t2 = loadbin(S2,Nz,'l');
    S2 = sprintf('%s%s%s',SS,SSs,'prandtl.bin');  Pr = loadbin(S2,Nz,'l');
    S2 = sprintf('%s%s%s',SS,SSs,'beta2.bin');    b2 = loadbin(S2,Nz,'l');
    S2 = sprintf('%s%s%s',SS,SSs,'ET.bin');       ET = loadbin(S2,Nz,'l');
end

width = 2.5; height = width / 1.5; 
%***********************************************************
%TIME SERIES PLOTS
%Figure: ustar/Ug vs time***********************************
if opt_t == 1
    uwT = auw(:,1)+tuw(:,1);
    vwT = avw(:,1)+tvw(:,1);
    u_sT= (uwT.^2+vwT.^2).^0.25;
    u_sT= u_sT*ustar;
    hfig = figure(1); set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 width height]);
    hold on;plot(T*f,u_sT,'.-','LineWidth',1, 'MarkerEdgeColor','r'); grid on;
    xlabel('tf','FontSize',10, 'interpreter', 'tex');
    ylabel('u_* (m/s)','FontSize',10, 'interpreter', 'tex');
    set(gca,'FontSize',10,'XLim',[0 max(T*f)],'YLim',[0.3 0.6]);
    ylim([0.2 0.5]);
    SPR = sprintf('%s%s',SS,'GoMRI_ustar'); print(Sprint,SPR);
   
end
%Figure: Cu vs time************************************
if opt_t == 1
    uwT   = (ustar^2)*(auw(:,1)+tuw(:,1));
    avT   = av*ustar;
    Cu    = -(f./uwT).*(mean(avT,2)-Vg)*L_z;
    hfig = figure(2); set(hfig, 'PaperUnits', 'inches', 'PaperPosition', [0 0  width height]);
    hold on;plot(T*f,Cu,'k','LineWidth',1);grid on
    xlabel('tf','FontSize',10, 'interpreter', 'tex');
    ylabel('C_u','FontSize',10, 'interpreter', 'tex');
    set(gca,'FontSize',10,'XLim',[0 max(T*f)],'YLim',[0 2]);
    SPR = sprintf('%s%s',SS,'GoMRI_Cu'); print(Sprint,SPR);
end
%Figure: Cv vs time************************************
if opt_t == 1
    vwT   = (ustar^2)*(avw(:,1)+tvw(:,1));
    auT   = au*ustar;
    Cv    = +(f./vwT).*(mean(auT,2)-Ug)*L_z;
    hfig = figure(3); set(hfig, 'PaperUnits', 'inches', 'PaperPosition', [0 0  width height]);
    hold on;plot(T*f,Cv,'k','LineWidth',1);grid on
    xlabel('tf','FontSize',10, 'interpreter', 'tex');
    ylabel('C_v','FontSize',10, 'interpreter', 'tex');
    set(gca,'FontSize',10,'XLim',[0 max(T*f)],'YLim',[0 3]);
    SPR = sprintf('%s%s',SS,'GoMRI_Cv'); print(Sprint,SPR);
end

%***********************************************************
mauw     = mean(auw(floor(Tstart/p_count):floor(Tend/p_count),:),1)';
mavw     = mean(avw(floor(Tstart/p_count):floor(Tend/p_count),:),1)';
mau      = mean(au(floor(Tstart/p_count):floor(Tend/p_count),:),1)';
mav      = mean(av(floor(Tstart/p_count):floor(Tend/p_count),:),1)';
maw      = mean(aw(floor(Tstart/p_count):floor(Tend/p_count),:),1)';
mtuw     = mean(tuw(floor(Tstart/p_count):floor(Tend/p_count),:),1)';
mtvw     = mean(tvw(floor(Tstart/p_count):floor(Tend/p_count),:),1)';
mu2      = mean(u2(floor(Tstart/p_count):floor(Tend/p_count),:),1)';
mv2      = mean(v2(floor(Tstart/p_count):floor(Tend/p_count),:),1)';
mw2      = mean(w2(floor(Tstart/p_count):floor(Tend/p_count),:),1)';
mdudz    = mean(dudz(floor(Tstart/p_count):floor(Tend/p_count),:),1)';
mdvdz    = mean(dvdz(floor(Tstart/p_count):floor(Tend/p_count),:),1)';
mCs      = mean(Cs(floor(Tstart/p_count):floor(Tend/p_count),:),1)';
mb1      = mean(b1(floor(Tstart/p_count):floor(Tend/p_count),:),1)';
mESGS    = mean(ESGS(floor(Tstart/p_count):floor(Tend/p_count),:),1)';

if(isempty(SSs) ~= 1)
    mawt     = mean(awt(floor(Tstart/p_count):floor(Tend/p_count),:),1)';
    maut     = mean(aut(floor(Tstart/p_count):floor(Tend/p_count),:),1)';
    mat      = mean(at(floor(Tstart/p_count):floor(Tend/p_count),:),1)';
    mtwt     = mean(twt(floor(Tstart/p_count):floor(Tend/p_count),:),1)';
    mt2      = mean(t2(floor(Tstart/p_count):floor(Tend/p_count),:),1)';
    mdtdz    = mean(dtdz(floor(Tstart/p_count):floor(Tend/p_count),:),1)';
    %Converting Inf to Nan
    Indx     = find(Pr == Inf); Pr(Indx) = 0/0;
    mPr      = nanmean(Pr(floor(Tstart/p_count):floor(Tend/p_count),:))';
    mb2      = mean(b2(floor(Tstart/p_count):floor(Tend/p_count),:),1)';
    mET      = mean(ET(floor(Tstart/p_count):floor(Tend/p_count),:),1)';
end

%***********************************************************
muw      = mauw+mtuw;
mvw      = mavw+mtvw;
u_sL     = (muw.^2+mvw.^2).^0.25;
u_s      = u_sL(1);

Ind = find(u_sL.^2 <= 0.05*(u_sL(1).^2));
    a1 = zH_2(Ind(1))*L_z; b1 = u_sL(Ind(1))^2; a2 = zH_2(Ind(1)-1)*L_z; b2 = u_sL(Ind(1)-1)^2;
    a = a2 + (a1-a2)*(0.05*(u_sL(1)^2)-b2)/(b1-b2); mBLA=a/0.95; clear Ind;
    display(['boundary layer height=',num2str(mBLA)])

if(isempty(SSs) ~= 1)
    mwt      = mawt+mtwt;
    t_sL     = -mwt./u_sL;
    t_s      = t_sL(1);
end
%Figure: phiM vs zf/u**********************************
if opt_m == 1
    X  = sqrt(mdudz.^2+mdvdz.^2)./u_s;
%     Y  = (z_1*f)./(ustar*u_s);
    Y  = (z_1)./(mBLA);
    hfig = figure(4); set(hfig, 'PaperUnits', 'inches', 'PaperPosition', [0 0  width height]);
    hold on;plot(0.4*(z_1/z_i).*X,Y,'o-k'); h = line([1 1],[min(Y) max(Y)]); hold off;
    set(h,'Color','k','LineStyle','--','LineWidth',2);
    xlabel('\Phi_M','FontSize',10, 'interpreter', 'tex');
    ylabel('zf/u_*','FontSize',10, 'interpreter', 'tex');
    set(gca,'FontSize',10,'XLim',[0 2],'YLim',[0 0.3]);
    SPR = sprintf('%s%s',SS,'GoMRI_phiM'); print(Sprint,SPR);
end
%Figure: phiC vs zf/u**********************************
if (opt_m == 1 && isempty(SSs) ~= 1)
    X  = mdtdz./t_s;
%     Y  = (z_1*f)./(ustar*u_s);
    Y  = (z_1)./(mBLA);
    hfig = figure(5); set(hfig, 'PaperUnits', 'inches', 'PaperPosition', [0 0  width height]);
    hold on;plot(0.4*(z_1/z_i).*X,Y,'o-k'); h = line([0.74 0.74],[min(Y) max(Y)]); hold off;
    set(h,'Color','k','LineStyle','--','LineWidth',2);
    xlabel('\Phi_C4','FontSize',10, 'interpreter', 'tex');
    ylabel('zf/u_*','FontSize',10, 'interpreter', 'tex');
    set(gca,'FontSize',10,'XLim',[0 2],'YLim',[0 0.1]);
    SPR = sprintf('%s%s',SS,'GoMRI_phiC'); print(Sprint,SPR);
end
%Figure: uw+vw vs zf/u*********************************
if opt_m == 1
%     Y  = (z_2*f)./(ustar*u_s);
    Y  = (z_2)./(mBLA);
    hfig = figure(6); set(hfig, 'PaperUnits', 'inches', 'PaperPosition', [0 0  width height]);
    hold on;plot((muw+mvw)/(u_s^2),Y,'-k',(mtuw+mtvw)/(u_s^2),Y,':k','LineWidth',2);
    P  = [0.241905 0.791389 0.203571 0.117738];
    LG = legend('Total','SGS',2); set(LG,'FontSize',10,'Position',P); legend boxoff;
    xlabel('<uw+vw>/u^2_*','FontSize',10, 'interpreter', 'tex');
    ylabel('zf/u_*','FontSize',10, 'interpreter', 'tex');
    set(gca,'FontSize',10,'XLim',[-1.0 0.2]);
    SPR = sprintf('%s%s',SS,'GoMRI_uw'); print(Sprint,SPR);
end
%Figure: vw vs zf/u**********************************
if opt_m == 1
%     Y  = (z_2*f)./(ustar*u_s);
    Y  = (z_2)./(mBLA);
    hfig = figure(7); set(hfig, 'PaperUnits', 'inches', 'PaperPosition', [0 0  width height]);
    hold on;plot(mvw/(u_s^2),Y,'-k',mtvw/(u_s^2),Y,':k','LineWidth',2);
    P  = [0.241905 0.791389 0.203571 0.117738];
    LG = legend('Total','SGS',2); set(LG,'FontSize',10,'Position',P); legend boxoff;
    xlabel('<vw>/u^2_*','FontSize',10, 'interpreter', 'tex');
    ylabel('zf/u_*','FontSize',10, 'interpreter', 'tex');
    set(gca,'FontSize',10,'XLim',[-0.7 0.3]);
    SPR = sprintf('%s%s',SS,'GoMRI_vw'); print(Sprint,SPR);
end
%Figure: wc vs zf/u**********************************
if (opt_m == 1 && isempty(SSs) ~= 1)
%     Y  = (z_2*f)./(ustar*u_s);
    Y  = (z_2)./(mBLA);
    hfig = figure(8); set(hfig, 'PaperUnits', 'inches', 'PaperPosition', [0 0  width height]);
    hold on;plot(-mwt/(u_s*t_s),Y,'-k',-mtwt/(u_s*t_s),Y,':k','LineWidth',2);
    P = [0.657381 0.822103 0.135714 0.0870238];
    LG = legend('Total','SGS',1); set(LG,'FontSize',10,'PaperPosition',P); legend boxoff;
    xlabel('-<wc>/u_*c_*','FontSize',10, 'interpreter', 'tex');
    ylabel('zf/u_*','FontSize',10, 'interpreter', 'tex');
    set(gca,'FontSize',10,'XLim',[0 1]);
    SPR = sprintf('%s%s',SS,'GoMRI_wc'); print(Sprint,SPR);
end
%Figure: SigU,SigV,SigW, SigC vs zf/u*******************************
if opt_m == 1
    i           = 1:Nz-1;
    Esgs_uvp    = mESGS;                %For stable remember to multiply by sqrt(1-Rf)
    Esgs_w      = zeros(size(mESGS));
    Esgs_w(i+1) = 0.5*(Esgs_uvp(i)+Esgs_uvp(i+1));

    varU        = mu2+(2/3)*Esgs_uvp;
    varV        = mv2+(2/3)*Esgs_uvp;
    varW        = mw2+(2/3)*Esgs_w;

%     Y_uvp       = (z_3*f)./(ustar*u_s);
%     Y_w         = (z_2*f)./(ustar*u_s);

    Y_uvp       = (z_3)./(mBLA);
    Y_w         = (z_2)./(mBLA);

    hfig = figure(9); set(hfig, 'PaperUnits', 'inches', 'PaperPosition', [0 0  width height]);
    hold on;plot(varU/u_s^2,Y_uvp,'-k',(2/3)*Esgs_uvp/u_s^2,Y_uvp,':k','LineWidth',2);
    P  = [0.641905 0.791389 0.203571 0.117738];
    LG = legend('Total','SGS',2); set(LG,'FontSize',10,'Position',P); legend boxoff;
    xlabel('<u^2>/u^2_*','FontSize',10, 'interpreter', 'tex');
    ylabel('zf/u_*','FontSize',10, 'interpreter', 'tex');
    set(gca,'FontSize',10,'XLim',[0 8]);
    SPR = sprintf('%s%s',SS,'GoMRI_sigU'); print(Sprint,SPR);

    hfig = figure(10); set(hfig, 'PaperUnits', 'inches', 'PaperPosition', [0 0  width height]);
    hold on;plot(varV/u_s^2,Y_uvp,'-k',(2/3)*Esgs_uvp/u_s^2,Y_uvp,':k','LineWidth',2);
    P  = [0.641905 0.791389 0.203571 0.117738];
    LG = legend('Total','SGS',1); set(LG,'FontSize',10,'Position',P); legend boxoff;
    xlabel('<v^2>/u^2_*','FontSize',10, 'interpreter', 'tex');
    ylabel('zf/u_*','FontSize',10, 'interpreter', 'tex');
    set(gca,'FontSize',10,'XLim',[0 4.2]);
    SPR = sprintf('%s%s',SS,'GoMRI_sigV'); print(Sprint,SPR);

    hfig = figure(11); set(hfig, 'PaperUnits', 'inches', 'PaperPosition', [0 0  width height]);
    hold on;plot(varW/u_s^2,Y_w,'-k',(2/3)*Esgs_w/u_s^2,Y_w,':k','LineWidth',2);
    P  = [0.641905 0.791389 0.203571 0.117738];
    LG = legend('Total','SGS',2); set(LG,'FontSize',10,'Position',P); legend boxoff;
    xlabel('<w^2>/u^2_*','FontSize',10, 'interpreter', 'tex');
    ylabel('zf/u_*','FontSize',10, 'interpreter', 'tex');
    set(gca,'FontSize',10,'XLim',[0 3]);
    SPR = sprintf('%s%s',SS,'GoMRI_sigW'); print(Sprint,SPR);

    if (isempty(SSs) ~= 1)
        ET_w        = mET;
        ET_uvp      = zeros(size(mET));
        ET_uvp(i)   = 0.5*(ET_w(i)+ET_w(i+1));
        ET_uvp      = ET_uvp.^2./(0.49^2*Esgs_uvp);
        varT        = mt2+ET_uvp;

        hfig = figure(12); set(hfig, 'PaperUnits', 'inches', 'PaperPosition', [0 0  width height]);
        hold on;plot(varT/t_s^2,Y_uvp,'-k',ET_uvp/t_s^2,Y_uvp,':k','LineWidth',2);
        P  = [0.641905 0.791389 0.203571 0.117738];
        LG = legend('Total','SGS',1); set(LG,'FontSize',10,'PaperPosition',P); legend boxoff;
        xlabel('<c^2>/c^2_*','FontSize',10, 'interpreter', 'tex');
        ylabel('zf/u_*','FontSize',10, 'interpreter', 'tex');
        set(gca,'FontSize',10,'XLim',[0 8]);
        SPR = sprintf('%s%s',SS,'GoMRI_sigC'); print(Sprint,SPR);
    end
end
% %Figure: M/ustar and Skewness of w**********************************
% opt = 0;
% if opt == 1
%     S2  = sprintf('%s%s',SS,'vel_frames.bin');
%     fid = fopen(S2,'r','b');
%     ts  = 0;
%     it  = 0;
%     Skw3 = zeros(Nz,1); Skwu2 = zeros(Nz,1); Skwv2 = zeros(Nz,1); SkwE = zeros(Nz,1);
%     M = zeros(Nz,1);
%     while ts == 0
%         dum = fread(fid,1,'int');
%         U1  = fread(fid,Nx*Ny*Nz,'double');
%         V1  = fread(fid,Nx*Ny*Nz,'double');
%         W1  = fread(fid,Nx*Ny*Nz,'double');
%         ts  = isempty(W1);
%         if ts == 1
%             break;
%         else
%             for k=1:Nz
%                 U2(:,:,k)= reshape(U1((k-1)*Nx*Ny+1:k*Nx*Ny),Nx,Ny);
%                 V2(:,:,k)= reshape(V1((k-1)*Nx*Ny+1:k*Nx*Ny),Nx,Ny);
%                 W2(:,:,k)= reshape(W1((k-1)*Nx*Ny+1:k*Nx*Ny),Nx,Ny);
%             end
%             U = U2*ustar + Ugal; V = V2*ustar; W = W2*ustar;
%
%             M3D = sqrt(U.^2+V.^2);
%             M   = M + squeeze(mean(mean(M3D)));
%
%             U_w = zeros(size(U));
%             V_w = zeros(size(V));
%             for k=2:Nz
%                 U_w(:,:,k) = 0.5*(U(:,:,k)+U(:,:,k-1));
%                 V_w(:,:,k) = 0.5*(V(:,:,k)+V(:,:,k-1));
%             end
%
%             for k=1:Nz
%                 U2D_w = squeeze(U_w(:,:,k));
%                 V2D_w = squeeze(V_w(:,:,k));
%                 W2D_w = squeeze(W(:,:,k));
%                 U2D_w = U2D_w - mean(mean(U2D_w));
%                 V2D_w = V2D_w - mean(mean(V2D_w));
%                 W2D_w = W2D_w - mean(mean(W2D_w));
%                 Skw3(k)=Skw3(k) + mean(mean(W2D_w.^3))/(ustar*u_s)^3;
%                 Skwu2(k)=Skwu2(k) + mean(mean(W2D_w.*U2D_w.*U2D_w))/(ustar*u_s)^3;
%                 Skwv2(k)=Skwv2(k) + mean(mean(W2D_w.*V2D_w.*V2D_w))/(ustar*u_s)^3;
%                 SkwE(k)=SkwE(k) + 0.5*mean(mean(W2D_w.*(U2D_w.^2+V2D_w.^2+W2D_w.^2)))/(ustar*u_s)^3;
%            end
%             dum = fread(fid,1,'int');
%             it = it + 1;
%         end
%     end
%     Skw3 = Skw3/it;
%     Skwu2= Skwu2/it;
%     Skwv2= Skwv2/it;
%     SkwE = SkwE/it;
%     X = M/it/(u_s*ustar);
%     figure(1);
%     semilogx(zH_3,X,'ok',zH_3,(1/vonk)*log(z_3/z0),'--k');
%     xlabel('z/H','FontSize',10); ylabel('M/u_*','FontSize',10);
%     set(gca,'FontSize',14,'XLim',[1e-2 1]); axis square;
%     SPR = sprintf('%s%s',SS,'GoMRI_Mbyustar'); print(Sprint,SPR);
%
%     Y_w         = (z_2*f)./(ustar*u_s);
%     figure(2);
%     plot(Skw3,Y_w,'--k',Skwu2,Y_w,'-.k',Skwv2,Y_w,':k',SkwE,Y_w,'-k','LineWidth',2);
%     P  = [0.641905 0.791389 0.203571 0.117738];
%     LG = legend('<w^3>/u_*^3','<wu^2>/u_*^3','<wv^2>/u_*^3','<wE>/u_*^3',1);
%     set(LG,'FontSize',14,'Position',P); legend boxoff;
%     xlabel('Normalized Velocity Variance Fluxes','FontSize',10); ylabel('zf/u_*','FontSize',10);
%     set(gca,'FontSize',14,'XLim',[-0.4 1.0]); axis square;
%     SPR = sprintf('%s%s',SS,'GoMRI_skW'); print(Sprint,SPR);
%     fclose(fid);
% end
%Figure: SGS parameters*******************************************
if opt_cs == 1
    Y_uvp       = (z_3*f)./(ustar*u_s);

    hfig = figure(13); set(hfig, 'PaperUnits', 'inches', 'PaperPosition', [0 0  width height]);
     set(hfig, 'visible', 'off'); plot(mCs,Y_uvp,'-k','LineWidth',2);
    P  = [0.641905 0.791389 0.203571 0.117738];
    xlabel('C_S','FontSize',10, 'interpreter', 'tex');
    ylabel('zf/u_*','FontSize',10, 'interpreter', 'tex');
    set(gca,'FontSize',10,'XLim',[0 0.2]);
    SPR = sprintf('%s%s',SS,'GoMRI_Cs'); print(Sprint,SPR);


    hfig = figure(14); set(hfig, 'PaperUnits', 'inches', 'PaperPosition', [0 0  width height]);
    set(hfig, 'visible', 'off');plot(mb1,Y_uvp,'-k','LineWidth',2);
    P  = [0.241905 0.791389 0.203571 0.117738];
    xlabel('\beta','FontSize',10, 'interpreter', 'tex');
    ylabel('zf/u_*','FontSize',10, 'interpreter', 'tex');
    set(gca,'FontSize',10,'XLim',[0 2.5]);
    SPR = sprintf('%s%s',SS,'GoMRI_betam'); print(Sprint,SPR);

    if (isempty(SSs) ~= 1)
        hfig = figure(15); set(hfig, 'PaperUnits', 'inches', 'PaperPosition', [0 0  width height]);
        set(hfig, 'visible', 'off');lot(mPr,Y_uvp,'-k','LineWidth',2);
        P  = [0.641905 0.791389 0.203571 0.117738];
        xlabel('Pr_{SGS}','FontSize',10, 'interpreter', 'tex');
        ylabel('zf/u_*','FontSize',10, 'interpreter', 'tex');
        set(gca,'FontSize',10,'XLim',[0 1.5]);
        SPR = sprintf('%s%s',SS,'GoMRI_Pr'); print(Sprint,SPR);

        hfig = figure(16); set(hfig, 'PaperUnits', 'inches', 'PaperPosition', [0 0  width height]);
        set(ffig, 'visible', 'off');lot(mb2,Y_uvp,'-k','LineWidth',2);
        P  = [0.241905 0.791389 0.203571 0.117738];
        xlabel('\beta_\theta','FontSize',10, 'interpreter', 'tex');
        ylabel('zf/u_*','FontSize',10, 'interpreter', 'tex');
        set(gca,'FontSize',10,'XLim',[0 1.5]);
        SPR = sprintf('%s%s',SS,'GoMRI_betah'); print(Sprint,SPR);
    end
end

%Figure: Spectra*******************************************
if opt_sp == 1

    %%%%%%%%%Spectra - u %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S = loadbin('spectru.bin',Nz,'l');
    R = Nx/2+1;
    C = length(S)/R;
    S = reshape(S,R,C);

    fx      = ([1:(Nx/2+1)]-1)/(2*pi*z_i);
    %Make sure to change the height when computing w-spectra
    hfig = figure(17); 
    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 width height]);
    pos = get(gca, 'Position');
    set(gca, 'Position', [pos(1)*1.2,.15,pos(3),0.8])
    set(hfig, 'visible', 'on');
    
    for h = 4:3:3*Nz/4
        Sp  = S(:,h:Nz:C);
        Sp  = mean(Sp(:,floor(Tstart/p_count):floor(Tend/p_count)),2);
%         f   = 2*pi*fx*(h-0.5)*L_z/(Nz-1);
        f   = fx*mBLA;
%         Sp  = Sp/((h-0.5)/(Nz-1));
%        Sp  = fx'*z_i.*Sp/(mu2(h));
        Sp  = (fx').*Sp/(mu2(h));
        semilogx(f(1:Nx/(2*fgr)),Sp(1:Nx/(2*fgr))'); hold all; %Explicit Filtering
%        loglog(fx(1:Nx/(2*fgr)-10).*z_i,Sp(1:Nx/(2*fgr)-10)'); hold all; %Explicit Filtering
    end
    %set(gca,'FontSize',14,'XLim',[1e-2 1e-2],'YLim',[3e-3 6e-2]);
    set(gca,'FontSize', 10);
    %hold on;X1=1.5;X2=50;H=line([X1 X2],0.02*[X1^(-2/3) X2^(-2/3)]);set(H,'color',[0 0 0],'LineWidth',2);
%     hold on; H=line([1*0.85 1*0.85],[1e-2*0.1 1e-2*3.0]);set(H,'color',[0 0 0],'LineWidth',1); hold on;
%     H=line([6 6],[1e-2*0.1 1e-2]);set(H,'color',[0 0 0],'LineWidth',1); hold on;
    % plot -2/3 slope line
%     xx = logspace(-0.699,1,10); cc =log10(0.08)-(-2/3.0)*log10(0.2); yy =  xx.^(-2/3.0)*10^(cc);
%     semilogx(xx,yy, '-k', 'markerfacecolor', 'None', 'markeredgecolor', 'r', 'markersize', 5, 'LineWidth',1); hold on;
    xlabel('$k_1\Delta$','FontSize',10,'interpreter','latex');
    ylabel('$k_1E_u(k_1)/\sigma_u^2$','FontSize',10, 'interpreter', 'latex');
%     xlim([0,10^3]); ylim([0.00001, 1.0]);
    SPR = sprintf('%s%s',SS,'GoMRI_specU_premultiplied'); print(Sprint,SPR);
    
    
    % regular spectra    
    hfig = figure(); 
    mrk={'o','s','*','v','+','^','x','.','<','>'};
    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 width height]);
    pos = get(gca, 'Position');
    set(hfig, 'visible', 'on');
    S = loadbin('spectru.bin',Nz,'l');
    R = Nx/2+1;
    C = length(S)/R;
    S = reshape(S,R,C);
    fx      = ([1:(Nx/2+1)]-1)/(2*pi*z_i);
    n = 1;
    for h = 2:2:18        
        Sp  = S(:,h:Nz:C);
        Sp  = mean(Sp(:,floor(Tstart/p_count):floor(Tend/p_count)),2);
        f   = fx*((h-0.5)*dz*z_i);
%         Sp  = Sp/((h-0.5)/(Nz-1));
        Sp = Sp/(u_s^2)/((h-0.5)*dz*z_i);
        loglog(f(1:Nx/(2*fgr)),Sp(1:Nx/(2*fgr))', lntyp{n},'LineWidth',0.25,'markerfacecolor', 'None', 'markersize', 2);
        hold all; %Explicit Filtering
        n=n+1; 
    end
    xx = linspace(0.5,2,10); cc =log10(10^-4*0.25)-(-5/3.0)*log10(0.35); yy =  xx.^(-5/3.0)*10^(cc);
    loglog(xx,yy, '-k', 'markerfacecolor', 'None', 'markersize', 2, 'LineWidth',0.5); hold on;
    xx=xlabel('$K_{x}z$','FontSize',9,'interpreter','latex');
    ylabel('$E_{u}/(z\ u_{*}^{2})$','FontSize',9, 'interpreter', 'latex');
    xlim([0.00005,100]); ylim([5*10^-8, 0.06]); 
    h=annotation('textarrow',[.4,.65],[.45,.7],'String',' ','LineWidth',0.25,'HeadWidth',5,'HeadStyle','vback3');
    set(gca, 'XTickMode', 'manual');
    set(gca, 'YTickMode', 'manual');   
    set(gca,'Fontsize',8);
    SPR = sprintf('%s%s',SS,'specU_ug2.eps');  
    set(gcf, 'Color', 'w');
    print(Sprint,SPR);


    
    
    
    
    %%%%%%%%%Spectra - v %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S = loadbin('spectrv.bin',Nz,'l');
    R = Nx/2+1;
    C = length(S)/R;
    S = reshape(S,R,C);

    fx      = ([1:(Nx/2+1)]-1)/(2*pi*z_i);

    %Make sure to change the height when computing w-spectra
    hfig = figure(18); set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 width height]);
    pos = get(gca, 'Position');
    set(gca, 'Position', [pos(1)*1.2,.15,pos(3),0.8])
    for h = 4:3:1*Nz/2
        Sp  = S(:,h:Nz:C);
        Sp  = mean(Sp(:,floor(Tstart/p_count):floor(Tend/p_count)),2);
%         Sp  = Sp/((h-0.5)/(Nz-1));
%         f   = 2*pi*fx*(h-0.5)*L_z/(Nz-1);
%        Sp  = fx'*z_i.*Sp/(mv2(h));
        Sp  = (fx'*z_i).^(1).*Sp/(mv2(h));
        f   = fx*mBLA;
        loglog(f(1:Nx/(2*fgr)-1),Sp(1:Nx/(2*fgr)-1)); hold all; %Explicit Filtering
    end
    set(gca,'FontSize', 10);
    %hold on;X1=1.5;X2=50;H=line([X1 X2],0.02*[X1^(-2/3) X2^(-2/3)]);set(H,'color',[0 0 0],'LineWidth',2);
    hold on; H=line([1*0.85 1*0.85],[1e-3*0.1 1e-2*3.0]);set(H,'color',[0 0 0],'LineWidth',1); hold on;
    H=line([3 3],[1e-3*0.1 1e-2]);set(H,'color',[0 0 0],'LineWidth',1); hold on;
    % plot -2/3 slope line
    xx = logspace(-0.699,1,10); cc =log10(0.08)-(-2/3.0)*log10(0.2); yy =  xx.^(-2/3.0)*10^(cc);
    loglog(xx,yy, '-k', 'markerfacecolor', 'None', 'markeredgecolor', 'r', 'markersize', 5, 'LineWidth',1); hold on;
    hold off;
    xlabel('k_1\delta','FontSize',10,'interpreter','tex');
    ylabel('k_1E_v(k_1)/\sigma_v^2','FontSize',10, 'interpreter', 'tex');
    xlim([0.3,40]); ylim([0.001, 0.10]);
    SPR = sprintf('%s%s',SS,'GoMRI_specV'); print(Sprint,SPR);

    %%%%%%%%%Spectra - w %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S = loadbin('spectrw.bin',Nz,'l');
    R = Nx/2+1;
    C = length(S)/R;
    S = reshape(S,R,C);

    fx      = ([1:(Nx/2+1)]-1)/(2*pi*z_i);

    %Make sure to change the height when computing w-spectra
    hfig = figure(19); set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0  width height]);
    pos = get(gca, 'Position');
    set(gca, 'Position', [pos(1)*1.2,.15,pos(3),0.8])
    for h = 4:3:1*Nz/2
        Sp  = S(:,h:Nz:C);
        Sp  = mean(Sp(:,floor(Tstart/p_count):floor(Tend/p_count)),2);
%         Sp  = Sp/((h-0.5)/(Nz-1));
%        Sp  = fx'*z_i.*Sp/(mw2(h));
        Sp  = (fx'*z_i).^(1).*Sp/(mw2(h));
%         f   = 2*pi*fx*(h-0.5)*L_z/(Nz-1);
        f   = fx*mBLA;
        loglog(f(1:Nx/(2*fgr)-1),Sp(1:Nx/(2*fgr)-1)); hold all; %Explicit Filtering
    end
    set(gca,'FontSize', 10);
    %hold on;X1=1.5;X2=50;H=line([X1 X2],0.02*[X1^(-2/3) X2^(-2/3)]);set(H,'color',[0 0 0],'LineWidth',2);
    hold on; H=line([1.7 1.7],[1e-3*0.1 1e-2*3.0]);set(H,'color',[0 0 0],'LineWidth',1); hold on;
    H=line([3 3],[1e-3*0.1 1e-2]);set(H,'color',[0 0 0],'LineWidth',1); hold on;
    % plot -2/3 slope line
    xx = logspace(-0.2,1,10); cc =log10(0.0536)-(-2/3.0)*log10(1.655); yy =  xx.^(-2/3.0)*10^(cc);
    loglog(xx,yy, '-k', 'markerfacecolor', 'None', 'markeredgecolor', 'r', 'markersize', 5, 'LineWidth',1); hold on;
    hold off;
    xlabel('k_1\Delta','FontSize',10,'interpreter','tex');
    ylabel('k_1E_w(k_1)/\sigma_w^2','FontSize',10, 'interpreter', 'tex');
    xlim([0.3,40]); ylim([0.001, 0.10]);
    SPR = sprintf('%s%s',SS,'GoMRI_specW'); print(Sprint,SPR); 

    %%%%%%%%%Spectra - T %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (isempty(SSs) ~= 1)
        cd(SSs);
        S = loadbin('spectrt.bin',Nz,'l');
        cd ..;
        R = Nx/2+1;
        C = length(S)/R;
        S = reshape(S,R,C);

        fx      = ([1:(Nx/2+1)]-1)/(2*pi*z_i);

        %Make sure to change the height when computing w-spectra
        hfig = figure(20); set(hfig, 'PaperUnits', 'inches', 'PaperPosition', [0 0  width height]);
        for h = 1:1:Nz/2
            Sp  = S(:,h:Nz:C);
            Sp  = mean(Sp(:,floor(Tstart/p_count):floor(Tend/p_count)),2);
            Sp  = Sp/((h-0.5)/(Nz-1));
            f   = 2*pi*fx*(h-0.5)*L_z/(Nz-1);
            loglog(f(1:Nx/(2*fgr)),Sp(1:Nx/(2*fgr)),'-k'); hold on; %Explicit Filtering
        end
        set(gca,'FontSize',10,'XLim',[1e-3 100],'YLim',[1e-8 1e-2]);
        hold on;X1=0.3;X2=3;H=line([X1 X2],2e-5*[X1^(-5/3) X2^(-5/3)]);set(H,'color',[0 0 0],'LineWidth',2);
        hold off;
        xlabel('k_1z','FontSize',10, 'interpreter', 'tex');
        ylabel('E_T(k)/z','FontSize',10, 'interpreter', 'tex');
        SPR = sprintf('%s%s',SS,'GoMRI_specT'); print(Sprint,SPR);
    end
end
disp(['Plotting done ...'])
