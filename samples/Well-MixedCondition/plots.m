clear,clc

SA='../../';

%add utilities to the path temorarily
addpath ../../utilities

%read inputs from LESinputs.txt
readinputs(SA)

%averaging period
t_avg=round([nsteps*3/4 nsteps]/p_count);

skipframe=1;

z_u=linspace(dz/2,l_z+dz/2,Nz)/l_z;
z_w=linspace(0,l_z,Nz)/l_z;

figure;hold on;

C=zeros(Nz-1,part_model,length(1:skipframe:nframe));
framecnt=0;
for t=1:skipframe:nframe
    framecnt=framecnt+1;
    break_flag=0;
    
    for m=1:part_model
    
        if(m==1)
            filename=sprintf('%s%s%04d%s',SA,'output/part_frame/part_frame_nomodel',t,'.bin');
        elseif(m==2)
             filename=sprintf('%s%s%04d%s',SA,'output/part_frame/part_frame_weiliso',t,'.bin');
        end
        if(~exist(filename,'file'))
            fprintf('%s does not exist...terminating.\n',filename)
            break_flag=1;
            break
        end
        fidpart=fopen(filename,'r','l');
        ipart=fread(fidpart,1,'int')/3/8;
        particle(:,:,m)=reshape(fread(fidpart,ipart*3,'double'),ipart,3);
        particle(:,:,m)=particle(:,:,m)*z_i;
        fclose('all');
    
        z=particle(:,3,m);
        indomain=find( (z>=0) .* (z<(Nz-1)*dz) );
%         fprintf('in domain = %d\n',length(indomain))
        
        k=ceil(z(indomain)/dz);
        
        for ii=1:length(indomain)
            C(k(ii),m,framecnt) = C(k(ii),m,framecnt) + 1;
        end
        
        RMSE(m)=sqrt(mean((C(:,m,framecnt)-npart/(Nz-1)).^2));
        if(t==1)
            RMSE0(m)=RMSE(m);
        end
        
    end
    
    if(break_flag==1)
        break
    end
        
    hold off
    plot(C(:,1,framecnt),z_u(1:end-1),'-ok')
    hold on;
    plot(npart/(Nz-1)*[1 1],[0 z_u(end-1)],'--k')
    if(part_model==1)
        legend('theoretical','nomodel')
    elseif(part_model>1)
        plot(C(:,2,framecnt),z_u(1:end-1),'-ob')
    end
    
    xlim=get(gca,'xlim');
    if(exist('ht1','var') && ishandle(ht1))
        delete(ht1)
    end
    ht1=text(0.8*diff(xlim)+xlim(1),0.9*z_u(end-1),sprintf('RMSE=%3.1f',RMSE(1)/RMSE0(1)));
    if(part_model>1)
        if(exist('ht2','var') && ishandle(ht2))
            delete(ht2)
        end
        ht2=text(0.8*diff(xlim)+xlim(1),0.8*z_u(end-1),sprintf('RMSE=%3.1f',RMSE(2)/RMSE0(2)),'color','b');
    end
    pause(1)
    
end
