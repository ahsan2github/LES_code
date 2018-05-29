%matlab function to read in fortan binary data%
%used as var = loadbin(filename,Nz,'ed')
%where var is the name of the variable being loaded
%and ed is little endian (for altix 'l') or big endian (for IBM 'b')
%example for loading au use: >>au=loadbin('au.bin',Nz,'r')
%to load au from the IBM P3 or P4 systems
%note that this function assumes Nx=Nz!!!

function [X] = loadbin(file,Nz,ed)
cnt=0;
ts=0; mns=0;
if(exist(file,'file')==0); error('%s%s',file,' not found');end
sz=getfield(dir(file),'bytes');
fid = fopen(file,'r',ed);
if(fid == -1); disp('file not found...'); return; end
while ts == 0
dum1 = fread(fid,1,'int');
    if (dum1 > 0) & (dum1 < 1000000)
        ts = isempty(dum1);
    else
        ts = 1;
    end
    if ts == 1
        break
    else
        X1 = fread(fid,dum1/8,'double');
        if (cnt == 0); NN=dum1/8/Nz; end
        if (length(X1) ~= NN*Nz); display('returning');return; end
        
%        if (strcmp(file(1:4),'spec') == 1)
            X1 = reshape(X1,NN,Nz);
%         elseif (strcmp(file(end-5:end),'fs.bin') == 1)
%             X1 = reshape(X1,NN,Nz);
%        end
        if cnt == 0
            zz = size(X1);
            X = zeros(floor(sz/(8+8*Nz*NN))*zz(1),zz(2));
            X(1:zz(1),:) = X1;
        else
            X(cnt*zz(1)+1:(cnt+1)*zz(1),:) = X1;
        end

        dum2 = fread(fid,1,'int');
        while(dum2 ~= dum1)
            dum2 = fread(fid,1,'int');
        end
        cnt=cnt+1;

    end
end

fclose(fid);