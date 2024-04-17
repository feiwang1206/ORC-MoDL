clear 
pathname=pwd;
addpath([pathname '/utils']);
addpath(genpath([pathname '/utils/mreg_recon_tool']));
addpath([pathname '/utils/image_reconstruction_toolbox']);
setup;
addpath(genpath([pathname '/utils/SphericalHarmonics']));
addpath([pathname '/utils/matmri-master']);
setPath;
    
debug=0;
if gpuDeviceCount>0
    gpu_acc=1;
end
scale=[1 1 1];

institu='_fbu';
folder=[pathname '/data/simu/'];
if ~exist(folder,'dir')
    mkdir(folder);
    mkdir([folder 'smaps/']);
end

subject{1}='subject04';

dim=[64 64 50];
trajectory=load([pathname '/data/trajectory.mat']);
traj = trajectory.trajectory;
traj_idx = trajectory.idx;
traj{1} = traj{1}(traj_idx{1},:);   
traj{2} = traj{2}(traj_idx{2},:); 
dis=traj{1}(:,1).^2+traj{1}(:,2).^2+traj{1}(:,3).^2;
[center]=find(dis==min(dis));
te=[0.0016 0.0016];
dt=5e-6;

traj1{1} = [traj{1}(:,2), -traj{1}(:,1) , traj{1}(:,3)];
traj1{2} = [traj{2}(:,2), -traj{2}(:,1) , traj{2}(:,3)];

xg=repmat((-dim(1)/2:(dim(1)/2-1))',[1 dim(2) dim(3)]);
yg=repmat((-dim(2)/2:(dim(2)/2-1)),[dim(1) 1 dim(3)]);
zg=repmat(permute((-dim(3)/2:(dim(3)/2-1))',[2 3 1]),[dim(1) dim(2) 1]);
if gpu_acc
    xg=gpuArray(xg);yg=gpuArray(yg);zg=gpuArray(zg);
end
%% simulate smaps
if ~exist([pathname '/data/simu/smaps/'],'dir')
    mkdir([pathname '/data/simu/smaps/']);
end
coil_num=16;
n=1;
for cc1=5:15
    for cc2=4:8
        fname=[pathname '/data/simu/smaps/' num2str(n) '.mat'];
        if ~exist(fname,'file')
            smaps=simRcvrSens(dim,cc1,cc2);
            smaps2d=reshape(smaps,[prod(data.dim) size(smaps,4)]);
            [Vt,~]=eig(smaps2d'*smaps2d);
            Vt=Vt(:,end:-1:(end-coil_num+1));
            Ut=smaps2d*Vt;    
            smaps=reshape(Ut,[size(data.wmap) coil_num]);
            save(fname,'smaps');
        end
        n=n+1;
    end
end

traj2=traj;
traj2{1}(:,1)=traj{1}(:,1)*scale(1);
traj2{1}(:,2)=traj{1}(:,2)*scale(2);
traj2{1}(:,3)=traj{1}(:,3)*scale(3);
subset1=(1:length(traj2{1}));
Tt{1}=(subset1+traj_idx{1}(1))*dt+te(1);


traj2{2}(:,1)=traj{2}(:,1)*scale(1);
traj2{2}(:,2)=traj{2}(:,2)*scale(2);
traj2{2}(:,3)=traj{2}(:,3)*scale(3);
subset2=(1:length(traj2{2}));
Tt{2}=(subset2+traj_idx{2}(1))*dt+te(2);
Fg0{1}=orc_segm_nuFTOperator_structure2(traj2(1),dim./scale,dt,10,Tt(1),0.01);
Fg0{2}=orc_segm_nuFTOperator_structure2(traj2(2),dim./scale,dt,10,Tt(2),0.01);
patch=1;
for subj=1:length(subject)
    image=load([pathname '/data/phantom/' subject{subj} '_image.mat']);
    bck=load([pathname '/data/phantom/' subject{subj} '_bck.mat']);
    bck=bck.bck+128;
    bck=255-bck;
    bck=imresize3D(bck,dim);

    for m=1:50
    tic
        if m<10
            disp([subject{subj} '_0' num2str(m)])
            fname=[folder subject{subj} '_0' num2str(m) '.mat'];
        else
            disp([subject{subj} '_' num2str(m)])
            fname=[folder subject{subj} '_' num2str(m) '.mat'];
        end
        if ~exist(fname,'file')
            recon=[];
            lim_l=0.8;
            lim_r=0.4;
            Tvalue=[0.9 0.8 1 0.8;
                1.3 0.9 4.1 1.2;
                0.06 0.05 0.7 0.025];
            cRho=zeros(1,4);
            cT1=zeros(1,4);
            cT2Star=zeros(1,4);
            for ii=1:4
                cRho(ii)=Tvalue(1,ii)*(rand*lim_r+lim_l);
                cT1(ii)=Tvalue(2,ii)*(rand*lim_r+lim_l);
                cT2Star(ii)=Tvalue(3,ii)*(rand*lim_r+lim_l);
            end
            T1=zeros(dim);
            T2Star=zeros(dim);
            rho=zeros(dim);
            brain=zeros(dim);
            for k=1:dim(3)
                for j=1:dim(2)
                    for i=1:dim(1)
                        tem=col(image.image(i,j,k,1:4));
                        tem=tem'/sum(tem);
                        tem2=tem;
                        tem2(isnan(tem))=0;
                        T1(i,j,k)=sum(cT1.*tem2);
                        T2Star(i,j,k)=sum(cT2Star.*tem2);
                        rho(i,j,k)=sum(cRho.*tem2);
                        tem=col(image.image(i,j,k,1:3));
                        tem=tem'/sum(tem);
                        tem2=tem;
                        tem2(isnan(tem))=0;
                        brain(i,j,k)=sum(tem2);
                    end
                end
            end

            TR=0.1;
            alphaa=25;
            T1a=ones(size(T1));
            for m1=1:(m*2)
                T1a=1-(1-T1a.*cos(alphaa/180*pi)).*(exp(-TR./T1));
            end
            % T1a=(1-exp(-TR./T1))./(1-cos(alpha/180*pi)*exp(-TR./T1));
            T1a((T1==0))=0;


            %% shift, rotation, deformation
            if ~debug
                N=64;
                dim0=[N N N];
                rho1=imresize3D(rho,dim0);
                T1a1=imresize3D(T1a,dim0);
                T11=imresize3D(T1,dim0);
                T2Star1=imresize3D(T2Star,dim0);
                bck1=imresize3D(bck,dim0);
                %% scaling
                tem=rand(2000,1);
                N1=floor((1-0.2*tem(m,1))*dim0);
                rho1s=imresize3D(rho1,N1);
                T1a1s=imresize3D(T1a1,N1);
                T2Star1s=imresize3D(T2Star1,N1);
                bck1s=imresize3D(bck1,N1);

                rho1=zeros(dim0);
                rho1(floor((N-N1)/2+1):floor((N-N1)/2+N1),floor((N-N1)/2+1):floor((N-N1)/2+N1),floor((N-N1)/2+1):floor((N-N1)/2+N1))=rho1s;
                T1a1=zeros(dim0);
                T1a1(floor((N-N1)/2+1):floor((N-N1)/2+N1),floor((N-N1)/2+1):floor((N-N1)/2+N1),floor((N-N1)/2+1):floor((N-N1)/2+N1))=T1a1s;
                T2Star1=zeros(dim0);
                T2Star1(floor((N-N1)/2+1):floor((N-N1)/2+N1),floor((N-N1)/2+1):floor((N-N1)/2+N1),floor((N-N1)/2+1):floor((N-N1)/2+N1))=T2Star1s;
                bck1=zeros(dim0);
                bck1(floor((N-N1)/2+1):floor((N-N1)/2+N1),floor((N-N1)/2+1):floor((N-N1)/2+N1),floor((N-N1)/2+1):floor((N-N1)/2+N1))=bck1s;

                %% rotation
                tem=(rand(3000,1)-0.5)*2;            
                rho1=imrotate(rho1,180*tem(m,1),'bilinear','crop');
                T1a1=imrotate(T1a1,180*tem(m,1),'bilinear','crop');
                T2Star1=imrotate(T2Star1,180*tem(m,1),'bilinear','crop');
                bck1=imrotate(bck1,180*tem(m,1),'bilinear','crop');

                rho1=permute(rho1,[3 2 1]);
                T1a1=permute(T1a1,[3 2 1]);
                T2Star1=permute(T2Star1,[3 2 1]);
                bck1=permute(bck1,[3 2 1]);

                tem=(rand(3000,3)-0.5)*2;            
                rho1=imrotate(rho1,180*(tem(m,1)-0.5),'bilinear','crop');
                T1a1=imrotate(T1a1,180*(tem(m,1)-0.5),'bilinear','crop');
                T2Star1=imrotate(T2Star1,180*(tem(m,1)-0.5),'bilinear','crop');
                bck1=imrotate(bck1,180*(tem(m,1)-0.5),'bilinear','crop');

                rho1=permute(rho1,[3 2 1]);
                T1a1=permute(T1a1,[3 2 1]);
                T2Star1=permute(T2Star1,[3 2 1]);
                bck1=permute(bck1,[3 2 1]);

                rho1=permute(rho1,[2 1 3]);
                T1a1=permute(T1a1,[2 1 3]);
                T2Star1=permute(T2Star1,[2 1 3]);
                bck1=permute(bck1,[2 1 3]);

                tem=(rand(3000,3)-0.5)*2;            
                rho1=imrotate(rho1,180*(tem(m,1)-0.5),'bilinear','crop');
                T1a1=imrotate(T1a1,180*(tem(m,1)-0.5),'bilinear','crop');
                T2Star1=imrotate(T2Star1,180*(tem(m,1)-0.5),'bilinear','crop');
                bck1=imrotate(bck1,180*(tem(m,1)-0.5),'bilinear','crop');

                rho1=permute(rho1,[2 1 3]);
                T1a1=permute(T1a1,[2 1 3]);
                T2Star1=permute(T2Star1,[2 1 3]);
                bck1=permute(bck1,[2 1 3]);
                %% deformation
                tem=rand([floor(dim/8),3000,3]);
                dx = -1+2*tem(:,:,:,m,1);
                dy = -1+2*tem(:,:,:,m,2);
                dz = -1+2*tem(:,:,:,m,3);

                tem=rand(3000,3);
                alpha=5;
                fdx=alpha*imresize3D(dx,dim0)*tem(m,1);
                fdy=alpha*imresize3D(dy,dim0)*tem(m,2);
                fdz=alpha*imresize3D(dz,dim0)*tem(m,3);
                [y, x, z]=ndgrid(1:dim0(1),1:dim0(2),1:dim0(3));
                rho1 = griddata(x-fdx,y-fdy,z-fdz,double(rho1),x,y,z);rho1(isnan(rho1))=0;
                T1a1 = griddata(x-fdx,y-fdy,z-fdz,double(T1a1),x,y,z);T1a1(isnan(T1a1))=0;
                T2Star1 = griddata(x-fdx,y-fdy,z-fdz,double(T2Star1),x,y,z);T2Star1(isnan(T2Star1))=0;
                bck1 = griddata(x-fdx,y-fdy,z-fdz,double(bck1),x,y,z);bck1(isnan(bck1))=0;
            else
                rho1=rho;
                T1a1=T1a;
                T2Star1=T2Star;
                bck1=bck;
            end
            rho1=gpuArray(imresize3D(rho1,dim));
            T1a1=gpuArray(imresize3D(T1a1,dim));
            T2Star1=gpuArray(imresize3D(T2Star1,dim));
            bck1=gpuArray(imresize3D(bck1,dim));

    %%                
            mask=zeros(dim);
            bck1s=smooth3(bck1,'box',3);
            mask(bck1s>0.2*max(bck1s(:)))=1;
            
        %%     wmap0    
            dim2=dim*5;
            mask2=zeros(dim2);
            mask2(dim(1)*2+(1:dim(1)),dim(2)*2+(1:dim(2)),dim(3)*2+(1:dim(3)))=mask;
            for ii=1:dim2(1)
                for jj=1:dim2(2)
                    mask2(ii,jj,1:(dim(3)*2))=mask2(ii,jj,2*dim(3)+1);
                end
            end
            px=repmat((-dim2(1)/2:(dim2(1)/2-1))',[1 dim2(2) dim2(3)]);
            py=repmat(-dim2(2)/2:(dim2(2)/2-1),[dim2(1) 1 dim2(3)]);
            pz=repmat(permute((-dim2(3)/2:(dim2(3)/2-1))',[3 2 1]),[dim2(1) dim2(2) 1]);
            pxyz=pz.^2./(px.^2+py.^2+pz.^2);
            B_temp=fftshift(1/3-pxyz).*fftn(((-9e-6)*mask2+0.36e-6*(1-mask2)));
            B_temp(isnan(B_temp))=0;
            wmap0=2*pi*3*42.58*10^6*ifftn(B_temp);
            wmap0=wmap0(dim(1)*2+(1:dim(1)),dim(2)*2+(1:dim(2)),dim(3)*2+(1:dim(3)));

            wmap0=shimming(wmap0,2,0.5,dim);
            wmap0=smooth3((wmap0.*mask),'box',5);
            for n=(m-1)*patch+(1:patch)
                if ~exist(fname,'file')
                    %% phase
                    map=smooth3(randn(dim*3),'box',15);
                    map=smooth3(map,'box',15);
                    map=map(dim(1)+(1:dim(1)),dim(2)+(1:dim(2)),dim(3)+(1:dim(3)));
                    map=map/max(abs(map(:)));
                    phase0=-1i*randn*map*2*pi;

                    %% wmap error
                    map=smooth3(randn(dim*3),'box',15);
                    map=smooth3(map,'box',15);
                    map=map(dim(1)+(1:dim(1)),dim(2)+(1:dim(2)),dim(3)+(1:dim(3)));
                    map=map/max(abs(map(:)));
                    wmap_err1=1.0*smooth3(2*(rand-0.5)*100*map.*mask+2*(rand-0.5)*0.1*wmap0,'box',5);
                    
                    %% wmap shift
                    
        %%
                    wmap_fake=imresize3D(wmap0,dim./scale);
                    wmap_real1=imresize3D(wmap0+wmap_err1,dim);
                    if gpuDeviceCount
                        wmap_fake=gpuArray(wmap_fake);
                        wmap_real1=gpuArray(wmap_real1);
                    end

        %% rawdata
                    smaps_num=floor(rand*54)+1;
                    smaps_l1=load([pathname '/data/simu/smaps/' num2str(smaps_num) '.mat']);
                    smaps_l1=smaps_l1.smaps;
                    image0=rho1.*T1a1.*sin(alphaa/180*pi).*exp(-(Tt{1}(center))./T2Star1+phase0);
                    image0(T2Star1<=0)=0;
                    recon{11}=gather(single(image0));

                    Fg{1}=orc_segm_nuFTOperator_multi_savetime_smaps(Fg0{1},gpuArray(smaps_l1),gpuArray(wmap_real1));
                    Fg{2}=orc_segm_nuFTOperator_multi_savetime_smaps(Fg0{2},gpuArray(smaps_l1),gpuArray(wmap_real1));
                    rawdata1=Fg{1}*image0;
                    rawdata2=Fg{2}*image0;
                    recon{7}{1}=gather(single(rawdata1));
                    recon{7}{2}=gather(single(rawdata2));

                    recon{3}=gather(single(wmap_fake));
                    recon{6}=gather(single(wmap_real1));
                    recon{12}=smaps_num;
                    save(fname,'recon');
                end
            end
        end
    toc                
    end
end


