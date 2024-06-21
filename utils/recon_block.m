function recon_block(trajectory,net,pname_source0,pname_source,pname_desti,lambda0,maxiter)
if ~exist(pname_desti,'dir')
    mkdir(pname_desti);
end
dim0=[64 64 50];
dim=[64 64 48];
dt=5e-6;
te=0.0016;

traj = trajectory.trajectory;
traj_idx = trajectory.idx;
for ii=1:3
    traj1{1}(:,ii) = traj{1}(traj_idx{1},ii)*(dim0(ii)/dim(ii));   
    traj1{2}(:,ii) = traj{2}(traj_idx{2},ii)*(dim0(ii)/dim(ii)); 
end
Tt{1}=dt*traj_idx{1}+te;
Tt{2}=dt*traj_idx{2}+te;
Fg0{1}=orc_segm_nuFTOperator_structure2(traj1(1),dim,dt,10,Tt(1),0.01);
Fg0{2}=orc_segm_nuFTOperator_structure2(traj1(2),dim,dt,10,Tt(2),0.01);

operator(1).handle = @finiteDifferenceOperator;
operator(1).args = {1};
operator(2).handle = @finiteDifferenceOperator;
operator(2).args = {2};
operator(3).handle = @finiteDifferenceOperator;
operator(3).args = {3};

counter = 1;
P{counter} = @L1Norm;
counter = counter + 1;
P{counter} = 1e-2;
counter = counter + 1;
for n=1:3
    P{counter} = operator(n).handle(operator(n).args{:});
    counter = counter + 1;
end


if ~exist([pname_desti '/inputsVal'],'dir')
    mkdir([pname_desti '/inputsTra']);
    mkdir([pname_desti '/inputsVal']);
    mkdir([pname_desti '/labelsTra']);
    mkdir([pname_desti '/labelsVal']);
    mkdir([pname_desti '/norm']);
end

for ii=1:2
    if ii==1
        fnames=dir([pname_source '/inputsVal']);
    else
        fnames=dir([pname_source '/inputsTra']);
    end
    for ff=3:length(fnames)
        if ii==1
            fname_input=[pname_desti '/inputsVal/' fnames(ff).name];
            fname_label=[pname_desti '/labelsVal/' fnames(ff).name];
        else
            fname_input=[pname_desti '/inputsTra/' fnames(ff).name];
            fname_label=[pname_desti '/labelsTra/' fnames(ff).name];
        end
        fname_norm_s=[pname_source '/norm/' fnames(ff).name];
        fname_norm_d=[pname_desti '/norm/' fnames(ff).name];
        if ~exist(fname_input,'file')
            image=load([fnames(ff).folder '/' fnames(ff).name]);
            image=image.image;
            recon=load([pname_source0 '/' fnames(ff).name]);
            recon=recon.recon;
            norm0=load(fname_norm_s);
            norm=norm0.norm;
            smaps=load([pname_source0 '/smaps/' num2str(recon{12}) '.mat']);
            smaps=imresize4D(smaps.smaps,dim);
            wi=load([pname_source0 '/smaps/wi_' num2str(recon{12}) '.mat']);
            wi=wi.wi;

            tmp=predict(net,image(:,:,:,1:5));
            recon{1}=(image(:,:,:,1)+1i*image(:,:,:,3))*norm;
            recon{2}=(image(:,:,:,2)+1i*image(:,:,:,4))*norm;
            image_c1=(tmp(:,:,:,2)+1i*tmp(:,:,:,3))*norm+(recon{1}+recon{2})/2;
            image_c2=(tmp(:,:,:,2)+1i*tmp(:,:,:,3))*norm+(recon{1}+recon{2})/2;
            field=(tmp(:,:,:,1)*100+image(:,:,:,5)*1000);
            if gpuDeviceCount>0
                Fg{1}=orc_segm_nuFTOperator_multi_savetime_smaps_wi(Fg0{1},gpuArray(double(smaps)),wi{1},gpuArray(double(field)));
                Fg{2}=orc_segm_nuFTOperator_multi_savetime_smaps_wi(Fg0{2},gpuArray(double(smaps)),wi{2},gpuArray(double(field)));
            else
                Fg{1}=orc_segm_nuFTOperator_multi_savetime_smaps_wi(Fg0{1},(double(smaps)),wi{1},(double(field)));
                Fg{2}=orc_segm_nuFTOperator_multi_savetime_smaps_wi(Fg0{2},(double(smaps)),wi{2},(double(field)));
            end
            P{2}=lambda0*max(abs(col(recon{7}{1})));
            recon{1}=single(gather(regularizedReconstruction(Fg{1},(double(recon{7}{1})),P{:},...
                'maxit',maxiter,'verbose_flag',0,'z0',image_c1)));
    
            P{2}=lambda0*max(abs(col(recon{7}{2})));
            recon{2}=single(gather(regularizedReconstruction(Fg{2},(double(recon{7}{2})),P{:},...
                'maxit',maxiter,'verbose_flag',0,'z0',image_c2)));

            image=zeros([dim 5]);
            image(:,:,:,1)=real(gather(recon{1}))/norm;
            image(:,:,:,2)=real(gather(recon{2}))/norm;
            image(:,:,:,3)=imag(gather(recon{1}))/norm;
            image(:,:,:,4)=imag(gather(recon{2}))/norm;
            image(:,:,:,5)=(gather(field))/1000;
            image=single(image);
            save(fname_input,'image');

            image=zeros([dim 3]);
            image(:,:,:,1)=gather(imresize3D(recon{6},dim)-field)/100;
            image(:,:,:,2)=gather(real(imresize3D(recon{11},dim)*sqrt(prod(dim0)/prod(dim))-(recon{1}+recon{2})/2))/norm;
            image(:,:,:,3)=gather(imag(imresize3D(recon{11},dim)*sqrt(prod(dim0)/prod(dim))-(recon{1}+recon{2})/2))/norm;
            image=single(image);
            save(fname_label,'image');

            save(fname_norm_d,'norm');
            if mod(ff-2,10)==0
                printf(fname_input);
            end


        end
    end
end

