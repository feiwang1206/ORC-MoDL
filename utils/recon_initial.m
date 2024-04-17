function recon_initial(trajectory,pname_source,pname_desti)
%% single image 256 

if ~exist(pname_desti,'dir')
    mkdir(pname_desti);
end
dt=5e-6;
te=0.0016;
dim0=[64 64 50];
dim=[64 64 48];

if ~exist(pname_desti,'dir')
    mkdir(pname_desti);
end

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
for i=1:55
    winame=[pname_source '/smaps/wi_' num2str(i) '.mat'];
    if ~exist(winame,'file')
        smaps=load([pname_source '/smaps/' num2str(i) '.mat']);
        smaps=imresize4D(smaps.smaps,dim);
        if gpuDeviceCount>0
            wi{1}=orc_segm_nuFTOperator_multi_savetime_smaps_structure(Fg0{1},gpuArray(double(smaps)));
            wi{2}=orc_segm_nuFTOperator_multi_savetime_smaps_structure(Fg0{2},gpuArray(double(smaps)));
        else
            wi{1}=orc_segm_nuFTOperator_multi_savetime_smaps_structure(Fg0{1},(double(smaps)));
            wi{2}=orc_segm_nuFTOperator_multi_savetime_smaps_structure(Fg0{2},(double(smaps)));
        end
        save(winame,"wi");
    end
end


if ~exist([pname_desti '/inputsTra'],'dir')
    mkdir([pname_desti '/inputsTra']);
    mkdir([pname_desti '/labelsTra']);
    mkdir([pname_desti '/inputsVal']);
    mkdir([pname_desti '/labelsVal']);
    mkdir([pname_desti '/norm']);
end
fnames=dir(pname_source);
for ff=4:length(fnames)
    if ff<54
        fname_input=[pname_desti '/inputsVal/' fnames(ff).name];
        fname_label=[pname_desti '/labelsVal/' fnames(ff).name];
    else
        fname_input=[pname_desti '/inputsTra/' fnames(ff).name];
        fname_label=[pname_desti '/labelsTra/' fnames(ff).name];
    end
    fname_norm=[pname_desti '/norm/' fnames(ff).name];
    
    if ~exist(fname_input,'file')
        recon=load([fnames(ff).folder '/' fnames(ff).name]);
        recon=recon.recon;
        field=imresize3D(recon{3},dim);
        smaps=load([pname_source '/smaps/' num2str(recon{12}) '.mat']);
        smaps=imresize4D(smaps.smaps,dim);
        wi=load([pname_source '/smaps/wi_' num2str(recon{12}) '.mat']);
        wi=wi.wi;
        if gpuDeviceCount>0
            Fg{1}=orc_segm_nuFTOperator_multi_savetime_smaps_wi(Fg0{1},gpuArray(double(smaps)),wi{1},gpuArray(double(field)));
            Fg{2}=orc_segm_nuFTOperator_multi_savetime_smaps_wi(Fg0{2},gpuArray(double(smaps)),wi{2},gpuArray(double(field)));
        else
            Fg{1}=orc_segm_nuFTOperator_multi_savetime_smaps_wi(Fg0{1},(double(smaps)),wi{1},(double(field)));
            Fg{2}=orc_segm_nuFTOperator_multi_savetime_smaps_wi(Fg0{2},(double(smaps)),wi{2},(double(field)));
        end
        P{2}=1e-2*max(abs(col(recon{7}{1})));
        recon{1}=single(gather(regularizedReconstruction(Fg{1},(double(recon{7}{1})),P{:},'maxit',10,'verbose_flag',0)));

        P{2}=1e-2*max(abs(col(recon{7}{2})));
        recon{2}=single(gather(regularizedReconstruction(Fg{2},(double(recon{7}{2})),P{:},'maxit',10,'verbose_flag',0)));
        
        norm=max(col(abs(recon{1}+recon{2})/2));

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

        save(fname_norm,'norm');
        if mod(ff-2,10)==0
            printf(fname_input);
        end
    
    end
    
end