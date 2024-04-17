clear
pathname=pwd;
pathname=pwd;
addpath([pathname '/utils']);
addpath(genpath([pathname '/utils/mreg_recon_tool']));
addpath([pathname '/utils/image_reconstruction_toolbox']);
setup;
addpath(genpath([pathname '/utils/SphericalHarmonics']));
addpath([pathname '/utils/matmri-master']);
setPath;
%%
pname_s=[pathname '/data/simu/'];
pname_d{1}=[pathname '/data/training_dataset/'];
traj=load([pathname '/data/trajectory.mat']);
traj=traj.traj;
recon_initial(traj,pname_s,pname_d{1});
%% itern field
    net_type='diff';
    n=1;
    maxiter=3;
    lambda=1e-2;
    for iter=1:maxiter
        dirname=[pathname '/data/net/'];

        fname{1}{1}=pname_d{iter};
        fname{1}{2}=pname_d{iter};

        if exist([dirname '/net_block' num2str(iter) '.mat'],'file')
            net_recon{iter}=load([dirname '/net_block' num2str(iter) '.mat'],'net');
        else
            if iter==1
                prepare_input_label(dirname,fname,50,net_type,1e-3);
            else
                prepare_input_label(dirname,fname,20,net_type,1e-4);
            end                
            copyfile([dirname '/net_block.mat'],[dirname '/net_block' num2str(iter) '.mat']);
            copyfile([dirname '/info.mat'],[dirname '/info' num2str(iter) '.mat']);
            net_recon{iter}=load([dirname '/net_block' num2str(iter) '.mat'],'net');
            
        end

        pname_d{iter+1}=[pname_d{1}(1:end-1) num2str(iter) '/'];
        if iter<maxiter
            recon_block(traj,net_recon{iter}.net,pname_s,pname_d{iter},...
                pname_d{iter+1},lambda,5); 
        end
        n=n+1;

    end   

