clear 
pathname=pwd;
addpath([pathname '/utils']);
addpath(genpath([pathname '/utils/mreg_recon_tool']));
addpath([pathname '/utils/image_reconstruction_toolbox']);
setup;
addpath(genpath([pathname '/utils/SphericalHarmonics']));
addpath([pathname '/utils/matmri-master']);
setPath;

data=load([pathname '/data/MREG_example_data/data.mat']);
data=data.data;
traj=data.trajectory.trajectory;
traj_idx=data.trajectory.idx;

param.dim0=data.dim;
param.dim=[64 64 48];
param.fov=data.fov;

% assemble arguments for regularizedReconstruction
operator_string='TV';
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

param.P=P;
param.traj=traj;
param.Vt=Vt;
param.GPU=1;
param.dt=5e-6;
param.te=0.0016;
param.tseg=10;
param.traj_idx=traj_idx;
param.smaps=data.smaps;
for i=1:length(traj)
    param.T{i}=param.dt*param.traj_idx{i}+param.te;
end
param.field0=imresize3D(data.wmap,param.dim);
param.traj1{1}=param.traj{1}(param.traj_idx{1},:);
param.traj1{2}=param.traj{2}(param.traj_idx{2},:);
for tt=1:size(param.traj1{1},2)
    param.traj2{1}(:,tt)=param.traj1{1}(:,tt)*(param.dim0(tt)/param.dim(tt));
    param.traj2{2}(:,tt)=param.traj1{2}(:,tt)*(param.dim0(tt)/param.dim(tt));
end

if param.GPU==1
    param.F0{1}=orc_segm_nuFTOperator_structure(param.traj2(1),param.dim,gpuArray(param.smaps),param.dt,param.tseg,param.T(1),0.01);
    param.F0{2}=orc_segm_nuFTOperator_structure(param.traj2(2),param.dim,gpuArray(param.smaps),param.dt,param.tseg,param.T(2),0.01);
    param.F0{4}=orc_segm_nuFTOperator_structure(param.traj2,param.dim,gpuArray(param.smaps),param.dt,param.tseg,param.T,0.01);
    
    param.F{1}=orc_segm_nuFTOperator_multi_savetime(param.F0{1},gpuArray(param.field0));
    param.F{2}=orc_segm_nuFTOperator_multi_savetime(param.F0{2},gpuArray(param.field0));
    param.F{4}=orc_segm_nuFTOperator_multi_savetime(param.F0{4},gpuArray(param.field0));
else
    param.F0{1}=orc_segm_nuFTOperator_structure(param.traj2(1),param.dim,param.smaps,param.dt,param.tseg,param.T(1),0.01);
    param.F0{2}=orc_segm_nuFTOperator_structure(param.traj2(2),param.dim,param.smaps,param.dt,param.tseg,param.T(2),0.01);
    param.F0{4}=orc_segm_nuFTOperator_structure(param.traj2,param.dim,param.smaps,param.dt,param.tseg,param.T,0.01);
    
    param.F{1}=orc_segm_nuFTOperator_multi_savetime(param.F0{1},param.field0);
    param.F{2}=orc_segm_nuFTOperator_multi_savetime(param.F0{2},param.field0);
    param.F{4}=orc_segm_nuFTOperator_multi_savetime(param.F0{4},param.field0);
end
rawdata1=load([pathname '/data/MREG_example_data/kspace_data_shot1.mat']);
rawdata2=load([pathname '/data/MREG_example_data/kspace_data_shot2.mat']);
rawdata{1}=double(rawdata1.rawdata);
rawdata{2}=double(rawdata2.rawdata);
%% recon MBIR
filename=[pathname '/data/MREG_example_data/recon_MBIR.mat'];
if ~exist(filename,'file')
    param.shot=1;
    recon_MBIR = reconstruction_MBIR(rawdata,param);
    save(filename,'recon_MBIR');
end
%% reconcor
filename=[pathname '/data/MREG_example_data/recon_MoDL.mat'];
if ~exist(filename,'file')
    net=cell(3);
    for n=1:3
        net{n}=load([pathname '/data/net/net_block' num2str(n) '.mat'],'net');
    end
    recon_MoDL=reconstruction_MoDL(rawdata,net,param);
    save(filename,'recon_MoDL');
end
%% display
figure,imagesc(array2mosaic(permute(abs(recon_MBIR{1}),[1 2 3])));axis equal;axis off;colormap gray;colorbar;title('Original MBIR')
figure,imagesc(array2mosaic(permute(abs(recon_MoDL{4}),[1 2 3])));axis equal;axis off;colormap gray;colorbar;title('Proposed ORC-MoDL')
figure,imagesc(array2mosaic(permute(recon_MoDL{3}-param.field0,[1 2 3])));axis equal;axis off;colormap parula;colorbar;title('Estimated field deviation')
