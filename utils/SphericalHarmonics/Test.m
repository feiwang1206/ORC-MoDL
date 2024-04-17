%% This files provides an example for the use of spherical harmonics

clear all
% close all

addpath(genpath(fullfile('.')))
load('data.mat');
%% we set first at which maximum degree and order we want to make the model
% and run the approximation
degreeMax = 1;
orderMax = 2*degreeMax; % We calculate the same number of order as degree.You do not have to
rhoReference = 0.7; % set the reference radius to 1 meter

% the point used for the description of the sphere are created. Note that
% they are specificaly choosen in order to make the numerical intergal
% more precise
rk = createTargetPointGaussLegendreAndRectangle7(rhoReference,degreeMax,orderMax);

% figure
% plot3(rk(:,1),rk(:,2),rk(:,3),'*');
% axis equal

%% the bc and bs coefficietn are inizialized for each field direction
bc(1).coefficient = zeros(degreeMax+1,orderMax+1);
bc(2).coefficient = zeros(degreeMax+1,orderMax+1);
bc(3).coefficient = zeros(degreeMax+1,orderMax+1);
bs(1).coefficient = zeros(degreeMax+1,orderMax+1);
bs(2).coefficient = zeros(degreeMax+1,orderMax+1);
bs(3).coefficient = zeros(degreeMax+1,orderMax+1);

% we decide to try with a field topology as:
% this is arbitrary
% bc(1).coefficient(1,1) = 2;
% bc(1).coefficient(2,1) = 1;
% bs(1).coefficient(1,1) = 0;
% bs(1).coefficient(2,1) = 3;
% bc(1).coefficient(3,3) = 1;

%% we calculate the value of the modeled field on the point of a sphere
% B  = RebuildField7bis(bc,bs,rhoReference,rk,'sch');
% displayFieldOnSphere(B,rk,'Field calculated from our topology')
dim=size(data.wmap);
for i=1:length(rk)
    po=floor((rk(i,1:3)+1).*dim/2)+1;
    B(1,i)=data.wmap(po(1),po(2),po(3));
    B(2,i)=0;
    B(3,i)=0;
end

%% Then we approximate the previously calculated field as a spherical harmonics expansion
% we should get exactly the same number, as we are making back an expansion
% with the correctorder
[bc2,bs2] = getSphericalHarmonicsCoefficientMeasure7(B(1,:,:),B(2,:,:),B(3,:,:),degreeMax,orderMax,rk,'sch');
% disp(sprintf('%g',bc2(1).coefficient(1,1)))
% disp(sprintf('%g',bs2(1).coefficient(2,2)))
% disp(sprintf('%g',bc2(1).coefficient(3,1)))
% disp(sprintf('%g',bc2(2).coefficient(5,1)))
% disp(sprintf('%g \n Done\n',bc2(3).coefficient(7,1)))

% we can then display them
% displaySHC(bc2,bs2,2,1)
%%
% and reconstruct the Bx field ampltiude in a 2D plan
% x = -rhoReference:rhoReference/10:rhoReference;
% y = -rhoReference:rhoReference/10:rhoReference;
% z= -rhoReference:rhoReference/10:rhoReference;;
x=((-dim(1)/2):(dim(1)/2-1))/(dim(1)/2);
y=((-dim(2)/2):(dim(2)/2-1))/(dim(2)/2);
z=((-dim(3)/2):(dim(3)/2-1))/(dim(3)/2);
bc3=bc2;
bc3(1).coefficient(2,[1 2])=0;
B1  = RebuildField7(bc3,bs2,rhoReference,x,y,z,'sch');
Bx1 = squeeze(B1(1,:,:,:));
mask=zeros(dim);
mask(data.anatomical>0.1*max(col(data.anatomical)))=1;
% %%
% figure
% imagesc(array2mosaic((data.wmap-Bx1).*mask),[-500 500]);axis equal;title('x component of the magentic field')
figure
imagesc(array2mosaic(Bx1.*mask),[-500 500]);axis equal;title('x component of the magentic field')