clear;close all;clc;
load CLprojection.mat
%% CLFDK重建
theta=-pi/4;
geo=defaultGeometry();
numProjs=90;
angles=[linspace(0,2*pi,numProjs);ones(1,numProjs).*theta;zeros(1,numProjs).*theta];
CLFDK=FDK(projections(:,:,1:2:180),geo,angles);
% clear projections 
%% 投影转换FDK重建
load ProjectionTransformation.mat
% plotImg(CTProjections,'Dim','Z')

alpha=pi/2-abs(theta);
geo1=geo;
geo1.DSD = geo.DSD*cos(alpha);                       % Distance Source Detector      (mm)
geo1.DSO = geo.DSO*cos(alpha);                        % Distance Source Origin        (mm)
% geo1.offDetector=[0;FDD*sin(alpha)]; 

% Detector parameters
geo1.nDetector=[294;367];					% number of pixels              (px)
geo1.sDetector=[294*0.8;367*0.8];            % total size of the detector    (mm)
geo1.dDetector=geo1.sDetector./geo1.nDetector;					% size of each pixel            (mm)
geo1.nVoxel=[256;256;256];                          % number of voxels              (vx)
geo1.sVoxel=[128;128;128];                   % total size of the image       (mm)
geo1.dVoxel=geo1.sVoxel./geo1.nVoxel;          % size of each voxel            (mm)

CTFDK=FDK(CTProjections(:,:,1:2:180),geo1,linspace(0,2*pi,numProjs));
% CTMax=max(CTFDK,[],'all');
% CTMin=min(CTFDK,[],'all');
% plotImg(CTFDK,'Dim','Z')

% figure
% imagesc(CLFDK(:,:,128));
% axis equal
% figure
% imagesc(CTFDK(:,:,128));
% axis equal

%% ASD-POCS投影转换重建
epsilon=im3Dnorm(Ax(FDK(CTProjections,geo1,linspace(0,2*pi,numProjs)),geo1,linspace(0,2*pi,numProjs))-CTProjections,'L2')*0.15;
alpha=0.002;
ng=25;
lambda=1;
lambdared=0.9999;
alpha_red=0.95;
ratio=0.94;
qualmeas={'RMSE'};
verb=true;
targetGpuName = 'NVIDIA GeForce RTX 3080';
gpuids = GpuIds(targetGpuName);

[CTASDPOCS,qualMeasOut]=ASD_POCS(CTProjections,geo1,linspace(0,2*pi,numProjs),50,...
                    'TViter',ng,'maxL2err',epsilon,'alpha',alpha,... % these are very important
                     'lambda',lambda,'lambda_red',lambdared,'Ratio',ratio,'Verbose',verb,'gpuids', gpuids,'QualMeas',qualmeas); % less important.
save CTASDPOCS.mat CTASDPOCS
 plotImg(CTASDPOCS,'Dim','Z')
% qualMeasossart=Measure_Quality(head,imgOSSART,{'RMSE','MSSIM','CC','UQI','error_norm'})
%% PTAR-POCS投影转换重建

[CTPTARPOCS,qualMeasOut]=ASD_POCS(CTProjections,geo1,linspace(0,2*pi,numProjs),50,...
                    'TViter',ng,'maxL2err',epsilon,'alpha',alpha,... % these are very important
                     'lambda',lambda,'lambda_red',lambdared,'Ratio',ratio,'Verbose',verb,'gpuids', gpuids,'QualMeas',qualmeas); % less important.
save CTPTARPOCS.mat CTPTARPOCS
 plotImg(CLFDK,'Dim','Z')


[CLPTARPOCS,qualMeasOut]=ASD_POCS(projections,geo,angles,50,...
                    'TViter',ng,'maxL2err',epsilon,'alpha',alpha,... % these are very important
                     'lambda',lambda,'lambda_red',lambdared,'Ratio',ratio,'Verbose',verb,'gpuids', gpuids,'QualMeas',qualmeas); % less important.
save CLPTARPOCS.mat CLPTARPOCS







shepp_type='Modified Shepp-Logan';  % (default).
shepp=sheppLogan3D([256,256,256],shepp_type); % Default are 128^3 and Modified shepp-logan
CLFDKMeasossart=Measure_Quality(shepp(:,:,128),CLFDK(:,:,128),{'RMSE','MSSIM','CC','UQI','error_norm'});
CTFDKMeasossart=Measure_Quality(shepp(:,:,128),CTFDK(:,:,128),{'RMSE','MSSIM','CC','UQI','error_norm'});
CTASDPOCSMeasossart=Measure_Quality(shepp(:,:,128),CTASDPOCS(:,:,128),{'RMSE','MSSIM','CC','UQI','error_norm'});
CTPTARPOCSMeasossart=Measure_Quality(shepp(:,:,128),CTPTARPOCS(:,:,128),{'RMSE','MSSIM','CC','UQI','error_norm'});
CLPTARPOCSSMeasossart=Measure_Quality(shepp(:,:,128),CLPTARPOCS(:,:,128),{'RMSE','MSSIM','CC','UQI','error_norm'});
Measossart=[CLFDKMeasossart;CTFDKMeasossart;CTASDPOCSMeasossart;CTPTARPOCSMeasossart;CLPTARPOCSSMeasossart];
% 
% [0.043239966,5.7951517e-08,0.97391754,0.97190255,177.11090]
% [0.37232536,1.2802381e-08,0.40903839,0.21237694,1525.0448]
% [0.41519678,1.0290297e-08,0.36071393,0.17049415,1700.6467]
% [0.36632201,1.2440876e-08,0.39973879,0.20624949,1500.4551]
% [0.13420469,4.2803784e-08,0.72072893,0.71410239,549.70239]
% 
% 
CTMax=max(CTFDK,[],'all');
CTMin=min(CTFDK,[],'all');
imwrite(uint16((CTFDK(:,:,1)-CTMin)./(CTMax-CTMin)*65535),'CTFDK.tif')
for i=2:256
    imwrite(uint16((CTFDK(:,:,i)-CTMin)./(CTMax-CTMin)*65535),'CTFDK.tif','WriteMode','append');
end

