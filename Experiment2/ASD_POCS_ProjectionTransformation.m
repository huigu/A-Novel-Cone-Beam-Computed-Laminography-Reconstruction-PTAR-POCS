%--------------------------------------------------------------------------
%% Initialize
clear;clc;
close all;
load('ASD_POCS_ProjectionTransformation20230906.mat')

theta=-60.0/180*pi;
alpha=pi/2-abs(theta);
FDD=SrcToDetector;

FDDCT=FDD*cos(alpha);
DetectorOffsetXCT=FDD*sin(alpha);
FODCT=SrcToObject*cos(alpha);

%% Define Geometry
geo=defaultGeometry('mode','cone');                     
    %% Example
    % VARIABLE                                   DESCRIPTION                    UNITS
    %-------------------------------------------------------------------------------------
    % Distances
    geo.DSD = FDDCT;                             % Distance Source Detector      (mm)
    geo.DSO = FODCT;                             % Distance Source Origin        (mm)
    % Detector parameters
    geo.nDetector=[436;333];					% number of pixels              (px)
    geo.sDetector=[436*DetectorPixelSizeX*4;333*DetectorPixelSizeY*4];            % total size of the detector    (mm)
    geo.dDetector=geo.sDetector./geo.nDetector;					% size of each pixel            (mm)

    % Image parameters
    geo.nVoxel=[VoxelsX/4;VoxelsY/4;VoxelsZ/4];                          % number of voxels              (vx)
    geo.sVoxel=[VoxelsX*VoxelSizeX;VoxelsY*VoxelSizeY;VoxelsZ*VoxelSizeZ*4];                      % total size of the image       (mm)
    geo.dVoxel=geo.sVoxel./geo.nVoxel;          % size of each voxel            (mm)
    % Offsets
    geo.offOrigin =[0;0;0];                     % Offset of image from origin   (mm)
    geo.offDetector=[0;DetectorOffsetXCT];                     % Offset of Detector            (mm)
angles1=angles(1:7:2513)./180*pi;
projections=zeros(333,436,359,'uint16');
for i=1:359
    projections(:,:,i)=imread('myMultipageFile-359-0.25.tif',i);
end
projections1=single(projections-min(projections,[],'all'))./single(max(projections,[],'all')-min(projections,[],'all'));


%% Load data and generate projections 
% see previous demo for explanation
% angles=linspace(0,2*pi,100);
% head=headPhantom(geo.nVoxel);
% angles1=angles(1,1:7:2513);
% projections1=single(projections(:,:,1:7:2513));
% 
% projections=Ax(head,geo,angles,'interpolated');
% noise_projections=addCTnoise(projections);
% noise_projections=projections;
% projections=zeros(1334,1747,2513,'uint16');
% for i=1:2513
%     projections(:,:,i)=NewProjections{i,1};
% end

%% Lets create a OS-SART test for comparison
% [imgOSSART,errL2OSSART]=OS_SART(noise_projections,geo,angles,60);

%% Total Variation algorithms
%
%  ASD-POCS: Adaptative Steeppest Descent-Projection On Convex Subsets
% Often called POCS-TV
%==========================================================================
%==========================================================================
%  ASD-POCS minimizes At-B and the TV norm separately in each iteration,
%  i.e. reconstructs the image first and then reduces the TV norm, every
%  iteration. As the other algorithms the mandatory inputs are projections,
%  geometry, angles and maximum iterations.
%
% ASD-POCS has a veriety of optional arguments, and some of them are crucial
% to determine the behaviour of the algorithm. The advantage of ASD-POCS is
% the power to create good images from bad data, but it needs a lot of
% tunning. 
%
%  Optional parameters that are very relevant:
% ----------------------------------------------
%    'maxL2err'    Maximum L2 error to accept an image as valid. This
%                  parameter is crucial for the algorithm, determines at
%                  what point an image should not be updated further.
%                  Default is 20% of the FDK L2 norm.
%                  
% its called epsilon in the paper
epsilon=im3Dnorm(Ax(FDK(projections1,geo,angles1),geo,angles1)-projections1,'L2')*0.15;
%   'alpha':       Defines the TV hyperparameter. default is 0.002. 
%                  However the paper mentions 0.2 as good choice
alpha=0.002;

%   'TViter':      Defines the amount of TV iterations performed per SART
%                  iteration. Default is 20

ng=25;

% Other optional parameters
% ----------------------------------------------
%   'lambda':      Sets the value of the hyperparameter for the SART iterations. 
%                  Default is 1
%
%   'lambdared':   Reduction of lambda Every iteration
%                  lambda=lambdared*lambda. Default is 0.99
%
lambda=1;
lambdared=0.9999;


%   'alpha_red':   Defines the reduction rate of the TV hyperparameter
alpha_red=0.95;

%   'Ratio':       The maximum allowed image/TV update ration. If the TV 
%                  update changes the image more than this, the parameter
%                  will be reduced.default is 0.95
ratio=0.94;

%   'Verbose'      1 or 0. Default is 1. Gives information about the
%                  progress of the algorithm.
qualmeas={'RMSE'};
verb=true;
targetGpuName = 'NVIDIA GeForce RTX 3080';
gpuids = GpuIds(targetGpuName);
% tic
% [imgASDPOCSxy1,qualMeasOut]=ASD_POCS(projections1,geo,angles1,50,...
%                     'TViter',ng,'maxL2err',epsilon,'alpha',alpha,... % these are very important
%                      'lambda',lambda,'lambda_red',lambdared,'Ratio',ratio,'Verbose',verb,'gpuids', gpuids,'QualMeas',qualmeas); % less important.
%                  toc
                 
imgASDPOCSxy1=FDK(projections1,geo,angles1);              
maxasd=max(imgASDPOCSxy1,[],'all');
minasd=min(imgASDPOCSxy1,[],'all');
imwrite(uint16((imgASDPOCSxy1(:,:,1)-minasd)./(maxasd-minasd)*65535),'imgASDPOCSxy1.tif')
for i=2:399
    imwrite(uint16((imgASDPOCSxy1(:,:,i)-minasd)./(maxasd-minasd)*65535),'imgASDPOCSxy1.tif','WriteMode','append')
end

imwrite(uint16(maxasd-(imgASDPOCSxy1(:,:,1))./(maxasd-minasd)*65535),'imgASDPOCSxy1.tif')
for i=2:399
    imwrite(uint16(maxasd-(imgASDPOCSxy1(:,:,i))./(maxasd-minasd)*65535),'imgASDPOCSxy1.tif','WriteMode','append')
end

qualMeasossart=Measure_Quality(head,imgOSSART,{'RMSE','MSSIM','CC','UQI','error_norm'})
qualMeasasdpocs=Measure_Quality(head,imgASDPOCS,{'RMSE','MSSIM','CC','UQI','error_norm'})
qualMeasasdpocsxy1=Measure_Quality(head,imgASDPOCSxy1,{'RMSE','MSSIM','CC','UQI','error_norm'})     
plotImg(abs([head, imgOSSART, imgASDPOCS, imgASDPOCSxy1]) ,'Dim','Z','Slice',10,'clims',[0 0.1])