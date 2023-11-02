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
    % VARIABLE                                   DESCRIPTION                    UNITS
    %-------------------------------------------------------------------------------------
    % Distances
    geo.DSD = FDDCT;                             % Distance Source Detector      (mm)
    geo.DSO = FODCT;                             % Distance Source Origin        (mm)
    % Detector parameters
    geo.nDetector=[RegionPixelsX;RegionPixelsY];					% number of pixels              (px)
    geo.sDetector=[RegionPixelsX*DetectorPixelSizeX;RegionPixelsY*DetectorPixelSizeY];            % total size of the detector    (mm)
    geo.dDetector=geo.sDetector./geo.nDetector;					% size of each pixel            (mm)

    % Image parameters
    geo.nVoxel=[VoxelsX;VoxelsY;VoxelsZ];                          % number of voxels              (vx)
    geo.sVoxel=[VoxelsX*VoxelSizeX;VoxelsY*VoxelSizeY;VoxelsZ*VoxelSizeZ*4];                      % total size of the image       (mm)
    geo.dVoxel=geo.sVoxel./geo.nVoxel;          % size of each voxel            (mm)
    % Offsets
    geo.offOrigin =[0;0;0];                     % Offset of image from origin   (mm)
    geo.offDetector=[0;DetectorOffsetXCT];                     % Offset of Detector            (mm)

    
    angles1=angles(1:7:2513)./180*pi;
projections=zeros(RegionPixelsY,RegionPixelsX,359,'uint16');
for i=1:359
    projections(:,:,i)=imread('myMultipageFile-359.tif',i);
end
projections1=single(projections-min(projections,[],'all'))./single(max(projections,[],'all')-min(projections,[],'all'));
clear projections angles;

targetGpuName = 'NVIDIA GeForce RTX 3080';
gpuids = GpuIds(targetGpuName);
imgASDPOCSxy1=FDK(projections1,geo,angles1,'gpuids',gpuids);    
plotImg(imgASDPOCSxy1(:,:,1000:1400),'Dim','Z','step',10)
maxasd=max(imgASDPOCSxy1,[],'all');
minasd=min(imgASDPOCSxy1,[],'all');
imwrite(uint16((imgASDPOCSxy1(:,:,1)-minasd)./(maxasd-minasd)*65535),'imgASDPOCSxy1.tif')
for i=500:VoxelsZ
    imwrite(uint16((imgASDPOCSxy1(:,:,i)-minasd)./(maxasd-minasd)*65535),'imgASDPOCSxy1.tif','WriteMode','append')
end

