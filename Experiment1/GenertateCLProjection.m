clear;clc;close all;
%% Initialize
clear;
close all;
%% Geometry
geo=defaultGeometry();
%% Define angles of projection and load phatom image


% shepp_type='yu-ye-wang'; 
% shepp_type='Shepp-Logan'; 
shepp_type='Modified Shepp-Logan';  % (default).
shepp=sheppLogan3D([256,256,256],shepp_type); % Default are 128^3 and Modified shepp-logan
head=headPhantom([256,256,256]);
% angles=[anglesZ1;anglesY;anglesZ2];
numProjs=180;
theta=pi/4;
angles=[linspace(0,2*pi,numProjs);zeros(1,numProjs).*theta;zeros(1,numProjs).*theta];
projections=Ax(shepp,geo,angles);
plotProj(projections,(1:numProjs)*2*pi/180); % angle information not right in the title

imwrite(projections(:,:,45),'pro45.png');
abc=projections(:,:,45);