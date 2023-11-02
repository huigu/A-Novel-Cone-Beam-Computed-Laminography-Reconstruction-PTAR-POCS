% 该脚本主要用来将CL投影转换为CT投影，参考文献《A reconstruction method for cone-beam　computed laminography based on　projection transformation》
clear;clc;close all
load('angles.mat')
fileList = dir('../Lego_Lamino30deg_XTH/*.tif');
VoxelsX=1596;
VoxelsY=1148;
VoxelsZ=1596;
VoxelSizeX=0.102950776887476;
VoxelSizeY=0.102950776887476;
VoxelSizeZ=0.102950776887476;
OffsetX=0;
OffsetY=0;
OffsetZ=0;
SrcToObject=392.072625160217;
SrcToDetector=967.320983886719;
RegionPixelsX=1596;
RegionPixelsY=1148;
DetectorPixelSizeX=0.254;
DetectorPixelSizeY=0.254;
DetectorOffsetX=0;
DetectorOffsetY=0;
ProjectionNum=2513;

theta=-60.0/180*pi;
alpha=pi/2-abs(theta);
FDD=SrcToDetector;

FDDCT=FDD*cos(alpha);
DetectorOffsetXCT=FDD*sin(alpha);

width=RegionPixelsX;
height=RegionPixelsY;
Psize=DetectorPixelSizeX;
Fi=width/2;
Fj=height/2;
NewProjections=zeros(ProjectionNum*1747*1334,'uint16');
NewP=zeros(ProjectionNum,3);
for iter=1:(ProjectionNum*1747*1334)
    [num,i,j]=ind2sub([ProjectionNum,1747,1334],iter);
    %     将CL投影坐标四角转为空间坐标
    XA=TransformToSpatialPosition(1,1,Fi,Fj,FDD,Psize);
    XB=TransformToSpatialPosition(width,1,Fi,Fj,FDD,Psize);
    XC=TransformToSpatialPosition(width,height,Fi,Fj,FDD,Psize);
    XD=TransformToSpatialPosition(1,height,Fi,Fj,FDD,Psize);
    %     将四角空间坐标转为CT投影坐标
    XA1=TransformToCTProjection(XA,FDD,theta);
    XB1=TransformToCTProjection(XB,FDD,theta);
    XC1=TransformToCTProjection(XC,FDD,theta);
    XD1=TransformToCTProjection(XD,FDD,theta);
    %     获取CT投影坐标范围
    xP=min(XA1(1),XD1(1));
    xQ=max(XB1(1),XC1(1));
    xt=linspace(xP,xQ,abs((xP-xQ)./Psize));
    yP=XA1(2);
    yS=XD1(2);
    yt=linspace(yP,yS,abs((yS-yP)./Psize./cos(theta)));
    zt=(yt.*sin(theta)-FDD.*sin(theta))./cos(theta);
    %     NewP(num,:)=[1,0,0;0,cos(alpha),sin(alpha);0,-sin(alpha),cos(alpha)]*[xt(1);yt(1);zt(1)];
    img=imread([fileList(num).folder '\' fileList(num).name]);
    
    %     for iter=1:(1747*1334)%size(xt,2)
    %         [i,j]=ind2sub([1747,1334],iter);
    %             获取CT投影坐标在CL投影上强度一致的坐标
    index=[FDD*xt(i)/yt(j),FDD,FDD*zt(j)/yt(j)];
    %             在CL投影上对应的像素坐标
    indexInt=[index(1)./Psize+Fi,Fj-index(3)./Psize];
    NewProjections(iter)=InstensityByInterpolation(indexInt,img);
    %     end
    
    
    
    
    %     parfor i=1:1747%size(xt,2)
    %         for j=1:1334%size(zt,2)
    %             %             获取CT投影坐标在CL投影上强度一致的坐标
    %             index=[FDD*xt(i)/yt(j),FDD,FDD*zt(j)/yt(j)];
    %             %             在CL投影上对应的像素坐标
    %             indexInt=[index(1)./Psize+Fi,Fj-index(3)./Psize];
    %             NewProjections(i,j,num)=InstensityByInterpolation(indexInt,img);
    %         end
    %     end
    %     subplot(1,2,1)
    %     imagesc(NewProjections(:,:,num))
    %    subplot(1,2,2)
    %     imagesc(img)
end

A=[1,1,0];B=[width,1,0];C=[width,height,0];D=[1,height,0];
abc=[XA1;XB1;XC1];
acd=[XA1;XC1;XD1];
hold on
% fill3(XA,XB,XC,'r')
% fill3(XA,XB,XD,'b')
fill3(abc(:,1),abc(:,2),abc(:,3),'r')
fill3(acd(1:3,1),acd(1:3,2),acd(1:3,3),'b')

[x,y]=meshgrid(xt,yt);
z=(y.*sin(theta)-FDD.*sin(theta))./cos(theta);
pcshow([x(:),y(:),z(:)]);

% img=imread([fileList(num).folder '\' fileList(num).name]);


function instensity=InstensityByInterpolation(indexInt,img)
indexFloor=floor(indexInt);
indexFloor(indexFloor<0)=0;
if indexFloor(1)>size(img,1)
    indexFloor(1)=size(img,1);
end
if indexFloor(2)>size(img,2)
    indexFloor(2)=size(img,2);
end

instensity=0;
if indexFloor(1)>0 && indexFloor(1)<=size(img,1) && indexFloor(2)>0 && indexFloor(2)<=size(img,2)
    instensity=instensity+prod(abs(indexInt-indexFloor-1)).*img(indexFloor(1),indexFloor(2));
end
if indexFloor(1)>0 && indexFloor(1)<=size(img,1) && indexFloor(2)>=0 && indexFloor(2)<size(img,2)
    instensity=instensity+prod(abs(indexInt-(indexFloor+[0,1])-1)).*img(indexFloor(1),indexFloor(2)+1);
end
if indexFloor(1)>=0 && indexFloor(1)<size(img,1) && indexFloor(2)>0 && indexFloor(2)<=size(img,2)
    instensity=instensity+prod(abs(indexInt-(indexFloor+[1,0])-1)).*img(indexFloor(1)+1,indexFloor(2));
end
if indexFloor(1)>=0 && indexFloor(1)<size(img,1) && indexFloor(2)>=0 && indexFloor(2)<size(img,2)
    instensity=instensity+prod(abs(indexInt-(indexFloor+[1,1])-1)).*img(indexFloor(1)+1,indexFloor(2)+1);
end
%     instensity=abs(indexInt-indexFloor-1).*img(indexFloor)+...
%         abs(indexInt-(indexFloor+[1,1])-1).*img(indexFloor+[1,1])+...
%         abs(indexInt-(indexFloor+[1,0])-1).*img(indexFloor+[1,0])+...
%         abs(indexInt-(indexFloor+[0,1])-1).*img(indexFloor+[0,1]);
return;
end

function Xw1=TransformToCTProjection(Xw,FDD,theta)
Xw1=FDD.*sin(theta)./(Xw(2).*sin(theta)-Xw(3).*cos(theta)).*Xw;
end
function Xw=TransformToSpatialPosition(i,j,Fi,Fj,FDD,Psize)
Xw=[(i-Fi).*Psize,FDD,(Fj-j).*Psize];
end
