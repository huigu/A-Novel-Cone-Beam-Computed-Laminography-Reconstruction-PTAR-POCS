% 该脚本主要用来将CL投影转换为CT投影，参考文献《A reconstruction method for cone-beam　computed laminography based on　projection transformation》
clear;clc;close all
%% 参数设置
load('angles.mat')
% 投影参数
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
RegionPixelsX=1747;%1596
RegionPixelsY=1334;%1148
DetectorPixelSizeX=0.254;
DetectorPixelSizeY=0.254;
DetectorOffsetX=0;
DetectorOffsetY=0;
ProjectionNum=2513;

% 原始投影路径
fileList = dir('../Lego_Lamino30deg_XTH/*.tif');
% imgs=cell(ProjectionNum);
% parfor num=1:ProjectionNum
% imgs{num}=imread([fileList(num).folder '\' fileList(num).name]);
% end


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
NewProjections=cell(ProjectionNum,1);
% NewProjections=zeros(ProjectionNum,1334,1747,'uint16');
%% 投影转换

for num=1:ProjectionNum
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
    img=imread([fileList(num).folder '\' fileList(num).name]);
    NewProjection=zeros(1555,1942,'uint16');%1334*1747
    for iter=1:(1555*1942)
        [j,i]=ind2sub([1555,1942],iter);
%     for i=1:size(xt,2)
%         for j=1:size(zt,2)
            %                         获取CT投影坐标在CL投影上强度一致的坐标
            index=[FDD*xt(i)/yt(j),FDD,FDD*zt(j)/yt(j)];
            %                         在CL投影上对应的像素坐标
            indexInt=[index(1)./Psize+Fi,Fj-index(3)./Psize];
%             indexInt=[index(1),index(3)]./Psize;
            NewProjection(iter)=InstensityByInterpolation(indexInt,img);
%         end
%     end
    end
    NewProjections{num,1}=NewProjection;
    %     subplot(1,2,1)
    %     imagesc(NewProjections{16,1})
    %     subplot(1,2,2)
    %     imagesc(img)
end

%% 保存投影
save('ProjectionTransformation.mat','NewProjections','FDDCT','DetectorOffsetXCT');


imwrite(uint16(projections1(:,:,1)),'myMultipageFile-359.tif')
for i=2:359
    imwrite(uint16(projections1(:,:,i)),'myMultipageFile-359.tif','WriteMode','append');
end
imwrite(projections1,'myMultipageFile-359.tif')


% parfor iter=1:(ProjectionNum*1334*1747)
%     [num,j,i]=ind2sub([ProjectionNum,1334,1747],iter);
%     %     将CL投影坐标四角转为空间坐标
%     XA=TransformToSpatialPosition(1,1,Fi,Fj,FDD,Psize);
%     XB=TransformToSpatialPosition(width,1,Fi,Fj,FDD,Psize);
%     XC=TransformToSpatialPosition(width,height,Fi,Fj,FDD,Psize);
%     XD=TransformToSpatialPosition(1,height,Fi,Fj,FDD,Psize);
%     %     将四角空间坐标转为CT投影坐标
%     XA1=TransformToCTProjection(XA,FDD,theta);
%     XB1=TransformToCTProjection(XB,FDD,theta);
%     XC1=TransformToCTProjection(XC,FDD,theta);
%     XD1=TransformToCTProjection(XD,FDD,theta);
%     %     获取CT投影坐标范围
%     xP=min(XA1(1),XD1(1));
%     xQ=max(XB1(1),XC1(1));
%     xt=linspace(xP,xQ,abs((xP-xQ)./Psize));
%     yP=XA1(2);
%     yS=XD1(2);
%     yt=linspace(yP,yS,abs((yS-yP)./Psize./cos(theta)));
%     zt=(yt.*sin(theta)-FDD.*sin(theta))./cos(theta);
%     
% %     for i=1:size(xt,2)
% %         for j=1:size(zt,2)
%             %             获取CT投影坐标在CL投影上强度一致的坐标
%             index=[FDD*xt(i)/yt(j),FDD,FDD*zt(j)/yt(j)];
%             %             在CL投影上对应的像素坐标
%                         indexInt=[index(1)./Psize+Fi,Fj-index(3)./Psize];
% %             indexInt=[index(1),index(3)]./Psize;
%             NewProjections(iter)=InstensityByInterpolation(indexInt,imgs{num});
% %         end
% %     end
% %     subplot(1,2,1)
% %     imagesc(NewProjections(:,:,num))
% %     subplot(1,2,2)
% %     imagesc(img)
% end







%% 子函数

function instensity=InstensityByInterpolation(indexInt,img)
indexFloor=floor(indexInt);
indexFloor(indexFloor<=0)=1;
if indexFloor(1)<2 ||indexFloor(2)<2 || indexFloor(1)>(size(img,2)-1)|| indexFloor(2)>(size(img,1)-1)
    instensity=0;
    return;
end

x=indexFloor(1);
z=indexFloor(2);
u=indexInt(1)-indexFloor(1);
v=indexInt(2)-indexFloor(2);
instensity=u*v*img(z,x)+(1-u)*v*img(z+1,x)+u*(1-v)*img(z,x+1)+(1-u)*(1-v)*img(z+1,x+1);
return;
end

function Xw1=TransformToCTProjection(Xw,FDD,theta)
Xw1=FDD.*sin(theta)./(Xw(2).*sin(theta)-Xw(3).*cos(theta)).*Xw;
end
function Xw=TransformToSpatialPosition(i,j,Fi,Fj,FDD,Psize)
Xw=[(i-Fi).*Psize,FDD,(Fj-j).*Psize];
end
