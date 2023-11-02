clear;close all;clc;
load result_90projections.mat

CLFDKMeasossart=Measure_Quality(shepp(:,:,128),CLFDK(:,:,128),{'RMSE','MSSIM','CC','UQI','error_norm'});
% CTFDKMeasossart=Measure_Quality(shepp(:,:,128),CTFDK(:,:,128),{'RMSE','MSSIM','CC','UQI','error_norm'});
% CTASDPOCSMeasossart=Measure_Quality(shepp(:,:,128),CTASDPOCS(:,:,128),{'RMSE','MSSIM','CC','UQI','error_norm'});
% CTPTARPOCSMeasossart=Measure_Quality(shepp(:,:,128),CTPTARPOCS(:,:,128),{'RMSE','MSSIM','CC','UQI','error_norm'});
% CLPTARPOCSSMeasossart=Measure_Quality(shepp(:,:,128),CLPTARPOCS(:,:,128),{'RMSE','MSSIM','CC','UQI','error_norm'});
% Measossart=[CLFDKMeasossart;CTFDKMeasossart;CTASDPOCSMeasossart;CTPTARPOCSMeasossart;CLPTARPOCSSMeasossart];

c1=Measure_Quality(abcd(1,:),abcd(4,:),{'RMSE','MSSIM','CC','UQI','error_norm'});
c2=Measure_Quality(abcd(2,:),abcd(4,:),{'RMSE','MSSIM','CC','UQI','error_norm'});
c3=Measure_Quality(abcd(3,:),abcd(4,:),{'RMSE','MSSIM','CC','UQI','error_norm'});
cc=[c1;c2;c3];
ASD_POCS_PT_Max=max(ASD_POCS_PT,[],'all');
ASD_POCS_PT_Min=min(ASD_POCS_PT,[],'all');
imwrite(uint16((ASD_POCS_PT(:,:,1)-ASD_POCS_PT_Min)./(ASD_POCS_PT_Max-ASD_POCS_PT_Min)*65535),'ASD_POCS_PT.tif')
for i=2:399
    imwrite(uint16((ASD_POCS_PT(:,:,i)-ASD_POCS_PT_Min)./(ASD_POCS_PT_Max-ASD_POCS_PT_Min)*65535),'ASD_POCS_PT.tif','WriteMode','append');
end

FDK_AR_Max=max(FDK_AR,[],'all');
FDK_AR_Min=min(FDK_AR,[],'all');
imwrite(uint16((FDK_AR(:,:,1)-FDK_AR_Min)./(FDK_AR_Max-FDK_AR_Min)*65535),'FDK_AR.tif')
for i=2:399
    imwrite(uint16((FDK_AR(:,:,i)-FDK_AR_Min)./(FDK_AR_Max-FDK_AR_Min)*65535),'FDK_AR.tif','WriteMode','append');
end

FDK_Orignal_Max=max(FDK_Orignal,[],'all');
FDK_Orignal_Min=min(FDK_Orignal,[],'all');
imwrite(uint16((FDK_Orignal(:,:,1)-FDK_Orignal_Min)./(FDK_Orignal_Max-FDK_Orignal_Min)*65535),'FDK_Orignal.tif')
for i=2:399
    imwrite(uint16((FDK_Orignal(:,:,i)-FDK_Orignal_Min)./(FDK_Orignal_Max-FDK_Orignal_Min)*65535),'FDK_Orignal.tif','WriteMode','append');
end

FDK_PT_Max=max(FDK_PT,[],'all');
FDK_PT_Min=min(FDK_PT,[],'all');
imwrite(uint16((FDK_PT(:,:,1)-FDK_PT_Min)./(FDK_PT_Max-FDK_PT_Min)*65535),'FDK_PT.tif')
for i=2:399
    imwrite(uint16((FDK_PT(:,:,i)-FDK_PT_Min)./(FDK_PT_Max-FDK_PT_Min)*65535),'FDK_PT.tif','WriteMode','append');
end

PATA_RPOCS_Max=max(PATA_RPOCS,[],'all');
PATA_RPOCS_Min=min(PATA_RPOCS,[],'all');
imwrite(uint16((PATA_RPOCS(:,:,1)-PATA_RPOCS_Min)./(PATA_RPOCS_Max-PATA_RPOCS_Min)*65535),'PATA_RPOCS.tif')
for i=2:399
    imwrite(uint16((PATA_RPOCS(:,:,i)-PATA_RPOCS_Min)./(PATA_RPOCS_Max-PATA_RPOCS_Min)*65535),'PATA_RPOCS.tif','WriteMode','append');
end

V_recon_Max=max(V_recon,[],'all');
V_recon_Min=min(V_recon,[],'all');
imwrite(uint16((V_recon(:,:,1)-V_recon_Min)./(V_recon_Max-V_recon_Min)*65535),'V_recon.tif')
for i=2:200
    imwrite(uint16((V_recon(:,:,i)-V_recon_Min)./(V_recon_Max-V_recon_Min)*65535),'V_recon.tif','WriteMode','append');
end
vvv=imread('FDK_0.103_1596_1148_359.tif',762);
% aaa=FDK_Orignal(180,:,200);
% bbb=FDK_PT(180,:,315);
% ccc=vvv(720,:);
% ddd=single(vvv(720,1:4:end));
% abcd=[(aaa-min(aaa))./(max(aaa)-min(aaa));(bbb-min(bbb))./(max(bbb)-min(bbb));(ddd-min(ddd))./(max(ddd)-min(ddd))];
% abcd(4,:)=mean(abcd);

minutes=abcd(1:3,:)-[abcd(4,:);abcd(4,:);abcd(4,:)];


V_recon(:,:,125)

defc=abcd(:,1:5:end);


