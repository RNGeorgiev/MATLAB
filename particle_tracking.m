%%%%%%%%%%%%%TO-DO LIST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%FIND MEMORY LEAK%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%PROPER FINAL FRAME ALGORITHM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%IMADJUST IN OUTPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%GUI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close
clc
set(0,'DefaultFigureWindowStyle','docked');
addpath pwd
addpath ./bfmatlab/
warning('off','all')
%%%%%%%%%%%%%Input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rratio='1.75';
Rdisk='50';
experiment='001';
shape='dumbbell';
colorFloor=11000;
colorCeiling=35000;
partvel=25;
ratio=0.65;
Cpoints=10;
fileStart=1;
fileEnd=650;
fileVis=780;
%%%%%%%%%%%%%Data reading%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_size=1024;
file1=bfopen(strcat('R_',Rratio,'_p_0.10_',experiment,'.nd2'));
%file1=bfopen('R_2.25_p_0.15_005.nd2');
file2=bfopen(strcat('R_',Rratio,'_p_0.10_back.nd2'));
%%%%%%%%%%%%%Metadata Generation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frames_file=length(file1{1,1});
microData=zeros(frames_file,5);
microData(:,2)=[1:frames_file];
frames_back=length(file1{1,2});
microBack=zeros(frames_back,5);
microBack(:,2)=[1:frames_back];
metadata=file1{1,2};
metaback=file2{1,2};
metadataKeys = metadata.keySet().iterator();
metabackKeys = metaback.keySet().iterator();
for i=1:metadata.size()
    key = metadataKeys.nextElement();
    value=metadata.get(key);
    if regexp(key,'Global Y position for position*')==1
        split=strsplit(char(key),'#');
        row=str2num(split{1,2});
        pos=strsplit(char(value),'[');
        pos=strsplit(pos{1,2},']');
        pos=str2double(pos{1,1});
        microData(row,4)=pos;
    elseif regexp(key,'Global X position for position*')==1
        split=strsplit(char(key),'#');
        row=str2num(split{1,2});
        pos=strsplit(char(value),'[');
        pos=strsplit(pos{1,2},']');
        pos=str2double(pos{1,1});
        microData(row,3)=pos;
    elseif regexp(key,'timestamp #*')==1
        split=strsplit(char(key),'#');
        row=str2num(split{1,2});
        microData(row,1)=value;
    elseif regexp(key,'Z position for position, plane #*')==1
        split=strsplit(char(key),'#');
        row=str2num(split{1,2});
        pos=strsplit(char(value),'[');
        pos=strsplit(pos{1,2},']');
        pos=str2double(pos{1,1});
        microData(row,5)=pos;
    end
end
for i=1:metaback.size()
    key = metabackKeys.nextElement();
    value=metaback.get(key);
    if regexp(key,'Global Y position for position*')==1
        split=strsplit(char(key),'#');
        row=str2num(split{1,2});
        pos=strsplit(char(value),'[');
        pos=strsplit(pos{1,2},']');
        pos=str2double(pos{1,1});
        microBack(row,4)=pos;
    elseif regexp(key,'Global X position for position*')==1
        split=strsplit(char(key),'#');
        row=str2num(split{1,2});
        pos=strsplit(char(value),'[');
        pos=strsplit(pos{1,2},']');
        pos=str2double(pos{1,1});
        microBack(row,3)=pos;
    elseif regexp(key,'timestamp #*')==1
        split=strsplit(char(key),'#');
        row=str2num(split{1,2});
        microBack(row,1)=value;
    elseif regexp(key,'Z position for position, plane #*')==1
        split=strsplit(char(key),'#');
        row=str2num(split{1,2});
        pos=strsplit(char(value),'[');
        pos=strsplit(pos{1,2},']');
        pos=str2double(pos{1,1});
        microBack(row,5)=pos;
    end
end
%%%%%%%%%%%%%Initial position of big disk%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
combos=nchoosek(1:Cpoints,3);
j=file1{1,1}{fileStart,1};
imshow(j,[]);hold on; 
for point=1:3 
    Mbigzoom(point,:)=ginput(1);
    plot(Mbigzoom(:,1),Mbigzoom(:,2),'+');
end
close;
clear point;    
[cbig0 Rbig0]=calcCircle(Mbigzoom);
imshow(j(cbig0(2)-Rbig0-10:cbig0(2)+Rbig0+10,cbig0(1)-Rbig0-10:cbig0(1)+Rbig0+10),[]);hold on;
for point=1:Cpoints 
    Mbig0(point,:)=ginput(1);
    plot(Mbig0(:,1),Mbig0(:,2),'+');
end
close;
clear point;
cBs=zeros(length(combos),2);
RBs=zeros(length(combos),1);
for N=1:length(combos)
    [cBs(N,:) RBs(N)]=calcCircle(Mbig0(combos(N,:),:));
end
centerYB=round(mean(cBs(:,2))+cbig0(2)-Rbig0-10);
centerXB=round(mean(cBs(:,1))+cbig0(1)-Rbig0-10);
widthB=2*round((mean(RBs)*2+10)/2);
%%%%%%%%%%%%%Initial position for small disks%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(shape,'disk')==1
    sprintf('Single disk needed.');
elseif strcmp(shape,'dumbbell')==1
    imshow(j,[]);hold on;
    for point=1:3 
        Msmallzoom(point,:)=ginput(1);
        plot(Msmallzoom(:,1),Msmallzoom(:,2),'+');
    end
    close;
    clear point;
    [csmall0 Rsmall0]=calcCircle(Msmallzoom);
    imshow(j(csmall0(2)-Rsmall0-10:csmall0(2)+Rsmall0+10,csmall0(1)-Rsmall0-10:csmall0(1)+Rsmall0+10),[]);hold on;
    for point=1:Cpoints 
        Msmall0(point,:)=ginput(1);
        plot(Msmall0(:,1),Msmall0(:,2),'+');
    end
    close;
    clear point;
    cSs=zeros(length(combos),2);
    RSs=zeros(length(combos),1);
    for N=1:length(combos)
        [cSs(N,:) RSs(N)]=calcCircle(Msmall0(combos(N,:),:));
    end
    centerYS=round(mean(cSs(:,2))+csmall0(2)-Rsmall0-10);
    centerXS=round(mean(cSs(:,1))+csmall0(1)-Rsmall0-10);
    widthS=2*round((mean(RSs)*2+10)/2);
elseif strcmp(shape,'tripod')==1
    for point=1:3 
        Msmallzoom(point,:)=ginput(1);
        plot(Msmallzoom(:,1),Msmallzoom(:,2),'+');
    end
    close;
    clear point;
    [csmall0 Rsmall0]=calcCircle(Msmallzoom);
    imshow(j(csmall0(2)-Rsmall0-10:csmall0(2)+Rsmall0+10,csmall0(1)-Rsmall0-10:csmall0(1)+Rsmall0+10),[]);hold on;
    for point=1:Cpoints 
        Msmall0(point,:)=ginput(1);
        plot(Msmall0(:,1),Msmall0(:,2),'+');
    end
    close;
    clear point;
    cSs1=zeros(length(combos),2);
    RSs1=zeros(length(combos),1);
    for N=1:length(combos)
        [cSs1(N,:) RSs1(N)]=calcCircle(Msmall0(combos(N,:),:));
    end
    centerYS1=round(mean(cSs1(:,2))+csmall0(2)-Rsmall0-10);
    centerXS1=round(mean(cSs1(:,1))+csmall0(1)-Rsmall0-10);
    widthS1=2*round((mean(RSs1)*2+10)/2);
    for point=1:3 
        Msmallzoom(point,:)=ginput(1);
        plot(Msmallzoom(:,1),Msmallzoom(:,2),'+');
    end
    close;
    clear point;
    [csmall0 Rsmall0]=calcCircle(Msmallzoom);
    imshow(j(csmall0(2)-Rsmall0-10:csmall0(2)+Rsmall0+10,csmall0(1)-Rsmall0-10:csmall0(1)+Rsmall0+10),[]);hold on;
    for point=1:Cpoints 
        Msmall0(point,:)=ginput(1);
        plot(Msmall0(:,1),Msmall0(:,2),'+');
    end
    close;
    clear point;
    cSs2=zeros(length(combos),2);
    RSs2=zeros(length(combos),1);
    for N=1:length(combos)
        [cSs2(N,:) RSs2(N)]=calcCircle(Msmall0(combos(N,:),:));
    end
    centerYS2=round(mean(cSs2(:,2))+csmall0(2)-Rsmall0-10);
    centerXS2=round(mean(cSs2(:,1))+csmall0(1)-Rsmall0-10);
    widthS2=2*round((mean(RSs2)*2+10)/2);
end
%%%%%%%%%%%%%Memory allocation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(shape,'disk')==1
    cBig=zeros(fileEnd,2);
    rBig=zeros(fileEnd,1);
    M=zeros(1024,1024,fileEnd-fileStart+1,'uint16');
    data=zeros(fileEnd-fileStart+1,5);
    dxs=zeros(length(microData)-1,1);
    dys=zeros(length(microData)-1,1);
    dxb=zeros(length(microBack)-1,1);
    dyb=zeros(length(microBack)-1,1);
    ts=microData(1:fileEnd,1);
    dts=cat(1,0,ts(2:length(ts))-ts(1:length(ts)-1));
    maxdx=ceil(dts*partvel/ratio);
elseif strcmp(shape,'dumbbell')==1
    cBig=zeros(fileEnd,2);
    cSmall=zeros(fileEnd,2);
    rBig=zeros(fileEnd,1);
    rSmall=zeros(fileEnd,1);
    M=zeros(1024,1024,fileEnd-fileStart+1,'uint16');
    data=zeros(fileEnd-fileStart+1,10);
    dxs=zeros(length(microData)-1,1);
    dys=zeros(length(microData)-1,1);
    dxb=zeros(length(microBack)-1,1);
    dyb=zeros(length(microBack)-1,1);
    ts=microData(1:fileEnd,1);
    dts=cat(1,0,ts(2:length(ts))-ts(1:length(ts)-1));
    maxdx=ceil(dts*partvel/ratio);
elseif strcmp(shape,'tripod')==1
    cBig=zeros(fileEnd,2);
    cSmall1=zeros(fileEnd,2);
    cSmall2=zeros(fileEnd,2);
    rBig=zeros(fileEnd,1);
    rSmall1=zeros(fileEnd,1);
    rSmall2=zeros(fileEnd,1);
    M=zeros(1024,1024,fileEnd-fileStart+1,'uint16');
    data=zeros(fileEnd-fileStart+1,15);
    dxs=zeros(length(microData)-1,1);
    dys=zeros(length(microData)-1,1);
    dxb=zeros(length(microBack)-1,1);
    dyb=zeros(length(microBack)-1,1);
    ts=microData(1:fileEnd,1);
    dts=cat(1,0,ts(2:length(ts))-ts(1:length(ts)-1));
    maxdx=ceil(dts*partvel/ratio);
end
%%%%%%%%%%%%%Image-to-Background Mapping%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for q=1:length(microData)-1
    dxs(q)=round((microData(q+1,3)-microData(q,3))/0.650);
    dys(q)=round((microData(q+1,4)-microData(q,4))/0.650);
end
fracsData=find(dxs);
clear q
for q=1:length(microBack)-1
    dxb(q)=round((microBack(q+1,3)-microBack(q,3))/0.650);
    dyb(q)=round((microBack(q+1,4)-microBack(q,4))/0.650);
end
fracsBack=find(dxb)+1;
map=zeros(fileEnd,2);
map(:,1)=[1:fileEnd];
map(1:fracsData(1),2)=1;
for q=1:length(fracsData)
    if q<length(fracsData)
        map(fracsData(q)+1:fracsData(q+1),2)=fracsBack(q);
    else
        map(fracsData(q):length(map),2)=fracsBack(q);
    end
end
%%%%%%%%%%%%%Image 3D Array Generation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for N=fileStart:fileEnd
    M(:,:,N)=file1{1,1}{map(N,1),1}-file2{1,1}{map(N,2),1};
end
clear N
%%%%%%%%%%%%%Finding Circles%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(shape,'disk')==1
    for N=fileStart:fileEnd
        [centersBig, radiiBig, metricBig]=imfindcircles(...
            imbinarize(...
            M(centerYB-widthB/2:centerYB+widthB/2,centerXB-round((widthB+maxdx(N)*0.2)/2)-maxdx(N):centerXB+round((widthB+maxdx(N)*0.2)/2)-maxdx(N),N)...
            ),...
            [round(mean(RBs))-5 round(mean(RBs))+5],'ObjectPolarity','dark','Sensitivity',0.97,'Method','TwoStage');
        if length(metricBig)>0
            cBig(N,:)=centersBig(find(max(metricBig)),:)+[centerXB-round((widthB+maxdx(N)*0.2)/2)-maxdx(N),centerYB-widthB/2];
            rBig(N)=radiiBig(find(max(metricBig)));
            centerXB=round(cBig(N,1)+dxs(N));
            centerYB=round(cBig(N,2)-dys(N));
        else
            cBig(N,:)=[NaN NaN];
            rBig(N)=NaN;
        end
        sprintf(strcat('Step ',num2str(N),'finished.'))
        jnew=file1{1,1}{N,1};
        if N>fileVis-1
            imshow(jnew(200:824,:),[colorFloor colorCeiling]); hold on;...
            viscircles(cBig(N,:)-[0,200],rBig(N)-2,'LineStyle',':','EnhanceVisibility',0,'LineWidth',1.5); pause(0.1);
        else
            continue
        end
    end
elseif strcmp(shape,'dumbbell')==1
    for N=fileStart:fileEnd
        [centersBig, radiiBig, metricBig]=imfindcircles(...
            imbinarize(...
            M(centerYB-widthB/2:centerYB+widthB/2,centerXB-round((widthB+maxdx(N)*0.2)/2)-maxdx(N):centerXB+round((widthB+maxdx(N)*0.2)/2)-maxdx(N),N)...
            ),...
            [round(mean(RBs))-5 round(mean(RBs))+5],'ObjectPolarity','dark','Sensitivity',0.97,'Method','TwoStage');
        if length(metricBig)>0
            cBig(N,:)=centersBig(find(max(metricBig)),:)+[centerXB-round((widthB+maxdx(N)*0.2)/2)-maxdx(N),centerYB-widthB/2];
            rBig(N)=radiiBig(find(max(metricBig)));
            centerXB=round(cBig(N,1)+dxs(N));
            centerYB=round(cBig(N,2)-dys(N));
        else
            cBig(N,:)=[NaN NaN];
            rBig(N)=NaN;
        end
        [centersSmall, radiiSmall, metricSmall]=imfindcircles(...
            imbinarize(...
            M(centerYS-widthS/2:centerYS+widthS/2,centerXS-round((widthS+maxdx(N)*0.2)/2)-maxdx(N):centerXS+round((widthS+maxdx(N)*0.2)/2)-maxdx(N),N)...
            ),...
            [round(mean(RSs))-5 round(mean(RSs))+5],'ObjectPolarity','dark','Sensitivity',0.98,'Method','TwoStage');
        if length(metricSmall)>0
            cSmall(N,:)=centersSmall(find(max(metricSmall)),:)+[centerXS-round((widthS+maxdx(N)*0.2)/2)-maxdx(N) centerYS-widthS/2];
            rSmall(N)=radiiSmall(find(max(metricSmall)));
            centerXS=round(cSmall(N,1)+dxs(N));
            centerYS=round(cSmall(N,2)-dys(N));
        else
            cSmall(N,:)=[NaN NaN];
            rSmall(N)=NaN;
        end
        sprintf(strcat('Step ',num2str(N),'finished.'))
        jnew=file1{1,1}{N,1};
        if N>fileVis-1
            imshow(jnew(200:824,:),[10129 34942]); hold on;...
                viscircles(cBig(N,:)-[0,200],rBig(N)-2,'LineStyle',':','EnhanceVisibility',0,'LineWidth',1.5); hold on;...
                viscircles(cSmall(N,:)-[0,200],rSmall(N)-2,'LineStyle',':','EnhanceVisibility',0,'LineWidth',1.5); pause(0.1);
        else
            continue
        end
    end
elseif strcmp(shape,'tripod')==1
    for N=fileStart:fileEnd
        [centersBig, radiiBig, metricBig]=imfindcircles(...
            imbinarize(...
            M(centerYB-widthB/2:centerYB+widthB/2,centerXB-round((widthB+maxdx(N)*0.2)/2)-maxdx(N):centerXB+round((widthB+maxdx(N)*0.2)/2)-maxdx(N),N)...
            ),...
            [round(mean(RBs))-5 round(mean(RBs))+5],'ObjectPolarity','dark','Sensitivity',0.97,'Method','TwoStage');
        if length(metricBig)>0
            cBig(N,:)=centersBig(find(max(metricBig)),:)+[centerXB-round((widthB+maxdx(N)*0.2)/2)-maxdx(N),centerYB-widthB/2];
            rBig(N)=radiiBig(find(max(metricBig)));
            centerXB=round(cBig(N,1)+dxs(N));
            centerYB=round(cBig(N,2)-dys(N));
        else
            cBig(N,:)=[NaN NaN];
            rBig(N)=NaN;
        end
        [centersSmall1, radiiSmall1, metricSmall1]=imfindcircles(...
            imbinarize(...
            M(centerYS1-widthS1/2:centerYS1+widthS1/2,...
            centerXS1-round((widthS1+maxdx(N)*0.2)/2)-maxdx(N):centerXS1+round((widthS1+maxdx(N)*0.2)/2)-maxdx(N),N)...
            ),...
            [round(mean(RSs1))-5 round(mean(RSs1))+5],'ObjectPolarity','dark','Sensitivity',0.98,'Method','TwoStage');
        if length(metricSmall1)>0
            cSmall1(N,:)=centersSmall1(find(max(metricSmall1)),:)+[centerXS1-round((widthS1+maxdx(N)*0.2)/2)-maxdx(N) centerYS1-widthS1/2];
            rSmall1(N)=radiiSmall1(find(max(metricSmall1)));
            centerXS1=round(cSmall1(N,1)+dxs(N));
            centerYS1=round(cSmall1(N,2)-dys(N));
        else
            cSmall1(N,:)=[NaN NaN];
            rSmall1(N)=NaN;
        end
        [centersSmall2, radiiSmall2, metricSmall2]=imfindcircles(...
            imbinarize(...
            M(centerYS2-widthS2/2:centerYS2+widthS2/2,...
            centerXS2-round((widthS2+maxdx(N)*0.2)/2)-maxdx(N):centerXS2+round((widthS2+maxdx(N)*0.2)/2)-maxdx(N),N)...
            ),...
            [round(mean(RSs2))-5 round(mean(RSs2))+5],'ObjectPolarity','dark','Sensitivity',0.98,'Method','TwoStage');
        if length(metricSmall2)>0
            cSmall2(N,:)=centersSmall1(find(max(metricSmall2)),:)+[centerXS1-round((widthS2+maxdx(N)*0.2)/2)-maxdx(N) centerYS1-widthS2/2];
            rSmall2(N)=radiiSmall2(find(max(metricSmall2)));
            centerXS2=round(cSmall2(N,1)+dxs(N));
            centerYS2=round(cSmall2(N,2)-dys(N));
        else
            cSmall2(N,:)=[NaN NaN];
            rSmall2(N)=NaN;
        end
        sprintf(strcat('Step ',num2str(N),'finished.'))
        jnew=file1{1,1}{N,1};
        if N>fileVis-1
            imshow(jnew(200:824,:),[10129 34942]); hold on;...
                viscircles(cBig(N,:)-[0,200],rBig(N)-2,'LineStyle',':','EnhanceVisibility',0,'LineWidth',1.5); hold on;...
                viscircles(cSmall1(N,:)-[0,200],rSmall1(N)-2,'LineStyle',':','EnhanceVisibility',0,'LineWidth',1.5); hold on;...
                viscircles(cSmall2(N,:)-[0,200],rSmall2(N)-2,'LineStyle',':','EnhanceVisibility',0,'LineWidth',1.5); pause(0.1);
        else
            continue
        end
    end
end
close all;
%%%%%%%%%%%%%Data Writting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(shape,'disk')==1
    data(:,1)=transpose([fileStart:fileEnd]);
    data(:,2)=ts(fileStart:fileEnd);
    data(:,3)=cBig(fileStart:fileEnd,1);
    data(:,4)=cBig(fileStart:fileEnd,2);
    data(:,5)=rBig(fileStart:fileEnd);
    data(:,6)=microData(fileStart:fileEnd,3);
    data(:,7)=microData(fileStart:fileEnd,4);
    data(:,8)=microData(fileStart:fileEnd,5);
    data(:,9)=data(:,3)*ratio-(microData(fileStart:fileEnd,3)-microData(1,3));
    dlmwrite(strcat('data_',Rdisk,'_',experiment,'_',num2str(fileStart),'to',num2str(fileEnd),'.dat'),data,'delimiter','\t','newline','pc');
elseif  shape=='dumbbell'
    data(:,1)=transpose([fileStart:fileEnd]);
    data(:,2)=ts(fileStart:fileEnd);
    data(:,3)=cBig(fileStart:fileEnd,1);
    data(:,4)=cBig(fileStart:fileEnd,2);
    data(:,5)=rBig(fileStart:fileEnd);
    data(:,6)=cSmall(fileStart:fileEnd,1);
    data(:,7)=cSmall(fileStart:fileEnd,2);
    data(:,8)=rSmall(fileStart:fileEnd);
    data(:,9)=microData(fileStart:fileEnd,3);
    data(:,10)=microData(fileStart:fileEnd,4);
    data(:,11)=microData(fileStart:fileEnd,5);
    data(:,12)=sqrt((data(:,6)-data(:,3)).^2+(data(:,7)-data(:,4)).^2);
    data(:,13)=acos(-(data(:,6)-data(:,3))./data(:,12));
    dlmwrite(strcat('data_',Rratio,'_',experiment,'_',num2str(fileStart),'to',num2str(fileEnd),'.dat'),data,'delimiter','\t','newline','pc');
elseif  shape=='tripod'
    data(:,1)=transpose([fileStart:fileEnd]);
    data(:,2)=ts(fileStart:fileEnd);
    data(:,3)=cBig(fileStart:fileEnd,1);
    data(:,4)=cBig(fileStart:fileEnd,2);
    data(:,5)=rBig(fileStart:fileEnd);
    data(:,6)=cSmall1(fileStart:fileEnd,1);
    data(:,7)=cSmall1(fileStart:fileEnd,2);
    data(:,8)=rSmall1(fileStart:fileEnd);
    data(:,9)=cSmall2(fileStart:fileEnd,1);
    data(:,10)=cSmall2(fileStart:fileEnd,2);
    data(:,11)=rSmall2(fileStart:fileEnd);
    data(:,12)=microData(fileStart:fileEnd,3);
    data(:,13)=microData(fileStart:fileEnd,4);
    data(:,14)=microData(fileStart:fileEnd,5);
    data(:,15)=sqrt((data(:,6)-data(:,3)).^2+(data(:,7)-data(:,4)).^2);
    data(:,16)=sqrt((data(:,9)-data(:,3)).^2+(data(:,10)-data(:,4)).^2);
    data(:,17)=acos(-(data(:,6)-data(:,3))./data(:,12));
    data(:,18)=acos(((data(:,6)-data(:,3))*(data(:,9)-data(:,3))+(data(:,7)-data(:,4))*(data(:,10)-data(:,4)))./(data(:,12).*data(:,13)));
    dlmwrite(strcat('data_',Rratio,'_',experiment,'_',num2str(fileStart),'to',num2str(fileEnd),'.dat'),data,'delimiter','\t','newline','pc');
end
%%%%%%%%%%%%%Data Plotting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
liny=log(tan(data(:,13)./2));
%p=polyfit(data(:,2),liny,1);
%set(0,'DefaultFigureWindowStyle','normal');
%scatter(data(:,2)-data(1,2),data(:,10)); hold on;...
%plot(data(:,2)-data(1,2),atan(exp((data(:,2)-data(1,2)).*p(1))*tan(data(1,10)/2))*2);
%figure
%scatter(-(data(:,2)-data(1,2)).*p(1),data(:,10)./(2*atan(exp(p(2))))); hold on; ...
%plot([0:0.1:ceil(-data(length(data(:,2)),2)*p(1))+1],atan(exp(-[0:0.1:ceil(-data(length(data(:,2)),2)*p(1))+1])*tan(data(1,10)/2))*2./(2*atan(exp(p(2)))));