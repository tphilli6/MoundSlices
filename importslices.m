% assuming that the pipe at the top is standard pvc schedule 40
% the wall thickness is about 5 pixels and the outer diameter is about 57
% pixels resulting in a wt/od = 0.0877
% comparing to PVC schedule size charts outer diameter of 1.25 inches
% actual outer diamter is 1.660, min wall thickness is 0.140 inches with a
% wt/od = 0.140/1.660 = 0.843
%
% https://formufit.com/pages/pvc-pipe-size-dimensions-chart
%
% The assumed pixel scaling measurement is 1.660 in per 57 pixels 
% 0.0296 in/pixel +- 0.00055 mm/pixel
% 0.753 mm/pixel +- 0.013 mm/pixel
%
% 1cm = 13.3 pixels


clc
clear all

createVideo=0;

data=importdata('traced/filelist.txt');
numfilelist=[1:5:100,105:5:800,811:10:1331];

mound.filelist=numfilelist;

imoundtop=1331;

[X, Y]=meshgrid(0:3264-1,0:2448-1);
X=X'; Y=Y';
X=X*0.753;
Y=Y*0.753;
pixelArea=0.753*0.753/1000^2;
mound.X=X;
mound.Y=Y;

MCscale=10;%mm
MCn=round(MCscale/.753,0);
nx=(floor(3264/MCn)+1)*MCn;
ny=(floor(2448/MCn)+1)*MCn;
MCIx= round(1:MCn:nx,0);
MCIy= round(1:MCn:ny,0);

nx=max(MCIx);
ny=max(MCIy);
MCnumfilelist=[1:10:100,110:10:800,811:10:1331];
Im_MC=zeros(numel(MCIx)-1,numel(MCIy)-1,numel(MCnumfilelist));
MCIz=MCnumfilelist;

nframe=numel(numfilelist);

if createVideo
    vidObj = VideoWriter('peaks.avi');
    vidObj.FrameRate=2;
    open(vidObj);

    % mod=zeros(1,nframe);

    figure(1)
    set(gca,'nextplot','replacechildren');
end

cnt=1;
MCcnt=1;
for i = 1:nframe
    ilabel=str2num(data{i}(6:8)); %pull number from file name
    iloc=imoundtop-ilabel+1; %location from the top in mm;
    
    %Import Image
    if numfilelist(i)<1000
        file=['traced/image',num2str(numfilelist(i),'%03d'),'.jpg'];
    else
      file=['traced/image',num2str(numfilelist(i),'%04d'),'.jpg'];
    end
    if ~exist(file,'file')
        continue
    end
    Ig=imread(file);
    
    
    [Im_out, Iborder, Irmean, Irad, Iprocessed] = process_trace_line(Ig,0);
    
    % Compute mound info
    Inotwhitewhite=uint32(find(Im_out(:,:,3)<255 & Im_out(:,:,2)<255 & Im_out(:,:,1)<255));
    for ii=1:3
        temp=Im_out(:,:,ii); Im_out_mean(ii)=mean(temp(Inotwhitewhite));
    end

    Iwhite=uint32(find(Im_out(:,:,3)>max(mean(Im_out_mean),90) ...
              & Im_out(:,:,2)>max(mean(Im_out_mean),90) ...
              & Im_out(:,:,1)>max(mean(Im_out_mean),90)));
    Iwhite=intersect(Iwhite,Inotwhitewhite);
    
    Im_bw=uint8( 255*ones(size(Im_out(:,:,1))));
    Im_bw(Inotwhitewhite)=0;
    Im_bw(Iwhite)=255;
    
    %mound area
    mound.totalArea(cnt)=numel(Inotwhitewhite)*pixelArea; %number of pixels times pixelArea
    mound.moundArea(cnt)=(numel(Inotwhitewhite)-numel(Iwhite))*pixelArea;
    mound.voidArea(cnt)=numel(Iwhite)*pixelArea;
    mound.voidRatio(cnt)=mound.voidArea(cnt)/mound.totalArea(cnt);
    
    mound.Z(cnt)=(imoundtop - numfilelist(i))/1000;
    mound.Center(cnt,1)=mean(X(Iprocessed))/1000;
    mound.Center(cnt,2)=mean(Y(Iprocessed))/1000;
    
    mound.Iborder{cnt}=uint32(Iborder);
    mound.Imound{cnt}=uint32( setdiff(Inotwhitewhite, Iwhite) );
    mound.Ivoid{cnt}=uint32( Iwhite );
    
    mound.Void(cnt) = findVoidSizes(mound.Imound{cnt},mound.Ivoid{cnt},mound.Iborder{cnt},size(X));
    
    if any(numfilelist(i)==MCnumfilelist)
      [Im_MC(:,:,MCcnt)] = minecraft_reduce(Im_bw(1:nx,1:ny),MCIx, MCIy);
      MCcnt=MCcnt+1;
    end
    
    
%     mound.Im_out(:,:,i)=rgb2gray(Im_out);
%     mound.Iborder{i}=uint32(Iborder);
%     mound.Irmean(i,:)=Irmean;
%     mound.Irad(i)=Irad;
%     mound.Imound{i}=uint32(Iprocessed);
    
    
    
%     Im_out=rgb2gray(Im_out);
%     Im_out(Iprocessed)=imadjust(Im_out(Iprocessed),[0.1,0.9],[]);

%     Im_out=uint8(ones(size(Im_out)));

    % Add cm measurement
    il=max(Irmean(1)-1,1); ih=min(Irmean(1)+1,size(Im_out,1));
    jl=max(Irmean(2)-6,1); jh=min(Irmean(2)+7,size(Im_out,2));
    
    Im_out(il:ih, jl:jh,:)=0;
    Im_out(il:ih, jl:jh,1)=255;
    
    
    
%     Im_out(Irmean(2)-6:Irmean(2)+7,Irmean(1)-1:Irmean(1)+1,2)=200;
%     Im_out(Irmean(2)-6:Irmean(2)+7,Irmean(1)-1:Irmean(1)+1,3)=150;

    ilow = max(Irmean(1)-Irad,1);
    ihigh = min(Irmean(1)+Irad,size(Im_out,1));
    jlow = max(Irmean(2)-Irad,1);
    jhigh = min(Irmean(2)+Irad,size(Im_out,2));
    
    Im_out=Im_out(ilow:ihigh,jlow:jhigh,:);
%     Im_out = imresize(Im_out,[600,600]);

    Im_outHSV=rgb2hsv(Im_out);
    
    Inotwhitewhite=uint32(find(Im_out(:,:,3)<255 & Im_out(:,:,2)<255 & Im_out(:,:,1)<255));
    for ii=1:3
        temp=Im_out(:,:,ii); Im_out_mean(ii)=mean(temp(Inotwhitewhite));
    end

%     Iwhite=find(Im_outHSV(:,:,3)>0.4 & Im_outHSV(:,:,2)<0.45);
    Iwhite=uint32(find(Im_out(:,:,3)>max(mean(Im_out_mean),90) ...
              & Im_out(:,:,2)>max(mean(Im_out_mean),90) ...
              & Im_out(:,:,1)>max(mean(Im_out_mean),90)));
%     Iwhite=find(Im_out(:,:,3)>100 | Im_out(:,:,2)>100 | Im_out(:,:,1)>100 ...
%                 & (Im_outHSV(:,:,3)>0.4 & Im_outHSV(:,:,2)<0.45) );
%         Iwhite=find(Im_outHSV(:,:,2)<0.45);
%     Iwhite=intersect(Iwhite,Inotwhitewhite);
    


    
    Im_outV=uint8(ones(size(Im_out(:,:,3))));
    Im_outV(Iwhite)=uint8(255);
%     Im_outV(Inotwhite)=uint(0);
%     mound.Inotsoil{i}=uint32(Iwhite);
%     mound.Isoil{i}=uint32(setdiff(Inotwhitewhite,Iwhite));

    if createVideo
        figure(1)
        subplot(1,2,1)
        Ig600 = imresize(Ig(ilow:ihigh,jlow:jhigh,:),[600,600]);
        imshow(Ig600);
        text(3,10,['Image ',num2str(numfilelist(i),'%03d')])

        subplot(1,2,2)
        Im_outV600 = imresize(Im_outV,[600,600]);
        imshow(Im_outV600);
        text(3,10,['Image ',num2str(numfilelist(i),'%03d')])

        currframe = getframe(1);
        writeVideo(vidObj,currframe);
    end
    
%     figure(2)
%     hold on
%     plot3(X(Iborder),Y(Iborder),iloc*ones(size(Iborder)),'r.')
    
    
    
%     Ig=imresize(Ig,0.25); %Resize
%     Ig=imgaussfilt(Ig,4); %Filter
% 
% 
% figure(1)
% Zr=Ig; Zr(:,:,2)=0;Zr(:,:,3)=0;
% Zg=Ig; Zg(:,:,1)=0;Zg(:,:,3)=0;
% Zb=Ig; Zb(:,:,1)=0;Zb(:,:,2)=0;
% subplot(2,3,1); imshow(Zr);
% subplot(2,3,2); imshow(Zg);
% subplot(2,3,3); imshow(Zb);
% subplot(2,3,4); imhist(Ig(:,:,1));
% subplot(2,3,5); imhist(Ig(:,:,2));
% subplot(2,3,6); imhist(Ig(:,:,3));
% 
% % Ig=min(Ig,[],3); %choose color closest to saturation
% % Ig=max(Ig,[],3); %choose color closest to white
% % Ig = rgb2gray(Ig); %grayscale
% Ig = Ig(:,:,3); %blue color
% % Sum of all colors
% % Ig=sum(Ig,3); 
% % Ig=uint8(255*Ig./max(max(Ig)));
% 
% %Find minimum in historagram and locate white areas
% [counts, binLocations] = imhist(Ig);
% [mincounts,cutoff]=min(counts(125:200));
% [Ilow]=find(counts(100:200)<5*mincounts,1);
% II=find(Ig>100+Ilow);
% 
% % Ig2=255*ones(size(Ig));%For black and white image
% Ig2=Ig;%for shaded image with highlighted voids
% Ig2(II)=0;
% 
% 
% % Ig2=imgaussfilt(Ig2,4); %Filter for smoothness
% 
% 
% figure(2)
% subplot(2,2,1)
% imshow(Ig)
% 
% subplot(2,2,2)
% imshow(Ig2)
% 
% subplot(2,2,3)
% imhist(Ig)
% 
% subplot(2,2,4)
% imhist(Ig2)
% 
% data(:,:,cnt) = Ig2;
% 
% area(cnt)=numel(II)/numel(Ig2);
cnt=cnt+1;
end

save('mound_slices.mat','mound');
save('MC.mat','Im_MC','MCIx','MCIy','MCIz');

% Close the file.
if createVideo
  close(vidObj);
end

% nx=size(data);
% nx(1:2)=nx(1:2)/max(nx(1:2));
% nx(3)=nx(3)*0.01;
% y=linspace(0,nx(1),size(data,1));
% x=linspace(0,nx(2),size(data,2));
% z=linspace(0,nx(3),size(data,3));
% [X,Y,Z]=meshgrid(x,y,z);


