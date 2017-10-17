%Temp

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
for i = 170:175+5
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
    
    imshow(Ig);
    pause(1);
    
    
end