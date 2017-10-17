%minecraft plot

clc
clear all
% close all
load('MC.mat');

MCIxc=(MCIx(1:end-1)+MCIx(2:end))/2;
MCIyc=(MCIy(1:end-1)+MCIy(2:end))/2;
[xx,yy]=meshgrid(MCIxc,MCIyc);
xx=xx'; yy=yy';zz(1,1,:)=MCIz;

xx=repmat(xx,[1,1,size(Im_MC,3)]);
yy=repmat(yy,[1,1,size(Im_MC,3)]);
zz=repmat(zz,[size(Im_MC,1),size(Im_MC,2),1]);


I=find(Im_MC==1);
plot3(xx(I),yy(I),zz(I),'k.')

% figure 
% hold off


% for kk=1:size(Im_MC,3)
%    for jj=1:numel(MCIyc)
%       for ii=1:numel(MCIxc)
%         if (Im_MC(ii,jj,kk)==1)
%             plot3( xx(ii,jj), yy(ii,jj), max(MCIz) - MCIz(kk),'ko' ) ;
%             hold on
%         end
%       end
%    end
%    
%    pause(0.1)
% end
