clear all
clc

load('mound_slices.mat');
viewAngle=[-45,0,15];

% =========================================================================
% Plot the void area, solid area, and total area ==========================
% =========================================================================
figure(4)

subplot(1,3,1)
plot(mound.totalArea,mound.Z,'k',mound.moundArea,mound.Z,'r',mound.voidArea,mound.Z,'b')
title('Area vs. Height')
ylabel('Height (m)')
xlabel('Area (m^2)')


subplot(1,3,2)
plot(sqrt(mound.totalArea/pi),mound.Z,'k',sqrt(mound.moundArea/pi),mound.Z,'r',sqrt(mound.voidArea/pi),mound.Z,'b')
title('"Radius" vs. Height')
ylabel('Height (m)')
xlabel('Radius (m)')

subplot(1,3,3)
plot(mound.voidArea./mound.totalArea,mound.Z,'g')
title('Area vs. Height')
ylabel('Height (m)')
xlabel('Area (m^2)')



radius=sqrt(mound.totalArea/pi);

%The camera was shift and this realligns everything as best as possible.
ishift=174;
xshift=mound.Center(ishift,:)-mound.Center(ishift-1,:);
mound.Center(ishift:end,1) = mound.Center(ishift:end,1)-xshift(1);
mound.Center(ishift:end,2) = mound.Center(ishift:end,2)-xshift(2);





% =========================================================================
% Plot the slice centroid versus height
% =========================================================================
figure(3)
subplot(2,2,1)
plot3(mound.Center(:,1)-mound.Center(end,1),mound.Center(:,2)-mound.Center(end,2),mound.Z,'k-o')
axis equal
xlabel('x')
ylabel('y')
zlabel('z')

skip=1;
cnt=1;
for i=1:skip:numel(mound.Z)-skip
%     dz=mound.Z(i)-mound.Z(i+skip);
%     zc(cnt)=(mound.Z(i)+mound.Z(i+skip))/2;
%     XX=mound.Center(i,1)-mound.Center(i+skip,1);%-mound.Center(i+1,1)/1000;
%     YY=mound.Center(i,2)-mound.Center(i+skip,2);%-mound.Center(i+1,2)/1000;
    
    dz=mound.Z(i);%-mound.Z(i+skip);
    zc(cnt)=mound.Z(i);
    XX=mound.Center(i,1)-mound.Center(end,1);%-mound.Center(i+1,1)/1000;
    YY=mound.Center(i,2)-mound.Center(end,2);%-mound.Center(i+1,2)/1000;
    
    ds=sqrt( XX.^2 + YY.^2 );
    theta(cnt)=atan(ds/dz)*180/pi;
    phi(cnt) = atan2(XX,YY)*180/pi;
    cnt=cnt+1;
end

subplot(2,2,2)
plot(theta,zc);

subplot(2,2,3)
plot(phi,zc);



% =========================================================================
% Plot the idealized mound
% =========================================================================
figure(2)
subplot(1,3,1)
hold off
th=linspace(0,2*pi,100);
for i=1:numel(mound.Z)
    plot3(radius(i)*cos(th),radius(i)*sin(th),mound.Z(i)*ones(1,numel(th)),'k-')
    hold on
end
axis([-1,1,-1,1,0,1.5])
axis equal
grid on
title('Idealized mound')
view(viewAngle)


% idealized mound with lean
subplot(1,3,2)
hold off
th=linspace(0,2*pi,100);
for i=1:numel(mound.Z)
    xsh=mound.Center(i,1)-mound.Center(end,1);
    ysh=mound.Center(i,2)-mound.Center(end,2);
    plot3(radius(i)*cos(th)+xsh,radius(i)*sin(th)+ysh,mound.Z(i)*ones(1,numel(th)),'k-')
    hold on
end
axis([-1,1,-1,1,0,1.5])
axis equal
title('Idealized mound with lean')
grid on
view(viewAngle)

% actual mound trace with lean (uses a lot of visual memory, uncomment to plot)
%
% subplot(1,3,3)
% hold off
% for i=1:3:numel(mound.Z)
%     I=mound.Iborder{i};
%     xc=mound.Center(end,1);
%     yc=mound.Center(end,2);
%     if i>=ishift
%        XX=mound.X(I)/1000-xshift(1);
%        YY=mound.Y(I)/1000-xshift(2);
%     else
%         XX=mound.X(I)/1000;
%         YY=mound.Y(I)/1000;
%     end
%    plot3(XX, YY, mound.Z(i)*ones(1,numel(1,I)),'k.');
%    hold on
%     
% end
% axis([-1,1,-1,1,0,1.5])
% axis equal
% grid on
% title('Extracted Border')
% view(viewAngle)


figure(1)
% set(gca,'nextplot','replacechildren');
% vidObj = VideoWriter('voids.avi');
% vidObj.FrameRate=2;
% open(vidObj);



colors='brcmk';
marker='ox+*sdv^<>ph';
for i=1:numel(mound.Z)
    plot(mound.X(mound.Imound{i}),mound.Y(mound.Imound{i}),'k.');
    axis equal
    hold on
    for j=1:min(numel(mound.Void(i).ID),numel(colors)*numel(marker))
      ind=mound.Void(i).pixels{j};
      col=rem(j-1,numel(colors))+1;
      mar=floor((j-1)/numel(colors))+1;
      plot(mound.X(ind),mound.Y(ind),[colors(col),marker(mar)]);
    end
    
%     title(['Image ',num2str(mound.filelist(i),'%03d')])
%     currframe = getframe(1);
%     writeVideo(vidObj,currframe);
    
    pause(0.1)
    hold off
end
% close(vidObj);