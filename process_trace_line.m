function [Ig_out, Iborder, Irmean, Irrad, Itest] = process_trace_line(Ig,MOD)

% Ig=imread(['traced/image1111.jpg']);
% MOD=0;

Ighsv=rgb2hsv(Ig);

rth=.98; gth=0.78; bth=0.78;
I=uint32(find(Ighsv(:,:,1)<rth & Ighsv(:,:,2)<gth & Ighsv(:,:,3)<bth));
Ir=uint32(find(Ighsv(:,:,1)>=rth & Ighsv(:,:,2)>=gth & Ighsv(:,:,3)>=bth)) ;


if isempty(Ir)
   Ig_out=Ig;
   Iborder=[];
   Irmean=[size(Ig,1)/2, size(Ig,2)/2];
   Irrad=[max(Irmean)];
   Itest=[];
   return
end

%Extract immediate region of trace
[Irx, Iry]=meshgrid(1:size(Ig,1),1:size(Ig,2));
Irx=Irx'; Iry=Iry';

Irmeanx=floor(mean(Irx(Ir))+1);
Irmeany=floor(mean(Iry(Ir))+1);
Irmean=[Irmeanx,Irmeany];

Irrad=floor(max(sqrt( (Irx(Ir)-Irmeanx).^2 + (Iry(Ir)-Irmeany).^2 ))+1)+10;

Irad=uint32(find( sqrt( (Irx-Irmeanx).^2 + (Iry-Irmeany).^2 ) <= Irrad ));


Ig2z=zeros(size(Ig(:,:,1)));
Ighsvr=Ighsv(:,:,1); Ig2r=Ig2z; Ig2r(Irad)=Ighsvr(Irad);
Ighsvg=Ighsv(:,:,2); Ig2g=Ig2z; Ig2g(Irad)=Ighsvg(Irad);
Ighsvb=Ighsv(:,:,3); Ig2b=Ig2z; Ig2b(Irad)=Ighsvb(Irad);
Ig2hsv(:,:,1)=Ig2r;
Ig2hsv(:,:,2)=Ig2g;
Ig2hsv(:,:,3)=Ig2b;

% temp=Ig(:,:,1); Ig2r=Ig2z; Ig2r(Irad)=temp(Irad);
% temp=Ig(:,:,2); Ig2r=Ig2z; Ig2r(Irad)=temp(Irad);
temp=Ig(:,:,3); Ig2r=Ig2z; Ig2r(Irad)=temp(Irad);
Ig2(:,:,1)=Ig2r;
Ig2(:,:,2)=Ig2g;
Ig2(:,:,3)=Ig2b;


% if MOD==1
%     Ir=find(Ig2hsv(:,:,1)>=rth | ( Ig2hsv(:,:,2)>=gth & Ig2hsv(:,:,3)>=bth)  );
% end
% Ir=find(Ig2(:,:,1)>=rth);
Ir=uint32(find(Ig2(:,:,1)>=rth & Ig2(:,:,2)>=gth & Ig2(:,:,3)>=bth)) ;


% Ir=find(Ig2hsv(:,:,1)>=rth | Ig2hsv(:,:,2)>=gth);
% if opt==1
% Ir=find( (Ig2hsv(:,:,1)>=rth & Ig2hsv(:,:,2)>=gth) | ( Ig2(:,:,1)>200 & Ig2(:,:,2)<100 & Ig2(:,:,3)<100) );
% end
% figure(1)
% imshow(Ig2hsv);


Iborder=Ir;%find(Ig2hsv(:,:,2)>=gth);

Igr=zeros(size(Ig(:,:,1)));
Igg=Igr; Igb=Igr;

Igr(Iborder)=255;
Igg(Iborder)=0; Igb(Iborder)=0;

Igwhite(:,:,1)=Igr; Igwhite(:,:,2)=Igg; Igwhite(:,:,3)=Igb;

% extract pixels within the trace (horizontal direction)
Ig_out=uint8(255*ones(size(Ig)));
Ig_out_test=zeros(size(Ig(:,:,1)));



%close small gaps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
smallGap=7;
isred=uint32(find(Igr==255));

for i = 1:numel(isred)
  ix=Irx(isred(i));
  iy=Iry(isred(i));
  il=max(ix-1,1); ih=min(ix+1,size(Igr,1));
  jl=max(iy-1,1); jh=min(iy+1,size(Igr,2));
  subArray=Igr( il:ih, jl:jh ); 

  % Adds horizontal and vertical connections (not just diagonal)
  %NW
  if subArray(1,1)==255 && subArray(2,1)~=255 && subArray(1,2)~=255
      subArray(2,1)=255;
      subArray(1,2)=255;
  end
  %NE
  if size(subArray,2)==3
      if subArray(1,3)==255 && subArray(2,3)~=255 && subArray(1,2)~=255
          subArray(1,2)=255;
          subArray(2,3)=255;
      end
  end
  %SE
  if size(subArray,2)==3 && size(subArray,1)==3
      if subArray(3,3)==255 && subArray(2,3)~=255 && subArray(3,2)~=255
          subArray(2,3)=255;
          subArray(3,2)=255;
      end
  end
  %SW
  if size(subArray,1)==3
      if subArray(3,1)==255 && subArray(2,1)~=255 && subArray(3,2)~=255
          subArray(2,1)=255;
          subArray(3,2)=255;
      end
  end

  deadEnd=1;
  %vertical cross
  if size(subArray,1)==3
      if any(subArray(1,:)==255) && any(subArray(3,:)==255)
        deadEnd=0;
      end
  end
  %horizontal cross
  if size(subArray,2)==3
      if any(subArray(:,1)==255) && any(subArray(:,3)==255)
        deadEnd=0;
      end
  end
  %horizontal cross
  if size(subArray,1)==3
      if any(subArray(1,:)==255) && any(subArray(3,:)==255)
        deadEnd=0;
      end
  end

  Igr( il:ih, jl:jh )=subArray;

%   movingRight=sum(subArray(:,3)==255)>=2;
%   movingLeft=sum(subArray(:,1)==255)>=2;
%   movingDown=sum(subArray(3,:)==255)>=2;
%   movingUp=sum(subArray(1,:)==255)>=2;


  if deadEnd
    % Not connected (maybe)
    exitCheck=0;
    loopCnt=0;

    nxm=2;
    nxp=2;
    nym=2;
    nyp=2;

    while exitCheck==0 
        clear colCheck rowCheck
        iil=max(ix-nxm,1); iih=min(ix+nxp,size(Igr,1));
        jjl=max(iy-nym,1); jjh=min(iy+nyp,size(Igr,2));
        nx0=nxm-(iil - (ix-nxm)   )+1;
        ny0=nym-(jjl - (iy-nym)   )+1;

        subArray2=Igr( iil:iih, jjl:jjh );
        for ii=1:size(subArray2,2)
           colCheck(1,ii)=any(subArray2(:,ii)==255);
        end
        for ii=1:size(subArray2,1)
           rowCheck(ii,1)=any(subArray2(ii,:)==255);
        end

        if sum(colCheck)~=size(subArray2,2) || sum(rowCheck)~=size(subArray2,1)
            cntC=0;
            cntR=0;
            if colCheck(1)==0
                nymNew=nym+1;
                cntC=cntC+1;
            else 
                nymNew=nym;
            end
            if colCheck(end)==0
                nypNew=nyp+1;
                cntC=cntC+1;
            else
                nypNew=nyp;
            end
            if rowCheck(1)==0
                nxmNew=nxm+1;
                cntR=cntR+1;
            else
                nxmNew=nxm;
            end
            if rowCheck(end)==0
                nxpNew=nxp+1;
                cntR=cntR+1;
            else
                nxpNew=nxp;
            end

            if cntC==0 %a gap is present
                sA2test=(subArray2==255);
                [~,IrowFilled]=max(sum(sA2test,2));
                IcolEmpty=find(colCheck==0);
                subArray2(nx0,IcolEmpty)=255;
                Igr( iil:iih, jjl:jjh  )=subArray2;
                exitCheck=1;
            end

            if cntR==0 %a gap is present
                sA2test=(subArray2==255);
                [~,IcolFilled]=max(sum(sA2test,1));
                IrowEmpty=find(rowCheck==0);
                subArray2(IrowEmpty,ny0)=255;
                Igr( iil:iih, jjl:jjh  )=subArray2;
                exitCheck=1;
            end

            nxp=nxpNew;
            nxm=nxmNew;
            nyp=nypNew;
            nym=nymNew;

        end %end if

        % Check if all pixels are connected %%%%%%%%%%%%%%%%%%%%%%%%%%%
        sA2connTest=zeros(size(subArray2));
        sA2connTest(nx0, ny0)=1;
        change=1;
        while change
            change=0;

            for jj=1:size(subArray2,2)
                for ii=1:size(subArray2,1)
                    il=max(ii-1,1); ih=min(ii+1,size(subArray2,1));
                    jl=max(jj-1,1); jh=min(jj+1,size(subArray2,2));
                    if subArray2(ii,jj)==255 ...
                       && sA2connTest(ii,jj)~=1 ...
                       && any(any(sA2connTest(il:ih,jl:jh)==1))
                        sA2connTest(ii,jj)=1;
                        change=1;
                    end

                end
            end
        end %end while

        notConnected=(sA2connTest==0 & subArray2==255);
        [sAx, sAy]=meshgrid(1:size(subArray2,1),1:size(subArray2,2));
        sAx=sAx'; sAy=sAy';
        Inc=uint32(find(notConnected));
        r=(sAx(Inc)-sAx(nx0,ny0)).^2 + (sAy(Inc)-sAy(nx0,ny0)).^2;

        if any(any(notConnected))
            for ii=1:numel(Inc)
                xfillDir=(sAx(Inc(ii))-(nx0))/abs(sAx(Inc(ii))-(nx0));
                if isnan(xfillDir); xfillDir=1; end % if on the same line
                Ixfill=nx0:xfillDir:sAx(Inc(ii));
                yfillDir=(sAy(Inc(ii))-(ny0))/abs(sAy(Inc(ii))-(ny0));
                if isnan(yfillDir); yfillDir=1; end %if on the same line
                Iyfill=ny0:yfillDir:sAy(Inc(ii));
                subArray2(Ixfill,Iyfill)=255;
            end
            Igr( iil:iih, jjl:jjh  )=subArray2;
            exitCheck=1;
        end


        if loopCnt>smallGap
            exitCheck=1;
        else
            loopCnt=loopCnt+1;
        end

    end %end while





  end
end

% Connect boundary if red is found on boundary %%%%%%%%%%%%%%%%%%%%%%%%%%%%
Itest=find(Igr(:,1)); Igr(min(Itest):max(Itest),1)=255;
Itest=find(Igr(:,end)); Igr(min(Itest):max(Itest),end)=255;
Itest=find(Igr(1,:)); Igr(1,min(Itest):max(Itest))=255;
Itest=find(Igr(end,:)); Igr(end,min(Itest):max(Itest))=255;


% Trace the connected path and exclude non-connected sections %%%%%%%%%%%%%
isred=uint32(find(Igr==255));
isredConnected=zeros(size(Igr));
conID=1;
xy2i=reshape( 1:numel(Irx),size(Irx));
lines=zeros(size(Igr));

while any(isredConnected(isred)==0)
  i=find(isredConnected(isred)==0,1,'first');
  
  queue=isred(i);
 
  while ~isempty(queue)
      
      queueNew=[];
      for jj=1:length(queue)
          pixel=queue(jj);
          isredConnected(pixel) = conID; %assign connection ID

          %pull neighbors
          iy=Iry(pixel);
          ix=Irx(pixel); %<===== changing this resulted in an infinite loop
          il=max(ix-1,1); ih=min(ix+1,size(Igr,1));
          jl=max(iy-1,1); jh=min(iy+1,size(Igr,2));
          subArray=Igr( il:ih, jl:jh );
          xy2iSA=xy2i( il:ih, jl:jh);
          
          Ineigh=find(subArray==255);
          Ineigh5=find(Ineigh~=5 & isredConnected(xy2iSA(Ineigh))==0  );
          Ineigh=Ineigh(Ineigh5);

          
          %Store connection ids for neighbors
          isredConnected(xy2iSA(Ineigh)) = conID;
          queueNew=[queueNew;xy2iSA(Ineigh)];

          
          %remove current pixel from queue
%           Iremove=find(queue~=pixel);
%           queue=queue(Iremove);



      end
      queue=queueNew;
      

  end
  
  Icon=find(isredConnected==conID);
  lines(Icon)=conID;
  
  conIDnum(conID)=numel(Icon);
  
  conID=conID+1;
  
end

plines= uint8( (lines/conID)*255);


[~,conIDmax]=max(conIDnum);
Ir=find(lines==conIDmax);
Iborder=Ir;

Igr=zeros(size(Ig(:,:,1)));
Igr(Iborder)=255;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



col0=0;
colL=0;

for i=1:size(Igr,1)
   Ired=find(Igr(i,:)==255);
   if ~isempty(Ired) || (col0==1 && colL==1)
   
       if ~isempty(Ired)
           if Ired(1)==1
               col0=1;
           elseif Ired(1)==2
               col0=0;
           end

           if Ired(end)==size(Igr,2)
               colL=1;
           elseif Ired(end)==size(Igr,2)-1
               colL=0;
           end
       end
           
       if col0==1
           Ired=[1,Ired];
       end
       if colL==1
           Ired=[Ired,size(Igr,2)];
       end
       

       for j=1:numel(Ired)-1
          if Ired(j)+1 ~= Ired(j+1)
%               Ig_out(i,Ired(j)+1:Ired(j+1)-1,:)=Ig(i,Ired(j)+1:Ired(j+1)-1,:);   
              Ig_out_test(i,Ired(j):Ired(j+1))=1;  
          end
       end
       
   end    
end

col0=0;
colL=0;

for i=1:size(Igr,2)
   Ired=find(Igr(:,i)==255);
   if ~isempty(Ired) || (col0==1 && colL==1)
   
       if ~isempty(Ired)
           if Ired(1)==1
               col0=1;
           elseif Ired(1)==2
               col0=0;
           end

           if Ired(end)==size(Igr,1)
               colL=1;
           elseif Ired(end)==size(Igr,1)-1
               colL=0;
           end
       end
           
       if col0==1
           Ired=[1;Ired];
       end
       if colL==1
           Ired=[Ired;size(Igr,1)];
       end
          
       for j=1:numel(Ired)-1
          if Ired(j)+1 ~= Ired(j+1)
%               Ig_out(Ired(j)+1:Ired(j+1)-1,i,:)=Ig(Ired(j)+1:Ired(j+1)-1,i,:);   
              Ig_out_test(Ired(j):Ired(j+1),i)=Ig_out_test(Ired(j):Ired(j+1),i)+1;
          end
       end
       
   end    
end





Itest=uint32(find(Ig_out_test>=2));
IgR=Ig(:,:,1); IgG=Ig(:,:,2); IgB=Ig(:,:,3);
tempR=Ig_out(:,:,1); tempR(Itest)=IgR(Itest);
tempG=Ig_out(:,:,2); tempG(Itest)=IgG(Itest);
tempB=Ig_out(:,:,3); tempB(Itest)=IgB(Itest);

% find the value closest to white and shift it to white
% [~,Imax]=max(max( double(tempR(Itest)).^2 + double(tempG(Itest)).^2 + double(tempB(Itest)).^2 ));
% tempR(Itest)=(255-tempR(Itest(Imax)))+tempR(Itest);
% tempG(Itest)=(255-tempG(Itest(Imax)))+tempG(Itest);
% tempB(Itest)=(255-tempB(Itest(Imax)))+tempB(Itest);



Ig_out(:,:,1)=tempR; Ig_out(:,:,2)=tempG; Ig_out(:,:,3)=tempB;
Ig_out=uint8(Ig_out);

Ig_outHSV=rgb2hsv(Ig_out);

InotWhite=uint32(find(Ig_outHSV(:,:,2)>0.35));

% tempR(InotWhite)=IgR(InotWhite);
% tempG(InotWhite)=IgG(InotWhite);
% tempB(InotWhite)=IgB(InotWhite);

% IWhite=find(Ig_outHSV(:,:,2)<=0.35);
% tempR(IWhite)=255;
% tempG(IWhite)=255;
% tempB(IWhite)=255;

% Ig_out(:,:,1)=tempR; Ig_out(:,:,2)=tempG; Ig_out(:,:,3)=tempB;
% Ig_out=uint8(Ig_out);



% grayImage = rgb2gray(Ig_out); % Convert to gray so we can get the mean luminance.
% % Extract the individual red, green, and blue color channels.
% redChannel = Ig_out(:, :, 1);
% greenChannel = Ig_out(:, :, 2);
% blueChannel = Ig_out(:, :, 3);
% meanR = mean2(redChannel(Itest));
% meanG = mean2(greenChannel(Itest));
% meanB = mean2(blueChannel(Itest));
% meanGray = mean2(grayImage(Itest));
% % Make all channels have the same mean
% redChannel(Itest) = uint8(double(redChannel(Itest)) * meanGray / meanR);
% greenChannel(Itest) = uint8(double(greenChannel(Itest)) * meanGray / meanG);
% blueChannel(Itest) = uint8(double(blueChannel(Itest)) * meanGray / meanB);
% % Recombine separate color channels into a single, true color RGB image.
% Ig_out = cat(3, redChannel, greenChannel, blueChannel);


% % 
% figure(2)
% imshow(Ig_out);
% 
% figure(3)
% imshow(Igr)
% axis equal

