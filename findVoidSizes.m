function Void = findVoidSizes(Imound,Ivoid,Iborder,Size)

A=255*ones(Size);
A(Imound)=0;
A(Iborder)=0;
A(Ivoid)=1;

Irx=repmat([1:Size(1)]',[1,Size(2)]);
Iry=repmat([1:Size(2)], [Size(1),1]);

% Trace the connected path and exclude non-connected sections %%%%%%%%%%%%%
isvoid=uint32(find(A==1));
isConnected=zeros(size(A));
conID=-1;
xy2i=reshape( 1:numel(Irx),size(Irx));
% lines=zeros(size(Igr));

nVoids=0;
while any(isConnected(isvoid)==0)
  i=find(isConnected(isvoid)==0,1,'first');
  
  queue=isvoid(i);
 
  while ~isempty(queue)
      
      queueNew=[];
      for jj=1:length(queue)
          pixel=queue(jj);
          isConnected(pixel) = conID; %assign connection ID

          %pull neighbors
          iy=Iry(pixel);
          ix=Irx(pixel); %<===== changing this resulted in an infinite loop
          il=max(ix-1,1); ih=min(ix+1,Size(1));
          jl=max(iy-1,1); jh=min(iy+1,Size(2));
          subArray=A( il:ih, jl:jh );
          xy2iSA=xy2i( il:ih, jl:jh);
          
          Ineigh=find(subArray==1);
          Ineigh5=find(Ineigh~=5 & isConnected(xy2iSA(Ineigh))==0  );
          Ineigh=Ineigh(Ineigh5);

          
          %Store connection ids for neighbors
          isConnected(xy2iSA(Ineigh)) = conID;
          queueNew=[queueNew;xy2iSA(Ineigh)];

          
          %remove current pixel from queue
%           Iremove=find(queue~=pixel);
%           queue=queue(Iremove);



      end
      queue=queueNew;
      

  end
  

  
  Icon=find(isConnected==conID);
  
  nVoids=nVoids+1;
  ID(nVoids)=conID;
  voidPixels{nVoids}=Icon;
  voidSize(nVoids)=numel(Icon);
  
%   lines(Icon)=conID;
  
%   conIDnum(conID)=numel(Icon);
  
  conID=conID-1;
  
end



[~,Isort]=sort(voidSize,'descend');

label=1;
for i=Isort
    Void.ID(label)=label;
    Void.size(label)=voidSize(i);
    Void.pixels{label}=voidPixels{i};
    
    A(voidPixels{label})=label;
    
    label=label+1;
end
% [minx] = min(Irx(Iborder));
% [maxx] = max(Irx(Iborder));
% [miny] = min(Iry(Iborder));
% [maxy] = max(Iry(Iborder));
% I255=find(A==255);
% % A(I255)=label;
% Argb(:,:,1)=A;
% Argb(:,:,2)=A;
% Argb(:,:,3)=A;
% % A(Imound)=label;
% A=A/label*255;
