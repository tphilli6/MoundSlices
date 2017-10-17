function [Im_out] = minecraft_reduce(Im,Ix,Iy)
%Minecraft 

[nx,ny]=size(Im);
if  size(Im,3)~=1
   error('Error! Image size not compatible with step size!'); 
end

Nx=numel(Ix)-1;
Ny=numel(Iy)-1;

% Im_out=sparse(zeros(Nx,Ny));
Im_out=(zeros(Nx,Ny));


maxPixel=max(max(Im));
minPixel=min(min(Im));
threshhold=0.5*(maxPixel-minPixel)+minPixel;

for jj=1:Ny-1
    for ii=1:Nx-1
        IIx=Ix(ii):Ix(ii+1)-1;
        IIy=Iy(jj):Iy(jj+1)-1;
        meanPixel=mean(mean( Im(IIx,IIy) ));
        if meanPixel<=threshhold
          Im_out(ii,jj) = 1;
        end
    end
end

% imshow(Im_out);