
clear all
clc
close all
%%
%Parameters of the simulation

Pxsize=0.7; %in nm/px
CellVolume=0.197; %in nm^3
NinCell=4; %for Cd atoms in CdS cubic cell
NumberOfAtomsPerCubicPixel=round(NinCell*(Pxsize^3)/CellVolume);

MaxParticleRadius=round(15*Pxsize); %in nm
MinParticleRadius=round(3*Pxsize); %in nm
AverageParticleRadius=round(9*Pxsize); %in nm
STD=3.8*Pxsize; %in nm

Dopant2Matrix=1/50; %ratio between dopant to matrix atoms (in fraction)

m = 400; n = 400; p=2*MaxParticleRadius+1; %Final size of the array 

%%
%initial parameters 
BulkMap = zeros(m,n,p);
SurfaceMap = zeros(m,n,p);
SumOfRadiiX=0;
SumOfRadiiY=0;
%%

while 1
        radius=round(STD*randn(1)+AverageParticleRadius);
        if ~(radius<=MaxParticleRadius && radius>=MinParticleRadius)
            continue
        elseif SumOfRadiiY>=m
            break
        end
 while 1
        radius=round(STD*randn(1)+AverageParticleRadius); %radius of the sphere
        if ~(radius<=MaxParticleRadius && radius>=MinParticleRadius)
            continue
        elseif SumOfRadiiX>=n
            break
        end
        Center=[radius+SumOfRadiiX+1 radius+SumOfRadiiY+1 MaxParticleRadius+1]; %calculate center
        [px,py,pz] = meshgrid(1:n, 1:m, 1:p);
        xc = Center(1); yc = Center(2); zc = Center(3); % the center of sphere
        
        BulkSphere = (px-xc).^2 + (py-yc).^2 + (pz-zc).^2 <=radius*radius;% create bulk map
        BulkMap(BulkSphere) = 1;

        SurfaceSphere=edge3(BulkSphere,'approxcanny',0.6);
        SurfaceMap(SurfaceSphere)=1;
       
        SumOfRadiiX=SumOfRadiiX+2*radius+3;
 end
  SumOfRadiiX=0;
  SumOfRadiiY=SumOfRadiiY+2*MaxParticleRadius+3;
end

%%

DopantMatrixBulk=zeros(size(BulkMap));%create the marix map 
MatrixAtomsLocation=find(BulkMap);%locate the location of each matrix atom


%dopant in Bulk

for i=1:size(MatrixAtomsLocation)
    for k=1:NumberOfAtomsPerCubicPixel
        RandomNumber=rand;
        if RandomNumber<=Dopant2Matrix
        DopantMatrixBulk(MatrixAtomsLocation(i))=1+DopantMatrixBulk(MatrixAtomsLocation(i));
        end
    end
end

%dopant in surface
DopantMatrixSurface=zeros(size(SurfaceMap));%create the dopant map
SurfaceAtomsLocation=find(SurfaceMap);%locate the location of each matrix atom on the surface

NumberOfMatrixAtoms=sum(BulkMap(:));
NumberOfSurfaceAtoms=size(SurfaceAtomsLocation,1);
MultiplicityFactor=NumberOfMatrixAtoms/NumberOfSurfaceAtoms;

for i=1:size(SurfaceAtomsLocation)
    RandomNumber=rand;
    for k=1:NumberOfAtomsPerCubicPixel
        if RandomNumber<=Dopant2Matrix*  MultiplicityFactor
         DopantMatrixSurface(SurfaceAtomsLocation(i))=1+DopantMatrixSurface(SurfaceAtomsLocation(i));
        end
    end
end
%%
%calculate Z for bulk and surface simulation
ZprojBulk=0;
ZprojSurf=0;
ZprojDopantBulk=0;
ZprojDopantSurface=0;


for z=1:p
   ZprojBulk=ZprojBulk+BulkMap(:,:,z);%Bulk simulation
   ZprojSurf=ZprojSurf+SurfaceMap(:,:,z);%Surface simulation

   ZprojDopantBulk=ZprojDopantBulk+DopantMatrixBulk(:,:,z);%
   ZprojDopantSurface=ZprojDopantSurface+DopantMatrixSurface(:,:,z);
end
%%

% Plot everything
Bulk=im2uint8(mat2gray(ZprojBulk,[1,100]));
Surface=im2uint8(mat2gray(ZprojSurf,[1,100]));
DopantBulk=im2uint8(mat2gray(ZprojDopantBulk,[1,100]));
Dopantsurface=im2uint8(mat2gray(ZprojDopantSurface,[1,100]));


% Bulk=ZprojBulk;
% Surface=ZprojSurf;
% DopantBulk=ZprojDopantBulk;
% Dopantsurface=ZprojDopantSurface;


%%
figure(1); imshow(Bulk);
figure (2); imshow(Surface);
figure (3); imshow(DopantBulk);
figure (4); imshow(Dopantsurface);






