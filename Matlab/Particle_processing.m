clear all
close all
clc
%%
%Plug in pixel size for radius calculation (in nm), file name and particle
%size (in px) for ploting
% xlsFileName='paritcal (simulation).xlsx'
% xlsFileName='paritcal 1_50_bulk.xlsx';
% xlsFileName='paritcal 1_50_surface.xlsx';
% xlsFileName='paritcal 1_1000 bulk.xlsx';
%xlsFileName='paritcal 1_1000 surface.xlsx';
xlsFileName='paritcalAll.xlsx';

PxSize=0.7;
ParticleSize=7;
%% define matrices
ExcelFile=xlsread(xlsFileName,'','','basic');

NumOfParticles=max(ExcelFile(:,1));
MaxRadius=max(ExcelFile(:,2));

CdCount=zeros(MaxRadius+3,MaxRadius+1);
CdIntensity=zeros(MaxRadius+3,MaxRadius+1);
MnCount=zeros(MaxRadius+3,MaxRadius+1);
MnIntensity=zeros(MaxRadius+3,MaxRadius+1);

ErrorMatrix=zeros(1,MaxRadius+1);

CdCountErrorMatrix=zeros(1,MaxRadius+2);
CdIntensityErrorMatrix=zeros(1,MaxRadius+2);
MnCountErrorMatrix=zeros(1,MaxRadius+2);
MnIntensityErrorMatrix=zeros(1,MaxRadius+2);

CdCountErrorMatrixSTD=zeros(MaxRadius,MaxRadius+1);
CdIntensityErrorMatrixSTD=zeros(MaxRadius,MaxRadius+1);
MnCountErrorMatrixSTD=zeros(MaxRadius,MaxRadius+1);
MnIntensityErrorMatrixSTD=zeros(MaxRadius,MaxRadius+1);

CdCountErrorMatrixMean=zeros(MaxRadius,MaxRadius+1);
CdIntensityErrorMatrixMean=zeros(MaxRadius,MaxRadius+1);
MnCountErrorMatrixMean=zeros(MaxRadius,MaxRadius+1);
MnIntensityErrorMatrixMean=zeros(MaxRadius,MaxRadius+1);

MaxNumOfDopants=300;
MnPoissonMatrix=zeros(MaxNumOfDopants,MaxRadius+3);
for i=0:MaxNumOfDopants
    MnPoissonMatrix(i+1,1)=i;
end


iParticle=zeros(MaxRadius+1,5);

for n=1:MaxRadius+1
    CdCount(1,n+1)=n;
    CdIntensity(1,n+1)=n;
    MnCount(1,n+1)=n;
    MnIntensity(1,n+1)=n;
    CdCountErrorMatrixSTD(1,n+1)=n;
    CdIntensityErrorMatrixSTD(1,n+1)=n;
    MnCountErrorMatrixSTD(1,n+1)=n;
    MnIntensityErrorMatrixSTD(1,n+1)=n;
    CdCountErrorMatrixMean(1,n+1)=n;
    CdIntensityErrorMatrixMean(1,n+1)=n;
    MnCountErrorMatrixMean(1,n+1)=n;
    MnIntensityErrorMatrixMean(1,n+1)=n;
end


%% analysis for Cd Signal -  Counts
%Average over particle size
for n=1:size(ExcelFile,1) 
    z=ExcelFile(n,2);
        if rem(ExcelFile(n,2),MaxRadius+1)==0 && n~=1 ;
        
        iParticleSize=nnz(iParticle(1:end,2)); %find size of particle
        [row,col]=find(ismember(CdCount(1,1:end),iParticleSize-1));%find where to add this particle
        
        %Add data to talbs
        CdCount(3:end,col)=CdCount(3:end,col)+iParticle(1:end,2);
        CdIntensity(3:end,col)=CdIntensity(3:end,col)+iParticle(1:end,3);
        MnCount(3:end,col)=MnCount(3:end,col)+iParticle(1:end,4);
        MnIntensity(3:end,col)=MnIntensity(3:end,col)+iParticle(1:end,5);

        % count particles
        CdCount(2,col)=CdCount(2,col)+1; 
        CdIntensity(2,col)=CdIntensity(2,col)+1;
        MnCount(2,col)=MnCount(2,col)+1;
        MnIntensity(2,col)=MnIntensity(2,col)+1;
        
        %add data for error calculation for particles with the same size
        TiParticle=transpose(iParticle);
        
        CdCountErrorMatrix(end+1,1:MaxRadius+1)=TiParticle(2,1:end);
        CdCountErrorMatrix(end,end)=iParticleSize;

        CdIntensityErrorMatrix(end+1,1:MaxRadius+1)=TiParticle(3,1:end);
        CdIntensityErrorMatrix(end,end)=iParticleSize;

        MnCountErrorMatrix(end+1,1:MaxRadius+1)=TiParticle(4,1:end);
        MnCountErrorMatrix(end,end)=iParticleSize;

        MnIntensityErrorMatrix(end+1,1:MaxRadius+1)=TiParticle(5,1:end);
        MnIntensityErrorMatrix(end,end)=iParticleSize;

        %for possion statistic extraction Here I count how many Mn in each
        %particles
        NumberOfDopants=sum(iParticle(:,4));
        MnPoissonMatrix(NumberOfDopants+1,iParticleSize)=MnPoissonMatrix(NumberOfDopants+1,iParticleSize)+1;
        

        %reset particle idx
        iParticle=zeros(MaxRadius+1,5);
        end

    iDistanceFromEdge=ExcelFile(n,2);
    iParticle(iDistanceFromEdge+1,1)=iDistanceFromEdge;
    iParticle(iDistanceFromEdge+1,2)=ExcelFile(n,3);%Cd Count
    iParticle(iDistanceFromEdge+1,3)=ExcelFile(n,4);%Cd Intensity
    iParticle(iDistanceFromEdge+1,4)=ExcelFile(n,5);%Mn Count
    iParticle(iDistanceFromEdge+1,5)=ExcelFile(n,6);%Mn Intensity

end


for k=1:size(MnPoissonMatrix,1)
     MnPoissonMatrix(k,size(MnPoissonMatrix,2))=sum(MnPoissonMatrix(k,2:size(MnPoissonMatrix,2)));
end
xlswrite('PoissonDist.xlsx',MnPoissonMatrix)

%%
        CdCountErrorMatrix=sortrows(CdCountErrorMatrix,size(CdCountErrorMatrix,2));
        CdIntensityErrorMatrix=sortrows(CdIntensityErrorMatrix,size(CdIntensityErrorMatrix,2));
        MnCountErrorMatrix=sortrows(MnCountErrorMatrix,size(MnCountErrorMatrix,2));
        MnIntensityErrorMatrix=sortrows(MnIntensityErrorMatrix,size(MnIntensityErrorMatrix,2));
        
%%
for n=1:MaxRadius+1%for each particle size
    [row,col]=find(CdCountErrorMatrix(:,MaxRadius+2)==n);
    if ~isempty(row) 
        for q=1:MaxRadius+1 %for each radius in particle
        M=mean(CdCountErrorMatrix(row,q));
        STDV=std(CdCountErrorMatrix(row,q));
        [row2,col2]=find(ismember(CdCountErrorMatrixMean(1,1:end),q));
        CdCountErrorMatrixMean(n,col2)=M;
        CdCountErrorMatrixSTD(n,col2)=STDV/sqrt(numel(CdCountErrorMatrix(row,q)));
        end
    end
end

for n=1:MaxRadius+1%for each particle size
    [row,col]=find(CdIntensityErrorMatrix(:,MaxRadius+2)==n);
    if ~isempty(row) 
        for q=1:MaxRadius+1 %for each radius in particle
        M=mean(CdIntensityErrorMatrix(row,q));
        STDV=std(CdIntensityErrorMatrix(row,q));
        [row2,col2]=find(ismember(CdIntensityErrorMatrixMean(1,1:end),q));
        CdIntensityErrorMatrixMean(n,col2)=M;
        CdIntensityErrorMatrixSTD(n,col2)=STDV/sqrt(numel(CdIntensityErrorMatrix(row,q)));
        end
    end
end
for n=1:MaxRadius+1%for each particle size
    [row,col]=find(MnCountErrorMatrix(:,MaxRadius+2)==n);
    if ~isempty(row) 
        for q=1:MaxRadius+1 %for each radius in particle
        M=mean(MnCountErrorMatrix(row,q));
        STDV=std(MnCountErrorMatrix(row,q));
        [row2,col2]=find(ismember(MnCountErrorMatrixMean(1,1:end),q));
        MnCountErrorMatrixMean(n,col2)=M;
        MnCountErrorMatrixSTD(n,col2)=STDV/sqrt(numel(MnCountErrorMatrix(row,q)));
        end
    end
end
for n=1:MaxRadius+1%for each particle size
    [row,col]=find(MnIntensityErrorMatrix(:,MaxRadius+2)==n);
    if ~isempty(row) 
        for q=1:MaxRadius+1 %for each radius in particle
        M=mean(MnIntensityErrorMatrix(row,q));%here i cal the average of Intensity/counts for a given distance of a given particle size
        STDV=std(MnIntensityErrorMatrix(row,q));
        [row2,col2]=find(ismember(MnIntensityErrorMatrixMean(1,1:end),q));
        MnIntensityErrorMatrixMean(n,col2)=M;
        MnIntensityErrorMatrixSTD(n,col2)=STDV/sqrt(numel(MnIntensityErrorMatrix(row,q)));
        end
    end
end


n=1;
CdCountErrorMatrixSTD(:,1+n:end+n)=CdCountErrorMatrixSTD(:,1:end);
CdCountErrorMatrixSTD(:,1:n)=0;
CdCountErrorMatrixSTD=CdCountErrorMatrixSTD.';
CdCountErrorMatrixSTD(1:2,2:end)=CdCount(1:2,2:end-1);

CdIntensityErrorMatrixSTD(:,1+n:end+n)=CdIntensityErrorMatrixSTD(:,1:end);
CdIntensityErrorMatrixSTD(:,1:n)=0;
CdIntensityErrorMatrixSTD=CdIntensityErrorMatrixSTD.';
CdIntensityErrorMatrixSTD(1:2,2:end)=CdCount(1:2,2:end-1);

MnCountErrorMatrixSTD(:,1+n:end+n)=MnCountErrorMatrixSTD(:,1:end);
MnCountErrorMatrixSTD(:,1:n)=0;
MnCountErrorMatrixSTD=MnCountErrorMatrixSTD.';
MnCountErrorMatrixSTD(1:2,2:end)=MnCount(1:2,2:end-1);

MnIntensityErrorMatrixSTD(:,1+n:end+n)=MnIntensityErrorMatrixSTD(:,1:end);
MnIntensityErrorMatrixSTD(:,1:n)=0;
MnIntensityErrorMatrixSTD=MnIntensityErrorMatrixSTD.';
MnIntensityErrorMatrixSTD(1:2,2:end)=MnCount(1:2,2:end-1);

CdCountErrorMatrixMean(:,1+n:end+n)=CdCountErrorMatrixMean(:,1:end);
CdCountErrorMatrixMean(:,1:n)=0;
CdCountErrorMatrixMean=CdCountErrorMatrixMean.';
CdCountErrorMatrixMean(1:2,2:end)=CdCount(1:2,2:end-1);

CdIntensityErrorMatrixMean(:,1+n:end+n)=CdIntensityErrorMatrixMean(:,1:end);
CdIntensityErrorMatrixMean(:,1:n)=0;
CdIntensityErrorMatrixMean=CdIntensityErrorMatrixMean.';
CdIntensityErrorMatrixMean(1:2,2:end)=CdCount(1:2,2:end-1);

MnCountErrorMatrixMean(:,1+n:end+n)=MnCountErrorMatrixMean(:,1:end);
MnCountErrorMatrixMean(:,1:n)=0;
MnCountErrorMatrixMean=MnCountErrorMatrixMean.';
MnCountErrorMatrixMean(1:2,2:end)=MnCount(1:2,2:end-1);

MnIntensityErrorMatrixMean(:,1+n:end+n)=MnIntensityErrorMatrixMean(:,1:end);
MnIntensityErrorMatrixMean(:,1:n)=0;
MnIntensityErrorMatrixMean=MnIntensityErrorMatrixMean.';
MnIntensityErrorMatrixMean(1:2,2:end)=MnCount(1:2,2:end-1);

%%
%Plot everything
close all
figure(1)
hold on

x=[CdCount(1,1:end-1)]*PxSize;

title('CdCount')
y=[CdCountErrorMatrixMean(3:end,ParticleSize+1)];
error=CdCountErrorMatrixSTD(3:end,ParticleSize+1);
y=transpose(y);
errorbar(x,y,error,'m','linewidth',2)
hold off


figure(2)
hold on
title('CdIntensity')
error=CdIntensityErrorMatrixSTD(3:end,ParticleSize+1);
y=[CdIntensityErrorMatrixMean(3:end,ParticleSize+1)];
y=transpose(y);
errorbar(x,y,error,'m','linewidth',2)
hold off

figure(3)
title('MnCount')
hold on
error=MnCountErrorMatrixSTD(3:end,ParticleSize+1);
y=[MnCountErrorMatrixMean(3:end,ParticleSize+1)];
y=transpose(y);
errorbar(x,y,error,'m','linewidth',2)
hold off

figure(4)
title('MnIntensity')
hold on
error=MnIntensityErrorMatrixSTD(3:end,ParticleSize+1);
y=[MnIntensityErrorMatrixMean(3:end,ParticleSize+1)];
y=transpose(y);
errorbar(x,y,error,'m','linewidth',2)
hold off

%%
%MLLS fitting

%X contain referenc espectrums, Y is the target spectrum, A contain the
%coeeficents of each refrence spectrum
%Y- Mn count signal
%X = x1-Cd count (surface contribution). x2-Cd intensity (bulk contribution)

x1=[CdCountErrorMatrixMean(3:end,ParticleSize+1)]/max(CdCountErrorMatrixMean(3:end,ParticleSize+1));
x2=[CdIntensityErrorMatrixMean(3:end,ParticleSize+1)]/max(CdIntensityErrorMatrixMean(3:end,ParticleSize+1));
Y=MnCountErrorMatrixMean(3:end,ParticleSize+1)/max(MnCountErrorMatrixMean(3:end,ParticleSize+1));
%


X=[ones(size(x1)) x1 x2];

A=inv(X'*X)*X'*Y;
SSE = (Y-X*A)'*(Y-X*A);
SStot = sum((Y-mean(Y)).^2);

R2 = 1-SSE/SStot;
sum(A);

MLLSretuls = ['MllS fitting results for particle size of ',num2str(ParticleSize),':',num2str(A(2)),' for surface contribution and ', num2str(A(3)),' for bulk contribution' ,' with R2 of ',num2str(R2)];
disp(MLLSretuls)

figure(5)
plot(Y,X*A, '.')
hold on
plot([min(Y),max(Y)],[min(Y),max(Y)])

figure(6)
plot(x,X*A,'m','linewidth',2)
hold on
errorMnCount=MnCountErrorMatrixSTD(3:end,ParticleSize+1);
plot(x,Y,'.')

%%
X=linspace(0,MaxRadius+1,MaxRadius+1)/(MaxRadius+1);
%consoladate all data into one particle
for i=2:MaxRadius+1%cols
    v1=nonzeros(CdCountErrorMatrixMean(3:end,i));
    v1=v1/max(v1);
    x1=nonzeros(CdCountErrorMatrixMean(3:2+size(v1,1),1))-1;
    x1=x1/max(x1);
    
    v2=nonzeros(CdIntensityErrorMatrixMean(3:end,i));
    v2=v2/max(v2);
    x2=nonzeros(CdIntensityErrorMatrixMean(3:2+size(v2,1),1))-1;
    x2=x2/max(x2);
    
    v3=nonzeros(MnCountErrorMatrixMean(3:end,i));
    v3=v3/max(v3);
    x3=nonzeros(MnCountErrorMatrixMean(3:2+size(v3,1),1))-1;
    x3=x3/max(x3);
    
    v4=nonzeros(MnIntensityErrorMatrixMean(3:end,i));
    v4=v4/max(v4);
    x4=nonzeros(MnIntensityErrorMatrixMean(3:2+size(v4,1),1))-1;
    x4=x4/max(x4);

if (~isempty(x1))&&(~isempty(v1))%rows
    for j=3:MaxRadius+3
        if (j-3)/(MaxRadius)==1
         CdCountErrorMatrixMean(j,i)=v1(end,1);
         CdIntensityErrorMatrixMean(j,i)=v2(end,1);
         MnCountErrorMatrixMean(j,i)=v3(end,1);
         MnIntensityErrorMatrixMean(j,i)=v4(end,1);
        else
         CdCountErrorMatrixMean(j,i)=interp((j-3)/(MaxRadius),x1,v1);
         CdIntensityErrorMatrixMean(j,i)=interp((j-3)/(MaxRadius),x2,v2);
         MnCountErrorMatrixMean(j,i)=interp((j-3)/(MaxRadius),x3,v3);
         MnIntensityErrorMatrixMean(j,i)=interp((j-3)/(MaxRadius),x4,v4);
        end
    end
end
end

for i=3:MaxRadius+3
CdCountErrorMatrixMean(i,MaxRadius+2)=mean(nonzeros(CdCountErrorMatrixMean(i,2:MaxRadius+1)));
CdCountErrorMatrixMean(i,MaxRadius+3)=std(nonzeros(CdCountErrorMatrixMean(i,2:MaxRadius+1)))/sqrt(length(nonzeros(CdCountErrorMatrixMean(i,2:MaxRadius+1))));

CdIntensityErrorMatrixMean(i,MaxRadius+2)=mean(nonzeros(CdIntensityErrorMatrixMean(i,2:MaxRadius+1)));
CdIntensityErrorMatrixMean(i,MaxRadius+3)=std(nonzeros(CdIntensityErrorMatrixMean(i,2:MaxRadius+1)))/sqrt(length(nonzeros(CdIntensityErrorMatrixMean(i,2:MaxRadius+1))));

MnCountErrorMatrixMean(i,MaxRadius+2)=mean(nonzeros(MnCountErrorMatrixMean(i,2:MaxRadius+1)));
MnCountErrorMatrixMean(i,MaxRadius+3)=std(nonzeros(MnCountErrorMatrixMean(i,2:MaxRadius+1)))/sqrt(length(nonzeros(MnCountErrorMatrixMean(i,2:MaxRadius+1))));

MnIntensityErrorMatrixMean(i,MaxRadius+2)=mean(nonzeros(MnIntensityErrorMatrixMean(i,2:MaxRadius+1)));
MnIntensityErrorMatrixMean(i,MaxRadius+3)=std(nonzeros(MnIntensityErrorMatrixMean(i,2:MaxRadius+1)))/sqrt(length(nonzeros(MnIntensityErrorMatrixMean(i,2:MaxRadius+1))));
end

figure(7)
x=X;
y=CdCountErrorMatrixMean(3:end,MaxRadius+2);
error=CdCountErrorMatrixMean(3:end,MaxRadius+3);
errorbar(x,y,error,'m','linewidth',2)
title('CdCountMean-All sizes')

figure(8)
x=X;
y=CdIntensityErrorMatrixMean(3:end,MaxRadius+2);
error=CdIntensityErrorMatrixMean(3:end,MaxRadius+3);
errorbar(x,y,error,'m','linewidth',2)
title('CdIntensityMean-All sizes')

figure(9)
x=X;
y=MnCountErrorMatrixMean(3:end,MaxRadius+2);
error=MnCountErrorMatrixMean(3:end,MaxRadius+3);
errorbar(x,y,error,'m','linewidth',2)
title('MnCountMean-All sizes')

figure(10)
x=X;
y=MnIntensityErrorMatrixMean(3:end,MaxRadius+2);
error=MnIntensityErrorMatrixMean(3:end,MaxRadius+3);
errorbar(x,y,error,'m','linewidth',2)
title('MnIntensityMean-All sizes')
%%
%MLLS fitting for 'reduced' particle

%X contain referenc espectrums, Y is the target spectrum, A contain the
%coeeficents of each refrence spectrum
%Y- Mn count signal
%X = x1-Cd count (surface contribution). x2-Cd intensity (bulk contribution)

x1 = CdCountErrorMatrixMean(3:end,MaxRadius+2);
x2=CdIntensityErrorMatrixMean(3:end,MaxRadius+2);
Y=MnCountErrorMatrixMean(3:end,MaxRadius+2);


X=[ones(size(x1)) x1 x2];

A=inv(X'*X)*X'*Y;
SSE = (Y-X*A)'*(Y-X*A);
SStot = sum((Y-mean(Y)).^2);

R2 = 1-SSE/SStot;
sum(A);

MLLSretuls = ['MllS fitting results for reduced particle (all particles) are ',num2str(A(2)),' for surface contribution and ',num2str(A(3)),' for bulk contribution' ,' with R2 of ',num2str(R2)];
disp(MLLSretuls)

figure(11)
plot(Y,X*A, '.')
hold on
plot([min(Y),max(Y)],[min(Y),max(Y)])
title('MLLS fitting-All sizes - residuals')

figure(12)
plot(x,X*A,'m','linewidth',2)
hold on
errorMnCount=MnCountErrorMatrixMean(3:end,MaxRadius+3);
plot(x,Y,'.')
title('MLLS fitting-All sizes')



%%
function ystar=interp(xstar,X,Y)

% find x0, x1,y0,y1
% x0 <= xstar <= x1
loc = find(xstar<X,1);
x1=X(loc);
x0=X(loc-1);
y1=Y(loc);
y0=Y(loc-1);

m=(y1-y0)/(x1-x0);
ystar=y0+m*(xstar-x0);
end

