%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Lifting through Radial Basis Functions
%Author: Kevin LoGrande
%Description: Takes discreet data from a high-dimensional manifold and its
             %lower-dimensional Diffusion Map representation and uses 
             %Radial Basis Functions to create a low-to-high dimension 
             %interpolation for new data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newCoord = lift(highDimData,lowDimData,newPoint,nn,p)
[numPoints,numHighDim]=size(highDimData);
numLowDim=min(size(lowDimData));
newCoord=zeros(1,numHighDim);
dist=zeros(numPoints,1);
Lambda=zeros(nn);
yvec=zeros(nn,1);

for a=1:numHighDim
    for i=1:numPoints
        %calculate squared distance between new point and existing data
        sqDist=0;
        for j=1:numLowDim
            sqDist=sqDist+(newPoint(j)-lowDimData(i,j))^2;
        end
        %put square distances in a vector
        dist(i)=sqDist;
    end
    [dist,index]=sort(dist);
    dist=dist(1:nn);
    dist=dist.^(p/2);
    index=index(1:nn);
    
    %calculating distance matrix for  nearest neighbors
    for i=1:nn
        for j=i+1:nn
            sqDist=0;
            for k=1:numLowDim
                sqDist=sqDist+(lowDimData(index(i),k)-lowDimData(index(j),k))^2;
            end
            Lambda(i,j)=sqDist^(p/2);
            Lambda(j,i)=Lambda(i,j);
        end
    end
    
    for i=1:nn
        yvec(i)=highDimData(index(i),a);
    end
    aCoeffs=Lambda\yvec;
    
    for i=1:nn
        newCoord(a)=newCoord(a)+aCoeffs(i)*dist(i);
    end
end