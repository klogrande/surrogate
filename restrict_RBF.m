%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Restriction through Radial Basis Functions
%Author: Kevin LoGrande
%Description: Takes discreet data from a low-dimensional manifold and its
             %higher-dimensional naturally embedded representation and uses 
             %Radial Basis Functions to create a high-to-low dimension 
             %interpolation for new data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newCoord=restrict_RBF(lowDimData,highDimData,newPoint,nn,p)
[numPoints,numHighDim]=size(highDimData);
numLowDim=min(size(lowDimData));
newCoord=zeros(1,numLowDim);
dist=zeros(numPoints,1);
Lambda=zeros(nn);
uvec=zeros(nn,1);

for a=1:numLowDim
    for i=1:numPoints
        %calculate squared distance between new point and existing data in
        %the high-dimensional space
        sqDist=0;
        for j=1:numHighDim
            sqDist=sqDist+(newPoint(j)-highDimData(i,j))^2;
        end
        %put square distances in a vector
        dist(i)=sqDist;
    end
    [dist,index]=sort(dist);
    dist=dist(1:nn);
    dist=dist.^(p/2);
    index=index(1:nn);
    
    %calculating distance matrix for  nearest neighbors in high-dim space
    for i=1:nn
        for j=i+1:nn
            sqDist=0;
            for k=1:numHighDim
                sqDist=sqDist+(highDimData(index(i),k)-highDimData(index(j),k))^2;
            end
            Lambda(i,j)=sqDist^(p/2);
            Lambda(j,i)=Lambda(i,j);
        end
    end
    
    %reading in low-dim coordinates
    for i=1:nn
        uvec(i)=lowDimData(index(i),a);
    end
    aCoeffs=Lambda\uvec;
    
    for i=1:nn
        newCoord(a)=newCoord(a)+aCoeffs(i)*dist(i);
    end
end