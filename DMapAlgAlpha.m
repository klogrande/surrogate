%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DMap Algorithm with alpha Adjustor%
%Author: Kevin LoGrande%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Data Sets%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Enter your input data here. The code opertates on data named "normData"
%The below readin files are data I have worked with, which you can delete
%
%{\
%Combustion Problem
parsedData=dlmread('ParsedNorm/2.0e_2.2dens_parse20.20.5.txt');
%
%get rid of some trajectories
[m,~]=size(parsedData);
lastPt=25*500;
parsedData(lastPt+1:m,:)=[];
%}
%
%Compute next data for old parse
nextData=parsedData;
[m,~]=size(parsedData);
for i=m:-1:1
    if mod(i,25)==0
        parsedData(i,:)=[];
    elseif mod(i,25)==1
        nextData(i,:)=[];
    end
end
%}
[m,n]=size(parsedData);
normData=zeros(m,n);
maxVal=max(parsedData);
multiplier=1./maxVal;
%
%
for i=1:m
    point=parsedData(i,:);
    shift=multiplier.*point;
    normData(i,:)=shift;
end
%}

%{
%Maryam DFT Data
normData=readmatrix('wfdata100000.csv','Range','B2:G10001');
[m,n]=size(normData);
%}

%{
%Sapta Defect Data
%normData=readmatrix('SaptaDefects/minushalf3D_DefectLayer.txt');
%normData=readmatrix('SaptaDefects/nanoparticle_inclusion.txt');
%normData=readmatrix('SaptaDefects/minushalf3d_coordinates.txt');
%normData=readmatrix('SaptaDefects/send_mp_2.txt');
normData=readmatrix('SaptaDefects/data.txt');
[m,n]=size(normData);
normData(:,3)=[];
%}

%{
%Cylinder Problem
normData=[xcyl,ycyl,zcyl];
[m,n]=size(normData);
%{
for i=m:-1:1
    if mod(i,2)==0
        normData(i,:)=[];
    end
end
 %}
%}

%{
%Swiss Roll Problem
normData=[xroll,yroll,zroll];
[m,n]=size(normData);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%Distance Matrix/Epsilon Choice%%%%%%%%%%%%%%%%%%%%%%
%Compute matrix of distance squared
dist2=zeros(m);
%{
%for PBCs
cellRad=18.01;
cellRad2=cellRad^2;
%}
for i=1:m
    for j=(i+1):m
        diff2=(normData(i,:)-normData(j,:)).^2;
        %{
        %periodic boundary conditions
        for k=1:n
            if diff2(k)>cellRad2
                diff2(k)=(abs(normData(i,k)-normData(j,k))-cellRad)^2;
            end
        end
        %}
        dist2(i,j)=sum(diff2);
        dist2(j,i)=dist2(i,j);
    end
end

%Calculating recommended epsilon value
dist4eps=dist2;
for i=1:m
    dist4eps(i,i)=NaN;
end
epsRec=sqrt(max(min(dist4eps)));
clear dist4eps %delete the matrix, since we no longer need it
fprintf('The recommended value for epsilon is %f\n', epsRec)
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Filter Points%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%An optional step for creating a more uniformly dense point cloud
%{
thresh=0.001;
track=0;
for i=1:m
    for j=i+1:m
        if dist2(i,j)<thresh && ismember(j,track)==logical(0) && ismember(i,track)==logical(0)
            track=[track;j];
        end
    end
end
track(1)=[];
del=sort(track);
normData(del,:)=[];
dist2(:,del)=[];
dist2(del,:)=[];
[m,n]=size(normData);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%Diffusion Map Algorithm%%%%%%%%%%%%%%%%%%%%%%%%%
eps=input('Input a value for epsilon.\n');
al=input('Input a value between 0 and 1 for alpha.\n');

%Create affinity matrix, diagonal matrix, and Markov matrix
W=exp(-dist2/(eps^2));
sumW=sum(W,2);
invD=spdiags((1./sumW),0,m,m);
if al==1
    L=(invD)*W*(invD);
    sumL=sum(L,2);
    DalInv=spdiags((1./sumL),0,m,m);
    K=DalInv*L;
elseif al==0
    K=invD*W;
else
    L=(invD^al)*W*(invD^al);
    sumL=sum(L,2);
    DalInv=spdiags((1./sumL),0,m,m);
    K=DalInv*L;
end

%Compute eigenvalues and eigenvectors of Markov matrix
[U,V]=eigs(K,11);
eigens=diag(V);
V=spdiags(eigens,0,11,11);
%index=[0:length(eigens)-1];
%scatter(index,eigens);
%hold on
%title('Scree Plot');
diffCoords=U*V;
diffCoords(:,1)=[];
%}
