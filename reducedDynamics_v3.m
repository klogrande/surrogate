%Reduced Dynamics v3
%Author: Kevin LoGrande
%Description: Uses discrete differentiation to define a discrete vector
%field on a low-dimensional manifold discovered in high dimensional data
%Before running this, complete the following:
    %1)save "parsedData" and "nextData" variables from parseData_v3
    %2)run DMapAlgAlpha with "parsedData" and DO NOT CLEAR VARIABLES
    %3)define your "lowData" coordinates and eigenvalues to proceed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%read in next point data for new parses
nextData=dlmread('ParsedNorm/2.0e_2.2dens_UNInext_500traj.txt');
nextPts=nextData(:,1:n);
timeSteps=nextData(:,n+1);
nLow=size(lowData,2);
%}
%
%read next data for old parses
nextPts=nextData;
timeSteps=0.005*ones(m,1);
nLow=size(lowData,2);
%}
%restrict next point data from ambient space to the manifold
%here restriction is done throguh radial basis functions with 50 nearest
    %neighbors and p=3
nextRed=zeros(m,nLow);
for i=1:m
    nextRed(i,:)=restrict_RBF(lowData,parsedData,nextPts(i,:),50,3);
end

%now produce the discrete vector field with simple differences
%note: dt is given in the timeSteps variable
byDt=spdiags((1./timeSteps),0,m,m);
dudt=byDt*(nextRed-lowData);