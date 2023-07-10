%Surrogate Model Function

function [traj,time] = surrogate(highData,lowData,dudt,newPoint,stepSize,numSteps)

%initializing variables
[~,lowDim]=size(lowData);
traj=zeros(numSteps+1,lowDim);
time=zeros(numSteps+1,1);

%restrict first point to the manifold
newPt_red=restrict_RBF(lowData,highData,newPoint,50,3);
traj(1,:)=newPt_red;

%calculate the rest of the trajectory by stepping along vector field
for i=2:numSteps+1
    %vector for new point approximated using radial basis functions
    vec=lift(dudt,lowData,traj(i-1,:),50,3);
    speed=norm(vec);
    
    %next step in trajectory and keeping track of time for dynamics
    traj(i,:)=traj(i-1,:)+stepSize*(vec/speed);
    time(i)=time(i-1)+stepSize/speed;
end
