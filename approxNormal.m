%create matrix of distance squared between each point of low dim data
[m,n]=size(lowData);
dist2low=zeros(m);
for i=1:m
    for j=(i+1):m
        diff2low=(lowData(i,:)-lowData(j,:)).^2;
        dist2low(i,j)=sum(diff2low);
        dist2low(j,i)=dist2low(i,j);
    end
end
lowData_red=lowData;
eps2=0.1;

%{
%filter points that are within a certain threshold
thresh=5e-8;
track=0;
for i=1:m
    for j=i+1:m
        if dist2low(i,j)<thresh && ismember(j,track)==logical(0) && ismember(i,track)==logical(0)
            track2=[track;j];
            clear track
            track=track2;
            clear track2
        end
    end
end
track(1)=[];
del=sort(track);
dist2low(:,del)=[];
dist2low(del,:)=[];
lowData_red(del,:)=[];
%}

[p,q]=size(dist2low);
indicator=zeros(p,q);
for i=1:p
    for j=i+1:p
        if dist2low(i,j)<eps2
            indicator(i,j)=1;
            indicator(j,i)=1;
        end
    end
end
rowSums=sum(indicator,2);
for i=1:p
    indicator(i,i)=-rowSums(i);
end
indicator=sparse(indicator);
rowDiag=spdiags(rowSums);
nVec=indicator*lowData_red;
modN=sqrt(sum(nVec.*nVec,2));