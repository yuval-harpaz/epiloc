function pairs=getErrors(distances)

distSum=0;
depthErr=0;
pairs=false(size(distances));
for pairi=1:min(size(distances))
    [minD,minDi]=sort(distances(:));
    x=find(sum(distances==minD(1),1));
    x=x(1); % when, say, there are two zero errors
    y=find(distances(:,x)==minD(1));
    y=y(1);
    %y=find(sum(distances==minD(1),2));
    %distSum=distSum+distances(y,x);
    %if size(distances,1)>size(distances,2) % more rows
    distances(:,x)=max(minD(end))+1;
    %else
    distances(y,:)=max(minD(end))+1;
    %end
    pairs(y,x)=true;
end