function marikVirtual29plot(input,pnt,results)

for linei=1:size(input,1)
    D(linei)=getDist(pnt,input(linei,:),input(linei,:));
end

binSum=[];
binNacc=[];
binDacc=[];
binD=[];
conds={'MED','AVG','R'};
fff=0;
for condi=1:3
    for bini=3:15
        rowi=logical((D>=(bini-10)*10).*(D<(bini*10)));
        binSum(bini-2)=sum(rowi);
        binNacc(bini-2)=mean(results(rowi,1+fff)); % ratio of number of dipoles (1 of 2 etc)
        binDacc(bini-2)=mean(results(rowi,2+fff)); % ratio of right location
        binD(bini-2)=mean(results(rowi,3+fff));
    end
    fff=fff+3;
    figure;bar([35:10:155],binNacc)
    xlabel('Distance between dipoles (mm)')
    ylabel('detected dipoles')
    title(['How many dipoles were found',' ',conds{condi}])
    figure;bar([35:10:155],binDacc)
    xlabel('Distance between dipoles (mm)')
    ylabel('correctly detected dipoles')
    title(['How many dipoles were correctly localized',' ',conds{condi}])
    figure;bar([35:10:155],binD)
    xlabel('Distance between dipoles (mm)')
    title(['Mean distance error',' ',conds{condi}])
    ylabel('Distance error per dipole (mm)')
end





function D=getDist(pnt,pnti,pntX)
D=[];
distances=[];
errVec(1)=length(pntX);
errVec(1,2)=sum(ismember(pntX,pnti));
% get all distances
for Xi=1:length(pntX)
    distances(Xi,1:length(pnti))=sqrt(sum(power(pnt(pnti,:)-repmat(pnt(pntX(Xi),:),length(pnti),1),2),2));
end
% find pairs of nearby sources and sum the error
distances(distances==0)=max(distances);
D=min(distances(:));



