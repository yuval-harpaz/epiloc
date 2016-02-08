function [miss,missR]=Simulate_tp_fn_table2(noiseFactor,Err)
    Rpower=100;
for Ndip=1:5;
    load(['results1_',num2str(Ndip),'_',num2str(Rpower),'_',num2str(noiseFactor),'.mat'])
    load (['Rcorr_Dist_',num2str(Ndip),'_',num2str(noiseFactor),'_',num2str(Rpower),'.mat'])
       
    load pnt
    load gain1
    load layer

    % simulate best fit, more than 2 dipoles
    dist=[];
    distR=[];
    for permi=1:1000
        pairs=getErrors(Dist{permi,2});
        dist=[dist;Dist{permi,2}(pairs)];
        pairs=getErrors(Dist{permi,3});
        distR=[distR;Dist{permi,3}(pairs)];
    end

    missAvgA=sum(uint8(results(:,5)-Ndip))./sum(results(:,5))*100;
    missAvgB=sum(dist>=Err)./sum(results(:,5))*100;
    missRA=sum(uint8(results(:,9)-Ndip))./sum(results(:,9))*100;
    missRB=sum(distR>=Err)./sum(results(:,9))*100;
    missR(Ndip)=missRB+missRA;
    miss(Ndip)=missAvgB+missAvgA;
    
    load(['resultsSeq1_',num2str(Ndip),'_',num2str(noiseFactor),'.mat'])
    load (['SEQ_Rcorr_Dist_',num2str(Ndip),'_',num2str(noiseFactor),'.mat'])
    
        % simulate best fit, more than 2 dipoles
    distS=[];
    for permi=1:1000
        pairs=getErrors(Dist{permi,1});
        distS=[distS;Dist{permi,1}(pairs)];
    end

    missSEQA=sum(uint8(results(:,1)-Ndip))./sum(results(:,1))*100;
    missSEQB=sum(dist>=Err)./sum(results(:,1))*100;
    missSEQ(Ndip)=missSEQB+missSEQA;

    figure;
    bar([missAvgB,missAvgA;missRB,missRA; missSEQB,missSEQA],'stacked')
    set(gca, 'xticklabel', [{'RIMDA'},{'BEST FIT'}, {'SEQ'}])
    xlim([0.5 3.5])
    legend('distant','superfluous')
    title([num2str(Ndip),' dipoles, false positive for ',num2str(Err),'mm or more'])
    ylabel('the retio of false positive dipoles (%)')
    ylim([0 25])

    
end



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
    if size(distances,1)>size(distances,2) % more rows
        distances(:,x)=max(minD(end));
    else
        distances(y,:)=max(minD(end));
    end
    pairs(y,x)=true;
end



% binSum=[];
% binNacc=[];
% binDacc=[];
% binD=[];
% conds={'MED','AVG','R'};
% fff=0;
% for condi=1:3
%     for bini=3:15
%         rowi=logical((results(:,3)>=(bini-10)).*(results(:,3)<(bini)));
%         binSum(bini-2)=sum(rowi);
%         binNacc(bini-2)=mean(results(rowi,1+fff)); % ratio of number of dipoles (1 of 2 etc)
%         binDacc(bini-2)=mean(results(rowi,2+fff)); % ratio of right location
%         binD(bini-2)=mean(results(rowi,3+fff));
%     end
%     figure;bar([35:10:155],binNacc)
%     xlabel('Distance between dipoles (mm)')
%     ylabel('detected dipoles')
%     title(['How many dipoles were found',' ',conds{condi}])
%     figure;bar([35:10:155],binDacc)
%     xlabel('Distance between dipoles (mm)')
%     ylabel('correctly detected dipoles')
%     title(['How many dipoles were correctly localized',' ',conds{condi}])
%     figure;bar([35:10:155],binD)
%     xlabel('Distance between dipoles (mm)')
%     title(['Mean distance error',' ',conds{condi}])
%     ylabel('Distance error per dipole (mm)')
% end