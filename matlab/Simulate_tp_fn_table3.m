function [miss,missR,missSEQ,mat_miss,sensi,sensiR,sensiSEQ,PVplus,PVplusR,PVplusSEQ]=Simulate_tp_fn_table3(noiseFactor,Err, isFig)
    Rpower=100;
    Ndip_options=1:5;
    mat_miss=cell(length(Ndip_options),1);
for Ndip=Ndip_options;
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
    missAvgB=sum(dist>=Err)./size(dist,1)*100;
    missRA=sum(uint8(results(:,9)-Ndip))./sum(results(:,9))*100;
    missRB=sum(distR>=Err)./size(distR,1)*100;
    missR(Ndip)=missRB+missRA;
    miss(Ndip)=missAvgB+missAvgA;
    
    % number of dipoles identified = results(1,:) or results(5,:);
    % results(9,:); SEQ -> results(1,:) 
    % false alarms (FP) = miss(Ndip)
    % hits = results(5,:) - miss(Ndip)
    
%     hits(Ndip)=(sum(uint8(results(:,5)))- (miss(Ndip)*(sum(results(:,5))*100)))./sum(results(:,5))*100;
%     hitsR(Ndip)=(sum(uint8(results(:,9)))- (missR(Ndip)*(sum(results(:,9))*100)))./sum(results(:,9))*100;
    hits(Ndip)= 100-miss(Ndip);
    hitsR(Ndip)= 100-missR(Ndip);
    
    % sensitivity = hits/(number of dipoles placed)
    % PV+ = hits/(number of dipoles identified)
    
    sensi(Ndip)=(hits(Ndip)*sum(results(:,5)))/(Ndip*length(results)*100);
    sensiR(Ndip)=(hitsR(Ndip)*sum(results(:,9)))/(Ndip*length(results)*100);
    PVplus(Ndip)=hits(Ndip)./100;
    PVplusR(Ndip)=hitsR(Ndip)./100;
 

    load(['resultsSeq1_',num2str(Ndip),'_',num2str(noiseFactor),'.mat'])
    load (['SEQ_Rcorr_Dist_',num2str(Ndip),'_',num2str(noiseFactor),'.mat'])
    
        % simulate best fit, more than 2 dipoles
    distS=[];
    for permi=1:1000
        pairs=getErrors(Dist{permi,1});
        distS=[distS;Dist{permi,1}(pairs)];
    end

    missSEQA=sum(uint8(results(:,1)-Ndip))./sum(results(:,1))*100;
    missSEQB=sum(distS>=Err)./size(distS,1)*100;
    missSEQ(Ndip)=missSEQB+missSEQA;
    
    %     hitsSEQ(Ndip)=(sum(uint8(results(:,1)))- (missSEQ(Ndip)*(sum(results(:,1))*100)))./sum(results(:,1))*100;
    hitsSEQ(Ndip)= 100 -missSEQ(Ndip);
    sensiSEQ(Ndip)=(hitsSEQ(Ndip)*sum(results(:,1)))/(Ndip*length(results)*100);
    PVplusSEQ(Ndip)=hitsSEQ(Ndip)./100;    

    mat_miss{Ndip}=[missAvgB,missAvgA;missRB,missRA; missSEQB,missSEQA];
    
    if isFig==1
        figure;
        bar(mat_miss{Ndip},'stacked')
        set(gca, 'xticklabel', [{'RIMDA'},{'BEST FIT'}, {'SEQ'}])
        xlim([0.5 3.5])
        legend('distant','superfluous')
        title([num2str(Ndip),' dipoles, false positive for ',num2str(Err),'mm or more'])
        ylabel('the ratio of false positive dipoles (%)')
        ylim([0 25])
    end
        
end


% % 
% % function pairs=getErrors(distances)
% % 
% % distSum=0;
% % depthErr=0;
% % pairs=false(size(distances));
% % for pairi=1:min(size(distances))
% %     [minD,minDi]=sort(distances(:));
% %     x=find(sum(distances==minD(1),1));
% %     x=x(1); % when, say, there are two zero errors
% %     y=find(distances(:,x)==minD(1));
% %     y=y(1);
% %     %y=find(sum(distances==minD(1),2));
% %     %distSum=distSum+distances(y,x);
% %     if size(distances,1)>size(distances,2) % more rows
% %         distances(:,x)=max(minD(end));
% %     else
% %         distances(y,:)=max(minD(end));
% %     end
% %     pairs(y,x)=true;
% % end



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