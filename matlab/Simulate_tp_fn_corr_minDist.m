function [pairCorr,pairCorrR, pairCorrS, his, hisR, hisS]=Simulate_tp_fn_corr_minDist(noiseFactor)
Rpower=100;
minDist=25;
for Ndip=2:5;
    
    load pnt
    load gain1
    load layer

    load(['results1_',num2str(Ndip),'_',num2str(Rpower),'_',num2str(noiseFactor),'.mat'])
    load (['Rcorr_Dist_',num2str(Ndip),'_',num2str(noiseFactor),'_',num2str(Rpower),'.mat'])
    % simulate best fit, more than 2 dipoles
    %dist=[];
    %distR=[];
    count=0;
    countR=0;
    pairCorr{Ndip}=[];
    pairCorrR{Ndip}=[];
    for permi=1:1000
        pairs=getErrors(Dist{permi,2});
        
        vec=Dist{permi,2}(pairs);
        vecDist=zeros(1,Ndip);
        pairs_ind=find(sum(pairs));
        vecDist(pairs_ind)=vec;
        vecDist=vecDist<=minDist; 
        vecDist=~vecDist;
        
        %dist=[dist;Dist{permi,2}(pairs)];
        rCorr=Rcorr(:,:,permi);
        rCorr(logical(eye(Ndip)))=nan;
        
        for i=1:Ndip
            
            j=Ndip;
            while j>i
                count=count+1;                
                pairCorr{Ndip}(count,1)=rCorr(i,j);
                pairCorr{Ndip}(count,2)=((sum(pairs(:,i))-vecDist(i))+(sum(pairs(:,j))-vecDist(j)))/2;               
                j=j-1;
            end
            
        end
        
        pairs=getErrors(Dist{permi,3});
        
        vec=Dist{permi,3}(pairs);
        vecDist=zeros(1,Ndip);
        pairs_ind=find(sum(pairs));
        vecDist(pairs_ind)=vec;
        vecDist=vecDist<=minDist; 
        vecDist=~vecDist;
        
        %dist=[dist;Dist{permi,2}(pairs)];
        rCorr=Rcorr(:,:,permi);
        rCorr(logical(eye(Ndip)))=nan;
        
        for i=1:Ndip
            j=Ndip;
            while j>i
                countR=countR+1;
                pairCorrR{Ndip}(countR,1)=rCorr(i,j);
                pairCorrR{Ndip}(countR,2)=((sum(pairs(:,i))-vecDist(i))+(sum(pairs(:,j))-vecDist(j)))/2;               
                j=j-1;
            end
        end
    end
    %     avg=histc(pairCorr{Ndip}(:,1),[0:0.05:1]);
    stepi=2/7;
    bins=[-1:stepi:1];
    his=NaN(length(bins)-1,1);
    for bini=1:length(his)
        linei=pairCorr{Ndip}(:,1)>=bins(bini) & pairCorr{Ndip}(:,1)<bins(bini+1);
        his(bini)=mean(pairCorr{Ndip}(linei,2));
    end
    %     for jj=1:length(pairCorr{Ndip}(:,1))
%     b=pairCorr{Ndip}(:,1)./stepi;
%     btemp=b;
%     btemp(b<0)=floor(b(b<0));
%     btemp(b>0)=ceil(b(b>0));
%         
%     for bi=1:length(bins)
%         bii=b==(bins(bi)./stepi);
%         his(bi)=mean(pairCorr{Ndip}(bii,2));
%     end
    
    hisR=NaN(length(bins)-1,1);
    %     for jj=1:length(pairCorr{Ndip}(:,1))
    for bini=1:length(his)
        linei=pairCorrR{Ndip}(:,1)>=bins(bini) & pairCorrR{Ndip}(:,1)<bins(bini+1);
        hisR(bini)=mean(pairCorrR{Ndip}(linei,2));
    end
  
    
    load(['resultsSeq1_',num2str(Ndip),'_',num2str(noiseFactor),'.mat'])
    load (['SEQ_Rcorr_Dist_',num2str(Ndip),'_',num2str(noiseFactor),'.mat'])
    countS=0;
    pairCorrS{Ndip}=[];
    for permi=1:1000
        
        pairs=getErrors(Dist{permi,1});
        
        vec=Dist{permi,1}(pairs);
        vecDist=zeros(1,Ndip);
        pairs_ind=find(sum(pairs));
        vecDist(pairs_ind)=vec;
        vecDist=vecDist<=minDist; 
        vecDist=~vecDist;
        
        %dist=[dist;Dist{permi,2}(pairs)];
        rCorr=Rcorr(:,:,permi);
        rCorr(logical(eye(Ndip)))=nan;
        
        for i=1:Ndip
            j=Ndip;
            while j>i
                countS=countS+1;
                pairCorrS{Ndip}(countS,1)=rCorr(i,j);
                pairCorrS{Ndip}(countS,2)=((sum(pairs(:,i))-vecDist(i))+(sum(pairs(:,j))-vecDist(j)))/2;               
                j=j-1;
            end
        end
        
    end
    %     avg=histc(pairCorr{Ndip}(:,1),[0:0.05:1]);
    hisS=NaN(length(bins)-1,1);
    for bini=1:length(his)
        linei=pairCorrS{Ndip}(:,1)>=bins(bini) & pairCorrS{Ndip}(:,1)<bins(bini+1);
        hisS(bini)=mean(pairCorrS{Ndip}(linei,2));
    end
     
    figure;
    plot(bins(1:end-1)+stepi/2,his, 'r*-');
    hold on;
    plot(bins(1:end-1)+stepi/2,hisR, 'bo-');     
    plot(bins(1:end-1)+stepi/2,hisS, 'gd-');
    ylim([0 1]);
    
%     sum(avg)
    %figure;hist(distR,100)
    %     figure;
    %     first=true;
    %     for err=15:2.5:40
    %         plot(err,sum(dist>=err)./length(dist),'og')
    %         hold on
    %         plot(err,sum(distR>=err)./length(distR),'^r')
    %         if first
    %             legend('RIMDA','Best Fit')
    %             first=false;
    %         end
    %     end
    %     xlim([13 42])
    %     ylim([-0.05 0.15])
    
    %err=25;
    % save(['results1_',num2str(Ndip),'_',num2str(Rpower),'_',num2str(noiseFactor)],'results','input')
    % marikVirtual29plot(input,pnt,results);
    %     missAvgA=sum(uint8(results(:,5)-Ndip))./sum(results(:,5))*100;
    %     missAvgB=sum(dist>=Err)./sum(results(:,5))*100;
    %     missRA=sum(uint8(results(:,9)-Ndip))./sum(results(:,9))*100;
    %     missRB=sum(distR>=Err)./sum(results(:,9))*100;
    %     missR(Ndip)=missRB+missRA;
    %     miss(Ndip)=missAvgB+missAvgA;
    %     figure;
    %     bar([missAvgB,missAvgA;missRB,missRA],'stacked')
    %     set(gca, 'xticklabel', [{'RIMDA'},{'BEST FIT'}])
    %     xlim([0.5 2.5])
    %     legend('distant','superfluous')
    %     title([num2str(Ndip),' dipoles, false positive for ',num2str(Err),'mm or more'])
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