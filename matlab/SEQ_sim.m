try 
    cd ('/home/yuval/Data/marik/som2/talk');
catch err
    cd('/home/oshrit/MyDocuments/DATA/Marik/epiloc/data/sim');
end
load pnt

try 
    cd('/home/oshrit/MyDocuments/DATA/Marik/epiloc/data/sim/SEQ/');
catch err
    cd('./');
end

%%
% % for noiseFactor=[0.1 0.3]
% %     for Ndip=1:5
% %         [results, Rcorr, Dist]=marikVirtual31(Ndip, noiseFactor);
% %         eval(['save SEQ_Rcorr_Dist_',num2str(Ndip),'_', num2str(noiseFactor),'.mat Dist Rcorr']);
% %     end
% %     clear all
% % %     cd ..
% %     load pnt
% % %     cd('./SEQ')
% % end

%%
cd ..
cd ./temp
colors=eye(3);
ls={'-', '--'};
noise_count=0;
h1=5; h2=6;
for noiseFactor=[0.1 0.3]
    [miss,missR, mat_miss]=Simulate_tp_fn_table2(noiseFactor,25, 1);
%     close all;
    
    m1=zeros(3,length(mat_miss));
    m2=zeros(3,length(mat_miss));
    for jj=1:length(mat_miss)
        m1(:,jj)=mat_miss{jj}(:,1);
        m2(:,jj)=mat_miss{jj}(:,2);
    end
    
%     figure;
%     bar(m1')
% %     hold on
% %     errorbar(m1', zeros(5,3));
%     xlabel('N dipoles');
%     legend([{'RIMDA'},{'BEST FIT'}, {'SEQ'}]);
%     ylabel('the ratio of false positive dipoles due to distance(%)')
%     title('Distant')
%     ylim([0 27])
%     
%     figure;
%     bar(m2')    
% %     hold on
% %     plot(m2')
%     xlabel('N dipoles');
%     legend([{'RIMDA'},{'BEST FIT'}, {'SEQ'}]);
%     ylabel('the ratio of false positive dipoles due to superfluous(%)')
%     title('Superfluous')
%     ylim([0 27])   
%     
%     noise_count=noise_count+1;
%     for jj=1:3
%         figure(h1)
%         hold on
%         plot(m1(jj,:)','color', colors(jj,:), 'LineStyle',ls{noise_count}, 'marker', '*');
%         figure(h2)
%         hold on
%         plot(m2(jj,:)','color', colors(jj,:), 'LineStyle',ls{noise_count}, 'marker', '*');
%     end
%     
    
end
figure(h1);
ylim([0 27])
xlim([0.75 5.25])
set(gca, 'xtick', 1:5); 
xlabel('N dipoles');
legend([{'RIMDA'},{'BEST FIT'}, {'SEQ'}]);
ylabel('the ratio of false positive dipoles due to distance(%)')
title('Distant')
figure(h2);
ylim([0 27])         
xlim([0.75 5.25])
set(gca, 'xtick', 1:5); 
xlabel('N dipoles');
legend([{'RIMDA'},{'BEST FIT'}, {'SEQ'}]);
ylabel('the ratio of false positive dipoles due to superfluous(%)')
title('Superfluous')

%% correlation
% [pairCorr,pairCorrR, pairCorrS, his, hisR, hisS]=Simulate_tp_fn_corr(0.3);
[pairCorr,pairCorrR, pairCorrS, his, hisR, hisS]=Simulate_tp_fn_corr_minDist(0.3);
[pairCorr,pairCorrR, pairCorrS, his, hisR, hisS]=Simulate_tp_fn_corr_minDist(0.1);


%% Specificity and Sensitivity
% TP - correctly identified dipoles (in a dist range of 25 mm) (hits)
% FN - uncorrectly unidentified dipoles (misses)
% FP - uncorrectly identified dipoles (distant + superflous) (false alarms)

% number of dipoles placed
% number of dipoles identified = false alarms + hits
% hits = (number of dipoles identified) - false alarms

% sensitivity = hits/(number of dipoles placed)
% PV+ = hits/(number of dipoles identified)

% number of dipoles identified = results(1,:) or results(5,:);
% results(9,:); SEQ -> results(1,:) 
% false alarms (FP) = miss(Ndip)
% hits = results(5,:) - miss(Ndip)


%%
% TN (PER NUMBER OF DIPOLES PLACED) - number of iterations all dipoles were correctly identified
% FP (PER NUMBER OF DIPOLES PLACED) - false alarms

% specificity = 
% (number of iterations all dipoles were correctly identified)/[(number of iterations all dipoles were correctly identified) + false alarms]
% PV- = 
% (number of iterations all dipoles were correctly identified)/[(number of iterations all dipoles were correctly identified) + misses]

% misses = (number of dipoles placed) - hits

% The effect of noise factor:
load pnt
load('resultsSeq1_2_0.1.mat')
[SEQ]=marikVirtual31plot(input,pnt,results)
load('resultsSeq1_2_0.3.mat')
[SEQ2]=marikVirtual31plot(input,pnt,results)
[SEQ' SEQ2']


%%
cd ..
cd ./temp
colors=eye(3);
noiseFactor=[0.1 0.3];
sensi=zeros(5,2);
sensiR=zeros(5,2);
sensiSEQ=zeros(5,2);
PVplus=zeros(5,2);
PVplusR=zeros(5,2);
PVplusSEQ=zeros(5,2);
for nF=1:2
    [miss,missR,missSEQ,mat_miss,sensi(:,nF),sensiR(:,nF),sensiSEQ(:,nF),PVplus(:,nF),PVplusR(:,nF),PVplusSEQ(:,nF)]=...
        Simulate_tp_fn_table3(noiseFactor(nF),25, 0);
%     close all;

    figure; 
    bar([sensi(:,nF)';sensiR(:,nF)';sensiSEQ(:,nF)']')
    xlabel('N dipoles');
    legend([{'RIMDA'},{'BEST FIT'}, {'SEQ'}]);
    ylabel('Sensitivity');
        
    figure; 
    bar([PVplus(:,nF)'; PVplusR(:,nF)'; PVplusSEQ(:,nF)']')
    xlabel('N dipoles');
    legend([{'RIMDA'},{'BEST FIT'}, {'SEQ'}]);
    ylabel('Predictive Value +');
    
end
figure; 
bar(sensi)
xlabel('N dipoles');
legend([{'0.1'},{'0.3'}]);
ylabel('Sensitivity');
title('RIMDA');
figure; 
bar(sensiR)
xlabel('N dipoles');
legend([{'0.1'},{'0.3'}]);
ylabel('Sensitivity');
title('BEST FIT');
figure; 
bar(sensiSEQ)
xlabel('N dipoles');
legend([{'0.1'},{'0.3'}]);
ylabel('Sensitivity');
title('SEQ');

figure; 
bar(PVplus)
xlabel('N dipoles');
legend([{'0.1'},{'0.3'}]);
ylabel('Predictive Value +');
title('RIMDA');
figure; 
bar(PVplusR)
xlabel('N dipoles');
legend([{'0.1'},{'0.3'}]);
ylabel('Predictive Value +');
title('BEST FIT');
figure; 
bar(PVplusSEQ)
xlabel('N dipoles');
legend([{'0.1'},{'0.3'}]);
ylabel('Predictive Value +');
title('SEQ');

%%
cd ..
cd ./temp
colors=eye(3);
noiseFactor=[0.1 0.3];
speci=zeros(5,2);
speciR=zeros(5,2);
speciSEQ=zeros(5,2);
PVminus=zeros(5,2);
PVminusR=zeros(5,2);
PVminusSEQ=zeros(5,2);
for nF=1:2
    [miss,missR,missSEQ,mat_miss,speci(:,nF),speciR(:,nF),speciSEQ(:,nF),PVminus(:,nF),PVminusR(:,nF),PVminusSEQ(:,nF)]=...
        Simulate_tp_fn_table4(noiseFactor(nF),25, 0);
%     close all;

    figure; 
    bar([speci(:,nF)';speciR(:,nF)';speciSEQ(:,nF)']')
    xlabel('N dipoles');
    legend([{'RIMDA'},{'BEST FIT'}, {'SEQ'}]);
    ylabel('Specificity');
        
    figure; 
    bar([PVminus(:,nF)'; PVminusR(:,nF)'; PVminusSEQ(:,nF)']')
    xlabel('N dipoles');
    legend([{'RIMDA'},{'BEST FIT'}, {'SEQ'}]);
    ylabel('Predictive Value -');
    
end
figure; 
bar(speci)
xlabel('N dipoles');
legend([{'0.1'},{'0.3'}]);
ylabel('Specificity');
title('RIMDA');
figure; 
bar(speciR)
xlabel('N dipoles');
legend([{'0.1'},{'0.3'}]);
ylabel('Specificity');
title('BEST FIT');
figure; 
bar(speciSEQ)
xlabel('N dipoles');
legend([{'0.1'},{'0.3'}]);
ylabel('Specificity');
title('SEQ');

figure; 
bar(PVminus)
xlabel('N dipoles');
legend([{'0.1'},{'0.3'}]);
ylabel('Predictive Value -');
title('RIMDA');
figure; 
bar(PVminusR)
xlabel('N dipoles');
legend([{'0.1'},{'0.3'}]);
ylabel('Predictive Value -');
title('BEST FIT');
figure; 
bar(PVminusSEQ)
xlabel('N dipoles');
legend([{'0.1'},{'0.3'}]);
ylabel('Predictive Value -');
title('SEQ');

%%
cd('/home/oshrit/Desktop/Marik_Fig/');
h = get(0,'children');
for i=1:length(h)
  saveas(h(i), ['figure' num2str(i)], 'png'); %+16
  saveas(h(i), ['figure' num2str(i)], 'fig');
end


%% Identified dipoles are with a distance error of:
load pnt
Rpower=100;
Ndip_options=1:5;
%noise(0.1,0.3)*Ndip(1:5)*distances(35mm,55mm,75mm)
numDip=zeros(1,6);
distErr=zeros(1,6);
deepErr=zeros(1,6);
noise_count=0;
for noiseFactor=[0.1 0.3]
    
    noise_count=noise_count+1;
    for Ndip=Ndip_options;
        % load('resultsSeq1_4_0.3.mat')
        load(['resultsSeq1_',num2str(Ndip),'_',num2str(noiseFactor),'.mat'])
        [SEQ]=marikVirtual31plot(input,pnt,results);
        
        % load('results1_4_100_0.3.mat')
        load(['results1_',num2str(Ndip),'_',num2str(Rpower),'_',num2str(noiseFactor),'.mat'])
        [R,MED]=marikVirtual29plot(input,pnt,results);
        
        mat=[R MED SEQ']';
        
%         if Ndip~=1

            numDip(end+1:end+3,:)=mat(:,[1:2, 7:8, 13:14]);
            distErr(end+1:end+3,:)=mat(:,[3:4, 9:10, 15:16]);
            deepErr(end+1:end+3,:)=mat(:,[5:6, 11:12, 17:18]);
%         else
%             numDip((1:3)*Ndip*noise_count,1:3)=mat(:,[1:2]);
%             distErr((1:3)*Ndip*noise_count,1:3)=mat(:,[3:4]);
%             deepErr((1:3)*Ndip*noise_count,1:3)=mat(:,[5:6]);
%         end
        
    end
    
    
end
numDip(1,:)=[];
distErr(1,:)=[];
deepErr(1,:)=[];

numDip_BF=numDip(1:3:30,:);
numDip_RIMDA=numDip(2:3:30,:);
numDip_SEQ=numDip(3:3:30,:);

distErr_BF=distErr(1:3:30,:);
distErr_RIMDA=distErr(2:3:30,:);
distErr_SEQ=distErr(3:3:30,:);

deepErr_BF=deepErr(1:3:30,:);
deepErr_RIMDA=deepErr(2:3:30,:);
deepErr_SEQ=deepErr(3:3:30,:);

save numDip_distErr_deepErr numDip_* distErr_* deepErr_*


%% For statistics

% Anova Number of dipoles and distance error 
load pnt
Rpower=100;
%noise(0.1,0.3)*Ndip(1:5)*distances(35mm,55mm,75mm)
noise_count=0;
R=cell(1,3); MED=cell(1,3); SEQ=cell(1,3);
Ndip=1;

matNum=zeros(1000,6);
matDistErr=zeros(1000,6);
for noiseFactor=[0.1 0.3]
    
    noise_count=noise_count+1;
    load(['resultsSeq1_',num2str(Ndip),'_',num2str(noiseFactor),'.mat'])
%     if Ndip==1
%         [SEQ(Ndip,1:3)]=marikVirtual31plot_stat(input,pnt,results);
%     else
        [SEQ(Ndip,:)]=marikVirtual31plot_stat(input,pnt,results);
%     end

    load(['results1_',num2str(Ndip),'_',num2str(Rpower),'_',num2str(noiseFactor),'.mat'])              
%     if Ndip==1
%         [R(Ndip,1:3),MED(Ndip,1:3)]=marikVirtual29plot_stat(input,pnt,results);
%     else
        [R(Ndip,:),MED(Ndip,:)]=marikVirtual29plot_stat(input,pnt,results);
%     end

    if noise_count==1;
        matNum(:, 1)=MED{Ndip,1};
        matNum(:, 2)=R{Ndip,1};
        matNum(:, 3)=SEQ{Ndip,1};
        matDistErr(:, 1)=MED{Ndip,2};
        matDistErr(:, 2)=R{Ndip,2};
        matDistErr(:, 3)=SEQ{Ndip,2};
    elseif noise_count==2;
        matNum(:, 4)=MED{Ndip,1};
        matNum(:, 5)=R{Ndip,1};
        matNum(:, 6)=SEQ{Ndip,1};
        matDistErr(:, 4)=MED{Ndip,2};
        matDistErr(:, 5)=R{Ndip,2};
        matDistErr(:, 6)=SEQ{Ndip,2};
    end
      
end


Ndip_options=2:5;
%noise(0.1,0.3)*Ndip(1:5)*distances(35mm,55mm,75mm)
noise_count=0;
R=cell(length(Ndip_options),9); MED=cell(length(Ndip_options),9); SEQ=cell(length(Ndip_options),9);
matNum35=NaN(1000,6*4);
matDistErr35=NaN(1000,6*4);
vec=[0:6:6*4];
for noiseFactor=[0.1 0.3]
    
    noise_count=noise_count+1;
    for Ndip=Ndip_options;
        load(['resultsSeq1_',num2str(Ndip),'_',num2str(noiseFactor),'.mat'])
        [SEQ(Ndip-1,:)]=marikVirtual31plot_stat(input,pnt,results);
        
        load(['results1_',num2str(Ndip),'_',num2str(Rpower),'_',num2str(noiseFactor),'.mat'])
        [R(Ndip-1,:),MED(Ndip-1,:)]=marikVirtual29plot_stat(input,pnt,results);
        
        if noise_count==1;
            % dipole, noise_level, method
            matNum35(1:length(MED{Ndip-1,1}{1}), vec(Ndip-1)+1)=MED{Ndip-1,1}{1};
            matNum35(1:length(R{Ndip-1,1}{1}), vec(Ndip-1)+2)=R{Ndip-1,1}{1};
            matNum35(1:length(SEQ{Ndip-1,1}{1}), vec(Ndip-1)+3)=SEQ{Ndip-1,1}{1};
            matDistErr35(1:length(MED{Ndip-1,2}{1}),vec(Ndip-1)+1)=MED{Ndip-1,2}{1};
            matDistErr35(1:length(R{Ndip-1,2}{1}), vec(Ndip-1)+2)=R{Ndip-1,2}{1};
            matDistErr35(1:length(SEQ{Ndip-1,2}{1}), vec(Ndip-1)+3)=SEQ{Ndip-1,2}{1};
        elseif noise_count==2;
            % dipole(4), noise_level(2), method(3)
            matNum35(1:length(MED{Ndip-1,1}{1}), vec(Ndip-1)+4)=MED{Ndip-1,1}{1};
            matNum35(1:length(R{Ndip-1,1}{1}), vec(Ndip-1)+5)=R{Ndip-1,1}{1};
            matNum35(1:length(SEQ{Ndip-1,1}{1}), vec(Ndip-1)+6)=SEQ{Ndip-1,1}{1};
            matDistErr35(1:length(MED{Ndip-1,2}{1}), vec(Ndip-1)+4)=MED{Ndip-1,2}{1};
            matDistErr35(1:length(R{Ndip-1,2}{1}), vec(Ndip-1)+5)=R{Ndip-1,2}{1};
            matDistErr35(1:length(SEQ{Ndip-1,2}{1}), vec(Ndip-1)+6)=SEQ{Ndip-1,2}{1};
        end
    end
end
 % SEQ R MED
 % 4 rows - levels of number of dipoles (2,3,4,5 dipoles) 
 % * 
 % 9 columns - 
 % ( [number of dipoles found] [distance error] [depth error] for distances 35(columns 1,2,3), 55(columns 4,5,6) and 75(columns 4,5,6)mm )
 %% I can do also matNum55 matNum75 matDistErr55 matDistErr75
 %% first ask **Maor** if I can run Anova ??

%%
% Chi-test
Err=25;
noiseFactor=[0.1 0.3];
for nF=1:2
    [miss,missR,missSEQ,mat_miss,mAA{nF},mAB{nF},mRA{nF},mRB{nF},mSA{nF},mSB{nF}]=Simulate_tp_fn_table5(noiseFactor(nF),Err, 0);  
end

Ndip_options=1:5;
mBNdipNoise=cell(length(Ndip_options),length(noiseFactor));
mANdipNoise=cell(length(Ndip_options),length(noiseFactor));
for nF=1:2
     for Ndip=Ndip_options;
         mANdipNoise{Ndip,nF}(1:3,1:2)=[mAA{nF}(Ndip,:);mRA{nF}(Ndip,:);mSA{nF}(Ndip,:)];
         mBNdipNoise{Ndip,nF}(1:3,1:2)=[mAB{nF}(Ndip,:);mRB{nF}(Ndip,:);mSB{nF}(Ndip,:)];         
     end
end

%% Observed data
% n1 = 51; N1 = 8193;
% n2 = 74; N2 = 8201;
% x1 = [repmat('a',N1,1); repmat('b',N2,1)];
% x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
% [tbl,chi2stat,pval] = crosstab(x1,x2)

% for nF=1:2
%     for Ndip=Ndip_options;
%     m=mANdipNoise{Ndip,nF};
%     x1 = [repmat('a',m(1,2),1); repmat('b',m(2,2),1); repmat('c',m(3,2),1)];
%     x2 = [repmat(1,m(1,1),1); repmat(2,m(1,2)-m(1,1),1); repmat(1,m(2,1),1); repmat(2,m(2,2)-m(2,1),1); repmat(1,m(3,1),1); repmat(2,m(3,2)-m(3,1),1)];
%     [tbl,chi2stat,pval] = crosstab(x1,x2)
%     end
% end

% ROWS in m: 1- RIMDA; 2- Best Fit; 3- SEQ
format shortEng
pvalA_A_R=zeros(length(noiseFactor),length(Ndip_options));
for nF=1:2
    for Ndip=Ndip_options;
    m=mANdipNoise{Ndip,nF};
    x1 = [repmat('a',m(1,2),1); repmat('b',m(2,2),1)];
    x2 = [repmat(1,m(1,1),1); repmat(2,m(1,2)-m(1,1),1); repmat(1,m(2,1),1); repmat(2,m(2,2)-m(2,1),1)];
    [~,~,pvalA_A_R(nF,Ndip)] = crosstab(x1,x2);
    end
end
pvalA_A_R./30
pvalA_A_S=zeros(length(noiseFactor),length(Ndip_options));
for nF=1:2
    for Ndip=Ndip_options;
    m=mANdipNoise{Ndip,nF};
    x1 = [repmat('a',m(1,2),1); repmat('b',m(3,2),1)];
    x2 = [repmat(1,m(1,1),1); repmat(2,m(1,2)-m(1,1),1); repmat(1,m(3,1),1); repmat(2,m(3,2)-m(3,1),1)];
    [~,~,pvalA_A_S(nF,Ndip)] = crosstab(x1,x2);
    end
end
pvalA_A_S./30
% %     m=mANdipNoise{1,1}; m=mANdipNoise{1,2};
%     x1 = [repmat('a',m(1,2),1); repmat('b',m(3,2),1)];
%     x2 = [repmat(1,m(1,1),1); repmat(2,m(1,2)-m(1,1),1); repmat(1,m(3,1),1); repmat(2,m(3,2)-m(3,1),1)];
%     [~,tbl,p] = crosstab(x1,x2)
pvalA_R_S=zeros(length(noiseFactor),length(Ndip_options));
for nF=1:2
    for Ndip=Ndip_options;
    m=mANdipNoise{Ndip,nF};
    x1 = [repmat('a',m(2,2),1); repmat('b',m(3,2),1)];
    x2 = [repmat(1,m(2,1),1); repmat(2,m(2,2)-m(2,1),1); repmat(1,m(3,1),1); repmat(2,m(3,2)-m(3,1),1)];
    [~,~,pvalA_R_S(nF,Ndip)] = crosstab(x1,x2);
    end
end
pvalA_R_S./30
%%
pvalB_A_R=zeros(length(noiseFactor),length(Ndip_options));
for nF=1:2
    for Ndip=Ndip_options;
    m=mBNdipNoise{Ndip,nF};
    x1 = [repmat('a',m(1,2),1); repmat('b',m(2,2),1)];
    x2 = [repmat(1,m(1,1),1); repmat(2,m(1,2)-m(1,1),1); repmat(1,m(2,1),1); repmat(2,m(2,2)-m(2,1),1)];
    [~,~,pvalB_A_R(nF,Ndip)] = crosstab(x1,x2);
    end
end
pvalB_A_R./30
% %     m=mBNdipNoise{1,1}; 
pvalB_A_S=zeros(length(noiseFactor),length(Ndip_options));
for nF=1:2
    for Ndip=Ndip_options;
    m=mBNdipNoise{Ndip,nF};
    x1 = [repmat('a',m(1,2),1); repmat('b',m(3,2),1)];
    x2 = [repmat(1,m(1,1),1); repmat(2,m(1,2)-m(1,1),1); repmat(1,m(3,1),1); repmat(2,m(3,2)-m(3,1),1)];
    [~,~,pvalB_A_S(nF,Ndip)] = crosstab(x1,x2);
    end
end
pvalB_A_S./30
% %     m=mBNdipNoise{1,1}; m=mBNdipNoise{1,2};
pvalB_R_S=zeros(length(noiseFactor),length(Ndip_options));
for nF=1:2
    for Ndip=Ndip_options;
    m=mBNdipNoise{Ndip,nF};
    x1 = [repmat('a',m(2,2),1); repmat('b',m(3,2),1)];
    x2 = [repmat(1,m(2,1),1); repmat(2,m(2,2)-m(2,1),1); repmat(1,m(3,1),1); repmat(2,m(3,2)-m(3,1),1)];
    [~,~,pvalB_R_S(nF,Ndip)] = crosstab(x1,x2);
    end
end
pvalB_R_S./30
% %     m=mBNdipNoise{1,1}; 

%% Comparing noise level:
pvalA_A=zeros(1,length(Ndip_options));
for Ndip=Ndip_options;
    m1=mANdipNoise{Ndip,1};m2=mANdipNoise{Ndip,2};
    x1 = [repmat('a',m1(1,2),1); repmat('b',m2(1,2),1)];
    x2 = [repmat(1,m1(1,1),1); repmat(2,m1(1,2)-m1(1,1),1); repmat(1,m2(1,1),1); repmat(2,m2(1,2)-m2(1,1),1)];
    [~,~,pvalA_A(Ndip)] = crosstab(x1,x2);
end
pvalA_A./15
pvalA_R=zeros(1,length(Ndip_options));
for Ndip=Ndip_options;
    m1=mANdipNoise{Ndip,1};m2=mANdipNoise{Ndip,2};
    x1 = [repmat('a',m1(2,2),1); repmat('b',m2(2,2),1)];
    x2 = [repmat(1,m1(2,1),1); repmat(2,m1(2,2)-m1(2,1),1); repmat(1,m2(2,1),1); repmat(2,m2(2,2)-m2(2,1),1)];
    [~,~,pvalA_R(Ndip)] = crosstab(x1,x2);
end
pvalA_R./15
pvalA_S=zeros(1,length(Ndip_options));
for Ndip=Ndip_options;
    m1=mANdipNoise{Ndip,1};m2=mANdipNoise{Ndip,2};
    x1 = [repmat('a',m1(3,2),1); repmat('b',m2(3,2),1)];
    x2 = [repmat(1,m1(3,1),1); repmat(2,m1(3,2)-m1(3,1),1); repmat(1,m2(3,1),1); repmat(2,m2(3,2)-m2(3,1),1)];
    [~,~,pvalA_S(Ndip)] = crosstab(x1,x2);
end
pvalA_S./15

pvalB_A=zeros(1,length(Ndip_options));
for Ndip=Ndip_options;
    m1=mBNdipNoise{Ndip,1};m2=mBNdipNoise{Ndip,2};
    x1 = [repmat('a',m1(1,2),1); repmat('b',m2(1,2),1)];
    x2 = [repmat(1,m1(1,1),1); repmat(2,m1(1,2)-m1(1,1),1); repmat(1,m2(1,1),1); repmat(2,m2(1,2)-m2(1,1),1)];
    [~,~,pvalB_A(Ndip)] = crosstab(x1,x2);
end
pvalB_A./15
pvalB_R=zeros(1,length(Ndip_options));
for Ndip=Ndip_options;
    m1=mBNdipNoise{Ndip,1};m2=mBNdipNoise{Ndip,2};
    x1 = [repmat('a',m1(2,2),1); repmat('b',m2(2,2),1)];
    x2 = [repmat(1,m1(2,1),1); repmat(2,m1(2,2)-m1(2,1),1); repmat(1,m2(2,1),1); repmat(2,m2(2,2)-m2(2,1),1)];
    [~,~,pvalB_R(Ndip)] = crosstab(x1,x2);
end
pvalB_R./15
pvalB_S=zeros(1,length(Ndip_options));
for Ndip=Ndip_options;
    m1=mBNdipNoise{Ndip,1};m2=mBNdipNoise{Ndip,2};
    x1 = [repmat('a',m1(3,2),1); repmat('b',m2(3,2),1)];
    x2 = [repmat(1,m1(3,1),1); repmat(2,m1(3,2)-m1(3,1),1); repmat(1,m2(3,1),1); repmat(2,m2(3,2)-m2(3,1),1)];
    [~,~,pvalB_S(Ndip)] = crosstab(x1,x2);
end
pvalB_S./15