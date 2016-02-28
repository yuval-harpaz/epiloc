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
for noiseFactor=[0.1 0.3]
    for Ndip=1:5
        [results, Rcorr, Dist]=marikVirtual31(Ndip, noiseFactor);
        eval(['save SEQ_Rcorr_Dist_',num2str(Ndip),'_', num2str(noiseFactor),'.mat Dist Rcorr']);
    end
    clear all
%     cd ..
    load pnt
%     cd('./SEQ')
end

cd ..
cd ./temp
colors=eye(3);
h1=figure;
h2=figure;
ls={'-', '--'};
noise_count=0;
for noiseFactor=[0.1 0.3]
    [miss,missR, mat_miss]=Simulate_tp_fn_table2(noiseFactor,25, 0);
%     close all;
    
    m1=zeros(3,length(mat_miss));
    m2=zeros(3,length(mat_miss));
    for jj=1:length(mat_miss)
        m1(:,jj)=mat_miss{jj}(:,1);
        m2(:,jj)=mat_miss{jj}(:,2);
    end
    
% %     figure;
% %     bar(m1')
% %     hold on
% %     errorbar(m1', zeros(5,3));
% %     xlabel('N dipoles');
% %     legend([{'RIMDA'},{'BEST FIT'}, {'SEQ'}]);
% %     ylabel('the ratio of false positive dipoles due to distance(%)')
% %     title('Distant')
% %     ylim([0 20])
% %     
% %     figure;
% %     bar(m2')    
% %     hold on
% %     plot(m2')
% %     xlabel('N dipoles');
% %     legend([{'RIMDA'},{'BEST FIT'}, {'SEQ'}]);
% %     ylabel('the ratio of false positive dipoles due to superfluous(%)')
% %     title('Superfluous')
% %     ylim([0 20])   
    
    noise_count=noise_count+1;
    for jj=1:3
        figure(h1)
        hold on
        plot(m1(jj,:)','color', colors(jj,:), 'LineStyle',ls{noise_count}, 'marker', '*');
        figure(h2)
        hold on
        plot(m2(jj,:)','color', colors(jj,:), 'LineStyle',ls{noise_count}, 'marker', '*');
    end
    
    
end
figure(h1);
ylim([0 17])
xlim([0.75 5.25])
set(gca, 'xtick', 1:5); 
xlabel('N dipoles');
legend([{'RIMDA'},{'BEST FIT'}, {'SEQ'}]);
ylabel('the ratio of false positive dipoles due to distance(%)')
title('Distant')
figure(h2);
ylim([0 17])         
xlim([0.75 5.25])
set(gca, 'xtick', 1:5); 
xlabel('N dipoles');
legend([{'RIMDA'},{'BEST FIT'}, {'SEQ'}]);
ylabel('the ratio of false positive dipoles due to superfluous(%)')
title('Superfluous')


load pnt
load('resultsSeq1_4_0.3.mat')
[SEQ]=marikVirtual31plot(input,pnt,results)
load('results1_4_100_0.3.mat')
[R,MED]=marikVirtual29plot(input,pnt,results)
[R MED SEQ']


load pnt
load('resultsSeq1_2_0.1.mat')
[SEQ]=marikVirtual31plot(input,pnt,results)
load('resultsSeq1_2_0.3.mat')
[SEQ2]=marikVirtual31plot(input,pnt,results)
[SEQ' SEQ2']

%%
% TP - correctly identified dipoles (hits)
% FN - uncorrectly unidentified dipoles (misses)
% FP - uncorrectly identified dipoles (distant + superflous) (false alarms)

% number of dipoles placed
% number of dipoles identified = false alarms + hits
% hits = (number of dipoles identified) - false alarms

% sensitivity = hits/(number of dipoles placed)
% PV+ = hits/(number of dipoles identified)

%%
% TN (PER NUMBER OF DIPOLES PLACED) - number of iterations all dipoles were correctly identified
% FP (PER NUMBER OF DIPOLES PLACED) - false alarms

% specificity = 
% (number of iterations all dipoles were correctly identified)/[(number of iterations all dipoles were correctly identified) + false alarms]
% PV- = 
% (number of iterations all dipoles were correctly identified)/[(number of iterations all dipoles were correctly identified) + misses]

% misses = (number of dipoles placed) - hits

