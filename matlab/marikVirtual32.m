%% plot figures for paper

%% Somatosensory Marik
% load /home/yuval/epiloc/data/avgFilt.mat
load /home/yuval/epiloc/data/MHR
cd /home/yuval/Data/marik/som2/1
% MHR=avg1_handR(:,138);
% avg1_handR.time(138)=0.0334;
[pow,pntR,current,fwdHR]=rimda(MHR);
load /home/yuval/epiloc/data/MHL
% MHL=avg1_handL(:,138);
% avg1_handL.time(138)=0.0334;
[pow,pntL,current,fwdHL]=rimda(MHL);
load /home/yuval/epiloc/data/MF
% MF=avg1_footL(:,180);
% avg1_footL.time(180)=0.0747;
[pow,pntF,current,fwdF]=rimda(MF);
load /home/yuval/epiloc/data/RF
[pow,pntRF,current,fwdRF]=rimda(RF);
load /home/yuval/epiloc/data/RFL
[pow,pntRFL,current,fwdRFL]=rimda(RFL);
save /home/yuval/epiloc/data/pntiAll pnt*

%% Correlations:
cFL=corr(fwdF,fwdHL);
cFR=corr(fwdF,fwdHR);
cLR=corr(fwdHL,fwdHR);

cMFL=corr(MF,MHL);
cMFR=corr(MF,MHR);
cMLR=corr(MHL,MHR);

%% Distances:
dist=norm(pnt(pntR,:)-pnt(pntRF(2),:),2);
% dist=norm(pnt(pntR,:)-pnt(pntRF(1),:),2)

%% Auditory Marik
cd /home/yuval/Data/marik/marikAud
load avgSt.mat
M=avgSt.avg(:,353);
% avgSt.time(353)=0.1465
[pow,pntAud,current]=rimda(M);
save /home/yuval/epiloc/data/pntAud pntAud
% If from Odd ball (Currently we decided to present only Marik - standart) 
% 1st componenent:
% avgOdd.time(343)=0.1366
% 2nd componenent:
% avgOdd.time(397)=0.1897

% cd '/home/oshrit/MyDocuments/DATA/Marik/epiloc/data/Figures/In_Afni/Marik/Aud'
% load('avgSt.mat')
% figure; plot(avgSt.time, avgSt.avg, 'b');
% figure; plot(avgSt.avg', 'b');
% load('avgOdd.mat')
% figure; plot(avgOdd.time, avgOdd.avg, 'r');
% figure; plot(avgOdd.avg', 'r');
% M=avgOdd.avg(:,397); figure; topoplot248(M);
% M=avgOdd.avg(:,343); figure; topoplot248(M);


%% Somatosensory Yuval
cd /home/yuval/epiloc/data/yuvSom
load avgFilt
RLF=[152,157,300];
% avg1_handR.time(RLF(1)) = 0.0472
% avg1_handL.time(RLF(2)) = 0.0521
% avg1_footL.time(RLF(3)) = 0.1927
% Hand R
% figure; plot(avg1_handR.avg');
MYHR=avg1_handR.avg(:,RLF(1));
% figure;topoplot248(MYHR);
[pow,pntYHR,current,fwdYHR]=rimda(MYHR);
% Hand L
MYHL=avg1_handL.avg(:,RLF(2));
% figure;topoplot248(MYHL);
[pow,pntYHL,current,fwdYHL]=rimda(MYHL);
% LEG  
MYF=avg1_footL.avg(:,RLF(3));
% figure;topoplot248(MYF);
[pow,pntYF,current,fwdYF]=rimda(MYF);
% Hand R + foot
MYRF=avg1_handR.avg(:,RLF(1))+avg1_footL.avg(:,RLF(3));
% figure;topoplot248(MYRF);
[pow,pntYRF,current]=rimda(MYRF);
% Hand R + Hand L + foot
MYRFL=avg1_handR.avg(:,RLF(1))+avg1_handL.avg(:,RLF(2))+avg1_footL.avg(:,RLF(3));
% figure;topoplot248(MYHR);
[pow,pntYRFL,current]=rimda(MYRFL);
% pnt(611,:)
%    -3.3202   36.0627   85.1601
% pnt(610,:)
%    -3.3202  -29.0331   93.1778
% pnt(624,:)
%    -8.1466   29.4913   88.9083
% pnt(308,:)
%   -14.3983  -13.4400  106.7133
% dist=norm(pnt(769,:)-pnt(611,:),2);
%%
% % Hand R + foot
% MYR2F=avg1_handR.avg(:,RLF(1))+2*avg1_footL.avg(:,RLF(3));
% % figure;topoplot248(MYRF);
% [pow,pntYR2F,current]=rimda(MYR2F);
% % Hand L + foot
% MYL2F=avg1_handL.avg(:,RLF(2))+2*avg1_footL.avg(:,RLF(3));
% % figure;topoplot248(MYL2F);
% [pow,pntYL2F,current]=rimda(MYL2F);
save /home/yuval/epiloc/data/pntYiAll pntY*

%% Correlations:
cYFL=corr(fwdYF,fwdYHL);
cYFR=corr(fwdYF,fwdYHR);
cYLR=corr(fwdYHL,fwdYHR);
% cYFL =  -0.9413
% cYFR =   0.0946
% cYLR =  -0.1798
save /home/yuval/epiloc/data/corrY cY* fwdY*
MYLF=avg1_handL.avg(:,RLF(2))+avg1_footL.avg(:,RLF(3));
% figure;topoplot248(MYLF);
[pow,pntYLF,current]=rimda(MYLF);
%% High diff in intensity:
figure;topoplot248(MYF);caxis([-0.26 0.44]*10^-12);
figure;topoplot248(MYHL);caxis([-0.26 0.44]*10^-12);
figure;topoplot248(MYHR);caxis([-0.26 0.44]*10^-12);


%% Auditory Yuval
cd /home/yuval/Data/marik/yuval/4
load avgSt.mat
% figure; plot(avgSt.avg');
% Yuval's first componenet does not include 2 dipoles - only 1. 
MY=avgSt.avg(:,385); 
% figure; topoplot248(MY);
[pow,pntYAud,current]=rimda(MY);
save /home/yuval/epiloc/data/pntYAud pntYAud
% MY=avgSt.avg(:,443); % 438
% avgSt.time(443)=0.2349
% We decided not to present Yuval's auditory 
% but if we did it will have been only the 2nd componenet
% Same oreintation as Marik's 2nd componenet (and anticrrelative to 1st
% componenet)
MY=avgOdd.avg(:,427); 
% figure; topoplot248(MY);
[pow,pntYAud,current]=rimda(MY);
% avgOdd.time(427)=0.2192

%% write sphere dipole to BRIK
cd /home/yuval/epiloc/data
%% Somatosensory
% Marik
load pntiAll
load pntSom
%AIR order 
x=num2str(-pnt(pntR,1));y=num2str(pnt(pntR,3));z=num2str(pnt(pntR,2));
[~,w]=afnix(['3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix R'])
afni

x=num2str(-pnt(pntL,1));y=num2str(pnt(pntL,3));z=num2str(pnt(pntL,2));
[~,w]=afnix(['3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix L'])

x=num2str(-pnt(pntF,1));y=num2str(pnt(pntF,3));z=num2str(pnt(pntF,2));
[~,w]=afnix(['3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix F'])
for ii=1:2
    x=num2str(-pnt(pntRF(ii),1));y=num2str(pnt(pntRF(ii),3));z=num2str(pnt(pntRF(ii),2));
    [~,w]=afnix(['3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix RF',num2str(ii)])
end
[~,w]=afnix('3dcalc -a RF1+orig -b RF2+orig -exp "a+b" -prefix RF')
for ii=1:3
    x=num2str(-pnt(pntRFL(ii),1));y=num2str(pnt(pntRFL(ii),3));z=num2str(pnt(pntRFL(ii),2));
    [~,w]=afnix(['3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix RFL',num2str(ii)])
end
[~,w]=afnix('3dcalc -a RFL1+orig -b RFL2+orig -c RFL3+orig -exp "a+b+c" -prefix RFL')
% [~,w]=afnix('3dcalc -a R+orig -b F+orig -c L+orig -d RFL+orig -exp "a+b+c+d" -prefix RFLRFL')

% Yuval
load('/home/yuval/epiloc/data/yuvSom/pnt.mat')
load pntYiAll
x=num2str(-pnt(pntYHR,1));y=num2str(pnt(pntYHR,3));z=num2str(pnt(pntYHR,2));
[~,w]=afnix(['3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix YR'])
% afni

x=num2str(-pnt(pntYHL,1));y=num2str(pnt(pntYHL,3));z=num2str(pnt(pntYHL,2));
[~,w]=afnix(['3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix YL'])

x=num2str(-pnt(pntYF,1));y=num2str(pnt(pntYF,3));z=num2str(pnt(pntYF,2));
[~,w]=afnix(['3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix YF'])

for ii=1:1 % not 2
    x=num2str(-pnt(pntYRF(ii),1));y=num2str(pnt(pntYRF(ii),3));z=num2str(pnt(pntYRF(ii),2));
    [~,w]=afnix(['3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix YRF',num2str(ii)])
end
% [~,w]=afnix('3dcalc -a YRF1+orig -b YRF2+orig -exp "a+b" -prefix YRF')

for ii=1:2 % not 3
    x=num2str(-pnt(pntYRFL(ii),1));y=num2str(pnt(pntYRFL(ii),3));z=num2str(pnt(pntYRFL(ii),2));
    [~,w]=afnix(['3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix YRFL',num2str(ii)])
end
% [~,w]=afnix('3dcalc -a YRFL1+orig -b YRFL2+orig -c YRFL3+orig -exp "a+b+c" -prefix YRFL')
[~,w]=afnix('3dcalc -a YRFL1+orig -b YRFL2+orig -exp "a+b" -prefix YRFL')


%% Auditory
% Marik
load /home/yuval/Data/marik/marikAud/pnt
cd /home/yuval/epiloc/data/
load pntAud
for ii=1:2
    x=num2str(-pnt(pntAud(ii),1));y=num2str(pnt(pntAud(ii),3));z=num2str(pnt(pntAud(ii),2));
    [~,w]=afnix(['3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix Aud',num2str(ii)])
end

% Yuval
load /home/yuval/Data/marik/yuval/4/pnt
cd /home/yuval/epiloc/data/
load pntYAud
% for ii=1:2
%     x=num2str(-pnt(pntYAud(ii),1));y=num2str(pnt(pntYAud(ii),3));z=num2str(pnt(pntYAud(ii),2));
%     [~,w]=afnix(['3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix YAud3_',num2str(ii)])
% end
% [~,w]=afnix('3dcalc -a YAud3_1+orig -b YAud3_2+orig -exp "a+b" -prefix YAud3')
% load pntYAud_temp
% pntYAud=pntYAud_385;
% pntYAud=pntYAud_385_2;
x=num2str(-pnt(pntYAud,1));y=num2str(pnt(pntYAud,3));z=num2str(pnt(pntYAud,2));
[~,w]=afnix(['3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix YAud']); % 2

   
%% 
% 3dcalc -a R+orig -b F+orig -c L+orig -d RFL+orig -exp "a+b+c+2*d" -prefix RFLRFL

cd /home/yuval/epiloc/data/
% cd /home/oshrit/MyDocuments/DATA/Marik/epiloc/data/drafts
load('MHR.mat')
figure; topoplot248(MHR);
caxis([-1.2*10^-13 1.0*10^-13])
load('MHL.mat')
figure; topoplot248(MHL);
caxis([-1.2*10^-13 1.0*10^-13])
load('MF.mat')
figure; topoplot248(MF);
caxis([-1.2*10^-13 1.0*10^-13])
load('RFL.mat')
figure; topoplot248(RFL);
caxis([-1.2*10^-13 1.0*10^-13])

cd /home/yuval/epiloc/data/yuvSom
load avgFilt
% load /home/oshrit/MyDocuments/DATA/Marik/epiloc/data/yuvSom/avgFilt
RLF=[152,157,300];
MYHR=avg1_handR.avg(:,RLF(1));
MYHL=avg1_handL.avg(:,RLF(2));
MYF=avg1_footL.avg(:,RLF(3));
MYRFL=avg1_handR.avg(:,RLF(1))+avg1_handL.avg(:,RLF(2))+avg1_footL.avg(:,RLF(3));
figure; topoplot248(MYHR);
caxis([-3*10^-13 3.5*10^-13])
figure; topoplot248(MYHL);
caxis([-3*10^-13 3.5*10^-13])
figure; topoplot248(MYF);
caxis([-3*10^-13 3.5*10^-13])
figure; topoplot248(MYRFL);
caxis([-3*10^-13 3.5*10^-13])

cd /home/yuval/Data/marik/marikAud
load avgSt.mat
% load /home/oshrit/MyDocuments/DATA/Marik/epiloc/data/In_Afni/Marik/Aud/avgSt.mat
M=avgSt.avg(:,353);
figure; topoplot248(M);
caxis([-1.2*10^-13 1.0*10^-13])

cd /home/yuval/Data/marik/yuval/4
% load avgSt.mat
% % load /home/oshrit/MyDocuments/DATA/Marik/epiloc/data/In_Afni/Yuval/Aud/avgSt.mat
% MY=avgSt.avg(:,385); %394
% figure; topoplot248(MY);
% caxis([-1.2.*10^-13 1.5*10^-13])
load avgOdd.mat
% load /home/oshrit/MyDocuments/DATA/Marik/epiloc/data/In_Afni/Yuval/Aud/avgOdd.mat
MY=avgOdd.avg(:,427); %394
figure; topoplot248(MY);
caxis([-3.75.*10^-13 3.5*10^-13])

% hs=ft_read_headshape('/home/oshrit/MyDocuments/DATA/Marik/epiloc/data/yuvSom/hs_file');
% pnti=726;
% figure;
% plot3pnt(hs,'.k');
% hold on
% scatter3pnt(pnt(pnti,:),25,[1 0 0])
% scatter3pnt(pnt,1,[1 1 1])
% colorbar off
% rotate3d on

%%
h = get(0,'children');
for i=1:length(h)
  saveas(h(i), ['figure' num2str(i)], 'fig');
  saveas(h(i), ['figure' num2str(i)], 'png');
end

%%
cd('/home/oshrit/MyDocuments/DATA/Marik/epiloc/data/Data_for_source/marikSom/');
% Marik somatosensory
load('gain.mat')
load('pnt.mat')
hs_path='/home/oshrit/MyDocuments/DATA/Marik/epiloc/data/Data_for_source/marikSom/hs_file';
load('avgFilt.mat');
% R L F and RLF
MHR=avg1_handR.avg(:,138);
MHL=avg1_handL.avg(:,138);
MF=avg1_footL.avg(:,180);
RFL=MHR+MF+MHL;
% figure; topoplot248(RFL);
[pntAvg_HR,pntR_HR,pntSeq_HR]=SEQ_and_Best_fit(MHR, hs_path, pnt, gain);
[pntAvg_HL,pntR_HL,pntSeq_HL]=SEQ_and_Best_fit(MHL, hs_path, pnt, gain);
[pntAvg_F,pntR_F,pntSeq_F]=SEQ_and_Best_fit(MF, hs_path, pnt, gain);
[pntAvg_RFL,pntR_RFL,pntSeq_RFL]=SEQ_and_Best_fit(RFL, hs_path, pnt, gain);
close all
save ./results/pnti_comparison pnt*
clear pntAvg* pntR* pntSeq*

clear all; 
% Marik Auditory
cd('/home/oshrit/MyDocuments/DATA/Marik/epiloc/data/Data_for_source/marikAud');
load('gain.mat')
load('pnt.mat')
hs_path='/home/oshrit/MyDocuments/DATA/Marik/epiloc/data/Data_for_source/marikAud/hs_file';
load('avgSt.mat');
M=avgSt.avg(:,353);
[pntAvg_Aud,pntR_Aud,pntSeq_Aud]=SEQ_and_Best_fit(M, hs_path, pnt, gain);
close all
save ./results/pnti_comparison pnt*
clear pntAvg* pntR* pntSeq*

clear all; 
% Yuval somatosensory
cd('/home/oshrit/MyDocuments/DATA/Marik/epiloc/data/Data_for_source/yuvSom');
load('gain.mat')
load('pnt.mat')
hs_path='/home/oshrit/MyDocuments/DATA/Marik/epiloc/data/Data_for_source/yuvSom/hs_file';
load('avgFilt.mat');
% R L F and RLF
RLF=[152,157,300];
MHR=avg1_handR.avg(:,RLF(1));
MHL=avg1_handL.avg(:,RLF(2));
MF=avg1_footL.avg(:,RLF(3));
RFL=MHR+MF+MHL;
% figure; topoplot248(RFL);
[pntAvg_HR,pntR_HR,pntSeq_HR]=SEQ_and_Best_fit(MHR, hs_path, pnt, gain);
[pntAvg_HL,pntR_HL,pntSeq_HL]=SEQ_and_Best_fit(MHL, hs_path, pnt, gain);
[pntAvg_F,pntR_F,pntSeq_F]=SEQ_and_Best_fit(MF, hs_path, pnt, gain);
[pntAvg_RFL,pntR_RFL,pntSeq_RFL]=SEQ_and_Best_fit(RFL, hs_path, pnt, gain);
close all
save ./results/pnti_comparison pnt*
clear pntAvg* pntR* pntSeq*


%% Creating Joint Afni files for all methods:
%% Marik somatosensory
cd('/home/oshrit/MyDocuments/DATA/Marik/epiloc/data/Data_for_source/marikSom/');
load('./results/pnti_comparison');
cd('/home/oshrit/MyDocuments/DATA/Marik/epiloc/data/Data_for_source/Figures/In_Afni/Marik/Som');
%% RIMDA
x=num2str(-pnt(pntAvg_HR,1));y=num2str(pnt(pntAvg_HR,3));z=num2str(pnt(pntAvg_HR,2));
[~,w]=afnix(['/home/megadmin/abin/3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix Avg_HR']);
x=num2str(-pnt(pntAvg_HL,1));y=num2str(pnt(pntAvg_HL,3));z=num2str(pnt(pntAvg_HL,2));
[~,w]=afnix(['/home/megadmin/abin/3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix Avg_HL']);
x=num2str(-pnt(pntAvg_F,1));y=num2str(pnt(pntAvg_F,3));z=num2str(pnt(pntAvg_F,2));
[~,w]=afnix(['/home/megadmin/abin/3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix Avg_F']);
length(pntAvg_RFL)
for ii=1:length(pntAvg_RFL)
    x=num2str(-pnt(pntAvg_RFL(ii),1));y=num2str(pnt(pntAvg_RFL(ii),3));z=num2str(pnt(pntAvg_RFL(ii),2));
    [~,w]=afnix(['/home/megadmin/abin/3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix Avg_RFL',num2str(ii)]);
end
[~,w]=afnix('/home/megadmin/abin/3dcalc -a Avg_RFL1+orig -b Avg_RFL2+orig -c Avg_RFL3+orig -exp "a+b+c" -prefix Avg_RFL');
[~,w]=afnix('/home/megadmin/abin/3dcalc -a Avg_HR+orig -b Avg_F+orig -c Avg_HL+orig -d Avg_RFL+orig -exp "a+b+c+2*d" -prefix Avg_RFLRFL');
%% Best Fit
length(pntR_HR)
for ii=1:length(pntR_HR)
    x=num2str(-pnt(pntR_HR(ii),1));y=num2str(pnt(pntR_HR(ii),3));z=num2str(pnt(pntR_HR(ii),2));
    [~,w]=afnix(['/home/megadmin/abin/3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix R_HR',num2str(ii)]);
end
[~,w]=afnix('/home/megadmin/abin/3dcalc -a R_HR1+orig -b R_HR2+orig -c R_HR3+orig -d R_HR4+orig -e R_HR5+orig -exp "a+b+c+d+e" -prefix R_HR');
length(pntR_HL)
for ii=1:length(pntR_HL)
    x=num2str(-pnt(pntR_HL(ii),1));y=num2str(pnt(pntR_HL(ii),3));z=num2str(pnt(pntR_HL(ii),2));
    [~,w]=afnix(['/home/megadmin/abin/3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix R_HL',num2str(ii)]);
end
[~,w]=afnix('/home/megadmin/abin/3dcalc -a R_HL1+orig -b R_HL2+orig -c R_HL3+orig -exp "a+b+c" -prefix R_HL');
length(pntR_F)
for ii=1:length(pntR_F)
    x=num2str(-pnt(pntR_F(ii),1));y=num2str(pnt(pntR_F(ii),3));z=num2str(pnt(pntR_F(ii),2));
    [~,w]=afnix(['/home/megadmin/abin/3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix R_F',num2str(ii)]);
end
[~,w]=afnix('/home/megadmin/abin/3dcalc -a R_F1+orig -b R_F2+orig -c R_F3+orig -exp "a+b+c" -prefix R_F');
length(pntR_RFL)
for ii=1:length(pntR_RFL)
    x=num2str(-pnt(pntR_RFL(ii),1));y=num2str(pnt(pntR_RFL(ii),3));z=num2str(pnt(pntR_RFL(ii),2));
    [~,w]=afnix(['/home/megadmin/abin/3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix R_RFL',num2str(ii)]);
end
[~,w]=afnix('/home/megadmin/abin/3dcalc -a R_RFL1+orig -b R_RFL2+orig -c R_RFL3+orig -exp "a+b+c" -prefix R_RFL');
[~,w]=afnix('/home/megadmin/abin/3dcalc -a R_HR+orig -b R_F+orig -c R_HL+orig -d R_RFL+orig -exp "a+b+c+2*d" -prefix R_RFLRFL');
%% SEQ
length(pntSeq_HR)
for ii=1:length(pntSeq_HR)
    x=num2str(-pnt(pntSeq_HR(ii),1));y=num2str(pnt(pntSeq_HR(ii),3));z=num2str(pnt(pntSeq_HR(ii),2));
    [~,w]=afnix(['/home/megadmin/abin/3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix Seq_HR',num2str(ii)]);
end
[~,w]=afnix('/home/megadmin/abin/3dcalc -a Seq_HR1+orig -b Seq_HR2+orig -exp "a+b" -prefix Seq_HR');
length(pntSeq_HL)
for ii=1:length(pntSeq_HL)
    x=num2str(-pnt(pntSeq_HL(ii),1));y=num2str(pnt(pntSeq_HL(ii),3));z=num2str(pnt(pntSeq_HL(ii),2));
    [~,w]=afnix(['/home/megadmin/abin/3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix Seq_HL',num2str(ii)]);
end
[~,w]=afnix('/home/megadmin/abin/3dcalc -a Seq_HL1+orig -b Seq_HL2+orig -exp "a+b" -prefix Seq_HL');
% length(pntSeq_F)
x=num2str(-pnt(pntSeq_F,1));y=num2str(pnt(pntSeq_F,3));z=num2str(pnt(pntSeq_F,2));
[~,w]=afnix(['/home/megadmin/abin/3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix Seq_F']);
length(pntSeq_RFL)
for ii=1:length(pntSeq_RFL)
    x=num2str(-pnt(pntSeq_RFL(ii),1));y=num2str(pnt(pntSeq_RFL(ii),3));z=num2str(pnt(pntSeq_RFL(ii),2));
    [~,w]=afnix(['/home/megadmin/abin/3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix Seq_RFL',num2str(ii)]);
end
[~,w]=afnix('/home/megadmin/abin/3dcalc -a Seq_RFL1+orig -b Seq_RFL2+orig -c Seq_RFL3+orig -d Seq_RFL4+orig -exp "a+b+c+d" -prefix Seq_RFL');
[~,w]=afnix('/home/megadmin/abin/3dcalc -a Seq_HR+orig -b Seq_F+orig -c Seq_HL+orig -d Seq_RFL+orig -exp "a+b+c+2*d" -prefix Seq_RFLRFL');


%% Yuval  somatosensory
cd('/home/oshrit/MyDocuments/DATA/Marik/epiloc/data/Data_for_source/yuvSom/');
load('./results/pnti_comparison');
cd('/home/oshrit/MyDocuments/DATA/Marik/epiloc/data/Data_for_source/Figures/In_Afni/Yuval/Som');
%% RIMDA
x=num2str(-pnt(pntAvg_HR,1));y=num2str(pnt(pntAvg_HR,3));z=num2str(pnt(pntAvg_HR,2));
[~,w]=afnix(['/home/megadmin/abin/3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix Avg_HR']);
x=num2str(-pnt(pntAvg_HL,1));y=num2str(pnt(pntAvg_HL,3));z=num2str(pnt(pntAvg_HL,2));
[~,w]=afnix(['/home/megadmin/abin/3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix Avg_HL']);
x=num2str(-pnt(pntAvg_F,1));y=num2str(pnt(pntAvg_F,3));z=num2str(pnt(pntAvg_F,2));
[~,w]=afnix(['/home/megadmin/abin/3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix Avg_F']);
length(pntAvg_RFL)
for ii=1:length(pntAvg_RFL)
    x=num2str(-pnt(pntAvg_RFL(ii),1));y=num2str(pnt(pntAvg_RFL(ii),3));z=num2str(pnt(pntAvg_RFL(ii),2));
    [~,w]=afnix(['/home/megadmin/abin/3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix Avg_RFL',num2str(ii)]);
end
[~,w]=afnix('/home/megadmin/abin/3dcalc -a Avg_RFL1+orig -b Avg_RFL2+orig -exp "a+b" -prefix Avg_RFL');
[~,w]=afnix('/home/megadmin/abin/3dcalc -a Avg_HR+orig -b Avg_F+orig -c Avg_HL+orig -d Avg_RFL+orig -exp "a+b+c+2*d" -prefix Avg_RFLRFL');
%% Best Fit
% length(pntR_HR)
x=num2str(-pnt(pntR_HR,1));y=num2str(pnt(pntR_HR,3));z=num2str(pnt(pntR_HR,2));
[~,w]=afnix(['/home/megadmin/abin/3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix R_HR']);
% length(pntR_HL)
x=num2str(-pnt(pntR_HL,1));y=num2str(pnt(pntR_HL,3));z=num2str(pnt(pntR_HL,2));
[~,w]=afnix(['/home/megadmin/abin/3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix R_HL']);
length(pntR_F)
for ii=1:length(pntR_F)
    x=num2str(-pnt(pntR_F(ii),1));y=num2str(pnt(pntR_F(ii),3));z=num2str(pnt(pntR_F(ii),2));
    [~,w]=afnix(['/home/megadmin/abin/3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix R_F',num2str(ii)]);
end
[~,w]=afnix('/home/megadmin/abin/3dcalc -a R_F1+orig -b R_F2+orig -c R_F3+orig -exp "a+b+c" -prefix R_F');
length(pntR_RFL)
for ii=1:length(pntR_RFL)
    x=num2str(-pnt(pntR_RFL(ii),1));y=num2str(pnt(pntR_RFL(ii),3));z=num2str(pnt(pntR_RFL(ii),2));
    [~,w]=afnix(['/home/megadmin/abin/3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix R_RFL',num2str(ii)]);
end
[~,w]=afnix('/home/megadmin/abin/3dcalc -a R_RFL1+orig -b R_RFL2+orig -c R_RFL3+orig -exp "a+b+c" -prefix R_RFL');
[~,w]=afnix('/home/megadmin/abin/3dcalc -a R_HR+orig -b R_F+orig -c R_HL+orig -d R_RFL+orig -exp "a+b+c+2*d" -prefix R_RFLRFL');
%% SEQ
% length(pntSeq_HR)
x=num2str(-pnt(pntSeq_HR,1));y=num2str(pnt(pntSeq_HR,3));z=num2str(pnt(pntSeq_HR,2));
[~,w]=afnix(['/home/megadmin/abin/3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix Seq_HR']);
% length(pntSeq_HL)
x=num2str(-pnt(pntSeq_HL,1));y=num2str(pnt(pntSeq_HL,3));z=num2str(pnt(pntSeq_HL,2));
[~,w]=afnix(['/home/megadmin/abin/3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix Seq_HL']);
% length(pntSeq_F)
x=num2str(-pnt(pntSeq_F,1));y=num2str(pnt(pntSeq_F,3));z=num2str(pnt(pntSeq_F,2));
[~,w]=afnix(['/home/megadmin/abin/3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix Seq_F']);
length(pntSeq_RFL)
for ii=1:length(pntSeq_RFL)
    x=num2str(-pnt(pntSeq_RFL(ii),1));y=num2str(pnt(pntSeq_RFL(ii),3));z=num2str(pnt(pntSeq_RFL(ii),2));
    [~,w]=afnix(['/home/megadmin/abin/3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix Seq_RFL',num2str(ii)]);
end
[~,w]=afnix('/home/megadmin/abin/3dcalc -a Seq_RFL1+orig -b Seq_RFL2+orig -c Seq_RFL3+orig -exp "a+b+c" -prefix Seq_RFL');
[~,w]=afnix('/home/megadmin/abin/3dcalc -a Seq_HR+orig -b Seq_F+orig -c Seq_HL+orig -d Seq_RFL+orig -exp "a+b+c+2*d" -prefix Seq_RFLRFL');

%% Marik Auditory
cd('/home/oshrit/MyDocuments/DATA/Marik/epiloc/data/Data_for_source/marikAud');
load('./results/pnti_comparison');
cd('/home/oshrit/MyDocuments/DATA/Marik/epiloc/data/Data_for_source/Figures/In_Afni/Marik/Aud');
%% RIMDA
length(pntAvg_Aud)
for ii=1:length(pntAvg_Aud)
    x=num2str(-pnt(pntAvg_Aud(ii),1));y=num2str(pnt(pntAvg_Aud(ii),3));z=num2str(pnt(pntAvg_Aud(ii),2));
    [~,w]=afnix(['/home/megadmin/abin/3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix Avg_Aud',num2str(ii)]);
end
[~,w]=afnix('/home/megadmin/abin/3dcalc -a Avg_Aud1+orig -b Avg_Aud2+orig -exp "a+b" -prefix Avg_Aud');
%% Best Fit
length(pntR_Aud)
for ii=1:length(pntR_Aud)
    x=num2str(-pnt(pntR_Aud(ii),1));y=num2str(pnt(pntR_Aud(ii),3));z=num2str(pnt(pntR_Aud(ii),2));
    [~,w]=afnix(['/home/megadmin/abin/3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix R_Aud',num2str(ii)]);
end
[~,w]=afnix('/home/megadmin/abin/3dcalc -a R_Aud1+orig -b R_Aud2+orig -c R_Aud3+orig -exp "a+b+c" -prefix R_Aud');
%% SEQ
length(pntSeq_Aud)
for ii=1:length(pntSeq_Aud)
    x=num2str(-pnt(pntSeq_Aud(ii),1));y=num2str(pnt(pntSeq_Aud(ii),3));z=num2str(pnt(pntSeq_Aud(ii),2));
    [~,w]=afnix(['/home/megadmin/abin/3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix Seq_Aud',num2str(ii)]);
end
[~,w]=afnix('/home/megadmin/abin/3dcalc -a Seq_Aud1+orig -b Seq_Aud2+orig -c Seq_Aud3+orig -d Seq_Aud4+orig -exp "a+b+c+d" -prefix Seq_Aud');