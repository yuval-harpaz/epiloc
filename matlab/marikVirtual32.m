%% plot figures for paper

%% Somatosensory Marik
load /home/yuval/epiloc/data/MHR
cd /home/yuval/Data/marik/som2/1
[pow,pntR]=rimda(MHR);
load /home/yuval/epiloc/data/MHL
[pow,pntL]=rimda(MHL);
load /home/yuval/epiloc/data/MF
[pow,pntF]=rimda(MF);
load /home/yuval/epiloc/data/RF
[pow,pntRF]=rimda(RF);
load /home/yuval/epiloc/data/RFL
[pow,pntRFL]=rimda(RFL);
save /home/yuval/epiloc/data/pntiAll pnt*

%% Auditory Marik
cd /home/yuval/Data/marik/marikAud
load avgSt.mat
M=avgSt.avg(:,353);
[pow,pntAud,current]=rimda(M);
save /home/yuval/epiloc/data/pntAud pntAud

%% Somatosensory Yuval
cd /home/yuval/epiloc/data/yuvSom
load avgFilt
RLF=[152,157,300];
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
% figure;topoplot248(MYHR);
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
save /home/yuval/epiloc/data/pntYiAll pntY*

%% What is the correlation between HR and F??
cYFL=corr(fwdYF,fwdYHL);
cYFR=corr(fwdYF,fwdYHR);
cYLR=corr(fwdYHL,fwdYHR);
% cYFL =  -0.9413
% cYFR =   0.0946
% cYLR =  -0.1798
save /home/yuval/epiloc/data/corrY cY* fwdY*


%% Auditory Yuval
cd /home/yuval/Data/marik/yuval/4
load avgSt.mat
% figure; plot(avgSt.avg');
MY=avgSt.avg(:,394);
% figure; topoplot248(MY);
[pow,pntYAud,current]=rimda(MY);
save /home/yuval/epiloc/data/pntYAud pntYAud

MY=avgSt.avg(:,438);


%%
cd /home/yuval/epiloc/data
load pntiAll
load pntSom
%AIR order 
x=num2str(-pnt(pntR,1));y=num2str(pnt(pntR,3));z=num2str(pnt(pntR,2));
[~,w]=afnix(['3dcalc -a ortho+orig -exp "step(9-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix R3'])
afni