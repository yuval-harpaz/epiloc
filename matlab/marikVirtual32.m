%% plot figures for paper

%% Somatosensory Marik
load /home/yuval/epiloc/data/MHR
cd /home/yuval/Data/marik/som2/1
[pow,pntR,current,fwdHR]=rimda(MHR);
load /home/yuval/epiloc/data/MHL
[pow,pntL,current,fwdHL]=rimda(MHL);
load /home/yuval/epiloc/data/MF
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
MY=avgSt.avg(:,394);
% figure; topoplot248(MY);
[pow,pntYAud,current]=rimda(MY);
save /home/yuval/epiloc/data/pntYAud pntYAud
%MY=avgSt.avg(:,438);


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
for ii=1:2
    x=num2str(-pnt(pntYAud(ii),1));y=num2str(pnt(pntYAud(ii),3));z=num2str(pnt(pntYAud(ii),2));
    [~,w]=afnix(['3dcalc -a ortho+orig -exp "step(25-(x-',x,')*(x-',x,')-(y-',y,')*(y-',y,')-(z-',z,')*(z-',z,'))" -prefix YAud',num2str(ii)])
end

%% 
% 3dcalc -a R+orig -b F+orig -c L+orig -d RFL+orig -exp "a+b+c+2*d" -prefix RFLRFL

cd /home/yuval/epiloc/data/
% cd /home/oshrit/MyDocuments/DATA/Marik/epiloc/data
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
caxis([-2.6*10^-13 4.3*10^-13])
figure; topoplot248(MYHL);
caxis([-2.6*10^-13 4.3*10^-13])
figure; topoplot248(MYF);
caxis([-2.6*10^-13 4.3*10^-13])
figure; topoplot248(MYRFL);
caxis([-2.6*10^-13 4.3*10^-13])

cd /home/yuval/Data/marik/marikAud
load avgSt.mat
% load /home/oshrit/MyDocuments/DATA/Marik/epiloc/data/In_Afni/Marik/Aud/avgSt.mat
M=avgSt.avg(:,353);
figure; topoplot248(M);
caxis([-1.2*10^-13 1.0*10^-13])

cd /home/yuval/Data/marik/yuval/4
load avgSt.mat
% load /home/oshrit/MyDocuments/DATA/Marik/epiloc/data/In_Afni/Yuval/Aud/avgSt.mat
MY=avgSt.avg(:,394);
figure; topoplot248(MY);
caxis([-1.*10^-13 1.2*10^-13])

hs=ft_read_headshape('/home/oshrit/MyDocuments/DATA/Marik/epiloc/data/yuvSom/hs_file');
pnti=726;
figure;
plot3pnt(hs,'.k');
hold on
scatter3pnt(pnt(pnti,:),25,[1 0 0])
scatter3pnt(pnt,1,[1 1 1])
colorbar off
rotate3d on