%% averaging odball
cd /home/yuval/Data/marik/yuval/4
[events,values]=trigOnset;
trl=events'-203;
trl(:,2)=events'+509;
trl(:,3)=-203;
trl(:,4)=values;
cfg=[];
cfg.dataset='xc,hb,lf_c,rfhp0.1Hz';
cfg.demean='yes';
cfg.blwindow=[-0.2 0];
cfg.bpfilter='yes';
cfg.bpfreq=[1 40];
cfg.padding=0.5;
cfg.channel='MEG';
cfg.trl=trl;
data=ft_preprocessing(cfg);

cfg=[];
cfg.method='summary';
cfg.preproc.hpfilter='yes';
cfg.preproc.hpfreq=60;
datacln=ft_rejectvisual(cfg,data);

%% ica
cfg            = [];
cfg.resamplefs = 300;
cfg.detrend    = 'no';
dummy           = ft_resampledata(cfg, datacln);

% run ica (it takes a long time have a break)
comp_dummy     = ft_componentanalysis([], dummy);

% see the components and find the artifacts
cfgb=[];
cfgb.layout='4D248.lay';
cfgb.channel = {comp_dummy.label{1:5}};
comppic=ft_databrowser(cfgb,comp_dummy);

% run the ICA on the original data
cfg = [];
cfg.topo      = comp_dummy.topo;
cfg.topolabel = comp_dummy.topolabel;
comp     = ft_componentanalysis(cfg, datacln);

% remove the artifact components
cfg = [];
cfg.component = 6; % change
dataica = ft_rejectcomponent(cfg, comp);

cfg=[];
cfg.method='summary';
datafinal=ft_rejectvisual(cfg,datacln);


trl(:,4)=values;
save trl trl

for trli=1:length(datafinal.trial)
    datafinal.trialinfo(trli,1)=trl(trl(:,1)==datafinal.sampleinfo(trli,1),4);
end
save datafinal datafinal

trlSt=find(datafinal.trialinfo==128);
trlOdd=find(datafinal.trialinfo==64);
choice=round(1:(485/99):485);
choice=round(1:(length(trlSt)/length(trlOdd)):length(trlSt));
trlSt=trlSt(choice);
cfg=[];
cfg.trials=trlSt;
avgSt=ft_timelockanalysis(cfg,datafinal);
cfg=[];
cfg.trials=trlOdd;
avgOdd=ft_timelockanalysis(cfg,datafinal);
save avgSt avgSt
save avgOdd avgOdd

%% head model etc

cd /home/yuval/Data/marik/yuval/4
load avgOdd


%% M150
M=avgOdd.avg(:,390);
N=10000;
Pow=zeros(length(gain),1);
tic
for permi=1:N
    Ran=[];
    
    [~,ran]=sort(rand(1,length(gain)/2));
    selected=ran(1:10);
    Ran=[Ran;selected];
    
    srcPerm=false(1,length(gain)/2);
    srcPerm(Ran)=true;
    Gain=gain(:,[srcPerm,srcPerm]);
    source=Gain\M;
    recon=Gain*source;
    R=corr(recon,M).^100;
    pow=zeros(size(Pow));
    pow([srcPerm,srcPerm])=source*R;
    Pow=Pow+pow;
    prog(permi)
end
toc
Pow2=sqrt(Pow(1:length(pnt)).^2+Pow(length(pnt)+1:length(pnt)*2).^2);
figure;
scatter3pnt(pnt,25,Pow2)
hold on
scatter3pnt(hs,5,'k')

%% somato
cd /home/yuval/Data/marik/yuval/som
fn='xc,hb,lf_c,rfhp0.1Hz';
cond={'handR','handL','footL'};
for runi=1
    %cd (num2str(runi))
%     trig=readTrig_BIU(fn);
%     trig=clearTrig(trig);
    [evt,val]=trigOnset;%(trig);
    
    trl=evt'-103;
    trl(:,2)=trl+410;
    trl(:,3)=-103;
    trl(:,4)=val;
    cfg.trl=trl;
    cfg.dataset=fn;
    cfg.channel='MEG';
    cfg.bpfilter='yes';
    cfg.bpfreq=[1 20];
    data=ft_preprocessing(cfg);
    %data.trialinfo=trig(evt)';
    %good=badTrials(data);
    for condi=1:3
        cfg=[];
        cfg.trials=find(data.trialinfo==condi*2);
        avg=ft_timelockanalysis(cfg,data);
        eval(['avg',num2str(runi),'_',cond{condi},'=correctBL(avg,[-0.1 0]);'])
    end
end
clear avg
save avgFilt avg*
%% dipoles
cd /home/yuval/Data/marik/yuval/3
load avgFilt
LRF=[157,150,214];

hs=ft_read_headshape('hs_file');
hs=hs.pnt*1000;
nPNT=642;
inward=[10 20 30];
[pnt1,lf1,vol]=sphereGrid(nPNT,inward(1)); % inner
pnt_complex1=findPerpendicular(pnt1);
[pnt2,lf2]=sphereGrid(nPNT,inward(2)); % outer
pnt_complex2=findPerpendicular(pnt2);
[pnt3,lf3,vol]=sphereGrid(nPNT,inward(3)); % inner
pnt_complex3=findPerpendicular(pnt3);
pnt=[pnt1(:,:,1);pnt2(:,:,1);pnt3(:,:,1)];
ini=true(length(pnt),1);
ini(pnt(:,3)<vol.o(3))=false; % keep only top half sphere
ini(pnt(:,2)<min(hs(:,2)'+10))=false;
ini(pnt(:,2)>max(hs(:,2)'-10))=false;
ini=find(ini);
pnt=pnt(ini,:);

pnt_complex=pnt_complex1;
pnt_complex(nPNT+1:nPNT*2,:,:)=pnt_complex2;
pnt_complex(nPNT*2+1:nPNT*3,:,:)=pnt_complex3;
pnt_complex=pnt_complex(ini,:,:);
layer(1:nPNT,1)=1;layer(nPNT+1:nPNT*2)=2;layer(nPNT*2+1:nPNT*3)=3;
layer=layer(ini);

LF1=lf1; %head pos 1
LF1.leadfield(nPNT+1:nPNT*2)=lf2.leadfield;
LF1.leadfield(nPNT*2+1:nPNT*3)=lf3.leadfield;
LF1.leadfield=LF1.leadfield(ini);
LF1.pos(nPNT+1:nPNT*2,:)=lf2.pos;
LF1.pos(nPNT*2+1:nPNT*3,:)=lf3.pos;
LF1.pos=LF1.pos(ini,:);
LF1.inside(nPNT+1:nPNT*2)=lf2.inside;
LF1.inside(nPNT*2+1:nPNT*3)=lf2.inside;
LF1.inside=LF1.inside(ini);


% close all
% figure
% plot3pnt(hs,'.')
% hold on
% plot3pnt(pnt,'.k')



gain=[];
for pnti=1:length(pnt)
    dip=(pnt_complex(pnti,:,2)-pnt_complex(pnti,:,1))*LF1.leadfield{pnti}';
    gain(1:248,pnti)=dip;
    dip=(pnt_complex(pnti,:,3)-pnt_complex(pnti,:,1))*LF1.leadfield{pnti}';
    gain(1:248,length(pnt)+pnti)=dip;
end



% MF=avg1_footL.avg(:,180);
% MHF=MH+MF;%*2.5;
% topoplot248(MHF(1:248))

%% right hand
M=avg1_handR.avg(:,LRF(1));
figure;topoplot248(M);
N=10000;
Pow=zeros(length(gain),1);
tic
for permi=1:N
    Ran=[];
    
    [~,ran]=sort(rand(1,length(gain)/2));
    selected=ran(1:10);
    Ran=[Ran;selected];
    
    srcPerm=false(1,length(gain)/2);
    srcPerm(Ran)=true;
    Gain=gain(:,[srcPerm,srcPerm]);
    source=Gain\M;
    recon=Gain*source;
    R=corr(recon,M).^100;
    pow=zeros(size(Pow));
    pow([srcPerm,srcPerm])=source*R;
    Pow=Pow+pow;
    prog(permi)
end
toc
Pow1=sqrt(Pow(1:length(pnt)).^2+Pow(length(pnt)+1:length(pnt)*2).^2);
figure;
scatter3pnt(pnt,25,Pow1)
hold on
scatter3pnt(hs,5,'k')

M=avg1_handL.avg(:,LRF(2));
figure;topoplot248(M);
N=10000;
Pow=zeros(length(gain),1);
tic
for permi=1:N
    Ran=[];
    
    [~,ran]=sort(rand(1,length(gain)/2));
    selected=ran(1:10);
    Ran=[Ran;selected];
    
    srcPerm=false(1,length(gain)/2);
    srcPerm(Ran)=true;
    Gain=gain(:,[srcPerm,srcPerm]);
    source=Gain\M;
    recon=Gain*source;
    R=corr(recon,M).^100;
    pow=zeros(size(Pow));
    pow([srcPerm,srcPerm])=source*R;
    Pow=Pow+pow;
    prog(permi)
end
toc
Pow1=sqrt(Pow(1:length(pnt)).^2+Pow(length(pnt)+1:length(pnt)*2).^2);
figure;
scatter3pnt(pnt,25,Pow1)
hold on
scatter3pnt(hs,5,'k')

M=avg1_footL.avg(:,LRF(3));
figure;topoplot248(M);
N=10000;
Pow=zeros(length(gain),1);
tic
for permi=1:N
    Ran=[];
    
    [~,ran]=sort(rand(1,length(gain)/2));
    selected=ran(1:10);
    Ran=[Ran;selected];
    
    srcPerm=false(1,length(gain)/2);
    srcPerm(Ran)=true;
    Gain=gain(:,[srcPerm,srcPerm]);
    source=Gain\M;
    recon=Gain*source;
    R=corr(recon,M).^100;
    pow=zeros(size(Pow));
    pow([srcPerm,srcPerm])=source*R;
    Pow=Pow+pow;
    prog(permi)
end
toc
Pow1=sqrt(Pow(1:length(pnt)).^2+Pow(length(pnt)+1:length(pnt)*2).^2);
figure;
scatter3pnt(pnt,25,Pow1)
hold on
scatter3pnt(hs,5,'k')

%% two hands
M=avg1_handR.avg(:,LRF(1))+avg1_handL.avg(:,LRF(2));
figure;
topoplot248(M);
N=10000;
Pow=zeros(length(gain),1);
tic
for permi=1:N
    Ran=[];
    
    [~,ran]=sort(rand(1,length(gain)/2));
    selected=ran(1:10);
    Ran=[Ran;selected];
    
    srcPerm=false(1,length(gain)/2);
    srcPerm(Ran)=true;
    Gain=gain(:,[srcPerm,srcPerm]);
    source=Gain\M;
    recon=Gain*source;
    R=corr(recon,M).^100;
    pow=zeros(size(Pow));
    pow([srcPerm,srcPerm])=source*R;
    Pow=Pow+pow;
    prog(permi)
end
toc
Pow1=sqrt(Pow(1:length(pnt)).^2+Pow(length(pnt)+1:length(pnt)*2).^2);
figure;
scatter3pnt(pnt,25,Pow1)
hold on
scatter3pnt(hs,5,'k')


%% aud right source?
cd /home/yuval/Data/marik/yuval/4
load avgOdd
load gain
load pnt
M=avgOdd.avg(:,390);
N=10000;
Pow=zeros(length(gain),1);
tic
for permi=1:N
    Ran=[];
    
    [~,ran]=sort(rand(1,length(gain)/2));
    selected=ran(1:10);
    Ran=[Ran;selected];
    
    srcPerm=false(1,length(gain)/2);
    srcPerm(Ran)=true;
    Gain=gain(:,[srcPerm,srcPerm]);
    source=Gain\M;
    recon=Gain*source;
    R=corr(recon,M).^10000;
    pow=zeros(size(Pow));
    pow([srcPerm,srcPerm])=source*R;
    Pow=Pow+pow;
    prog(permi)
end
toc
Pow2=sqrt(Pow(1:length(pnt)).^2+Pow(length(pnt)+1:length(pnt)*2).^2);
figure;
scatter3pnt(pnt,25,Pow2)
% hold on
% scatter3pnt(hs,5,'k')

%% median
Pow=zeros(length(gain),10);
for NN=1:10
    N=10000;
    
    disp(' ')
    disp(['round ',num2str(NN)])
    disp(' ')
    for permi=1:N
        Ran=[];
        [~,ran]=sort(rand(1,length(gain)/2));
        selected=ran(1:10);
        Ran=[Ran;selected];
        srcPerm=false(1,length(gain)/2);
        srcPerm(Ran)=true;
        Gain=gain(:,[srcPerm,srcPerm]);
        source=Gain\M;
        recon=Gain*source;
        R=corr(recon,M).^10000;
        pow=zeros(size(Pow,1),1);
        pow([srcPerm,srcPerm])=source*R;
        Pow(1:length(pow),NN)=Pow(1:length(pow),NN)+pow;
        prog(permi)
    end
end
Pow1=median(Pow,2);
Pow2=sqrt(Pow1(1:length(pnt)).^2+Pow1(length(pnt)+1:length(pnt)*2).^2);
figure;
scatter3pnt(pnt,25,Pow2)


% hold on
% scatter3pnt(hs,5,'k')

%% median for somatosensory
cd /home/yuval/Data/marik/yuval/3
load avgFilt
LRF=[157,150,214];
load gain
load pnt
M=avg1_handR.avg(:,LRF(1))+avg1_handL.avg(:,LRF(2));

Pow=zeros(length(gain),10);
for NN=1:10
    N=10000;
    
    disp(' ')
    disp(['round ',num2str(NN)])
    disp(' ')
    for permi=1:N
        Ran=[];
        [~,ran]=sort(rand(1,length(gain)/2));
        selected=ran(1:10);
        Ran=[Ran;selected];
        srcPerm=false(1,length(gain)/2);
        srcPerm(Ran)=true;
        Gain=gain(:,[srcPerm,srcPerm]);
        source=Gain\M;
        recon=Gain*source;
        R=corr(recon,M).^1000;
        pow=zeros(size(Pow,1),1);
        pow([srcPerm,srcPerm])=source*R;
        Pow(1:length(pow),NN)=Pow(1:length(pow),NN)+pow;
        prog(permi)
    end
end
Pow1=median(Pow,2);
Pow2=sqrt(Pow1(1:length(pnt)).^2+Pow1(length(pnt)+1:length(pnt)*2).^2);
figure;
scatter3pnt(pnt,25,Pow2)

%% forward 

fwd=gain*Pow1;
figure; topoplot248(fwd);

[~,bigi]=sort(Pow2,'descend')
maxPow=zeros(size(Pow2));
maxPow(bigi(1:2))=1;
figure;
scatter3pnt(pnt,25,maxPow)
