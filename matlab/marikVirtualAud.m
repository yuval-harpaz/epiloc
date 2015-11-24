%% averaging odball
cd /home/yuval/Data/marik/marikAud
events=trigOnset;
trl=events'-203;
trl(:,2)=events'+509;
trl(:,3)=-203;
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
cfg.component = [1 2]; % change
dataica = ft_rejectcomponent(cfg, comp);

cfg=[];
cfg.method='summary';
datafinal=ft_rejectvisual(cfg,datacln);

save datafinal datafinal

for trli=1:length(datafinal.trial)
    datafinal.trialinfo(trli,1)=trl(trl(:,1)==datafinal.sampleinfo(trli,1),4);
end
trlSt=find(datafinal.trialinfo==128);
trlOdd=find(datafinal.trialinfo==64);
choice=round(1:(485/99):485);
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

cd /home/yuval/Data/marik/marikAud
load avg
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

%% M100
M=avg.avg(:,353);
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

%% M150
M=avg.avg(:,390);
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

M=avg.avg(:,500);
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
Pow3=sqrt(Pow(1:length(pnt)).^2+Pow(length(pnt)+1:length(pnt)*2).^2);
figure;
scatter3pnt(pnt,25,Pow3)
hold on
scatter3pnt(hs,5,'k')