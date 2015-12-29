%% averaging odball
cd /home/yuval/Data/marik/yuval/som



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

%% one hand
cd /home/yuval/Data/marik/yuval/som
load gain
load pnt
load avgFilt
hs=ft_read_headshape('hs_file')
hs=hs.pnt*1000;
RLF=[152,157,300];
M=avg1_handR.avg(:,RLF(1));
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

M=avg1_handL.avg(:,RLF(2));
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

M=avg1_footL.avg(:,RLF(3));
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
M=avg1_handR.avg(:,RLF(1))+avg1_handL.avg(:,RLF(2));
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


%% two hands + foot
M=avg1_handR.avg(:,RLF(1))+avg1_handL.avg(:,RLF(2))+avg1_footL.avg(:,RLF(3));
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


