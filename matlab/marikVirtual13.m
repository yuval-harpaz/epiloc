
nPNT=642;
try
    cd /home/yuval/Data/marik/som2/1
catch
    cd /home/oshrit/MyDocuments/DATA/epiloc/data
end
hs=ft_read_headshape('hs_file');
hs=hs.pnt*1000;

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

cd /home/yuval/Data/marik/som2
load avgFilt avg1_handR

MH=avg1_handR.avg(:,138);
% MF=avg1_footL.avg(:,180);
% MHF=MH+MF;%*2.5;
% topoplot248(MHF(1:248))


%% hand
Pow=zeros(length(gain),1);
tic
for permi=1:100000
    Ran=[];
    
    [~,ran]=sort(rand(1,length(gain)/2));
    selected=ran(1:10);
    Ran=[Ran;selected];
    
    srcPerm=false(1,length(gain)/2);
    srcPerm(Ran)=true;
    Gain=gain(:,[srcPerm,srcPerm]);
    source=Gain\MH;
    recon=Gain*source;
    R=corr(recon,MH).^100;
    pow=zeros(size(Pow));
    pow([srcPerm,srcPerm])=source*R;
    Pow=Pow+pow;
    prog(permi)
end
toc
Pow1=sqrt(Pow(1:920).^2+Pow(921:1840).^2);
figure;
scatter3pnt(pnt,25,Pow1)
[~,maxPNT]=max(Pow1);
hold on
scatter3(pnt(maxPNT,1),pnt(maxPNT,2),pnt(maxPNT,3),30,0)

[~,sortPNT]=sort(Pow1,'descend');
sortPNT(1:5)

%% anatomy - Hand
pos=pnt(sortPNT(1),:);
% ori=pnt_complex(sortPNT(1),:,2);
% pnt_complex(
cd 1
if exist('ori.mat','file')
    load ori
    load LF
else
    
    fwd=mne_read_forward_solution('/home/yuval/Data/marik/som2/MNE/pos1_raw-oct-6-fwd.fif');
    cd 1
    cfg=[];
    cfg.fileName='/home/yuval/Data/marik/som2/MNE/pos1_raw-oct-6-fwd.fif';
    cfg.subset=false;
    [Grid,ftLeadfield,ori]=fs2grid(cfg);
    save ori ori
    LF=ft_convert_units(ftLeadfield,'mm')
    save LF LF
end

% figure;
% plot3pnt(hs,'.k')
% hold on
% plot3pnt(LF.pos,'.g')
% plot3pnt(pos,'.r')
% view([-90 90])
distnc=sqrt(sum((LF.pos-repmat(pos,length(LF.pos),1)).^2,2));
srci=find(distnc<20);
srcPos=LF.pos(srci,:);
gain=[];
for si=1:length(srci)
    gain(1:length(MH),si)=LF.leadfield{srci(si)}*ori(srci(si),:)';
end
srcN=size(gain,2);
Pow=zeros(srcN,1);
tic
for permi=1:100000
    Ran=[];
    [~,Ran]=sort(rand(1,srcN));
    Ran=Ran(1:3);
    srcPerm=false(1,srcN);
    srcPerm(Ran)=true;
    Gain=gain(:,srcPerm);
    source=Gain\MH;
    recon=Gain*source;
    R=corr(recon,MH).^100;
    pow=zeros(size(Pow));
    pow(srcPerm)=source*R;
    Pow=Pow+pow;
    prog(permi)
end
toc

figure;
scatter3pnt(srcPos,10,abs(Pow))
[~,maxPNT]=max(Pow1);
hold on
scatter3pnt(hs,1,'k')

mm=abs(gain'*MH);
figure;
scatter3pnt(srcPos,10,mm)
hold on
scatter3pnt(hs,1,'k')



% Pow1=zeros(length(LF.pos),1);
% Pow1(srci)=abs(Pow);
% figure;
% scatter3pnt(LF.pos,25,Pow1)
% [~,maxPNT]=max(Pow1);
% hold on
% scatter3(LF.pos(maxPNT,1),LF.pos(maxPNT,2),LF.pos(maxPNT,3),30,0)
% 
% [~,sortPNT]=sort(Pow1,'descend');
% sortPNT(1:5)
fwd=mne_read_forward_solution('/home/yuval/Data/marik/som2/MNE/pos1_raw-oct-6-fwd.fif');
nL=fwd.src(1).np;
nR=fwd.src(2).np;
src=zeros(nL+nR,1);
src(srci)=mm;
srcMultL(1:nL)=src(1:nL);
srcMultR(1:nR)=src(nL+1:end);
src=zeros(nL+nR,1);
src(srci)=Pow;
srcL(1:nL)=src(1:nL);
srcR(1:nR)=src(nL+1:end);
save src srcL srcR srcMu*







