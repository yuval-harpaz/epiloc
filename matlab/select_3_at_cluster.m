nPNT=642;
try
    cd /home/yuval/Data/marik/som2/1
catch err
    cd /home/oshrit/MyDocuments/DATA/som2/1
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

gain1=[];
for pnti=1:length(pnt)
    dip=(pnt_complex(pnti,:,2)-pnt_complex(pnti,:,1))*LF1.leadfield{pnti}';
    gain1(1:248,pnti)=dip;
    dip=(pnt_complex(pnti,:,3)-pnt_complex(pnti,:,1))*LF1.leadfield{pnti}';
    gain1(1:248,length(pnt)+pnti)=dip;
end

try 
    cd /home/yuval/Data/marik/som2
catch err
    cd /home/oshrit/MyDocuments/DATA/som2
end
load avgFilt avg1_footL avg1_handR

MH=avg1_handR.avg(:,138);

%% hand
Pow=zeros(length(gain1),1);
tic
for permi=1:100000
    Ran=[];
    
    [~,ran]=sort(rand(1,length(gain1)/2));
    selected=ran(1:10);
    Ran=[Ran;selected];
    
    srcPerm=false(1,length(gain1)/2);
    srcPerm(Ran)=true;
    Gain=gain1(:,[srcPerm,srcPerm]);
    source=Gain\MH;
    recon=Gain*source;
    R=corr(recon,MH).^100;
    pow=zeros(size(Pow));
    pow([srcPerm,srcPerm])=source*R;
    Pow=Pow+pow;
    prog(permi)
end
toc
Pow1=sqrt(Pow(1:length(gain1)/2).^2+Pow((length(gain1)/2+1):length(gain1)).^2);
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
if exist('./1/ori.mat','file')
    cd 1
    load ori
    load LF
else
    fwd=mne_read_forward_solution('/home/oshrit/MyDocuments/DATA/som2/MNE/pos1_raw-oct-6-fwd.fif');
    cd 1
    cfg=[];
    cfg.fileName='/home/oshrit/MyDocuments/DATA/som2/MNE/pos1_raw-oct-6-fwd.fif'; %'/home/yuval/Data/marik/som2/MNE/pos1_raw-oct-6-fwd.fif';
    cfg.subset=false;
    [Grid,ftLeadfield,ori]=fs2grid(cfg);
    save ori ori
    LF=ft_convert_units(ftLeadfield,'mm');
    save LF LF
end

distnc=sqrt(sum((LF.pos-repmat(pos,length(LF.pos),1)).^2,2));
srci=find(distnc<20);
srcPos=LF.pos(srci,:);
gain2=[];
for si=1:length(srci)
    gain2(1:length(MH),si)=LF.leadfield{srci(si)}*ori(srci(si),:)';
end

srcN=size(gain2,2);
Pow=zeros(srcN,1);
tic
for permi=1:100000
    Ran=[];
    [~,Ran]=sort(rand(1,srcN));
    Ran=Ran(1:3);
    srcPerm=false(1,srcN);
    srcPerm(Ran)=true;
    Gain=gain2(:,srcPerm);
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

mm=gain2'*MH;
figure;
scatter3pnt(srcPos,10,abs(mm))
hold on
scatter3pnt(hs,1,'k')


%%
srcN=size(gain2,2);
Pow1=zeros(length(gain1)+srcN,1);
tic
for permi=1:100000
    
    [~,ran]=sort(rand(1,length(gain1)/2));
    Ran=ran(1:7);   
    srcPerm1=false(1,length(gain1)/2);
    srcPerm1(Ran)=true;
    Gain=gain1(:,[srcPerm1,srcPerm1]);
    
    [~,Ran]=sort(rand(1,srcN));
    Ran=Ran(1:3);
    srcPerm2=false(1,srcN);
    srcPerm2(Ran)=true;
    Gain=[Gain, gain2(:,srcPerm2)];
    
    source=Gain\MH;
    recon=Gain*source;
    R=corr(recon,MH).^100;
       
    pow=zeros(size(Pow1,1),1);
    pow([srcPerm1,srcPerm1,srcPerm2])=source*R;
    Pow1=Pow1+pow;    
        
    prog(permi)
end
toc
Pow11=sqrt(Pow1(1:length(gain1)/2).^2+Pow1((length(gain1)/2+1):length(gain1)).^2);
Pow2=[Pow11; Pow1(length(gain1)+1:end)];

figure;
scatter3pnt(pnt,25,Pow11)
[~,maxPNT]=max(Pow11);
hold on
scatter3(pnt(maxPNT,1),pnt(maxPNT,2),pnt(maxPNT,3),30,0)

[~,sortPNT]=sort(Pow11,'descend');
sortPNT(1:5)
    
% figure;
scatter3pnt(srcPos,10,abs(Pow1(length(gain1)+1:end)))
[~,maxPNT]=max(Pow1(length(gain1)+1:end));
hold on
scatter3pnt(hs,1,'k')

% caxis([0 0.2]*10^-7)
% caxis([0 0.05]*10^-7)
% caxis([0 0.005]*10^-7)
% caxis([0 0.0005]*10^-7)
% caxis([0 0.00005]*10^-7)
% caxis([0 0.000005]*10^-7)
% caxis([0 0.00004]*10^-7)
% caxis([0 0.00002]*10^-7)

% mm=gain2'*MH;
% figure;
% scatter3pnt(srcPos,10,abs(mm))
% hold on
% scatter3pnt(hs,1,'k')

