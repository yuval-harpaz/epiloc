
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

%% hand R
load avgFilt avg1_handR
MHR=avg1_handR.avg(:,138);
% MF=avg1_footL.avg(:,180);
% MHF=MH+MF;%*2.5;
% topoplot248(MHF(1:248))
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
    source=Gain\MHR;
    recon=Gain*source;
    R=corr(recon,MHR).^100;
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

save talk/MHR MHR
save talk/Pow1 Pow1
%% hand L
load avgFilt avg1_handL
MHL=avg1_handL.avg(:,138);
% MF=avg1_footL.avg(:,180);
% MHF=MH+MF;%*2.5;
% topoplot248(MHF(1:248))
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
    source=Gain\MHL;
    recon=Gain*source;
    R=corr(recon,MHL).^100;
    pow=zeros(size(Pow));
    pow([srcPerm,srcPerm])=source*R;
    Pow=Pow+pow;
    prog(permi)
end
toc
Pow2=sqrt(Pow(1:920).^2+Pow(921:1840).^2);
figure;
scatter3pnt(pnt,25,Pow2)
[~,maxPNT]=max(Pow2);
hold on
scatter3(pnt(maxPNT,1),pnt(maxPNT,2),pnt(maxPNT,3),30,0)

save talk/MHL MHL
save talk/Pow2 Pow2

%% foot
load avgFilt avg1_footL

MF=avg1_footL.avg(:,180);
% MF=avg1_footL.avg(:,180);
% MHF=MH+MF;%*2.5;
% topoplot248(MHF(1:248))

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
    source=Gain\MF;
    recon=Gain*source;
    R=corr(recon,MF).^100;
    pow=zeros(size(Pow));
    pow([srcPerm,srcPerm])=source*R;
    Pow=Pow+pow;
    prog(permi)
end
toc
Pow3=sqrt(Pow(1:920).^2+Pow(921:1840).^2);
figure;
scatter3pnt(pnt,25,Pow3)
[~,maxPNT]=max(Pow3);
hold on
scatter3(pnt(maxPNT,1),pnt(maxPNT,2),pnt(maxPNT,3),30,0)
figure;
scatter3pnt(pnt,25,Pow3)

save talk/MF MF
save talk/Pow3 Pow3

%% hand R + foot

RF=MHR+MF;
% MF=avg1_footL.avg(:,180);
% MHF=MH+MF;%*2.5;
% topoplot248(MHF(1:248))

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
    source=Gain\RF;
    recon=Gain*source;
    R=corr(recon,RF).^100;
    pow=zeros(size(Pow));
    pow([srcPerm,srcPerm])=source*R;
    Pow=Pow+pow;
    prog(permi)
end
toc
Pow4=sqrt(Pow(1:920).^2+Pow(921:1840).^2);
figure;
scatter3pnt(pnt,25,Pow4)
[~,maxPNT]=max(Pow4);
hold on
scatter3(pnt(maxPNT,1),pnt(maxPNT,2),pnt(maxPNT,3),30,0)
figure;
scatter3pnt(pnt,25,Pow4)

save talk/RF RF
save talk/Pow4 Pow4

%% foot

RFL=RF+MHL;
% MF=avg1_footL.avg(:,180);
% MHF=MH+MF;%*2.5;
% topoplot248(MHF(1:248))

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
    source=Gain\RFL;
    recon=Gain*source;
    R=corr(recon,RFL).^100;
    pow=zeros(size(Pow));
    pow([srcPerm,srcPerm])=source*R;
    Pow=Pow+pow;
    prog(permi)
end
toc
Pow5=sqrt(Pow(1:920).^2+Pow(921:1840).^2);
figure;
scatter3pnt(pnt,25,Pow5)
[~,maxPNT]=max(Pow5);
hold on
scatter3(pnt(maxPNT,1),pnt(maxPNT,2),pnt(maxPNT,3),30,0)
figure;
scatter3pnt(pnt,25,Pow5)

save talk/RFL RFL
save talk/Pow5 Pow5

%% make BRIK
cd /home/yuval/Data/marik/som2/talk
load pnt
load Pow4
load Pow5

cd /home/yuval/Data/marik/som2/1
pnt2vol(pnt,Pow4,'RF');
pnt2vol(pnt,Pow5,'RFL');
pnt2vol(pnt,Pow1,'R');
pnt2vol(pnt,Pow2,'L');
pnt2vol(pnt,Pow3,'F');