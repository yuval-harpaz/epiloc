
nPNT=642;

cd /home/yuval/Copy/MEGdata/alice/idan
hs=ft_read_headshape('hs_file');
hs=hs.pnt*1000;
inward=[10 20 30 40];
pnt=[];
pnt_complex=[];
for shifti=1:length(inward)
    [pntX,lfX,volX]=sphereGrid(nPNT,inward(shifti)); 
    pnt_complexX=findPerpendicular(pntX);
    pnt=[pnt;pntX(:,:,1)];
    if shifti==1
        pnt_complex=pnt_complexX;
        LF1=lfX;
    else
        pnt_complex(nPNT*(shifti-1)+1:nPNT*shifti,:,:)=pnt_complexX;
        LF1.leadfield(nPNT*(shifti-1)+1:nPNT*shifti)=lfX.leadfield;
        LF1.pos(nPNT*(shifti-1)+1:nPNT*shifti,:)=lfX.pos;
        LF1.inside(nPNT*(shifti-1)+1:nPNT*shifti)=lfX.inside;
    end
    layer(nPNT*(shifti-1)+1:nPNT*shifti)=shifti;
end
% [pnt2,lf2]=sphereGrid(nPNT,inward(2)); 
% pnt_complex2=findPerpendicular(pnt2);
% [pnt3,lf3,vol]=sphereGrid(nPNT,inward(3)); 
% pnt_complex3=findPerpendicular(pnt3);
% 
% pnt=[pnt1(:,:,1);pnt2(:,:,1);pnt3(:,:,1)];
ini=true(length(pnt),1);
ini(pnt(:,3)<(volX.o(3)-20))=false; % keep only top half sphere
ini(pnt(:,2)<min(hs(:,2)'+10))=false;
ini(pnt(:,2)>max(hs(:,2)'-10))=false;
ini=find(ini);

pnt=pnt(ini,:);
pnt_complex=pnt_complex(ini,:,:);
layer=layer(ini);
LF1.pos=LF1.pos(ini,:);
LF1.leadfield=LF1.leadfield(ini);
LF1.inside=LF1.inside(ini);

plot3pnt(hs,'.k')
hold on
plot3pnt(pnt,'.b')

gain=[];
for pnti=1:length(pnt)
    dip=(pnt_complex(pnti,:,2)-pnt_complex(pnti,:,1))*LF1.leadfield{pnti}';
    gain(1:248,pnti)=dip;
    dip=(pnt_complex(pnti,:,3)-pnt_complex(pnti,:,1))*LF1.leadfield{pnti}';
    gain(1:248,length(pnt)+pnti)=dip;
end

load /home/yuval/Copy/MEGdata/alice/ga2015/GavgEqTrlsubs Mr1
t1=nearest(Mr1.time,0.09);
t2=nearest(Mr1.time,0.11);
M=mean(Mr1.avg(1:248,t1:t2),2);

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
mid=length(pnt);
Pow1=sqrt(Pow(1:mid).^2+Pow((mid+1):mid.*2).^2);
figure;
scatter3pnt(pnt,25,Pow1)
[~,maxPNT]=max(Pow1);
hold on
scatter3(pnt(maxPNT,1),pnt(maxPNT,2),pnt(maxPNT,3),30,0)

%% depth bias

Ninv=1000;
Nfwd=1000;
Ndip=300;
NdipInv=10;%Ndip; % 10
Pow=zeros(length(gain),1);
tic
for fwdi=1:Nfwd
    Mrand=zeros(248,1);
    for dipi=1:Ndip
        Ran=[];
        [~,Ran]=max(rand(1,length(gain)/2));
        %Ran=ran(1);
        srcPerm=false(1,length(gain)/2);
        srcPerm(Ran)=true;
        GainFwd=gain(:,[srcPerm,srcPerm]);
        randOri=rand(1);
        randOri(2,1)=sqrt(1-randOri^2);
        randOri=randOri.*((rand(2,1)<0.5)-0.5)*2;
        Mrand=Mrand+GainFwd*randOri;
    end
    %disp('done fwd')
    for invi=1:Ninv
        Ran=[];
        [~,ran]=sort(rand(1,length(gain)/2));
        Ran=ran(1:NdipInv);
        srcPerm=false(1,length(gain)/2);
        srcPerm(Ran)=true;
        Gain=gain(:,[srcPerm,srcPerm]);
        source=Gain\Mrand;
        recon=Gain*source;
        R=corr(recon,Mrand).^100;
        pow=zeros(size(Pow));
        pow([srcPerm,srcPerm])=source*R;
        Pow=Pow+pow;
    end
    prog(fwdi)
end
toc


    
PowRand=sqrt(Pow(1:920).^2+Pow(921:1840).^2);
%save (['bias_f',num2str(Nfwd),'_i',num2str(Ninv),'_d',num2str(Ndip)],'PowRand')
%eval(['PowRand',num2str(Ndip),'=PowRand;']);
figure;
scatter3pnt(pnt,25,PowRand)

figure;
scatter3pnt(pnt,25,Pow1./PowRand)

