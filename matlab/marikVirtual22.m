
cd /home/yuval/Data/marik/som2/talk
load pnt
load gain1
load layer
% hand R + foot
load RF
load Pow4
% figure;
% scatter3pnt(pnt,25,Pow4)
[current,ori,pnti,PowClust]=getCurrent(Pow,pnt,RF,gain,30,0.3);

%% depth normalize remote sources


% Ninv=100000;
Ninv=100; %10000
Nfwd=1000;
Ndip=3;

PowN=zeros(length(gain),1);
Pow1=zeros(length(gain),1);
Pow10=zeros(length(gain),10);
tic
for fwdi=1:Nfwd
    Mrand=zeros(248,1);
    usedDips=false(length(gain)/2,1);
    dips=[];
    for dipi=1:Ndip
        ran=rand(1,length(gain)/2);
        ran(usedDips)=0;
        [~,Ran]=max(ran);
        %Ran=ran(1);
        srcPerm=false(1,length(gain)/2);
        srcPerm(Ran)=true;
        GainFwd=gain(:,[srcPerm,srcPerm]);
        randOri=rand(1);
        randOri(2,1)=sqrt(1-randOri^2);
        randOri=randOri.*((rand(2,1)<0.5)-0.5)*2;
        Mrand=Mrand+GainFwd*randOri;
        distnc=sqrt(sum((pnt-repmat(pnt(Ran,:),length(pnt),1)).^2,2));
        usedDips(distnc<30)=true;
        dips(dipi)=Ran;
    end
    for NN=1:50   
%         disp(' ')
%         disp(['round ',num2str(NN)])
%         disp(' ')
        for invi=1:Ninv
            Ran=[];
            [~,ran]=sort(rand(1,length(gain)/2));
            Ran=ran(1:10);
            srcPerm=false(1,length(gain)/2);
            srcPerm(Ran)=true;
            Gain=gain(:,[srcPerm,srcPerm]);
            source=Gain\Mrand;
            recon=Gain*source;
            R=corr(recon,Mrand).^100;
            pow=zeros(size(PowN));
            pow([srcPerm,srcPerm])=source*R;
            PowN=PowN+pow;
        end
        Pow10(1:length(pow),NN)=PowN;
    end
    Pow1=Pow1+median(Pow10,2);
    prog(fwdi)
end
toc   
PowRand=sqrt(Pow1(1:920).^2+Pow1(921:1840).^2);
%save (['noise/bias_f',num2str(Nfwd),'_i',num2str(Ninv),'_d',num2str(Ndip)],'PowRand')
figure;
scatter3pnt(pnt,25,PowRand)
figure;
scatter3pnt(pnt,25,Pow4./PowRand)



%load ../bias_f10000_i1000_d33.mat
PowRand_s=zeros(size(pnt,1),1);
Th=mean(PowRand)+3*std(PowRand);
Th_inv=mean(1./PowRand)+3*std(1./PowRand);
% diff_neig=[];
for ly=1:3
    index=find(layer==ly);
    for v=index'
        i_neig=[];
        for p=index'
            if sqrt((pnt(v,1)-pnt(p,1))^2+(pnt(v,2)-pnt(p,2))^2+(pnt(v,3)-pnt(p,3))^2)<20
                p_n=p;
                i_neig=[i_neig,p_n]; % i_neig include the point itself (as the distance to itself is zero)
            end
        end
        %             i_neig_t=i_neig(PowRand(i_neig)<=Th & ((1./PowRand(i_neig))<=Th_inv));
        PowRand_s(v)=median(PowRand(i_neig));
    end
end
figure;
scatter3pnt(pnt,25,PowRand_s)
figure;
scatter3pnt(pnt,25,1./PowRand_s)
figure;
scatter3pnt(pnt,25,Pow4./PowRand_s)
[current,ori,pnti,Pow3]=getCurrent(Pow./[PowRand_s;PowRand_s],pnt,RF,gain,30,0.4);
[current,ori,pnti,Pow3]=getCurrent(Pow,pnt,RF,gain,30,0.4);
% figure;
% scatter3pnt(pnt,25,Pow3)
% cd /home/yuval/Data/marik/som2/talk
cd /home/yuval/Data/marik/som2/1
pnt2vol(pnt,Pow3,'RF');


% [~,maxPNT]=max(Pow4);
% hold on
% scatter3(pnt(maxPNT,1),pnt(maxPNT,2),pnt(maxPNT,3),30,0)
% figure;
% scatter3pnt(pnt,25,Pow4)
% 
% save talk/RF RF
% save talk/Pow4 Pow4

%% hand R
load MHR
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
[current,ori,pnti,Pow3]=getCurrent(Pow,pnt,MHR,gain,30,0.3)
pnt2vol(pnt,Pow3,'MHR');

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
pnt2vol(pnt,Pow4,'RF');

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