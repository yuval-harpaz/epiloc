% sequential noise norm, no rand
cd /home/yuval/Data/marik/som2/talk
load pnt
load gain1
load layer

%% hand R + foot

load RF
load Pow4

%[current,ori,pnti,PowClust]=getCurrent(Pow,pnt,RF,gain,30,0.3);



% Ninv=100000;
Ninv=100; %10000

PowN=zeros(length(gain),1);
Pow1=zeros(length(gain),1);
Pow10=zeros(length(gain),10);
tic
for fwdi=1:length(gain)
    Mrand=zeros(248,1); % not random really
    GainFwd=gain(:,fwdi);
    Mrand=Mrand+GainFwd;
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

figure;scatter3pnt(pnt,25,Pow4./PowRand)


load('/home/yuval/Data/marik/som2/talk/noise/bias_f1000_i100_d1.mat')



%[current,ori,pnti,Pow3]=getCurrent(Pow./[PowRand;PowRand],pnt,RF,gain,30,0.4);
[current,ori,pnti,PowClust]=getCurrent(Pow,pnt,RF,gain,30,0.4);
mask=PowClust>0;
figure;scatter3pnt(pnt,25,Pow4.*mask./PowRand)
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