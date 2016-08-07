
%% aud right source
cd /home/yuval/Data/marik/yuval/4
load avgOdd
load gain
load pnt
M=avgOdd.avg(:,390);

%% median
Pow=zeros(length(gain),10);
N=10000;
for NN=1:10
disp(' ')
    disp(['round ',num2str(NN)])
    disp(' ')
    for permi=1:N
        Ran=[];
        [~,ran]=sort(rand(1,length(gain)/2));
        selected=ran(1:100);
        Ran=[Ran;selected];
        srcPerm=false(1,length(gain)/2);
        srcPerm(Ran)=true;
        Gain=gain(:,[srcPerm,srcPerm]);
        source=Gain\M;
        %source=pinv(Gain)*M;
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


% hold on
% scatter3pnt(hs,5,'k')


%% forward 

fwd=gain*Pow1;
figure; topoplot248(fwd);

[~,bigi]=sort(Pow2,'descend')
maxPow=zeros(size(Pow2));
maxPow(bigi(1:2))=1;
figure;
scatter3pnt(pnt,25,maxPow)

gain1=gain(:,bigi(1))*Pow1(bigi(1))+gain(:,bigi(1)+length(Pow2))*Pow1(bigi(1)+length(Pow2));
figure;topoplot248(gain1);
gain2=gain(:,bigi(2))*Pow1(bigi(2))+gain(:,bigi(2)+length(Pow2))*Pow1(bigi(2)+length(Pow2));
% FIXME - normalise gain1 and gain2
figure;topoplot248(gain2);
source=[gain1,gain2]\M;
fwd=[gain1,gain2]*source;
figure; topoplot248(fwd);

