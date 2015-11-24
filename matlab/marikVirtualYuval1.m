
%% aud right source
cd /home/yuval/Data/marik/yuval/4
load avgOdd
load gain
load pnt
M=avgOdd.avg(:,390);

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


%% forward 

fwd=gain*Pow1;
figure; topoplot248(fwd);

[~,bigi]=sort(Pow2,'descend')
maxPow=zeros(size(Pow2));
maxPow(bigi(1:2))=1;
figure;
scatter3pnt(pnt,25,maxPow)
