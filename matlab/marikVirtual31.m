function results=marikVirtual31(Ndip, noiseFactor)

% simulate sequential dipole fitting

%cd /home/yuval/Data/marik/som2/talk
load pnt
load gain1
load layer
N=10000;
% Ndip=1;
%% simulate
results=[];
input=[];
% noiseFactor=0.3;
for dermi=1:1000
    %% make simulated field
    
    ran=rand(1,length(gain)/2);
    [~,Ran]=sort(ran,'descend');
    Ran=Ran(1:Ndip);
    input(dermi,1:Ndip)=Ran;
    pnti=Ran;
    Mrand=zeros(248,1);
    for dipi=1:Ndip
        %Ran=ran(1);
        srcPerm=false(1,length(gain)/2);
        srcPerm(Ran(dipi))=true;
        GainFwd=gain(:,[srcPerm,srcPerm]);
        randOri=rand(1);
        randOri(2,1)=sqrt(1-randOri^2);
        randOri=randOri.*((rand(2,1)<0.5)-0.5)*2;
        Mrand=Mrand+GainFwd*randOri;
    end
    noise=randn(248,1);
    noise=noise./std(noise)*noiseFactor/3;
    interf=gain*randn(length(layer)*2,1);
    interf=interf./std(interf)*(2*noiseFactor/3);
    Mrand=Mrand./std(Mrand);
%     noise=noise./max(abs(noise)).*max(abs(Mrand'))*(noiseFactor/3); 
%     interf=interf./max(abs(interf)).*max(abs(Mrand'))*(2*noiseFactor/3);
    Mrand=Mrand+noise+interf;
    %     Pow3=zeros(size(PowSim));
    %     Pow3(pnti)=1;
    %     figure;scatter3pnt(pnt,25,Pow3)
    %tic
    Pow=zeros(length(gain),1);
    Mresid=Mrand;
    
    pntAll=[];
    for dipi=1:5
        Rmax=0;
        rp=[];
        for loci=1:length(Pow)./2;
            srcPerm=false(1,length(gain)/2);
            srcPerm(loci)=true;
            Gain=gain(:,[srcPerm,srcPerm]);
            source=Gain\Mresid;
            recon=Gain*source;
            R=corr(recon,Mresid).^2;
            rp(loci,2:3)=source;
            if R>Rmax
                Rmax=R; % we want R^2
                pntMaxi=loci;
                powMax=zeros(size(Pow));
                powMax([srcPerm,srcPerm])=source;
                reconMax=recon;
            end
        end
        pntAll=[pntAll,pntMaxi];
        Mresid=Mresid-reconMax;
        Pow=Pow+powMax;
    end
    [~,~,pntSeq,~]=getCurrent(Pow,pnt,Mrand,gain,30,0.3,false); % FIXME - check ori error?
    results(dermi,1:4)=getErrors(pnt,pnti,pntSeq,layer);
    prog(dermi)
end

save(['resultsSeq_',num2str(Ndip),'_',num2str(noiseFactor),'.mat'],'results','input')
disp('done');
%marikVirtual29plot(input,pnt,results);
function errVec=getErrors(pnt,pnti,pntX,layer)
errVec=[];
distances=[];
errVec(1)=length(pntX);
errVec(1,2)=sum(ismember(pntX,pnti));
% get all distances
for Xi=1:length(pntX)
    distances(Xi,1:length(pnti))=sqrt(sum(power(pnt(pnti,:)-repmat(pnt(pntX(Xi),:),length(pnti),1),2),2));
end
% find pairs of nearby sources and sum the error
distSum=0;
depthErr=0;
for pairi=1:min([length(pnti),length(pntX)])
    [minD,minDi]=sort(distances(:));
    x=find(sum(distances==minD(1),1));
    x=x(1); % when, say, there are two zero errors
    y=find(distances(:,x)==minD(1));
    y=y(1);
    %y=find(sum(distances==minD(1),2));
    distSum=distSum+distances(y,x);
    if length(pnti)>length(pntX) % more rows
        distances(:,x)=max(minD(end));
    else
        distances(y,:)=max(minD(end));
    end
    depthErr=depthErr+layer(pnti(x))-layer(pntX(y));
end
depthErr=depthErr./pairi;
errVec(1,3)=distSum./pairi;
errVec(1,4)=depthErr;
        

