function results=marikVirtual29(Ndip, noiseFactor)

% simulate best fit, more than 2 dipoles

%cd /home/yuval/Data/marik/som2/talk
load pnt
load gain1
load layer
N=10000;
% Ndip=1;
%% simulate
inputLabel='simulated points index';
resultsLabel={'N dipoles found','N correct location','distance'};
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
    %% solve with RIMDA
    %tic
    Rmax=0;
    pntMaxi=zeros(1,10);
    pntMaxi=zeros(1,10);
    allThere=false;
    Pow10=zeros(length(gain),10);
    for medi=1:10
        Pow=zeros(length(gain),1);
        for permi=1:N
            Ran=[];
            [~,ran]=sort(rand(1,length(gain)/2));
            selected=ran(1:10);
            Ran=[Ran;selected];
            srcPerm=false(1,length(gain)/2);
            srcPerm(Ran)=true;
            Gain=gain(:,[srcPerm,srcPerm]);
            source=Gain\Mrand;
            recon=Gain*source;
            R=corr(recon,Mrand).^100;
            pow=zeros(size(Pow));
            pow([srcPerm,srcPerm])=source*R;
            Pow=Pow+pow;
            %prog(permi)
            if R^0.02>Rmax
                Rmax=R^0.02; % we want R^2
                pntMaxi=Ran;
                powMax=zeros(size(Pow));
                powMax([srcPerm,srcPerm])=source;
            end
            if sum(ismember(Ran,input(dermi,:)))==Ndip
                allThere=true;
            end
        end
        Pow10(:,medi)=Pow;
    end
    PowMed=median(Pow10,2);
    PowAvg=mean(Pow10,2);
    %toc
    %PowSim=sqrt(Pow(1:920).^2+Pow(921:1840).^2);
    %figure;scatter3pnt(pnt,25,PowSim)
    [~,~,pntMed,~]=getCurrent(PowMed,pnt,Mrand,gain,30,0.3,false); % FIXME - check ori error?
    %PowSimR=sqrt(powMax(1:920).^2+powMax(921:1840).^2);
    [~,~,pntAvg,~]=getCurrent(PowAvg,pnt,Mrand,gain,30,0.3,false); % FIXME - check ori error?
    [~,~,pntR,~]=getCurrent(powMax,pnt,Mrand,gain,30,0.3,false); %
%     results(dermi,1:6)=[length(pntMed),sum(ismember(pntMed,input(dermi,:))),...
%         length(pntAvg),sum(ismember(pntAvg,input(dermi,:))),...
%         length(pntR),sum(ismember(pntinvR,input(dermi,:)))];
    
    results(dermi,1:4)=getErrors(pnt,pnti,pntMed,layer);
    results(dermi,5:8)=getErrors(pnt,pnti,pntAvg,layer);
    results(dermi,9:12)=getErrors(pnt,pnti,pntR,layer);
    %     resultsR(dermi,1)=length(pntinvR);
    %     resultsR(dermi,2)=sum(ismember(pntMaxi,input(dermi,:)));
    %     resultsR(dermi,3)=Rmax;
    %     resultsR(dermi,4)=allThere;
    %
    %
    %     goodIteration(dermi)=allThere;
    %
    %     %%
    %     results(dermi,1)=length(pntinv); % should be 2 dipoles found
    %     [goodfit,dipi]=ismember(pntinv,pnti);
    %     results(dermi,2)=sum(goodfit); % should be two dipoles in exact location
    %     results(dermi,3)=sqrt(sum((pnt(pnti(1),:)-pnt(pnti(2),:)).^2)); %distance
    %     %distnc=sqrt(sum((pnt-repmat(pnt(Ran,:),length(pnt),1)).^2,2));
    %
    prog(dermi)
end

save(['results_',num2str(Ndip),'_',num2str(noiseFactor)],'results','input')
disp('done');
marikVirtual29plot(input,pnt,results);
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