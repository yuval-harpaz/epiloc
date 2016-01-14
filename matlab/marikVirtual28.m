% simulate best fit, all dipoles + noise added
cd /home/yuval/Data/marik/som2/talk
load pnt
load gain1
load layer
N=10000;
Ndip=4;
%% simulate
inputLabel={'pnt1','pnt2','ori1x','ori1y','ori2x','ori2y'};
resultsLabel={'N dipoles found','N correct location','distance','dist error 1','dist error 2','level error 1','level error 2'};
resultsR=[];
results=[];
input=[];
goodIteration=false;
for dermi=1:1000
    %% make simulated field
    Pow=zeros(length(gain),1);
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
    %% solve with RIMDA
    for permi=1:N
        noise=randn(248,1);
        noise=noise./max(abs(noise)).*max(abs(Mrand'));
        noise=noise*0.3;
        Nrand=Mrand+noise;
        source=gain\Nrand;
        Pow=Pow+source;
        prog(permi)
    end

    PowSim=sqrt(Pow(1:920).^2+Pow(921:1840).^2);
    [~,~,pntinv,~]=getCurrent(Pow,pnt,Mrand,gain,30,0.3,false); % FIXME - check ori error?
    PowSimR=sqrt(powMax(1:920).^2+powMax(921:1840).^2);
    resultsR(dermi,1)=Rmax;
    resultsR(dermi,2)=sum(ismember(pntMaxi,input(dermi,1:2)));
    [~,~,pntinvR,~]=getCurrent(powMax,pnt,Mrand,gain,30,0.3,false); %
    goodIteration(dermi)=bothThere;
    
    %%
    results(dermi,1)=length(pntinv); % should be 2 dipoles found
    [goodfit,dipi]=ismember(pntinv,pnti);
    results(dermi,2)=sum(goodfit); % should be two dipoles in exact location
    results(dermi,3)=sqrt(sum((pnt(pnti(1),:)-pnt(pnti(2),:)).^2)); %distance
    %distnc=sqrt(sum((pnt-repmat(pnt(Ran,:),length(pnt),1)).^2,2));
    if length(pntinv)==2
        if sum(goodfit)==2
            results(dermi,4:7)=0; %no or layer distance error
        elseif sum(goodfit)==1
            gooddip=dipi(find(goodfit));
            %gooddip=find(goodfit);
            baddip=3-gooddip;
            results(dermi,3+gooddip)=0; %distance
            results(dermi,5+gooddip)=0; %layer error
            pntimissed=pnti;
            pntimissed(dipi(goodfit))=[];
            dist=sqrt(sum((pnt(pntinv(baddip),:)-pnt(pntimissed,:)).^2,2));
            results(dermi,3+baddip)=dist;
            results(dermi,5+gooddip)=layer(pntimissed)-layer(baddip); % inside layer is 1. pos diff is bias to the center
        elseif sum(goodfit)==0
            % got to match them
            distA=sqrt(sum((pnt(pntinv(1),:)-pnt(pnti(1),:)).^2,2));
            distB=sqrt(sum((pnt(pntinv(2),:)-pnt(pnti(2),:)).^2,2));
            distC=sqrt(sum((pnt(pntinv(1),:)-pnt(pnti(2),:)).^2,2));
            distD=sqrt(sum((pnt(pntinv(2),:)-pnt(pnti(1),:)).^2,2));
            if distA+distB<distC+distD
                results(dermi,3+1)=distA;
                results(dermi,3+2)=distB;
                results(dermi,5+1)=layer(pnti(1))-layer(pntinv(1));
                results(dermi,5+2)=layer(pnti(2))-layer(pntinv(2));
            else
                results(dermi,3+1)=distD;
                results(dermi,3+2)=distC;
                results(dermi,5+1)=layer(pnti(1))-layer(pntinv(2));
                results(dermi,5+2)=layer(pnti(2))-layer(pntinv(1));
            end
        end
    elseif length(pntinv)==1
        
        if sum(goodfit)==0 % one wrong dipole
            distA=sqrt(sum((pnt(pntinv,:)-pnt(pnti(1),:)).^2,2));
            distC=sqrt(sum((pnt(pntinv,:)-pnt(pnti(2),:)).^2,2));
            if distA<distC
                results(dermi,3+1)=distA;
                results(dermi,3+2)=nan;
                results(dermi,5+1)=layer(pnti(1))-layer(pntinv);
                results(dermi,5+2)=nan;
            else
                results(dermi,3+2)=distC;
                results(dermi,3+1)=nan;
                results(dermi,5+2)=layer(pnti(2))-layer(pntinv);
                results(dermi,5+1)=nan;
            end
        else % one good dipole
            results(dermi,3+dipi)=0; %distance
            results(dermi,5+dipi)=0; %layer error
            results(dermi,3+3-dipi)=nan;
            results(dermi,5+3-dipi)=nan;
        end
    else % more than 2 dipoles found
        if sum(goodfit)==2
            results(dermi,4:7)=0; %no layer or distance error
        elseif sum(goodfit)==1
            gooddip=dipi(find(goodfit));
            baddip=3-gooddip;
            results(dermi,3+gooddip)=0; %distance
            results(dermi,5+gooddip)=0; %layer error
            pntimissed=pnti;
            pntimissed(dipi(goodfit))=[];
            dist=sqrt(sum((pnt(pntinv(baddip),:)-pnt(pntimissed,:)).^2,2));
            distA=sqrt(sum((pnt(pntinv,:)-repmat(pnt(pntimissed,:),size(pntinv,1),1)).^2,2));
            distA(find(goodfit))=1000;
            [distA,Ai]=min(distA);
            results(dermi,3+baddip)=distA;
            results(dermi,5+baddip)=layer(pntimissed)-layer(pntinv(Ai)); % inside layer is 1. pos diff is bias to the center
        elseif sum(goodfit)==0
            % got to match them
            distA=sqrt(sum((pnt(pntinv,:)-repmat(pnt(pnti(1),:),size(pntinv,1),1)).^2,2));
            distB=sqrt(sum((pnt(pntinv,:)-repmat(pnt(pnti(2),:),size(pntinv,1),1)).^2,2));
            if min(distA)<min(distB)
                [distA,Ai]=min(distA);
                pntiA=pntinv(Ai);
                %distB(Ai,:)=[];
                [distB,Bi]=sort(distB);
                if Bi(1)==Ai
                    Bi=Bi(2);
                    distB=distB(2);
                else
                    Bi=Bi(1);
                    distB=distB(1);
                end
                results(dermi,4)=distA;
                results(dermi,5)=distB;
                results(dermi,6)=layer(pnti(1))-layer(pntinv(Ai));
                results(dermi,7)=layer(pnti(2))-layer(pntinv(Bi));
            elseif min(distA)>min(distB)
                [distB,Bi]=min(distB);
                pntiB=pntinv(Bi);
                %distB(Ai,:)=[];
                [distA,Ai]=sort(distA);
                if Ai(1)==Bi
                    Ai=Ai(2);
                    distA=distA(2);
                else
                    Ai=Ai(1);
                    distA=distA(1);
                end
                results(dermi,4)=distA;
                results(dermi,5)=distB;
                results(dermi,6)=layer(pnti(1))-layer(pntinv(Ai));
                results(dermi,7)=layer(pnti(2))-layer(pntinv(Bi)); %#ok<*SAGROW>
            else
                %disp('same distance?!')
                results(dermi,4:7)=nan;
            end
        end
        %disp('I say!');
    end
    
    prog(dermi)
end
disp('done');
%save resultsRnr results input resultsR

figure;scatter(results(:,3),results(:,2))

binSum=[];
binNacc=[];
binDacc=[];
binD=[];
for bini=3:15
    rowi=logical((results(:,3)>=(10*bini-10)).*(results(:,3)<(bini*10)));
    binSum(bini-2)=sum(rowi);
    binNacc(bini-2)=mean(results(rowi,1)); % ratio of number of dipoles (1 of 2 etc)
    binDacc(bini-2)=mean(results(rowi,2)); % ratio of right location
    binD(bini-2)=sum(nansum(results(rowi,4:5),2))./sum(results(rowi,1));
end
figure;bar([35:10:155],binNacc)
xlabel('Distance between dipoles (mm)')
ylabel('detected dipoles')
title('How many dipoles were found')
figure;bar([35:10:155],binDacc)
xlabel('Distance between dipoles (mm)')
ylabel('correctly detected dipoles')
title('How many dipoles were correctly localized')
figure;bar([35:10:155],binD)
xlabel('Distance between dipoles (mm)')
title('Mean distance error')
ylabel('Distance error per dipole (mm)')
disp(['min R = ',num2str(min(resultsR))])

layA=[results(:,6);results(:,7)];
layA(isnan(layA))=[];
layA(layA==0)=[];
[~,p]=ttest(layA)
mean(layA)
% distErr=results(results(:,1)==2,4:5);
% distErr=max(distErr,[],2);
% figure;scatter(results(results(:,1)==2,3),distErr)
% title('distance error for 2 dipoles')
% xlabel('Distance between 2 dipoles (mm)')
% ylabel('Localization error (mm)')

