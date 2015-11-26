
%% aud local maxima
cd /home/yuval/Data/marik/yuval/4
load Pow1000
load pnt
load avgOdd
load gain
M=avgOdd.avg(:,390);
[current,ori,pnti]=getCurrent(Pow1,pnt,M,gain)

%% older method
maxdist=30;
maxima=false(size(Pow2));
for pnti=1:length(maxima)
    pos=pnt(pnti,:);
    distnc=sqrt(sum((pnt-repmat(pos,length(pnt),1)).^2,2));
    neighb=distnc<maxdist;
    if Pow2(pnti)==max(Pow2(neighb))
        maxima(pnti)=true;
    end

end
maxima = find(maxima);
[~,maxOrder]=sort(Pow2(maxima),'descend');
maxima=maxima(maxOrder);
src=zeros(size(Pow2));
src(maxima)=length(maxima)*2-[1:length(maxima)];
figure;
scatter3pnt(pnt,25,src)
% right is 777, left is 757
Gain=[];
for maxi=1:length(maxima)
    maxj=maxima(maxi);
    maxk=maxima(maxi)+length(Pow2);
    normFac=sqrt(1/((Pow1(maxj)).^2+(Pow1(maxk)).^2));
    ori=[Pow1(maxj)*normFac Pow1(maxk)*normFac];
    Gain(1:248,maxi)=gain(:,maxj)*ori(1)+gain(:,maxk)*ori(2);
end
%source=Gain\M;
R=[];
for maxi=1:length(maxima)
    source=Gain(:,1:maxi)\M;
    R(maxi)=(corr(Gain(:,1:maxi)*source,M)).^2;
end
R
% figure;
% scatter3pnt(pnt,25,Pow2)

src=zeros(size(Pow2));
src(maxima)=source;
figure;
scatter3pnt(pnt,25,src)
fwd=Gain*source;
figure; topoplot248(fwd);

