function [current,ori,pnti]=getCurrent(pow,pnt,M,gain)

% pow is Pow1, 2 dipoles per location
maxdist=30;
threshold=0.5;

current=[];
%% find local maxima
Pow2=sqrt(pow(1:length(pnt)).^2+pow(length(pnt)+1:length(pnt)*2).^2);
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

pp=Pow2(maxima);
pp=pp/max(pp);
pnti=maxima(pp>=0.5);

% right is 777, left is 757
Gain=[];
for maxi=1:length(pnti)
    maxj=pnti(maxi);
    maxk=pnti(maxi)+length(Pow2);
    normFac=sqrt(1/((pow(maxj)).^2+(pow(maxk)).^2));
    ori(maxi,1:2)=[pow(maxj)*normFac pow(maxk)*normFac];
    Gain(1:248,maxi)=gain(:,maxj)*ori(maxi,1)+gain(:,maxk)*ori(maxi,2);
end
%source=Gain\M;
% R=[];
% for maxi=1:length(maxima)
%     source=Gain(:,1:maxi)\M;
%     R(maxi)=(corr(Gain(:,1:maxi)*source,M)).^2;
% end
% R
% figure;
% scatter3pnt(pnt,25,Pow2)
current=Gain\M;
R=(corr(Gain*current,M)).^2;
src=zeros(size(Pow2));
src(pnti)=current;
figure;
scatter3pnt(pnt,25,src)
fwd=Gain*current;
figure; topoplot248(fwd);
title([num2str(size(pnti,1)),' dipoles explain ',num2str(round(R*100)),'%'])

