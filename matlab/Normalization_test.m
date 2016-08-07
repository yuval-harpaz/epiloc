%%
% bias_f1000_i100000_d* contains the calculation of noise
% * indicate the number of dipoles used
% I can use the diffrent numbers of dipoles for normalizations to calculate
% the effect on the normalization
% o is the cordinate of the origin of the sphere
% pnt in position 1
% Pow1 is data of the hand in position 1

nPNT=642;
try
    cd /home/yuval/Data/marik/som2/2
catch err
    cd /home/oshrit/MyDocuments/DATA/som2/2
    %     cd /home/oshrit/MyDocuments/DATA/epiloc/data
end
hs=ft_read_headshape('hs_file');
hs=hs.pnt*1000;

inward=[10 20 30];

[pnt1,lf1,vol]=sphereGrid(nPNT,inward(1)); % inner
pnt_complex1=findPerpendicular(pnt1);
[pnt2,lf2]=sphereGrid(nPNT,inward(2)); % outer
pnt_complex2=findPerpendicular(pnt2);
[pnt3,lf3,vol]=sphereGrid(nPNT,inward(3)); % inner
pnt_complex3=findPerpendicular(pnt3);
pnt=[pnt1(:,:,1);pnt2(:,:,1);pnt3(:,:,1)];
ini=true(length(pnt),1);
ini(pnt(:,3)<vol.o(3))=false; % keep only top half sphere
ini(pnt(:,2)<min(hs(:,2)'+10))=false;
ini(pnt(:,2)>max(hs(:,2)'-10))=false;
ini=find(ini);
pnt=pnt(ini,:);

pnt_complex=pnt_complex1;
pnt_complex(nPNT+1:nPNT*2,:,:)=pnt_complex2;
pnt_complex(nPNT*2+1:nPNT*3,:,:)=pnt_complex3;
pnt_complex=pnt_complex(ini,:,:);
layer(1:nPNT,1)=1;layer(nPNT+1:nPNT*2)=2;layer(nPNT*2+1:nPNT*3)=3;
layer=layer(ini);
save /home/oshrit/MyDocuments/DATA/epiloc/data/For_normalization/layer layer

figure;
% scatter3pnt(pnt_,25,'r')
% hold on
scatter3pnt(pnt,25,'k')
figure;
scatter3pnt(pnt,25,Pow1)
figure;
scatter3pnt(pnt,25,layer)
figure;
scatter3pnt(pnt,25,PowRand)

% [lat,lon,~] = xyz2ell(pnt(:,1)/1000,pnt(:,2)/1000,pnt(:,3)/1000);
% [arclen,~] = distance(lat1,lon1,lat2,lon2);


%% Normalozation Smoothing
% go over all points
% define neighbours if closer than 4 cm
% replace each point by the mean of its neighbours
% except if it is a large or a small one (that becomes large at 1/.), than
% leave out of the calculations of the mean
% than plot the smoothed noise as well as 1/noise to control for red points (so also smooth)

PowRand_s=zeros(size(pnt,1),1);
Th=mean(PowRand)+3*std(PowRand);
Th_inv=mean(1./PowRand)+3*std(1./PowRand);
for ly=1:3
    index=find(layer==ly);
    for v=index'
        i_neig=[];
        for p=index'
            if sqrt((pnt(v,1)-pnt(p,1))^2+(pnt(v,2)-pnt(p,2))^2+(pnt(v,3)-pnt(p,3))^2)<40
                p_n=p;
                i_neig=[i_neig,p_n]; % i_neig include the point itself (as the distance to itself is zero)
            end
        end
        i_neig_t=i_neig(PowRand(i_neig)<=Th & ((1./PowRand(i_neig))<=Th_inv));
        PowRand_s(v)=mean(PowRand(i_neig_t));
        % c=zeros(size(pnt,1),1);
        % c(i_neig)=1;
        % c(v)=2;
        % figure;
        % scatter3pnt(pnt,25,c)
    end
end
figure;
scatter3pnt(pnt,25,PowRand_s)
figure;
scatter3pnt(pnt,25,1./PowRand_s)

%% Normalization smoothing by 1./noise
PowRand_inv_s=zeros(size(pnt,1),1);
Th=mean(PowRand)+3*std(PowRand);
Th_inv=mean(1./PowRand)+3*std(1./PowRand);
% diff_neig=[];
for ly=1:3
    index=find(layer==ly);
    for v=index'
        i_neig=[];
        for p=index'
            if sqrt((pnt(v,1)-pnt(p,1))^2+(pnt(v,2)-pnt(p,2))^2+(pnt(v,3)-pnt(p,3))^2)<40
                p_n=p;
                i_neig=[i_neig,p_n]; % i_neig include the point itself (as the distance to itself is zero)
            end
        end
        i_neig_t=i_neig(PowRand(i_neig)<=Th & ((1./PowRand(i_neig))<=Th_inv));
        PowRand_inv_s(v)=mean(1./PowRand(i_neig_t));
        %             diff_neig=length(i_neig)-length(i_neig_t);
        %             if diff_neig>10
        %                 disp([num2str(diff_neig), ' ', num2str(length(i_neig))])
        %             end
        
    end
end
figure;
scatter3pnt(pnt,25,PowRand_inv_s)
figure;
scatter3pnt(pnt,25,1./PowRand_inv_s)

%% Normalization smoothing by median
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


%% For comparison
figure;
scatter3pnt(pnt,25,PowRand)
figure;
scatter3pnt(pnt,25,Pow1)
figure;
scatter3pnt(pnt,25,Pow1./PowRand_s)