function dipoleCrowds(Data_path, data, sampInt, nPNT, nHeadPos, nPerm, N_sources, R_pow, normalization_flag)
% Inputs:
% Data_path, in this directory we need to have folders 1,2,3.. (as the number of head positions). e.g., '/home/yuval/Data/marik/som2/1'
% data, we need an .avg field. Can be data of a single head position, or
% data unite for 3 head positions. e.g., [data, ~]=marikViUnite(avg1_footL,avg2_footL,avg3_footL)
% sampInt, the time point of intrest that we're solving for
% nPNT, nHeadPos, nPerm, N_sources, R_pow, normalization_flag
% nHeadPos can be either 1 or 3
% nPerm is number of permutations
% N_sources is number of sources to be randomly selected
% R_pow for the goodness of fit: the power of the correlation
% normalization_flag

% Examples:
% For 3 Head positions:
% load avgFilt % filtered data 1-20 Hz
% [dataF, ~]=marikViUnite(avg1_footL,avg2_footL,avg3_footL);
% sampInt=180;
% nPNT=642;
% nHeadPos=3;
% nPerm=1000;
% N_sources=5;
% R_pow=10;
% normalization_flag=1;
% dipoleCrowds('/home/oshrit/MyDocuments/DATA/som2/', dataF, sampInt, nPNT, nHeadPos, nPerm, N_sources, R_pow, normalization_flag)
% 
% For 1 Head position:
% load avgFilt % filtered data 1-20 Hz
% sampInt=180;
% nPNT=642;
% nHeadPos=1;
% nPerm=1000;
% N_sources=5;
% R_pow=10;
% normalization_flag=1;
% dipoleCrowds('/home/oshrit/MyDocuments/DATA/som2/', avg1_footL, sampInt, nPNT, nHeadPos, nPerm, N_sources, R_pow, normalization_flag)

eval(['cd ', Data_path]);

cd ('./1');
hs=ft_read_headshape('hs_file');
hs=hs.pnt*1000;

inward=[10 20 30];

%% Creating sphearical grid with 3 layers
% Layer 1
[pnt1,lf1_1,vol]=sphereGrid(nPNT,inward(1)); % inner
% Unit vectors at each point
pnt_complex1=findPerpendicular(pnt1);
% Layer 2
[pnt2,lf2_1]=sphereGrid(nPNT,inward(2)); % outer
pnt_complex2=findPerpendicular(pnt2);
% Layer 3
[pnt3,lf3_1,vol]=sphereGrid(nPNT,inward(3)); % inner
pnt_complex3=findPerpendicular(pnt3);

pnt=[pnt1(:,:,1);pnt2(:,:,1);pnt3(:,:,1)];
% To cheack if inside head (ini = inside)
ini=true(length(pnt),1);
ini(pnt(:,3)<vol.o(3))=false; % keep only top half sphere
ini(pnt(:,2)<min(hs(:,2)'+10))=false;
ini(pnt(:,2)>max(hs(:,2)'-10))=false;
ini=find(ini);
pnt=pnt(ini,:);

% pnt_complex contains xyz of origin and tangent1(anterior) and tangent2 and radial
pnt_complex=pnt_complex1;
pnt_complex(nPNT+1:nPNT*2,:,:)=pnt_complex2;
pnt_complex(nPNT*2+1:nPNT*3,:,:)=pnt_complex3;
pnt_complex=pnt_complex(ini,:,:);
layer(1:nPNT,1)=1; layer(nPNT+1:nPNT*2)=2; layer(nPNT*2+1:nPNT*3)=3;
layer=layer(ini);

LF1=lf1_1; %head pos 1
LF1.leadfield(nPNT+1:nPNT*2)=lf2_1.leadfield;
LF1.leadfield(nPNT*2+1:nPNT*3)=lf3_1.leadfield;
LF1.leadfield=LF1.leadfield(ini);
LF1.pos(nPNT+1:nPNT*2,:)=lf2_1.pos;
LF1.pos(nPNT*2+1:nPNT*3,:)=lf3_1.pos;
LF1.pos=LF1.pos(ini,:);
LF1.inside(nPNT+1:nPNT*2)=lf2_1.inside;
LF1.inside(nPNT*2+1:nPNT*3)=lf2_1.inside;
LF1.inside=LF1.inside(ini);

cd ..

if nHeadPos~=1
    for ii=2:nHeadPos
        eval(['cd  ./',num2str(ii)]);
        eval(['[~,lf1_',num2str(ii),']=sphereGrid(nPNT,inward(1));']);
        eval(['[~,lf2_',num2str(ii),']=sphereGrid(nPNT,inward(2));']);
        eval(['[~,lf3_',num2str(ii),']=sphereGrid(nPNT,inward(3));']);
        
        eval(['LF',num2str(ii),'=lf1_',num2str(ii),';']); %head pos 2
        eval(['LF',num2str(ii),'.leadfield(nPNT+1:nPNT*2)=lf2_',num2str(ii),'.leadfield;']);
        eval(['LF',num2str(ii),'.leadfield(nPNT*2+1:nPNT*3)=lf3_',num2str(ii),'.leadfield;']);
        eval(['LF',num2str(ii),'.leadfield=LF',num2str(ii),'.leadfield(ini);']);
        eval(['LF',num2str(ii),'.pos(nPNT+1:nPNT*2,:)=lf2_',num2str(ii),'.pos;']);
        eval(['LF',num2str(ii),'.pos(nPNT*2+1:nPNT*3,:)=lf3_',num2str(ii),'.pos;']);
        eval(['LF',num2str(ii),'.pos=LF',num2str(ii),'.pos(ini,:);']);
        eval(['LF',num2str(ii),'.inside(nPNT+1:nPNT*2)=lf2_',num2str(ii),'.inside;']);
        eval(['LF',num2str(ii),'.inside(nPNT*2+1:nPNT*3)=lf3_',num2str(ii),'.inside;']);
        eval(['LF',num2str(ii),'.inside=LF',num2str(ii),'.inside(ini);']);
        
        cd ..
    end
end

figure;
plot3pnt(hs,'.')
hold on
plot3pnt(pnt,'.k')

gain=[];
for pnti=1:length(pnt)
    dip=[];
    for ii=1:nHeadPos
        eval(['dip=[dip,(pnt_complex(pnti,:,2)-pnt_complex(pnti,:,1))*LF',num2str(ii),'.leadfield{pnti}''];']);
    end
    gain(:,pnti)=dip;
    dip=[];
    for ii=1:nHeadPos
        eval(['dip=[dip,(pnt_complex(pnti,:,3)-pnt_complex(pnti,:,1))*LF',num2str(ii),'.leadfield{pnti}''];']);
    end
    gain(:,length(pnt)+pnti)=dip;
end
% figure; topoplot248(gain(1:248,1))
% figure; topoplot248(gain(1:248,921))
disp(['Size of Gain Matrix: ', num2str(size(gain))]);

% figure;
% plot(data.data1.time,data.data1.avg)

if nHeadPos>1 %(==3)
    M=data.dataU.avg(:,sampInt); % Meassuared (Original Data)
elseif nHeadPos==1
    M=data.avg(:,sampInt);
end

figure;
topoplot248(M(1:248));

Pow=zeros(length(gain),1);
tic
for permi=1:nPerm
    [~,ran]=sort(rand(1,length(gain)/2));
    Ran=ran(1:N_sources);
    srcPerm=false(length(gain)/2,1);
    srcPerm(Ran)=true;
    Gain=gain(:,[srcPerm;srcPerm]);
    source=Gain\M; % left divide source localization
    recon=Gain*source; % recounstruction
    R=corr(recon,M).^R_pow; % Goodness of fit of reciunstracted vs. Meassured (Original - unite)
    pow=zeros(size(Pow));
    pow([srcPerm;srcPerm])=source*R;
    Pow=Pow+pow;
    prog(permi)
end
toc

Pow1=sqrt(Pow(1:(length(gain)/2)).^2+Pow((length(gain)/2+1):length(gain)).^2); % Power (vector magnitude)
figure;
scatter3pnt(pnt,25,Pow1)
[~,maxPNT]=max(Pow1);
hold on
scatter3(pnt(maxPNT,1),pnt(maxPNT,2),pnt(maxPNT,3),30,0)

[~,sortPNT]=sort(Pow1,'descend');
sortPNT(1:5)