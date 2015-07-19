%% 
nPNT=642;
try
    cd /home/yuval/Data/marik/som2/1
catch err
    cd /home/oshrit/MyDocuments/DATA/som2/1
end
hs=ft_read_headshape('hs_file');
hs=hs.pnt*1000;

inward=[10 20 30];

%% Creating sphearical grid with 3 layers - The grid has 920 points
%% Layer 1
[pnt1,lf1,vol]=sphereGrid(nPNT,inward(1)); % inner
%% Unit vectors at each point
pnt_complex1=findPerpendicular(pnt1);
%% Layer 2
[pnt2,lf2]=sphereGrid(nPNT,inward(2)); % outer
pnt_complex2=findPerpendicular(pnt2);
%% Layer 3
[pnt3,lf3,vol]=sphereGrid(nPNT,inward(3)); % inner
pnt_complex3=findPerpendicular(pnt3);
pnt=[pnt1(:,:,1);pnt2(:,:,1);pnt3(:,:,1)];
%% To cheack if inside head (ini = inside)
ini=true(length(pnt),1);
ini(pnt(:,3)<vol.o(3))=false; % keep only top half sphere
ini(pnt(:,2)<min(hs(:,2)'+10))=false;
ini(pnt(:,2)>max(hs(:,2)'-10))=false;
ini=find(ini);
pnt=pnt(ini,:);

%% xyz tangent1(anterior) tangent2 radial
pnt_complex=pnt_complex1;
pnt_complex(nPNT+1:nPNT*2,:,:)=pnt_complex2;
pnt_complex(nPNT*2+1:nPNT*3,:,:)=pnt_complex3;
pnt_complex=pnt_complex(ini,:,:);
layer(1:nPNT,1)=1;layer(nPNT+1:nPNT*2)=2;layer(nPNT*2+1:nPNT*3)=3;
layer=layer(ini);

%% Lead field (Foward solution of new vectors) 
try
    cd /home/yuval/Data/marik/som2/2
catch err
    cd /home/oshrit/MyDocuments/DATA/som2/2
end
[~,lf4]=sphereGrid(nPNT,inward(1));
[~,lf5]=sphereGrid(nPNT,inward(2));
[~,lf6]=sphereGrid(nPNT,inward(3));

try
    cd /home/yuval/Data/marik/som2/3
catch err
    cd /home/oshrit/MyDocuments/DATA/som2/3
end
[~,lf7]=sphereGrid(nPNT,inward(1));
[~,lf8]=sphereGrid(nPNT,inward(2));
[~,lf9]=sphereGrid(nPNT,inward(3));

%% Lead field of 3 Head positions (sensor positions)
LF1=lf1; %head pos 1
LF1.leadfield(nPNT+1:nPNT*2)=lf2.leadfield;
LF1.leadfield(nPNT*2+1:nPNT*3)=lf3.leadfield;
LF1.leadfield=LF1.leadfield(ini);
LF1.pos(nPNT+1:nPNT*2,:)=lf2.pos;
LF1.pos(nPNT*2+1:nPNT*3,:)=lf3.pos;
LF1.pos=LF1.pos(ini,:);
LF1.inside(nPNT+1:nPNT*2)=lf2.inside;
LF1.inside(nPNT*2+1:nPNT*3)=lf2.inside;
LF1.inside=LF1.inside(ini);
LF2=lf4; %head pos 2
LF2.leadfield(nPNT+1:nPNT*2)=lf5.leadfield;
LF2.leadfield(nPNT*2+1:nPNT*3)=lf6.leadfield;
LF2.leadfield=LF2.leadfield(ini);
LF2.pos(nPNT+1:nPNT*2,:)=lf5.pos;
LF2.pos(nPNT*2+1:nPNT*3,:)=lf6.pos;
LF2.pos=LF2.pos(ini,:);
LF2.inside(nPNT+1:nPNT*2)=lf5.inside;
LF2.inside(nPNT*2+1:nPNT*3)=lf6.inside;
LF2.inside=LF2.inside(ini);
LF3=lf7; %pos3 
LF3.leadfield(nPNT+1:nPNT*2)=lf8.leadfield;
LF3.leadfield(nPNT*2+1:nPNT*3)=lf9.leadfield;
LF3.leadfield=LF3.leadfield(ini);
LF3.pos(nPNT+1:nPNT*2,:)=lf8.pos;
LF3.pos(nPNT*2+1:nPNT*3,:)=lf9.pos;
LF3.pos=LF3.pos(ini,:);
LF3.inside(nPNT+1:nPNT*2)=lf8.inside;
LF3.inside(nPNT*2+1:nPNT*3)=lf9.inside;
LF3.inside=LF3.inside(ini);

close all
figure
plot3pnt(hs,'.')
hold on
plot3pnt(pnt,'.k')

%% Data with 3N channels (N=248) - Virtual Helmeat
try
    cd /home/yuval/Data/marik/som2/
catch err
    cd /home/oshrit/MyDocuments/DATA/som2/
end
load avgFilt % filtered data 1-20 Hz
[data, cfg]=marikViUnite(avg1_handR,avg2_handR,avg3_handR);
%LF.label=data.dataU.label

%% gain matrix (2 sources for each point 920*2=1840), so we get the FW
%% solution first for tangent1 and then for tanget2 - y
%% in the dimenssions of 3N (virtual helmet) - x
gain=[];
for pnti=1:length(pnt)
    dip=(pnt_complex(pnti,:,2)-pnt_complex(pnti,:,1))*LF1.leadfield{pnti}';
    dip=[dip,(pnt_complex(pnti,:,2)-pnt_complex(pnti,:,1))*LF2.leadfield{pnti}'];
    dip=[dip,(pnt_complex(pnti,:,2)-pnt_complex(pnti,:,1))*LF3.leadfield{pnti}'];
    gain(1:744,pnti)=dip;
    dip=(pnt_complex(pnti,:,3)-pnt_complex(pnti,:,1))*LF1.leadfield{pnti}';
    dip=[dip,(pnt_complex(pnti,:,3)-pnt_complex(pnti,:,1))*LF2.leadfield{pnti}'];
    dip=[dip,(pnt_complex(pnti,:,3)-pnt_complex(pnti,:,1))*LF3.leadfield{pnti}'];
    gain(1:744,length(pnt)+pnti)=dip;
end
% figure; topoplot248(gain(1:248,1))
% figure; topoplot248(gain(1:248,921))


%% foot
samp=180 ; % 145  %% specific time point
[dataF, ~]=marikViUnite(avg1_footL,avg2_footL,avg3_footL);
figure;plot(data.data1.time,data.data1.avg)
MF=dataF.dataU.avg(:,samp); % Meassuared Foot (Origina Data unite)

topoplot248(MF(1:248))

Pow=zeros(length(gain),1);
tic
for permi=1:1000
    Ran=[];
    
    [~,ran]=sort(rand(1,length(gain)/2));
    selected=ran(1:5);
    Ran=[Ran;selected];
    
    srcPerm=false(length(gain)/2,1);
    srcPerm(Ran)=true;
    Gain=gain(:,[srcPerm;srcPerm]);
    source=Gain\MF; % left divide source localization 
    recon=Gain*source; % recounstruction
    R=corr(recon,MF).^100; % Goodness of fit of reciunstracted vs. Meassured (Original - unite) 
    pow=zeros(size(Pow));
    pow([srcPerm;srcPerm])=source*R;
    Pow=Pow+pow;
    prog(permi)
end
toc
Pow1=sqrt(Pow(1:920).^2+Pow(921:1840).^2); % Power (vector magnitude)
figure;
scatter3pnt(pnt,25,Pow1)
[~,maxPNT]=max(Pow1);
hold on
scatter3(pnt(maxPNT,1),pnt(maxPNT,2),pnt(maxPNT,3),30,0)

[~,sortPNT]=sort(Pow1,'descend');
sortPNT(1:5)


%% hand+foot
MHF=M+MF;%*2.5;
Pow=zeros(length(srci)*2,1);
tic
for permi=1:10000
    [~,ran]=sort(rand(1,length(srci)));
    Ran=ran(1:5);
    srcPerm=false(size(srci));
    srcPerm(Ran)=true;
    Gain=gain(:,[srcPerm;srcPerm]);
    source=Gain\MHF;
    recon=Gain*source;
    R=corr(recon,MHF).^100;
    pow=zeros(size(Pow));
    pow([srcPerm;srcPerm])=source*R;
    Pow=Pow+pow;
    prog(permi)
end
toc
Pow1=sqrt(Pow(1:920).^2+Pow(921:1840).^2);
figure;
scatter3pnt(pnt,25,Pow1)
[~,maxPNT]=max(Pow1);
hold on
scatter3(pnt(maxPNT,1),pnt(maxPNT,2),pnt(maxPNT,3),30,0)

[~,sortPNT]=sort(Pow1,'descend');
sortPNT(1:5)

%% one layers
MHF248=MHF(1:248);%*2.5;
Pow=zeros(length(srci)*2,1);
tic
for permi=1:10000
    [~,ran]=sort(rand(1,length(srci)));
    Ran=ran(1:10);
    srcPerm=false(size(srci));
    srcPerm(Ran)=true;
    Gain=gain(1:248,[srcPerm;srcPerm]);
    source=Gain\MHF248;
    recon=Gain*source;
    R=corr(recon,MHF248).^100;
    pow=zeros(size(Pow));
    pow([srcPerm;srcPerm])=source*R;
    Pow=Pow+pow;
    prog(permi)
end
toc
Pow1=sqrt(Pow(1:920).^2+Pow(921:1840).^2);
figure;
scatter3pnt(pnt,25,Pow1)
[~,maxPNT]=max(Pow1);
hold on
scatter3(pnt(maxPNT,1),pnt(maxPNT,2),pnt(maxPNT,3),30,0)

[~,sortPNT]=sort(Pow1,'descend');
sortPNT(1:5)