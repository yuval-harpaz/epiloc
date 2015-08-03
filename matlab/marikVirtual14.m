
nPNT=642;
try
    cd /home/yuval/Data/marik/som2/1
catch
    cd /home/oshrit/MyDocuments/DATA/epiloc/data
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


% close all
% figure
% plot3pnt(hs,'.')
% hold on
% plot3pnt(pnt,'.k')



gain=[];
for pnti=1:length(pnt)
    dip=(pnt_complex(pnti,:,2)-pnt_complex(pnti,:,1))*LF1.leadfield{pnti}';
    gain(1:248,pnti)=dip;
    dip=(pnt_complex(pnti,:,3)-pnt_complex(pnti,:,1))*LF1.leadfield{pnti}';
    gain(1:248,length(pnt)+pnti)=dip;
end

cd /home/yuval/Data/marik/som2
load avgFilt avg1_handR

MH=avg1_handR.avg(:,138);
% MF=avg1_footL.avg(:,180);
% MHF=MH+MF;%*2.5;
% topoplot248(MHF(1:248))


%% hand
Pow=zeros(length(gain),1);
tic
for permi=1:100000
    Ran=[];
    
    [~,ran]=sort(rand(1,length(gain)/2));
    selected=ran(1:10);
    Ran=[Ran;selected];
    
    srcPerm=false(1,length(gain)/2);
    srcPerm(Ran)=true;
    Gain=gain(:,[srcPerm,srcPerm]);
    source=Gain\MH;
    recon=Gain*source;
    R=corr(recon,MH).^100;
    pow=zeros(size(Pow));
    pow([srcPerm,srcPerm])=source*R;
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

%% anatomy - Hand
pos=pnt(sortPNT(1),:);
% ori=pnt_complex(sortPNT(1),:,2);
% pnt_complex(
cd 1
if exist('ori.mat','file')
    load ori
    load LF
else
    
    fwd=mne_read_forward_solution('/home/yuval/Data/marik/som2/MNE/pos1_raw-oct-6-fwd.fif');
    cd 1
    cfg=[];
    cfg.fileName='/home/yuval/Data/marik/som2/MNE/pos1_raw-oct-6-fwd.fif';
    cfg.subset=false;
    [Grid,ftLeadfield,ori]=fs2grid(cfg);
    save ori ori
    LF=ft_convert_units(ftLeadfield,'mm')
    save LF LF
end

% figure;
% plot3pnt(hs,'.k')
% hold on
% plot3pnt(LF.pos,'.g')
% plot3pnt(pos,'.r')
% view([-90 90])
distnc=sqrt(sum((LF.pos-repmat(pos,length(LF.pos),1)).^2,2));
srci=find(distnc<20);
srcPos=LF.pos(srci,:);
gain=[];
for si=1:length(srci)
    gain(1:length(MH),si)=LF.leadfield{srci(si)}*ori(srci(si),:)';
end
srcN=size(gain,2);
Pow=zeros(srcN,1);
tic
for permi=1:100000
    Ran=[];
    [~,Ran]=sort(rand(1,srcN));
    Ran=Ran(1:3);
    srcPerm=false(1,srcN);
    srcPerm(Ran)=true;
    Gain=gain(:,srcPerm);
    source=Gain\MH;
    recon=Gain*source;
    R=corr(recon,MH).^100;
    pow=zeros(size(Pow));
    pow(srcPerm)=source*R;
    Pow=Pow+pow;
    prog(permi)
end
toc

figure;
scatter3pnt(srcPos,10,abs(Pow))
[~,maxPNT]=max(Pow1);
hold on
scatter3pnt(hs,1,'k')

mm=abs(gain'*MH);
figure;
scatter3pnt(srcPos,10,mm)
hold on
scatter3pnt(pnt(sortPNT(1:2),:),25,'c')
scatter3pnt(hs,1,'k')

%% make patches of neigbouring nodes

[~,pntA]=min(sqrt(sum(power((srcPos-repmat(pnt(maxPNT,:),length(srcPos),1)),2),2)));

figure;
plot3pnt(srcPos(pntA,:),'oc')
hold on
scatter3pnt(srcPos,10,'k')
scatter3pnt(pnt(sortPNT(1),:),25,'r')
%scatter3pnt(srcPos(pntA,:),25,'c')
scatter3pnt(hs,1,'k')


%nodesi(pntA)=true;

areas=[100 200 300 400];
radii=sqrt(areas./pi);

posA=srcPos(pntA,:);
radius=radii(1);
nodesi=false(length(srcPos),1);
increment=2;
for step=1:increment:radius;
    for posi=1:size(posA,1);
        distnc=sqrt(sum(power((srcPos-repmat(posA(posi,:),length(srcPos),1)),2),2));
        nodesi(distnc<increment)=true;
    end
    posA=srcPos(nodesi,:);
end
    

figure;
scatter3pnt(posA,5,'r')
hold on
%scatter3pnt(srcPos(pntA,:),25,'c')
plot3pnt(srcPos(pntA,:),'oc')
scatter3pnt(srcPos,1,'k')
%scatter3pnt(hs,1,'k')


for pnti=1:length(srcPos)
    posA=srcPos(pnti,:);
    nodesi=false(length(srcPos),1);
    for step=1:increment:radius;
        for posi=1:size(posA,1);
            distnc=sqrt(sum(power((srcPos-repmat(posA(posi,:),length(srcPos),1)),2),2));
            nodesi(distnc<increment)=true;
        end
        posA=srcPos(nodesi,:);
    end
    recon=sum(gain(:,nodesi),2);
    R(pnti)=corr(recon,MH).^2;
    prog(pnti)
end
figure;
%scatter3pnt(posA,5,'r')
scatter3pnt(srcPos,5,R)
hold on
[~,maxi]=max(R);
%scatter3pnt(srcPos(pntA,:),25,'c')
plot3pnt(srcPos(maxi,:),'ok')


radius=radii(4);


for pnti=1:length(srcPos)
    posA=srcPos(pnti,:);
    nodesi=false(length(srcPos),1);
    for step=1:increment:radius;
        for posi=1:size(posA,1);
            distnc=sqrt(sum(power((srcPos-repmat(posA(posi,:),length(srcPos),1)),2),2));
            nodesi(distnc<increment)=true;
        end
        posA=srcPos(nodesi,:);
    end
    recon=sum(gain(:,nodesi),2);
    R(pnti)=corr(recon,MH).^2;
    prog(pnti)
end
figure;
%scatter3pnt(posA,5,'r')
scatter3pnt(srcPos,5,R)
hold on
[~,maxi]=max(R);
%scatter3pnt(srcPos(pntA,:),25,'c')
plot3pnt(srcPos(maxi,:),'ok')
%% dipole fit
% cfg=[];
% cfg.latency=avg1_handR.time(138);
% cfg.method='fieldtrip';
% dipFT=dipolefitBIU(cfg,avg1_handR);
% plot3pnt(dipFT.dip.pos,'og')
% cfg.symmetry='y';
% dipFT=dipolefitBIU(cfg,avg1_handR);
% plot3pnt(dipFT.dip.pos(1,:),'om')
% 
% figure;
% scatter3pnt(posA,5,'r')
% hold on
% %scatter3pnt(srcPos(pntA,:),25,'c')
% plot3pnt(srcPos(pntA,:),'oc')
% scatter3pnt(srcPos,1,'k')
% %scatter3pnt(hs,1,'k')
% 
% cfg=[];
% cfg.latency=[avg1_handR.time(138) avg1_handR.time(138)];
% cfg.method='pinv';
% dip=dipolefitBIU(cfg,avg1_handR);
% plot3pnt(dip.dip.pos,'og')
% cfg.symmetry='y';
% dip=dipolefitBIU(cfg,avg1_handR);
% plot3pnt(dip.dip.pos(1,:),'om')






