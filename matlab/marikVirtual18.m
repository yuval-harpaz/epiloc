%% depth normalization on sparse grid
nPNT=42;
try
    cd /home/yuval/Data/marik/som2/2
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

Ninv=10000;
Nfwd=10000;
Ndip=33;
NdipInv=33;%Ndip; % 10
Pow=zeros(size(gain,2),1);
tic
for fwdi=1:Nfwd
    Mrand=zeros(248,1);
    for dipi=1:Ndip
        %Ran=[];
        [~,Ran]=max(rand(1,size(gain,2)/2));
        %Ran=ran(1);
        srcPerm=false(1,size(gain,2)/2);
        srcPerm(Ran)=true;
        GainFwd=gain(:,[srcPerm,srcPerm]);
        %randOri=rand(1);
        %randOri(2,1)=sqrt(1-randOri^2);
        randOri=((rand(2,1)<0.5)-0.5)*2;
        Mrand=Mrand+GainFwd*randOri;
    end
    %disp('done fwd')
    for invi=1:Ninv
        Ran=[];
        [~,ran]=sort(rand(1,size(gain,2)/2));
        Ran=ran(1:NdipInv);
        srcPerm=false(1,size(gain,2)/2);
        srcPerm(Ran)=true;
        Gain=gain(:,[srcPerm,srcPerm]);
        source=Gain\Mrand;
        recon=Gain*source;
        R=corr(recon,Mrand).^100;
        pow=zeros(size(Pow));
        pow([srcPerm,srcPerm])=source*R;
        Pow=Pow+pow;
    end
    prog(fwdi)
end
toc


    
PowRand=sqrt(Pow(1:size(gain,2)./2).^2+Pow((size(gain,2)./2+1):(size(gain,2))).^2);
save (['bias_f',num2str(Nfwd),'_i',num2str(Ninv),'_marik'],'PowRand')
%eval(['PowRand',num2str(Ndip),'=PowRand;']);
figure;
scatter3pnt(pnt,25,PowRand)

figure;
scatter3pnt(pnt,25,Pow1./PowRand)





