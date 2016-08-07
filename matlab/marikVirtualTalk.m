cd /home/yuval/Data/marik/som2/talk
load pnt
load gain1
load MHR

hs=ft_read_headshape('/home/yuval/Data/marik/som2/1/hs_file');
hs=hs.pnt*1000;
%talk3
topoplot248(MHR)


srcN=size(gain,2);
Pow=zeros(srcN,1);

permi=1;
% Ran=[];
% [~,Ran]=sort(rand(1,srcN./2));
% Ran=Ran(1:10);
Ran=[328   259   549   423   115   654   588   580   330   269];
srcPerm=false(1,srcN./2);
srcPerm(Ran)=true;
figure1=figure('position',[0,0,700,700]);
scatter3pnt(pnt,25,double(srcPerm))
hold on
scatter3pnt(hs,5,'k')
view([-177 32])
grid off
axis off
set(gcf,'color','white')
colorbar off

Gain=gain(:,[srcPerm,srcPerm]);
source=Gain\MHR;
pow=zeros(size(Pow));
pow([srcPerm,srcPerm])=source;
Pow=sqrt(pow(1:920).^2+pow(921:1840).^2);




figure1=figure('position',[0,0,700,700]);
scatter3pnt(pnt,25,double(Pow))
hold on
dipi=7;
%ft_plot_dipole(pnt(Ran(dipi),:),[source([dipi;dipi+10]);0],'units','mm','color','red');
colorbar off
set(gcf,'color','white')
scatter3pnt(hs,5,'k')
view([-177 32])
grid off
axis off

figure2=figure('position',[0,0,700,700]);
scatter3pnt(pnt,25,double(Pow))
hold on
dipi=7;
ft_plot_dipole(pnt(Ran(dipi),:),[source([dipi;dipi+10]);0],'units','mm','color','red');
colorbar off
set(gcf,'color','white')
scatter3pnt(hs,5,'k')
view([-177 32])
%grid off



recon=Gain*source;

R=corr(recon,MHR).^100;
pow=zeros(size(Pow));
pow(srcPerm)=source*R;
Pow=Pow+pow;
prog(permi)


    
Pow=zeros(srcN,1);    
for permi=1:100000
    Ran=[];
    [~,Ran]=sort(rand(1,srcN./2));
    Ran=Ran(1:10);
    srcPerm=false(1,srcN./2);
    srcPerm(Ran)=true;
    Gain=gain(:,[srcPerm,srcPerm]);
    source=Gain\MHR;
    recon=Gain*source;
    R=corr(recon,MHR).^100;
    pow=zeros(size(Pow));
    pow([srcPerm,srcPerm])=source*R;
    Pow=Pow+pow;
    prog(permi)
end
Pow1=sqrt(Pow(1:920).^2+pow(921:1840).^2);



figure2=figure('position',[0,0,700,700]);
scatter3pnt(pnt,25,Pow1)
hold on
dipi=7;
% ft_plot_dipole(pnt(Ran(dipi),:),[source([dipi;dipi+10]);0],'units','mm','color','red');
% colorbar off
set(gcf,'color','white')
scatter3pnt(hs,5,'k')
view([-177 32])
grid off
axis off
colorbar off

%% 
load RFL
Pow=zeros(srcN,1);    
for permi=1:100000
    Ran=[];
    [~,Ran]=sort(rand(1,srcN./2));
    Ran=Ran(1:10);
    srcPerm=false(1,srcN./2);
    srcPerm(Ran)=true;
    Gain=gain(:,[srcPerm,srcPerm]);
    source=Gain\RFL;
    recon=Gain*source;
    R=corr(recon,RFL).^100;
    pow=zeros(size(Pow));
    pow([srcPerm,srcPerm])=source*R;
    Pow=Pow+pow;
    prog(permi)
end
Pow1=sqrt(Pow(1:920).^2+pow(921:1840).^2);



figure2=figure('position',[0,0,700,700]);
scatter3pnt(pnt,25,Pow1)
hold on
dipi=7;
% ft_plot_dipole(pnt(Ran(dipi),:),[source([dipi;dipi+10]);0],'units','mm','color','red');
% colorbar off
set(gcf,'color','white')
scatter3pnt(hs,5,'k')
view([-177 32])
grid off
axis off
colorbar off
