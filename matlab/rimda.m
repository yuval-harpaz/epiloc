function PowCur=rimda(M)

N=10000;

nPNT=642;
hs=ft_read_headshape('hs_file');
hs=hs.pnt*1000;
inward=[10 20 30];
if exist('gain.mat','file')
    load gain
    load pnt
else
    [pnt1,lf1,vol]=sphereGrid(nPNT,inward(1)); % inner
    pnt_complex1=findPerpendicular(pnt1);
    [pnt2,lf2]=sphereGrid(nPNT,inward(2)); % outer
    pnt_complex2=findPerpendicular(pnt2);
    [pnt3,lf3,vol]=sphereGrid(nPNT,inward(3)); % inner
    pnt_complex3=findPerpendicular(pnt3);
    pnt=[pnt1(:,:,1);pnt2(:,:,1);pnt3(:,:,1)];
    ini=true(length(pnt),1);
    ini(pnt(:,3)<vol.o(3)-20)=false; % keep only top half sphere
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
    LF=lf1;
    LF.leadfield(nPNT+1:nPNT*2)=lf2.leadfield;
    LF.leadfield(nPNT*2+1:nPNT*3)=lf3.leadfield;
    LF.leadfield=LF.leadfield(ini);
    LF.pos(nPNT+1:nPNT*2,:)=lf2.pos;
    LF.pos(nPNT*2+1:nPNT*3,:)=lf3.pos;
    LF.pos=LF.pos(ini,:);
    LF.inside(nPNT+1:nPNT*2)=lf2.inside;
    LF.inside(nPNT*2+1:nPNT*3)=lf2.inside;
    LF.inside=LF.inside(ini);
    gain=[];
    for pnti=1:length(pnt)
        dip=(pnt_complex(pnti,:,2)-pnt_complex(pnti,:,1))*LF.leadfield{pnti}';
        gain(1:length(dip),pnti)=dip;
        dip=(pnt_complex(pnti,:,3)-pnt_complex(pnti,:,1))*LF.leadfield{pnti}';
        gain(1:length(dip),length(pnt)+pnti)=dip;
    end
    save pnt pnt
    save layer layer
    save gain gain
end
%% 
Npnt=length(pnt);
Pow=zeros(2*Npnt,size(M,2));
for permi=1:N
    Ran=[];
    [~,ran]=sort(rand(1,Npnt));
    selected=ran(1:10);
    Ran=[Ran;selected];
    srcPerm=false(1,Npnt);
    srcPerm(Ran)=true;
    Gain=gain(:,[srcPerm,srcPerm]);
    source=Gain\M;
    recon=Gain*source;
    R=corr(recon(:),M(:)).^100;
    pow=zeros(size(Pow,1),size(M,2));
    pow([srcPerm,srcPerm],:)=source.*R;
    Pow=Pow+pow;
    prog(permi)
end


%PowSim=sqrt(Pow(1:Npnt).^2+Pow(Npnt+1:2*Npnt).^2);
%figure;scatter3pnt(pnt,25,PowSim)
[current,~,pnti,fwd]=getCurrent2(Pow,pnt,M,gain,30,0.3,false);
%[current,~,pnti,~,fwd]=getCurrent(Pow,pnt,M,gain,30,0.3,false);
if size(M,2)==1
    PowCur=zeros(Npnt,1);
    PowCur(pnti)=mean(current,2);
    cfg=[];
    cfg.interactive='no';
    cfg.comment='Measured';
    figure;
    subplot(2,2,1)
    topoplot248(mean(M,2),cfg);
    subplot(2,2,2)
    cfg.comment='Reconstructed';
    topoplot248(mean(fwd,2),cfg);
    subplot(2,2,3)
    cfg.comment='Residual';
    topoplot248(mean(M,2)-mean(fwd,2),cfg);
    subplot(2,2,4)
    plot3pnt(hs,'.k');
    hold on
    scatter3pnt(pnt(pnti,:),25,mean(current,2))
    scatter3pnt(pnt,1,PowCur)
    colorbar off
    rotate3d on
else
%     FIXME - make more plots
%     perhaps make in getCurrent2 a source estimation based on single time points, concatenated
    [~,sizei]=sort(max(current,[],2),'descend');
    figure;
    plot(current','color',[0.9 0.9 0.9])
    hold on
    for linei=1:5
        plot(current(sizei(linei),:),'color',(linei-1).*[0.15 0.15 0.15],'linewidth',6-linei)
    end
    figure;
    hold on
    for linei=1:5
        scatter3pnt(pnt(pnti(sizei(linei)),:),(6-linei)*5,(linei-1).*[0.15 0.15 0.15]);
    end
    scatter3pnt(pnt(pnti(sizei(6:end)),:),5,[0.9 0.9 0.9]);
    scatter3pnt(pnt,3,[0.9 0.9 0.9]);
end
