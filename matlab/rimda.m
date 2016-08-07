function [PowCur,pnti,current,fwd]=rimda(M)
% This function runs RIMDA source localization, finding dipoles for a
% measured field M, a vector of 248 values in Tesla, the number of channels
% for Magnes WH 3600 system.
% We run the function in subject directory where the headshape can be found
% (hs_file), along with data (c,rfhp0.1Hz) and the config file.

% Output
% current is the current at the chosen dipoles. pnti is the index
% of the chosen points in the spheric grid. fwd is the reconstructed field based on these
% sources.
% PowCur needs fixing, it is for multiple timepoints, not applicable for
% one time point.
PowCur=[];
% number of iterations
N=10000;
% set a spherical grid of equally distributed points (642). here we rely on
% GridSphere, an external package, although FieldTrip creates such a grid
% by default when no grid is specified for creating leadfield matrices.
nPNT=642;
hs=ft_read_headshape('hs_file');
hs=hs.pnt*1000;
inward=[10 20 30];
if exist('gain.mat','file')
    load gain
    load pnt
else
    % create 3 layers of 642 points spheres
    % compute leadfields for all the points
    [pnt1,lf1,vol]=sphereGridEL(nPNT,inward(1)); % outer
    pnt_complex1=findPerpendicular(pnt1);
    [pnt2,lf2]=sphereGridEL(nPNT,inward(2)); 
    pnt_complex2=findPerpendicular(pnt2);
    [pnt3,lf3,vol]=sphereGridEL(nPNT,inward(3)); % inner
    pnt_complex3=findPerpendicular(pnt3); % use two vectors per dipole location
    pnt=[pnt1(:,:,1);pnt2(:,:,1);pnt3(:,:,1)];
    % find indices for points inside the headshape, top half of the
    % spheres. cut points outside both ears.
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
    % index for depth of points
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
    % making a gain matrix for all the points
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
%% RIMDA stage
% run permutations: select 10 locations in each run of the loop and check
% their possible contribution to the measured field M.
Npnt=length(pnt);
% prepare an empty vector Pow to store the results of all the permutations
Pow=zeros(2*Npnt,size(M,2));
% loop fo N times (10000)
for permi=1:N
    % choose 10 locations
    Ran=[];
    [~,ran]=sort(rand(1,Npnt));
    selected=ran(1:10);
    Ran=[Ran;selected];
    srcPerm=false(1,Npnt);
    srcPerm(Ran)=true;
    % use the gain matrix of the chosen 10 points only
    Gain=gain(:,[srcPerm,srcPerm]);
    source=Gain\M; % left divide (similar to pinv)
    recon=Gain*source; % reconstruct the field from hypothesised sources
    R=corr(recon(:),M(:)).^100; % compute correlation between measured and reconstructed field
                                % use power of 100 to supress solutions
                                % that have less than excellent fit
    pow=zeros(size(Pow,1),size(M,2)); 
    pow([srcPerm,srcPerm],:)=source.*R; % weight the sources by r^100
    Pow=Pow+pow;    % sum the result with previous ones
    prog(permi) % feedback on progress
end

%% Choose sources
% this function evaluates which sources are 'local maxima' within a radius of  
% 30mm, and pass threshold of at least 0.3 strength compared to the strongest source.
% it returns current for sources that are strong maxima. pnti is the index
% of the chosen points. fwd is the reconstructed field based on these
% sources.
[current,~,pnti,fwd]=getCurrent2(Pow,pnt,M,gain,30,0.3,false);

% make some plots and prepare output
if size(M,2)==1 % for localizing more than one time sample
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
