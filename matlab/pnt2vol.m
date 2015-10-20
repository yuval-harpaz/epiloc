cd /home/yuval/Data/marik/som2/1
[V,info]=BrikLoad('ortho+orig');

% 0 0 0 = 127 180 127
% info.ORIGIN [-127,180,127]
% info.DELTA 1    -1    -1
% give x y z in RAI order and get ijk (in ASL order)

[~,w]=unix('~/SAM_BIU/docs/coordstoijk.csh ortho+orig 1 2 3 > zero'); % 129 177 126
vox=zeros(length(pnt),3);
for pnti=1:length(pnt)
    [~,w]=unix(['~/SAM_BIU/docs/coordstoijk.csh ortho+orig ',num2str(pnt(pnti,2)),' ',num2str(-pnt(pnti,1)),' ',num2str(pnt(pnti,3)),' > vox']); % 129 177 126
    vox(pnti,1:3)=importdata('vox');
    prog(pnti)
end
% M=importdata('ijkmat.1D');
% x=1;
% y=2;
% z=3;
% 
% X=sum(M(:,1)*x)+


VV=zeros(256,256,256);
for pnti=1:920
    VV(vox(pnti,1),vox(pnti,2),vox(pnti,3))=Pow1(pnti); % ASL order, see output of 3dinfo ortho+orig
end
OptTSOut.Scale = 1;
OptTSOut.Prefix = 'test';
OptTSOut.verbose = 1;
infoNew = info;
infoNew.TAGSET_LABELS='Nasion~Left Ear~Right Ear';
if exist('test+orig.BRIK','file')
    eval('!rm test+orig*')
end
WriteBrik (VV, infoNew, OptTSOut);
VV=zeros(256,256,256);
dist=20;
for i=1:256
    for j=1:256
        for k=1:256
            neib=find(sqrt(sum(power(repmat([i,j,k],920,1)-vox,2),2))<dist); % look for neighbouring pnt
            if isempty(neib) % no pnt near by
                VV(i,j,k)=0;
            else
                neibdist=sqrt(sum(power(repmat([i,j,k],size(neib,1),1)-vox(neib,:),2),2));
                gotpnt=find(neibdist==0);
                if ~isempty(gotpnt) % one of the voxels has pnt in it
                    VV(i,j,k)=Pow1(neib(gotpnt));
                else
                    if size(neib,1)==2
                        neibv=Pow1(neib);
                        %disp(['[i j k size] = ',num2str([i,j,k,length(neibv)])])
                        neibv(3)=0;
                        neibdist(3)=dist;
                    elseif size(neib,1)==1
                        neibv=Pow1(neib);
                        %disp(['[i j k size] = ',num2str([i,j,k,length(neibv)])])
                        neibv(2:3)=0;
                        neibdist(2:3)=dist;
                    else
                        neibv=Pow1(neib);
                    end
                    ratio=neibdist./sum(neibdist);
                    ratio=1./ratio;
                    ratio=ratio.^3;
                    ratio=ratio./sum(ratio);
                    VV(i,j,k)=sum(ratio.*neibv);
                end
            end
        end
    end
    disp(['done slice ',num2str(i),])
end
                