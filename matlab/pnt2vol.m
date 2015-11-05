function pnt2vol(pnt,Pow,prefix)
% takes pnt (n by 3) + Pow (n by 1) and saves in BRIK format, got to have
% ortho+orig.BRIK in currecnt directory

if exist('ortho+orig.BRIK','file')
    [~,info]=BrikLoad('ortho+orig');
else
    error('no ortho file here')
end
if ~exist('prefix','var')
    prefix='test';
end

vox=zeros(length(pnt),3);
for pnti=1:length(pnt)
    [~,w]=unix(['~/SAM_BIU/docs/coordstoijk.csh ortho+orig ',num2str(pnt(pnti,2)),' ',num2str(-pnt(pnti,1)),' ',num2str(pnt(pnti,3)),' > vox']); % 129 177 126
    vox(pnti,1:3)=importdata('vox');
    prog(pnti)
end

OptTSOut.Scale = 1;
OptTSOut.Prefix = prefix;
OptTSOut.verbose = 1;
infoNew = info;
infoNew.TAGSET_LABELS='Nasion~Left Ear~Right Ear';

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
                    VV(i,j,k)=Pow(neib(gotpnt));
                else
                    if size(neib,1)==2
                        neibv=Pow(neib);
                        %disp(['[i j k size] = ',num2str([i,j,k,length(neibv)])])
                        neibv(3)=0;
                        neibdist(3)=dist;
                    elseif size(neib,1)==1
                        neibv=Pow(neib);
                        %disp(['[i j k size] = ',num2str([i,j,k,length(neibv)])])
                        neibv(2:3)=0;
                        neibdist(2:3)=dist;
                    else
                        neibv=Pow(neib);
                    end
                    ratio=neibdist./sum(neibdist);
                    ratio=1./ratio;
                    ratio=ratio.^2;
                    ratio=ratio./sum(ratio);
                    VV(i,j,k)=sum(ratio.*neibv);
                end
            end
        end
    end
    disp(['done slice ',num2str(i),])
end
VV=VV.*10^13;
if exist([prefix,'+orig.BRIK'],'file')
    unix(['rm ',prefix,'+orig*'])
end
WriteBrik (VV, infoNew, OptTSOut);