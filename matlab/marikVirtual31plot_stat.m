function [SEQ]=marikVirtual31plot_stat(input,pnt,results)

for linei=1:size(input,1)
    D(linei)=getDist(pnt,input(linei,:),input(linei,:));
end

% binSum=[];
binNacc=[];
% binNaccSTD=[];
% binDacc=[];
binD=[];
% binDSTD=[];
binDe=[];
% binDeSTD=[];
conds={'SEQ','AVG','R'};
%R=zeros(18,1);
SEQ=[];
if size(input,2)>1
    condi=1;
        fff=4*condi-4;
        for bini=3:15
            rowi=logical((D>=(bini*10-10)).*(D<(bini*10)));
%             binSum(bini-2)=sum(rowi);
            binNacc{bini-2}=results(rowi,1+fff); % ratio of number of dipoles (1 of 2 etc)
%             binNaccSTD(bini-2)=std(results(rowi,1+fff)); % ratio of number of dipoles (1 of 2 etc)
%             binDacc(bini-2)=mean(results(rowi,2+fff)); % ratio of right location
            binD{bini-2}=results(rowi,3+fff);
%             binDSTD(bini-2)=std(results(rowi,3+fff));
            binDe{bini-2}=results(rowi,4+fff);
%             binDeSTD(bini-2)=std(results(rowi,4+fff));
        end
        %     figure;bar([35:10:155],binNacc)
        %     xlabel('Distance between dipoles (mm)')
        %     ylabel('detected dipoles')
        %     title(['How many dipoles were found',' ',conds{condi}])
        %     figure;bar([35:10:155],binDacc)
        %     xlabel('Distance between dipoles (mm)')
        %     ylabel('correctly detected dipoles')
        %     title(['How many dipoles were correctly localized',' ',conds{condi}])
        %     figure;bar([35:10:155],binD)
        %     xlabel('Distance between dipoles (mm)')
        %     title(['Mean distance error',' ',conds{condi}])
        %     ylabel('Distance error per dipole (mm)')
        %     figure;bar([35:10:155],binDe)
        %     xlabel('Distance between dipoles (mm)')
        %     title(['Mean depth distance error',' ',conds{condi}])
        %     ylabel('Depth distance error per dipole (mm)')
        
        % the 3 parameters: id_dip,dist_err, dist_depth_err)
        % distances 35mm 85mm 135mm
        
        %         if condi==3
        %             count_dist=1;
        %             for dist=find(ismember([35:10:155],[35, 55, 75]))
        %                 R(count_dist)=binNacc(dist);
        %                 R(count_dist+2)=binD(dist);
        %                 R(count_dist+4)=binDe(dist);
        %                 R(count_dist+1)=binNaccSTD(dist);
        %                 R(count_dist+3)=binDSTD(dist);
        %                 R(count_dist+5)=binDeSTD(dist);
        %                 count_dist=count_dist+6;
        %             end
        %         elseif condi==1
        count_dist=1;
        for dist=find(ismember([35:10:155],[35, 55, 75]))
            SEQ{count_dist}=binNacc(dist);
            SEQ{count_dist+1}=binD(dist);
            SEQ{count_dist+2}=binDe(dist);
%             SEQ(count_dist+1)=binNaccSTD(dist);
%             SEQ(count_dist+3)=binDSTD(dist);
%             SEQ(count_dist+5)=binDeSTD(dist);
            count_dist=count_dist+3;
        end
        %         end
        
else
    %     R(1)=mean(results(:,9));
    %     R(2)=std(results(:,9));
    %     R(3)=mean(results(:,11));
    %     R(4)=std(results(:,11));
    %     R(5)=mean(results(:,12));
    %     R(6)=std(results(:,12));
    SEQ{1}=results(:,1);
%     SEQ(2)=std(results(:,1));
    SEQ{2}=results(:,3);
%     SEQ(4)=std(results(:,3));
    SEQ{3}=results(:,4);
%     SEQ(6)=std(results(:,4));
end





function D=getDist(pnt,pnti,pntX)
D=[];
distances=[];
errVec(1)=length(pntX);
errVec(1,2)=sum(ismember(pntX,pnti));
% get all distances
for Xi=1:length(pntX)
    distances(Xi,1:length(pnti))=sqrt(sum(power(pnt(pnti,:)-repmat(pnt(pntX(Xi),:),length(pnti),1),2),2));
end
% find pairs of nearby sources and sum the error
distances(distances==0)=max(distances);
D=min(distances(:));



