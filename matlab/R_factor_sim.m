cd ('/home/yuval/Data/marik/som2/talk');
load('layer.mat')
load('gain1.mat')
load('pnt.mat')
cd('./R_temp')
for noiseFactor=[0.1 0.3]
    for RFactor=[1, 100, 1000]
        for Ndip=1:5
            load(['results_',num2str(Ndip),'_',num2str(RFactor),'_',num2str(noiseFactor),'.mat']);
            eval(['[R',num2str(Ndip),'_',num2str(RFactor),'_',num2str(noiseFactor*10),',MED',num2str(Ndip),'_',num2str(RFactor),...
                '_',num2str(noiseFactor*10),']=marikVirtual29plot(input,pnt,results);']);
        end
    end
end

% For comparison:
mat=[MED3_1_1, MED3_100_1, MED3_1000_1];

% As R power increases we find less solutions, but with a lower distance
% error. 
