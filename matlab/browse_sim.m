try 
    cd ('/home/yuval/Data/marik/som2/talk');
catch err
    cd('/home/oshrit/MyDocuments/DATA/Marik/epiloc/data/sim');
end
load('layer.mat')
load('gain1.mat')
load('pnt.mat')

try
    cd('/home/oshrit/MyDocuments/DATA/Marik/epiloc/data/sim/R_100');
catch err
    cd .
end

for noiseFactor=[0.1 0.3]
    for Ndip=1:5
        load(['results1_',num2str(Ndip),'_100_',num2str(noiseFactor),'.mat']);
        eval(['[R',num2str(Ndip),'_100_',num2str(noiseFactor*10),',MED',num2str(Ndip),'_100_',num2str(noiseFactor*10),']=marikVirtual29plot(input,pnt,results);']);
    end
end

% For comparison:
mat=[MED3_100_1, MED3_100_3];

% As R power increases we find less solutions, but with a lower distance
% error. 
