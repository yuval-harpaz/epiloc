
load('layer.mat')
load('gain1.mat')
load('pnt.mat')

load('results_3_0.1.mat')
[R,MED]=marikVirtual29plot(input,pnt,results)

load('results_3_0.1_0.5.mat')
[R_5,MED_5]=marikVirtual29plot(input,pnt,results)

load('results_3_0.1_0.75.mat')
[R_75,MED_75]=marikVirtual29plot(input,pnt,results)

MAT1=[MED_5, MED_75, MED];
MAT2=[R_5, R_75, R];

% In unpublished results, we saw that by decraesing the magnitude of one of
% 3 dipoles the sensitivity is reduced between 10 to 20 %. 