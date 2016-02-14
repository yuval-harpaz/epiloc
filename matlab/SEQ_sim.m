try 
    cd ('/home/yuval/Data/marik/som2/talk');
catch err
    cd('/home/oshrit/MyDocuments/DATA/Marik/epiloc/data/sim');
end
load pnt

try 
    cd('/home/oshrit/MyDocuments/DATA/Marik/epiloc/data/sim/SEQ/');
catch err
    cd('./');
end
for noiseFactor=[0.1 0.3]
    for Ndip=1:5
        [results, Rcorr, Dist]=marikVirtual31(Ndip, noiseFactor);
        eval(['save SEQ_Rcorr_Dist_',num2str(Ndip),'_', num2str(noiseFactor),'.mat Dist Rcorr']);
    end
    clear all
%     cd ..
    load pnt
    cd('./SEQ')
end


for noiseFactor=[0.1 0.3]
    [miss,missR]=Simulate_tp_fn_table2(noiseFactor,25);
    
end

load pnt
load('resultsSeq1_4_0.3.mat')
[SEQ]=marikVirtual31plot(input,pnt,results)
load('results1_4_100_0.3.mat')
[R,MED]=marikVirtual29plot(input,pnt,results)
[R MED SEQ']


load pnt
load('resultsSeq1_2_0.1.mat')
[SEQ]=marikVirtual31plot(input,pnt,results)
load('resultsSeq1_2_0.3.mat')
[SEQ2]=marikVirtual31plot(input,pnt,results)
[SEQ' SEQ2']

