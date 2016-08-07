%% make subject list
cd /home/oshrit/MyDocuments/DATA/som2

%% save unaveraged data
fn='hb,xc,lf_c,rfhp0.1Hz';
cond={'handR','handL','footL'};
for runi=1:3
    cd (num2str(runi))
    trig=readTrig_BIU(fn);
    trig=clearTrig(trig);
    evt=trigOnset(trig);
    
    trl=evt'-103;
    trl(:,2)=trl+410;
    trl(:,3)=-103;
    cfg=[];
    cfg.trl=trl;
    cfg.dataset=fn;
    cfg.channel='MEG';
    cfg.hpfilter='yes';
    cfg.hpfreq=60;%[1 20];
    data=ft_preprocessing(cfg);
    
    trialinfo=trig(evt)';
    cfg=[];
    cfg.method='var';
    criterion='median';
    cfg.critval=0.5;
    [good,bad]=badTrials(cfg,data,1);
    title(['run ',num2str(runi)])
    
    trialinfo=trialinfo(good,:);
    cfg=[];
    cfg.trl=trl(good,:);
    cfg.dataset=fn;
    cfg.channel='MEG';
    cfg.bpfilter='yes';
    cfg.bpfreq=[1 20];
    data=ft_preprocessing(cfg);
    
    cfg=[];
    cfg.method='var';
    criterion='median';
    cfg.critval=3;
    [good,bad]=badTrials(cfg,data,1);

    trials=find(trialinfo);
    data.trial=data.trial(ismember(trials,good));
    data.trialinfo=trialinfo(ismember(trials,good));
    data.time=data.time(ismember(trials,good));
    data.sampleinfo=data.sampleinfo(ismember(trials,good),:);  
   
    %good=badTrials(data);
%     data=correctBL(data,[-0.1 0]);
    for condi=1:3
        eval(['data',num2str(runi),'_',cond{condi},'=data;'])
        eval(['data',num2str(runi),'_',cond{condi},'.trial=data',num2str(runi),'_',cond{condi},'.trial(data.trialinfo==condi*2);'])
        eval(['data',num2str(runi),'_',cond{condi},'.time=data',num2str(runi),'_',cond{condi},'.time(data.trialinfo==condi*2);'])
        eval(['data',num2str(runi),'_',cond{condi},'.sampleinfo=data',num2str(runi),'_',cond{condi},'.sampleinfo(data.trialinfo==condi*2,:);'])
        eval(['data',num2str(runi),'_',cond{condi},'.trialinfo=data',num2str(runi),'_',cond{condi},'.trialinfo(data.trialinfo==condi*2);'])
    end
    cd ../
end
clear data
save dataFilt data*
close all


%% averaging 1
load dataFilt
cond={'handR','handL','footL'};
runi=2;
step=[75, 150]; 
for condi=1:3
    for s=step
        eval(['vec=1:s:size(data',num2str(runi),'_',cond{condi},'.trialinfo,1);']);
        count_s2=0;
        for s2=vec
            cfg=[];
            if s2~=vec(end)
                cfg.trials=s2:(s2+s-1);
            else
               eval(['cfg.trials=s2:size(data',num2str(runi),'_',cond{condi},'.trialinfo,1);']); 
            end
            count_s2=count_s2+1;
            eval(['avg=ft_timelockanalysis(cfg,data',num2str(runi),'_',cond{condi},');']);
            eval(['avg',num2str(runi),'_',cond{condi},'_',num2str(s),'_',num2str(count_s2),'=correctBL(avg,[-0.1 0]);']);
        end
    end
    
end
clear avg
save avgFilt_group_of_trials avg*


%% averaging 2
load dataFilt
cond={'handR','handL','footL'};
runi=2;
random=[37, 75, 150];
for condi=1:3
    for r=random
        cfg=[];
        eval(['cfg.trials=randperm(size(data',num2str(runi),'_',cond{condi},'.trialinfo,1), r);']);
        eval(['avg=ft_timelockanalysis(cfg,data',num2str(runi),'_',cond{condi},');']);
        eval(['avg',num2str(runi),'_',cond{condi},'_',num2str(r),'=correctBL(avg,[-0.1 0]);']);
    end
end
clear avg
save avgFilt_random_trials avg*


%%
nPNT=642;
cd /home/oshrit/MyDocuments/DATA/som2/2
hs=ft_read_headshape('hs_file');
hs=hs.pnt*1000;

inward=[10 20 30];

[pnt1,lf1,vol]=sphereGrid(nPNT,inward(1)); % inner
pnt_complex1=findPerpendicular(pnt1);
[pnt2,lf2]=sphereGrid(nPNT,inward(2)); % outer
pnt_complex2=findPerpendicular(pnt2);
[pnt3,lf3,vol]=sphereGrid(nPNT,inward(3)); % inner
pnt_complex3=findPerpendicular(pnt3);
pnt=[pnt1(:,:,1);pnt2(:,:,1);pnt3(:,:,1)];
ini=true(length(pnt),1);
ini(pnt(:,3)<vol.o(3))=false; % keep only top half sphere
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

LF1=lf1; %head pos 1
LF1.leadfield(nPNT+1:nPNT*2)=lf2.leadfield;
LF1.leadfield(nPNT*2+1:nPNT*3)=lf3.leadfield;
LF1.leadfield=LF1.leadfield(ini);
LF1.pos(nPNT+1:nPNT*2,:)=lf2.pos;
LF1.pos(nPNT*2+1:nPNT*3,:)=lf3.pos;
LF1.pos=LF1.pos(ini,:);
LF1.inside(nPNT+1:nPNT*2)=lf2.inside;
LF1.inside(nPNT*2+1:nPNT*3)=lf2.inside;
LF1.inside=LF1.inside(ini);

gain=[];
for pnti=1:length(pnt)
    dip=(pnt_complex(pnti,:,2)-pnt_complex(pnti,:,1))*LF1.leadfield{pnti}';
    gain(1:248,pnti)=dip;
    dip=(pnt_complex(pnti,:,3)-pnt_complex(pnti,:,1))*LF1.leadfield{pnti}';
    gain(1:248,length(pnt)+pnti)=dip;
end

cd /home/oshrit/MyDocuments/DATA/som2
load avgFilt
M=avg2_handR.avg(:,138)+avg2_handL.avg(:,138)+avg2_footL.avg(:,180);

load avgFilt_random_trials
load avgFilt_group_of_trials

M=avg2_handR_150.avg(:,138);
M=avg2_handL_150.avg(:,138);
M=avg2_footL_150.avg(:,180);
M=avg2_handR_150.avg(:,138)+avg2_handL_150.avg(:,138)+avg2_footL_150.avg(:,180);

M=avg2_handR_75.avg(:,138);
M=avg2_handL_75.avg(:,138);
M=avg2_footL_75.avg(:,180);
M=avg2_handR_75.avg(:,138)+avg2_handL_75.avg(:,138)+avg2_footL_75.avg(:,180);

M=avg2_handR_37.avg(:,138);
M=avg2_handL_37.avg(:,138);
M=avg2_footL_37.avg(:,180);
M=avg2_handR_37.avg(:,138)+avg2_handL_37.avg(:,138)+avg2_footL_37.avg(:,180);

M=avg2_handR_75_1.avg(:,138);
M=avg2_handL_75_1.avg(:,138);
M=avg2_footL_75_1.avg(:,180);
M=avg2_handR_75_1.avg(:,138)+avg2_handL_75_1.avg(:,138)+avg2_footL_75_1.avg(:,180);

M=avg2_handR_75_2.avg(:,138);
M=avg2_handL_75_2.avg(:,138);
M=avg2_footL_75_2.avg(:,180);
M=avg2_handR_75_2.avg(:,138)+avg2_handL_75_2.avg(:,138)+avg2_footL_75_2.avg(:,180);

M=avg2_handR_75_3.avg(:,138);
M=avg2_handL_75_3.avg(:,138);
M=avg2_footL_75_3.avg(:,180);
M=avg2_handR_75_3.avg(:,138)+avg2_handL_75_3.avg(:,138)+avg2_footL_75_3.avg(:,180);


%% median of 10 iteration of 10000 permutations 
Pow=zeros(length(gain),10);
for NN=1:10
    N=10000;
    
    disp(' ')
    disp(['round ',num2str(NN)])
    disp(' ')
    for permi=1:N
        Ran=[];
        [~,ran]=sort(rand(1,length(gain)/2));
        selected=ran(1:10);
        Ran=[Ran;selected];
        srcPerm=false(1,length(gain)/2);
        srcPerm(Ran)=true;
        Gain=gain(:,[srcPerm,srcPerm]);
        source=Gain\M;
        recon=Gain*source;
        R=corr(recon,M).^10;
        pow=zeros(size(Pow,1),1);
        pow([srcPerm,srcPerm])=source*R;
        Pow(1:length(pow),NN)=Pow(1:length(pow),NN)+pow;
        prog(permi)
    end
end
Pow1=median(Pow,2);
Pow2=sqrt(Pow1(1:length(pnt)).^2+Pow1(length(pnt)+1:length(pnt)*2).^2);
figure;
scatter3pnt(pnt,25,Pow2)


% Pow1_R_100_combined3_300trials=Pow1; 
% save Pow1_R_100_combined3_300trials Pow1_R_100_combined3_300trials

% Pow1_R_1000_combined3_300trials=Pow1; 
% save Pow1_R_1000_combined3_300trials Pow1_R_1000_combined3_300trials

% Pow1_R_10_combined3_300trials=Pow1; 
% save Pow1_R_10_combined3_300trials Pow1_R_10_combined3_300trials

%% local maxima
% [current,ori,pnti]=getCurrent(Pow1,pnt,M,gain, maxdist, threshold);
[current,ori,pnti]=getCurrent(Pow1_R_100_combined3_300trials,pnt,M,gain, 30, 0.5);
[current,ori,pnti]=getCurrent(Pow1_R_100_combined3_300trials,pnt,M,gain, 30, 0.4); % !!!
[current,ori,pnti]=getCurrent(Pow1_R_100_combined3_300trials,pnt,M,gain, 20, 0.5); 
[current,ori,pnti]=getCurrent(Pow1_R_100_combined3_300trials,pnt,M,gain, 20, 0.4); 

[current,ori,pnti]=getCurrent(Pow1_R_1000_combined3_300trials,pnt,M,gain, 30, 0.5);
[current,ori,pnti]=getCurrent(Pow1_R_1000_combined3_300trials,pnt,M,gain, 30, 0.3); % !!!
[current,ori,pnti]=getCurrent(Pow1_R_1000_combined3_300trials,pnt,M,gain, 20, 0.5); 
[current,ori,pnti]=getCurrent(Pow1_R_1000_combined3_300trials,pnt,M,gain, 20, 0.3); % !!!

[current,ori,pnti]=getCurrent(Pow1_R_10_combined3_300trials,pnt,M,gain, 30, 0.5); % !!!
[current,ori,pnti]=getCurrent(Pow1_R_10_combined3_300trials,pnt,M,gain, 30, 0.4); 
[current,ori,pnti]=getCurrent(Pow1_R_10_combined3_300trials,pnt,M,gain, 20, 0.5); 
[current,ori,pnti]=getCurrent(Pow1_R_10_combined3_300trials,pnt,M,gain, 20, 0.3); 


%% 1 iteration of 10000 permutations 
% % N=10000;
% % Pow=zeros(length(gain),1);
% % tic
% % for permi=1:N
% %     Ran=[];
% %     [~,ran]=sort(rand(1,length(gain)/2));
% %     selected=ran(1:10);
% %     Ran=[Ran;selected];
% %     
% %     srcPerm=false(1,length(gain)/2);
% %     srcPerm(Ran)=true;
% %     Gain=gain(:,[srcPerm,srcPerm]);
% %     source=Gain\M;
% %     recon=Gain*source;
% %     R=corr(recon,M).^100;
% %     pow=zeros(size(Pow));
% %     pow([srcPerm,srcPerm])=source*R;
% %     Pow=Pow+pow;
% %     prog(permi)
% % end
% % toc
% % Pow1=sqrt(Pow(1:920).^2+Pow(921:1840).^2);
% % figure;
% % scatter3pnt(pnt,25,Pow1)
% % [~,maxPNT]=max(Pow1);
% % hold on
% % scatter3(pnt(maxPNT,1),pnt(maxPNT,2),pnt(maxPNT,3),30,0)

%%     
% % figure;
% % plot(avg1_foot1.time,avg1_hand2.avg,'g')
% % hold on
% % plot(avg1_foot1.time,avg1_foot2.avg,'r')
% % plot(avg1_foot1.time,avg1_foot1.avg,'k')
% % legend(cond)
% % title('position1')
% % figure;
% % plot(avg1_foot1.time,avg2_hand2.avg,'g')
% % hold on
% % plot(avg1_foot1.time,avg2_foot2.avg,'r')
% % plot(avg1_foot1.time,avg2_foot1.avg,'k')
% % legend(cond)
% % title('position2')
% % figure;
% % plot(avg1_foot1.time,avg3_hand2.avg,'g')
% % hold on
% % plot(avg1_foot1.time,avg3_foot2.avg,'r')
% % plot(avg1_foot1.time,avg3_foot1.avg,'k')
% % legend(cond)
% % title('position3')


%% Yuval auditory

cd /home/oshrit/MyDocuments/DATA/Marik/yuval/4
% load Pow1000
load pnt
load avgOdd
load gain
M=avgOdd.avg(:,390);
figure; topoplot248(M);

%% median of 10 iteration of 10000 permutations 
Pow=zeros(length(gain),10);
for NN=1:10
    N=10000;
    
    disp(' ')
    disp(['round ',num2str(NN)])
    disp(' ')
    for permi=1:N
        Ran=[];
        [~,ran]=sort(rand(1,length(gain)/2));
        selected=ran(1:10);
        Ran=[Ran;selected];
        srcPerm=false(1,length(gain)/2);
        srcPerm(Ran)=true;
        Gain=gain(:,[srcPerm,srcPerm]);
        source=Gain\M;
        recon=Gain*source;
        R=corr(recon,M).^100;
        pow=zeros(size(Pow,1),1);
        pow([srcPerm,srcPerm])=source*R;
        Pow(1:length(pow),NN)=Pow(1:length(pow),NN)+pow;
        prog(permi)
    end
end
Pow1=median(Pow,2);
Pow2=sqrt(Pow1(1:length(pnt)).^2+Pow1(length(pnt)+1:length(pnt)*2).^2);
figure;
scatter3pnt(pnt,25,Pow2)

[current,ori,pnti]=getCurrent(Pow1,pnt,M,gain, 30, 0.3);

Pow1_R_1000=Pow1;
Pow1_R_100=Pow1;

save Pow1_R_1000 Pow1_R_1000
save Pow1_R_100 Pow1_R_100

