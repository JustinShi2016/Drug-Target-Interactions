function RunCVWithRepetition(Task,nComp)
% % The default values of nComp are as follows,  tuned in different datasets in different scenarios, such as 
% % nComp: 8 -NR, 18-GPCR,  18-IC, 50-EN
% % 
% % REMARK: We just roughly tune it. The user may tune to get better results.
clc
load('..\NR.mat')
% nComp = 8;
% Task = 'S4'; 
CV = 10; %
Repetition =50;

AUC_=zeros(Repetition,1);
AUPR_=zeros(Repetition,1);
for r=1 : Repetition
    switch upper(Task)
        case 'S1'
    [ AUC_(r,1), AUPR_(r,1)]  =...
        DrugRepositioningForS1_TMF(DTI, {d_s},{t_s}, CV, nComp,false);
        case 'S2'
    [ AUC_(r,1), AUPR_(r,1)]  =...
        DrugRepositioningForS2AndS3_TMF(DTI,{d_s}, {t_s}, CV, nComp,false);
        case 'S3'
    [ AUC_(r,1), AUPR_(r,1)]  =...
        DrugRepositioningForS2AndS3_TMF(DTI', {t_s},{d_s} ,CV, nComp,false);
        case 'S4'
            CV=5;
    [ AUC_(r,1), AUPR_(r,1)]  =...
        DrugRepositioningForS4_TMF(DTI,{d_s}, {t_s}, CV, nComp,false);
    end
    disp(sprintf('%d-round CV-is done...',r) );

end

disp([Task,'-----------------------------------------------------'])
disp([mean(AUC_);mean(AUPR_)]);
disp([std(AUC_);std(AUPR_)]);