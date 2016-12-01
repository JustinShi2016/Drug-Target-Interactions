% The following demo codes were written by MATLAB. The users can apply them to run four different scenarios of predicting DTIs.
% DTI is the m X n matrix accounting for DTIs between m drugs and n targets
% d_s is the similarity matrix between m drugs
% t_s is the similarity matrix between n targets
% Four benchmark datasets are EN, IC, GPCR, NR.
% Four scenarios are S1, S2, S3 and S4
% nComp is the tunable parameters for our TMF. Its best value depends on the size and the complexity of datasets.
% The default values of nComp for the datasets are 50, 18, 18, 8 respectively.(in additon 4 for S1 only in NR since it is too small) 
% We just tuned it roughly, but the users may tune it to achieve better results.

load('..\NR.mat')
nComp = 8;
Task = 'S4'; 
CV = 10; %
    switch upper(Task)
        case 'S1' % for existing drugs and exising targets
        DrugRepositioningForS1_SVD_Regress(DTI, {d_s},{t_s}, CV, nComp,false);
        case 'S2' % for new drugs and exising targets
        DrugRepositioningForS2AndS3_SVD_Regress(DTI,{d_s}, {t_s}, CV, nComp,false);
        case 'S3' % for existing drugs and new targets
        DrugRepositioningForS2AndS3_SVD_Regress(DTI', {t_s},{d_s} ,CV, nComp,false);
        case 'S4' % for new drugs and new targets
            CV=5;
        DrugRepositioningForS4_SVD_Regress(DTI,{d_s}, {t_s}, CV, nComp,false);
    end
    