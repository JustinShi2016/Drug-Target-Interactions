nComp = 8;
Task = 'S4'; 
CV = 10; %
Repetition =50;

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
    