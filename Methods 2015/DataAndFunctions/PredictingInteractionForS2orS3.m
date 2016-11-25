function OutputMat= PredictingInteractionForS2orS3(DTI,SimD_cell, SimT_cell, nCV)
% In the scenario S2 or S3 (dx-ty), generate samples indices for k-fold cross validation
%Example:
%  PredictingInteractionForS2orS3(DTI,{d_s,d_ATC}, {t_s,t_s_Class}); %for predicting interactions for new drugs,  S2;
%  PredictingInteractionForS2orS3(DTI', {t_s,t_s_Class},{d_s,d_ATC}); %for predicting interactions for new targets, S3;
%
% INPUT:
% SimD_cell, SimT_cell  -- the cells which store several kinds of
% drug similarities and targets simillarities respectively.
% nCV                           -- the number of fold in cross validation
%
% OUTPUT: 
% OutputMat -- all predicted scores
%
% NOTE:
% This function had been proven that  performing on DTI by PerformS2 is
% equivalent to performing on DTI' by PerformS3. In addition, PerformS3 is
% provided to validate PerformS2 in a symmetric way.
%
%Writtern by J-Y. Shi, jianyushi@nwpu.edu.cn

if nargin <4
    nCV=5;
end

%make similarity matrices symmetric
for s=1:length(SimD_cell)
    SimD_cell{s}=  (SimD_cell{s}+ SimD_cell{s} )/2;
end

for s=1:length(SimT_cell)
    SimT_cell{s}=  (SimT_cell{s}+ SimT_cell{s} )/2;
end

% nCV = 5;
[nDrug, ~]= size(DTI);

rand('state',1234567890); % fix random seed for debug, S2
CV_Idx_Row= GenerateIdxForCV(nDrug,nCV);


% 1 using similarity combination
w = 0.5;
Kd =w*SimD_cell{1}+(1-w)* SimD_cell{2} ;
Kt =w*SimT_cell{1}+(1-w)*SimT_cell{2} ;

%  2 Predicting
[Outputs_ , OriginalOutput]=PerformS2(DTI,Kd,Kt,CV_Idx_Row); %


% 3 ASSERT
[~,IX] = sort( cell2mat(CV_Idx_Row) );

OriginalMat = cell2mat(OriginalOutput) ;     OriginalMat = OriginalMat (IX,:);
% % ASSERT OriginalMat==DTI
OriginalMat (OriginalMat==-1)=0;
if any ( any(OriginalMat-DTI) )
    error('Unmatched');
end

% % % put them together to assess
OutputMat = cell2mat(Outputs_) ;     OutputMat = OutputMat (IX,:);

%% Subfunctions
function [Outputs_ , OriginalOutput_, AUC_S2, AUPR_S2] = PerformS2(DTI,Kd,Kt,CV_Idx_Row)
% Input: CV_Idx_Row= GenerateIdxForCV(nDrug,nCV); %split rows
%     disp('--> supertarget: S2/S3 in S4')

cutoff_Super = 1.1; % default 1.1

nCV = length(CV_Idx_Row);
DTI_Test = DTI;
SuperTargetDTI = GetSuperTargetDTI(DTI,Kt,cutoff_Super);

AUC_S2 = zeros(nCV,1);
AUPR_S2 = zeros(nCV,1);
Outputs_= cell(nCV,1); OriginalOutput_= cell(nCV,1);

useKNN = true;
if ~useKNN
    disp('Other Classifier')
end

for k_fold_r = 1:nCV
    
    % find training entries of both interactions and non-interactions
    % Split training set and testing set.
    CV_Temp= CV_Idx_Row;    CV_Temp(k_fold_r) = [];
    ID_TrainDrug =  cell2mat(CV_Temp) ;
    clear CV_Temp;
    ID_TestDrug =  CV_Idx_Row{k_fold_r} ;
    % ID_TrainDrug = 6:nDrug;%## Unit test
    % ID_TestDrug = 1:3;%##Unit test
    
    %Output individual
    TrnDTI_S2 = DTI_Test(ID_TrainDrug,:);
    TstDTI_S2 =  DTI_Test (ID_TestDrug , :) ; % find testing entries of both interactions and non-interactions.
    Out_ = PredictMonopartite(TrnDTI_S2,TstDTI_S2, Kd(ID_TrainDrug,ID_TrainDrug), Kd(ID_TestDrug,ID_TrainDrug), useKNN);
    
   
    %Output super
    TrnDTI_S2_s = SuperTargetDTI(ID_TrainDrug,:);
    TstDTI_S2_s =  SuperTargetDTI (ID_TestDrug , :) ; % find testing entries of both interactions and non-interactions.
    Out_s = PredictMonopartite(TrnDTI_S2_s,TstDTI_S2_s, Kd(ID_TrainDrug,ID_TrainDrug), Kd(ID_TestDrug,ID_TrainDrug),useKNN);
    
    %Final scores and original labels
    Outputs_{k_fold_r} =(Out_s .* Out_); % must be values >0
    OriginalOutput_{k_fold_r}= TstDTI_S2;
    
    TrueScore = Outputs_{k_fold_r}(TstDTI_S2==1);
    FalseScore= Outputs_{k_fold_r}(TstDTI_S2~=1);
    
   
    % Assess the performance in the current CV
    [AUC_S2(k_fold_r), AUPR_S2(k_fold_r) ]=EstimationAUC(TrueScore,FalseScore,2000,0);
    disp('-');
    
end
disp([ mean(AUC_S2),  mean(AUPR_S2) ])


function [Outputs_, areaROC_ , areaPR_ ] = PredictMonopartite(TrnDTI,TstDTI,TrnSimlarity,TstSimlarity, useKNN)
% averaged assessment from drug or target at one time
if nargin<5
    useKNN = true;
end

if useKNN
    Num=3;Smooth=1; % default parameters and requirements of MLKNN
    TrnDTI = TrnDTI'; %transpose for MLKNN
    TstDTI= TstDTI';%transpose for MLKNN
    
    TrnDTI( TrnDTI~=1)  =-1; %required by MLKNN
    TstDTI( TstDTI~=1)  =-1; %required by MLKNN
    
    %KNN Classifier
    [Prior,PriorN,Cond,CondN]=MLKNN_trainWithKernel(TrnSimlarity,TrnDTI,Num,Smooth);
    [~,~,~,~,~,Outputs_,~]...
        =MLKNN_testWithKernel(TrnDTI,TstSimlarity,TstDTI,Num,Prior,PriorN,Cond,CondN);
    
    Outputs_= Outputs_'; %transpose
    %
    TstDTI=TstDTI'; %transpose again
    
else
    % We can use other classfiers here. %
    
end
%6 split scores
TrueScore = Outputs_(TstDTI==1);
FalseScore= Outputs_(TstDTI~=1);
% disp('S2/S3 in S4')
[areaROC_ , areaPR_ ]=EstimationAUC(TrueScore,FalseScore,2000,0);


%%
%%
function SuperTargetDTI= GetSuperTargetDTI(DTI,Kt, cutoff)
% T: clustering results containing cluster labels of columns
DTI(DTI ~=1)= 0;

% # --- clustering
dissimilarity = 1- Kt;
Y=squareform(dissimilarity-diag(diag(dissimilarity)) );
Z=linkage(Y,'ward');
if nargin <3
    cutoff = 0.7*(max(Z(:,3)));
end

T = cluster(Z,'cutoff',cutoff,'criterion','distance');

figure,[H,~] = dendrogram(Z,size(Kt,1),'colorthreshold',cutoff);'default';
set(H,'LineWidth',2)

% # --- Shrinking
DTI_shrink = DTI;

uT=unique(T);
MergedColIdx  =[];
t=0;
for u=1:length(uT)
    idx = find(T==uT(u));
    if length(idx)<1 %
        continue;
    else
        t=t+1;
        MergedColIdx{t,1} = idx;
    end
end

% # ---Get the interaction matrix between drugs and Super-targets
MergedCol = zeros(size(DTI_shrink,1),length(MergedColIdx));
MergedColCell = cell(1,length(MergedColIdx));

for t=1:length(MergedColIdx)
    colIdx = MergedColIdx{t};
    for c =1: length(colIdx)
        MergedCol(:,t) =  MergedCol(:,t) | DTI_shrink(:, colIdx(c) );
    end
    MergedColCell{t}= repmat(MergedCol(:,t),1,length(colIdx));
end

allMergedCol = cell2mat( MergedColIdx );

%
DTI_shrink(:,allMergedCol) = cell2mat( MergedColCell);

SuperTargetDTI = DTI_shrink;

