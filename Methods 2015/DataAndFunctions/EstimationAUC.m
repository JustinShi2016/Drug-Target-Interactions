function [areaROC,areaPR,prec, tpr, fpr]=...
    EstimationAUC(TrueScore,FalseScore,count,plotFlag,runAUPR)
if nargin <5
    runAUPR = true;
     if nargin< 4
    plotFlag=0;
    if nargin<3
        count =2000;
    end
     end
end

if size(TrueScore,2)~=1
    TrueScore = TrueScore';
end
if size(FalseScore,2)~=1
        FalseScore = FalseScore';
end

if isempty(TrueScore)
    TrueScore= 0; disp('Empty score: No postive') ; % SJY 2015-01-06
end

% 
% [r_t,c_t]=size(TrueScore);
% 
% [r_f,c_f]=size(FalseScore);


ConcerningScores = [TrueScore;FalseScore]; 

[prec, tpr, fpr, thresh] =...
    prec_rec(ConcerningScores,...
    [ones(length(TrueScore),1);-ones(length(FalseScore),1)],...
    'style','.-','plotPR',plotFlag,'plotROC',plotFlag,...
    'plotBaseline',1,'numThresh',count);

if tpr ==0
    s=0;;
end
prec = [1;prec]; % append the first point
tpr = [0;tpr];
fpr = [0;fpr];



areaPR=0;
if runAUPR  
    for q=1: ( length(prec)-1 )
        areaPR=areaPR+( prec(q+1)+prec(q) )*abs( tpr(q+1) -tpr(q) )/2;
    end
end
areaROC=0;
for q=1: ( length(tpr)-1 )
    areaROC=areaROC+( tpr(q+1)+tpr(q) )*abs( fpr(q+1) -fpr(q) )/2;
end
disp(sprintf('AUC/AUPR:\t%.3f/%.3f \t #Positive = %d ',areaROC,areaPR,length(TrueScore)))

%%