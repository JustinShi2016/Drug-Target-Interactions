function [DoubleSingle,SingleTargetIntOnly,SingleDrugIntOnly,Nonsingle]=...
    SplitInteractionIdx(DTI)
% Split interaction indices in terms of four scenarios of degrees
%NOTE: 
% Cannot tell where is an isolated node in DTI matrix containing isolated nodes

Degree_d= sum(DTI,2); singleDrugsRow=find(Degree_d==1);
Degree_t= sum(DTI,1); singleTargetsCol=find(Degree_t==1)';

SingleTargetInteraction = zeros(length(singleTargetsCol),1);
for col=1:length(singleTargetsCol)
    s_r(col,1) = find(DTI(:,singleTargetsCol(col))==1);
    SingleTargetInteraction(col,1)=...
        sub2ind(size(DTI),s_r(col),singleTargetsCol(col));
end

SingleDrugInteraction = zeros(length(singleDrugsRow),1);
for row=1:length(singleDrugsRow)
    s_c(row,1) = find(DTI(singleDrugsRow(row),:)==1);
    SingleDrugInteraction(row,1) = ...
        sub2ind(size(DTI),singleDrugsRow(row),s_c(row,1));
end

%common entries between SingleTargetInteraction, SingleDrugInteraction
% is the single pair
InteractionIdx = find(DTI);
DoubleSingle = intersect(SingleTargetInteraction, SingleDrugInteraction);
SingleTargetIntOnly=setdiff(SingleTargetInteraction,DoubleSingle);
SingleDrugIntOnly=setdiff(SingleDrugInteraction,DoubleSingle);
Nonsingle = setdiff(InteractionIdx,unique([SingleTargetInteraction; SingleDrugInteraction]));
