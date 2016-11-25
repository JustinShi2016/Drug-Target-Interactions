%%
function CV_Idx= GenerateIdxForCV(nCount,nCV)
%split nCount numbers in to non-overlapped nCV subsets
%CV_Idx contains nCV subsets of which each contains indices for each round
%of CV.
Division = repmat(1:nCV,1,ceil( nCount / nCV ) );
elements = Division(randperm( nCount));
elements = reshape(elements,1,nCount);
CV_Idx = cell(nCV,1); % initialize subset
for c=1:nCV
    CV_Idx{c} = find(elements==c)';
end
