function [stars]=pval2stars(pvalm,StarsOrNum)

if strcmp(StarsOrNum,'num')
    if pvalm<=1E-3
        stars='p<0.001';
    elseif pvalm<=1E-2
        stars='p<0.01';
    elseif pvalm<=0.05
        stars='p<0.05';
         else
        stars='p=n.s.';
    end
elseif strcmp(StarsOrNum,'stars')
    if pvalm<=1E-3
        stars='***';
    elseif pvalm<=1E-2
        stars='**';
    elseif pvalm<=0.05
        stars='*';
    else
        stars='n.s.';
        
    end
end
end