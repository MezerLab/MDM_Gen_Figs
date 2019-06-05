function [hip,FDRval,R2mat]=gen_young_old_matrix(young_ind,old_ind,fit_info,vec,l,r,comrl,str_vec)

% this function arranges the date of young and old subjects for each ROI
% and preforms a two-sample t-test to evaluate the significance of aging
% related changes.
% hip is a matrix of [ROIs X age group X [slope,intersection,mean MTV, mean
% qMRI parameter] X qMRI parameters]
% FDRval is a matrix of FDR corrected p-values with dimentions of [slope,
% mean MTV, mean qMRI parameter] X ROIs X qMRI parameters.
% R2mat is a matrix of R2* values for all ROIs- it is created only if there
% is R2* data in the fit_info matrix


fit_Coeff=[1 3 4]; % [slope,MTV,qMRI par,Int]
hip=nan(2*max([length(young_ind),length(old_ind)]),2,4,length(vec),length(str_vec));
for j=1:length(str_vec)
    qMRIfit=fit_info.data{str_vec(j)};
    for ii=1:length(vec)
        part=vec(ii); % ROI
        for n=1:length(fit_Coeff)
            hip(1:length(old_ind),2,n,ii,j)=squeeze(qMRIfit(part,fit_Coeff(n),old_ind)); %slope old
            hip(1:length(young_ind),1,n,ii,j)=squeeze(qMRIfit(part,fit_Coeff(n),young_ind));%slope young
            if ismember(part,l) && comrl
                ind=find(l==part);
                hip(length(old_ind)+1:2*length(old_ind),2,n,ii,j)=squeeze(qMRIfit(r(ind),fit_Coeff(n),old_ind));
                hip(length(young_ind)+1:2*length(young_ind),1,n,ii,j)=squeeze(qMRIfit(r(ind),fit_Coeff(n),young_ind));
            end
            for paremeter=1:3 % statistics for slope,MTV, qMRI par
                [hALL(paremeter,ii,j),pALL(paremeter,ii,j)] = ttest2(hip(:,1,paremeter,ii,j),hip(:,2,paremeter,ii,j));
            end
        end
    end
end

R2mat=[];

if length(fit_info.str)>4
    str_vec=[1 2 3 4]; % which qMRI parameters to take
    R2mat=hip(:,:,:,:,5);
    hip=hip(:,:,:,:,1:4);
    pALL=pALL(:,:,:);
    hALL=hALL(:,:,:);
end

pALLtest=pALL;
pALLtest(2,:,[1 3 4])=nan; % take mean MTV only for R1
[FDRval] = mafdr(pALLtest(:),'BHFDR', true);
FDRval=reshape(FDRval,size(pALL));
Hval=FDRval;
Hval(FDRval>(5./100))=nan;

end