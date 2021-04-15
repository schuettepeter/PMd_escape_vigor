function rsq = varexp(y,yhat)
% y=pathway_data{21}.dat';
%yhat = yhat_pbn_and_vlth_dpIns;
SS=sum(y(:).^2);
SE=(yhat-y).^2;
SSE=sum(SE(:));
rsq=1-SSE/SS;