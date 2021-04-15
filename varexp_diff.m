function rsq_d = varexp_diff(y,yhat_i,yhat_j)
% y=pathway_data{21}.dat';
%yhat = yhat_pbn_and_vlth_dpIns;
SS=sum(y(:).^2);
SE_i=(yhat_i-y).^2;
SE_j=(yhat_j-y).^2;

SSE_i=sum(SE_i(:));
rsq_i=1-SSE_i/SS;

SSE_j=sum(SE_j(:));
rsq_j=1-SSE_j/SS;

rsq_d=rsq_i-rsq_j;