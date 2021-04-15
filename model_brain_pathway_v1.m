function stats = model_brain_pathway_v1(obj,source_one,source_two,target_one,target_two,varargin)
% Models connections between a pair of brain regions using partial least squares regression
%
% Usage:
% ::
%
%    stats = model_brain_pathway(obj,source_one,source_two,target_one,target_two,varargin);
%
% This is a method for an fmri_data object that attempts to identify a set of
% weights in one brain region (the source) that can be used to explain activity
% in another brain region (the target region). This many to many mapping is
% performed using Partial Least Squares regression, and is based on the
% idea that connected neural populations in the two brain regions produce
% correlated signals (e.g., fMRI activation). For now this function uses
% two different sources that should project to the same brain region (e.g., PBN
% and Thalamus projecting to the amygdala) The function returns a
% set of statistics that characterize the strength of each pathway, the
% pattern of weights in the source region that reliably predict activity
% in the target region, feature-wise variance explained in the target region, test
% that show which source best predicts activity in each voxel and plots that
% shows these values.
%
% ..
%     Author and copyright information:
%
%     Copyright (C) 2019 Phil Kragel
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
%
% :Inputs:
%
%   **obj:**
%        An image object with one or more images loaded (e.g., single-trial
%        or preprocessed time-series)
%
% :Optional inputs:
%
%   **
%
%
% :Outputs:
%
%   **stats:**
%        Structure including:
%           - .PLS_betas, fmri_data object with PLS regression coefficients
%           - .PLS_bootstrap_stats, statistic image for PLS regression coefficients
%           - .target_varexp, fmri data object with indicating voxel-wise variance explained
%           -
%
% :Examples:
% ::
%
%
% :See also:

% ..
%    Programmers' notes:
% 10/10/2019 Wrote funcion (Phil Kragel)
% ..


do_alignment=true;
do_boot=true;
write_output_objects=false;
nboot = 1000;
indices=crossvalind('Kfold',size(obj.dat,2),10); %do 10-fold by default.
if any(strcmp(varargin,'Indices'))
    indices = varargin{find(strcmp(varargin,'Indices'))+1};
end
names={'Source_One' 'Source_Two' 'Target_One' 'Target_Two'}; %TODO: add option to specify these

source_one_obj=apply_mask(obj,source_one);
source_two_obj=apply_mask(obj,source_two);
target_one_obj=apply_mask(obj,target_one);
target_two_obj=apply_mask(obj,target_two);

ndim=1; %number of latent dims for PLS. TODO: add option to specify
% ndim=5; AVI paper %number of latent dims for PLS. TODO: add option to specify


%% estimate simple correlations using k-fold CV
for k=1:max(indices)
    
    train = indices~=k;
    test = ~train;
    
    %simple linear association of means for comparison
    beta_pathway_one_mean = glmfit(mean(source_one_obj.dat(:,train))',mean(target_one_obj.dat(:,train))');
    beta_pathway_two_mean = glmfit(mean(source_two_obj.dat(:,train))',mean(target_one_obj.dat(:,train))');
    beta_pathway_three_mean = glmfit(mean(source_one_obj.dat(:,train))',mean(target_two_obj.dat(:,train))');
    beta_pathway_four_mean = glmfit(mean(source_two_obj.dat(:,train))',mean(target_two_obj.dat(:,train))');
    
    
    %simple regression model predictions
    yhat_pathway_one_mean(test,:)=[ones(size(mean(source_one_obj.dat(:,test))',1),1) mean(source_one_obj.dat(:,test))'] * beta_pathway_one_mean;
    yhat_pathway_two_mean(test,:)=[ones(size(mean(source_two_obj.dat(:,test))',1),1) mean(source_two_obj.dat(:,test))'] * beta_pathway_two_mean;
    yhat_pathway_three_mean(test,:)=[ones(size(mean(source_one_obj.dat(:,test))',1),1) mean(source_one_obj.dat(:,test))'] * beta_pathway_three_mean;
    yhat_pathway_four_mean(test,:)=[ones(size(mean(source_two_obj.dat(:,test))',1),1) mean(source_two_obj.dat(:,test))'] * beta_pathway_four_mean;
    
    stats.yhat_mean(test,1)=mean(yhat_pathway_one_mean(test,:),2);
    stats.yhat_mean(test,2)=mean(yhat_pathway_two_mean(test,:),2);
    stats.yhat_mean(test,3)=mean(yhat_pathway_three_mean(test,:),2);
    stats.yhat_mean(test,4)=mean(yhat_pathway_four_mean(test,:),2);

    
    %simple correlations of means
    
    
    stats.simple_correlations(k,:)=[corr(yhat_pathway_one_mean(test),mean(target_one_obj.dat(:,test))') corr(yhat_pathway_one_mean(test),mean(target_two_obj.dat(:,test))') corr(yhat_pathway_four_mean(test),mean(target_one_obj.dat(:,test))') corr(yhat_pathway_four_mean(test),mean(target_two_obj.dat(:,test))')];
    
end
%%


for k=1:max(indices)
    to_align_dat_one{k}=source_one_obj.dat(:,indices==k);
    to_align_dat_two{k}=source_two_obj.dat(:,indices==k);
    to_align_dat_three{k}=target_one_obj.dat(:,indices==k);
    to_align_dat_four{k}=target_two_obj.dat(:,indices==k);
    
end

if do_alignment
    
    aligned_dat_one = hyperalign(to_align_dat_one{:});
    aligned_dat_two = hyperalign(to_align_dat_two{:});
    aligned_dat_three = hyperalign(to_align_dat_three{:});
    aligned_dat_four = hyperalign(to_align_dat_four{:});
    
    for k=1:max(indices)
        source_one_obj.dat(:,indices==k)=aligned_dat_one{k};
        source_two_obj.dat(:,indices==k)=aligned_dat_two{k};
        target_one_obj.dat(:,indices==k)=aligned_dat_three{k};
        target_two_obj.dat(:,indices==k)=aligned_dat_four{k};
    end
end

%% estimate pattern-based connectivity performance with kfold cval

for k=1:max(indices)
    
    train = indices~=k;
    test = ~train;
    
    %pls models to estimate correlation between latent sources
    [xl_pathway_one,yl_pathway_one,xs_pathway_one,ys_pathway_one,beta_pathway_one,~,~,st_one] = plsregress(source_one_obj.dat(:,train)',target_one_obj.dat(:,train)',ndim);
    [xl_pathway_two,yl_pathway_two,xs_pathway_two,ys_pathway_two,beta_pathway_two,~,~,st_two] = plsregress(source_two_obj.dat(:,train)',target_one_obj.dat(:,train)',ndim);
    [xl_pathway_three,yl_pathway_three,xs_pathway_three,ys_pathway_three,beta_pathway_three,~,~,st_three] = plsregress(source_one_obj.dat(:,train)',target_two_obj.dat(:,train)',ndim);
    [xl_pathway_four,yl_pathway_four,xs_pathway_four,ys_pathway_four,beta_pathway_four,~,~,st_four] = plsregress(source_two_obj.dat(:,train)',target_two_obj.dat(:,train)',ndim);
    
    meanX = mean(source_one_obj.dat(:,train)',1);
    meanY = mean(target_one_obj.dat(:,train)',1);
    beta_pathway_one_crossed=st_one.W*yl_pathway_two';
    beta_pathway_one_crossed = [meanY - meanX*beta_pathway_one_crossed; beta_pathway_one_crossed];
    
    meanX = mean(source_two_obj.dat(:,train)',1);
    meanY = mean(target_two_obj.dat(:,train)',1);
    beta_pathway_four_crossed=st_four.W*yl_pathway_three';
    beta_pathway_four_crossed = [meanY - meanX*beta_pathway_four_crossed; beta_pathway_four_crossed];
    
    
    %combined pls models to see if more variance is explained when jointly modeled
    [~,~,~,~,beta_pathways_one_and_two] = plsregress([source_one_obj.dat(:,train)' source_two_obj.dat(:,train)'],target_one_obj.dat(:,train)',ndim);
    [~,~,~,~,beta_pathways_three_and_four] = plsregress([source_one_obj.dat(:,train)' source_two_obj.dat(:,train)'],target_two_obj.dat(:,train)',ndim);
    
    
    %pls model predictions
    yhat_pathway_one(test,:)=[ones(size(source_one_obj.dat(:,test)',1),1) source_one_obj.dat(:,test)'] * beta_pathway_one;
    yhat_pathway_two(test,:)=[ones(size(source_two_obj.dat(:,test)',1),1) source_two_obj.dat(:,test)'] * beta_pathway_two;
    yhat_pathway_three(test,:)=[ones(size(source_one_obj.dat(:,test)',1),1) source_one_obj.dat(:,test)'] * beta_pathway_three;
    yhat_pathway_four(test,:)=[ones(size(source_two_obj.dat(:,test)',1),1) source_two_obj.dat(:,test)'] * beta_pathway_four;
    
    yhat_pathway_one_crossed(test,:)=[ones(size(source_one_obj.dat(:,test)',1),1) source_one_obj.dat(:,test)'] * beta_pathway_one_crossed;
    yhat_pathway_four_crossed(test,:)=[ones(size(source_two_obj.dat(:,test)',1),1) source_two_obj.dat(:,test)'] * beta_pathway_four_crossed;    
    
    stats.yhat(test,1)=mean(yhat_pathway_one(test,:),2);
    stats.yhat(test,4)=mean(yhat_pathway_four(test,:),2);
    
    stats.yhat(test,2)=mean(yhat_pathway_one_crossed(test,:),2);
    stats.yhat(test,3)=mean(yhat_pathway_four_crossed(test,:),2);  
        
    stats.yhat_allvox{1}(test,:)=yhat_pathway_one(test,:);
    stats.yhat_allvox{2}(test,:)=yhat_pathway_one_crossed(test,:);
    stats.yhat_allvox{3}(test,:)=yhat_pathway_four_crossed(test,:);
    stats.yhat_allvox{4}(test,:)=yhat_pathway_four(test,:);
    
    
    ysub_one=ys_pathway_one*yl_pathway_one';
    ysub_two=ys_pathway_two*yl_pathway_two';
    ysub_three=ys_pathway_three*yl_pathway_three';
    ysub_four=ys_pathway_four*yl_pathway_four';
   
    stats.yhat_corr(k,1)=mean(diag(corr(yhat_pathway_one(test,:),target_one_obj.dat(:,test)')));
    stats.yhat_corr(k,2)=mean(diag(corr(yhat_pathway_one_crossed(test,:),target_one_obj.dat(:,test)')));
    stats.yhat_corr(k,3)=mean(diag(corr(yhat_pathway_four_crossed(test,:),target_two_obj.dat(:,test)')));
    stats.yhat_corr(k,4)=mean(diag(corr(yhat_pathway_four(test,:),target_two_obj.dat(:,test)')));

    
    
    stats.yhat_weighted(test,1)=yhat_pathway_one(test,:)*diag(corr(ys_pathway_one*yl_pathway_one',target_one_obj.dat(:,train)')).^2;
    stats.yhat_weighted(test,2)=yhat_pathway_one_crossed(test,:)*diag(corr(ys_pathway_one*yl_pathway_two',target_one_obj.dat(:,train)')).^2;
    stats.yhat_weighted(test,3)=yhat_pathway_four_crossed(test,:)*diag(corr(ys_pathway_four*yl_pathway_three',target_two_obj.dat(:,train)')).^2;
    stats.yhat_weighted(test,4)=yhat_pathway_four(test,:)*diag(corr(ys_pathway_four*yl_pathway_four',target_two_obj.dat(:,train)')).^2;  
   
    %combined pls model predictions
    yhat_pathways_one_and_two(test,:)=[ones(size(source_one_obj.dat(:,test)',1),1) [source_one_obj.dat(:,test)' source_two_obj.dat(:,test)']] * beta_pathways_one_and_two;
    yhat_pathways_three_and_four(test,:)=[ones(size(source_one_obj.dat(:,test)',1),1) [source_one_obj.dat(:,test)' source_two_obj.dat(:,test)']] * beta_pathways_three_and_four;
    
    
    %correlation of latent dimensions
    Ytest_target_one=(target_one_obj.dat(:,test))';
    for i=1:size(Ytest_target_one,2)
        Ytest_target_one(:,i)=Ytest_target_one(:,i)-mean(Ytest_target_one(:,i));
    end
    
    Ytest_target_two=(target_two_obj.dat(:,test))';
    for i=1:size(Ytest_target_two,2)
        Ytest_target_two(:,i)=Ytest_target_two(:,i)-mean(Ytest_target_two(:,i));
    end
    
    Xtest_source_one=(source_one_obj.dat(:,test))';
    for i=1:size(Xtest_source_one,2)
        Xtest_source_one(:,i)=Xtest_source_one(:,i)-mean(Xtest_source_one(:,i));
    end
    
    
    Xtest_source_two=(source_two_obj.dat(:,test))';
    for i=1:size(Xtest_source_two,2)
        Xtest_source_two(:,i)=Xtest_source_two(:,i)-mean(Xtest_source_two(:,i));
    end
    
    %optimized models
    YS_Test_target_one_pathway_one = Ytest_target_one*yl_pathway_one;
    YS_Test_target_one_pathway_two = Ytest_target_one*yl_pathway_two;
    YS_Test_target_two_pathway_three = Ytest_target_two*yl_pathway_three;
    YS_Test_target_two_pathway_four = Ytest_target_two*yl_pathway_four;
    
    XS_Test_source_one_pathway_one = Xtest_source_one*xl_pathway_one;
    XS_Test_source_two_pathway_two = Xtest_source_two*xl_pathway_two;
    XS_Test_source_one_pathway_three = Xtest_source_one*xl_pathway_three;
    XS_Test_source_two_pathway_four = Xtest_source_two*xl_pathway_four;
    
    stats.latent_timeseries(indices==k,1)=YS_Test_target_one_pathway_one(:,1);
    stats.latent_timeseries(indices==k,4)=YS_Test_target_two_pathway_four(:,1);
    
    %optimized pathways
    latent_correlation_pathway_one(k,:)=diag(corr(YS_Test_target_one_pathway_one,XS_Test_source_one_pathway_one));
    latent_correlation_pathway_four(k,:)=diag(corr(YS_Test_target_two_pathway_four,XS_Test_source_two_pathway_four));
    
    %switched sources
    latent_correlation_pathway_one_crossed(k,:)=diag(corr(YS_Test_target_one_pathway_two,XS_Test_source_one_pathway_one));%;XS_Test_source_one_pathway_three
    latent_correlation_pathway_four_crossed(k,:)=diag(corr(YS_Test_target_two_pathway_three,XS_Test_source_two_pathway_four));%XS_Test_source_two_pathway_two
    
     stats.latent_timeseries(indices==k,2)=YS_Test_target_one_pathway_two(:,1);
    stats.latent_timeseries(indices==k,3)=YS_Test_target_two_pathway_three(:,1);
    
    
end


stats.latent_correlations=[latent_correlation_pathway_one(:,1) latent_correlation_pathway_one_crossed(:,1) latent_correlation_pathway_four_crossed(:,1) latent_correlation_pathway_four(:,1) ];

% on average, are 'on target' pathways more functionally connected than 'off target' pathways
[~,~,~,stats.latent_correlation_interaction_ttest]=ttest(atanh(stats.latent_correlations(:,1))-atanh(stats.latent_correlations(:,2))-atanh(stats.latent_correlations(:,3))+atanh(stats.latent_correlations(:,4)));
[~,~,~,stats.simple_correlation_interaction_ttest]=ttest(atanh(stats.simple_correlations(:,1))-atanh(stats.simple_correlations(:,2))-atanh(stats.simple_correlations(:,3))+atanh(stats.simple_correlations(:,4)));


[~,~,~,stats.latent_correlation_pathway_one_ttest]=ttest(atanh(stats.latent_correlations(:,1))-atanh(stats.latent_correlations(:,2)));
[~,~,~,stats.latent_correlation_pathway_two_ttest]=ttest(atanh(stats.latent_correlations(:,4))-atanh(stats.latent_correlations(:,3)));
[~,~,~,stats.simple_correlation_pathway_one_ttest]=ttest(atanh(stats.simple_correlations(:,1))-atanh(stats.simple_correlations(:,2)));
[~,~,~,stats.simple_correlation_pathway_two_ttest]=ttest(atanh(stats.simple_correlations(:,4))-atanh(stats.simple_correlations(:,3)));


%% partial correlations in target for each pathway

[r]=partialcorr(yhat_pathway_one,target_one_obj.dat',[yhat_pathway_two]);
bs_R=bootstrp(1000,@partialcorr,yhat_pathway_one,target_one_obj.dat',[yhat_pathway_two]);
bs_R=reshape(bs_R,1000,size(beta_pathway_one,2),size(beta_pathway_one,2));
bs_Z=nanmean(bs_R)./nanstd(bs_R);
p = squeeze(2*normcdf(-1*abs(bs_Z),0,1));
unique_pathway_one_r=diag(r);
simple_pathway_one_r=corr(yhat_pathway_one,target_one_obj.dat');
unique_pathway_one_p=diag(p);

[r]=partialcorr(yhat_pathway_three,target_two_obj.dat',yhat_pathway_four);
bs_R=bootstrp(1000,@partialcorr,yhat_pathway_three,target_two_obj.dat',yhat_pathway_four);
bs_R=reshape(bs_R,1000,size(beta_pathway_one,2),size(beta_pathway_one,2));
bs_Z=nanmean(bs_R)./nanstd(bs_R);
p = squeeze(2*normcdf(-1*abs(bs_Z),0,1));
unique_pathway_three_r=diag(r);
simple_pathway_three_r=corr(yhat_pathway_one,target_one_obj.dat');
unique_pathway_three_p=diag(p);



[r]=partialcorr(yhat_pathway_two,target_one_obj.dat',[yhat_pathway_one]);
unique_pathway_two_r=diag(r);
simple_pathway_two_r=corr(yhat_pathway_two,target_one_obj.dat');
bs_R=bootstrp(1000,@partialcorr,yhat_pathway_two,target_one_obj.dat',[yhat_pathway_one]);
bs_R=reshape(bs_R,1000,size(beta_pathway_one,2),size(beta_pathway_one,2));
bs_Z=nanmean(bs_R)./nanstd(bs_R);
p = squeeze(2*normcdf(-1*abs(bs_Z),0,1));
unique_pathway_two_p=diag(p);

[r, p]=partialcorr(yhat_pathway_four,target_two_obj.dat',yhat_pathway_three);
bs_R=bootstrp(1000,@partialcorr,yhat_pathway_four,target_two_obj.dat',yhat_pathway_three);
bs_R=reshape(bs_R,1000,size(beta_pathway_one,2),size(beta_pathway_one,2));
bs_Z=nanmean(bs_R)./nanstd(bs_R);
p = squeeze(2*normcdf(-1*abs(bs_Z),0,1));
unique_pathway_four_r=diag(r);
simple_pathway_four_r=corr(yhat_pathway_two,target_one_obj.dat');
unique_pathway_four_p=diag(p);

%% objects with partial correlations in target one

pathway_one_partial_corr=target_one_obj;
pathway_one_partial_corr.dat=unique_pathway_one_r;
pathway_one_partial_corr.dat=pathway_one_partial_corr.dat .* double(unique_pathway_one_p<FDR(unique_pathway_one_p+1e-12,.05));
pathway_one_partial_corr.fullpath='pathway_one_partial_correlation.nii';



pathway_two_partial_corr=target_one_obj;
pathway_two_partial_corr.dat=unique_pathway_two_r;
pathway_two_partial_corr.dat=pathway_two_partial_corr.dat .* double(unique_pathway_two_p<FDR(unique_pathway_two_p+1e-12,.05));
pathway_two_partial_corr.fullpath='pathway_two_partial_correlation.nii';


if write_output_objects
    write(pathway_one_partial_corr,'overwrite')
    write(pathway_two_partial_corr,'overwrite')
end


%% objects with partial correlations in target two

pathway_three_partial_corr=target_two_obj;
pathway_three_partial_corr.dat=unique_pathway_three_r;
pathway_three_partial_corr.dat=pathway_three_partial_corr.dat .* double(unique_pathway_three_p<FDR(unique_pathway_three_p+1e-12,.05));
pathway_three_partial_corr.fullpath='pathway_three_partial_correlation.nii';



pathway_four_partial_corr=target_two_obj;
pathway_four_partial_corr.dat=unique_pathway_four_r;
pathway_four_partial_corr.dat=pathway_four_partial_corr.dat .* double(unique_pathway_four_p<FDR(unique_pathway_four_p+1e-12,.05));
pathway_four_partial_corr.fullpath='pathway_four_partial_correlation.nii';


if write_output_objects
    write(pathway_three_partial_corr,'overwrite')
    write(pathway_four_partial_corr,'overwrite')
end

%% comparisons of correlations in target one

target_one_stats=statistic_image;
for v=1:size(yhat_pathway_one,2)
    [target_one_stats.p(v,1), target_one_stats.dat(v,1)]=r_test_paired(target_one_obj.dat(v,:)',yhat_pathway_one(:,v),yhat_pathway_two(:,v),0);
end
target_one_stats.volInfo=target_one_obj.volInfo;
target_one_stats.removed_voxels=target_one_obj.removed_voxels;
% target_one_stats=threshold(target_one_stats,.05,'FDR');
% target_one_stats.dat(~target_one_stats.sig)=NaN;
target_one_stats.fullpath='correlation_difference_pathway_one_vs_pathway_two.nii';
if write_output_objects
    write(target_one_stats,'overwrite')
end

stats.target_one_stats=target_one_stats;
%% comparisons of correlations in target two

target_two_stats=statistic_image;
target_two_stats.volInfo=target_two_obj.volInfo;
for v=1:size(yhat_pathway_three,2)
    [target_two_stats.p(v,1), target_two_stats.dat(v,1)]=r_test_paired(target_two_obj.dat(v,:)',yhat_pathway_four(:,v),yhat_pathway_three(:,v),0);
end
target_two_stats.removed_voxels=target_two_obj.removed_voxels;

% target_two_stats=threshold(target_two_stats,.05,'FDR');
% target_two_stats.dat(~target_two_stats.sig)=NaN;
target_two_stats.fullpath='correlation_difference_pathway_three_vs_pathway_four.nii';
if write_output_objects
    write(target_two_stats,'overwrite')
end

stats.target_two_stats=target_one_stats;


%% Model parameters based on whole sample


[xl_pathway_one,yl_pathway_one,xs_pathway_one,ys_pathway_one,beta_pathway_one,~,~,stats.st{1}] = plsregress(source_one_obj.dat',target_one_obj.dat',ndim);
[xl_pathway_two,yl_pathway_two,xs_pathway_two,ys_pathway_two,beta_pathway_two,~,~,stats.st{2}] = plsregress(source_two_obj.dat',target_one_obj.dat',ndim);

[xl_pathway_three,yl_pathway_three,xs_pathway_three,ys_pathway_three,beta_pathway_three,~,~,stats.st{3}] = plsregress(source_one_obj.dat',(target_two_obj.dat)',ndim);
[xl_pathway_four,yl_pathway_four,xs_pathway_four,ys_pathway_four,beta_pathway_four,~,~,stats.st{4}] = plsregress(source_two_obj.dat',(target_two_obj.dat)',ndim);



stats.yl{1}=yl_pathway_one(:,1) ;
stats.yl{2}=yl_pathway_two(:,1) ;
stats.yl{3}=yl_pathway_three(:,1);
stats.yl{4}=yl_pathway_four(:,1);

stats.xl{1}=xl_pathway_one(:,1) ;
stats.xl{2}=xl_pathway_two(:,1) ;
stats.xl{3}=xl_pathway_three(:,1);
stats.xl{4}=xl_pathway_four(:,1);


if do_boot
    bs_beta_pathway_three=bootstrp(nboot,@bootpls,source_one_obj.dat',target_two_obj.dat');
    bs_beta_pathway_four=bootstrp(nboot,@bootpls,source_two_obj.dat',target_two_obj.dat');
    
    bs_beta_pathway_one=bootstrp(nboot,@bootpls,source_one_obj.dat',target_one_obj.dat');
    bs_beta_pathway_two=bootstrp(nboot,@bootpls,source_two_obj.dat',target_one_obj.dat');
    
    
    bs_xl_pathway_three=bootstrp(nboot,@bootpls_loadings,source_one_obj.dat',target_two_obj.dat');
    bs_xl_pathway_four=bootstrp(nboot,@bootpls_loadings,source_two_obj.dat',target_two_obj.dat');
    
    bs_xl_pathway_one=bootstrp(nboot,@bootpls_loadings,source_one_obj.dat',target_one_obj.dat');
    bs_xl_pathway_two=bootstrp(nboot,@bootpls_loadings,source_two_obj.dat',target_one_obj.dat');
    
    bs_yl_pathway_three=bootstrp(nboot,@bootpls_yloadings,source_one_obj.dat',target_two_obj.dat');
    bs_yl_pathway_four=bootstrp(nboot,@bootpls_yloadings,source_two_obj.dat',target_two_obj.dat');
    
    bs_yl_pathway_one=bootstrp(nboot,@bootpls_yloadings,source_one_obj.dat',target_one_obj.dat');
    bs_yl_pathway_two=bootstrp(nboot,@bootpls_yloadings,source_two_obj.dat',target_one_obj.dat');
end

%% mean betas for each pathway
beta_obj_pathway_one=source_one_obj;
beta_obj_pathway_one.dat=mean(beta_pathway_one(2:end,:),2);
stats.PLS_objects(1)=beta_obj_pathway_one;

beta_obj_pathway_two=source_one_obj;
beta_obj_pathway_two.dat=mean(beta_pathway_two(2:end,:),2);
stats.PLS_objects(2)=beta_obj_pathway_two;

beta_obj_pathway_three=source_two_obj;
beta_obj_pathway_three.dat=mean(beta_pathway_three(2:end,:),2);
stats.PLS_objects(3)=beta_obj_pathway_three;

beta_obj_pathway_four=source_two_obj;
beta_obj_pathway_four.dat=mean(beta_pathway_four(2:end,:),2);
stats.PLS_objects(4)=beta_obj_pathway_four;

stats.PLS_betas{1}=beta_pathway_one;
stats.PLS_betas{2}=beta_pathway_two;
stats.PLS_betas{3}=beta_pathway_three;
stats.PLS_betas{4}=beta_pathway_four;


%% bootstrap stats pathway one
if do_boot
    bs_Z=mean(bs_beta_pathway_one)./std(bs_beta_pathway_one);
    bs_Z=reshape(bs_Z,size(bs_Z,2)/size(beta_pathway_one,2),size(beta_pathway_one,2));
    bs_Z=mean(bs_Z,2);
    bs_P = 2*normcdf(-1*abs(bs_Z),0,1);
    
    
    bs_stat=statistic_image;
    bs_stat.volInfo=source_one_obj.volInfo;
    bs_stat.dat=bs_Z(2:end,1);
    bs_stat.p=bs_P(2:end,1);
    bs_stat.removed_voxels=source_one_obj.removed_voxels;
    stats.PLS_bootstrap_stats(1)=bs_stat;
    
    bs_Z=mean(bs_xl_pathway_one)./std(bs_xl_pathway_one);
    bs_P = 2*normcdf(-1*abs(bs_Z),0,1);
    
    
    bs_stat=statistic_image;
    bs_stat.volInfo=source_one_obj.volInfo;
    bs_stat.dat=bs_Z';
    bs_stat.p=bs_P';
    bs_stat.removed_voxels=source_one_obj.removed_voxels;
    stats.PLS_bootstrap_stats_xl(1)=bs_stat;
    
    
    bs_Z=mean(bs_yl_pathway_one)./std(bs_yl_pathway_one);
    bs_P = 2*normcdf(-1*abs(bs_Z),0,1);
    
    
    bs_stat=statistic_image;
    bs_stat.volInfo=target_one_obj.volInfo;
    bs_stat.dat=bs_Z';
    bs_stat.p=bs_P';
    bs_stat.removed_voxels=target_one_obj.removed_voxels;
    stats.PLS_bootstrap_stats_yl(1)=bs_stat;
    %% bootstrap stats pathway 3
    bs_Z=mean(bs_beta_pathway_three)./std(bs_beta_pathway_three);
    bs_Z=reshape(bs_Z,size(bs_Z,2)/size(beta_pathway_three,2),size(beta_pathway_three,2));
    bs_Z=mean(bs_Z,2);
    bs_P = 2*normcdf(-1*abs(bs_Z),0,1);
    
    
    bs_stat=statistic_image;
    bs_stat.volInfo=source_one_obj.volInfo;
    bs_stat.dat=bs_Z(2:end,1);
    bs_stat.p=bs_P(2:end,1);
    bs_stat.removed_voxels=source_one_obj.removed_voxels;
    stats.PLS_bootstrap_stats(3)=bs_stat;
    
    
    
    bs_Z=mean(bs_xl_pathway_three)./std(bs_xl_pathway_three);
    bs_P = 2*normcdf(-1*abs(bs_Z),0,1);
    
    
    bs_stat=statistic_image;
    bs_stat.volInfo=source_one_obj.volInfo;
    bs_stat.dat=bs_Z';
    bs_stat.p=bs_P';
    bs_stat.removed_voxels=source_one_obj.removed_voxels;
    stats.PLS_bootstrap_stats_xl(3)=bs_stat;
    
    
    bs_Z=mean(bs_yl_pathway_three)./std(bs_yl_pathway_three);
    bs_P = 2*normcdf(-1*abs(bs_Z),0,1);
    
    
    bs_stat=statistic_image;
    bs_stat.volInfo=target_two_obj.volInfo;
    bs_stat.dat=bs_Z';
    bs_stat.p=bs_P';
    bs_stat.removed_voxels=target_two_obj.removed_voxels;
    stats.PLS_bootstrap_stats_yl(3)=bs_stat;
    
    
    %% bootstrap stats pathway 4
    
    bs_Z=mean(bs_beta_pathway_four)./std(bs_beta_pathway_four);
    bs_Z=reshape(bs_Z,size(bs_Z,2)/size(beta_pathway_four,2),size(beta_pathway_four,2));
    bs_Z=mean(bs_Z,2);
    bs_P = 2*normcdf(-1*abs(bs_Z),0,1);
    
    
    bs_stat=statistic_image;
    bs_stat.volInfo=source_two_obj.volInfo;
    bs_stat.dat=bs_Z(2:end,1);
    bs_stat.p=bs_P(2:end,1);
    bs_stat.removed_voxels=source_two_obj.removed_voxels;
    stats.PLS_bootstrap_stats(4)=bs_stat;
    
    
    bs_Z=mean(bs_xl_pathway_four)./std(bs_xl_pathway_four);
    bs_P = 2*normcdf(-1*abs(bs_Z),0,1);
    
    
    bs_stat=statistic_image;
    bs_stat.volInfo=source_two_obj.volInfo;
    bs_stat.dat=bs_Z';
    bs_stat.p=bs_P';
    bs_stat.removed_voxels=source_two_obj.removed_voxels;
    stats.PLS_bootstrap_stats_xl(4)=bs_stat;
    
    
    bs_Z=mean(bs_yl_pathway_four)./std(bs_yl_pathway_four);
    bs_P = 2*normcdf(-1*abs(bs_Z),0,1);
    
    
    bs_stat=statistic_image;
    bs_stat.volInfo=target_two_obj.volInfo;
    bs_stat.dat=bs_Z';
    bs_stat.p=bs_P';
    bs_stat.removed_voxels=target_two_obj.removed_voxels;
    stats.PLS_bootstrap_stats_yl(4)=bs_stat;
    
    
    %% bootstrap stats pathway 2
    bs_Z=mean(bs_beta_pathway_two)./std(bs_beta_pathway_two);
    bs_Z=reshape(bs_Z,size(bs_Z,2)/size(beta_pathway_two,2),size(beta_pathway_two,2));
    bs_Z=mean(bs_Z,2);
    bs_P = 2*normcdf(-1*abs(bs_Z),0,1);
    
    
    bs_stat=statistic_image;
    bs_stat.volInfo=source_two_obj.volInfo;
    bs_stat.dat=bs_Z(2:end,1);
    bs_stat.p=bs_P(2:end,1);
    bs_stat.removed_voxels=source_two_obj.removed_voxels;
    stats.PLS_bootstrap_stats(2)=bs_stat;
    
    bs_Z=mean(bs_xl_pathway_two)./std(bs_xl_pathway_two);
    bs_P = 2*normcdf(-1*abs(bs_Z),0,1);
    
    bs_stat=statistic_image;
    bs_stat.volInfo=source_two_obj.volInfo;
    bs_stat.dat=bs_Z';
    bs_stat.p=bs_P';
    bs_stat.removed_voxels=source_two_obj.removed_voxels;
    stats.PLS_bootstrap_stats_xl(2)=bs_stat;
    
    
    bs_Z=mean(bs_yl_pathway_two)./std(bs_yl_pathway_two);
    bs_P = 2*normcdf(-1*abs(bs_Z),0,1);
    
    
    bs_stat=statistic_image;
    bs_stat.volInfo=target_one_obj.volInfo;
    bs_stat.dat=bs_Z';
    bs_stat.p=bs_P';
    bs_stat.removed_voxels=target_one_obj.removed_voxels;
    stats.PLS_bootstrap_stats_yl(2)=bs_stat;
    
    %% variance explained in target one
    
    y=target_one_obj.dat';
    yhat=yhat_pathway_one;
    rsq_pathway_one=varexp(y,yhat);
    rsq_ci_pathway_one=bootci(nboot,@varexp,y,yhat);
    
    
    yhat=yhat_pathway_two;
    rsq_pathway_two=varexp(y,yhat);
    rsq_ci_pathway_two=bootci(nboot,@varexp,y,yhat);
    
    %% variance explained in target two
    
    y=target_two_obj.dat';
    yhat=yhat_pathway_three;
    rsq_pathway_three=varexp(y,yhat);
    rsq_ci_pathway_three=bootci(nboot,@varexp,y,yhat);
    
    
    yhat=yhat_pathway_four;
    rsq_pathway_four=varexp(y,yhat);
    rsq_ci_pathway_four=bootci(nboot,@varexp,y,yhat);
    
    
    stats.varexp(1,:)=[rsq_ci_pathway_one(1) rsq_pathway_one rsq_ci_pathway_one(2)];
    stats.varexp(2,:)=[rsq_ci_pathway_two(1) rsq_pathway_two rsq_ci_pathway_two(2)];
    stats.varexp(3,:)=[rsq_ci_pathway_three(1) rsq_pathway_three rsq_ci_pathway_three(2)];
    stats.varexp(4,:)=[rsq_ci_pathway_four(1) rsq_pathway_four rsq_ci_pathway_four(2)];
    
end
%%

stats.target_partialcorr(1)=pathway_one_partial_corr;
stats.target_partialcorr(2)=pathway_two_partial_corr;
stats.target_partialcorr(3)=pathway_three_partial_corr;
stats.target_partialcorr(4)=pathway_four_partial_corr;


%%

stats.compare_correlations(1)=target_one_stats;
stats.compare_correlations(2)=target_two_stats;
