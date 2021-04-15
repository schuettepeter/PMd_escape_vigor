%% load masks and fMRI data for analysis
% load(which('MPA_dataset_visual_pathway.mat'),'visual_pathway','S','COND','VAL','DAT','CeM', 'PAG','hythal');
load MPA_data_pHythal_PAG.mat

%% create ROIs - two sources and the same target
% PAG=load_atlas('Kragel2019PAG');
% PAG.dat=double(PAG.dat>0);
% 
% PAG=load_atlas('Kragel2019PAG');
% PAG.dat=double(PAG.dat>0);
% 
% subcort=load_atlas('CIT168');
% hythal = select_atlas_subset(subcort, {'Hythal'});

%% estimate models

stats = model_brain_pathway_v1(masked_dat,hythal,CeM,PAG,PAG,'Indices',S); 

%% full supplemental figure

create_figure('Pathway Specificity - PLS');
barplot_columns((stats.latent_correlations),'dolines','nofig','names',{'Hythal <-> PAG' 'Hythal <-> CeA Model' 'CeA <-> PAG Null' 'CeA <-> PAG'},'nostars','nofig','color',{[.1 0 .5],[.2 .2 .2]+[.1 0 .5],[.2 .2 .2]+[.7 .7 .7],[.7 .7 .7]},'noviolin');
ylabel 'Pearson Correlation'
xlabel ''
title ''
set(gca,'XTick',[])

legend_names={'HTH_H_T_H_-_P_A_G,PAG_H_T_H_-_P_A_G',...
'HTH_H_T_H_-_P_A_G,PAG_C_e_A_-_P_A_G',...
'CeA_C_e_A_-_P_A_G,PAG_H_T_H_-_P_A_G',...
'CeA_C_e_A_-_P_A_G,PAG_C_e_A_-_P_A_G'};

h=findobj(gca,'Type','Bar');
legend(h(1:4),fliplr(legend_names));


%% compute average expression for different conditions

for s=1:max(S) %for each subject
    for c=1:max(COND)  % for each condition
        mean_pathway_expression(s,c,:)=mean(stats.latent_timeseries(S==s & COND==c,:)); 
    end
    
     for v=1:max(VAL) %for each intensity level for negative images
        mean_pathway_expression_avoidance(s,v,:)=mean(stats.latent_timeseries(S==s & VAL==v & COND==4,:)); 
     end 
    
end

%% examine relationship between pathway expression and intensity of images

create_figure('Functional Relevance');

barplot_columns(mean_pathway_expression_avoidance(:,:,1),'dolines','nofig','nostars','noviolin','color',[.1 0 .8]); 
xlabel 'Stimulus Intensity'
ylabel 'Pathway Expression'
[~, p(1) ci st(1)]=ttest(mean_pathway_expression_avoidance(:,1,1),mean_pathway_expression_avoidance(:,2,1));
[~, p(2) ci st(2)]=ttest(mean_pathway_expression_avoidance(:,2,1),mean_pathway_expression_avoidance(:,3,1));
[~, p(3) ci st(3)]=ttest(mean_pathway_expression_avoidance(:,3,1),mean_pathway_expression_avoidance(:,4,1));

%% evaluate response for positive vs. negative images
create_figure('Functional Specificity - PLS ');
barplot_columns(squeeze(mean_pathway_expression(:,4:5,1)),'dolines','names',{'Negative Images' 'Positive Images'},'nofig','nostars','noviolin','color',{[.1 0 .8],[.5 .5 .5]}); 
ylabel 'Pathway Expression'
set(gca,'XTick',[])
xlabel ''
h=findobj(gca,'Type','Bar');
legend(h(1:2),'Positive Images','Negative Images');
[h p ci st]=ttest(squeeze(mean_pathway_expression(:,4,1)),squeeze(mean_pathway_expression(:,5,1)));

%% show model weights in hythal
zmap = threshold(stats.PLS_bootstrap_stats(1),.05,'FDR','k',10);
zmap.dat(~zmap.sig)=0;
zmap.fullpath='pathway_one_bpls_zmap.nii';
% write(zmap)
orthviews(zmap)

