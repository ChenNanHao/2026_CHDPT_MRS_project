clear all;clc;close all
addpath('F:\Brain Connectivity Toolbox'); % Image load code, or use what ever you want
%% Read the master files (Single-voxel version)
T_summary{1} = readtable('E:\CHD_MRS\NH_MRS_results_summary_V2.xlsx','sheet','cchd_results'); 
T_summary{1} = T_summary{1}(2:end,:); 
T_summary{2} = readtable('E:\CHD_MRS\NH_MRS_results_summary_V2.xlsx','sheet','achd_results'); 
T_summary{2} = T_summary{2}(2:end,:); 
T_summary{3} = readtable('E:\CHD_MRS\NH_MRS_results_summary_V2.xlsx','sheet','avpt_results'); 
T_summary{3} = T_summary{3}(2:end,:); 
T_cchd_coor = T_summary{1}(:,[1 10:13]); 
T_cchd = T_summary{1}(:,[1 14:16]); 
T_cchd.Properties.VariableNames = {'Subject_name', 'TNAA absolute conc (mM)', 'TCho absolute conc (mM)', 'TCr absolute conc (mM)'};
T_cchd = addvars(T_cchd, table2array(T_cchd(:,2))./table2array(T_cchd(:,3)), 'NewVariableNames' , 'TNAA/TCho ratio', 'After', 'TCr absolute conc (mM)'); 
T_cchd = addvars(T_cchd, table2array(T_cchd(:,2))./table2array(T_cchd(:,4)), 'NewVariableNames' , 'TNAA/TCr ratio', 'After', 'TNAA/TCho ratio'); 
T_cchd = addvars(T_cchd, table2array(T_cchd(:,3))./table2array(T_cchd(:,4)), 'NewVariableNames' , 'TCho/TCr ratio', 'After', 'TNAA/TCr ratio'); 
T_achd_coor = T_summary{2}(:,[1 10:13]); 
T_achd = T_summary{2}(:,[1 14:16]); 
T_achd.Properties.VariableNames = {'Subject_name', 'TNAA absolute conc (mM)', 'TCho absolute conc (mM)', 'TCr absolute conc (mM)'};
T_achd = addvars(T_achd, table2array(T_achd(:,2))./table2array(T_achd(:,3)), 'NewVariableNames' , 'TNAA/TCho ratio', 'After', 'TCr absolute conc (mM)'); 
T_achd = addvars(T_achd, table2array(T_achd(:,2))./table2array(T_achd(:,4)), 'NewVariableNames' , 'TNAA/TCr ratio', 'After', 'TNAA/TCho ratio'); 
T_achd = addvars(T_achd, table2array(T_achd(:,3))./table2array(T_achd(:,4)), 'NewVariableNames' , 'TCho/TCr ratio', 'After', 'TNAA/TCr ratio'); 
T_avpt_coor = T_summary{3}(:,[1 10:13]); 
T_avpt = T_summary{3}(:,[1 14:16]); 
T_avpt.Properties.VariableNames = {'Subject_name', 'TNAA absolute conc (mM)', 'TCho absolute conc (mM)', 'TCr absolute conc (mM)'};
T_avpt = addvars(T_avpt, table2array(T_avpt(:,2))./table2array(T_avpt(:,3)), 'NewVariableNames' , 'TNAA/TCho ratio', 'After', 'TCr absolute conc (mM)'); 
T_avpt = addvars(T_avpt, table2array(T_avpt(:,2))./table2array(T_avpt(:,4)), 'NewVariableNames' , 'TNAA/TCr ratio', 'After', 'TNAA/TCho ratio'); 
T_avpt = addvars(T_avpt, table2array(T_avpt(:,3))./table2array(T_avpt(:,4)), 'NewVariableNames' , 'TCho/TCr ratio', 'After', 'TNAA/TCr ratio'); 

%*******************************************************************************
% Exclusion cirteria: metabolites with consistently lower   
%*******************************************************************************
% Mean and std of concentrations 
idx_cchd = find(~isnan(table2array(T_cchd_coor(:,2))) & ~isnan(table2array(T_cchd_coor(:,4)))); 
idx_achd = find(~isnan(table2array(T_achd_coor(:,2))) & ~isnan(table2array(T_achd_coor(:,4)))); 
idx_avpt = find(~isnan(table2array(T_avpt_coor(:,2))) & ~isnan(table2array(T_avpt_coor(:,4)))); 

T_cchd_include = T_cchd(idx_cchd,:); 
T_achd_include = T_achd(idx_achd,:); 
T_avpt_include = T_avpt(idx_avpt,:); 

for ii = 2:7
    mean_meta(ii-1,1) = mean(table2array(T_cchd(idx_cchd,ii)),'omitnan');
    mean_meta(ii-1,2) = mean(table2array(T_achd(idx_achd,ii)),'omitnan');
    mean_meta(ii-1,3) = mean(table2array(T_avpt(idx_avpt,ii)),'omitnan');
    std_meta(ii-1,1) = std(table2array(T_cchd(idx_cchd,ii)),'omitnan');
    std_meta(ii-1,2) = std(table2array(T_achd(idx_achd,ii)),'omitnan');
    std_meta(ii-1,3) = std(table2array(T_avpt(idx_avpt,ii)),'omitnan');
end

group_labels = [repmat({'CTL'}, length(idx_cchd), 1); repmat({'CHD'}, length(idx_achd), 1); repmat({'Preterm'}, length(idx_avpt), 1)];
idx_groups = [ones(length(idx_cchd), 1); ones(length(idx_achd), 1)*2; ones(length(idx_avpt), 1)*3];

age = [table2array(T_summary{1}(idx_cchd,2));table2array(T_summary{2}(idx_achd,2));table2array(T_summary{3}(idx_avpt,2))]; 
sex = cell2mat([table2array(T_summary{1}(idx_cchd,3));table2array(T_summary{2}(idx_achd,3));table2array(T_summary{3}(idx_avpt,3))]); 
BMI = [table2array(T_summary{1}(idx_cchd,4));table2array(T_summary{2}(idx_achd,4));table2array(T_summary{3}(idx_avpt,4))]; 
SES = [table2array(T_summary{1}(idx_cchd,5));table2array(T_summary{2}(idx_achd,5));table2array(T_summary{3}(idx_avpt,5))]; 
Brain_abnorm = [table2array(T_summary{1}(idx_cchd,18));table2array(T_summary{2}(idx_achd,18));table2array(T_summary{3}(idx_avpt,18))]; 
TBV = [table2array(T_summary{1}(idx_cchd,8));table2array(T_summary{2}(idx_achd,8));table2array(T_summary{3}(idx_avpt,8))]; 
WMV = [table2array(T_summary{1}(idx_cchd,9));table2array(T_summary{2}(idx_achd,9));table2array(T_summary{3}(idx_avpt,9))]; 
sex_idx = ones(length(sex),1); sex_idx(sex == 'F') = 2; 
sex_idx_cat = categorical(sex_idx);
Brain_abnorm_cat = categorical(Brain_abnorm);

%*******************************************************************************
% Clinical variables 
%*******************************************************************************
CHD_variables = T_summary{2}(idx_achd,31:45); 
PT_variables = T_summary{3}(idx_avpt,31:35); 



color_groups = [0 0.4470 0.7410;0.4660 0.6740 0.1880;0.6350 0.0780 0.1840]; 

%% Statistical analyses
%**************************************************************************
% Normality test (Shapiro-Wilk test)   
%**************************************************************************
for ii = 2:7
    [norm_H_meta(1,ii-1),norm_p_meta(1,ii-1)] = swtest(table2array(T_cchd(idx_cchd,ii)));
    [norm_H_meta(2,ii-1),norm_p_meta(2,ii-1)] = swtest(table2array(T_achd(idx_achd,ii)));
    [norm_H_meta(3,ii-1),norm_p_meta(3,ii-1)] = swtest(table2array(T_avpt(idx_avpt,ii)));
end

%********************************************************************************
% We start with a linear model to test differences in metabolite concentration 
%********************************************************************************
for ii = 2:7
    meta_conc = [table2array(T_cchd_include(:,ii));table2array(T_achd_include(:,ii));table2array(T_avpt_include(:,ii))];
    varNames = {'Metabolite_concentration','Groups','Sex','SES'};
    tbl = table(meta_conc, categorical(idx_groups), sex_idx_cat, SES,'VariableNames',varNames);
    T = [0 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1; 0 1 1 0]; 
    lm{ii-1} = fitlm(tbl, T);
end
disp(lm{5})

%**************************************************************************
% Sex-specific analyses (3 groups combined, ANCOVA controlling SES)
%**************************************************************************
male_idx_cchd = find(sex(idx_groups == 1) == 'M'); 
female_idx_cchd = find(sex(idx_groups == 1) == 'F'); 
male_idx_achd = find(sex(idx_groups == 2) == 'M'); 
female_idx_achd = find(sex(idx_groups == 2) == 'F'); 
male_idx_avpt = find(sex(idx_groups == 3) == 'M'); 
female_idx_avpt = find(sex(idx_groups == 3) == 'F'); 
SES_male = [SES(sex == 'M' & idx_groups == 1);SES(sex == 'M' & idx_groups == 2);SES(sex == 'M' & idx_groups == 3)]; 
SES_female = [SES(sex == 'F' & idx_groups == 1);SES(sex == 'F' & idx_groups == 2);SES(sex == 'F' & idx_groups == 3)]; 

idx_sex_cchd = [ones(length(male_idx_cchd), 1); ones(length(female_idx_cchd), 1)*2];
SES_sex_cchd = [SES(sex == 'M' & idx_groups == 1);SES(sex == 'F' & idx_groups == 1)]; 
idx_sex_achd = [ones(length(male_idx_achd), 1); ones(length(female_idx_achd), 1)*2];
SES_sex_achd = [SES(sex == 'M' & idx_groups == 2);SES(sex == 'F' & idx_groups == 2)]; 
idx_sex_avpt = [ones(length(male_idx_avpt), 1); ones(length(female_idx_avpt), 1)*2];
SES_sex_avpt = [SES(sex == 'M' & idx_groups == 3);SES(sex == 'F' & idx_groups == 3)]; 

% Metabolites: Male vs. Female 
clearvars lm_sex
for ii = 2:7
    meta_conc = [table2array(T_cchd_include(male_idx_cchd,ii)); table2array(T_cchd_include(female_idx_cchd,ii))];
    varNames = {'Metabolite_concentration','Groups', 'SES'};
    tbl = table(meta_conc, idx_sex_cchd, SES_sex_cchd,'VariableNames',varNames);
    lm_sex{ii-1} = fitlm(tbl, 'Metabolite_concentration ~ Groups + SES');
end

clearvars lm_sex
for ii = 2:7
    meta_conc = [table2array(T_achd_include(male_idx_achd,ii)); table2array(T_achd_include(female_idx_achd,ii))];
    varNames = {'Metabolite_concentration','Groups', 'SES'};
    tbl = table(meta_conc, idx_sex_achd, SES_sex_achd,'VariableNames',varNames);
    lm_sex{ii-1} = fitlm(tbl, 'Metabolite_concentration ~ Groups + SES');
end

clearvars lm_sex cohen_d_sex
for ii = 2:7
    meta_conc = [table2array(T_avpt_include(male_idx_avpt,ii)); table2array(T_avpt_include(female_idx_avpt,ii))];
    varNames = {'Metabolite_concentration','Groups', 'SES'};
    tbl = table(meta_conc, idx_sex_avpt, SES_sex_avpt,'VariableNames',varNames);
    lm_sex{ii-1} = fitlm(tbl, 'Metabolite_concentration ~ Groups + SES');
    cohen_d_sex(ii-1,:) = (mean(table2array(T_avpt_include(male_idx_avpt,ii)))-mean(table2array(T_avpt_include(female_idx_avpt,ii))))/std([table2array(T_avpt_include(male_idx_avpt,ii));table2array(T_avpt_include(female_idx_avpt,ii))]); 
end

disp(lm_sex{4})


% Male group comparisons
idx_Groups_male = [ones(length(male_idx_cchd), 1); ones(length(male_idx_achd), 1)*2; ones(length(male_idx_avpt), 1)*3];
idx_Groups_male = categorical(idx_Groups_male); 

for ii = 2:7
    meta_conc = [table2array(T_cchd_include(male_idx_cchd,ii));table2array(T_achd_include(male_idx_achd,ii));table2array(T_avpt_include(male_idx_avpt,ii))];
    varNames = {'Metabolite_concentration','Groups', 'SES'};
    tbl = table(meta_conc, idx_Groups_male, SES_male,'VariableNames',varNames);
    lm_male{ii-1} = fitlm(tbl, 'Metabolite_concentration ~ Groups + SES');
end

disp(lm_male{4})


% Female group comparisons
idx_Groups_female = [ones(length(female_idx_cchd), 1); ones(length(female_idx_achd), 1)*2; ones(length(female_idx_avpt), 1)*3];
idx_Groups_female = categorical(idx_Groups_female); 
for ii = 2:7
    meta_conc = [table2array(T_cchd_include(female_idx_cchd,ii));table2array(T_achd_include(female_idx_achd,ii));table2array(T_avpt_include(female_idx_avpt,ii))];
    varNames = {'Metabolite_concentration','Groups', 'SES'};
    tbl = table(meta_conc, idx_Groups_female, SES_female,'VariableNames',varNames);
    lm_female{ii-1} = fitlm(tbl, 'Metabolite_concentration ~ Groups + SES');
end
disp(lm_female{4})


%****************************************************************************************
% Associations with clinical variables (PT), Pearson/Spearman correlation with SES correction 
%****************************************************************************************
idx_PT_variables = ~isnan(table2array(PT_variables(:,1))); 
[r,p] = corr(table2array(T_avpt_include(idx_PT_variables,5)),table2array(PT_variables(idx_PT_variables,1)))

idx_PT_variables = ~isnan(table2array(PT_variables(:,2))); 
[r,p] = corr(table2array(T_avpt_include(idx_PT_variables,5)),table2array(PT_variables(idx_PT_variables,2)))

idx_PT_variables = ~isnan(table2array(PT_variables(:,3))); 
[r,p] = corr(table2array(T_avpt_include(idx_PT_variables,5)),table2array(PT_variables(idx_PT_variables,3)),'type','spearman')

idx_PT_variables = ~isnan(table2array(PT_variables(:,4))); 
[r,p] = corr(table2array(T_avpt_include(idx_PT_variables,5)),table2array(PT_variables(idx_PT_variables,4)))


%****************************************************************************************
% Comparisons of clinical variables (CHD), Mann-Whitney U-test  
%****************************************************************************************
for iii = 1:length(table2array(CHD_variables(:,14)))
    aaa = table2array(CHD_variables(iii,14)); 
    if strcmp(aaa{:},'no')
        BAS_before_sx_cat(iii,1) = 0; 
    elseif strcmp(aaa{:},'yes')
        BAS_before_sx_cat(iii,1) = 1; 
    else
        BAS_before_sx_cat(iii,1) = nan;
    end
end
CHD_variables = addvars(CHD_variables, BAS_before_sx_cat, 'NewVariableNames' , 'BAS_before_sx_cat', 'After', 'nb_BAS'); 

idx_CHD_variables = ~isnan(table2array(CHD_variables(:,7))); 
[r,p] = corr(table2array(T_achd_include(idx_CHD_variables,5)),table2array(CHD_variables(idx_CHD_variables,7)),'type','spearman')

idx_CHD_variables = ~isnan(table2array(CHD_variables(:,8))); 
[r,p] = corr(table2array(T_achd_include(idx_CHD_variables,5)),table2array(CHD_variables(idx_CHD_variables,8)),'type','spearman')

idx_CHD_variables = ~isnan(table2array(CHD_variables(:,9))); 
[r,p] = corr(table2array(T_achd_include(idx_CHD_variables,5)),table2array(CHD_variables(idx_CHD_variables,9)))

idx_CHD_variables = ~isnan(table2array(CHD_variables(:,10))); 
[r,p] = corr(table2array(T_achd_include(idx_CHD_variables,5)),table2array(CHD_variables(idx_CHD_variables,10)))

idx_CHD_variables = ~isnan(table2array(CHD_variables(:,11))); 
[r,p] = corr(table2array(T_achd_include(idx_CHD_variables,5)),table2array(CHD_variables(idx_CHD_variables,11)))

idx_CHD_variables = ~isnan(table2array(CHD_variables(:,13))); 
[r,p] = corr(table2array(T_achd_include(idx_CHD_variables,5)),table2array(CHD_variables(idx_CHD_variables,13)),'type','spearman')

idx_CHD_variables = ~isnan(table2array(CHD_variables(:,16))); 
[r,p] = corr(table2array(T_achd_include(idx_CHD_variables,5)),table2array(CHD_variables(idx_CHD_variables,16)),'type','spearman')



