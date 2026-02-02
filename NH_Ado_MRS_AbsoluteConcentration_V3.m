%% Output section (Single-voxel version, two voxels per participants)
%**************************************************************************
% We extract the outputs of MRS-DP for estimating tissue fractions and the
% outputs of Tarquin for unprocessed metabolites concentrations
% 
% The main purpose: With the selection of "nice" white matter voxels within
% the VOI, we ensure that our estimation of metabolites concentrations is
% in good quality. Considering partial volume effect contributed by grey
% matter and CSF may support us to more precisely measure the "Absolute"
% concentration representing white matter. 
% 
% Things need to be done beforehand: 
%       1. Tissue fractions estimated by MRS-DP
%                                      (TissueTypesForEachVoxel.txt)
% 
%       2. Metabolites concentrations estimated by Tarquin
%                                      (exported the results in txt format)
%
%**************************************************************************

%**************************************************************************
% We take the tissue fractions as the estimation of partial volume effect 
%**************************************************************************
clear all
Parent_folder = 'E:\CHD_MRS'; % Hard-coded now -- need to be changed 

subject_names = dir(fullfile(append(Parent_folder,'\Outputs\Tarquin_Concentration_Results\'),'achd_*')); 

for kk = 1:length(subject_names)
% clearvars -except T kk subject_names Parent_folder

TNAA_absolute_GG = []; 
TCho_absolute_GG = []; 
TCr_absolute_GG = []; 
WM_vm = [];

subject_name = subject_names(kk).name; 
subject_folder = append(Parent_folder,'\dcm\',subject_name(1:4),'\',subject_name); 
MRSDP_output_voxmap_dir = append(Parent_folder,'\Outputs\Scripts_Output\',subject_name); 
dicom_dir = append(subject_folder,'\DICOM'); 

if ~isfolder(MRSDP_output_voxmap_dir)
    error('No voxel maps generated from MRS-DP.py - Run the script first!')
end

if isempty(fullfile(dicom_dir,'*IM_*'))
    error('No voxel maps generated from MRS-DP.py - Run the script first!')
else
    dicomFiles = dir(fullfile(dicom_dir,'*IM_*'));
end

%**************************************************************************
% Check if the output folder is existed or not 
%**************************************************************************
MRSDP_output_voxmap_dir_new = append(Parent_folder,'\Outputs');
if ~isfolder(MRSDP_output_voxmap_dir_new)
    mkdir(MRSDP_output_voxmap_dir_new)
end

if ~isfolder(append(MRSDP_output_voxmap_dir_new,'\Scripts_Output'))
    mkdir(append(MRSDP_output_voxmap_dir_new,'\Scripts_Output'))
end
if ~isfolder(append(MRSDP_output_voxmap_dir_new,'\Scripts_Output\',subject_name))
    mkdir(append(MRSDP_output_voxmap_dir_new,'\Scripts_Output\',subject_name))
end

if ~isfolder(append(MRSDP_output_voxmap_dir_new,'\Tarquin_Concentration_Results'))
    mkdir(append(MRSDP_output_voxmap_dir_new,'\Tarquin_Concentration_Results'))
end
if ~isfolder(append(MRSDP_output_voxmap_dir_new,'\Tarquin_Concentration_Results\',subject_name))
    mkdir(append(MRSDP_output_voxmap_dir_new,'\Tarquin_Concentration_Results\',subject_name))
end

if ~isfolder(append(MRSDP_output_voxmap_dir_new,'\Absolute_Concentration_Results'))
    mkdir(append(MRSDP_output_voxmap_dir_new,'\Absolute_Concentration_Results'))
end
if ~isfolder(append(MRSDP_output_voxmap_dir_new,'\Absolute_Concentration_Results\',subject_name))
    mkdir(append(MRSDP_output_voxmap_dir_new,'\Absolute_Concentration_Results\',subject_name))
end

%**************************************************************************
% Read inputs (Tissue Fraction and Tarquin metabolite concentrations) 
%**************************************************************************
TissueFraction = fread(fopen(append(MRSDP_output_voxmap_dir_new,'\Scripts_Output\',subject_name,'\TissueTypesForEachVoxel.txt'), 'r')...
    , '*char')';

Metabolites_Tarquin = fread(fopen(append(MRSDP_output_voxmap_dir_new,'\Tarquin_Concentration_Results\',subject_name,'\',subject_name,'.txt'), 'r')...
    , '*char')';

% Tissue Fractions and voxel locations 
pattern = ['Voxels at position\s*\((\d+),(\d+)\)\s*' ...
           'Grey Matter Percentage:\s*([\d.]+)\s*' ...
           'White Matter Percentage:\s*([\d.]+)\s*' ...
           'CSF Percentage:\s*([\d.]+)'];
tokens = regexp(TissueFraction, pattern, 'tokens');

n = length(tokens);
TissueFraction_MRSDP = zeros(n, 5); % columns: Row, Col, GM, WM, CSF
for i = 1:n
    vals = str2double(tokens{i});
    TissueFraction_MRSDP(i, :) = vals;
end
voxel_position_MRSDP = TissueFraction_MRSDP(:,1:2); 
% Find corner coordinates (Flip the row and column so it matches with Tarquin outputs
voxel_corner_MRSDP = [max(voxel_position_MRSDP(:,2)) max(voxel_position_MRSDP(:,1)); ...
    min(voxel_position_MRSDP(:,2)) max(voxel_position_MRSDP(:,1)); ...
    max(voxel_position_MRSDP(:,2)) min(voxel_position_MRSDP(:,1)); ...
    min(voxel_position_MRSDP(:,2)) min(voxel_position_MRSDP(:,1))]; 

% Tarquin 
% Find all voxel coordinate
coordExpr = 'Row\s*:\s*(\d+),\s*Col\s*:\s*(\d+),\s*Slice\s*:\s*(\d+)';
[starts, ends, tokens] = regexp(Metabolites_Tarquin, coordExpr, 'start', 'end', 'tokens');
voxel_position_tarquin = []; 
for zz = 1:size(tokens,2)
    voxel_position_tarquin(zz,:) = str2double(tokens{zz}(1:2)); 
end
% Find corner coordinates 
voxel_corner_tarquin = [min(voxel_position_tarquin(:,1)) min(voxel_position_tarquin(:,2)); ...
    max(voxel_position_tarquin(:,1)) min(voxel_position_tarquin(:,2)); ...
    min(voxel_position_tarquin(:,1)) max(voxel_position_tarquin(:,2)); ...
    max(voxel_position_tarquin(:,1)) max(voxel_position_tarquin(:,2))]; 

% Relationship between MRSDP and Tarquin: 
voxel_position_tarquin_trans = [max(voxel_position_tarquin(:,1))-voxel_position_tarquin(:,1) max(voxel_position_tarquin(:,2))-voxel_position_tarquin(:,2)];  

% Find corresponding voxels and their tissue fractions
TissueFraction_MRSDP_selected = [];
for zz = 1:size(voxel_position_tarquin_trans,1)
    x = find(TissueFraction_MRSDP(:,2) == voxel_position_tarquin_trans(zz,1) & TissueFraction_MRSDP(:,1) == voxel_position_tarquin_trans(zz,2)); 
    TissueFraction_MRSDP_selected(zz,:) = TissueFraction_MRSDP(x,3:end); 
end

%**************************************************************************
% Correction of Tissue volume fractions into mole tissue fraction
%**************************************************************************
% Density of GM, WM and CSF from Gasparovic 2006 and Ernst 1993
TissueFraction_MRSDP_selected_corr = []; 
TissueFraction_MRSDP_selected_corr(:,1) = 0.78*TissueFraction_MRSDP_selected(:,1); 
TissueFraction_MRSDP_selected_corr(:,2) = 0.65*TissueFraction_MRSDP_selected(:,2); 
TissueFraction_MRSDP_selected_corr(:,3) = 0.97*TissueFraction_MRSDP_selected(:,3); 
TissueFraction_final = TissueFraction_MRSDP_selected_corr./sum(TissueFraction_MRSDP_selected_corr,2); 

%**************************************************************************
% We extract the metabolites concentration (mM) from Tarquin (UNPROCESSED!)
%**************************************************************************
nVoxels = length(tokens);
TNAA_data = [];
TCho_data = [];
TCr_data = [];

FWHM_data = [];
SNR_data = [];

for i = 1:nVoxels
    row = str2double(tokens{i}{1});
    col = str2double(tokens{i}{2});
    slice = str2double(tokens{i}{3});
    % Define block boundaries
    startIdx = ends(i);
    if i < nVoxels
        endIdx = starts(i+1) - 1;
    else
        endIdx = length(Metabolites_Tarquin);
    end
    block = Metabolites_Tarquin(startIdx:endIdx);
    % Extract TNAA data using regex
    tnaaExpr = 'TNAA\s+([\d\.e\+\-]+)\s+([\d\.#IOe\+\-]+)\s+([\d\.e\+\-]+)';
    match = regexp(block, tnaaExpr, 'tokens');
    if ~isempty(match)
        conc = str2double(match{1}{1});
        perc_sd = str2double(match{1}{2});
        sd = str2double(match{1}{3});
        TNAA_data = [TNAA_data; conc, perc_sd, sd];
    end
    % Extract TCho data using regex
    tchoExpr = 'TCho\s+([\d\.e\+\-]+)\s+([\d\.#IOe\+\-]+)\s+([\d\.e\+\-]+)';
    match = regexp(block, tchoExpr, 'tokens');
    if ~isempty(match)
        conc = str2double(match{1}{1});
        perc_sd = str2double(match{1}{2});
        sd = str2double(match{1}{3});
        TCho_data = [TCho_data; conc, perc_sd, sd];
    end
    % Extract TCr data using regex
    tcrExpr = 'TCr\s+([\d\.e\+\-]+)\s+([\d\.#IOe\+\-]+)\s+([\d\.e\+\-]+)';
    match = regexp(block, tcrExpr, 'tokens');
    if ~isempty(match)
        conc = str2double(match{1}{1});
        perc_sd = str2double(match{1}{2});
        sd = str2double(match{1}{3});
        TCr_data = [TCr_data; conc, perc_sd, sd];
    end
    % Extract data quality (FWHM)
    FWHMExpr = 'Metab FWHM \(PPM\)\s*:\s*([\d\.]+)';
    match = regexp(block, FWHMExpr, 'tokens');
    if ~isempty(match)
        fwhm = str2double(match{1}{1});
        FWHM_data = [FWHM_data; fwhm];
    end
    % Extract data quality (SNR)
    SNRExpr = 'SNR residual\s*:\s*([\d\.]+)';
    match = regexp(block, SNRExpr, 'tokens');
    if ~isempty(match)
        snr = str2double(match{1}{1});
        SNR_data = [SNR_data; snr];
    end
end

%**************************************************************************
% We calculate the absolute concentration of each metabolites 
%**************************************************************************
% Units of TE, TR, T1 and T2 relaxation: ms 
% Units of water concentration: mmol/L 
Parameters.conc_H20 = 35880;
Parameters.att_H20 = 0.7; 
Parameters.ProtonWater = 55556; 
% Parameters.H20_WM = 36100; Parameters.H20_GM = 43300; Parameters.H20_CSF = 53800; 
Parameters.T2_WM = 79.2; Parameters.T2_GM = 110; Parameters.T2_CSF = 503; 
Parameters.T1_WM = 832; Parameters.T1_GM = 1331; Parameters.T1_CSF = 3817; 
Parameters.TE = 144; Parameters.TR = 2000; 

% T1 and T2 relaxation times were based on this study: Posse et al. 2007
Parameters.T2_TNAA_GM = 247; Parameters.T2_TCho_GM = 222; Parameters.T2_TCr_GM = 162; 
Parameters.T2_TNAA_WM = 301; Parameters.T2_TCho_WM = 222; Parameters.T2_TCr_WM = 178; 
Parameters.T1_TNAA_GM = 1470; Parameters.T1_TCho_GM = 1250; Parameters.T1_TCr_GM = 1339; 
Parameters.T1_TNAA_WM = 1560; Parameters.T1_TCho_WM = 1210; Parameters.T1_TCr_WM = 1400; 

idx_WM = find(TissueFraction_MRSDP_selected(:,2) >= 0.5);
idx_MRS = find(FWHM_data(idx_WM) < 0.1 & SNR_data(idx_WM) > 5);
nVoxels_final = size(idx_WM(idx_MRS),1); 
for i = 1:size(idx_WM(idx_MRS),1)
    ii = idx_WM(idx_MRS(i)); 
    AttH20GM = exp(-Parameters.TE/Parameters.T2_GM) * (1-exp(-Parameters.TR/Parameters.T1_GM)); 
    AttH20WM = exp(-Parameters.TE/Parameters.T2_WM) * (1-exp(-Parameters.TR/Parameters.T1_WM)); 
    AttH20CSF = exp(-Parameters.TE/Parameters.T2_CSF) * (1-exp(-Parameters.TR/Parameters.T1_CSF)); 
    RelaxH20 = (TissueFraction_MRSDP_selected(ii,1)*AttH20GM)+(TissueFraction_MRSDP_selected(ii,2)*AttH20WM)+(TissueFraction_MRSDP_selected(ii,3)*AttH20CSF);
    FracH20 = (0.81*TissueFraction_MRSDP_selected(ii,1))+(0.71*TissueFraction_MRSDP_selected(ii,2))+(TissueFraction_MRSDP_selected(ii,3));

    AttNAAGM=(1-exp(-Parameters.TR/Parameters.T1_TNAA_GM))*exp(-Parameters.TE/Parameters.T2_TNAA_GM);
    AttNAAWM=(1-exp(-Parameters.TR/Parameters.T1_TNAA_WM))*exp(-Parameters.TE/Parameters.T2_TNAA_WM);
    RelaxNAA=((TissueFraction_MRSDP_selected(ii,1)*AttNAAGM)+(TissueFraction_MRSDP_selected(ii,2)*AttNAAWM))/(1-TissueFraction_MRSDP_selected(ii,3));
    FracNAA = 1-TissueFraction_MRSDP_selected(ii,3); 
    TNAA_absolute_GG(i) = TNAA_data(ii,1)/Parameters.att_H20*(FracH20/FracNAA)*(Parameters.ProtonWater/Parameters.conc_H20)*(RelaxH20/RelaxNAA);

    AttChoGM=(1-exp(-Parameters.TR/Parameters.T1_TCho_GM))*exp(-Parameters.TE/Parameters.T2_TCho_GM);
    AttChoWM=(1-exp(-Parameters.TR/Parameters.T1_TCho_WM))*exp(-Parameters.TE/Parameters.T2_TCho_WM);
    RelaxCho=((TissueFraction_MRSDP_selected(ii,1)*AttChoGM)+(TissueFraction_MRSDP_selected(ii,2)*AttChoWM))/(1-TissueFraction_MRSDP_selected(ii,3));
    FracCho = 1-TissueFraction_MRSDP_selected(ii,3); 
    TCho_absolute_GG(i) = TCho_data(ii,1)/Parameters.att_H20*(RelaxH20/RelaxCho)*((FracH20/FracCho)*Parameters.ProtonWater)/Parameters.conc_H20;

    AttCrGM=(1-exp(-Parameters.TR/Parameters.T1_TCr_GM))*exp(-Parameters.TE/Parameters.T2_TCr_GM);
    AttCrWM=(1-exp(-Parameters.TR/Parameters.T1_TCr_WM))*exp(-Parameters.TE/Parameters.T2_TCr_WM);
    RelaxCr=((TissueFraction_MRSDP_selected(ii,1)*AttCrGM)+(TissueFraction_MRSDP_selected(ii,2)*AttCrWM))/(1-TissueFraction_MRSDP_selected(ii,3));
    FracCr= 1-TissueFraction_MRSDP_selected(ii,3); 
    TCr_absolute_GG(i) = TCr_data(ii,1)/Parameters.att_H20*(RelaxH20/RelaxCr)*((FracH20/FracCr)*Parameters.ProtonWater)/Parameters.conc_H20;

    WM_vm(i) = TissueFraction_MRSDP_selected(ii,2); 
end

% Select a voxel with the highest WM probability in each brain hemisphere (2 voxels in total) 
idx_WM_vm = zeros(length(TissueFraction_MRSDP_selected),1);
idx_WM_vm(idx_WM(idx_MRS)) = 1; 

idx_midline = length(unique(TissueFraction_MRSDP(:,1)))/2; 
rh_include = voxel_position_tarquin_trans(:,2)>idx_midline; 
lh_include = voxel_position_tarquin_trans(:,2)<=idx_midline; 
TissueFraction_MRSDP_selected_rh = TissueFraction_MRSDP_selected((rh_include & idx_WM_vm),2); 
TissueFraction_MRSDP_selected_lh = TissueFraction_MRSDP_selected((lh_include & idx_WM_vm),2); 
if isempty(TissueFraction_MRSDP_selected_rh)
    rh_WMmax_voxel = [];
    rh_WMmax_voxel_conc_idx = [];
    aa = nan; 
    aa1 = nan; 
else
    rh_WMmax_voxel = find((TissueFraction_MRSDP(:,4) == max(TissueFraction_MRSDP_selected_rh)) == 1);
    rh_WMmax_voxel_conc_idx = find(WM_vm == max(TissueFraction_MRSDP_selected_rh)); 
    aa = TissueFraction_MRSDP(rh_WMmax_voxel,1); 
    aa1 = TissueFraction_MRSDP(rh_WMmax_voxel,2); 
end

if isempty(TissueFraction_MRSDP_selected_lh)
    lh_WMmax_voxel = [];
    lh_WMmax_voxel_conc_idx = [];
    bb = nan; 
    bb1 = nan; 
else
    lh_WMmax_voxel = find((TissueFraction_MRSDP(:,4) == max(TissueFraction_MRSDP_selected_lh)) == 1);
    lh_WMmax_voxel_conc_idx = find(WM_vm == max(TissueFraction_MRSDP_selected_lh)); 
    bb = TissueFraction_MRSDP(lh_WMmax_voxel,1); 
    bb1 = TissueFraction_MRSDP(lh_WMmax_voxel,2); 
end

% T(kk,:) = table({subject_name},nVoxels_final,mean(WM_vm),mean(TNAA_absolute_GG),mean(TCho_absolute_GG),mean(TCr_absolute_GG),...
%     'VariableNames', {'Subject_name','Number of voxels', 'WM volume fraction','TNAA absolute conc (mM)', 'TCho absolute conc (mM)', 'TCr absolute conc (mM)'}); 

T(kk,:) = table({subject_name},aa,aa1,bb,bb1,...
    mean(TNAA_absolute_GG([rh_WMmax_voxel_conc_idx lh_WMmax_voxel_conc_idx])),mean(TCho_absolute_GG([rh_WMmax_voxel_conc_idx lh_WMmax_voxel_conc_idx])),mean(TCr_absolute_GG([rh_WMmax_voxel_conc_idx lh_WMmax_voxel_conc_idx])),...
    'VariableNames', {'Subject_name','Voxel_location_rh_x','Voxel_location_rh_y','Voxel_location_lh_x','Voxel_location_lh_y', ...
    'TNAA absolute conc (mM)', 'TCho absolute conc (mM)', 'TCr absolute conc (mM)'}); 

% if nVoxels_final ~= 0
%     T(kk,:) = table({subject_name},nVoxels_final,max(WM_vm),TNAA_absolute_GG(WM_vm == max(WM_vm)),TCho_absolute_GG(WM_vm == max(WM_vm)),TCr_absolute_GG(WM_vm == max(WM_vm)),...
%          'VariableNames', {'Subject_name','Number of voxels', 'WM volume fraction','TNAA absolute conc (mM)', 'TCho absolute conc (mM)', 'TCr absolute conc (mM)'}); 
% else
%     T(kk,:) = table({subject_name},nVoxels_final,nan,nan,nan,nan,...
%          'VariableNames', {'Subject_name','Number of voxels', 'WM volume fraction','TNAA absolute conc (mM)', 'TCho absolute conc (mM)', 'TCr absolute conc (mM)'}); 
% end

disp(append('Finished subject ',subject_name))
end

figure
violinplot([table2array(T(:,6)),table2array(T(:,7)),table2array(T(:,8))],{'TNAA','TCho','TCr'})
ylabel('Concentration (mM)','fontweight','bold')
title('CHD','fontweight','bold')
set(gca,'linewidth',2,'fontname','Times New Romans','fontsize',16,'FontWeight','bold')

T_achd = T; 

T_achd = addvars(T_achd, table2array(T_achd(:,6))./table2array(T_achd(:,7)), 'NewVariableNames' , 'TNAA/TCho ratio', 'After', 'TCr absolute conc (mM)'); 
T_achd = addvars(T_achd, table2array(T_achd(:,6))./table2array(T_achd(:,8)), 'NewVariableNames' , 'TNAA/TCr ratio', 'After', 'TNAA/TCho ratio'); 
T_achd = addvars(T_achd, table2array(T_achd(:,7))./table2array(T_achd(:,8)), 'NewVariableNames' , 'TCho/TCr ratio', 'After', 'TNAA/TCr ratio'); 

writetable(T_achd,'E:\CHD_MRS\Outputs\Absolute_Concentration_Results\achd_results_2voxels.xlsx')

%% Controls 
clearvars T
subject_names = dir(fullfile(append(Parent_folder,'\Outputs\Tarquin_Concentration_Results\'),'cchd_*')); 

for kk = 1:length(subject_names)
% clearvars -except T kk subject_names Parent_folder

TNAA_absolute_GG = []; 
TCho_absolute_GG = []; 
TCr_absolute_GG = []; 
WM_vm = [];

subject_name = subject_names(kk).name; 
subject_folder = append(Parent_folder,'\dcm\',subject_name(1:4),'\',subject_name); 
MRSDP_output_voxmap_dir = append(Parent_folder,'\Outputs\Scripts_Output\',subject_name); 
dicom_dir = append(subject_folder,'\DICOM'); 

if ~isfolder(MRSDP_output_voxmap_dir)
    error('No voxel maps generated from MRS-DP.py - Run the script first!')
end

if isempty(fullfile(dicom_dir,'*IM_*'))
    error('No voxel maps generated from MRS-DP.py - Run the script first!')
else
    dicomFiles = dir(fullfile(dicom_dir,'*IM_*'));
end

%**************************************************************************
% Check if the output folder is existed or not 
%**************************************************************************
MRSDP_output_voxmap_dir_new = append(Parent_folder,'\Outputs');
if ~isfolder(MRSDP_output_voxmap_dir_new)
    mkdir(MRSDP_output_voxmap_dir_new)
end

if ~isfolder(append(MRSDP_output_voxmap_dir_new,'\Scripts_Output'))
    mkdir(append(MRSDP_output_voxmap_dir_new,'\Scripts_Output'))
end
if ~isfolder(append(MRSDP_output_voxmap_dir_new,'\Scripts_Output\',subject_name))
    mkdir(append(MRSDP_output_voxmap_dir_new,'\Scripts_Output\',subject_name))
end

if ~isfolder(append(MRSDP_output_voxmap_dir_new,'\Tarquin_Concentration_Results'))
    mkdir(append(MRSDP_output_voxmap_dir_new,'\Tarquin_Concentration_Results'))
end
if ~isfolder(append(MRSDP_output_voxmap_dir_new,'\Tarquin_Concentration_Results\',subject_name))
    mkdir(append(MRSDP_output_voxmap_dir_new,'\Tarquin_Concentration_Results\',subject_name))
end

if ~isfolder(append(MRSDP_output_voxmap_dir_new,'\Absolute_Concentration_Results'))
    mkdir(append(MRSDP_output_voxmap_dir_new,'\Absolute_Concentration_Results'))
end
if ~isfolder(append(MRSDP_output_voxmap_dir_new,'\Absolute_Concentration_Results\',subject_name))
    mkdir(append(MRSDP_output_voxmap_dir_new,'\Absolute_Concentration_Results\',subject_name))
end

%**************************************************************************
% Read inputs (Tissue Fraction and Tarquin metabolite concentrations) 
%**************************************************************************
TissueFraction = fread(fopen(append(MRSDP_output_voxmap_dir_new,'\Scripts_Output\',subject_name,'\TissueTypesForEachVoxel.txt'), 'r')...
    , '*char')';

Metabolites_Tarquin = fread(fopen(append(MRSDP_output_voxmap_dir_new,'\Tarquin_Concentration_Results\',subject_name,'\',subject_name,'.txt'), 'r')...
    , '*char')';

% Tissue Fractions and voxel locations 
pattern = ['Voxels at position\s*\((\d+),(\d+)\)\s*' ...
           'Grey Matter Percentage:\s*([\d.]+)\s*' ...
           'White Matter Percentage:\s*([\d.]+)\s*' ...
           'CSF Percentage:\s*([\d.]+)'];
tokens = regexp(TissueFraction, pattern, 'tokens');

n = length(tokens);
TissueFraction_MRSDP = zeros(n, 5); % columns: Row, Col, GM, WM, CSF
for i = 1:n
    vals = str2double(tokens{i});
    TissueFraction_MRSDP(i, :) = vals;
end
voxel_position_MRSDP = TissueFraction_MRSDP(:,1:2); 
% Find corner coordinates (Flip the row and column so it matches with Tarquin outputs
voxel_corner_MRSDP = [max(voxel_position_MRSDP(:,2)) max(voxel_position_MRSDP(:,1)); ...
    min(voxel_position_MRSDP(:,2)) max(voxel_position_MRSDP(:,1)); ...
    max(voxel_position_MRSDP(:,2)) min(voxel_position_MRSDP(:,1)); ...
    min(voxel_position_MRSDP(:,2)) min(voxel_position_MRSDP(:,1))]; 

% Tarquin 
% Find all voxel coordinate
coordExpr = 'Row\s*:\s*(\d+),\s*Col\s*:\s*(\d+),\s*Slice\s*:\s*(\d+)';
[starts, ends, tokens] = regexp(Metabolites_Tarquin, coordExpr, 'start', 'end', 'tokens');
voxel_position_tarquin = []; 
for zz = 1:size(tokens,2)
    voxel_position_tarquin(zz,:) = str2double(tokens{zz}(1:2)); 
end
% Find corner coordinates 
voxel_corner_tarquin = [min(voxel_position_tarquin(:,1)) min(voxel_position_tarquin(:,2)); ...
    max(voxel_position_tarquin(:,1)) min(voxel_position_tarquin(:,2)); ...
    min(voxel_position_tarquin(:,1)) max(voxel_position_tarquin(:,2)); ...
    max(voxel_position_tarquin(:,1)) max(voxel_position_tarquin(:,2))]; 

% Relationship between MRSDP and Tarquin: 
voxel_position_tarquin_trans = [max(voxel_position_tarquin(:,1))-voxel_position_tarquin(:,1) max(voxel_position_tarquin(:,2))-voxel_position_tarquin(:,2)];  

% Find corresponding voxels and their tissue fractions
TissueFraction_MRSDP_selected = [];
for zz = 1:size(voxel_position_tarquin_trans,1)
    x = find(TissueFraction_MRSDP(:,2) == voxel_position_tarquin_trans(zz,1) & TissueFraction_MRSDP(:,1) == voxel_position_tarquin_trans(zz,2)); 
    TissueFraction_MRSDP_selected(zz,:) = TissueFraction_MRSDP(x,3:end); 
end

%**************************************************************************
% Correction of Tissue volume fractions into mole tissue fraction
%**************************************************************************
% Density of GM, WM and CSF from Gasparovic 2006 and Ernst 1993
TissueFraction_MRSDP_selected_corr = []; 
TissueFraction_MRSDP_selected_corr(:,1) = 0.78*TissueFraction_MRSDP_selected(:,1); 
TissueFraction_MRSDP_selected_corr(:,2) = 0.65*TissueFraction_MRSDP_selected(:,2); 
TissueFraction_MRSDP_selected_corr(:,3) = 0.97*TissueFraction_MRSDP_selected(:,3); 
TissueFraction_final = TissueFraction_MRSDP_selected_corr./sum(TissueFraction_MRSDP_selected_corr,2); 

%**************************************************************************
% We extract the metabolites concentration (mM) from Tarquin (UNPROCESSED!)
%**************************************************************************
nVoxels = length(tokens);
TNAA_data = [];
TCho_data = [];
TCr_data = [];

FWHM_data = [];
SNR_data = [];

for i = 1:nVoxels
    row = str2double(tokens{i}{1});
    col = str2double(tokens{i}{2});
    slice = str2double(tokens{i}{3});
    % Define block boundaries
    startIdx = ends(i);
    if i < nVoxels
        endIdx = starts(i+1) - 1;
    else
        endIdx = length(Metabolites_Tarquin);
    end
    block = Metabolites_Tarquin(startIdx:endIdx);
    % Extract TNAA data using regex
    tnaaExpr = 'TNAA\s+([\d\.e\+\-]+)\s+([\d\.#IOe\+\-]+)\s+([\d\.e\+\-]+)';
    match = regexp(block, tnaaExpr, 'tokens');
    if ~isempty(match)
        conc = str2double(match{1}{1});
        perc_sd = str2double(match{1}{2});
        sd = str2double(match{1}{3});
        TNAA_data = [TNAA_data; conc, perc_sd, sd];
    end
    % Extract TCho data using regex
    tchoExpr = 'TCho\s+([\d\.e\+\-]+)\s+([\d\.#IOe\+\-]+)\s+([\d\.e\+\-]+)';
    match = regexp(block, tchoExpr, 'tokens');
    if ~isempty(match)
        conc = str2double(match{1}{1});
        perc_sd = str2double(match{1}{2});
        sd = str2double(match{1}{3});
        TCho_data = [TCho_data; conc, perc_sd, sd];
    end
    % Extract TCr data using regex
    tcrExpr = 'TCr\s+([\d\.e\+\-]+)\s+([\d\.#IOe\+\-]+)\s+([\d\.e\+\-]+)';
    match = regexp(block, tcrExpr, 'tokens');
    if ~isempty(match)
        conc = str2double(match{1}{1});
        perc_sd = str2double(match{1}{2});
        sd = str2double(match{1}{3});
        TCr_data = [TCr_data; conc, perc_sd, sd];
    end
    % Extract data quality (FWHM)
    FWHMExpr = 'Metab FWHM \(PPM\)\s*:\s*([\d\.]+)';
    match = regexp(block, FWHMExpr, 'tokens');
    if ~isempty(match)
        fwhm = str2double(match{1}{1});
        FWHM_data = [FWHM_data; fwhm];
    end
    % Extract data quality (SNR)
    SNRExpr = 'SNR residual\s*:\s*([\d\.]+)';
    match = regexp(block, SNRExpr, 'tokens');
    if ~isempty(match)
        snr = str2double(match{1}{1});
        SNR_data = [SNR_data; snr];
    end
end

%**************************************************************************
% We calculate the absolute concentration of each metabolites 
%**************************************************************************
% Units of TE, TR, T1 and T2 relaxation: ms 
% Units of water concentration: mmol/L 
Parameters.conc_H20 = 35880;
Parameters.att_H20 = 0.7; 
Parameters.ProtonWater = 55556; 
% Parameters.H20_WM = 36100; Parameters.H20_GM = 43300; Parameters.H20_CSF = 53800; 
Parameters.T2_WM = 79.2; Parameters.T2_GM = 110; Parameters.T2_CSF = 503; 
Parameters.T1_WM = 832; Parameters.T1_GM = 1331; Parameters.T1_CSF = 3817; 
Parameters.TE = 144; Parameters.TR = 2000; 

% T1 and T2 relaxation times were based on this study: Posse et al. 2007
Parameters.T2_TNAA_GM = 247; Parameters.T2_TCho_GM = 222; Parameters.T2_TCr_GM = 162; 
Parameters.T2_TNAA_WM = 301; Parameters.T2_TCho_WM = 222; Parameters.T2_TCr_WM = 178; 
Parameters.T1_TNAA_GM = 1470; Parameters.T1_TCho_GM = 1250; Parameters.T1_TCr_GM = 1339; 
Parameters.T1_TNAA_WM = 1560; Parameters.T1_TCho_WM = 1210; Parameters.T1_TCr_WM = 1400; 

idx_WM = find(TissueFraction_MRSDP_selected(:,2) >= 0.5);
idx_MRS = find(FWHM_data(idx_WM) < 0.1 & SNR_data(idx_WM) > 5);
nVoxels_final = size(idx_WM(idx_MRS),1); 
for i = 1:size(idx_WM(idx_MRS),1)
    ii = idx_WM(idx_MRS(i)); 
    AttH20GM = exp(-Parameters.TE/Parameters.T2_GM) * (1-exp(-Parameters.TR/Parameters.T1_GM)); 
    AttH20WM = exp(-Parameters.TE/Parameters.T2_WM) * (1-exp(-Parameters.TR/Parameters.T1_WM)); 
    AttH20CSF = exp(-Parameters.TE/Parameters.T2_CSF) * (1-exp(-Parameters.TR/Parameters.T1_CSF)); 
    RelaxH20 = (TissueFraction_MRSDP_selected(ii,1)*AttH20GM)+(TissueFraction_MRSDP_selected(ii,2)*AttH20WM)+(TissueFraction_MRSDP_selected(ii,3)*AttH20CSF);
    FracH20 = (0.81*TissueFraction_MRSDP_selected(ii,1))+(0.71*TissueFraction_MRSDP_selected(ii,2))+(TissueFraction_MRSDP_selected(ii,3));

    AttNAAGM=(1-exp(-Parameters.TR/Parameters.T1_TNAA_GM))*exp(-Parameters.TE/Parameters.T2_TNAA_GM);
    AttNAAWM=(1-exp(-Parameters.TR/Parameters.T1_TNAA_WM))*exp(-Parameters.TE/Parameters.T2_TNAA_WM);
    RelaxNAA=((TissueFraction_MRSDP_selected(ii,1)*AttNAAGM)+(TissueFraction_MRSDP_selected(ii,2)*AttNAAWM))/(1-TissueFraction_MRSDP_selected(ii,3));
    FracNAA = 1-TissueFraction_MRSDP_selected(ii,3); 
    TNAA_absolute_GG(i) = TNAA_data(ii,1)/Parameters.att_H20*(FracH20/FracNAA)*(Parameters.ProtonWater/Parameters.conc_H20)*(RelaxH20/RelaxNAA);

    AttChoGM=(1-exp(-Parameters.TR/Parameters.T1_TCho_GM))*exp(-Parameters.TE/Parameters.T2_TCho_GM);
    AttChoWM=(1-exp(-Parameters.TR/Parameters.T1_TCho_WM))*exp(-Parameters.TE/Parameters.T2_TCho_WM);
    RelaxCho=((TissueFraction_MRSDP_selected(ii,1)*AttChoGM)+(TissueFraction_MRSDP_selected(ii,2)*AttChoWM))/(1-TissueFraction_MRSDP_selected(ii,3));
    FracCho = 1-TissueFraction_MRSDP_selected(ii,3); 
    TCho_absolute_GG(i) = TCho_data(ii,1)/Parameters.att_H20*(RelaxH20/RelaxCho)*((FracH20/FracCho)*Parameters.ProtonWater)/Parameters.conc_H20;

    AttCrGM=(1-exp(-Parameters.TR/Parameters.T1_TCr_GM))*exp(-Parameters.TE/Parameters.T2_TCr_GM);
    AttCrWM=(1-exp(-Parameters.TR/Parameters.T1_TCr_WM))*exp(-Parameters.TE/Parameters.T2_TCr_WM);
    RelaxCr=((TissueFraction_MRSDP_selected(ii,1)*AttCrGM)+(TissueFraction_MRSDP_selected(ii,2)*AttCrWM))/(1-TissueFraction_MRSDP_selected(ii,3));
    FracCr= 1-TissueFraction_MRSDP_selected(ii,3); 
    TCr_absolute_GG(i) = TCr_data(ii,1)/Parameters.att_H20*(RelaxH20/RelaxCr)*((FracH20/FracCr)*Parameters.ProtonWater)/Parameters.conc_H20;

    WM_vm(i) = TissueFraction_MRSDP_selected(ii,2); 
end

% Select a voxel with the highest WM probability in each brain hemisphere (2 voxels in total) 
idx_WM_vm = zeros(length(TissueFraction_MRSDP_selected),1);
idx_WM_vm(idx_WM(idx_MRS)) = 1; 

idx_midline = length(unique(TissueFraction_MRSDP(:,1)))/2; 
rh_include = voxel_position_tarquin_trans(:,2)>idx_midline; 
lh_include = voxel_position_tarquin_trans(:,2)<=idx_midline; 
TissueFraction_MRSDP_selected_rh = TissueFraction_MRSDP_selected((rh_include & idx_WM_vm),2); 
TissueFraction_MRSDP_selected_lh = TissueFraction_MRSDP_selected((lh_include & idx_WM_vm),2); 
if isempty(TissueFraction_MRSDP_selected_rh)
    rh_WMmax_voxel = [];
    rh_WMmax_voxel_conc_idx = [];
    aa = nan; 
    aa1 = nan; 
else
    rh_WMmax_voxel = find((TissueFraction_MRSDP(:,4) == max(TissueFraction_MRSDP_selected_rh)) == 1);
    rh_WMmax_voxel_conc_idx = find(WM_vm == max(TissueFraction_MRSDP_selected_rh)); 
    aa = TissueFraction_MRSDP(rh_WMmax_voxel,1); 
    aa1 = TissueFraction_MRSDP(rh_WMmax_voxel,2); 
end

if isempty(TissueFraction_MRSDP_selected_lh)
    lh_WMmax_voxel = [];
    lh_WMmax_voxel_conc_idx = [];
    bb = nan; 
    bb1 = nan; 
else
    lh_WMmax_voxel = find((TissueFraction_MRSDP(:,4) == max(TissueFraction_MRSDP_selected_lh)) == 1);
    lh_WMmax_voxel_conc_idx = find(WM_vm == max(TissueFraction_MRSDP_selected_lh)); 
    bb = TissueFraction_MRSDP(lh_WMmax_voxel,1); 
    bb1 = TissueFraction_MRSDP(lh_WMmax_voxel,2); 
end

% T(kk,:) = table({subject_name},nVoxels_final,mean(WM_vm),mean(TNAA_absolute_GG),mean(TCho_absolute_GG),mean(TCr_absolute_GG),...
%     'VariableNames', {'Subject_name','Number of voxels', 'WM volume fraction','TNAA absolute conc (mM)', 'TCho absolute conc (mM)', 'TCr absolute conc (mM)'}); 

T(kk,:) = table({subject_name},aa,aa1,bb,bb1,...
    mean(TNAA_absolute_GG([rh_WMmax_voxel_conc_idx lh_WMmax_voxel_conc_idx])),mean(TCho_absolute_GG([rh_WMmax_voxel_conc_idx lh_WMmax_voxel_conc_idx])),mean(TCr_absolute_GG([rh_WMmax_voxel_conc_idx lh_WMmax_voxel_conc_idx])),...
    'VariableNames', {'Subject_name','Voxel_location_rh_x','Voxel_location_rh_y','Voxel_location_lh_x','Voxel_location_lh_y', ...
    'TNAA absolute conc (mM)', 'TCho absolute conc (mM)', 'TCr absolute conc (mM)'}); 

% if nVoxels_final ~= 0
%     T(kk,:) = table({subject_name},nVoxels_final,max(WM_vm),TNAA_absolute_GG(WM_vm == max(WM_vm)),TCho_absolute_GG(WM_vm == max(WM_vm)),TCr_absolute_GG(WM_vm == max(WM_vm)),...
%          'VariableNames', {'Subject_name','Number of voxels', 'WM volume fraction','TNAA absolute conc (mM)', 'TCho absolute conc (mM)', 'TCr absolute conc (mM)'}); 
% else
%     T(kk,:) = table({subject_name},nVoxels_final,nan,nan,nan,nan,...
%          'VariableNames', {'Subject_name','Number of voxels', 'WM volume fraction','TNAA absolute conc (mM)', 'TCho absolute conc (mM)', 'TCr absolute conc (mM)'}); 
% end

disp(append('Finished subject ',subject_name))
end

figure
violinplot([table2array(T(:,6)),table2array(T(:,7)),table2array(T(:,8))],{'TNAA','TCho','TCr'})
ylabel('Concentration (mM)','fontweight','bold')
title('Controls','fontweight','bold')
set(gca,'linewidth',2,'fontname','Times New Romans','fontsize',16,'FontWeight','bold')

T_cchd = T; 

T_cchd = addvars(T_cchd, table2array(T_cchd(:,6))./table2array(T_cchd(:,7)), 'NewVariableNames' , 'TNAA/TCho ratio', 'After', 'TCr absolute conc (mM)'); 
T_cchd = addvars(T_cchd, table2array(T_cchd(:,6))./table2array(T_cchd(:,8)), 'NewVariableNames' , 'TNAA/TCr ratio', 'After', 'TNAA/TCho ratio'); 
T_cchd = addvars(T_cchd, table2array(T_cchd(:,7))./table2array(T_cchd(:,8)), 'NewVariableNames' , 'TCho/TCr ratio', 'After', 'TNAA/TCr ratio'); 

writetable(T_cchd,'E:\CHD_MRS\Outputs\Absolute_Concentration_Results\cchd_results_2voxels.xlsx')

%% PT 
clearvars T
subject_names = dir(fullfile(append(Parent_folder,'\Outputs\Tarquin_Concentration_Results\'),'avpt_*')); 

for kk = 1:length(subject_names)
% clearvars -except T kk subject_names Parent_folder

TNAA_absolute_GG = []; 
TCho_absolute_GG = []; 
TCr_absolute_GG = []; 
WM_vm = [];

subject_name = subject_names(kk).name; 
subject_folder = append(Parent_folder,'\dcm\',subject_name(1:4),'\',subject_name); 
MRSDP_output_voxmap_dir = append(Parent_folder,'\Outputs\Scripts_Output\',subject_name); 
dicom_dir = append(subject_folder,'\DICOM'); 

if ~isfolder(MRSDP_output_voxmap_dir)
    error('No voxel maps generated from MRS-DP.py - Run the script first!')
end

if isempty(fullfile(dicom_dir,'*IM_*'))
    error('No voxel maps generated from MRS-DP.py - Run the script first!')
else
    dicomFiles = dir(fullfile(dicom_dir,'*IM_*'));
end

%**************************************************************************
% Check if the output folder is existed or not 
%**************************************************************************
MRSDP_output_voxmap_dir_new = append(Parent_folder,'\Outputs');
if ~isfolder(MRSDP_output_voxmap_dir_new)
    mkdir(MRSDP_output_voxmap_dir_new)
end

if ~isfolder(append(MRSDP_output_voxmap_dir_new,'\Scripts_Output'))
    mkdir(append(MRSDP_output_voxmap_dir_new,'\Scripts_Output'))
end
if ~isfolder(append(MRSDP_output_voxmap_dir_new,'\Scripts_Output\',subject_name))
    mkdir(append(MRSDP_output_voxmap_dir_new,'\Scripts_Output\',subject_name))
end

if ~isfolder(append(MRSDP_output_voxmap_dir_new,'\Tarquin_Concentration_Results'))
    mkdir(append(MRSDP_output_voxmap_dir_new,'\Tarquin_Concentration_Results'))
end
if ~isfolder(append(MRSDP_output_voxmap_dir_new,'\Tarquin_Concentration_Results\',subject_name))
    mkdir(append(MRSDP_output_voxmap_dir_new,'\Tarquin_Concentration_Results\',subject_name))
end

if ~isfolder(append(MRSDP_output_voxmap_dir_new,'\Absolute_Concentration_Results'))
    mkdir(append(MRSDP_output_voxmap_dir_new,'\Absolute_Concentration_Results'))
end
if ~isfolder(append(MRSDP_output_voxmap_dir_new,'\Absolute_Concentration_Results\',subject_name))
    mkdir(append(MRSDP_output_voxmap_dir_new,'\Absolute_Concentration_Results\',subject_name))
end

%**************************************************************************
% Read inputs (Tissue Fraction and Tarquin metabolite concentrations) 
%**************************************************************************
TissueFraction = fread(fopen(append(MRSDP_output_voxmap_dir_new,'\Scripts_Output\',subject_name,'\TissueTypesForEachVoxel.txt'), 'r')...
    , '*char')';

Metabolites_Tarquin = fread(fopen(append(MRSDP_output_voxmap_dir_new,'\Tarquin_Concentration_Results\',subject_name,'\',subject_name,'.txt'), 'r')...
    , '*char')';

% Tissue Fractions and voxel locations 
pattern = ['Voxels at position\s*\((\d+),(\d+)\)\s*' ...
           'Grey Matter Percentage:\s*([\d.]+)\s*' ...
           'White Matter Percentage:\s*([\d.]+)\s*' ...
           'CSF Percentage:\s*([\d.]+)'];
tokens = regexp(TissueFraction, pattern, 'tokens');

n = length(tokens);
TissueFraction_MRSDP = zeros(n, 5); % columns: Row, Col, GM, WM, CSF
for i = 1:n
    vals = str2double(tokens{i});
    TissueFraction_MRSDP(i, :) = vals;
end
voxel_position_MRSDP = TissueFraction_MRSDP(:,1:2); 
% Find corner coordinates (Flip the row and column so it matches with Tarquin outputs
voxel_corner_MRSDP = [max(voxel_position_MRSDP(:,2)) max(voxel_position_MRSDP(:,1)); ...
    min(voxel_position_MRSDP(:,2)) max(voxel_position_MRSDP(:,1)); ...
    max(voxel_position_MRSDP(:,2)) min(voxel_position_MRSDP(:,1)); ...
    min(voxel_position_MRSDP(:,2)) min(voxel_position_MRSDP(:,1))]; 

% Tarquin 
% Find all voxel coordinate
coordExpr = 'Row\s*:\s*(\d+),\s*Col\s*:\s*(\d+),\s*Slice\s*:\s*(\d+)';
[starts, ends, tokens] = regexp(Metabolites_Tarquin, coordExpr, 'start', 'end', 'tokens');
voxel_position_tarquin = []; 
for zz = 1:size(tokens,2)
    voxel_position_tarquin(zz,:) = str2double(tokens{zz}(1:2)); 
end
% Find corner coordinates 
voxel_corner_tarquin = [min(voxel_position_tarquin(:,1)) min(voxel_position_tarquin(:,2)); ...
    max(voxel_position_tarquin(:,1)) min(voxel_position_tarquin(:,2)); ...
    min(voxel_position_tarquin(:,1)) max(voxel_position_tarquin(:,2)); ...
    max(voxel_position_tarquin(:,1)) max(voxel_position_tarquin(:,2))]; 

% Relationship between MRSDP and Tarquin: 
voxel_position_tarquin_trans = [max(voxel_position_tarquin(:,1))-voxel_position_tarquin(:,1) max(voxel_position_tarquin(:,2))-voxel_position_tarquin(:,2)];  

% Find corresponding voxels and their tissue fractions
TissueFraction_MRSDP_selected = [];
for zz = 1:size(voxel_position_tarquin_trans,1)
    x = find(TissueFraction_MRSDP(:,2) == voxel_position_tarquin_trans(zz,1) & TissueFraction_MRSDP(:,1) == voxel_position_tarquin_trans(zz,2)); 
    TissueFraction_MRSDP_selected(zz,:) = TissueFraction_MRSDP(x,3:end); 
end

%**************************************************************************
% Correction of Tissue volume fractions into mole tissue fraction
%**************************************************************************
% Density of GM, WM and CSF from Gasparovic 2006 and Ernst 1993
TissueFraction_MRSDP_selected_corr = []; 
TissueFraction_MRSDP_selected_corr(:,1) = 0.78*TissueFraction_MRSDP_selected(:,1); 
TissueFraction_MRSDP_selected_corr(:,2) = 0.65*TissueFraction_MRSDP_selected(:,2); 
TissueFraction_MRSDP_selected_corr(:,3) = 0.97*TissueFraction_MRSDP_selected(:,3); 
TissueFraction_final = TissueFraction_MRSDP_selected_corr./sum(TissueFraction_MRSDP_selected_corr,2); 

%**************************************************************************
% We extract the metabolites concentration (mM) from Tarquin (UNPROCESSED!)
%**************************************************************************
nVoxels = length(tokens);
TNAA_data = [];
TCho_data = [];
TCr_data = [];

FWHM_data = [];
SNR_data = [];

for i = 1:nVoxels
    row = str2double(tokens{i}{1});
    col = str2double(tokens{i}{2});
    slice = str2double(tokens{i}{3});
    % Define block boundaries
    startIdx = ends(i);
    if i < nVoxels
        endIdx = starts(i+1) - 1;
    else
        endIdx = length(Metabolites_Tarquin);
    end
    block = Metabolites_Tarquin(startIdx:endIdx);
    % Extract TNAA data using regex
    tnaaExpr = 'TNAA\s+([\d\.e\+\-]+)\s+([\d\.#IOe\+\-]+)\s+([\d\.e\+\-]+)';
    match = regexp(block, tnaaExpr, 'tokens');
    if ~isempty(match)
        conc = str2double(match{1}{1});
        perc_sd = str2double(match{1}{2});
        sd = str2double(match{1}{3});
        TNAA_data = [TNAA_data; conc, perc_sd, sd];
    end
    % Extract TCho data using regex
    tchoExpr = 'TCho\s+([\d\.e\+\-]+)\s+([\d\.#IOe\+\-]+)\s+([\d\.e\+\-]+)';
    match = regexp(block, tchoExpr, 'tokens');
    if ~isempty(match)
        conc = str2double(match{1}{1});
        perc_sd = str2double(match{1}{2});
        sd = str2double(match{1}{3});
        TCho_data = [TCho_data; conc, perc_sd, sd];
    end
    % Extract TCr data using regex
    tcrExpr = 'TCr\s+([\d\.e\+\-]+)\s+([\d\.#IOe\+\-]+)\s+([\d\.e\+\-]+)';
    match = regexp(block, tcrExpr, 'tokens');
    if ~isempty(match)
        conc = str2double(match{1}{1});
        perc_sd = str2double(match{1}{2});
        sd = str2double(match{1}{3});
        TCr_data = [TCr_data; conc, perc_sd, sd];
    end
    % Extract data quality (FWHM)
    FWHMExpr = 'Metab FWHM \(PPM\)\s*:\s*([\d\.]+)';
    match = regexp(block, FWHMExpr, 'tokens');
    if ~isempty(match)
        fwhm = str2double(match{1}{1});
        FWHM_data = [FWHM_data; fwhm];
    end
    % Extract data quality (SNR)
    SNRExpr = 'SNR residual\s*:\s*([\d\.]+)';
    match = regexp(block, SNRExpr, 'tokens');
    if ~isempty(match)
        snr = str2double(match{1}{1});
        SNR_data = [SNR_data; snr];
    end
end

%**************************************************************************
% We calculate the absolute concentration of each metabolites 
%**************************************************************************
% Units of TE, TR, T1 and T2 relaxation: ms 
% Units of water concentration: mmol/L 
Parameters.conc_H20 = 35880;
Parameters.att_H20 = 0.7; 
Parameters.ProtonWater = 55556; 
% Parameters.H20_WM = 36100; Parameters.H20_GM = 43300; Parameters.H20_CSF = 53800; 
Parameters.T2_WM = 79.2; Parameters.T2_GM = 110; Parameters.T2_CSF = 503; 
Parameters.T1_WM = 832; Parameters.T1_GM = 1331; Parameters.T1_CSF = 3817; 
Parameters.TE = 144; Parameters.TR = 2000; 

% T1 and T2 relaxation times were based on this study: Posse et al. 2007
Parameters.T2_TNAA_GM = 247; Parameters.T2_TCho_GM = 222; Parameters.T2_TCr_GM = 162; 
Parameters.T2_TNAA_WM = 301; Parameters.T2_TCho_WM = 222; Parameters.T2_TCr_WM = 178; 
Parameters.T1_TNAA_GM = 1470; Parameters.T1_TCho_GM = 1250; Parameters.T1_TCr_GM = 1339; 
Parameters.T1_TNAA_WM = 1560; Parameters.T1_TCho_WM = 1210; Parameters.T1_TCr_WM = 1400; 

idx_WM = find(TissueFraction_MRSDP_selected(:,2) >= 0.5);
idx_MRS = find(FWHM_data(idx_WM) < 0.1 & SNR_data(idx_WM) > 5);
nVoxels_final = size(idx_WM(idx_MRS),1); 
for i = 1:size(idx_WM(idx_MRS),1)
    ii = idx_WM(idx_MRS(i)); 
    AttH20GM = exp(-Parameters.TE/Parameters.T2_GM) * (1-exp(-Parameters.TR/Parameters.T1_GM)); 
    AttH20WM = exp(-Parameters.TE/Parameters.T2_WM) * (1-exp(-Parameters.TR/Parameters.T1_WM)); 
    AttH20CSF = exp(-Parameters.TE/Parameters.T2_CSF) * (1-exp(-Parameters.TR/Parameters.T1_CSF)); 
    RelaxH20 = (TissueFraction_MRSDP_selected(ii,1)*AttH20GM)+(TissueFraction_MRSDP_selected(ii,2)*AttH20WM)+(TissueFraction_MRSDP_selected(ii,3)*AttH20CSF);
    FracH20 = (0.81*TissueFraction_MRSDP_selected(ii,1))+(0.71*TissueFraction_MRSDP_selected(ii,2))+(TissueFraction_MRSDP_selected(ii,3));

    AttNAAGM=(1-exp(-Parameters.TR/Parameters.T1_TNAA_GM))*exp(-Parameters.TE/Parameters.T2_TNAA_GM);
    AttNAAWM=(1-exp(-Parameters.TR/Parameters.T1_TNAA_WM))*exp(-Parameters.TE/Parameters.T2_TNAA_WM);
    RelaxNAA=((TissueFraction_MRSDP_selected(ii,1)*AttNAAGM)+(TissueFraction_MRSDP_selected(ii,2)*AttNAAWM))/(1-TissueFraction_MRSDP_selected(ii,3));
    FracNAA = 1-TissueFraction_MRSDP_selected(ii,3); 
    TNAA_absolute_GG(i) = TNAA_data(ii,1)/Parameters.att_H20*(FracH20/FracNAA)*(Parameters.ProtonWater/Parameters.conc_H20)*(RelaxH20/RelaxNAA);

    AttChoGM=(1-exp(-Parameters.TR/Parameters.T1_TCho_GM))*exp(-Parameters.TE/Parameters.T2_TCho_GM);
    AttChoWM=(1-exp(-Parameters.TR/Parameters.T1_TCho_WM))*exp(-Parameters.TE/Parameters.T2_TCho_WM);
    RelaxCho=((TissueFraction_MRSDP_selected(ii,1)*AttChoGM)+(TissueFraction_MRSDP_selected(ii,2)*AttChoWM))/(1-TissueFraction_MRSDP_selected(ii,3));
    FracCho = 1-TissueFraction_MRSDP_selected(ii,3); 
    TCho_absolute_GG(i) = TCho_data(ii,1)/Parameters.att_H20*(RelaxH20/RelaxCho)*((FracH20/FracCho)*Parameters.ProtonWater)/Parameters.conc_H20;

    AttCrGM=(1-exp(-Parameters.TR/Parameters.T1_TCr_GM))*exp(-Parameters.TE/Parameters.T2_TCr_GM);
    AttCrWM=(1-exp(-Parameters.TR/Parameters.T1_TCr_WM))*exp(-Parameters.TE/Parameters.T2_TCr_WM);
    RelaxCr=((TissueFraction_MRSDP_selected(ii,1)*AttCrGM)+(TissueFraction_MRSDP_selected(ii,2)*AttCrWM))/(1-TissueFraction_MRSDP_selected(ii,3));
    FracCr= 1-TissueFraction_MRSDP_selected(ii,3); 
    TCr_absolute_GG(i) = TCr_data(ii,1)/Parameters.att_H20*(RelaxH20/RelaxCr)*((FracH20/FracCr)*Parameters.ProtonWater)/Parameters.conc_H20;

    WM_vm(i) = TissueFraction_MRSDP_selected(ii,2); 
end

% Select a voxel with the highest WM probability in each brain hemisphere (2 voxels in total) 
idx_WM_vm = zeros(length(TissueFraction_MRSDP_selected),1);
idx_WM_vm(idx_WM(idx_MRS)) = 1; 

idx_midline = length(unique(TissueFraction_MRSDP(:,1)))/2; 
rh_include = voxel_position_tarquin_trans(:,2)>idx_midline; 
lh_include = voxel_position_tarquin_trans(:,2)<=idx_midline; 
TissueFraction_MRSDP_selected_rh = TissueFraction_MRSDP_selected((rh_include & idx_WM_vm),2); 
TissueFraction_MRSDP_selected_lh = TissueFraction_MRSDP_selected((lh_include & idx_WM_vm),2); 
if isempty(TissueFraction_MRSDP_selected_rh)
    rh_WMmax_voxel = [];
    rh_WMmax_voxel_conc_idx = [];
    aa = nan; 
    aa1 = nan; 
else
    rh_WMmax_voxel = find((TissueFraction_MRSDP(:,4) == max(TissueFraction_MRSDP_selected_rh)) == 1);
    rh_WMmax_voxel_conc_idx = find(WM_vm == max(TissueFraction_MRSDP_selected_rh)); 
    aa = TissueFraction_MRSDP(rh_WMmax_voxel,1); 
    aa1 = TissueFraction_MRSDP(rh_WMmax_voxel,2); 
end

if isempty(TissueFraction_MRSDP_selected_lh)
    lh_WMmax_voxel = [];
    lh_WMmax_voxel_conc_idx = [];
    bb = nan; 
    bb1 = nan; 
else
    lh_WMmax_voxel = find((TissueFraction_MRSDP(:,4) == max(TissueFraction_MRSDP_selected_lh)) == 1);
    lh_WMmax_voxel_conc_idx = find(WM_vm == max(TissueFraction_MRSDP_selected_lh)); 
    bb = TissueFraction_MRSDP(lh_WMmax_voxel,1); 
    bb1 = TissueFraction_MRSDP(lh_WMmax_voxel,2); 
end

% T(kk,:) = table({subject_name},nVoxels_final,mean(WM_vm),mean(TNAA_absolute_GG),mean(TCho_absolute_GG),mean(TCr_absolute_GG),...
%     'VariableNames', {'Subject_name','Number of voxels', 'WM volume fraction','TNAA absolute conc (mM)', 'TCho absolute conc (mM)', 'TCr absolute conc (mM)'}); 

T(kk,:) = table({subject_name},aa,aa1,bb,bb1,...
    mean(TNAA_absolute_GG([rh_WMmax_voxel_conc_idx lh_WMmax_voxel_conc_idx])),mean(TCho_absolute_GG([rh_WMmax_voxel_conc_idx lh_WMmax_voxel_conc_idx])),mean(TCr_absolute_GG([rh_WMmax_voxel_conc_idx lh_WMmax_voxel_conc_idx])),...
    'VariableNames', {'Subject_name','Voxel_location_rh_x','Voxel_location_rh_y','Voxel_location_lh_x','Voxel_location_lh_y', ...
    'TNAA absolute conc (mM)', 'TCho absolute conc (mM)', 'TCr absolute conc (mM)'}); 

% if nVoxels_final ~= 0
%     T(kk,:) = table({subject_name},nVoxels_final,max(WM_vm),TNAA_absolute_GG(WM_vm == max(WM_vm)),TCho_absolute_GG(WM_vm == max(WM_vm)),TCr_absolute_GG(WM_vm == max(WM_vm)),...
%          'VariableNames', {'Subject_name','Number of voxels', 'WM volume fraction','TNAA absolute conc (mM)', 'TCho absolute conc (mM)', 'TCr absolute conc (mM)'}); 
% else
%     T(kk,:) = table({subject_name},nVoxels_final,nan,nan,nan,nan,...
%          'VariableNames', {'Subject_name','Number of voxels', 'WM volume fraction','TNAA absolute conc (mM)', 'TCho absolute conc (mM)', 'TCr absolute conc (mM)'}); 
% end

disp(append('Finished subject ',subject_name))
end

figure
violinplot([table2array(T(:,6)),table2array(T(:,7)),table2array(T(:,8))],{'TNAA','TCho','TCr'})
ylabel('Concentration (mM)','fontweight','bold')
title('PT','fontweight','bold')
set(gca,'linewidth',2,'fontname','Times New Romans','fontsize',16,'FontWeight','bold')

T_avpt = T; 
T_avpt = addvars(T_avpt, table2array(T_avpt(:,6))./table2array(T_avpt(:,7)), 'NewVariableNames' , 'TNAA/TCho ratio', 'After', 'TCr absolute conc (mM)'); 
T_avpt = addvars(T_avpt, table2array(T_avpt(:,6))./table2array(T_avpt(:,8)), 'NewVariableNames' , 'TNAA/TCr ratio', 'After', 'TNAA/TCho ratio'); 
T_avpt = addvars(T_avpt, table2array(T_avpt(:,7))./table2array(T_avpt(:,8)), 'NewVariableNames' , 'TCho/TCr ratio', 'After', 'TNAA/TCr ratio'); 

writetable(T_avpt,'E:\CHD_MRS\Outputs\Absolute_Concentration_Results\avpt_results_2voxels.xlsx')



%% Overview stats
idx_groups = [ones(length(table2array(T_cchd(:,4))), 1); ones(length(table2array(T_achd(:,4))), 1)*2; ones(length(table2array(T_avpt(:,4))), 1)*3];

figure
violinplot([table2array(T_cchd(:,6));table2array(T_achd(:,6));table2array(T_avpt(:,6))],idx_groups)
ylabel('Concentration (mM)','fontweight','bold'); xticklabels({'CTL','CHD','PT'})
title('tNAA','fontweight','bold')
set(gca,'linewidth',2,'fontname','Times New Romans','fontsize',16,'FontWeight','bold')
ylim([0 30]); box off 
figure
violinplot([table2array(T_cchd(:,7));table2array(T_achd(:,7));table2array(T_avpt(:,7))],idx_groups)
ylabel('Concentration (mM)','fontweight','bold'); xticklabels({'CTL','CHD','PT'})
title('tCho','fontweight','bold')
set(gca,'linewidth',2,'fontname','Times New Romans','fontsize',16,'FontWeight','bold')
ylim([0 5]); box off 
figure
violinplot([table2array(T_cchd(:,8));table2array(T_achd(:,8));table2array(T_avpt(:,8))],idx_groups)
ylabel('Concentration (mM)','fontweight','bold'); xticklabels({'CTL','CHD','PT'})
title('tCr','fontweight','bold')
set(gca,'linewidth',2,'fontname','Times New Romans','fontsize',16,'FontWeight','bold')
ylim([4 16]); box off 

% Ratio plot
figure
violinplot([table2array(T_cchd(:,9));table2array(T_achd(:,9));table2array(T_avpt(:,9))],idx_groups)
ylabel('TNAA/TCho ratio','fontweight','bold'); xticklabels({'CTL','CHD','PT'})
set(gca,'linewidth',2,'fontname','Times New Romans','fontsize',16,'FontWeight','bold')
ylim([2 10]); box off 
figure
violinplot([table2array(T_cchd(:,10));table2array(T_achd(:,10));table2array(T_avpt(:,10))],idx_groups)
ylabel('TNAA/TCr ratio','fontweight','bold'); xticklabels({'CTL','CHD','PT'})
set(gca,'linewidth',2,'fontname','Times New Romans','fontsize',16,'FontWeight','bold')
ylim([0.5 2.5]); box off 
figure
violinplot([table2array(T_cchd(:,11));table2array(T_achd(:,11));table2array(T_avpt(:,11))],idx_groups)
ylabel('TCho/TCr ratio','fontweight','bold'); xticklabels({'CTL','CHD','PT'})
set(gca,'linewidth',2,'fontname','Times New Romans','fontsize',16,'FontWeight','bold')
ylim([0.2 0.45]); box off 

% Mean and std of concentrations 
for ii = 6:11
    mean_meta(ii-2,1) = mean(table2array(T_cchd(:,ii)),'omitnan');
    mean_meta(ii-2,2) = mean(table2array(T_achd(:,ii)),'omitnan');
    mean_meta(ii-2,3) = mean(table2array(T_avpt(:,ii)),'omitnan');
    std_meta(ii-2,1) = std(table2array(T_cchd(:,ii)),'omitnan');
    std_meta(ii-2,2) = std(table2array(T_achd(:,ii)),'omitnan');
    std_meta(ii-2,3) = std(table2array(T_avpt(:,ii)),'omitnan');
end



