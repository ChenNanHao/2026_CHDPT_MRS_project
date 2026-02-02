%**************************************************************************
% We generate the DICOM outputs for the reference of voxel selection in the
% Tarquin software
% 
% Things need to be done beforehand: 
%       1. DICOM T1w images should be in Philips Classic format (Tarquin
%       readable) 
% 
%       2. MRS-DP.py should be run first to obtain tissue segmentation and
%       the estimated coverage of good MRS-acquired voxel (acqusition box
%       overlapped with Volume of interest on the scanner) 
%
%**************************************************************************
clear all
Parent_folder = 'E:\CHD_MRS'; % Hard-coded now -- need to be changed 

%subject_names = dir(fullfile(append(Parent_folder,'\Outputs\Scripts_Output'),'achd_*')); 
subject_names = dir(fullfile(append(Parent_folder,'\Outputs\Tarquin_whole_fitting\'),'avpt_*')); 

for ii = 1:length(subject_names)

subject_name = subject_names(ii).name; 
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
    mkdir(append(MRSDP_output_voxmap_dir_new,'\Script_Output'))
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
% Take the necessary files and remove those redundant files
%**************************************************************************
voxmap = niftiread(append(MRSDP_output_voxmap_dir,'\WhiteMatterVoxel\constructedWhiteMatterVoxel.nii.gz')); 

% Rotate to the same directionality as DICOM files 
for z = 1:size(voxmap,3)
    voxmap(:,:,z) = imrotate(voxmap(:,:,z),90); 
end

% Write total voxel coverage in the DICOM format 
voxmask = voxmap; 
mkdir(append(subject_folder,'/DICOM_voxel_mask'));
for z = 1:size(voxmap,3)
    % Read original DICOM header for slice z
    origHeader = dicominfo(fullfile(dicom_dir, dicomFiles(z).name));

    % Replace image with masked data but keep header
    newImage = int16(voxmask(:,:,z)*100);  % Ensure data type compatibility

    % Create unique UID for output
    origHeader.SOPInstanceUID = dicomuid;

    % Write the DICOM slice
    outputFile = fullfile(append(subject_folder,'/DICOM_voxel_mask'), sprintf('masked_slice_%03d', z));
    dicomwrite(newImage, outputFile, origHeader, 'CreateMode', 'Copy');
end

% Write qualified white matter voxels in the DICOM format 
voxmap(voxmap < 0.7) = 0; 
mkdir(append(subject_folder,'/DICOM_WMvoxels_70WM'));
for z = 1:size(voxmap,3)
    % Read original DICOM header for slice z
    origHeader = dicominfo(fullfile(dicom_dir, dicomFiles(z).name));

    % Replace image with masked data but keep header
    newImage = int16(voxmap(:,:,z)*100);  % Ensure data type compatibility

    % Create unique UID for output
    origHeader.SOPInstanceUID = dicomuid;

    % Write the DICOM slice
    outputFile = fullfile(append(subject_folder,'/DICOM_WMvoxels_70WM'), sprintf('masked_slice_%03d', z));
    dicomwrite(newImage, outputFile, origHeader, 'CreateMode', 'Copy');
end

disp(append('Finished subject ',subject_name))
end
