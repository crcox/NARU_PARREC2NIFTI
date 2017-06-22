function par2nii = spm_cfg_par2nii
% SPM Configuration file for Coregister
%__________________________________________________________________________
% Copyright (C) 2005-2014 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_par2nii.m 6148 2014-09-03 15:49:04Z guillaume $
if ~isdeployed, addpath(fullfile(spm('dir'),'toolbox','PARREC2NIFTI')); end

%--------------------------------------------------------------------------
% par PAR Headers
%--------------------------------------------------------------------------
par         = cfg_files;
par.tag     = 'parfiles';
par.name    = 'PAR header files';
par.help    = {'These are the Phillips PAR images that are to be converted to NIFTI-1 format (techinically, "n+1", meaning that a single .nii file with header and data combined will be produced).'};
par.filter  = 'any';
par.ufilter = {'.PAR','.par'};
par.num     = [1 inf];

%--------------------------------------------------------------------------
% epi_branch Single Echo EPI Scan Protocol
%--------------------------------------------------------------------------
epi_protocol_label = cfg_const();
epi_protocol_label.tag = 'label';
epi_protocol_label.name = 'Protocol Label';
epi_protocol_label.help = {'A constant which uniquely identifies this scan protocol.'};
epi_protocol_label.val  = {'EPI'};

epi_branch = cfg_branch();
epi_branch.tag = 'epi';
epi_branch.name = 'Single Echo EPI';
epi_branch.help = {'Single Echo EPI'};
epi_branch.val = {epi_protocol_label};
%--------------------------------------------------------------------------
% epide_branch Dual Echo EPI Scan Protocol
%--------------------------------------------------------------------------
epide_protocol_label = cfg_const();
epide_protocol_label.tag = 'label';
epide_protocol_label.name = 'Protocol Label';
epide_protocol_label.help = {'A constant while uniquely identifies this scan protocol.'};
epide_protocol_label.val  = {'EPI_DE'};

epide_echo_labels = cfg_entry();
epide_echo_labels.tag    = 'echo_labels';
epide_echo_labels.name   = 'Echo Labels';
epide_echo_labels.help   = {
    'Provide a list of labels for the echoes within the files as a cell array. For example: {''short'',''long''}.'
    ''
    'These will be appended to the NIFTI filenames. You must enter a label for every echo in the file---no more, no less.'
    ''
    'Providing the wrong number of labels will result in an error. The n-th echo label you provide will be applied to the n-th longest echo contained within the PAR/REC.'
    ''
    'Order matters! Make sure to list labels in order of ascending length. In other words, short comes before long.'
    ''};
epide_echo_labels.val    = {{'short', 'long'}};

epide_branch = cfg_branch();
epide_branch.tag = 'epi_de';
epide_branch.name = 'Dual Echo EPI';
epide_branch.help = {'Dual Echo EPI'};
epide_branch.val = {epide_protocol_label, epide_echo_labels};
%--------------------------------------------------------------------------
% b0_branch B0 Scan Protocol
%--------------------------------------------------------------------------
b0_protocol_label = cfg_const();
b0_protocol_label.tag = 'label';
b0_protocol_label.name = 'Protocol Label';
b0_protocol_label.help = {'A constant while uniquely identifies this scan protocol.'};
b0_protocol_label.val  = {'B0'};

b0_branch = cfg_branch();
b0_branch.tag = 'b0';
b0_branch.name = 'B0';
b0_branch.help = {'B0 (phase/magnitude map)'};
b0_branch.val = {b0_protocol_label};
%--------------------------------------------------------------------------
% t1_branch T1 Scan Protocol
%--------------------------------------------------------------------------
t1_protocol_label = cfg_const();
t1_protocol_label.tag = 'label';
t1_protocol_label.name = 'Protocol Label';
t1_protocol_label.help = {'A constant while uniquely identifies this scan protocol.'};
t1_protocol_label.val  = {'T1'};

t1_branch = cfg_branch();
t1_branch.tag = 't1';
t1_branch.name = 'T1';
t1_branch.help = {'T1 anatomical'};
t1_branch.val = {t1_protocol_label};
%--------------------------------------------------------------------------
% stype Scan Protocol
%--------------------------------------------------------------------------
stype         = cfg_choice;
stype.tag     = 'stype';
stype.name    = 'Scan Protocol';
stype.help    = {'Scans are collected under a variety of protocols which require different conversion routines. Task and resting state EPI scans are handled in the same way, and are saved as 4D NIFTI files. T1 and B0 files will be converted into 3D NIFTI files.'};
stype.values  = {epi_branch, epide_branch, t1_branch, b0_branch};
stype.val     = {epi_branch};

%%--------------------------------------------------------------------------
%% stype Scan Type
%%--------------------------------------------------------------------------
%stype         = cfg_menu;
%stype.tag     = 'stype';
%stype.name    = 'Scan Protocol';
%stype.help    = {'Scans are collected under a variety of protocols which require different conversion routines. Task and resting state EPI scans are handled in the same way, and are saved as 4D NIFTI files. T1 and B0 files will be converted into 3D NIFTI files.'};
%stype.labels  = {
%                'Single Echo EPI'
%                'Dual Echo EPI'
%                'T1 Structural'
%                'B0 Magnitude/Phase Map'
%}';
%stype.values  = {'EPI','EPI_DE', 'T1', 'B0'};
%stype.val     = {0};

%--------------------------------------------------------------------------
% dataset Data set
%--------------------------------------------------------------------------
dataset      = cfg_branch;
dataset.tag  = 'dataset';
dataset.name = 'Data set';
dataset.help = {'Collections of PAR files to convert, grouped by protocol.'};
dataset.val  = {par stype};

%--------------------------------------------------------------------------
% dtc Data to convert
%--------------------------------------------------------------------------
dtc        = cfg_repeat;
dtc.tag    = 'dtc';
dtc.name   = 'Data to convert';
dtc.help   = {'Collections of PAR files to convert, grouped by protocol.'};
dtc.values = {dataset};
dtc.num    = [1 Inf];

%--------------------------------------------------------------------------
% outdir NIFTI output directory
%--------------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'NIFTI output directory';
outdir.help    = {'Specify the directory where the .nii file should be written.'};
outdir.filter  = 'dir';
outdir.ufilter = {'.*'};
outdir.num     = [1 1];

%--------------------------------------------------------------------------
% par2nii Convert PAR to NIFTI
%--------------------------------------------------------------------------
par2nii       = cfg_exbranch;
par2nii.tag   = 'par2nii';
par2nii.name  = 'Convert PAR to NIFTI';
par2nii.val   = {dtc outdir};
par2nii.help  = {'This is for converting Phillips PAR/REC files to NIFTI.'};
par2nii.prog  = @spm_run_par2nii;
par2nii.vout  = @vout;

%==========================================================================
function dep = vout(job)
dep(1)            = cfg_dep;
dep(1).sname      = 'Converted Images';
dep(1).src_output = substruct('.','cfiles');
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
