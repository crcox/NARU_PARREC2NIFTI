function [dyn, image_meta, image_info] = convert_par_to_nii( par_fullpaths, outdir, scan_type, echo_labels )
%CONVERT_PAR_TO_NII SPM interface to conversion functions
%   Four broad classes of fMRI PAR/REC files are supported:
%    EPI, Dual Echo EPI, B0 (for phase mapping), and T1.
%
%   All file types are read the same way, and contain the same range of
%   information. That is to say, the PAR/REC format is consistent across these
%   file types. The distinction between types has to do with what information
%   will be checked for/expected within the PAR header, and how to split a
%   single PAR/REC into multiple NIFTI outputs. For example, a B0 protocol will
%   have involved both phase and magnitude scans, and these will need to be
%   written to separate NIFTI files. Similarly, protocols involving multiple
%   echoes will need to be parsed into different files.
%
%   In truth, this function is little more than a wrapper around load_parrec.
%   The convert_ functions are all quite similar. Differences come down to
%   expectations (a SingleEcho conversion will throw an error if applied to a
%   file with multiple echoes; a DualEcho conversion will throw an error if the
%   number of provided labels does not match the number of echoes in the file).
%   At some point, it may be a good idea to refactor the code to cut down on
%   redundancy.
%
%   It is worth noting that the the DualEcho function will work for any number
%   of echoes, even for a single echo. The only condition is that echo labels
%   must be provided, and the number of labels must match the number of echoes
%   referenced within the PAR header.
%
%   Dependencies:
%       SPM 8+
    if iscell(outdir) && numel(outdir) == 1
        outdir = outdir{1};
    end

    switch scan_type
        case 'EPI' % Including task and resting state, gradient or spin echo (I think!)
            [dyn, image_meta, image_info] = convert_SingleEcho_4D(par_fullpaths, outdir);
        case 'EPI_DE' % Including task and resting state, gradient or spin echo (I think!)
            [dyn, image_meta, image_info] = convert_DualEcho_4D(par_fullpaths, echo_labels, outdir);
        case 'B0'
            [dyn,image_meta, image_info]  = convert_PhaseMagnitude_3D(par_fullpaths, outdir);
        case 'T1'
            [dyn, image_meta, image_info] = convert_T1_3D(par_fullpaths, outdir);
        case 'survey'
            % This converstion is not implemented yet. Data are
            % acquired from all 3 orientations.
            % Although, this is not a useful scan for any research purpose,
            % so it probably will never need to be converted.
        case 'MPR_RPS'
            % This converstion is not implemented yet. Data are
            % acquired from 2 orientations.
    end
end

function [dyn, image_meta, image_info] = convert_SingleEcho_4D(par_fullpaths, outdir)
    par_basenames = cell(1, numel(par_fullpaths));
    for i = 1:numel(par_fullpaths)
        [~, par_basenames{i}] = fileparts(par_fullpaths{i});
    end

    if ~exist(outdir, 'dir');
        mkdir(outdir);
    end
    dyn = cell(1, numel(par_fullpaths));
    image_meta = cell(1, numel(par_fullpaths));
    image_info = cell(1,numel(par_fullpaths));
    for i=1:numel(par_fullpaths);
        [dyn{i}, image_meta{i}, image_info{i}] = load_parrec(par_fullpaths{i}, ...
            'dataformat_source', 'int16', ...
            'dataformat_target', 'int16', ...
            'fliplr', true);
        ntr = image_meta{i}.Max_number_of_dynamics;
        % Check that the number of echoes in the PAR header matches expectations
        if numel(dyn{i}) > 1
            error('Reference to %d echoes were identified in the PAR file, but expected a single echo. Aborting ...', numel(dyn{i}));
        end
        % This is done so that labels (which should be ordered from shortest to
        % longest) match with the proper echo.
        [~,ix_sort] = sort([dyn{i}.echo_time]);
        dyn{i} = dyn{i}(ix_sort);

        RI = dyn{i}.rescale_intercept;
        RS = dyn{i}.rescale_slope;
        SS = dyn{i}.scale_slope;
        slope = 1/SS;
        intercept = RI/(RS*SS);
        r = scaleoffset_nii(dyn{i}, image_meta{i});
        dim = size(dyn{i}.data);
        dyn{i}.run = i;
        dyn{i}.data = (dyn{i}.data .* slope) + intercept;
        dyn{i}.hdr = struct(                   ...
            'fname',   '',                     ...
            'dim',     dim(1:3),               ...
            'pinfo',   [slope,intercept,352]', ...
            'dt',      [spm_type('int16'),1], ...
            'mat',     r.mat,                  ...
            'n',       [],                     ...
            'descrip', '',                     ...
            'private', []);
        % Write out one volume per TR
        tmp_root = tempname;
        mkdir(tmp_root);
        for t = 1:size(dyn{i}.data, 4)
            basename = sprintf('%s_%03d.nii', par_basenames{i}, t);
            cur_nii_vol = fullfile(tmp_root, basename);
            dyn{i}.hdr.fname = cur_nii_vol;
            spm_write_vol(dyn{i}.hdr, squeeze(dyn{i}.data(:,:,:,t)));
        end
        % Compose a list off all the file names (probably should just
        % populate this as they are made... but this is quick enough.)
        fmt_tr = strjoin({par_basenames{i}, '%03d.nii,1'}, '_');
        basenames_tr = arrayfun(@(i) sprintf(fmt_tr,i), 1:ntr, 'Unif', 0);
        paths_tr = cellfun(@(b) fullfile(tmp_root, b), basenames_tr, 'Unif', 0);
        basename_4D_nii = sprintf('%s.nii', par_basenames{i});
        basename_4D_mat = sprintf('%s_parhdr.mat', par_basenames{i});
        cur_nii_vol_4D = fullfile(outdir, basename_4D_nii);
        cur_nii_mat_4D = fullfile(outdir, basename_4D_mat);
        spm_file_merge(paths_tr(:), cur_nii_vol_4D, spm_type('int16'));
        dyn{i}.data = [];
        naru_write_par_auxilary(cur_nii_mat_4D, dyn{i}, image_meta{i}, image_info{i});
        if exist(cur_nii_vol_4D,'file')
            rmdir(tmp_root,'s')
        end
    end
end
function [dyn, image_meta, image_info] = convert_DualEcho_4D(par_fullpaths, echo_labels, outdir)
    par_basenames = cell(1, numel(par_fullpaths));
    for i = 1:numel(par_fullpaths)
        [~, par_basenames{i}] = fileparts(par_fullpaths{i});
    end

    if ~exist(outdir, 'dir');
        mkdir(outdir);
    end
    dyn = cell(1, numel(par_fullpaths));
    image_meta = cell(1, numel(par_fullpaths));
    image_info = cell(1,numel(par_fullpaths));
    for i=1:numel(par_fullpaths);
        [dyn{i}, image_meta{i}, image_info{i}] = load_parrec(par_fullpaths{i}, ...
            'dataformat_source', 'int16', ...
            'dataformat_target', 'int16', ...
            'fliplr', true);
        ntr = image_meta{i}.Max_number_of_dynamics;
        % Check that the number of echoes in the PAR header matches expectations
        if numel(echo_labels) ~= numel(dyn{i})
            error('Reference to %d echoes were identified in the PAR file, but %d labels were provided. Aborting ...', numel(dyn{i}), numel(echo_labels));
        end
        % This is done so that labels (which should be ordered from shortest to
        % longest) match with the proper echo.
        [~,ix_sort] = sort([dyn{i}.echo_time]);
        dyn{i} = dyn{i}(ix_sort);

        for j = 1:numel(dyn{i});
            RI = dyn{i}(j).rescale_intercept;
            RS = dyn{i}(j).rescale_slope;
            SS = dyn{i}(j).scale_slope;
            slope = 1/SS;
            intercept = RI/(RS*SS);
            r = scaleoffset_nii(dyn{i}(j), image_meta{i});
            dim = size(dyn{i}(j).data);
            dyn{i}(j).run = i;
            dyn{i}(j).hdr = struct(                   ...
                'fname',   '',                     ...
                'dim',     dim(1:3),               ...
                'pinfo',   [slope/1000,intercept,352]', ...
                'dt',      [spm_type('int16'),1], ...
                'mat',     r.mat,                  ...
                'n',       [],                     ...
                'descrip', '',                     ...
                'private', []);
            % Write out one volume per TR
            tmp_root = tempname;
            mkdir(tmp_root);
            for t = 1:size(dyn{i}(j).data, 4)
                basename = sprintf('%s_%s_%03d.nii', par_basenames{i}, echo_labels{j}, t);
                cur_nii_vol = fullfile(tmp_root, basename);
                dyn{i}(j).hdr.fname = cur_nii_vol;
                spm_write_vol(dyn{i}(j).hdr, squeeze(dyn{i}(j).data(:,:,:,t)));
            end
            % Compose a list off all the file names
            fmt_tr = strjoin({par_basenames{i}, echo_labels{j}, '%03d.nii,1'}, '_');
            basenames_tr = arrayfun(@(i) sprintf(fmt_tr,i), 1:ntr, 'Unif', 0);
            paths_tr = cellfun(@(b) fullfile(tmp_root, b), basenames_tr, 'Unif', 0);
            basename_4D_nii = sprintf('%s_%s.nii', par_basenames{i}, echo_labels{j});
            basename_4D_mat = sprintf('%s_%s_parhdr.mat', par_basenames{i}, echo_labels{j});
            cur_nii_vol_4D = fullfile(outdir, basename_4D_nii);
            cur_nii_mat_4D = fullfile(outdir, basename_4D_mat);
            spm_file_merge(paths_tr(:), cur_nii_vol_4D, spm_type('int16'));
            dyn{i}(j).data = [];
            naru_write_par_auxilary(cur_nii_mat_4D, dyn{i}(j), image_meta{i}, image_info{i}(j));
            if exist(cur_nii_vol_4D,'file')
                rmdir(tmp_root,'s')
            end
        end
    end
end
function [dyn, image_meta, image_info] = convert_PhaseMagnitude_3D(par_fullpaths, outdir)
    par_basenames = cell(1, numel(par_fullpaths));
    for i = 1:numel(par_fullpaths)
        [~, par_basenames{i}] = fileparts(par_fullpaths{i});
    end

    if ~exist(outdir, 'dir');
        mkdir(outdir);
    end
    dyn = cell(1, numel(par_fullpaths));
    image_meta = cell(1, numel(par_fullpaths));
    image_info = cell(1,numel(par_fullpaths));
    for i=1:numel(par_fullpaths);
        [dyn{i}, image_meta{i}, image_info{i}] = load_parrec(par_fullpaths{i}, ...
            'dataformat_source', 'int16', ...
            'dataformat_target', 'int16', ...
            'fliplr', true);

        for j = 1:numel(dyn{i});
            RI = dyn{i}(j).rescale_intercept;
            RS = dyn{i}(j).rescale_slope;
            SS = dyn{i}(j).scale_slope;
            slope = 1/SS;
            intercept = RI/(RS*SS);
            r = scaleoffset_nii(dyn{i}(j), image_meta{i});
            dim = size(dyn{i}(j).data);
            dyn{i}(j).run = i;
            dyn{i}(j).hdr = struct(                   ...
                'fname',   '',                     ...
                'dim',     dim(1:3),               ...
                'pinfo',   [slope,intercept,352]', ...
                'dt',      [spm_type('int16'),1], ...
                'mat',     r.mat,                  ...
                'n',       [],                     ...
                'descrip', '',                     ...
                'private', []);

            switch dyn{i}(j).image_type_mr
                case 0
                    basename_nii = sprintf('%s_%s.nii', par_basenames{i}, 'magnitude');
                    basename_mat = sprintf('%s_%s_parhdr.mat', par_basenames{i}, 'magnitude');
                case 3
                    basename_nii = sprintf('%s_%s.nii', par_basenames{i}, 'phase');
                    basename_mat = sprintf('%s_%s_parhdr.mat', par_basenames{i}, 'phase');
                otherwise
                    error('Programming error---unexpected image_type_mr code.');
            end
            cur_nii_vol = fullfile(outdir, basename_nii);
            cur_nii_mat = fullfile(outdir, basename_mat);
            dyn{i}(j).hdr.fname = cur_nii_vol;
            spm_write_vol(dyn{i}(j).hdr, squeeze(dyn{i}(j).data));
            dyn{i}(j).data = [];
            naru_write_par_auxilary(cur_nii_mat, dyn{i}(j), image_meta{i}, image_info{i}(j));
        end
    end
end
function [dyn, image_meta, image_info] = convert_T1_3D(par_fullpaths, outdir)
    par_basenames = cell(1, numel(par_fullpaths));
    for i = 1:numel(par_fullpaths)
        [~, par_basenames{i}] = fileparts(par_fullpaths{i});
    end

    if ~exist(outdir, 'dir');
        mkdir(outdir);
    end
    dyn = cell(1, numel(par_fullpaths));
    image_meta = cell(1, numel(par_fullpaths));
    image_info = cell(1,numel(par_fullpaths));
    for i=1:numel(par_fullpaths);
        [dyn{i}, image_meta{i}, image_info{i}] = load_parrec(par_fullpaths{i}, ...
            'dataformat_source', 'int16', ...
            'dataformat_target', 'int16', ...
            'fliplr', true);
        for j = 1:numel(dyn{i});
            RI = dyn{i}(j).rescale_intercept;
            RS = dyn{i}(j).rescale_slope;
            SS = dyn{i}(j).scale_slope;
            slope = 1/SS;
            intercept = RI/(RS*SS);
            r = scaleoffset_nii(dyn{i}(j), image_meta{i});
            dim = size(dyn{i}(j).data);
            dyn{i}(j).run = i;
            dyn{i}(j).hdr = struct(                   ...
                'fname',   '',                     ...
                'dim',     dim(1:3),               ...
                'pinfo',   [slope,intercept,352]', ...
                'dt',      [spm_type('int16'),1], ...
                'mat',     r.mat,                  ...
                'n',       [],                     ...
                'descrip', '',                     ...
                'private', []);
            basename_nii = sprintf('%s.nii', par_basenames{i});
            basename_mat = sprintf('%s_parhdr.mat', par_basenames{i});
            cur_nii_vol = fullfile(outdir, basename_nii);
            cur_nii_mat = fullfile(outdir, basename_mat);
            dyn{i}(j).hdr.fname = cur_nii_vol;
            spm_write_vol(dyn{i}(j).hdr, squeeze(dyn{i}(j).data));
            dyn{i}(j).data = [];
            naru_write_par_auxilary(cur_nii_mat, dyn{i}(j), image_meta{i}, image_info{i}(j));
        end
    end
end

function naru_write_par_auxilary(fname, dyn, image_meta, image_info) %#ok<INUSD>
    % Allows saving files with the intended variable names
    save(fname, 'dyn','image_meta','image_info');
end
%% Detritus ...
% This is a more manual way to do it with the nifti tools ...
% but I don't know how to define the quaternions.
%             nii = nifti.make_nii(dyn{i}(j).data);
%             nii.hdr.dime.pixdim = [1, r.voxel_dim, 0, 0, 0, 0];
%             nii.hdr.dime.scl_slope = slope/1000;
%             nii.hdr.dime.scl_inter = intercept;
%             nii.hdr.dime.vox_offset = 352;
%             nii.hdr.dime.xyzt_units = 10;
%             nii.hdr.hist.descrip = '4D image';
%             nii.hdr.hist.aux_file = '';
%             nii.hdr.hist.qform_code = 2;
%             nii.hdr.hist.sform_code = 2;
%             nii.hdr.hist.qoffset_x = r.mat(1,4);
%             nii.hdr.hist.qoffset_y = r.mat(2,4);
%             nii.hdr.hist.qoffset_z = r.mat(3,4);
%             nii.hdr.hist.srow_x = r.mat(1,:);
%             nii.hdr.hist.srow_y = r.mat(2,:);
%             nii.hdr.hist.srow_z = r.mat(3,:);
%             nii.hdr.hist.magic = 'n+1';
%             basename_4D_nii = sprintf('%s_%s2.nii', par_basenames{i}, echo_labels{j});
%             nifti.save_nii(nii, basename_4D_nii)

% function s = init_image_meta(n)
%     s(n) = struct(...
%         'echo_number', [], ...
%         'image_type_mr', [], ...
%         'scanning_sequence', [], ...
%         'image_pixel_size', [], ...
%         'recon_resolution', [], ...
%         'rescale_intercept', [], ...
%         'rescale_slope', [], ...
%         'scale_slope', [], ...
%         'image_angulation', [], ...
%         'slice_thickness', [], ...
%         'slice_gap', [], ...
%         'image_display_orientation', [], ...
%         'slice_orientation', [], ...
%         'image_type_ed_es', [], ...
%         'pixel_spacing', [], ...
%         'echo_time', [], ...
%         'trigger_time', [], ...
%         'diffusion_b_factor', [], ...
%         'number_of_averages', [], ...
%         'image_flip_angle', [], ...
%         'minimum_RR', [], ...
%         'maximum_RR', [], ...
%         'TURBO_factor', [], ...
%         'Inversion_delay', [], ...
%         'diffusion_b', [], ...
%         'gradient_orientation_number', [], ...
%         'contrast_type', [], ...
%         'diffusion_anisotropy_type', [], ...
%         'label_type', [], ...
%         'data', [], ...
%         'run', [], ...
%         'volume_index', []);
% end
% function s = init_image_info(n)
%     s(n) = struct(...
%         'Patient_name', [], ...
%         'Examination_name', [], ...
%         'Protocol_name', [], ...
%         'Examination_date', [], ...
%         'Series_Type', [], ...
%         'Acquisition_nr', [], ...
%         'Reconstruction_nr', [], ...
%         'Scan_Duration', [], ...
%         'Max_number_of_cardiac_phases', [], ...
%         'Max_number_of_echoes', [], ...
%         'Max_number_of_slices', [], ...
%         'Max_number_of_dynamics', [], ...
%         'Max_number_of_mixes', [], ...
%         'Patient_position', [], ...
%         'Preparation_direction', [], ...
%         'Technique', [], ...
%         'Scan_resolution', [], ...
%         'Scan_mode', [], ...
%         'Repetition_time', [], ...
%         'FOV', [], ...
%         'Water_Fat_shift', [], ...
%         'Angulation_midslice', [], ...
%         'Off_Centre_midslice', [], ...
%         'Flow_compensation', [], ...
%         'Presaturation', [], ...
%         'Phase_encoding_velocity', [], ...
%         'MTC', [], ...
%         'SPIR', [], ...
%         'EPI_factor', [], ...
%         'Dynamic_scan', [], ...
%         'Diffusion', [], ...
%         'Diffusion_echo_time', [], ...
%         'Max_number_of_diffusion_values', [], ...
%         'Max_number_of_gradient_orients', [], ...
%         'Number_of_label_types', []);
% end
