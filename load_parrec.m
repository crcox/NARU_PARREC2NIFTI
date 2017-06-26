function [dynamics, image_meta, image_info] = ...
    load_parrec(filename, version, dataformat_source, bits_source, dataformat_target, bits_target)
% dualecho.LOAD_PARREC Load binary image data from PAR/REC file pair.
% 
%  import(dualecho)
%  [dynamics] = load_parrec(filename) Read dynamic volumes from
%      Phillips PAR/REC into a structured array. The structure will contain
%      information about the acquisition protocol. If there were multiple
%      protocols or echo types, this structure will be multi-dimensional.
%  [dynamics, meta, info] = load_parrec(filename) Will also return
%      metadata from the PAR header file.
%  [dynamics] = load_parrec(filename,version,dataformat_source,bits_source)
%      These options affect how the files are read and interpretted.
%      Currently, the function is only known to recognize data written by
%      the Phillips Research image export tool V4.2. The dataformat and
%      bits options pertain to how data are serialized in the REC file.
%      Because PAR/REC is a closed format, the community can only guess as
%      to how data are encoded in this file. The general consensus is that
%      REC files store uint16 values that are little Endian.
%        - Re: little Endian, see https://github.com/nipy/nibabel/issues/274
%        - Re: uint16, see https://github.com/nipy/nibabel/issues/275
%        (although technically it suggests uint16).
%        - CRC: I have loaded one fMRI PAR/REC file as both uint16 and int16
%        and the resulting matrices are identical.
%        [superceded] CRC: For historical reasons (i.e., this is how it is
%        done in dual_echo_analyse.m), the default assumed format for this
%        function is int16. Because we cannot know for sure, it might be
%        worth loading data as int16 and uint16 (particularly for DTI data,
%        see link above for anecdotal evidence) and proving that they are
%        the same.
%        !!! The previous statement isn't true anymore... I am running with
%        uint16 as the default since that is what Nibable and many others
%        are assuming.
%
%  CRC: I have checked the procedure below against the procedure included
%  as part of dual_echo_analyse.m and they produce identical matrices.
%  CRC: Update... this is not true anymore, because I changed how data are
%  loaded (as double rather than int) and scaled (I scale using the
%  equations specified in the PAR header without KARLS_RESCALEFACTOR. But
%  nothing else was changed, so they still match up to a linear
%  transformation.
%
%  CRC: A note on computational efficiency. This implementation reads the
%  whole PAR file and the whole REC file into memory and holds it there.
%  Searching for text fields in the PAR file involve scanning all the text
%  several times. Clearly, this could be improved upon, but PAR files are
%  not huge so we get away with it. The package ReadData3D (posted to
%  Matlab Central by Dirk-Jan Kroon (c) 2010) includes a far more efficient
%  functoin for reading from the PAR file, if anyone ever feels like
%  refactoring this code to be more like that.
%
% <Chris Cox 23/03/2017>
% <Adapted from work by Ajay Halai>

    % KARLS_RESCALEFACTOR = 1000;
    % CRC: For the time being, I am operating without this as a default.
    DEFAULT_TARGET_DATAFORMAT = 'double';
    % N.B. if left empty, will also default to source type.
    % Since this function loads the data into memory and returns it for
    % further analysis we probably want to load as double. If loaded as
    % double, and the data were stored as integers in the REC file, the
    % returned values will be scaled.
    if iscell(filename) 
        if numel(filename) == 1
            filename = filename{1};
        else
            error('Expected a single filename.');
        end
    end

    [pr_path,pr_name,~] = fileparts(filename);
    par_fullfile = [fullfile(pr_path,pr_name),'.PAR'];
    rec_fullfile = [fullfile(pr_path,pr_name),'.REC'];

    %% Load text data from PAR file.
    fid = fopen(par_fullfile);
    data_description = textscan(fid,'%s','Delimiter','\n');
    data_description = data_description{1};
    fclose(fid);

    % if version number is not provided, look it up.
    if nargin < 2 || isempty(version)
        version = get_RIET_version(data_description);
    end
    
    %% Initialize data structures
    image_meta = init_meta_struct(version);
    image_info = init_info_struct(version);
    
    %% Parse metadata
    image_meta = parse_par_structured_text(image_meta, data_description);
    
    %% Parse image info
    image_info = parse_par_tabular_text(image_info, data_description);
    
    %% Parse data and information by protocol/echo type (and return)
    [dynamics, dyn_ind] = info_by_configuration(image_info, version);
    % dyn_ind gives a selector for each slice. Since we need to select by
    % volume, and protocol information is not varying by slice, we can
    % reshape so that each column corresponds to a volume containing slicen
    % slices, and skim off the first row to index into the volumes.
    
    % extract number of slices and volumes (i.e., dynamics; determined by
    % gradient)
    slicen = image_meta.Max_number_of_slices;
    gradientn = image_meta.Max_number_of_dynamics * numel(dynamics);
    %gradientn = 2 * gradient;
    
    dyn_ind = reshape(dyn_ind, slicen, gradientn);

    % extract the pixel resolution of each slice
    x = image_info(1).recon_resolution(1,1); % (right-left) (check these comments...)
    y = image_info(1).recon_resolution(1,2); % (anterior-posterior)

    % extract voxel dimensions (in mm)
    outdim = zeros(1,3);
    outdim(1:2) = image_info(1).pixel_spacing(1,1:2);
    outdim(3) = image_info(1).slice_thickness(1,1);

    % extract scaling factor, which is applied to the functional data.
    RI = image_info(1).rescale_intercept(1,1);
    RS = image_info(1).rescale_slope(1,1);
    SS = image_info(1).scale_slope(1,1);
    
    scalingfactor = 1/SS;
    scalingintercept = RI/(RS*SS);

    % gradientnhalf = round(gradientn/2);
    num_images = gradientn * slicen;
    num_elements = num_images * x * y;
    if length(image_info(1).slice_number) ~= num_images
        error('load_parrec:badMetaRead', 'The number of slices/images in the REC file implied by the max number of slices per volume and the number of volumes acquired differ from the number of slices listed in the info table. Either metadata is being misread or the source data may be corrupted.');
    end
    
    %% Note echo times for each dynamic
    echo_number = image_info(1).echo_number;
    TE = image_info(1).echo_time;
    tmp = unique(sortrows([echo_number,TE]),'rows');
    % echo_index = tmp(:,1);
    TE_index = tmp(:,2);
    
    %% Compose the data specification
    % If number of bits per value is not specified, read from file.
    % The data is most likely stored as a 16 bit integer. There is a
    % question as to whether it is signed or unsigned. We might want to
    % represent the data as 
    if nargin < 3 || isempty(dataformat_source)
      dataformat_source = 'uint';
    end
    if nargin < 4 || isempty(bits_source)
        bits_source = image_info(1).image_pixel_size(1);
        if ~ismember(bits_source,[8, 16, 32, 64]);
            error('load_parrec:badDataSpec', 'Attempted to read number of bits per value (in REC file) from PAR file, but got a value other than 8, 16, 32, or 64.');
        end
    end
    if nargin < 5 || isempty(dataformat_target)
        if isempty(DEFAULT_TARGET_DATAFORMAT)
            dataformat_target = dataformat_source;
        else
            dataformat_target = DEFAULT_TARGET_DATAFORMAT;
        end
    end
    if nargin < 6 || isempty(bits_target)
        if strcmp(dataformat_target, dataformat_source)
            bits_target = bits_source;
        elseif strcmp(dataformat_target, 'float')
            bits_target = 32;
        elseif strcmp(dataformat_target, 'double')
            bits_target = 64;
        end
    end
    
    switch dataformat_source
    case {'float','double'}
        dataformat_source_code = dataformat_source;
    case {'int','uint'}
        dataformat_source_code = sprintf('%s%d', dataformat_source, bits_source);
    end
    
    switch dataformat_target
    case {'float','double'}
        dataformat_target_code = dataformat_target;
    case {'int','uint'}
        dataformat_target_code = sprintf('%s%d', dataformat_target, bits_target);
    end
    dataformatcode = sprintf('%s=>%s', dataformat_source_code, dataformat_target_code);
    
    %% Report data info before reading
    fprintf('Started %s\n', datestr(now, 'dd mmmm yyyy, HH:MM:SS.FFF'));
    fprintf('PAR/REC export version: %s\n', version);
    fprintf('Data format code: %s\n', dataformatcode);
    fprintf('Recon Resolution (Field Of View): %d x %d (pixels)\n', x, y);
    fprintf('Number of slices: %d\n', slicen);
    fprintf('Max number of dynamics (functional volumes/dynamics): %.2f\n', gradientn);
    fprintf('Voxel size: %.2f x %.2f x %.2f (mm)\n', outdim);
    fprintf('Total number of images (i.e., slices) to read from REC: %d\n', num_images);
    fprintf('Total number of data points (i.e., voxels) to read from REC: %d\n', num_elements);
    fprintf('Echo time(s):'); fprintf(' %.2f', TE_index); fprintf(' (ms)\n');
    
    %% Read data from REC file
    fid = fopen(rec_fullfile); %open data stream for rec file
    img = fread(fid,num_elements,dataformatcode)'; % read stream as specified in dataformatcode
    fclose(fid);
    
    %% Reshape and rescale vector into a recognizable format (3D+time).
    % CRC: When it comes to rescaling, there are two distinct formulas.
    % From the PAR header:
    %  # === PIXEL VALUES =============================================================
    %  #  PV = pixel value in REC file, FP = floating point value, DV = displayed value on console
    %  #  RS = rescale slope,           RI = rescale intercept,    SS = scale slope
    %  #  DV = PV * RS + RI             FP = DV / (RS * SS)
    %
    % The "displayed value" is different than the "floating point value"
    % for some reason. Generally, people don't seem to think it matters
    % much but that the floating point (FP) value to be the most
    % definiative. To move directly from the stored pixel value (PV) to the
    % FP value:
    %   FP = PV * (1/SS) + RI/(RS*SS)
    % So if we were to provide a slope and intercept to the Nifti format,
    % slope = 1/SS and intercept = RI/(RS*SS). If we were going to store
    % the data again to disk as a nifti file we would reverse the
    % operation:
    %   PV = uint16((FP - intercept) / slope);
    if strcmpi(dataformat_target, 'double') && strcmpi(dataformat_source, 'int')
        dual4d = (reshape(img,[x y slicen gradientn]) * scalingfactor) + scalingintercept; % ./ KARLS_RESCALEFACTOR;
    else
        dual4d = (reshape(img,[x y slicen gradientn]));
    end
    
    % loop over the different types of dynamics and plug in the data.
    for i = 1:numel(dynamics)
        dynamics(i).data = dual4d(:,:,:,dyn_ind(1,:)==i);
    end
end

function image_meta = init_meta_struct(version)
    switch lower(version)
    case 'v4.2'
        image_meta = struct( ...
            'Patient_name', [], ...
            'Examination_name', [], ...
            'Protocol_name', [], ...
            'Examination_date', [], ... /time
            'Series_Type', [], ...
            'Acquisition_nr', [], ...
            'Reconstruction_nr', [], ...
            'Scan_Duration', [], ... _[sec]
            'Max_number_of_cardiac_phases', [], ...
            'Max_number_of_echoes', [], ...
            'Max_number_of_slices', [], ... /locations
            'Max_number_of_dynamics', [], ...
            'Max_number_of_mixes', [], ...
            'Patient_position', [], ...
            'Preparation_direction', [], ...
            'Technique', [], ...
            'Scan_resolution', [], ...  (x,y)
            'Scan_mode', [], ...
            'Repetition_time', [], ... _[ms]
            'FOV', [], ... _(ap,fh,rl)_[mm]
            'Water_Fat_shift', [], ... _[pixels]
            'Angulation_midslice', [], ... (ap,fh,rl)[degr]
            'Off_Centre_midslice', [], ... (ap,fh,rl)_[mm]
            'Flow_compensation', [], ...  _<0=no_1=yes>_?
            'Presaturation', [], ...     <0=no_1=yes>_?
            'Phase_encoding_velocity', [], ... _[cm/sec]
            'MTC', [], ...               <0=no_1=yes>_?
            'SPIR', [], ...              <0=no_1=yes>_?
            'EPI_factor', [], ...        <0,1=no_EPI>
            'Dynamic_scan', [], ...      <0=no_1=yes>_?
            'Diffusion', [], ...         <0=no_1=yes>_?
            'Diffusion_echo_time', [], ... _[ms]
            'Max_number_of_diffusion_values', [], ...
            'Max_number_of_gradient_orients', [], ...
            'Number_of_label_types', [] ...   <0=no_ASL>
        );
    otherwise
        error('load_parrec:init_meta_struct:unknownVersion', 'Unknown PAR/REC version.');
    end
end

function image_info = init_info_struct(version)
% Will return a 2-D structured array. The first element is an empty array, the
% second is the number of columns that should be filled into the empty array
% from the par-file image data.
    switch lower(version)
    case 'v4.2'
      image_info = struct(...
          'slice_number', {[], 1}, ...                             (integer) *
          'echo_number', {[], 1}, ...                              (integer) *
          'dynamic_scan_number', {[], 1}, ...                      (integer) *
          'cardiac_phase_number', {[], 1}, ...                     (integer) *
          'image_type_mr', {[], 1}, ...                            (integer) x
          'scanning_sequence', {[], 1}, ...                        (integer) x
          'index_in_REC_file',{[], 1}, ...   (in images)           (integer) *
          'image_pixel_size', {[], 1}, ...   (in bits)             (integer) x
          'scan_percentage', {[], 1}, ...                          (integer) *
          'recon_resolution', {[], 2}, ...   (x,y)                 (2*integer)x
          'rescale_intercept', {[], 1}, ...                        (float) x
          'rescale_slope', {[], 1}, ...                            (float) x
          'scale_slope', {[], 1}, ...                              (float) x
          'window_center', {[], 1}, ...                            (integer) *
          'window_width', {[], 1}, ...                             (integer) *
          'image_angulation', {[], 3}, ...  (ap,fh,rl in degrees)  (3*float) x
          'image_offcentre', {[], 3}, ...   (ap,fh,rl in mm)       (3*float) *
          'slice_thickness', {[], 1}, ...   (in mm)                (float) x
          'slice_gap', {[], 1}, ...         (in mm)                (float) x
          'image_display_orientation', {[], 1}, ...                (integer) x
          'slice_orientation', {[], 1}, ... (TRA/SAG/COR)          (integer) x
          'fmri_status_indication', {[], 1}, ...                   (integer) *
          'image_type_ed_es', {[], 1}, ...  (end_diast/end_syst)   (integer) x
          'pixel_spacing', {[], 2}, ...     (x,y) (in mm)          (2*float) x
          'echo_time', {[], 1}, ...                                (float) x
          'dyn_scan_begin_time', {[], 1}, ...                      (float) *
          'trigger_time', {[], 1}, ...                             (float) x
          'diffusion_b_factor', {[], 1}, ...                       (float) x
          'number_of_averages', {[], 1}, ...                       (integer) x
          'image_flip_angle', {[], 1}, ...  (in degrees)           (float) x
          'cardiac_frequency', {[], 1}, ...   (bpm)                (integer) *
          'minimum_RR', {[], 1}, ...  (interval in ms)             (integer) x
          'maximum_RR', {[], 1}, ...  (interval in ms)             (integer) x
          'TURBO_factor', {[], 1}, ...  <0=no_turbo>               (integer) x
          'Inversion_delay', {[], 1}, ...   (in ms)                (float) x
          'diffusion_b', {[], 1}, ... value_number    (imagekey!)  (integer) x
          'gradient_orientation_number', {[], 1}, ... (imagekey!)  (integer) x
          'contrast_type', {[], 1}, ...                            (string) x
          'diffusion_anisotropy_type', {[], 1}, ...                (string) x
          'diffusion', {[], 3}, ...   (ap, fh, rl)                 (3*float) *
          'label_type', {[], 1} ...   (ASL)           (imagekey!)  (integer) x
    );
    end
end

function targetstruct = parse_par_structured_text(targetstruct, partext)
    selection_beg = rowfind(partext, '= GENERAL INFORMATION =') + 1;
    selection_end = rowfind(partext, '= PIXEL VALUES =') - 1;
    parstructtext = partext(selection_beg:selection_end);
    
    if exist('strsplit', 'file') ~= 2
        strsplit = @strsplit_;
    end
    
    % Clean the data of troublesome characters
    parstructtext = strrep(parstructtext, 'Max.', 'Max');

    fn = fieldnames(targetstruct);
    for i = 1:length(fn)
        row_id = rowfind(parstructtext, strrep(fn{i},'_',' '));
        row_str = parstructtext{row_id};
        delim_pos = strfind(row_str, ':');
        row_val = strtrim(row_str((delim_pos(1)+1):end));
        row_vals = strsplit(row_val);
        row_nums = str2double(row_vals); % This will yield [] unless the content of row_val is exclusively numeric
        if any(isnan(row_nums));
            targetstruct.(fn{i}) = row_vals;
        else
            targetstruct.(fn{i}) = row_nums;
        end
    end
    targetstruct.Angulation_midslice = targetstruct.Angulation_midslice([3,1,2]);
    targetstruct.Off_Centre_midslice = targetstruct.Off_Centre_midslice([3,1,2]);
    
    function c = strsplit_(s, pattern)
        if nargin == 1
            pattern = '\s+';
        end 
        c = regexp(s,pattern,'split');
    end
end

function row_id = rowfind(TEXT, PATTERN)
    tmp = strfind(TEXT, PATTERN);
    row_id = find(~cellfun('isempty', tmp));
end

function targetstruct = parse_par_tabular_text(targetstruct, partext)
    selection_beg = rowfind(partext, '= IMAGE INFORMATION =') + 1;
    selection_end = rowfind(partext, '= END OF DATA DESCRIPTION FILE =') - 1;
    partabtext = partext(selection_beg:selection_end);
    fn = fieldnames(targetstruct);
    z = ~(cellfun('isempty', partabtext) | cellfun(@(x) strncmp(x,'#', 1), partabtext));
    M = cell2mat(cellfun(@(x) str2double(strsplit(x)), partabtext(z), 'unif', 0));
    cursor = 0;
    for i = 1:length(fn)
        n = targetstruct(2).(fn{i});
        a = cursor + 1;
        b = cursor + n;
        targetstruct(1).(fn{i}) = M(:,a:b);
        cursor = cursor + n;
    end
    targetstruct(1).image_angulation = targetstruct(1).image_angulation(:,[2,3,1]);
    targetstruct(1).image_offcentre = targetstruct(1).image_offcentre(:,[2,3,1]);
    
    if cursor < size(M, 2);
      warning('load_parrec:parse_par_tabular_text:dataMismatch', 'Some columns of the image data were not parsed. The assignment of columns to fields may be incorrect.');
      choice = input('Continue anyway? y/[n]: ');
      if ~strcmpi(choice, 'y')
          fprintf('\n\nBailing out... check that definitions in load_parrec:init_info_struct() is compatible with your PAR file format.\n\n');
      end
    end
end

function version = get_RIET_version(partext)
% GET_RIET_VERSION Parse PAR for Research Image Export Tool version.
    selection_beg = rowfind(partext, '= DATA DESCRIPTION FILE =') + 1;
    selection_end = rowfind(partext, '= GENERAL INFORMATION =') - 1;
    parheadertext = partext(selection_beg:selection_end);
    index = regexp(parheadertext, 'V[0-9]+\.[0-2]+$');
    z = ~cellfun('isempty', index);
    if nnz(z) > 1;
        error('load_parrec:get_REIT_version:tooManyMatches', 'Could not identify version number (too many strings match regexp).');
    elseif nnz(z) == 0
        error('load_parrec:get_REIT_version:noMatch', 'Could not identify version number (no strings match regexp).');
    end
    i = index{z};
    version = parheadertext{z}(i:end);
end

function [dynamics, ic] = info_by_configuration(image_info, version)
    switch lower(version)
    case 'v4.2'
        vary_by_cfg = {...
            'echo_number','image_type_mr','scanning_sequence','image_pixel_size', ...
            'recon_resolution','rescale_intercept','rescale_slope', ...
            'scale_slope','image_angulation','slice_thickness', ...
            'slice_gap','image_display_orientation','slice_orientation', ...
            'image_type_ed_es','pixel_spacing','echo_time', ...
            'trigger_time','diffusion_b_factor','number_of_averages', ...
            'image_flip_angle','minimum_RR','maximum_RR', ...
            'TURBO_factor','Inversion_delay','diffusion_b', ...
            'gradient_orientation_number','contrast_type','diffusion_anisotropy_type', ...
            'label_type'};
        vary_by_img = {...
            'slice_number','dynamic_scan_number', ...
            'cardiac_phase_number','index_in_REC_file','scan_percentage', ...
            'window_center','window_width','image_offcentre', ...
            'fmri_status_indication','dyn_scan_begin_time','cardiac_frequency', ...
            'diffusion'};
        image_info_vary_by_cfg = rmfield(image_info, vary_by_img);
        fn = fieldnames(image_info_vary_by_cfg);
        if numel(fn) ~= numel(vary_by_cfg) && ~all(strcmp(sort(fn(:)), sort(vary_by_cfg(:))))
            error('load_parrec:info_by_configuration:programmingError','The set of factors that define different scanner configurations does not match expectations. Inspect the code.');
        end
    otherwise
        error('load_parrec:info_by_configuration:unknownVersion', 'Unknown PAR/REC version.');
    end
    
    M = struct2array(image_info_vary_by_cfg(1));
    n = struct2array(image_info_vary_by_cfg(2));
    [U, ~, ic] = unique(M, 'rows');
    dynamics = cell2struct(mat2cell(U, ones(size(U,1),1), n), fn, 2);
    n_slice = max(image_info(1).slice_number);
    n_volume = nnz(image_info(1).slice_number == n_slice);
    slice_ind = reshape(ic, n_slice, n_volume);
    volume_ind = slice_ind(1,:);
    for i = 1:numel(dynamics)
        dynamics(i).data = [];
        dynamics(i).run = [];
        dynamics(i).volume_index = find(volume_ind==i);
    end
end