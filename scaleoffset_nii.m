function r = scaleoffset_nii(protocol, meta)
    voxel_dim = [protocol.pixel_spacing,protocol.slice_thickness + protocol.slice_gap];
    volume_dim = [protocol.recon_resolution(1), protocol.recon_resolution(2), meta.Max_number_of_slices];
    stackOffcentre = meta.Off_Centre_midslice;

    %% SETUP THE INITIAL ROTATION MATRIX
    % NB. The following are adapted from dcm2niix:nii_dicom.cpp#1523:1585
    %
    % First, convert degrees to radians and compute sine and cosine of the
    % angles.
    ca = [...
        cos(deg2rad(meta.Angulation_midslice(1))),...
        cos(deg2rad(meta.Angulation_midslice(2))),...
        cos(deg2rad(meta.Angulation_midslice(3)))...
    ];
    sa = [...
        sin(deg2rad(meta.Angulation_midslice(1))),...
        sin(deg2rad(meta.Angulation_midslice(2))),...
        sin(deg2rad(meta.Angulation_midslice(3)))...
    ];
    rx = [1.0, 0.0, 0.0; 0.0, ca(1), -sa(1); 0.0, sa(1), ca(1)];
    ry = [ca(2), 0.0, sa(2); 0.0, 1.0, 0.0; -sa(2), 0.0, ca(2)];
    rz = [ca(3), -sa(3), 0.0; sa(3), ca(3), 0.0; 0.0, 0.0, 1.0];
    R = (rx * ry) * rz;

    kSliceOrientTra = 1;
    kSliceOrientSag = 2;
    kSliceOrientCor = 3;
    sliceOrient = protocol.slice_orientation;

    ixyz = [1,2,3];
    if sliceOrient == kSliceOrientSag
        ixyz = [2,3,1];
        for r = 1:3
            for c = 1:3
                if (c ~= 1), R(r,c) = -R(r,c), end; % invert first and final columns
            end
        end
    elseif sliceOrient == kSliceOrientCor
        ixyz = [1,3,2];
        for r = 1:3
            R(r,3) = -R(r,3); % invert rows of final column
        end
    end

    R = R(:,ixyz); % dicom rotation matrix

    orient = zeros(1,7);
    orient(2) = R(1,1);
    orient(3) = R(2,1);
    orient(4) = R(3,1);
    orient(5) = R(1,2);
    orient(6) = R(2,2);
    orient(7) = R(3,2);

    dim_mat = zeros(3,3);
    dim_mat([1,5,9]) = voxel_dim;
    R = R * dim_mat;
    R44 = eye(4);
    R44(1:3,1:4) = [R, stackOffcentre'];

    offset_mat = eye(4);
    offset_mat(1:3,4) = (volume_dim - 1) / 2;
    % offset_mat_inv = inv(offset_mat);
    R44 = R44 / offset_mat; % same as R44 * inv(offset_mat)
    y_v = [0,0,volume_dim(3) - 1, 1];
    y_v = R44 * y_v';

    iOri = 3; % for axial, slices are 3rd dimenson (indexed from 0) (k)
    if (sliceOrient == kSliceOrientSag), iOri = 1, end; % for sagittal, slices are 1st dimension (i)
    if (sliceOrient == kSliceOrientCor), iOri = 2, end; % for coronal, slices are 2nd dimension (j)
    if  (( (y_v(iOri)-R44(iOri,4))>0 ) == ( (y_v(iOri)-stackOffcentre(iOri)) > 0 ) )
        patientPosition = R44(1:3,4);
        patientPositionLast = y_v(1:3);
    else
        patientPosition = y_v(1:3);
        patientPositionLast = R44.m(1:3,4);
    end

    r = struct(...
        'voxel_dim', voxel_dim, ...
        'volume_dim', volume_dim, ...
        'mat', R44, ...
        'patientPosition', patientPosition, ...
        'patientPositionLast', patientPositionLast);
end