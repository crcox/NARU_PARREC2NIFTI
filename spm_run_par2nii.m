function out = spm_run_par2nii(job)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2005-2014 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_coreg.m 5956 2014-04-16 14:34:25Z guillaume $

%-Converty PAR to NII
%--------------------------------------------------------------------------

outdir = job.outdir{1};
for j = 1:numel(job.dataset)
    protocol_code = fieldnames(job.dataset(j).stype);
    p = protocol_code{1};
    switch upper(p)
        case 'EPI_DE' % Only case where job.dataset.echo_labels will be defined
            parfiles = job.dataset(j).parfiles;
            protocol_label = job.dataset(j).stype.(p).label;
            echo_labels = job.dataset(j).stype.(p).echo_labels;
            convert_par_to_nii( parfiles, outdir, protocol_label, echo_labels);
        otherwise
            parfiles = job.dataset(j).parfiles;
            protocol_label = job.dataset(j).stype.(p).label;
            convert_par_to_nii( job.dataset(j).parfiles, outdir, job.dataset(j).stype.(p));
    end
end
out = [];
