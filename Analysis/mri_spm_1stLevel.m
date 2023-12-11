%__________________________________________________________________________
function mri_spm_1stLevel(mypath,mri,splt,runs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-IMPORT FILES AND INFORMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop across subjects
nsubjs = length(mri.subject.list);
for cs = 1:nsubjs
    
    %-Create 1st-level and temporary folders
    subj = mri.subject.list{cs};
    funcdir = eval(mypath.folder.subj.mri.functional);
    firstdir = eval(mypath.folder.subj.mri.firstLevel);
    mytempdir = eval(mypath.folder.temp.subj);
    mkdir(firstdir);
    mkdir(mytempdir);
    
    %-Read the Functional and Denoising Files
    funcfilesList = cellstr(spm_select('List',funcdir,mri.file.preprocessed(splt)));
    mrfilesList = cellstr(spm_select('List',funcdir,mri.file.denoise(splt)));
    nruns = size(funcfilesList,1);
    if isequal(funcfilesList{1},'')
        warning('No NIFTI file found')
        return
    end
    if isequal(mrfilesList{1},'')
        warning('No realignment parameter found')
        return
    end
    for r = 1:nruns
        funcfiles{r} = cellstr(spm_select('ExtFPList',funcdir,funcfilesList{r},Inf));
        mrfiles{r} = cellstr(spm_select('FPList',funcdir,mrfilesList{r}));
    end

    %-Load Conditions and Contrasts
    [design.conds,design.conts]=mri_design(mypath,mri,subj,runs);
    
    clear matlabbatch;
    disp '***** Creating first-level job *****'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-1ST-LEVEL ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-Model Specification
%==========================================================================
stage = 1;
stage_modelspec = stage;
matlabbatch{stage}.spm.stats.fmri_spec.dir = {mytempdir};
matlabbatch{stage}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{stage}.spm.stats.fmri_spec.timing.RT = mri.spm.TR;
matlabbatch{stage}.spm.stats.fmri_spec.timing.fmri_t = mri.spm.fmri_t;
matlabbatch{stage}.spm.stats.fmri_spec.timing.fmri_t0 = mri.spm.fmri_t0;
%-Input
%..........................................................................
for r=1:nruns
    matlabbatch{stage}.spm.stats.fmri_spec.sess(r).scans = funcfiles{r};
    % ::: Conditions ::: %
    for c=1:design.conds(r).num
        % Condition name
        matlabbatch{stage}.spm.stats.fmri_spec.sess(r).cond(c).name = design.conds(r).name{c};
        % Condition onsets in a run
        matlabbatch{stage}.spm.stats.fmri_spec.sess(r).cond(c).onset = design.conds(r).onset{c};
        % Condition duration in a run (it could be a number of an array based on the desing)
        matlabbatch{stage}.spm.stats.fmri_spec.sess(r).cond(c).duration = design.conds(r).dur{c};
        matlabbatch{stage}.spm.stats.fmri_spec.sess(r).cond(c).tmod = 0;
        matlabbatch{stage}.spm.stats.fmri_spec.sess(r).cond(c).pmod = struct('name',{},'param',{},'poly',{});
        matlabbatch{stage}.spm.stats.fmri_spec.sess(r).cond(c).orth = 1;
    end
    matlabbatch{stage}.spm.stats.fmri_spec.sess(r).multi = {''};
    matlabbatch{stage}.spm.stats.fmri_spec.sess(r).regress = struct('name', {}, 'val', {});
    matlabbatch{stage}.spm.stats.fmri_spec.sess(r).multi_reg = mrfiles{r};
    matlabbatch{stage}.spm.stats.fmri_spec.sess(r).hpf = 128;
end
%..........................................................................
matlabbatch{stage}.spm.stats.fmri_spec.fact = struct('name',{},'levels',{});
matlabbatch{stage}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{stage}.spm.stats.fmri_spec.volt = 1;
matlabbatch{stage}.spm.stats.fmri_spec.global = mri.spm.global;
matlabbatch{stage}.spm.stats.fmri_spec.mthresh = mri.spm.mthresh;
matlabbatch{stage}.spm.stats.fmri_spec.mask = {''};
matlabbatch{stage}.spm.stats.fmri_spec.cvi = mri.spm.cvi;

%-Model Estimation
%==========================================================================
stage = stage + 1;
stage_modelest = stage;
matlabbatch{stage}.spm.stats.fmri_est.spmmat(1) = ...
    cfg_dep('fMRI model specification: SPM.mat File',substruct('.','val','{}',{stage_modelspec},'.','val','{}',{1},'.','val','{}',{1}),substruct('.','spmmat'));
matlabbatch{stage}.spm.stats.fmri_est.write_residuals = 1;
matlabbatch{stage}.spm.stats.fmri_est.method.Classical = 1;

%-Contrast Manager
%==========================================================================
stage = stage + 1;
stage_contrast = stage;
matlabbatch{stage}.spm.stats.con.spmmat(1) = ...
    cfg_dep('Model estimation: SPM.mat File',substruct('.','val','{}',{stage_modelest},'.','val','{}',{1},'.','val','{}',{1}),substruct('.','spmmat'));
%-Contrasts
%..........................................................................
for c=1:design.conts.num
    % Contrast name
    matlabbatch{stage}.spm.stats.con.consess{c}.tcon.name = design.conts.name{c};
    % Contrast weight
    matlabbatch{stage}.spm.stats.con.consess{c}.tcon.weights = design.conts.weight{c};
    matlabbatch{stage}.spm.stats.con.consess{c}.tcon.sessrep = 'repl';
end
%..........................................................................
matlabbatch{stage}.spm.stats.con.delete = 0;

%-Results Report
%==========================================================================
stage = stage + 1;
stage_result = stage;
matlabbatch{stage}.spm.stats.results.spmmat(1) = ...
    cfg_dep('Contrast Manager: SPM.mat File',substruct('.','val','{}',{stage_contrast},'.','val','{}',{1},'.','val','{}',{1}),substruct('.','spmmat'));
%-Contrasts
%..........................................................................
for c=1:design.conts.num
    matlabbatch{stage}.spm.stats.results.conspec(c).titlestr = design.conts.name{c};
    % Contrast index
    matlabbatch{stage}.spm.stats.results.conspec(c).contrasts = c;
    matlabbatch{stage}.spm.stats.results.conspec(c).threshdesc = mri.spm.result.thresholdType;
    matlabbatch{stage}.spm.stats.results.conspec(c).thresh = mri.spm.result.threshold;
    matlabbatch{stage}.spm.stats.results.conspec(c).extent = mri.spm.result.extent;
    if strcmp(mri.spm.result.mask,'none')
        matlabbatch{stage}.spm.stats.results.conspec(c).mask.none = 1;
    end
end
%..........................................................................
matlabbatch{stage}.spm.stats.results.units = 1;
matlabbatch{stage}.spm.stats.results.print = 'ps';
matlabbatch{stage}.spm.stats.results.write.none = 1;

%-Moving the PostScript File
%==========================================================================
stage = stage + 1;
stage_move = stage;
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.files = {strcat(mytempdir,'/','spm_',char(datetime('today','Format','yyyyMMMdd')),'.ps')};
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.action.moveto = {firstdir};

%-Moving the 1st-Level Files (con files)
%==========================================================================
stage = stage + 1;
stage_move = stage;
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.files = cellstr(strcat(mytempdir,'/',mri.file.(mri.task.name).con.nii.values))';
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.action.moveto = {firstdir};

%-Moving the 1st-Level Files (spmT files)
%==========================================================================
stage = stage + 1;
stage_move = stage;
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.files = cellstr(strcat(mytempdir,'/',mri.file.(mri.task.name).spmt.nii.values)');
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.action.moveto = {firstdir};

%-Change Directory
%==========================================================================
stage = stage + 1;
stage_change = stage;
matlabbatch{stage}.cfg_basicio.file_dir.dir_ops.cfg_cd.dir = {mypath.folder.code.root};

%-Delete Temporary Folder
%==========================================================================
stage = stage + 1;
stage_delete = stage;
matlabbatch{stage}.cfg_basicio.file_dir.dir_ops.dir_move.dir = {mytempdir};
matlabbatch{stage}.cfg_basicio.file_dir.dir_ops.dir_move.action.delete = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-SAVE ONE JOB PER SUBJECT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str.folder = mytempdir;
str.file = eval(mypath.file.subj.jobfile.firstLevel);
save(strcat(str.folder,'/',str.file),'matlabbatch');

end

end
