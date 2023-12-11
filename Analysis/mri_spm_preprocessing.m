%__________________________________________________________________________
function mri_spm_preprocessing(mypath,mri)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-INITIALIZE FILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Loop across subjects
nsubjs = length(mri.subject.list);
for cs = 1:nsubjs
    subj = mri.subject.list{cs};
    
    %-Create anatomical, functional, and temporary folders
    %-The temporary folder is for saving batch, functional and distortion
    %correction files. They will be deleted after the analysis.
    rawdir = eval(mypath.folder.subj.mri.raw);
    anatdir = eval(mypath.folder.subj.mri.anatomical);
    funcdir = eval(mypath.folder.subj.mri.functional);
    mytempdir = eval(mypath.folder.temp.subj);
    mkdir(anatdir);
    mkdir(funcdir);
    mkdir(mytempdir);

    if strcmp(mri.spm.distortion,'fieldmap')
        %-Copy  distortion correction files to mytempdir
        str.folder = rawdir;
        str.file = mri.file.fieldmap;
        fn.orig = strcat(str.folder,'/',str.file);
        fn.dest = repelem({mytempdir},length(mri.file.fieldmap));
        cellfun(@copyfile,fn.orig,fn.dest);
    end
    
    %-Copy anatomical file to anatdir
    str.folder = rawdir;
    str.file = mri.file.anatomical;
    fn.orig = strcat(str.folder,'/',str.file);
    fn.dest = repelem({anatdir},length(mri.file.anatomical));
    cellfun(@copyfile,fn.orig,fn.dest);
    
    %-Copy  functional to mytempdir
    str.folder = rawdir;
    str.file = mri.file.functional;
    fn.orig = strcat(str.folder,'/',str.file);
    fn.dest = repelem({mytempdir},length(mri.file.functional));
    cellfun(@copyfile,fn.orig,fn.dest);

    if strcmp(mri.spm.distortion,'fieldmap')
        %-Read the distortion correction files
        magnitudefile = cellstr(spm_select('ExtFPList',mytempdir,'magnitude.nii',Inf));
        phasefile = cellstr(spm_select('ExtFPList',mytempdir,'phase.nii',Inf));
        if isequal(magnitudefile{1},'')|isequal(phasefile{1},'')
            warning('No distortion correction file (phase or magnitude) found!')
            return
        end
    end
    
    %-Read the anatomical file
    anatfile = cellstr(spm_select('FPList',anatdir,'\.nii'));
    if isequal(anatfile{1},'')
        warning('No anatomical file found!')
        return
    end

    %-Read the functional files
    nruns = length(mri.file.functional);
    funcfiles = cell(nruns,1);
    for r = 1:nruns
        %funcfiles{r} = cellstr(spm_select('ExtFPList',mytempdir,mri.file.functional{r},Inf));
        funcfiles{r} = cellstr(spm_select('ExtFPList',mytempdir,mri.file.functional{r},mri.spm.start_vol(r):mri.spm.end_vol(r)));
    end
    if isequal(funcfiles{1},'')
        warning('No functional file found!')
        return
    end
    
    %-Set the ROIs for denoising
    roifiles = {
                strcat(anatdir,'/','rc2anat.nii');
                strcat(anatdir,'/','rc3anat.nii')
                };

    clear matlabbatch;
    disp '********** Creating preprocessing job **********'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-PREPROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%==========================================================================
%-With FIELDMAP distortion correction
%==========================================================================
if strcmp(mri.spm.distortion,'fieldmap')

%-Preprocessing Steps
%--------------------------------------------------------------------------
%https://andysbrainbook.readthedocs.io/en/latest/SPM/SPM_Short_Course/SPM_04_Preprocessing.html

%-1. slice timing
%-2. distortion correction
%-3. realignement
%-4. coregistration: structural image is co-registered with EPI images
%-5. segmentation: the normalized structural image is segmented into
%different classes of tissue using probabilities maps
%-6. normalization: the parameters from segmentation are applied to the
%EPI images for normalization
%-7. smoothing

%-SPM Image File Prefixes
%--------------------------------------------------------------------------
%-a - slice timing correction
%-r - resliced (this can be from coregistration or realignment)
%-u - undistorted (from Realign unwarp - which requires reslicing)
%-w - warped (typically this is done by normilization)
%-s - smoothed
%-c - a tissue class created by segment (c1 is Gray matter, for example,
%mostly for anatomy)

%-Slice-Timing Correction (Temporal Preprocessing)
%==========================================================================
%-Goal: Correct for different acquisition time of each slice within an
%image volume
%-How: All voxel time series are aligned to acquisition time of 1 slice
%via Sinc-interpolation of each voxels time series
%-Input: fmri.nii
%-Output: afmri.nii
stage = 1;
stage_slicetiming = stage;
matlabbatch{stage}.spm.temporal.st.scans = {funcfiles{:}};
matlabbatch{stage}.spm.temporal.st.nslices = mri.spm.nslices;
matlabbatch{stage}.spm.temporal.st.tr = mri.spm.TR;
matlabbatch{stage}.spm.temporal.st.ta = mri.spm.TA;
matlabbatch{stage}.spm.temporal.st.so = mri.spm.slice_order.time;
matlabbatch{stage}.spm.temporal.st.refslice = mri.spm.refslice.time;
matlabbatch{stage}.spm.temporal.st.prefix = 'a_';

%-Calculate VDM
%==========================================================================
%-Input: fmri.nii, phase.nii, magnitude.nii
%-Output: ufmri.nii, wfmag_fmri.nii, afmri_uw.nii, scphase.nii,
%fpm_scphase.nii, vdm5_scphase.nii
stage = stage + 1;
stage_vdm = stage;
matlabbatch{stage}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase = phasefile;
matlabbatch{stage}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = magnitudefile;
matlabbatch{stage}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = [mri.spm.short_TE mri.spm.long_TE];
matlabbatch{stage}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 0;
matlabbatch{stage}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = -1;
matlabbatch{stage}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = mri.spm.total_readout_time;
matlabbatch{stage}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 0;
matlabbatch{stage}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm = 0;
matlabbatch{stage}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = 'Mark3D';
matlabbatch{stage}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm = 10;
matlabbatch{stage}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad = 0;
matlabbatch{stage}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws = 1;
matlabbatch{stage}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = {strcat(mri.spm.path,'/toolbox/FieldMap/T1.nii')};
matlabbatch{stage}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm = 5;
matlabbatch{stage}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode = 2;
matlabbatch{stage}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate = 4;
matlabbatch{stage}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh = 0.5;
matlabbatch{stage}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg = 0.02;
%-Input
%..........................................................................
for r=1:nruns
    matlabbatch{stage}.spm.tools.fieldmap.calculatevdm.subj.session(r).epi = {funcfiles{r,1}{1,1}};
end
%..........................................................................
matlabbatch{stage}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;
matlabbatch{stage}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'session';
matlabbatch{stage}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 0;
matlabbatch{stage}.spm.tools.fieldmap.calculatevdm.subj.anat = '';
matlabbatch{stage}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 0;

%-Realignment and Unwarp (Spatial Preprocessing)
%==========================================================================
%-Goal: Correct for subject motion between volumes by minimising
%mean-squared difference
%-How: Rigid-body transformation
%-Input: afmri.nii
%-Output: uafmri.nii, meanuafmri.nii, rp_afmri.txt, ranat.nii
stage = stage + 1;
stage_realign = stage;
%-Input
%..........................................................................
for r=1:nruns
    matlabbatch{stage}.spm.spatial.realignunwarp.data(r).scans = ...
        cfg_dep(strcat('Slice Timing: Slice Timing Corr. Images (Sess ',num2str(r),')'), substruct('.','val','{}',{stage_slicetiming},'.','val','{}',{1},'.','val','{}',{1}),substruct('()',{r},'.','files'));
    matlabbatch{stage}.spm.spatial.realignunwarp.data(r).pmscan = ...
        cfg_dep(strcat('Calculate VDM: Voxel displacement map (Subj 1, Session ',num2str(r),')'), substruct('.','val', '{}',{stage_vdm}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','vdmfile', '{}',{r}));
end
%..........................................................................
matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.quality = 1;
matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.sep = 2;
matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.rtm = 1;
matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.einterp = 2;
matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.weight = '';
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.jm = 0;
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.sot = [];
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.rem = 1;
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.noi = 5;
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
matlabbatch{stage}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
matlabbatch{stage}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
matlabbatch{stage}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
matlabbatch{stage}.spm.spatial.realignunwarp.uwroptions.mask = 1;
matlabbatch{stage}.spm.spatial.realignunwarp.uwroptions.prefix = 'u_';

%-Co-registration - Estimate and Reslice (Spatial Preprocessing)
%==========================================================================
%-Goal: Match geometry of functional and structural images from same subject
%-How: Find affine transformation (rotation/translation/shear/scaling)
%that maximizes mutual information (similarity) between both images
%-Input: meanuafmri.nii (Reference), struct.nii (Source)
%-Output
stage = stage + 1;
stage_coregister = stage;
%-If you are normalising the data you don't need to reslice as this
%writing will be done later
% Fixed image (Reference) --functional image--
matlabbatch{stage}.spm.spatial.coreg.estwrite.ref(1) = ...
    cfg_dep('Realign & Unwarp: Unwarped Mean Image',substruct('.','val','{}',{stage_realign},'.','val','{}',{1},'.','val','{}',{1},'.','val','{}',{1}),substruct('.','meanuwr'));
% Image that is transformed (Source) --structural image--
matlabbatch{stage}.spm.spatial.coreg.estwrite.source = anatfile;
matlabbatch{stage}.spm.spatial.coreg.estwrite.other = {''};
matlabbatch{stage}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{stage}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{stage}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{stage}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{stage}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{stage}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{stage}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{stage}.spm.spatial.coreg.estwrite.roptions.prefix = 'r_';

%-Segmentation (Spatial Preprocessing)
%==========================================================================
%-Goal: Match geometry of subject brain to standard space (for group studies)
%-How: Find non-linear transformation (deformation field) that makes
%tissue class distribution in structural image most plausible
%-Input
%-Output: mstruct.nii, c1-c6struct.nii, mwc1-mwc6struct.nii,
%rc1-rc6struct.nii, wc1-wc6struct.nii, y_struct.nii
stage = stage + 1;
stage_segmentation = stage;
matlabbatch{stage}.spm.spatial.preproc.channel.vols(1) = ...
    cfg_dep('Coregister: Estimate: Coregistered Images',substruct('.','val','{}',{stage_coregister},'.','val','{}',{1},'.','val','{}',{1},'.','val','{}',{1}),substruct('.','cfiles'));
matlabbatch{stage}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{stage}.spm.spatial.preproc.channel.biasfwhm = 60;
%matlabbatch{stage}.spm.spatial.preproc.channel.write = [0 0];   %-Save Bias Corrected: Save Nothing
matlabbatch{stage}.spm.spatial.preproc.channel.write = [0 1];   %-Save Bias Corrected: Save Bias Corrected
matlabbatch{stage}.spm.spatial.preproc.tissue(1).tpm = {strcat(mri.spm.path,'/tpm/TPM.nii,1')};
matlabbatch{stage}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{stage}.spm.spatial.preproc.tissue(1).native = [1 1];
matlabbatch{stage}.spm.spatial.preproc.tissue(1).warped = [1 1];
matlabbatch{stage}.spm.spatial.preproc.tissue(2).tpm = {strcat(mri.spm.path,'/tpm/TPM.nii,2')};
matlabbatch{stage}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{stage}.spm.spatial.preproc.tissue(2).native = [1 1];
matlabbatch{stage}.spm.spatial.preproc.tissue(2).warped = [1 1];
matlabbatch{stage}.spm.spatial.preproc.tissue(3).tpm = {strcat(mri.spm.path,'/tpm/TPM.nii,3')};
matlabbatch{stage}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{stage}.spm.spatial.preproc.tissue(3).native = [1 1];
matlabbatch{stage}.spm.spatial.preproc.tissue(3).warped = [1 1];
matlabbatch{stage}.spm.spatial.preproc.tissue(4).tpm = {strcat(mri.spm.path,'/tpm/TPM.nii,4')};
matlabbatch{stage}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{stage}.spm.spatial.preproc.tissue(4).native = [1 1];
matlabbatch{stage}.spm.spatial.preproc.tissue(4).warped = [1 1];
matlabbatch{stage}.spm.spatial.preproc.tissue(5).tpm = {strcat(mri.spm.path,'/tpm/TPM.nii,5')};
matlabbatch{stage}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{stage}.spm.spatial.preproc.tissue(5).native = [1 1];
matlabbatch{stage}.spm.spatial.preproc.tissue(5).warped = [1 1];
matlabbatch{stage}.spm.spatial.preproc.tissue(6).tpm = {strcat(mri.spm.path,'/tpm/TPM.nii,6')};
matlabbatch{stage}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{stage}.spm.spatial.preproc.tissue(6).native = [1 1];
matlabbatch{stage}.spm.spatial.preproc.tissue(6).warped = [1 1];
matlabbatch{stage}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{stage}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{stage}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{stage}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{stage}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{stage}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{stage}.spm.spatial.preproc.warp.write = [0 1];   %-Deformation Fields: Forward
%matlabbatch{stage}.spm.spatial.preproc.warp.write = [1 1];   %-Deformation Fields: Inverse + Forward

%-Normalisation - Write (Spatial Preprocessing)
%==========================================================================
%-Goal: Write out functional/structural image in standard space for
%multi-subject statistical analysis
%-How: Applies estimated deformation fields from Unified Segmentation to
%all functional and structural data
%-Input: uafmri.nii, y_struct.nii (Deformation Field)
%-Output: wuafmri.nii
stage = stage + 1;
stage_normalisation = stage;
matlabbatch{stage}.spm.spatial.normalise.write.subj.def(1) = ...
    cfg_dep('Segment: Forward Deformations', substruct('.','val','{}',{stage_segmentation},'.','val','{}',{1},'.','val','{}',{1}),substruct('.','fordef','()',{':'}));
%-Input
%..........................................................................
for r=1:nruns
    matlabbatch{stage}.spm.spatial.normalise.write.subj.resample(r) = ...
        cfg_dep(strcat('Realign & Unwarp: Unwarped Images (Sess ',num2str(r),')'),substruct('.','val','{}',{stage_realign},'.','val','{}',{1},'.','val','{}',{1},'.','val','{}',{1}),substruct('.','sess','()',{r},'.','uwrfiles'));
end
%..........................................................................
%matlabbatch{stage}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
matlabbatch{stage}.spm.spatial.normalise.write.woptions.bb = NaN(2,3);
matlabbatch{stage}.spm.spatial.normalise.write.woptions.vox = mri.spm.voxel_size;
matlabbatch{stage}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{stage}.spm.spatial.normalise.write.woptions.prefix = 'w_';

%-Smoothing (Spatial Preprocessing)
%==========================================================================
%-Goal: Increase sensitivity by reducing thermal noise and inter-subject
%variability in functional images
%-How: Convolution ("Blurring") with a 3D Gaussian kernel
%-Input: wrafmri.nii
%-Output: swuafmri.nii

%stage = stage + 1;
%stage_smoothing = stage;
%matlabbatch{stage}.spm.spatial.smooth.data(1) = ...
%    cfg_dep(strcat('Normalise: Write: Normalised Images (',subj,')'),substruct('.','val','{}',{stage_normalisation},'.','val','{}',{1},'.','val','{}',{1},'.','val','{}',{1}),substruct('()',{1},'.','files'));
%matlabbatch{stage}.spm.spatial.smooth.fwhm = [mri.spm.fwhm mri.spm.fwhm mri.spm.fwhm];
%matlabbatch{stage}.spm.spatial.smooth.dtype = 0;
%matlabbatch{stage}.spm.spatial.smooth.im = 0;
%matlabbatch{stage}.spm.spatial.smooth.prefix = 's_';

%-Denoise with PhysIO Toolbox
%==========================================================================
for r=1:nruns
    stage = stage + 1;
    stage_denoise = stage;
    matlabbatch{stage}.spm.tools.physio.save_dir = {mytempdir};
    matlabbatch{stage}.spm.tools.physio.log_files.vendor = mri.spm.vendor;
    matlabbatch{stage}.spm.tools.physio.log_files.cardiac = {''};
    matlabbatch{stage}.spm.tools.physio.log_files.respiration = {''};
    matlabbatch{stage}.spm.tools.physio.log_files.scan_timing = {''};
    matlabbatch{stage}.spm.tools.physio.log_files.sampling_interval = [];
    matlabbatch{stage}.spm.tools.physio.log_files.relative_start_acquisition = 0;
    matlabbatch{stage}.spm.tools.physio.log_files.align_scan = 'last';
    matlabbatch{stage}.spm.tools.physio.scan_timing.sqpar.Nslices = mri.spm.nslices;
    matlabbatch{stage}.spm.tools.physio.scan_timing.sqpar.NslicesPerBeat = [];
    matlabbatch{stage}.spm.tools.physio.scan_timing.sqpar.TR = mri.spm.TR;
    matlabbatch{stage}.spm.tools.physio.scan_timing.sqpar.Ndummies = 0;
    matlabbatch{stage}.spm.tools.physio.scan_timing.sqpar.Nscans = mri.spm.nvols(r);
    matlabbatch{stage}.spm.tools.physio.scan_timing.sqpar.onset_slice = mri.spm.refslice.index;
    matlabbatch{stage}.spm.tools.physio.scan_timing.sqpar.time_slice_to_slice = [];
    matlabbatch{stage}.spm.tools.physio.scan_timing.sqpar.Nprep = [];
    matlabbatch{stage}.spm.tools.physio.scan_timing.sync.nominal = struct([]);
    matlabbatch{stage}.spm.tools.physio.preproc.cardiac.modality = 'ECG';
    matlabbatch{stage}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.min = 0.4;
    matlabbatch{stage}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.file = 'initial_cpulse_kRpeakfile.mat';
    matlabbatch{stage}.spm.tools.physio.preproc.cardiac.posthoc_cpulse_select.off = struct([]);
    matlabbatch{stage}.spm.tools.physio.model.output_multiple_regressors = mri.file.mr{r};
    matlabbatch{stage}.spm.tools.physio.model.output_physio = 'physio.mat';
    matlabbatch{stage}.spm.tools.physio.model.orthogonalise = 'none';
    matlabbatch{stage}.spm.tools.physio.model.censor_unreliable_recording_intervals = false;
    matlabbatch{stage}.spm.tools.physio.model.retroicor.no = struct([]);
    matlabbatch{stage}.spm.tools.physio.model.rvt.no = struct([]);
    matlabbatch{stage}.spm.tools.physio.model.hrv.no = struct([]);
    matlabbatch{stage}.spm.tools.physio.model.noise_rois.yes.fmri_files = cellstr(strcat(mytempdir,'/','u_a_',mri.file.functional{r}))';
    matlabbatch{stage}.spm.tools.physio.model.noise_rois.yes.roi_files = roifiles;
    matlabbatch{stage}.spm.tools.physio.model.noise_rois.yes.force_coregister = 'No';
    matlabbatch{stage}.spm.tools.physio.model.noise_rois.yes.thresholds = 0.9;
    matlabbatch{stage}.spm.tools.physio.model.noise_rois.yes.n_voxel_crop = 1;
    matlabbatch{stage}.spm.tools.physio.model.noise_rois.yes.n_components = 5;
    matlabbatch{stage}.spm.tools.physio.model.movement.yes.file_realignment_parameters = cellstr(strcat(mytempdir,'/',mri.file.rp{r}))';
    matlabbatch{stage}.spm.tools.physio.model.movement.yes.order = 6;
    matlabbatch{stage}.spm.tools.physio.model.movement.yes.censoring_method = 'none';
    matlabbatch{stage}.spm.tools.physio.model.movement.yes.censoring_threshold = 0.5;
    matlabbatch{stage}.spm.tools.physio.model.other.no = struct([]);
    matlabbatch{stage}.spm.tools.physio.verbose.level = 2;
    matlabbatch{stage}.spm.tools.physio.verbose.fig_output_file = 'physio.fig';
    matlabbatch{stage}.spm.tools.physio.verbose.use_tabs = false;
end

%-Moving the PostScript File
%==========================================================================
stage = stage + 1;
stage_move = stage;
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.files = {strcat(mypath.folder.code.root,'/','spm_',char(datetime('today','Format','yyyyMMMdd')),'.ps')};
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.action.moveto = {funcdir};

%-Moving and Renaming the Functional Files (u_a_)
%==========================================================================
stage = stage + 1;
stage_move = stage;
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.files = cellstr(strcat(mytempdir,'/','u_a_',mri.file.functional))';
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.action.copyren.copyto = {funcdir};
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.action.copyren.patrep.pattern = 'u_a_';
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.action.copyren.patrep.repl = 'fwhm0_u_a_';
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.action.copyren.unique = false;

%-Moving and Renaming the Functional Files (w_u_a_)
%==========================================================================
stage = stage + 1;
stage_move = stage;
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.files = cellstr(strcat(mytempdir,'/','w_u_a_',mri.file.functional))';
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.action.copyren.copyto = {funcdir};
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.action.copyren.patrep.pattern = 'w_u_a_';
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.action.copyren.patrep.repl = 'fwhm0_w_u_a_';
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.action.copyren.unique = false;

%-Moving the Denoising (Realignment Parameters) Files
%==========================================================================
stage = stage + 1;
stage_move = stage;
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.files = cellstr(strcat(mytempdir,'/',mri.file.rp))';
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.action.copyto = {funcdir};

%-Moving the Denoising (Multiple Regression) Files
%==========================================================================
stage = stage + 1;
stage_move = stage;
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.files = cellstr(strcat(mytempdir,'/',mri.file.mr))';
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.action.copyto = {funcdir};

%-Delete MAT Files
%==========================================================================
stage = stage + 1;
stage_delete = stage;
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.files = cellstr(strcat(funcdir,'/','u_a_',strcat(erase(mri.file.functional,'.nii'),'.mat')))';
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.action.delete = true;

%-Delete Temporary Folder
%==========================================================================
stage = stage + 1;
stage_delete = stage;
matlabbatch{stage}.cfg_basicio.file_dir.dir_ops.dir_move.dir = {mytempdir};
matlabbatch{stage}.cfg_basicio.file_dir.dir_ops.dir_move.action.delete = true;

%==========================================================================
%-Without FIELDMAP distortion correction
%==========================================================================
elseif strcmp(mri.spm.distortion,'none')

%-Preprocessing Steps
%--------------------------------------------------------------------------
%https://andysbrainbook.readthedocs.io/en/latest/SPM/SPM_Short_Course/SPM_04_Preprocessing.html

%-1. slice timing
%-2. realignement
%-3. coregistration: structural image is co-registered registered with EPI
%images
%-4. segmentation: the normalized structural image is segmented into
%different classes of tissue using probabilities maps
%-5. normalization: the parameters from segmentation are applied to the
%EPI images for normalization
%-6. smoothing

%-SPM Image File Prefixes
%--------------------------------------------------------------------------
%-a - slice timing correction
%-r - resliced (this can be from coregistration or realignment)
%-u - undistorted (from Realign unwarp - which requires reslicing)
%-w - warped (typically this is done by normilization)
%-s - smoothed
%-c - a tissue class created by segment (c1 is Gray matter, for example,
%mostly for anatomy)

%-Slice-Timing Correction (Temporal Preprocessing)
%==========================================================================
%-Goal: Correct for different acquisition time of each slice within an
%image volume
%-How: All voxel time series are aligned to acquisition time of 1 slice
%via Sinc-interpolation of each voxel’s time series
%-Input: fmri.nii
%-Output: afmri.nii
stage = 1;
stage_slicetiming = stage;
matlabbatch{stage}.spm.temporal.st.scans = {funcfiles{:}};
matlabbatch{stage}.spm.temporal.st.nslices = mri.spm.nslices;
matlabbatch{stage}.spm.temporal.st.tr = mri.spm.TR;
matlabbatch{stage}.spm.temporal.st.ta = mri.spm.TA;
matlabbatch{stage}.spm.temporal.st.so = mri.spm.slice_order.time;
matlabbatch{stage}.spm.temporal.st.refslice = mri.spm.refslice.time;
matlabbatch{stage}.spm.temporal.st.prefix = 'a_';

%-Realignment and Unwarp (Spatial Preprocessing)
%==========================================================================
%-Goal: Correct for subject motion between volumes by minimising
%mean-squared difference
%-How: Rigid-body transformation
%-Input: afmri.nii
%-Output: uafmri.nii, meanuafmri.nii, rp_afmri.txt, ranat.nii
stage = stage + 1;
stage_realign = stage;
%-Input
%..........................................................................
for r=1:nruns
    matlabbatch{stage}.spm.spatial.realignunwarp.data(r).scans = ...
        cfg_dep(strcat('Slice Timing: Slice Timing Corr. Images (Sess ',num2str(r),')'), substruct('.','val','{}',{stage_slicetiming},'.','val','{}',{1},'.','val','{}',{1}),substruct('()',{r},'.','files'));
    matlabbatch{stage}.spm.spatial.realignunwarp.data(r).pmscan = {};
end
%..........................................................................
matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.quality = 1;
matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.sep = 2;
matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.rtm = 1;
matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.einterp = 2;
matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
matlabbatch{stage}.spm.spatial.realignunwarp.eoptions.weight = '';
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.jm = 0;
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.sot = [];
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.rem = 1;
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.noi = 5;
matlabbatch{stage}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
matlabbatch{stage}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
matlabbatch{stage}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
matlabbatch{stage}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
matlabbatch{stage}.spm.spatial.realignunwarp.uwroptions.mask = 1;
matlabbatch{stage}.spm.spatial.realignunwarp.uwroptions.prefix = 'u_';

%-Co-registration - Estimate and Reslice (Spatial Preprocessing)
%==========================================================================
%-Goal: Match geometry of functional and structural images from same subject
%-How: Find affine transformation (rotation/translation/shear/scaling)
%that maximizes mutual information (similarity) between both images
%-Input: meanuafmri.nii (Reference), struct.nii (Source)
%-Output
stage = stage + 1;
stage_coregister = stage;
%-If you are normalising the data you don't need to reslice as this
%writing will be done later
% Fixed image (Reference) --functional image--
matlabbatch{stage}.spm.spatial.coreg.estwrite.ref(1) = ...
    cfg_dep('Realign & Unwarp: Unwarped Mean Image',substruct('.','val','{}',{stage_realign},'.','val','{}',{1},'.','val','{}',{1},'.','val','{}',{1}),substruct('.','meanuwr'));
% Image that is transformed (Source) --structural image--
matlabbatch{stage}.spm.spatial.coreg.estwrite.source = anatfile;
matlabbatch{stage}.spm.spatial.coreg.estwrite.other = {''};
matlabbatch{stage}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{stage}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{stage}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{stage}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{stage}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{stage}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{stage}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{stage}.spm.spatial.coreg.estwrite.roptions.prefix = 'r_';

%-Segmentation (Spatial Preprocessing)
%==========================================================================
%-Goal: Match geometry of subject brain to standard space (for group
%studies)
%-How: Find non-linear transformation (deformation field) that makes
%tissue class distribution in structural image most plausible
%-Input: rstruct.nii
%-Output: c1-c5rstruct.nii, mrstruct.nii, y_rstruct.nii
stage = stage + 1;
stage_segmentation = stage;
matlabbatch{stage}.spm.spatial.preproc.channel.vols(1) = ...
    cfg_dep('Coregister: Estimate: Coregistered Images',substruct('.','val','{}',{stage_coregister},'.','val','{}',{1},'.','val','{}',{1},'.','val','{}',{1}),substruct('.','cfiles'));
matlabbatch{stage}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{stage}.spm.spatial.preproc.channel.biasfwhm = 60;
%matlabbatch{stage}.spm.spatial.preproc.channel.write = [0 0];   %-Save Bias Corrected: Save Nothing
matlabbatch{stage}.spm.spatial.preproc.channel.write = [0 1];   %-Save Bias Corrected: Save Bias Corrected
matlabbatch{stage}.spm.spatial.preproc.tissue(1).tpm = {strcat(mri.spm.path,'/tpm/TPM.nii,1')};
matlabbatch{stage}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{stage}.spm.spatial.preproc.tissue(1).native = [1 1];
matlabbatch{stage}.spm.spatial.preproc.tissue(1).warped = [1 1];
matlabbatch{stage}.spm.spatial.preproc.tissue(2).tpm = {strcat(mri.spm.path,'/tpm/TPM.nii,2')};
matlabbatch{stage}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{stage}.spm.spatial.preproc.tissue(2).native = [1 1];
matlabbatch{stage}.spm.spatial.preproc.tissue(2).warped = [1 1];
matlabbatch{stage}.spm.spatial.preproc.tissue(3).tpm = {strcat(mri.spm.path,'/tpm/TPM.nii,3')};
matlabbatch{stage}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{stage}.spm.spatial.preproc.tissue(3).native = [1 1];
matlabbatch{stage}.spm.spatial.preproc.tissue(3).warped = [1 1];
matlabbatch{stage}.spm.spatial.preproc.tissue(4).tpm = {strcat(mri.spm.path,'/tpm/TPM.nii,4')};
matlabbatch{stage}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{stage}.spm.spatial.preproc.tissue(4).native = [1 1];
matlabbatch{stage}.spm.spatial.preproc.tissue(4).warped = [1 1];
matlabbatch{stage}.spm.spatial.preproc.tissue(5).tpm = {strcat(mri.spm.path,'/tpm/TPM.nii,5')};
matlabbatch{stage}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{stage}.spm.spatial.preproc.tissue(5).native = [1 1];
matlabbatch{stage}.spm.spatial.preproc.tissue(5).warped = [1 1];
matlabbatch{stage}.spm.spatial.preproc.tissue(6).tpm = {strcat(mri.spm.path,'/tpm/TPM.nii,6')};
matlabbatch{stage}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{stage}.spm.spatial.preproc.tissue(6).native = [1 1];
matlabbatch{stage}.spm.spatial.preproc.tissue(6).warped = [1 1];
matlabbatch{stage}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{stage}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{stage}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{stage}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{stage}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{stage}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{stage}.spm.spatial.preproc.warp.write = [0 1];   %-Deformation Fields: Forward
%matlabbatch{stage}.spm.spatial.preproc.warp.write = [1 1];   %-Deformation Fields: Inverse + Forward

%-Normalisation - Write (Spatial Preprocessing)
%==========================================================================
%-Goal: Write out functional/structural image in standard space for
%multi-subject statistical analysis
%-How: Applies estimated deformation fields from Unified Segmentation to
%all functional and structural data
%-Input: rafmri.nii, y_rstruct.nii (Deformation Field)
%-Output: wrafmri.nii
stage = stage + 1;
stage_normalisation = stage;
matlabbatch{stage}.spm.spatial.normalise.write.subj.def(1) = ...
    cfg_dep('Segment: Forward Deformations', substruct('.','val','{}',{stage_segmentation},'.','val','{}',{1},'.','val','{}',{1}),substruct('.','fordef','()',{':'}));
%-Input
%..........................................................................
for r=1:nruns
    matlabbatch{stage}.spm.spatial.normalise.write.subj.resample(r) = ...
        cfg_dep(strcat('Realign & Unwarp: Unwarped Images (Sess ',num2str(r),')'),substruct('.','val','{}',{stage_realign},'.','val','{}',{1},'.','val','{}',{1},'.','val','{}',{1}),substruct('.','sess','()',{r},'.','uwrfiles'));
end
%..........................................................................
%matlabbatch{stage}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
matlabbatch{stage}.spm.spatial.normalise.write.woptions.bb = NaN(2,3);
matlabbatch{stage}.spm.spatial.normalise.write.woptions.vox = mri.spm.voxel_size;
matlabbatch{stage}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{stage}.spm.spatial.normalise.write.woptions.prefix = 'w_';

%-Smoothing (Spatial Preprocessing)
%==========================================================================
%-Goal: Increase sensitivity by reducing thermal noise and inter-subject
%variability in functional images
%-How: Convolution ("Blurring") with a 3D Gaussian kernel
%-Input: wrafmri.nii
%-Output: swrafmri.nii

%stage = stage + 1;
%stage_smoothing = stage;
%matlabbatch{stage}.spm.spatial.smooth.data(1) = ...
%    cfg_dep(strcat('Normalise: Write: Normalised Images (',subj,')'),substruct('.','val','{}',{stage_normalisation},'.','val','{}',{1},'.','val','{}',{1},'.','val','{}',{1}),substruct('()',{1},'.','files'));
%matlabbatch{stage}.spm.spatial.smooth.fwhm = [mri.spm.fwhm mri.spm.fwhm mri.spm.fwhm];
%matlabbatch{stage}.spm.spatial.smooth.dtype = 0;
%matlabbatch{stage}.spm.spatial.smooth.im = 0;
%matlabbatch{stage}.spm.spatial.smooth.prefix = 's_';

%-Denoise with PhysIO Toolbox
%==========================================================================
for r=1:nruns
    stage = stage + 1;
    stage_denoise = stage;
    matlabbatch{stage}.spm.tools.physio.save_dir = {mytempdir};
    matlabbatch{stage}.spm.tools.physio.log_files.vendor = mri.spm.vendor;
    matlabbatch{stage}.spm.tools.physio.log_files.cardiac = {''};
    matlabbatch{stage}.spm.tools.physio.log_files.respiration = {''};
    matlabbatch{stage}.spm.tools.physio.log_files.scan_timing = {''};
    matlabbatch{stage}.spm.tools.physio.log_files.sampling_interval = [];
    matlabbatch{stage}.spm.tools.physio.log_files.relative_start_acquisition = 0;
    matlabbatch{stage}.spm.tools.physio.log_files.align_scan = 'last';
    matlabbatch{stage}.spm.tools.physio.scan_timing.sqpar.Nslices = mri.spm.nslices;
    matlabbatch{stage}.spm.tools.physio.scan_timing.sqpar.NslicesPerBeat = [];
    matlabbatch{stage}.spm.tools.physio.scan_timing.sqpar.TR = mri.spm.TR;
    matlabbatch{stage}.spm.tools.physio.scan_timing.sqpar.Ndummies = 0;
    matlabbatch{stage}.spm.tools.physio.scan_timing.sqpar.Nscans = mri.spm.nvols(r);
    matlabbatch{stage}.spm.tools.physio.scan_timing.sqpar.onset_slice = mri.spm.refslice.index;
    matlabbatch{stage}.spm.tools.physio.scan_timing.sqpar.time_slice_to_slice = [];
    matlabbatch{stage}.spm.tools.physio.scan_timing.sqpar.Nprep = [];
    matlabbatch{stage}.spm.tools.physio.scan_timing.sync.nominal = struct([]);
    matlabbatch{stage}.spm.tools.physio.preproc.cardiac.modality = 'ECG';
    matlabbatch{stage}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.min = 0.4;
    matlabbatch{stage}.spm.tools.physio.preproc.cardiac.initial_cpulse_select.auto_matched.file = 'initial_cpulse_kRpeakfile.mat';
    matlabbatch{stage}.spm.tools.physio.preproc.cardiac.posthoc_cpulse_select.off = struct([]);
    matlabbatch{stage}.spm.tools.physio.model.output_multiple_regressors = mri.file.mr{r};
    matlabbatch{stage}.spm.tools.physio.model.output_physio = 'physio.mat';
    matlabbatch{stage}.spm.tools.physio.model.orthogonalise = 'none';
    matlabbatch{stage}.spm.tools.physio.model.censor_unreliable_recording_intervals = false;
    matlabbatch{stage}.spm.tools.physio.model.retroicor.no = struct([]);
    matlabbatch{stage}.spm.tools.physio.model.rvt.no = struct([]);
    matlabbatch{stage}.spm.tools.physio.model.hrv.no = struct([]);
    matlabbatch{stage}.spm.tools.physio.model.noise_rois.yes.fmri_files = cellstr(strcat(mytempdir,'/','u_a_',mri.file.functional{r}))';
    matlabbatch{stage}.spm.tools.physio.model.noise_rois.yes.roi_files = roifiles;
    matlabbatch{stage}.spm.tools.physio.model.noise_rois.yes.force_coregister = 'No';
    matlabbatch{stage}.spm.tools.physio.model.noise_rois.yes.thresholds = 0.9;
    matlabbatch{stage}.spm.tools.physio.model.noise_rois.yes.n_voxel_crop = 1;
    matlabbatch{stage}.spm.tools.physio.model.noise_rois.yes.n_components = 5;
    matlabbatch{stage}.spm.tools.physio.model.movement.yes.file_realignment_parameters = cellstr(strcat(mytempdir,'/',mri.file.rp{r}))';
    matlabbatch{stage}.spm.tools.physio.model.movement.yes.order = 6;
    matlabbatch{stage}.spm.tools.physio.model.movement.yes.censoring_method = 'none';
    matlabbatch{stage}.spm.tools.physio.model.movement.yes.censoring_threshold = 0.5;
    matlabbatch{stage}.spm.tools.physio.model.other.no = struct([]);
    matlabbatch{stage}.spm.tools.physio.verbose.level = 2;
    matlabbatch{stage}.spm.tools.physio.verbose.fig_output_file = 'physio.fig';
    matlabbatch{stage}.spm.tools.physio.verbose.use_tabs = false;
end

%-Moving the PostScript File
%==========================================================================
stage = stage + 1;
stage_move = stage;
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.files = {strcat(mypath.folder.code.root,'/','spm_',char(datetime('today','Format','yyyyMMMdd')),'.ps')};
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.action.moveto = {funcdir};

%-Moving and Renaming the Functional Files (u_a_)
%==========================================================================
stage = stage + 1;
stage_move = stage;
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.files = cellstr(strcat(mytempdir,'/','u_a_',mri.file.functional))';
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.action.copyren.copyto = {funcdir};
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.action.copyren.patrep.pattern = 'u_a_';
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.action.copyren.patrep.repl = 'fwhm0_u_a_';
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.action.copyren.unique = false;

%-Moving and Renaming the Functional Files (w_u_a_)
%==========================================================================
stage = stage + 1;
stage_move = stage;
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.files = cellstr(strcat(mytempdir,'/','w_u_a_',mri.file.functional))';
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.action.copyren.copyto = {funcdir};
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.action.copyren.patrep.pattern = 'w_u_a_';
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.action.copyren.patrep.repl = 'fwhm0_w_u_a_';
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.action.copyren.unique = false;

%-Moving the Denoising (Realignment Parameters) Files
%==========================================================================
stage = stage + 1;
stage_move = stage;
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.files = cellstr(strcat(mytempdir,'/',mri.file.rp))';
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.action.copyto = {funcdir};

%-Moving the Denoising (Multiple Regression) Files
%==========================================================================
stage = stage + 1;
stage_move = stage;
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.files = cellstr(strcat(mytempdir,'/',mri.file.mr))';
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.action.copyto = {funcdir};

%-Delete MAT Files
%==========================================================================
stage = stage + 1;
stage_delete = stage;
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.files = cellstr(strcat(funcdir,'/','u_a_',strcat(erase(mri.file.functional,'.nii'),'.mat')))';
matlabbatch{stage}.cfg_basicio.file_dir.file_ops.file_move.action.delete = true;

%-Delete Temporary Folder
%==========================================================================
stage = stage + 1;
stage_delete = stage;
matlabbatch{stage}.cfg_basicio.file_dir.dir_ops.dir_move.dir = {mytempdir};
matlabbatch{stage}.cfg_basicio.file_dir.dir_ops.dir_move.action.delete = true;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-SAVE ONE JOB PER SUBJECT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str.folder = mytempdir;
str.file = eval(mypath.file.subj.jobfile.preprocess);
save(strcat(str.folder,'/',str.file),'matlabbatch');

end

end
