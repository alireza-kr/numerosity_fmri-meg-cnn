%__________________________________________________________________________
function mri_cosmo_group(mypath,mri,mvpa,mode)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-ANALYSIS MODE - GLM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode,'glm-surface-tfce')

%-IMPORT INFORMATION
%==========================================================================
nsubjs = length(mri.subject.list);

targets=ones(1,nsubjs)';
chunks=(1:nsubjs)';

%-MERGE FILES AND ANALYSIS
%==========================================================================
for p=1:length(mvpa.rdm.predictor.name)
    
    %-LOAD DATA SET
    %======================================================================
    ds_cell = {};
    for s=1:nsubjs
        subj = mri.subject.list{s};
        
        str.file = eval(mypath.file.subj.mri.glm.searchlight);
        str.folder = eval(mypath.folder.result.mri.searchlight);
        fn.in = strcat(str.folder,'/',str.file);
        ds_cell{s} = cosmo_surface_dataset(fn.in,'targets',1,'chunks',s);
    end

    %https://gitlab.pavlovia.org/tgro/2020_ab/blob/0d4320498d6f675e86fd4048b06e5e99a84748f1/get_acti.m
    ds_stacked = cosmo_stack(ds_cell);
    ds_stacked.sa.targets = targets;
    ds_stacked.sa.chunks = chunks;
    
    %-Load standard surface (MNI or Talairach)
    str.file = eval(mypath.file.surface.intermediate);
    str.folder = eval(mypath.folder.data.surface);
    fn.intermediate = strcat(str.folder,'/',str.file);
    [vertices,faces] = surfing_read(fn.intermediate);

    %-ANALYSIS
    %======================================================================
    %-Fisher transform values
    if strcmp(mvpa.cosmo.glm.analysis.mri,'corr') || strcmp(mvpa.cosmo.glm.analysis.mri,'semipartialcorr')
        ds_stacked.samples = atanh(ds_stacked.samples);
    end
    
    %-Define neighborhood for each feature
    cluster_nbrhood=cosmo_cluster_neighborhood(ds_stacked,'vertices',vertices,'faces',faces);

    fprintf('Cluster neighborhood:\n');
    cosmo_disp(cluster_nbrhood);

    opt=struct();
    opt.niter=mvpa.cosmo.group.niter;
    opt.h0_mean=mvpa.cosmo.group.h0_mean;
    opt.null=[];

    fprintf('Running multiple-comparison correction with these options:\n');
    cosmo_disp(opt);

    tfce_z_ds = cosmo_montecarlo_cluster_stat(ds_stacked,cluster_nbrhood,opt);
    
    %-Convert Z-Score to P-Value
    %......................................................................
    %tfce_p_ds = tfce_z_ds;
    
    %-Right-tailed
    %tfce_p_ds.samples = 1-normcdf(tfce_z_ds.samples);
    %-Left-tailed
    %tfce_p_ds.samples = normcdf(tfce_z_ds.samples);
    %-Two-tailed
    %tfce_p_ds.samples = 2*(1-normcdf(tfce_z_ds.samples));
    %......................................................................

    %-SAVE RESULTS
    %======================================================================
    fprintf('TFCE p-value dataset\n');
    cosmo_disp(tfce_z_ds);

    str.file = eval(mypath.file.result.tfce.mri.glm);
    str.folder = eval(mypath.folder.result.mri.group);
    fn.out.dset = strcat(str.folder,'/',str.file);
    cosmo_map2surface(tfce_z_ds,fn.out.dset);
    %cosmo_map2surface(tfce_z_ds,fn.out.gifti);
    %cosmo_map2surface(tfce_z_ds,fn.out.smp);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-ANALYSIS MODE - FIRST-LEVEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode,'second-level-surface-tfce')

%-IMPORT INFORMATION
%==========================================================================
nsubjs = length(mri.subject.list);

targets=ones(1,nsubjs)';
chunks=(1:nsubjs)';
    
%-LOAD DATA SET
%======================================================================
ds_cell = {};
for s=1:nsubjs
    subj = mri.subject.list{s};

    str.file = eval(mypath.file.subj.mri.secondLevel.dset);
    str.folder = eval(mypath.folder.subj.mri.firstLevel);
    fn.in = strcat(str.folder,'/',str.file);
    ds_cell{s} = cosmo_surface_dataset(fn.in,'targets',1,'chunks',s);
end

%https://gitlab.pavlovia.org/tgro/2020_ab/blob/0d4320498d6f675e86fd4048b06e5e99a84748f1/get_acti.m
ds_stacked = cosmo_stack(ds_cell);
ds_stacked.sa.targets = targets;
ds_stacked.sa.chunks = chunks;

%-Load standard surface (MNI or Talairach)
str.file = eval(mypath.file.surface.intermediate);
str.folder = eval(mypath.folder.data.surface);
fn.intermediate = strcat(str.folder,'/',str.file);
[vertices,faces] = surfing_read(fn.intermediate);

%-ANALYSIS
%======================================================================
%-Define neighborhood for each feature
cluster_nbrhood=cosmo_cluster_neighborhood(ds_stacked,'vertices',vertices,'faces',faces);

fprintf('Cluster neighborhood:\n');
cosmo_disp(cluster_nbrhood);

opt=struct();
opt.niter=mvpa.cosmo.group.niter;
opt.h0_mean=mvpa.cosmo.group.h0_mean;
opt.null=[];

fprintf('Running multiple-comparison correction with these options:\n');
cosmo_disp(opt);

tfce_z_ds = cosmo_montecarlo_cluster_stat(ds_stacked,cluster_nbrhood,opt);

%-Convert Z-Score to P-Value
%......................................................................
%tfce_p_ds = tfce_z_ds;

%-Right-tailed
%tfce_p_ds.samples = 1-normcdf(tfce_z_ds.samples);
%-Left-tailed
%tfce_p_ds.samples = normcdf(tfce_z_ds.samples);
%-Two-tailed
%tfce_p_ds.samples = 2*(1-normcdf(tfce_z_ds.samples));
%......................................................................

%-SAVE RESULTS
%======================================================================
fprintf('TFCE p-value dataset\n');
cosmo_disp(tfce_z_ds);

str.file = eval(mypath.file.result.tfce.mri.secondLevel);
str.folder = eval(mypath.folder.result.mri.group);
fn.out.dset = strcat(str.folder,'/',str.file);
cosmo_map2surface(tfce_z_ds,fn.out.dset);
%cosmo_map2surface(tfce_z_ds,fn.out.gifti);
%cosmo_map2surface(tfce_z_ds,fn.out.smp);
    
end

end
