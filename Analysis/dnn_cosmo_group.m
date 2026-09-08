%__________________________________________________________________________
function dnn_cosmo_group(mypath,mri,meg,mvpa,type,layer)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-DNN @ MRI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(type,strcat(mvpa.dnn.architecture,'-',mvpa.dnn.database,'@mri'))

%-IMPORT INFORMATION
%==========================================================================
nsubjs = length(mri.subject.list);

targets=ones(1,nsubjs)';
chunks=(1:nsubjs)';

mvpa.dnn.predictor.name = mvpa.dnn.layer.name;

%-MERGE FILES AND ANALYSIS
%==========================================================================
for p=1:length(mvpa.dnn.predictor.name)
    
    %-LOAD DATA SET
    %======================================================================
    ds_cell = {};
    for s=1:nsubjs
        subj = mri.subject.list{s};
        
        str.file = eval(mypath.file.subj.dnn.glm.searchlight);
        str.folder = eval(mypath.folder.result.dnn.searchlight);
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
    ds_stacked.samples = atanh(ds_stacked.samples);
    
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

    str.file = eval(mypath.file.result.tfce.dnn.glm.mri);
    str.folder = eval(mypath.folder.result.dnn.group);
    fn.out.dset = strcat(str.folder,'/',str.file);
    cosmo_map2surface(tfce_z_ds,fn.out.dset);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Model & DNN @ MRI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexpcmp(type,strcat('model-',mvpa.dnn.architecture,'-',mvpa.dnn.database,'-.*@mri'))

%-IMPORT INFORMATION
%==========================================================================
nsubjs = length(mri.subject.list);

targets=ones(1,nsubjs)';
chunks=(1:nsubjs)';

mvpa.dnn.predictor.name = [mvpa.rdm.predictor.name,mvpa.dnn.model.name(layer)];

%-MERGE FILES AND ANALYSIS
%==========================================================================
for p=1:length(mvpa.dnn.predictor.name)
    
    %-LOAD DATA SET
    %======================================================================
    ds_cell = {};
    for s=1:nsubjs
        subj = mri.subject.list{s};
        
        str.file = eval(mypath.file.subj.dnn.glm.searchlight);
        str.folder = eval(mypath.folder.result.dnn.searchlight);
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
    ds_stacked.samples = atanh(ds_stacked.samples);
    
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

    str.file = eval(mypath.file.result.tfce.dnn.glm.mri);
    str.folder = eval(mypath.folder.result.dnn.group);
    fn.out.dset = strcat(str.folder,'/',str.file);
    cosmo_map2surface(tfce_z_ds,fn.out.dset);
end

end

end
