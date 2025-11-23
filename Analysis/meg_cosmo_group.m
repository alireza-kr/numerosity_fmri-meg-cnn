%__________________________________________________________________________
function result = meg_cosmo_group(mypath,meg,mvpa,data,mode)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-ANALYSIS MODE - Sensor-level GLM RSA (Time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode,'glm-sensor-glm-time-rsa-tfce')

%-IMPORT FILES AND INFORMATION
%==========================================================================
nsubjs = length(meg.subject.list);

%-ANALYSIS
%==========================================================================
time_axis = linspace(1000*meg.epoch.tmin,1000*meg.epoch.tmax,size(data,2));

%https://www.cosmomvpa.org/faq.html#get-temporal-data-in-a-cosmomvpa-struct
ds=cosmo_flatten(data(:,:),{'time'},{time_axis},2);

%-cluster neighborhood over time points
cl_nh=cosmo_cluster_neighborhood(ds);

%-cluster neighborhood not connecting time points
%cl_nh_not_over_time=cosmo_cluster_neighborhood(ds,'time',false);

ds.sa.targets=ones(nsubjs,1);
ds.sa.chunks=(1:nsubjs)';

%......................................................................
opt=struct();
opt.niter=mvpa.cosmo.group.niter;
opt.h0_mean=mvpa.cosmo.group.h0_mean;
%......................................................................

tfce_z_ds=cosmo_montecarlo_cluster_stat(ds,cl_nh,opt);
%tfce_z_ds=cosmo_montecarlo_cluster_stat(ds,cl_nh_not_over_time,opt);

result = tfce_z_ds.samples;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-ANALYSIS MODE - Sensor-level GLM RSA (Time-Frequency)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode,'glm-sensor-glm-time-freq-rsa-tfce')

%-IMPORT FILES AND INFORMATION
%==========================================================================
nsubjs = length(meg.subject.list);

%-ANALYSIS
%==========================================================================
group_cell=cell(nsubjs,1);
for n=1:nsubjs
    dgm_ds=cosmo_flatten(data{n},{'freq','time'},{1:size(data{n},1),1:size(data{n},2)},1);

    % make freq and time a feature dimension
    ds=cosmo_dim_transpose(dgm_ds,{'freq','time'},2);

    % for one-sample t-test
    ds.sa.targets=1;

    % each participant is independent
    ds.sa.chunks=n;

    group_cell{n}=ds;
end

group_ds=cosmo_stack(group_cell);

% define the clustering neighborhood. By default, features next to each
% other in time (train time or test time) are considered neighbors.
nbrhood=cosmo_cluster_neighborhood(group_ds);

%..........................................................................
opt=struct();
opt.cluster_stat='tfce';
opt.niter=mvpa.cosmo.group.niter;
opt.h0_mean=mvpa.cosmo.group.h0_mean;
%..........................................................................

% run multiple comparison correction
ds_result=cosmo_montecarlo_cluster_stat(group_ds, nbrhood, opt);

% for easier unflattenening, the time dimensions are moved back
% from feature to sample dimensions
ds_result_time = cosmo_dim_transpose(ds_result,{'freq','time'});

% flatten into a 2D array
[result.z_score,result.dim_labels,result.dim_values]=cosmo_unflatten(ds_result_time,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-ANALYSIS MODE - Sensor-level Time Generalization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode,'glm-sensor-time-generalization-tfce')

%-IMPORT FILES AND INFORMATION
%==========================================================================
nsubjs = length(meg.subject.list);

%-ANALYSIS
%==========================================================================
group_cell=cell(nsubjs,1);
for n=1:nsubjs
    % convert data to data structure
    %https://groups.google.com/g/cosmomvpa/c/GWepADY2M_g
    dgm_ds=cosmo_flatten(data{n},{'train_time','test_time'},{1:meg.epoch.ntimes,1:meg.epoch.ntimes},1);

    % make train_time and test_time a feature dimension
    ds=cosmo_dim_transpose(dgm_ds,{'train_time','test_time'},2);

    % for one-sample t-test
    ds.sa.targets=1;

    % each participant is independent
    ds.sa.chunks=n;

    group_cell{n}=ds;
end

group_ds=cosmo_stack(group_cell);

% define the clustering neighborhood. By default, features next to each
% other in time (train time or test time) are considered neighbors.
nbrhood=cosmo_cluster_neighborhood(group_ds);

%..........................................................................
opt=struct();
opt.cluster_stat='tfce';
opt.niter=mvpa.cosmo.group.niter;
opt.h0_mean=mvpa.cosmo.group.h0_mean;
%..........................................................................

% run multiple comparison correction
ds_result=cosmo_montecarlo_cluster_stat(group_ds,nbrhood,opt);

% for easier unflattenening, the time dimensions are moved back
% from feature to sample dimensions
ds_result_time = cosmo_dim_transpose(ds_result,{'train_time','test_time'});

% flatten into a 2D array
[result.z_score,result.dim_labels,result.dim_values]=cosmo_unflatten(ds_result_time,1);

end

end
