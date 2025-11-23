%__________________________________________________________________________
function meg_cosmo_glm_sensor_searchlight(mypath,meg,mvpa,sensor)

%-IMPORT FILES AND INFORMATION
%==========================================================================
%-Load predictor RDMs
str.file = eval(mypath.file.data.rdm.predictor);
str.folder = eval(mypath.folder.data.rdm);
load(strcat(str.folder,'\',str.file));
pred_cell = cell(2,length(mvpa.rdm.predictor.name));
for i=1:length(mvpa.rdm.predictor.name)
   pred_cell{1,i} = mvpa.rdm.predictor.name{i};
   pred_cell{2,i} = eval(sprintf('RDM.%s.%s',...
       mvpa.rdm.predictor.type,...
       mvpa.rdm.predictor.name{i}));
end

nsubjs = length(meg.subject.list);
npreds = length(mvpa.rdm.predictor.name);

for s = 1:nsubjs
    subj = meg.subject.list{s};

    %-LOAD DATA
    %======================================================================
    str.file = eval(mypath.file.subj.meg.preprocessed);
    str.folder = eval(mypath.folder.subj.meg);
    fn.in = strcat(str.folder,'/',str.file);
    ds = load(fn.in);
    ds = convert_cosmo_ds(ds);
    
    %-Load data of Gradiometer or Magnetometer based on sensor parameter.
    if strcmp(sensor,'grad')
        chanx='meg_axial';
    elseif strcmp(sensor,'mag')
        chanx='meg_planar';
    end
    
    layout=cosmo_meeg_find_layout(ds,'chantype',chanx);
    layout_label = layout.label(1:end-2,:);
    feature_msk = cosmo_dim_match(ds,'chan',layout_label);
    ds = cosmo_slice(ds,feature_msk,2);
    ds = cosmo_fx(ds,@(x)mean(x,1),'targets',1);
    
    %-ANALYSIS
    %======================================================================
    %-Regression
    measure = @MeasureRDMCorr;

    %-Measure arguments
    %......................................................................
    measure_args = struct();
    measure_args.mode = 'mri';
    measure_args.analysis = mvpa.cosmo.glm.analysis.meg;
    measure_args.metric_dsm = mvpa.cosmo.glm.metric.meg;
    measure_args.type = mvpa.cosmo.glm.type.meg;
    measure_args.frrsa = mvpa.cosmo.glm.frrsa;
    measure_args.center_data = mvpa.cosmo.glm.chanscaling;
    measure_args.pseudo = mvpa.cosmo.glm.npseudo;
    measure_args.perm = mvpa.cosmo.glm.nperm;
    measure_args.cv = mvpa.cosmo.glm.cv;
    measure_args.mnn = mvpa.cosmo.glm.mnn;
    measure_args.smooth = mvpa.cosmo.glm.smooth;
    measure_args.vp = {};
    measure_args.target_dsm = pred_cell(2,:);
    measure_args.labels = pred_cell(1,:)';
    %......................................................................

    chan_nbrhood = cosmo_meeg_chan_neighborhood(ds,'count',mvpa.cosmo.glm.searchlight.nchan,'chantype',chanx);
    time_nbrhood = cosmo_interval_neighborhood(ds,'time','radius',0);
    nbrhood = cosmo_cross_neighborhood(ds,{chan_nbrhood,time_nbrhood});
    
    %-Run searchlight
    ds_sl = cosmo_searchlight(ds,nbrhood,measure,measure_args);
    
    %-SAVE THE RESULTS
    %======================================================================
    for p=1:npreds
        %-Map to FieldTrip struct for visualization
        tf_map(p) = cosmo_map2meeg(cosmo_slice(ds_sl,p));
    end

    str.file = eval(mypath.file.subj.meg.glm.searchlight);
    str.folder = eval(mypath.folder.result.meg.searchlight);
    fn.out = strcat(str.folder,'/',str.file);
    save(fn.out,'tf_map');

end

end

%__________________________________________________________________________
function ds_cosmo = convert_cosmo_ds(ds)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%https://www.cosmomvpa.org/faq.html#faq-get-ecog-data-in-cosmomvpa-struct
%-This function is input structure dependent!
%-It means that 'trial', 'channel_name', 'time', and 'index' depends on
%naming convention of 'ds'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ds_cosmo = cosmo_flatten(ds.trial,{'chan','time'},{cellstr(ds.channel_name),ds.time});
    ds_cosmo.a.meeg.samples_field = 'trial';
    ds_cosmo.sa.targets = ds.index';
end
