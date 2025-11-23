%__________________________________________________________________________
function meg_cosmo_vp_sensor(mypath,meg,mvpa,sensor)

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
ntimes = meg.epoch.ntimes;

%-LOAD DATA
%==========================================================================
data = load_data(mypath,mvpa,meg,sensor,'rdm');

%-ANALYSIS
%==========================================================================
for s = 1:nsubjs
    subj = meg.subject.list{s};
    
    %-LOAD SUBJECT DATA
    %----------------------------------------------------------------------
    %-Raw
    %......................................................................
    %ds = data.(subj).ds;
    %......................................................................

    %-RDM
    %......................................................................
    ds = data.(subj).ds;
    %......................................................................
    
    %-ANALYSIS
    %----------------------------------------------------------------------
    measure = @MeasureRDMCorr;

    %-Measure arguments
    %......................................................................
    measure_args = struct();
    measure_args.mode = 'meg';
    measure_args.analysis = mvpa.cosmo.glm.analysis.meg;
    measure_args.metric = mvpa.cosmo.glm.metric.meg;
    measure_args.type = mvpa.cosmo.glm.type.meg;
    measure_args.center_data = mvpa.cosmo.glm.chanscaling;
    measure_args.pseudo = mvpa.cosmo.glm.npseudo;
    measure_args.perm = mvpa.cosmo.glm.nperm;
    measure_args.cv = mvpa.cosmo.glm.cv;
    measure_args.mnn = mvpa.cosmo.glm.mnn;
    measure_args.smooth = mvpa.cosmo.glm.smooth;
    measure_args.vp = mvpa.vp;
    measure_args.target_dsm = pred_cell(2,:);
    measure_args.labels = pred_cell(1,:)';
    %......................................................................
    
    %-Apply the measure
    ds_vp = measure(ds,measure_args);

    %-PLOT THE RESULTS
    %----------------------------------------------------------------------
    x = linspace(1000*meg.epoch.tmin,1000*meg.epoch.tmax,ntimes);
    for n = 1:mvpa.vp.nmodels
        for m = 1:size(ds_vp,2)
            y(n,m) = ds_vp(m).samples(n);
        end
    end
    
    vis_args = struct();
    vis_args.type = 'meg-subject-vp';
    vis_args.legend_name = mvpa.vp.visualize;
    vis_args.label_name.title = mvpa.vis.meg.glm.title;
    vis_args.label_name.xlabel = mvpa.vis.meg.glm.xlabel;
    vis_args.label_name.ylabel = mvpa.vis.meg.glm.ylabel;
    vis_args.color = cell2mat(values(mvpa.vis.color.glm,vis_args.legend_name)');
    DrawScatter(x,y(find(ismember(ds_vp(1).sa.labels,mvpa.vp.visualize)),:),vis_args);
    
    str.file = eval(mypath.file.figure.meg.vp.subj);
    str.folder = eval(mypath.folder.result.meg.figure);
    fn.figure.scatter = strcat(str.folder,'/',str.file);
    saveas(gcf,fn.figure.scatter);
    clf;
    
    %-SAVE THE RESULTS
    %----------------------------------------------------------------------
    %-Save each subject result in a seperate mat file
    vp_table_subj = array2table(y);
    vp_table_subj.Properties.VariableNames = compose('time_%d', 1:length(y));
    vp_table_subj.Properties.RowNames = ds_vp(1).sa.labels;
    str.file = eval(mypath.file.subj.meg.vp.sensor);
    str.folder = eval(mypath.folder.result.meg.vp);
    fn.result.vp = strcat(str.folder,'/',str.file);
    save(fn.result.vp,'vp_table_subj','pred_cell');
    
end

end

%__________________________________________________________________________
function data = load_data(mypath,mvpa,meg,sensor,type)

nsubjs = length(meg.subject.list);

if strcmp(type,'raw')
    for s = 1:nsubjs
        subj = meg.subject.list{s};
        
        str.file = eval(mypath.file.subj.meg.preprocessed);
        str.folder = eval(mypath.folder.subj.meg);
        fn.in = strcat(str.folder,'/',str.file);
        data.(subj).ds = load(fn.in);

        %-Load data of Gradiometer or Magnetometer based on sensor parameter.
        if strcmp(sensor,'mag')
            data.(subj).ds.trial = data.(subj).ds.trial(:,data.(subj).ds.channel_type(:,1)=='m',:);
        elseif strcmp(sensor,'grad')
            data.(subj).ds.trial = data.(subj).ds.trial(:,data.(subj).ds.channel_type(:,1)=='g',:);
        else
            warning("Please choose sensor type (mag or grad)!");
        end
    end
elseif strcmp(type,'rdm')
    for s = 1:nsubjs
        subj = meg.subject.list{s};
        
        str.file = eval(mypath.file.subj.meg.rdm);
        str.folder = eval(mypath.folder.result.meg.rdm);
        fn.in = strcat(str.folder,'/',str.file);
        load(fn.in);
        data.(subj).ds = RDM.(sensor);
    end
end

end
