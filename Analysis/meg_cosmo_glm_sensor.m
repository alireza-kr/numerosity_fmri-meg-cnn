%__________________________________________________________________________
function meg_cosmo_glm_sensor(mypath,meg,mvpa,sensor,mode)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-ANALYSIS MODE - Sensor-level GLM RSA (Time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode,'time')
    
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
data = load_data(mypath,mvpa,meg,sensor,'rdm-time');

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
    measure_args.input_dsm = true;
    measure_args.metric_dsm = mvpa.cosmo.glm.metric.meg;
    measure_args.type = mvpa.cosmo.glm.type.meg;
    measure_args.frrsa = mvpa.cosmo.glm.frrsa;
    measure_args.center_data = mvpa.cosmo.glm.chanscaling;
    measure_args.pseudo = mvpa.cosmo.glm.npseudo;
    measure_args.perm = mvpa.cosmo.glm.nperm;
    measure_args.cv = mvpa.cosmo.glm.cv;
    measure_args.mnn = mvpa.cosmo.glm.mnn;
    measure_args.smooth = mvpa.cosmo.glm.smooth;
    measure_args.target_dsm = pred_cell(2,:);
    measure_args.labels = pred_cell(1,:)';
    %......................................................................
    
    %-Apply the measure
    ds_glm = measure(ds,measure_args);

    %-PLOT THE RESULTS
    %----------------------------------------------------------------------
    x = linspace(1000*meg.epoch.tmin,1000*meg.epoch.tmax,ntimes);
    for n = 1:length(ds_glm(1).sa.labels)
        for m = 1:size(ds_glm,2)
            y(n,m) = ds_glm(m).samples(n);
        end
    end
    
    vis_args = struct();
    vis_args.type = 'meg-subject-glm';
    vis_args.legend_name = pred_cell(1,:);
    vis_args.label_name.title = mvpa.vis.meg.glm.title;
    vis_args.label_name.xlabel = mvpa.vis.meg.glm.xlabel;
    vis_args.label_name.ylabel = mvpa.vis.meg.glm.ylabel;
    vis_args.color = cell2mat(values(mvpa.vis.color.glm,vis_args.legend_name)');
    DrawScatter(x,y,vis_args);
    
    str.file = eval(mypath.file.figure.meg.glm.subj);
    str.folder = eval(mypath.folder.result.meg.figure);
    fn.figure.scatter = strcat(str.folder,'/',str.file);
    saveas(gcf,fn.figure.scatter);
    clf;
    
    %-SAVE THE RESULTS
    %----------------------------------------------------------------------
    %-Save each subject result in a seperate mat file
    glm_table_subj = array2table(y);
    glm_table_subj.Properties.VariableNames = compose('time_%d', 1:length(y));
    glm_table_subj.Properties.RowNames = ds_glm(1).sa.labels;
    str.file = eval(mypath.file.subj.meg.glm.sensor);
    str.folder = eval(mypath.folder.result.meg.glm);
    fn.result.glm = strcat(str.folder,'/',str.file);
    save(fn.result.glm,'glm_table_subj','pred_cell');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-ANALYSIS MODE - Sensor-level GLM RSA (Time-Frequency)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode,'time-frequency')

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
nfreq = meg.convo.frex_num;

%-LOAD DATA
%==========================================================================
data = load_data(mypath,mvpa,meg,sensor,'rdm-time-frequency');

%-ANALYSIS
%==========================================================================
for s = 1:nsubjs
    subj = meg.subject.list{s};
    
    for f = 1:nfreq
    
        %-LOAD SUBJECT DATA
        %------------------------------------------------------------------
        %-RDM
        %..................................................................
        ds = data.(subj).ds{1,f};
        %..................................................................

        %-ANALYSIS
        %------------------------------------------------------------------
        measure = @MeasureRDMCorr;

        %-Measure arguments
        %..................................................................
        measure_args = struct();
        measure_args.mode = 'meg';
        measure_args.analysis = mvpa.cosmo.glm.analysis.meg;
        measure_args.input_dsm = true;
        measure_args.metric_dsm = mvpa.cosmo.glm.metric.meg;
        measure_args.type = mvpa.cosmo.glm.type.meg;
        measure_args.frrsa = mvpa.cosmo.glm.frrsa;
        measure_args.center_data = mvpa.cosmo.glm.chanscaling;
        measure_args.pseudo = mvpa.cosmo.glm.npseudo;
        measure_args.perm = mvpa.cosmo.glm.nperm;
        measure_args.cv = mvpa.cosmo.glm.cv;
        measure_args.mnn = mvpa.cosmo.glm.mnn;
        measure_args.smooth = mvpa.cosmo.glm.smooth;
        measure_args.target_dsm = pred_cell(2,:);
        measure_args.labels = pred_cell(1,:)';
        %..................................................................

        %-Apply the measure
        ds_glm{f} = measure(ds,measure_args);
    
    end

    %-PLOT THE RESULTS
    %----------------------------------------------------------------------
    glm_cell_subj = {};
    
    for p=1:length(ds_glm{f}(1).sa.labels)
        for f=1:nfreq
            for t = 1:size(ds_glm{f},2)
                mat(f,t) = ds_glm{f}(t).samples(p);
            end
        end
        
        glm_cell_subj{p} = mat;
        
        vis_args = struct();
        vis_args.xlabel = 'Time';
        vis_args.ylabel = 'Frequency';
        vis_args.xtick = linspace(1,ntimes,10);
        vis_args.ytick = linspace(1,nfreq,nfreq);
        vis_args.xticklabel = linspace(1000*meg.epoch.tmin,1000*meg.epoch.tmax,10);
        vis_args.yticklabel = linspace(meg.convo.frex_min,meg.convo.frex_max,nfreq);
        vis_args.range = [min(min(mat)) max(max(mat))];
        vis_args.title = ds_glm{f}(1).sa.labels{p};
        DrawMatrix(mat,vis_args);
        mat_plot{p} = gcf;
        figure;
    end

    for p=1:length(ds_glm{f}(1).sa.labels)
        mat_subplot{p} = subplot(1,length(ds_glm{f}(1).sa.labels),p);
        title(ds_glm{f}(1).sa.labels{p}) 
        copyobj(allchild(get(mat_plot{p},'CurrentAxes')),mat_subplot{p});
    end
    set(gcf,'Position',[100 100 1700 250]);

    str.file = eval(mypath.file.figure.meg.glm.tfreq.subj);
    str.folder = eval(mypath.folder.result.meg.figure);
    fn.figure.matrix = strcat(str.folder,'/',str.file);
    saveas(gcf,fn.figure.matrix);
    close all;

    %-SAVE THE RESULTS
    %----------------------------------------------------------------------
    %-Save each subject result in a seperate mat file
    str.file = eval(mypath.file.subj.meg.glm.tfreq.sensor);
    str.folder = eval(mypath.folder.result.meg.glm);
    fn.result.glm = strcat(str.folder,'/',str.file);
    save(fn.result.glm,'glm_cell_subj','pred_cell');
end

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
elseif strcmp(type,'rdm-time')
    for s = 1:nsubjs
        subj = meg.subject.list{s};
        
        str.file = eval(mypath.file.subj.meg.rdm);
        str.folder = eval(mypath.folder.result.meg.rdm);
        fn.in = strcat(str.folder,'/',str.file);
        load(fn.in);
        data.(subj).ds = RDM.(sensor);
    end
elseif strcmp(type,'rdm-time-frequency')
    for s = 1:nsubjs
        subj = meg.subject.list{s};
        
        str.file = eval(mypath.file.subj.meg.freq.rdm);
        str.folder = eval(mypath.folder.result.meg.rdm);
        fn.in = strcat(str.folder,'/',str.file);
        load(fn.in);
        data.(subj).ds = RDM.(sensor);
    end
end

end
