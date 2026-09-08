%__________________________________________________________________________
function dnn_cosmo_glm(mypath,mri,meg,mvpa,type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Model @ DNN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(type,strcat('model@',mvpa.dnn.architecture,'-',mvpa.dnn.database))

%-IMPORT FILES AND INFORMATION
%==========================================================================
%-Target and Chunk
targets = mvpa.cosmo.glm.targets;

%-Load Model RDMs
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

%-Load DNN RDMs
str.file = eval(mypath.file.data.rdm.dnn);
str.folder = eval(mypath.folder.data.rdm);
DNN = load(strcat(str.folder,'\',str.file));

nlayers = length(mvpa.dnn.layer.name);

%-ANALYSIS
%==========================================================================
for l=1:nlayers
   corr(l) = MeasureRDMCorr(DNN.RDM.(matlab.lang.makeValidName(mvpa.dnn.layer.name{l})),...
       'mode','mri','target_dsm',pred_cell(2,:),'analysis',mvpa.cosmo.glm.analysis.dnn,'input_dsm',true); 
end

%-PLOT THE RESULTS
%==========================================================================
%-Create X and Y
part = length(mvpa.dnn.layer.name);
x1 = categorical(mvpa.dnn.layer.name);
x1 = reordercats(x1,mvpa.dnn.layer.name);
x2 = categorical();
for p = 1:length(pred_cell(1,:))
    for l = 1:nlayers
        y(p,l) = corr(l).samples(p);
    end
end

vis_args = struct();
vis_args.part = part;
vis_args.type = 'mri-subject-glm';
vis_args.legend_name = pred_cell(1,:);
vis_args.label_name.title = mvpa.vis.dnn.glm.title;
vis_args.label_name.xlabel = mvpa.vis.dnn.glm.xlabel;
vis_args.label_name.ylabel = mvpa.vis.dnn.glm.ylabel;
vis_args.color = cell2mat(values(mvpa.vis.color.glm,vis_args.legend_name)');
DrawScatterSplit(x1,x2,y,vis_args);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-DNN @ MRI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(type,strcat(mvpa.dnn.architecture,'-',mvpa.dnn.database,'@mri'))

%-IMPORT FILES AND INFORMATION
%==========================================================================
%-Target and Chunk
targets = mvpa.cosmo.glm.targets;

%-Load DNN RDMs
str.file = eval(mypath.file.data.rdm.dnn);
str.folder = eval(mypath.folder.data.rdm);
load(strcat(str.folder,'\',str.file));
pred_cell = cell(2,length(mvpa.dnn.layer.name));
for l=1:length(mvpa.dnn.layer.name)
   pred_cell{1,l} = mvpa.dnn.layer.name{l};
   pred_cell{2,l} = eval(sprintf('RDM.%s',matlab.lang.makeValidName(mvpa.dnn.layer.name{l})));
end

nsubjs = length(mri.subject.list);
nmasks = length(mvpa.cosmo.glm.mask.name);

%-Create output structure
table_size = [nsubjs nmasks];
table_type = repmat({'cell'},1,nmasks);
glm_table = table('Size',table_size,...
    'VariableTypes',table_type,...
    'VariableNames',mvpa.cosmo.glm.mask.name,...
    'RowNames',mri.subject.list);

%-LOAD DATA
%==========================================================================
data = load_mri_rdm(mypath,mvpa,mri);

%-ANALYSIS
%==========================================================================
for m = 1:nmasks
    for s = 1:nsubjs
        subj = mri.subject.list{s};

        %-RDM
        %..................................................................
        ds = data.(subj).ds.(mvpa.cosmo.glm.mask.name{m});
        %..................................................................
        
        %-ANALYSIS
        %------------------------------------------------------------------
        %-Regression
        measure = @MeasureRDMCorr;
        
        %-Measure arguments
        %..................................................................
        measure_args = struct();
        measure_args.mode = 'mri';
        measure_args.analysis = mvpa.cosmo.glm.analysis.mri;
        measure_args.input_dsm = true;
        measure_args.metric_dsm = mvpa.cosmo.glm.metric.mri;
        measure_args.type = mvpa.cosmo.glm.type.mri;
        measure_args.frrsa = mvpa.cosmo.glm.frrsa;
        measure_args.center_data = mvpa.cosmo.glm.voxscaling;
        measure_args.target_dsm = pred_cell(2,:);
        measure_args.labels = pred_cell(1,:)';
        %..................................................................
        
        %-Apply the measure
        ds_glm = measure(ds,measure_args);
        
        %-Beta values of each subject
        glm_table(s,m) = {num2str(ds_glm.samples)};
    end
end

%-PLOT THE RESULTS
%==========================================================================
for s = 1:nsubjs
    %-Create X and Y
    part = mvpa.cosmo.glm.mask.partition;
    x1 = categorical(mvpa.cosmo.glm.mask.name(1:part));
    x1 = reordercats(x1,mvpa.cosmo.glm.mask.name(1:part));
    x2 = categorical(mvpa.cosmo.glm.mask.name(part+1:end));
    x2 = reordercats(x2,mvpa.cosmo.glm.mask.name(part+1:end));
    for n = 1:length(ds_glm.sa.labels)
        for m = 1:length(mvpa.cosmo.glm.mask.name)
            sm_pred = str2num(glm_table{s,m}{:});
            y(n,m) = sm_pred(n);
        end
    end

    vis_args = struct();
    vis_args.part = part;
    vis_args.type = 'mri-subject-glm';
    vis_args.legend_name = pred_cell(1,:);
    vis_args.label_name.title = mvpa.vis.mri.glm.title;
    vis_args.label_name.xlabel = mvpa.vis.mri.glm.xlabel;
    vis_args.label_name.ylabel = mvpa.vis.mri.glm.ylabel;
    vis_args.color = cell2mat(values(mvpa.vis.color.layer,vis_args.legend_name)');
    DrawScatterSplit(x1,x2,y,vis_args);
    
    %-Save figure
    str.file = eval(mypath.file.figure.dnn.glm.mri.subj);
    str.folder = eval(mypath.folder.result.dnn.figure);
    fn.figure.scatter = strcat(str.folder,'/',str.file);
    saveas(gcf,fn.figure.scatter);
    clf;
end

%-SAVE THE RESULTS
%==========================================================================
%-Save each subject result in a seperate mat file
for s = 1:nsubjs
    glm_table_subj = glm_table(s,:);
    str.file = eval(mypath.file.subj.dnn.glm.roi);
    str.folder = eval(mypath.folder.result.dnn.glm);
    fn.result.glm = strcat(str.folder,'/',str.file);
    save(fn.result.glm,'glm_table_subj','pred_cell');
end

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Model & DNN @ MRI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexpcmp(type,strcat('model-',mvpa.dnn.architecture,'-',mvpa.dnn.database,'-.*@mri'))

%-IMPORT FILES AND INFORMATION
%==========================================================================
%-Target and Chunk
targets = mvpa.cosmo.glm.targets;

%-Load Model RDMs
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

%-Load DNN RDMs
str.file = eval(mypath.file.data.rdm.dnn);
str.folder = eval(mypath.folder.data.rdm);
load(strcat(str.folder,'\',str.file));
layers = mvpa.dnn.model.name(char(extractBetween(type,strcat('model-',mvpa.dnn.architecture,'-',mvpa.dnn.database,'-'),'@')));
for l=1:numel(layers)
        pred_cell{1,end+1} = layers{l};
        pred_cell{2,end} = eval(sprintf('RDM.%s',matlab.lang.makeValidName(layers{l})));
end

nsubjs = length(mri.subject.list);
nmasks = length(mvpa.cosmo.glm.mask.name);

%-Create output structure
table_size = [nsubjs nmasks];
table_type = repmat({'cell'},1,nmasks);
glm_table = table('Size',table_size,...
    'VariableTypes',table_type,...
    'VariableNames',mvpa.cosmo.glm.mask.name,...
    'RowNames',mri.subject.list);

%-LOAD DATA
%==========================================================================
data = load_mri_rdm(mypath,mvpa,mri);

%-ANALYSIS
%==========================================================================
for m = 1:nmasks
    for s = 1:nsubjs
        subj = mri.subject.list{s};

        %-RDM
        %..................................................................
        ds = data.(subj).ds.(mvpa.cosmo.glm.mask.name{m});
        %..................................................................
        
        %-ANALYSIS
        %------------------------------------------------------------------
        %-Regression
        measure = @MeasureRDMCorr;
        
        %-Measure arguments
        %..................................................................
        measure_args = struct();
        measure_args.mode = 'mri';
        measure_args.analysis = mvpa.cosmo.glm.analysis.mri;
        measure_args.input_dsm = true;
        measure_args.metric_dsm = mvpa.cosmo.glm.metric.mri;
        measure_args.type = mvpa.cosmo.glm.type.mri;
        measure_args.frrsa = mvpa.cosmo.glm.frrsa;
        measure_args.center_data = mvpa.cosmo.glm.voxscaling;
        measure_args.target_dsm = pred_cell(2,:);
        measure_args.labels = pred_cell(1,:)';
        %..................................................................
        
        %-Apply the measure
        ds_glm = measure(ds,measure_args);
        
        %-Beta values of each subject
        glm_table(s,m) = {num2str(ds_glm.samples)};
    end
end

%-PLOT THE RESULTS
%==========================================================================
for s = 1:nsubjs
    %-Create X and Y
    part = mvpa.cosmo.glm.mask.partition;
    x1 = categorical(mvpa.cosmo.glm.mask.name(1:part));
    x1 = reordercats(x1,mvpa.cosmo.glm.mask.name(1:part));
    x2 = categorical(mvpa.cosmo.glm.mask.name(part+1:end));
    x2 = reordercats(x2,mvpa.cosmo.glm.mask.name(part+1:end));
    for n = 1:length(ds_glm.sa.labels)
        for m = 1:length(mvpa.cosmo.glm.mask.name)
            sm_pred = str2num(glm_table{s,m}{:});
            y(n,m) = sm_pred(n);
        end
    end

    vis_args = struct();
    vis_args.part = part;
    vis_args.type = 'mri-subject-glm';
    vis_args.legend_name = pred_cell(1,:);
    vis_args.label_name.title = mvpa.vis.mri.glm.title;
    vis_args.label_name.xlabel = mvpa.vis.mri.glm.xlabel;
    vis_args.label_name.ylabel = mvpa.vis.mri.glm.ylabel;
    vis_args.color = cell2mat(values(mvpa.vis.color.model,vis_args.legend_name)');
    DrawScatterSplit(x1,x2,y,vis_args);
    
    %-Save figure
    str.file = eval(mypath.file.figure.dnn.glm.mri.subj);
    str.folder = eval(mypath.folder.result.dnn.figure);
    fn.figure.scatter = strcat(str.folder,'/',str.file);
    saveas(gcf,fn.figure.scatter);
    clf;
end

%-SAVE THE RESULTS
%==========================================================================
%-Save each subject result in a seperate mat file
for s = 1:nsubjs
    glm_table_subj = glm_table(s,:);
    str.file = eval(mypath.file.subj.dnn.glm.roi);
    str.folder = eval(mypath.folder.result.dnn.glm);
    fn.result.glm = strcat(str.folder,'/',str.file);
    save(fn.result.glm,'glm_table_subj','pred_cell');
end

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-DNN @ MEG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(type,strcat(mvpa.dnn.architecture,'-',mvpa.dnn.database,'@meg'))

%-IMPORT FILES AND INFORMATION
%==========================================================================
%-Target and Chunk
targets = mvpa.cosmo.glm.targets;

%-Load DNN RDMs
str.file = eval(mypath.file.data.rdm.dnn);
str.folder = eval(mypath.folder.data.rdm);
load(strcat(str.folder,'\',str.file));
pred_cell = cell(2,length(mvpa.dnn.layer.name));
for l=1:length(mvpa.dnn.layer.name)
   pred_cell{1,l} = mvpa.dnn.layer.name{l};
   pred_cell{2,l} = eval(sprintf('RDM.%s',matlab.lang.makeValidName(mvpa.dnn.layer.name{l})));
end

nsubjs = length(meg.subject.list);
ntimes = meg.epoch.ntimes;
sensor = 'grad';

%-LOAD DATA
%==========================================================================
data = load_meg_rdm(mypath,mvpa,meg,sensor);

%-ANALYSIS
%==========================================================================
for s = 1:nsubjs
    subj = meg.subject.list{s};
    
    %-LOAD SUBJECT DATA
    %----------------------------------------------------------------------
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
    vis_args.color = cell2mat(values(mvpa.vis.color.layer,vis_args.legend_name)');
    DrawScatter(x,y,vis_args);
    
    str.file = eval(mypath.file.figure.meg.glm.subj);
    str.folder = eval(mypath.folder.result.dnn.figure);
    fn.figure.scatter = strcat(str.folder,'/',str.file);
    saveas(gcf,fn.figure.scatter);
    clf;
    
    %-SAVE THE RESULTS
    %----------------------------------------------------------------------
    %-Save each subject result in a seperate mat file
    glm_table_subj = array2table(y);
    glm_table_subj.Properties.VariableNames = compose('time_%d', 1:length(y));
    glm_table_subj.Properties.RowNames = ds_glm(1).sa.labels;
    str.file = eval(mypath.file.subj.dnn.glm.sensor);
    str.folder = eval(mypath.folder.result.dnn.glm);
    fn.result.glm = strcat(str.folder,'/',str.file);
    save(fn.result.glm,'glm_table_subj','pred_cell');
    
end

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Model & DNN @ MEG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexpcmp(type,strcat('model-',mvpa.dnn.architecture,'-',mvpa.dnn.database,'-.*@meg'))

%-IMPORT FILES AND INFORMATION
%==========================================================================
%-Target and Chunk
targets = mvpa.cosmo.glm.targets;

%-Load Model RDMs
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

%-Load DNN RDMs
str.file = eval(mypath.file.data.rdm.dnn);
str.folder = eval(mypath.folder.data.rdm);
load(strcat(str.folder,'\',str.file));
layers = mvpa.dnn.model.name(char(extractBetween(type,strcat('model-',mvpa.dnn.architecture,'-',mvpa.dnn.database,'-'),'@')));
for l=1:numel(layers)
        pred_cell{1,end+1} = layers{l};
        pred_cell{2,end} = eval(sprintf('RDM.%s',matlab.lang.makeValidName(layers{l})));
end

nsubjs = length(meg.subject.list);
ntimes = meg.epoch.ntimes;
sensor = 'grad';

%-LOAD DATA
%==========================================================================
data = load_meg_rdm(mypath,mvpa,meg,sensor);

%-ANALYSIS
%==========================================================================
for s = 1:nsubjs
    subj = meg.subject.list{s};
    
    %-LOAD SUBJECT DATA
    %----------------------------------------------------------------------
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
    vis_args.color = cell2mat(values(mvpa.vis.color.model,vis_args.legend_name)');
    DrawScatter(x,y,vis_args);
    
    str.file = eval(mypath.file.figure.meg.glm.subj);
    str.folder = eval(mypath.folder.result.dnn.figure);
    fn.figure.scatter = strcat(str.folder,'/',str.file);
    saveas(gcf,fn.figure.scatter);
    clf;
    
    %-SAVE THE RESULTS
    %----------------------------------------------------------------------
    %-Save each subject result in a seperate mat file
    glm_table_subj = array2table(y);
    glm_table_subj.Properties.VariableNames = compose('time_%d', 1:length(y));
    glm_table_subj.Properties.RowNames = ds_glm(1).sa.labels;
    str.file = eval(mypath.file.subj.dnn.glm.sensor);
    str.folder = eval(mypath.folder.result.dnn.glm);
    fn.result.glm = strcat(str.folder,'/',str.file);
    save(fn.result.glm,'glm_table_subj','pred_cell');
    
end

close all;

end

end

%__________________________________________________________________________
function data = load_mri_rdm(mypath,mvpa,mri)

nsubjs = length(mri.subject.list);
nmasks = length(mvpa.cosmo.glm.mask.name);

for s = 1:nsubjs
    subj = mri.subject.list{s};

    str.file = eval(mypath.file.subj.mri.rdm);
    str.folder = eval(mypath.folder.result.mri.rdm);
    fn.in = strcat(str.folder,'/',str.file);
    load(fn.in);
    data.(subj).ds = RDM;
end

end

%__________________________________________________________________________
function data = load_meg_rdm(mypath,mvpa,meg,sensor)

nsubjs = length(meg.subject.list);

    for s = 1:nsubjs
        subj = meg.subject.list{s};
        
        str.file = eval(mypath.file.subj.meg.rdm);
        str.folder = eval(mypath.folder.result.meg.rdm);
        fn.in = strcat(str.folder,'/',str.file);
        load(fn.in);
        data.(subj).ds = RDM.(sensor);
    end

end
