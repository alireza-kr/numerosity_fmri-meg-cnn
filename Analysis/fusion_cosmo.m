function fusion_cosmo(mypath,mri,meg,mvpa,sensor,mode)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-ANALYSIS MODE - INDIVIDUAL MEG & GROUP-AVERAGE MRI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode,'meg')

%-IMPORT FILES AND INFORMATION
%==========================================================================
%-Target and Chunk
targets = mvpa.cosmo.glm.targets;

%-Load Predictor RDMs (Model RDM(s) and Group-Average MRI RDM)
str.file.model = eval(mypath.file.data.rdm.predictor);
str.file.mri = eval(mypath.file.data.rdm.avg.mri);
str.folder = eval(mypath.folder.data.rdm);

%-Load model RDM(s) in 'preds'
load(strcat(str.folder,'\',str.file.model));
for i=1:length(mvpa.rdm.predictor.name)-1
   preds{i} = ...
       eval(sprintf('RDM.%s.%s',...
       mvpa.rdm.predictor.type,...
       mvpa.rdm.predictor.name{i}));
end

%-Load Group-Average MRI RDM
load(strcat(str.folder,'\',str.file.mri),'RDM');

nmasks = length(mvpa.cosmo.glm.mask.name);
ntimes = meg.epoch.ntimes;
nsubjs = length(meg.subject.list);

%-LOAD DATA
%==========================================================================
data = load_data_meg(mypath,mvpa,meg,sensor,'rdm');

%-ANALYSIS
%==========================================================================
for m = 1:nmasks
    for s = 1:nsubjs
        subj = meg.subject.list{s};

        %-LOAD SUBJECT DATA
        %------------------------------------------------------------------
        %-Raw
        %..................................................................
        %ds = data.(subj).ds;
        %..................................................................
        
        %-RDM
        %..................................................................
        ds = data.(subj).ds;
        %..................................................................

        %-ANALYSIS
        %------------------------------------------------------------------
        %-Add Group-Average MRI RDM of each ROI to 'preds'
        preds{length(mvpa.rdm.predictor.name)} = RDM.(mvpa.cosmo.glm.mask.name{m});

        %-Regression
        measure = @MeasureRDMCorr;

        %-Measure arguments
        %..................................................................
        measure_args = struct();
        measure_args.mode = 'meg';
        measure_args.analysis = mvpa.cosmo.glm.analysis.meg;
        measure_args.input_dsm = true;
        measure_args.metric = mvpa.cosmo.glm.metric.meg;
        measure_args.type = mvpa.cosmo.glm.type.meg;
        measure_args.frrsa = mvpa.cosmo.glm.frrsa;
        measure_args.frrsa_dsm = preds{end};   % MRI at region m
        measure_args.center_data = mvpa.cosmo.glm.chanscaling;
        measure_args.pseudo = mvpa.cosmo.glm.npseudo;
        measure_args.perm = mvpa.cosmo.glm.nperm;
        measure_args.cv = mvpa.cosmo.glm.cv;
        measure_args.mnn = mvpa.cosmo.glm.mnn;
        measure_args.smooth = mvpa.cosmo.glm.smooth;
        measure_args.vp = mvpa.vp;
        measure_args.target_dsm = preds;
        measure_args.labels = mvpa.rdm.predictor.name';
        %..................................................................

        %-Apply the measure
        ds_fusion = measure(ds,measure_args);
        
        %-R2 values of each subject
        for p = 1:length(ds_fusion(1).sa.labels)
            for t = 1:size(ds_fusion,2)
                %data(MASK,PREDICTORS,SUBJECTS,TIME)
                result(m,p,s,t) = ds_fusion(t).samples(p);
            end
        end
    end
end

%-PLOT THE RESULTS
%==========================================================================
x = linspace(1000*meg.epoch.tmin,1000*meg.epoch.tmax,ntimes);
for s = 1:nsubjs
    figure('WindowState','maximized'); pause(0.1);
    for m = 1:nmasks
        vis_args = struct();
        vis_args.type = 'subject-fusion';
        vis_args.legend_name = mvpa.vp.visualize;
        vis_args.label_name.title.name = mvpa.cosmo.glm.mask.name{m};
        vis_args.label_name.title.list = mvpa.cosmo.glm.mask.name;
        vis_args.label_name.xlabel = mvpa.vis.meg.fusion.xlabel;
        vis_args.label_name.ylabel = mvpa.vis.meg.fusion.ylabel;
        vis_args.color = cell2mat(values(mvpa.vis.color.fusion,vis_args.legend_name)');
        vis_args.subplot_size = [ceil(length(mvpa.cosmo.glm.mask.name)/mvpa.cosmo.glm.mask.partition) mvpa.cosmo.glm.mask.partition];
        Y = squeeze(result(m,find(ismember(mvpa.vp.xname,mvpa.vp.visualize)),s,:));
        DrawScatter(x,Y,vis_args);
    end
    %-Save all figure in one file
    str.file = eval(mypath.file.figure.fusion.subj);
    str.folder = eval(mypath.folder.result.fusion.figure);
    fn.figure.scatter = strcat(str.folder,'/',str.file);
    saveas(gcf,fn.figure.scatter);
    close all;
end

%-SAVE THE RESULTS
%==========================================================================
%-Save each subject result in a seperate mat file
for s = 1:nsubjs
    subj = meg.subject.list{s};
    for m = 1:nmasks
        roi_name = mvpa.cosmo.glm.mask.name{m};
        if size(squeeze(result(m,:,s,:)),2)==1
            fusion_table_subj.(roi_name) = array2table(squeeze(result(m,:,s,:))');
        else
            fusion_table_subj.(roi_name) = array2table(squeeze(result(m,:,s,:)));
        end
        fusion_table_subj.(roi_name).Properties.VariableNames = compose('time_%d', 1:meg.epoch.ntimes);
        fusion_table_subj.(roi_name).Properties.RowNames = mvpa.vp.xname;
    end
    str.file = eval(mypath.file.subj.fusion);
    str.folder = eval(mypath.folder.result.fusion.file);
    fn.result.fusion = strcat(str.folder,'/',str.file);
    save(fn.result.fusion,'fusion_table_subj');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-ANALYSIS MODE - INDIVIDUAL MRI & GROUP-AVERAGE MEG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode,'mri')

%-CHECK EXTERNALS
%==========================================================================
cosmo_check_external('surfing');
cosmo_check_external('afni');

%-IMPORT FILES AND INFORMATION
%==========================================================================
%-Target and Chunk
targets = mvpa.cosmo.glm.targets;

%-Load Predictor RDMs (Model RDM(s) and Group-Average MEG RDM)
str.file.model = eval(mypath.file.data.rdm.predictor);
str.file.meg = eval(mypath.file.data.rdm.avg.meg);
str.folder = eval(mypath.folder.data.rdm);

%-Load model RDM(s) in 'preds'
load(strcat(str.folder,'\',str.file.model));
for i=1:length(mvpa.rdm.predictor.name)-1
   preds{i} = ...
       eval(sprintf('RDM.%s.%s',...
       mvpa.rdm.predictor.type,...
       mvpa.rdm.predictor.name{i}));
end

%-Load Group-Average MEG RDM
load(strcat(str.folder,'\',str.file.meg),'RDM');
%-Smooth Group-Average MEG RDM
RDM.(sensor) = smooth_rdm(RDM.(sensor));

nsubjs = length(mri.subject.list);
nmasks = length(mvpa.cosmo.glm.mask.name);
ntimes = meg.epoch.ntimes;

%-Create output structure
table_size = [nsubjs nmasks];
table_type = repmat({'cell'},1,nmasks);
fusion_table = table('Size',table_size,...
    'VariableTypes',table_type,...
    'VariableNames',mvpa.cosmo.glm.mask.name,...
    'RowNames',mri.subject.list);

%-LOAD DATA
%==========================================================================
data = load_data_mri(mypath,mvpa,mri,'rdm');

%-ANALYSIS
%==========================================================================
for t = 1:ntimes
    for m = 1:nmasks
        for s = 1:nsubjs
            subj = mri.subject.list{s};

            %-LOAD SUBJECT DATA
            %--------------------------------------------------------------
            %-Raw
            %..............................................................
            %ds_full = read_surf(mypath,mri,mvpa,subj,m,data);
            %ds_full.sa.targets = targets;
            %-Remove constant features
            %ds = cosmo_remove_useless_data(ds_full);
            %..............................................................
            
            %-RDM
            %..............................................................
            ds = data.(subj).ds.(mvpa.cosmo.glm.mask.name{m});
            %..............................................................
            
            %-ANALYSIS
            %--------------------------------------------------------------
            %-Regression
            measure = @MeasureRDMCorr;

            %-Add Group-Average MEG RDM of each sensor and time to 'preds'
            preds{length(mvpa.rdm.predictor.name)} = RDM.(sensor)(:,:,t);
            
            %-Measure arguments
            %..............................................................
            measure_args = struct();
            measure_args.mode = 'mri';
            measure_args.analysis = mvpa.cosmo.glm.analysis.mri;
            measure_args.input_dsm = true;
            measure_args.metric = mvpa.cosmo.glm.metric.mri;
            measure_args.type = mvpa.cosmo.glm.type.mri;
            measure_args.frrsa = mvpa.cosmo.glm.frrsa;
            measure_args.frrsa_dsm = preds{end};   % MEG at time t
            measure_args.center_data = mvpa.cosmo.glm.voxscaling;
            measure_args.vp = mvpa.vp;
            measure_args.target_dsm = preds;
            measure_args.labels = mvpa.rdm.predictor.name';
            %..............................................................

            %-Apply the measure
            ds_fusion = measure(ds,measure_args);

            %-R2 values of each subject
            result(m,:,s,t) = ds_fusion.samples;
        end
    end
end

%-PLOT THE RESULTS
%==========================================================================
x = linspace(1000*meg.epoch.tmin,1000*meg.epoch.tmax,ntimes);
for s = 1:nsubjs
    figure('WindowState','maximized'); pause(0.1);
    for m = 1:nmasks
        vis_args = struct();
        vis_args.type = 'subject-fusion';
        vis_args.legend_name = mvpa.vp.visualize;
        vis_args.label_name.title.name = mvpa.cosmo.glm.mask.name{m};
        vis_args.label_name.title.list = mvpa.cosmo.glm.mask.name;
        vis_args.label_name.xlabel = mvpa.vis.mri.fusion.xlabel;
        vis_args.label_name.ylabel = mvpa.vis.mri.fusion.ylabel;
        vis_args.color = cell2mat(values(mvpa.vis.color.fusion,vis_args.legend_name)');
        vis_args.subplot_size = [ceil(length(mvpa.cosmo.glm.mask.name)/mvpa.cosmo.glm.mask.partition) mvpa.cosmo.glm.mask.partition];
        Y = squeeze(result(m,find(ismember(mvpa.vp.xname,mvpa.vp.visualize)),s,:));
        DrawScatter(x,Y,vis_args);
    end
    %-Save all figure in one file
    str.file = eval(mypath.file.figure.fusion.subj);
    str.folder = eval(mypath.folder.result.fusion.figure);
    fn.figure.scatter = strcat(str.folder,'/',str.file);
    saveas(gcf,fn.figure.scatter);
    close all;
end

%-SAVE THE RESULTS
%==========================================================================
%-Save each subject result in a seperate mat file
for s = 1:nsubjs
    subj = mri.subject.list{s};
    for m = 1:nmasks
        roi_name = mvpa.cosmo.glm.mask.name{m};
        fusion_table_subj.(roi_name) = array2table(squeeze(result(m,:,s,:)));
        fusion_table_subj.(roi_name).Properties.VariableNames = compose('time_%d', 1:meg.epoch.ntimes);
        fusion_table_subj.(roi_name).Properties.RowNames = mvpa.vp.xname;
    end
    str.file = eval(mypath.file.subj.fusion);
    str.folder = eval(mypath.folder.result.fusion.file);
    fn.result.fusion = strcat(str.folder,'/',str.file);
    save(fn.result.fusion,'fusion_table_subj');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-ANALYSIS MODE - INDIVIDUAL FR-MRI & GROUP-AVERAGE MEG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode,'fr-mri')

%-CHECK EXTERNALS
%==========================================================================
cosmo_check_external('surfing');
cosmo_check_external('afni');

%-IMPORT FILES AND INFORMATION
%==========================================================================
%-Target and Chunk
targets = mvpa.cosmo.glm.targets;

%-Load Predictor RDMs (Model RDM(s) and Group-Average MEG RDM)
str.file.model = eval(mypath.file.data.rdm.predictor);
str.file.meg = eval(mypath.file.data.rdm.avg.meg);
str.folder = eval(mypath.folder.data.rdm);

%-Load model RDM(s) in 'preds'
load(strcat(str.folder,'\',str.file.model));
for i=1:length(mvpa.rdm.predictor.name)-1
   preds{i} = ...
       eval(sprintf('RDM.%s.%s',...
       mvpa.rdm.predictor.type,...
       mvpa.rdm.predictor.name{i}));
end

%-Load Group-Average MEG RDM
load(strcat(str.folder,'\',str.file.meg),'RDM');
%-Smooth Group-Average MEG RDM
RDM.(sensor) = smooth_rdm(RDM.(sensor));

nsubjs = length(mri.subject.list);
nmasks = length(mvpa.cosmo.glm.mask.name);
ntimes = meg.epoch.ntimes;

%-Create output structure
table_size = [nsubjs nmasks];
table_type = repmat({'cell'},1,nmasks);
fusion_table = table('Size',table_size,...
    'VariableTypes',table_type,...
    'VariableNames',mvpa.cosmo.glm.mask.name,...
    'RowNames',mri.subject.list);

%-LOAD DATA
%==========================================================================
data = load_data_mri(mypath,mvpa,mri,'fr-rdm');

%-SMOOTH DATA
%==========================================================================
% for m = 1:nmasks
%     for s = 1:nsubjs
%         subj = mri.subject.list{s};
%                 
%         data.(subj).ds.(sensor).(mvpa.cosmo.glm.mask.name{m}) = ...
%             smooth_rdm(data.(subj).ds.(sensor).(mvpa.cosmo.glm.mask.name{m}));
%     end
% end

%-ANALYSIS
%==========================================================================
for t = 1:ntimes
    for m = 1:nmasks
        for s = 1:nsubjs
            subj = mri.subject.list{s};

            %-LOAD SUBJECT DATA
            %--------------------------------------------------------------            
            %-RDM
            %..............................................................
            %ds = data.(subj).ds.(sensor).(mvpa.cosmo.glm.mask.name{m})(:,:,t);
            ds = data.(subj).ds.(sensor).(mvpa.cosmo.glm.mask.name{m});
            %..............................................................
            
            %-ANALYSIS
            %--------------------------------------------------------------
            %-Regression
            measure = @MeasureRDMCorr;

            %-Add Group-Average MEG RDM of each sensor and time to 'preds'
            preds{length(mvpa.rdm.predictor.name)} = RDM.(sensor)(:,:,t);
            
            %-Measure arguments
            %..............................................................
            measure_args = struct();
            measure_args.mode = 'mri';
            measure_args.analysis = mvpa.cosmo.glm.analysis.mri;
            measure_args.input_dsm = true;
            measure_args.metric = mvpa.cosmo.glm.metric.mri;
            measure_args.type = mvpa.cosmo.glm.type.mri;
            measure_args.frrsa = mvpa.cosmo.glm.frrsa;
            measure_args.frrsa_dsm = preds{end};   % MEG at time t
            measure_args.center_data = mvpa.cosmo.glm.voxscaling;
            measure_args.vp = mvpa.vp;
            measure_args.target_dsm = preds;
            measure_args.labels = mvpa.rdm.predictor.name';
            %..............................................................

            %-Apply the measure
            ds_fusion = measure(ds,measure_args);

            %-R2 values of each subject
            result(m,:,s,t) = ds_fusion.samples;
        end
    end
end

%-PLOT THE RESULTS
%==========================================================================
x = linspace(1000*meg.epoch.tmin,1000*meg.epoch.tmax,ntimes);
for s = 1:nsubjs
    figure('WindowState','maximized'); pause(0.1);
    for m = 1:nmasks
        vis_args = struct();
        vis_args.type = 'subject-fusion';
        vis_args.legend_name = mvpa.vp.visualize;
        vis_args.label_name.title.name = mvpa.cosmo.glm.mask.name{m};
        vis_args.label_name.title.list = mvpa.cosmo.glm.mask.name;
        vis_args.label_name.xlabel = mvpa.vis.mri.fusion.xlabel;
        vis_args.label_name.ylabel = mvpa.vis.mri.fusion.ylabel;
        vis_args.color = cell2mat(values(mvpa.vis.color.fusion,vis_args.legend_name)');
        vis_args.subplot_size = [ceil(length(mvpa.cosmo.glm.mask.name)/mvpa.cosmo.glm.mask.partition) mvpa.cosmo.glm.mask.partition];
        Y = squeeze(result(m,find(ismember(mvpa.vp.xname,mvpa.vp.visualize)),s,:));
        DrawScatter(x,Y,vis_args);
    end
    %-Save all figure in one file
    str.file = eval(mypath.file.figure.fusion.subj);
    str.folder = eval(mypath.folder.result.fusion.figure);
    fn.figure.scatter = strcat(str.folder,'/',str.file);
    saveas(gcf,fn.figure.scatter);
    close all;
end

%-SAVE THE RESULTS
%==========================================================================
%-Save each subject result in a seperate mat file
for s = 1:nsubjs
    subj = mri.subject.list{s};
    for m = 1:nmasks
        roi_name = mvpa.cosmo.glm.mask.name{m};
        fusion_table_subj.(roi_name) = array2table(squeeze(result(m,:,s,:)));
        fusion_table_subj.(roi_name).Properties.VariableNames = compose('time_%d', 1:meg.epoch.ntimes);
        fusion_table_subj.(roi_name).Properties.RowNames = mvpa.vp.xname;
    end
    str.file = eval(mypath.file.subj.fusion);
    str.folder = eval(mypath.folder.result.fusion.file);
    fn.result.fusion = strcat(str.folder,'/',str.file);
    save(fn.result.fusion,'fusion_table_subj');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-ANALYSIS MODE - GROUP-AVERAGE MEG & GROUP-AVERAGE MRI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode,'none')

%-LOAD DATA
%==========================================================================
%-Load group-averaged RDM (MRI)
str.file = eval(mypath.file.data.rdm.avg.mri);
str.folder = eval(mypath.folder.data.rdm);
MRI = load(strcat(str.folder,'/',str.file));

%-Load group-averaged RDM (MEG)
str.file = eval(mypath.file.data.rdm.avg.meg);
str.folder = eval(mypath.folder.data.rdm);
MEG = load(strcat(str.folder,'/',str.file));

%-Load predictor RDMs
str.file = eval(mypath.file.data.rdm.predictor);
str.folder = eval(mypath.folder.data.rdm);
load(strcat(str.folder,'\',str.file));
for i=1:length(mvpa.rdm.predictor.name)-1
   predictors{i} = ...
       eval(sprintf('RDM.%s.%s',...
       mvpa.rdm.predictor.type,...
       mvpa.rdm.predictor.name{i}));
end

nmasks = length(mvpa.cosmo.glm.mask.name);
ntimes = size(MEG.RDM.grad,3);

%-ANALYSIS
%==========================================================================
for t = 1:ntimes
    for m = 1:nmasks
        %-Add Group-Average MEG RDM of each sensor and time to 'preds'
        predictors{length(mvpa.rdm.predictor.name)} = MEG.RDM.(sensor)(:,:,t);
        %-Compute R2
        r2(t,m) = MeasureRDMCorr(MRI.RDM.(mvpa.cosmo.glm.mask.name{m}),...
            'mode','mri','target_dsm',predictors,'analysis','vp','vp',mvpa.vp,'input_dsm',true);
    end
end

%-PLOT THE RESULTS
%==========================================================================
%-Plot the Result in One Figure by ROI
%--------------------------------------------------------------------------
x = linspace(1000*meg.epoch.tmin,1000*meg.epoch.tmax,ntimes);
figure('WindowState','maximized'); pause(0.1);
for m = 1:nmasks
    roi_name = mvpa.cosmo.glm.mask.name{m};
    for n = 1:mvpa.vp.nmodels
        for t = 1:ntimes
            r2_temp = r2(t,m).samples;
            y.(roi_name)(n,t) = r2_temp(n);
        end
    end
    
    for model=1:mvpa.vp.nmodels
        model_name = replace(mvpa.vp.xname{model},'-','_');
        y.(model_name)(m,:) = ...
            y.(roi_name)(find(strcmp(mvpa.vp.xname,mvpa.vp.xname{model})),:);
        err.(model_name)(m,:) = 0;
    end
    
    vis_args = struct();
    vis_args.type = 'subject-fusion';
    vis_args.legend_name = mvpa.vp.visualize;
    vis_args.label_name.title.name = mvpa.cosmo.glm.mask.name{m};
    vis_args.label_name.title.list = mvpa.cosmo.glm.mask.name;
    vis_args.label_name.xlabel = 'Time';
    vis_args.label_name.ylabel = 'Correlation Coefficient';
    vis_args.color = cell2mat(values(mvpa.vis.color.fusion,vis_args.legend_name)');
    vis_args.subplot_size = [ceil(length(mvpa.cosmo.glm.mask.name)/mvpa.cosmo.glm.mask.partition) mvpa.cosmo.glm.mask.partition];
    Y = y.(roi_name)(find(ismember(mvpa.vp.xname,mvpa.vp.visualize)),:);
    DrawScatter(x,Y,vis_args);
end

%-Plot the Result in One Figure by Model
%--------------------------------------------------------------------------
figure('WindowState','maximized'); pause(2);
for model=1:length(mvpa.vp.visualize)
     model_name = replace(mvpa.vp.visualize{model},'-','_');
     
    vis_args = struct();
    vis_args.type = 'subject-fusion';
    vis_args.legend_name = mvpa.mask.name(1:mvpa.mask.partition);
    vis_args.label_name.title.name = mvpa.vp.visualize{model};
    vis_args.label_name.title.list =  mvpa.vp.visualize;
    vis_args.label_name.xlabel = 'Time';
    vis_args.label_name.ylabel = 'Correlation Coefficient';
    vis_args.color = cell2mat(values(mvpa.vis.color.roi,vis_args.legend_name)');
    vis_args.subplot_size = [length(mvpa.vp.visualize)+1 1];
    Y = y.(model_name)(1:mvpa.mask.partition,:);
    DrawScatter(x,Y,vis_args);
end
%..........................................................................
figure('WindowState','maximized'); pause(2);
for model=1:length(mvpa.vp.visualize)
    model_name = replace(mvpa.vp.visualize{model},'-','_');
     
    vis_args = struct();
    vis_args.type = 'subject-fusion';
    vis_args.legend_name = mvpa.mask.name(mvpa.mask.partition+1:end);
    vis_args.label_name.title.name = mvpa.vp.visualize{model};
    vis_args.label_name.title.list =  mvpa.vp.visualize;
    vis_args.color = cell2mat(values(mvpa.vis.color.roi,vis_args.legend_name)');
    vis_args.subplot_size = [length(mvpa.vp.visualize)+1 1];
    Y = y.(model_name)(mvpa.mask.partition+1:end,:);
    DrawScatter(x,Y,vis_args);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-ANALYSIS MODE - INDIVIDUAL MEG & INDIVIDUAL MRI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode,'both')
    
%-IMPORT FILES AND INFORMATION
%==========================================================================
%-Target and Chunk
targets = mvpa.cosmo.glm.targets;

%-Load Predictor RDMs (Model RDMs)
str.file.model = eval(mypath.file.data.rdm.predictor);
str.file.mri = eval(mypath.file.data.rdm.avg.mri);
str.folder = eval(mypath.folder.data.rdm);

%-Load model RDM(s) in 'preds'
load(strcat(str.folder,'\',str.file.model));
for i=1:length(mvpa.rdm.predictor.name)-1
   preds{i} = ...
       eval(sprintf('RDM.%s.%s',...
       mvpa.rdm.predictor.type,...
       mvpa.rdm.predictor.name{i}));
end

nmasks = length(mvpa.cosmo.glm.mask.name);
ntimes = meg.epoch.ntimes;
nsubjs = length(meg.subject.list);

%-LOAD DATA
%==========================================================================
data = load_data_meg(mypath,mvpa,meg,sensor,'rdm');

%-ANALYSIS
%==========================================================================
for m = 1:nmasks
    for s = 1:nsubjs
        subj = meg.subject.list{s};

        %-LOAD SUBJECT DATA
        %------------------------------------------------------------------
        %-Raw
        %..................................................................
        %ds = data.(subj).ds;
        %..................................................................
        
        %-RDM
        %..................................................................
        ds = data.(subj).ds;
        %..................................................................

        %-ANALYSIS
        %------------------------------------------------------------------
        %-Add Subject MRI RDM of each ROI to 'preds'
        str.file = eval(mypath.file.subj.mri.rdm);
        str.folder = eval(mypath.folder.result.mri.rdm);
        fn.in = strcat(str.folder,'/',str.file); load(fn.in);
        preds{length(mvpa.rdm.predictor.name)} = RDM.(mvpa.cosmo.glm.mask.name{m});

        %-Regression
        measure = @MeasureRDMCorr;

        %-Measure arguments
        %..................................................................
        measure_args = struct();
        measure_args.mode = 'meg';
        measure_args.analysis = mvpa.cosmo.glm.analysis.meg;
        measure_args.input_dsm = true;
        measure_args.metric = mvpa.cosmo.glm.metric.meg;
        measure_args.type = mvpa.cosmo.glm.type.meg;
        measure_args.frrsa = mvpa.cosmo.glm.frrsa;
        measure_args.frrsa_dsm = preds{end};   % MRI at region m
        measure_args.center_data = mvpa.cosmo.glm.chanscaling;
        measure_args.pseudo = mvpa.cosmo.glm.npseudo;
        measure_args.perm = mvpa.cosmo.glm.nperm;
        measure_args.cv = mvpa.cosmo.glm.cv;
        measure_args.mnn = mvpa.cosmo.glm.mnn;
        measure_args.smooth = mvpa.cosmo.glm.smooth;
        measure_args.vp = mvpa.vp;
        measure_args.target_dsm = preds;
        measure_args.labels = mvpa.rdm.predictor.name';
        %..................................................................

        %-Apply the measure
        ds_fusion = measure(ds,measure_args);
        
        %-R2 values of each subject
        for p = 1:length(ds_fusion(1).sa.labels)
            for t = 1:size(ds_fusion,2)
                %data(MASK,PREDICTORS,SUBJECTS,TIME)
                result(m,p,s,t) = ds_fusion(t).samples(p);
            end
        end
    end
end

%-PLOT THE RESULTS
%==========================================================================
x = linspace(1000*meg.epoch.tmin,1000*meg.epoch.tmax,ntimes);
for s = 1:nsubjs
    figure('WindowState','maximized'); pause(0.1);
    for m = 1:nmasks
        vis_args = struct();
        vis_args.type = 'subject-fusion';
        vis_args.legend_name = mvpa.vp.visualize;
        vis_args.label_name.title.name = mvpa.cosmo.glm.mask.name{m};
        vis_args.label_name.title.list = mvpa.cosmo.glm.mask.name;
        vis_args.label_name.xlabel = mvpa.vis.meg.fusion.xlabel;
        vis_args.label_name.ylabel = mvpa.vis.meg.fusion.ylabel;
        vis_args.color = cell2mat(values(mvpa.vis.color.fusion,vis_args.legend_name)');
        vis_args.subplot_size = [ceil(length(mvpa.cosmo.glm.mask.name)/mvpa.cosmo.glm.mask.partition) mvpa.cosmo.glm.mask.partition];
        Y = squeeze(result(m,find(ismember(mvpa.vp.xname,mvpa.vp.visualize)),s,:));
        DrawScatter(x,Y,vis_args);
    end
    %-Save all figure in one file
    str.file = eval(mypath.file.figure.fusion.subj);
    str.folder = eval(mypath.folder.result.fusion.figure);
    fn.figure.scatter = strcat(str.folder,'/',str.file);
    saveas(gcf,fn.figure.scatter);
    close all;
end

%-SAVE THE RESULTS
%==========================================================================
%-Save each subject result in a seperate mat file
for s = 1:nsubjs
    subj = meg.subject.list{s};
    for m = 1:nmasks
        roi_name = mvpa.cosmo.glm.mask.name{m};
        if size(squeeze(result(m,:,s,:)),2)==1
            fusion_table_subj.(roi_name) = array2table(squeeze(result(m,:,s,:))');
        else
            fusion_table_subj.(roi_name) = array2table(squeeze(result(m,:,s,:)));
        end
        fusion_table_subj.(roi_name).Properties.VariableNames = compose('time_%d', 1:meg.epoch.ntimes);
        fusion_table_subj.(roi_name).Properties.RowNames = mvpa.vp.xname;
    end
    str.file = eval(mypath.file.subj.fusion);
    str.folder = eval(mypath.folder.result.fusion.file);
    fn.result.fusion = strcat(str.folder,'/',str.file);
    save(fn.result.fusion,'fusion_table_subj');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-ANALYSIS MODE - GROUP-AVERAGE RANDOMIZED MEG & GROUP-AVERAGE MRI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mode,'null')

%-LOAD DATA
%==========================================================================
%-Load group-averaged RDM (MRI)
str.file = eval(mypath.file.data.rdm.avg.mri);
str.folder = eval(mypath.folder.data.rdm);
MRI = load(strcat(str.folder,'/',str.file));

%-Load group-averaged RDM (MEG)
str.file = eval(mypath.file.data.rdm.avg.meg);
str.folder = eval(mypath.folder.data.rdm);
MEG = load(strcat(str.folder,'/',str.file));

%-Load predictor RDMs
str.file = eval(mypath.file.data.rdm.predictor);
str.folder = eval(mypath.folder.data.rdm);
load(strcat(str.folder,'\',str.file));
for i=1:length(mvpa.rdm.predictor.name)-1
   predictors{i} = ...
       eval(sprintf('RDM.%s.%s',...
       mvpa.rdm.predictor.type,...
       mvpa.rdm.predictor.name{i}));
end

nmasks = length(mvpa.cosmo.glm.mask.name);
ntimes = size(MEG.RDM.grad,3);

%-ANALYSIS
%==========================================================================
for n = 1:mvpa.vp.null.size
    for m = 1:nmasks
        %-Add Group-Average MEG RDM of each sensor and time to 'preds'
        predictors{length(mvpa.rdm.predictor.name)} = MRI.RDM.(mvpa.cosmo.glm.mask.name{m});
        
        %-Randomize MEG RDM
        meg_rdm = MEG.RDM.(sensor);
        for t = 1:size(meg_rdm,3)
            meg_rdm_vec(:,t) = cosmo_squareform(meg_rdm(:,:,t),'tovector');
        end
        idx = randperm(size(meg_rdm_vec,1));
        meg_rdm_vec_random = meg_rdm_vec(idx,:);
        for t = 1:size(meg_rdm,3)
            meg_rdm_random(:,:,t) = cosmo_squareform(meg_rdm_vec_random(:,t),'tomatrix');
        end
        
        %-Compute R2
        r2 = MeasureRDMCorr(meg_rdm_random,...
            'mode','meg','target_dsm',predictors,'analysis','vp','vp',mvpa.vp,'input_dsm',true);
          
        for p = 1:length(r2(1).sa.labels)
            for t = 1:size(r2,2)
                %data(MASK,PREDICTORS,SUBJECTS,TIME)
                null(m,p,n,t) = r2(t).samples(p);
            end
        end
    end
end

%-SAVE THE RESULTS
%==========================================================================
str.file = eval(mypath.file.data.null.fusion);
str.folder = eval(mypath.folder.result.fusion.null);
fn.result.null = strcat(str.folder,'/',str.file);
description = 'MASK,PREDICTORS,SUBJECTS,TIME';
save(fn.result.null,'null','description');

end

end

%__________________________________________________________________________
function data = load_data_meg(mypath,mvpa,meg,sensor,type)

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

%__________________________________________________________________________
function data = load_data_mri(mypath,mvpa,mri,type)

nsubjs = length(mri.subject.list);
nmasks = length(mvpa.cosmo.glm.mask.name);

if strcmp(type,'raw')
    %-Read Atlas
    %----------------------------------------------------------------------
    str.file = eval(mypath.file.data.atlas);
    str.folder = eval(mypath.folder.data.atlas);
    fn.atlas = strcat(str.folder,'/',str.file);
    [err,msk,Info,Com] = Read_1D(fn.atlas);

    for m = 1:length(mvpa.cosmo.glm.mask.name) 
        data.mask{m}(:,1) = msk(:,1);
        %-Create the mask by putting to one the index of the ROIs of interest
        %and putting to zero the values of other regions.
        data.mask{m}(:,2) = ismember(msk(:,2),mvpa.cosmo.glm.mask.map(mvpa.cosmo.glm.mask.name{m}));
    end

    %-Read Sources
    %----------------------------------------------------------------------
    for s = 1:nsubjs
        subj = mri.subject.list{s};
        
        str.file = eval(mypath.file.subj.mri.source.(mri.mask.voxel.type).dset);
        str.folder = eval(mypath.folder.result.mri.source.(mri.mask.voxel.type));
        fn.source = strcat(str.folder,'/',str.file);
        if strcmp(mri.mask.voxel.type,'firstLevel')
            data.source{s} = cosmo_surface_dataset(fn.source);
        elseif strcmp(mri.mask.voxel.type,'secondLevel')
            data.source{s} = cosmo_surface_dataset(fn.source);
            data.source{s}.samples = data.source{s}.samples(1,:);
        end
    end

    %-Read Subjects' Data (1st-Level T-Map)
    %----------------------------------------------------------------------
    for sbj = 1:nsubjs
        subj = mri.subject.list{sbj};
        for s = mri.analysis.split.num
            for d=1:length(mvpa.cosmo.glm.file.(mri.task.name))
                c = mvpa.cosmo.glm.file.(mri.task.name){d};
                str.file = eval(mypath.file.subj.mri.firstLevel.dset);
                str.folder = eval(mypath.folder.subj.mri.firstLevel);
                fn.data = strcat(str.folder,'/',str.file);
                data.(subj){s,d} = cosmo_surface_dataset(fn.data);
            end
        end
    end
    
elseif strcmp(type,'rdm')
    for s = 1:nsubjs
        subj = mri.subject.list{s};
        
        str.file = eval(mypath.file.subj.mri.rdm);
        str.folder = eval(mypath.folder.result.mri.rdm);
        fn.in = strcat(str.folder,'/',str.file);
        load(fn.in);
        data.(subj).ds = RDM;
    end

elseif strcmp(type,'fr-rdm')
    for s = 1:nsubjs
        subj = mri.subject.list{s};
        
        str.file = eval(mypath.file.subj.frrsa.mri.rdm);
        str.folder = eval(mypath.folder.result.mri.rdm);
        fn.in = strcat(str.folder,'/',str.file);
        load(fn.in);
        data.(subj).ds = RDM;
    end
end

end

%__________________________________________________________________________
function ds_full = read_surf(mypath,mri,mvpa,subj,m,data)

    ds_full = struct();
    
    activeVoxels = read_voxel(mypath,mri,mvpa,subj,m,data);
    
    for s = mri.analysis.split.num
        for d=1:length(mvpa.cosmo.glm.file.(mri.task.name))
            %-Select active voxels
            ds = cosmo_slice(data.(subj){s,d},logical(activeVoxels),2);

            if(isequal(ds_full,struct()))
                ds_full=ds;
            else
                ds_full = cosmo_stack({ds_full,ds});
            end
        end
    end
    
end

%__________________________________________________________________________
function activeVoxels = read_voxel(mypath,mri,mvpa,subj,m,data)

    %-If the desired number of voxels within each ROIs is set to be 0,
    %then select all possible voxels within the ROI.
    if mri.mask.nvox==0
        %-Read the atlas file
        activeVoxels = data.mask{m};
        activeVoxels = activeVoxels(:,2)';
    else
        %-Read the atlas file
        mask = data.mask{m};
        
        %-Read the source file
        source = data.source{find(strcmp(mri.subject.list,subj))};
        
        %-Looking for voxels that are in the mask and have a nonzero t-value
        inmask=find(source.samples~=0 & mask(:,2)'>0);
        mask_vol_nonzeroVox=zeros(1,length(source.samples));
        mask_vol_nonzeroVox(inmask)=source.samples(inmask);
        [sorted,inds]=sort(mask_vol_nonzeroVox(inmask));
        if length(sorted) < mri.mask.nvox
            disp('Insufficient number of voxels in mask!');
        end
        mask_vol_activeVox=zeros(1,length(source.samples));
        mask_vol_activeVox(inmask(inds(length(sorted)-(mri.mask.nvox-1):length(sorted))))=1;
        source.samples = mask_vol_activeVox;
            
        activeVoxels = source.samples;
        
    end
end

%__________________________________________________________________________
function smoothed_rdm=smooth_rdm(rdm)
    for i=1:size(rdm,1)
        for j=1:size(rdm,2)
            if i>j
                for t=1:size(rdm,3)
                    rdm_vector(t)=rdm(i,j,t);
                end
                
                %-Now we average the vector across 3 timepoints
                tp2avg = 3; % smooth over n timepoints
                smoothed_rdm_vector=conv([rdm_vector(1:2),rdm_vector,rdm_vector(end-1:end)],ones(1,tp2avg).*(1/tp2avg),'same');
                smoothed_rdm_vector=smoothed_rdm_vector(3:end-2);
                rdm_vector=smoothed_rdm_vector;
                
                %-Now we put it back into matrix for each t
                for t=1:size(rdm,3)
                    smoothed_rdm(i,j,t)=rdm_vector(t);
                    smoothed_rdm(j,i,t)=rdm_vector(t);
                end
            end
        end
    end
end
