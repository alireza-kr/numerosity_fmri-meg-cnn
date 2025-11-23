function meg_glm_source(mypath,mri,meg,mvpa)

%-IMPORT FILES AND INFORMATION
%==========================================================================
%-Target and Chunk
targets = mvpa.cosmo.glm.targets;

%-Load model RDMs
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

nmasks = length(mvpa.cosmo.glm.mask.name);
ntimes = meg.epoch.ntimes;
nsubjs = length(meg.subject.list);

%-LOAD DATA
%==========================================================================
data = load_data(mypath,mvpa,mri,meg);

%-ANALYSIS
%==========================================================================
for m = 1:nmasks
    for s = 1:nsubjs
        subj = meg.subject.list{s};

        %-LOAD SUBJECT DATA (RDM)
        %------------------------------------------------------------------
        ds = data.(subj).ds.(mvpa.cosmo.glm.mask.name{m});

        %-ANALYSIS
        %------------------------------------------------------------------
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
        measure_args.center_data = mvpa.cosmo.glm.chanscaling;
        measure_args.pseudo = mvpa.cosmo.glm.npseudo;
        measure_args.perm = mvpa.cosmo.glm.nperm;
        measure_args.cv = mvpa.cosmo.glm.cv;
        measure_args.mnn = mvpa.cosmo.glm.mnn;
        measure_args.smooth = mvpa.cosmo.glm.smooth;
        measure_args.vp = mvpa.vp;
        measure_args.target_dsm = pred_cell(2,:);
        measure_args.labels = pred_cell(1,:)';
        %..................................................................

        %-Apply the measure
        ds_source = measure(ds,measure_args);
        
        %-R2 values of each subject
        for p = 1:length(ds_source(1).sa.labels)
            for t = 1:size(ds_source,2)
                %data(MASK,PREDICTORS,SUBJECTS,TIME)
                result(m,p,s,t) = ds_source(t).samples(p);
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
        vis_args.legend_name = pred_cell(1,:);
        vis_args.label_name.title.name = mvpa.cosmo.glm.mask.name{m};
        vis_args.label_name.title.list = mvpa.cosmo.glm.mask.name;
        vis_args.label_name.xlabel = mvpa.vis.meg.fusion.xlabel;
        vis_args.label_name.ylabel = mvpa.vis.meg.fusion.ylabel;
        vis_args.color = cell2mat(values(mvpa.vis.color.glm,vis_args.legend_name)');
        vis_args.subplot_size = [ceil(length(mvpa.cosmo.glm.mask.name)/mvpa.cosmo.glm.mask.partition) mvpa.cosmo.glm.mask.partition];
        Y = squeeze(result(m,:,s,:));
        DrawScatter(x,Y,vis_args);
    end
    %-Save all figure in one file
    str.file = eval(mypath.file.figure.meg.glm.source.subj);
    str.folder = eval(mypath.folder.result.meg.figure);
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
            source_table_subj.(roi_name) = array2table(squeeze(result(m,:,s,:))');
        else
            source_table_subj.(roi_name) = array2table(squeeze(result(m,:,s,:)));
        end
        source_table_subj.(roi_name).Properties.VariableNames = compose('time_%d', 1:meg.epoch.ntimes);
        source_table_subj.(roi_name).Properties.RowNames = pred_cell(1,:);
    end
    str.file = eval(mypath.file.subj.meg.glm.source);
    str.folder = eval(mypath.folder.result.meg.glm);
    fn.result.fusion = strcat(str.folder,'/',str.file);
    save(fn.result.fusion,'source_table_subj');
end

end

%__________________________________________________________________________
function data = load_data(mypath,mvpa,mri,meg)

nsubjs = length(meg.subject.list);

for s = 1:nsubjs
    subj = meg.subject.list{s};

    str.file = eval(mypath.file.subj.meg.sources.rdm);
    str.folder = eval(mypath.folder.result.meg.rdm);
    fn.in = strcat(str.folder,'/',str.file);
    load(fn.in);
    data.(subj).ds = RDM;
end

end
