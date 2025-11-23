%__________________________________________________________________________
function meg_cosmo_decoding_sensor(mypath,meg,mvpa,sensor,dim)

%-IMPORT FILES AND INFORMATION
%==========================================================================
nsubjs = length(meg.subject.list);

for s = 1:nsubjs
    subj = meg.subject.list{s};
    
    %-LOAD DATA
    %======================================================================
    str.file = eval(mypath.file.subj.meg.preprocessed);
    str.folder = eval(mypath.folder.subj.meg);
    fn.in = strcat(str.folder,'/',str.file);
    ds=load(fn.in);
    ds = convert_cosmo_ds(ds,dim,meg);
    time = ds.fa.time;
    
    if strcmp(sensor,'grad')
        chanx='meg_axial';
    elseif strcmp(sensor,'mag')
        chanx='meg_planar';
    end
    
    layout=cosmo_meeg_find_layout(ds,'chantype',chanx);
    layout_label = layout.label(1:end-2,:);
    feature_msk = cosmo_dim_match(ds,'chan',layout_label);
    ds = cosmo_slice(ds,feature_msk,2);
    
    %-ANALYSIS
    %======================================================================
    %-Find the number of trials for each samples
    n_unique = length(unique(ds.sa.targets));
    [cnt_unique,~] = hist(ds.sa.targets,unique(ds.sa.targets));
    %-As the number of trials for each samples is equal, choose the first one
    cnt_unique = cnt_unique(1);
    %-Initialize zero array of accuracy across time
    acc_sum = zeros(1,length(unique(time)));
    
    for r=1:mvpa.cosmo.decoding.nperm
        %-Classifier
        measure = @cosmo_crossvalidation_measure;
        
        %..................................................................
        measure_args = struct();
        measure_args.ratio = 1/mvpa.cosmo.decoding.npseudo;
        measure_args.resamplings = 1;
        measure_args.split_by = {'targets','chunks'};
        ds_avg = cosmo_average_samples(ds,measure_args);
        %..................................................................
        
        %-Measure arguments (without feature selection)
        %..................................................................
        %measure_args = struct();
        %measure_args.partitions = @cosmo_nfold_partitioner;
        %measure_args.output = 'accuracy';
        %measure_args.classifier = mvpa.cosmo.decoding.clf.type;
        %measure_args.pca_explained_ratio = 0.95;
        %..................................................................
        
        %-Measure arguments (with feature selection)
        %..................................................................
        %measure_args = struct();
        %measure_args.partitions = @cosmo_nfold_partitioner;
        %measure_args.output = 'accuracy';
        %measure_args.classifier = @cosmo_classify_meta_feature_selection;
        %measure_args.child_classifier = mvpa.cosmo.decoding.clf.type;
        %measure_args.feature_selector = @cosmo_anova_feature_selector;
        %measure_args.feature_selection_ratio_to_keep=0.6;
        %..................................................................
        
        chan_nbrhood = cosmo_meeg_chan_neighborhood(ds_avg,'count',length(layout.label)-2,'chantype',chanx);
        time_nbrhood = cosmo_interval_neighborhood(ds_avg,'time','radius',0);
        nbrhood = cosmo_cross_neighborhood(ds_avg,{chan_nbrhood,time_nbrhood});

        %-Run searchlight
        ds_sl = cosmo_searchlight(ds_avg,nbrhood,measure,measure_args);

        %-All channels have the same accuracy due to the searchligh
        %on all channels
        tf_map = cosmo_map2meeg(ds_sl);
        %-It would be fine if one choose one channel's accuracy instead
        %of averaging as all channels have the same accuracy
        acc = mean(tf_map.avg,1);
        acc_sum = acc_sum + acc;
        
    end

    %-Find the average accuracy across permutation
    acc_avg = acc_sum / mvpa.cosmo.decoding.nperm;

    %-PLOT THE RESULTS
    %======================================================================
    vis_args = struct();
    vis_args.line = 1/n_unique;
    vis_args.type = 'meg-subject-decoding';
    vis_args.label_name.xlabel = mvpa.vis.meg.decoding.xlabel;
    vis_args.label_name.ylabel = mvpa.vis.meg.decoding.ylabel;
    vis_args.label_name.title = mvpa.vis.meg.decoding.title;
    DrawScatter(time,acc_avg,vis_args);

    str.file = eval(mypath.file.figure.meg.decoding.subj);
    str.folder = eval(mypath.folder.result.meg.figure);
    fn.figure.scatter = strcat(str.folder,'/',str.file);
    saveas(gcf,fn.figure.scatter);
    clf;
    
    %-SAVE THE RESULTS
    %======================================================================
    %-Save each subject result in a seperate mat file
    acc_table_subj = array2table(acc_avg);
    acc_table_subj.Properties.VariableNames = compose('time_%d', 1:length(acc_avg));
    acc_table_subj.Properties.RowNames = {'SVM'};
    
    str.file = eval(mypath.file.subj.meg.decoding.sensor);
    str.folder = eval(mypath.folder.result.meg.decoding);
    fn.result.decoding = strcat(str.folder,'/',str.file);
    save(fn.result.decoding,'acc_table_subj');
end

end

%__________________________________________________________________________
function ds_cosmo = convert_cosmo_ds(ds,dim,meg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%https://www.cosmomvpa.org/faq.html#faq-get-ecog-data-in-cosmomvpa-struct
%-This function is input structure dependent!
%-It means that 'trial', 'channel_name', 'time', and 'index' depends on
%naming convention of 'ds'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ds_cosmo = cosmo_flatten(ds.trial,{'chan','time'},{cellstr(ds.channel_name),ds.time});
    ds_cosmo.a.meeg.samples_field = 'trial';
    ds_cosmo.sa.targets = double(ds.index)';
    
    %-Set the label for the specified dimension (N/S/TFA) 
    switch dim
        case 'Number'
            for i=1:length(ds_cosmo.sa.targets)
                if any(ds_cosmo.sa.targets(i)==meg.stim.n.N1(:))
                    ds_cosmo.sa.targets(i)=1;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.n.N2(:))
                    ds_cosmo.sa.targets(i)=2;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.n.N3(:))
                    ds_cosmo.sa.targets(i)=3;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.n.N4(:))
                    ds_cosmo.sa.targets(i)=4;
                else
                    warning("No Matching Parameter!");
                end
            end
        case 'Size'
            for i=1:length(ds_cosmo.sa.targets)
                if any(ds_cosmo.sa.targets(i)==meg.stim.s.S1(:))
                    ds_cosmo.sa.targets(i)=1;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.s.S2(:))
                    ds_cosmo.sa.targets(i)=2;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.s.S3(:))
                    ds_cosmo.sa.targets(i)=3;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.s.S4(:))
                    ds_cosmo.sa.targets(i)=4;
                else
                    warning("No Matching Parameter!");
                end
            end
        case 'TFA'
            for i=1:length(ds_cosmo.sa.targets)
                if any(ds_cosmo.sa.targets(i)==meg.stim.tfa.TFA1(:))
                    ds_cosmo.sa.targets(i)=1;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.tfa.TFA2(:))
                    ds_cosmo.sa.targets(i)=2;
                else
                    warning("No Matching Parameter!");
                end
            end
        otherwise
            warning(strcat('Decoding for ',dim,' is not implemented yet!'));
    end
    
    %-Chunkize the dataset
    ds_cosmo.sa.chunks = (1:size(ds_cosmo.samples,1))';
    switch dim
        case 'Number'
            ds_cosmo.sa.chunks = cosmo_chunkize(ds_cosmo,4);
        case 'Size'
            ds_cosmo.sa.chunks = cosmo_chunkize(ds_cosmo,4);
        case 'TFA'
            ds_cosmo.sa.chunks = cosmo_chunkize(ds_cosmo,2);
    end
end
