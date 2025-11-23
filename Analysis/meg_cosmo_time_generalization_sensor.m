%__________________________________________________________________________
function meg_cosmo_time_generalization_sensor(mypath,meg,mvpa,sensor,dim)

%-IMPORT FILES AND INFORMATION
%==========================================================================
nsubjs = length(meg.subject.list);
nclfs = length(mvpa.cosmo.decoding.clf.name);

for c=1:nclfs
    for s = 1:nsubjs
        subj = meg.subject.list{s};

        %-LOAD DATA
        %======================================================================
        str.file = eval(mypath.file.subj.meg.preprocessed);
        str.folder = eval(mypath.folder.subj.meg);
        fn.in = strcat(str.folder,'/',str.file);
        ds = load(fn.in);
        ds = convert_cosmo_ds(ds,dim,meg);
        time = ds.fa.time;

        if strcmp(sensor,'grad')
            chanx='meg_axial';
        elseif strcmp(sensor,'mag')
            chanx='meg_planar';
        end

        layout = cosmo_meeg_find_layout(ds,'chantype',chanx);
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
        acc_sum = zeros(length(unique(time)),length(unique(time)));

        for r=1:mvpa.cosmo.decoding.nperm
            %-Classifier
            measure = @cosmo_dim_generalization_measure;

            %..................................................................
            measure_args = struct();
            measure_args.ratio = 1/mvpa.cosmo.decoding.npseudo;
            measure_args.resamplings = 1;
            measure_args.split_by = {'targets','chunks'};
            %..................................................................
            ds_avg = cosmo_average_samples(ds,measure_args);

            %..................................................................
            % Now use cosmo_dim_transpose to make 'time' a sample dimension,
            % and assign to ds_tr
            ds_tr = cosmo_dim_transpose(ds_avg,'time',1);
            %..................................................................

            %-Measure arguments (without feature selection)
            %..................................................................
            measure_args = struct();
            measure_args.measure = @cosmo_crossvalidation_measure;
            measure_args.output = 'accuracy';
            measure_args.dimension = 'time';
            measure_args.classifier = mvpa.cosmo.decoding.clf.type{c};
            %measure_args.pca_explained_ratio = 0.95;
            %..................................................................

            %-Measure arguments (with feature selection)
            %..................................................................
            %measure_args = struct();
            %measure_args.measure = @cosmo_crossvalidation_measure;
            %measure_args.output = 'accuracy';
            %measure_args.dimension = 'time';
            %measure_args.classifier = @cosmo_classify_meta_feature_selection;
            %measure_args.child_classifier = mvpa.cosmo.decoding.clf.type{c};
            %measure_args.feature_selector = @cosmo_anova_feature_selector;
            %measure_args.feature_selection_ratio_to_keep=0.6;
            %..................................................................

            %-Now apply the measure to the dataset
            ds_sl = measure(ds_tr,measure_args);
            [acc,labels,values] = cosmo_unflatten(ds_sl,1);

            acc_sum = acc_sum + acc;

        end

        %-Find the average accuracy across permutation
        acc_subj = acc_sum / mvpa.cosmo.decoding.nperm;

        %-PLOT THE RESULTS
        %======================================================================
        vis_args = struct();
        vis_args.xlabel = strrep(labels{2},'_',' ');
        vis_args.ylabel = strrep(labels{1},'_',' ');
        vis_args.xtick = round(linspace(1, numel(values{2}), 10));
        vis_args.ytick = round(linspace(1, numel(values{1}), 10));
        vis_args.xticklabel = values{2}(vis_args.xtick);
        vis_args.yticklabel = values{1}(vis_args.ytick);
        vis_args.range = [min(min(acc_subj)) max(max(acc_subj))];
        vis_args.title = sprintf('%s -- Accuracy%',upper(mvpa.cosmo.decoding.clf.name{1}));
        DrawMatrix(acc_subj,vis_args);

        str.file = eval(mypath.file.figure.meg.tg.subj);
        str.folder = eval(mypath.folder.result.meg.figure);
        fn.figure.scatter = strcat(str.folder,'/',str.file);
        saveas(gcf,fn.figure.scatter);
        clf;

        %-SAVE THE RESULTS
        %======================================================================
        %-Save each subject result in a seperate mat file
        str.file = eval(mypath.file.subj.meg.tg.sensor);
        str.folder = eval(mypath.folder.result.meg.decoding);
        fn.result.decoding = strcat(str.folder,'/',str.file);
        save(fn.result.decoding,'acc_subj');
    end
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
        case 'TSA'
            for i=1:length(ds_cosmo.sa.targets)
                if any(ds_cosmo.sa.targets(i)==meg.stim.n.N1(:)) && any(ds_cosmo.sa.targets(i)==meg.stim.s.S1(:))
                    ds_cosmo.sa.targets(i)=1;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.n.N2(:)) && any(ds_cosmo.sa.targets(i)==meg.stim.s.S1(:))
                    ds_cosmo.sa.targets(i)=2;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.n.N3(:)) && any(ds_cosmo.sa.targets(i)==meg.stim.s.S1(:))
                    ds_cosmo.sa.targets(i)=3;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.n.N4(:)) && any(ds_cosmo.sa.targets(i)==meg.stim.s.S1(:))
                    ds_cosmo.sa.targets(i)=4;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.n.N1(:)) && any(ds_cosmo.sa.targets(i)==meg.stim.s.S2(:))
                    ds_cosmo.sa.targets(i)=5;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.n.N2(:)) && any(ds_cosmo.sa.targets(i)==meg.stim.s.S2(:))
                    ds_cosmo.sa.targets(i)=6;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.n.N3(:)) && any(ds_cosmo.sa.targets(i)==meg.stim.s.S2(:))
                    ds_cosmo.sa.targets(i)=7;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.n.N4(:)) && any(ds_cosmo.sa.targets(i)==meg.stim.s.S2(:))
                    ds_cosmo.sa.targets(i)=8;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.n.N1(:)) && any(ds_cosmo.sa.targets(i)==meg.stim.s.S3(:))
                    ds_cosmo.sa.targets(i)=9;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.n.N2(:)) && any(ds_cosmo.sa.targets(i)==meg.stim.s.S3(:))
                    ds_cosmo.sa.targets(i)=10;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.n.N3(:)) && any(ds_cosmo.sa.targets(i)==meg.stim.s.S3(:))
                    ds_cosmo.sa.targets(i)=11;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.n.N4(:)) && any(ds_cosmo.sa.targets(i)==meg.stim.s.S3(:))
                    ds_cosmo.sa.targets(i)=12;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.n.N1(:)) && any(ds_cosmo.sa.targets(i)==meg.stim.s.S4(:))
                    ds_cosmo.sa.targets(i)=13;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.n.N2(:)) && any(ds_cosmo.sa.targets(i)==meg.stim.s.S4(:))
                    ds_cosmo.sa.targets(i)=14;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.n.N3(:)) && any(ds_cosmo.sa.targets(i)==meg.stim.s.S4(:))
                    ds_cosmo.sa.targets(i)=15;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.n.N4(:)) && any(ds_cosmo.sa.targets(i)==meg.stim.s.S4(:))
                    ds_cosmo.sa.targets(i)=16;
                else
                    warning("No Matching Parameter!");
                end
            end
        case 'Density'
            for i=1:length(ds_cosmo.sa.targets)
                if any(ds_cosmo.sa.targets(i)==meg.stim.n.N1(:)) && any(ds_cosmo.sa.targets(i)==meg.stim.tfa.TFA1(:))
                    ds_cosmo.sa.targets(i)=1;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.n.N2(:)) && any(ds_cosmo.sa.targets(i)==meg.stim.tfa.TFA1(:))
                    ds_cosmo.sa.targets(i)=2;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.n.N3(:)) && any(ds_cosmo.sa.targets(i)==meg.stim.tfa.TFA1(:))
                    ds_cosmo.sa.targets(i)=3;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.n.N4(:)) && any(ds_cosmo.sa.targets(i)==meg.stim.tfa.TFA1(:))
                    ds_cosmo.sa.targets(i)=4;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.n.N1(:)) && any(ds_cosmo.sa.targets(i)==meg.stim.tfa.TFA2(:))
                    ds_cosmo.sa.targets(i)=5;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.n.N2(:)) && any(ds_cosmo.sa.targets(i)==meg.stim.tfa.TFA2(:))
                    ds_cosmo.sa.targets(i)=6;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.n.N3(:)) && any(ds_cosmo.sa.targets(i)==meg.stim.tfa.TFA2(:))
                    ds_cosmo.sa.targets(i)=7;
                elseif any(ds_cosmo.sa.targets(i)==meg.stim.n.N4(:)) && any(ds_cosmo.sa.targets(i)==meg.stim.tfa.TFA2(:))
                    ds_cosmo.sa.targets(i)=8;
                else
                    warning("No Matching Parameter!");
                end
            end
        otherwise
            warning(strcat('Decoding for ',dim,' is not implemented yet!'));
    end
    
    %-Chunkize the dataset
    % Reduce the number of chunks to have only two chunks
    ds_cosmo.sa.chunks = (1:size(ds_cosmo.samples,1))';
    switch dim
        case 'Number'
            ds_cosmo.sa.chunks = cosmo_chunkize(ds_cosmo,2);
        case 'Size'
            ds_cosmo.sa.chunks = cosmo_chunkize(ds_cosmo,2);
        case 'TFA'
            ds_cosmo.sa.chunks = cosmo_chunkize(ds_cosmo,2);
        case 'TSA'
            ds_cosmo.sa.chunks = cosmo_chunkize(ds_cosmo,2);
        case 'Density'
            ds_cosmo.sa.chunks = cosmo_chunkize(ds_cosmo,2);
    end
end
