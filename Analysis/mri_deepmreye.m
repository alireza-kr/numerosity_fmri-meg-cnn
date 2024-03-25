params.tr = 2; % 2.03
params.num_run = 6;
params.stim_dur = 0.5; % stimulus duration in second
params.gaze_subSample = 10;
params.gaze_len = 11760/params.num_run;
params.gaze_len_subSample = params.gaze_len/params.gaze_subSample;
params.exp_len = round(params.gaze_len_subSample*params.tr);

params.num.key = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 ...
    25 26 27 28 29 30 31 32];
params.num.val = [6 6 6 6 6 6 6 6 10 10 10 10 10 10 10 10 17 17 17 17 17 17 17 17 ...
    29 29 29 29 29 29 29 29];
params.num.dict = containers.Map(params.num.key,params.num.val);
params.num.uval = unique(params.num.val);
params.num.label = ["Six";"Ten";"Seventeen";"Twenty-Nine"];

subjs = mri.subject.list;
subjs = setdiff(subjs,{'Subj38'});

str.file = eval(mypath.file.data.gaze);
str.folder = eval(mypath.folder.data.deepmreye);
gaze = load(strcat(str.folder,'/',str.file));

%-Create output structure
varNames = {'MeanHorizental','MeanVertical','STDHorizental','STDVertical'};
result.withOutlier = table('Size',[numel(subjs)+1 numel(varNames)],...
    'VariableTypes',repmat({'cell'},1,numel(varNames)),...
    'VariableNames',varNames,...
    'RowNames',[subjs,'Mean']);
result.withoutOutlier = table('Size',[numel(subjs)+1 numel(varNames)],...
    'VariableTypes',repmat({'cell'},1,numel(varNames)),...
    'VariableNames',varNames,...
    'RowNames',[subjs,'Mean']);

for p = 1:numel(subjs)
subj = subjs{p};

str.file = eval(mypath.file.subj.response.summary);
str.folder = eval(mypath.folder.subj.response);
response = load(strcat(str.folder,'/',str.file));

%-Make stimulus information (onset and values)
%--------------------------------------------------------------------------
for r = 1:params.num_run
   Stim(r).Info(1,:) = round(response.Matrix.Run(r).OnsetAll);

   for s = 1:length(response.Matrix.Run(r).OnsetAll)-1
      if any(response.Matrix.Run(r).OnsetMatch == response.Matrix.Run(r).OnsetAll(s))
          Stim(r).Info(2,s) = 0;
      else
          Stim(r).Info(2,s) = params.num.dict(response.Matrix.Run(r).CatAll1D(s));
      end
   end
   Stim(r).Info(2,length(response.Matrix.Run(r).OnsetAll)) = 0;
end

%-Make stimulus time series
%--------------------------------------------------------------------------
for r = 1:params.num_run
    stim = 0;
    for t = 1:params.exp_len
        if any(Stim(r).Info(1,:) == t)
            stim = Stim(r).Info(2,Stim(r).Info(1,:) == t);
        end
        Stim(r).series(1,t) = stim;
        % If you want to only have the result during stimulus presentation,
        % uncomment the following line
        stim = 0;
    end
end

%-Subsample gaze data by taking the median of ten elements
%--------------------------------------------------------------------------
gaze_h = gaze.(subj)(:,1);
gaze_v = gaze.(subj)(:,2);
m = numel(gaze_h);
gaze_h_subsample = nanmedian(reshape([gaze_h(:);nan(mod(-m,params.gaze_subSample),1)],params.gaze_subSample,[]));
gaze_v_subsample = nanmedian(reshape([gaze_v(:);nan(mod(-m,params.gaze_subSample),1)],params.gaze_subSample,[]));

%-Upsample gaze data by TR
%--------------------------------------------------------------------------
x_coords = 1:length(gaze_h_subsample);
new_x_coords = linspace(1, length(gaze_h_subsample), 2*length(gaze_h_subsample));
gaze_position(:,1) = interp1(x_coords, gaze_h_subsample, new_x_coords, 'linear');
gaze_position(:,2) = interp1(x_coords, gaze_v_subsample, new_x_coords, 'linear');

%-Add gaze data to stimulus time series
%--------------------------------------------------------------------------
gaze_len_upSample = 2*params.gaze_len_subSample;

Stim(1).series(2,:) = gaze_position(1:gaze_len_upSample,1);
Stim(2).series(2,:) = gaze_position(gaze_len_upSample+1:2*gaze_len_upSample,1);
Stim(3).series(2,:) = gaze_position(2*gaze_len_upSample+1:3*gaze_len_upSample,1);
Stim(4).series(2,:) = gaze_position(3*gaze_len_upSample+1:4*gaze_len_upSample,1);
Stim(5).series(2,:) = gaze_position(4*gaze_len_upSample+1:5*gaze_len_upSample,1);
Stim(6).series(2,:) = gaze_position(5*gaze_len_upSample+1:6*gaze_len_upSample,1);

Stim(1).series(3,:) = gaze_position(1:gaze_len_upSample,2);
Stim(2).series(3,:) = gaze_position(gaze_len_upSample+1:2*gaze_len_upSample,2);
Stim(3).series(3,:) = gaze_position(2*gaze_len_upSample+1:3*gaze_len_upSample,2);
Stim(4).series(3,:) = gaze_position(3*gaze_len_upSample+1:4*gaze_len_upSample,2);
Stim(5).series(3,:) = gaze_position(4*gaze_len_upSample+1:5*gaze_len_upSample,2);
Stim(6).series(3,:) = gaze_position(5*gaze_len_upSample+1:6*gaze_len_upSample,2);

%==========================================================================
%-With outliers positions
%==========================================================================

%-Draw the gaze position
%--------------------------------------------------------------------------
% for r = 1:params.num_run
%     for g = 1:gaze_len_upSample
%         if Stim(r).series(1,g)~=0
%             if Stim(r).series(1,g)==6
%                 plot(Stim(r).series(2,g),Stim(r).series(3,g),'.','Color','green');
%                 hold on;
%             elseif Stim(r).series(1,g)==10
%                 plot(Stim(r).series(2,g),Stim(r).series(3,g),'.','Color','red');
%                 hold on;
%             elseif Stim(r).series(1,g)==17
%                 plot(Stim(r).series(2,g),Stim(r).series(3,g),'.','Color','blue');
%                 hold on;
%             elseif Stim(r).series(1,g)==29
%                 plot(Stim(r).series(2,g),Stim(r).series(3,g),'.','Color','black');
%                 hold on;
%             end
%         end
%     end
% end

%-Find average and STD
%--------------------------------------------------------------------------
% Concatenate runs
series = [];
for r = 1:params.num_run
    series = [series Stim(r).series()];
end

MeanHorizontal = rescale([mean(series(2,series(1,:)==6));mean(series(2,series(1,:)==10));...
    mean(series(2,series(1,:)==17));mean(series(2,series(1,:)==29))]);
MeanVertical = rescale([mean(series(3,series(1,:)==6));mean(series(3,series(1,:)==10));...
    mean(series(3,series(1,:)==17));mean(series(3,series(1,:)==29))]);
STDHorizontal = [std(series(2,series(1,:)==6));std(series(2,series(1,:)==10));...
    std(series(2,series(1,:)==17));std(series(2,series(1,:)==29))];
STDVertical = [std(series(3,series(1,:)==6));std(series(3,series(1,:)==10));...
    std(series(3,series(1,:)==17));std(series(3,series(1,:)==29))];

result.withOutlier(p,:) = {MeanHorizontal,MeanVertical,STDHorizontal,STDVertical};

%==========================================================================
%-Without outliers positions
%==========================================================================

outlier_h = isoutlier(series(2,:));
outlier_v = isoutlier(series(3,:));
outlier = or(outlier_h,outlier_v);

series(:,outlier) = [];

%-Draw the gaze position
%--------------------------------------------------------------------------
% figure;
% for g = 1:length(series)
%     if series(1,g)~=0
%         if series(1,g)==6
%             plot(series(2,g),series(3,g),'.','Color','green');
%             hold on;
%         elseif series(1,g)==10
%             plot(series(2,g),series(3,g),'.','Color','red');
%             hold on;
%         elseif series(1,g)==17
%             plot(series(2,g),series(3,g),'.','Color','blue');
%             hold on;
%         elseif series(1,g)==29
%             plot(series(2,g),series(3,g),'.','Color','black');
%             hold on;
%         end
%     end
% end

%-Find average and STD
%--------------------------------------------------------------------------
MeanHorizental = rescale([mean(series(2,series(1,:)==6));mean(series(2,series(1,:)==10));...
    mean(series(2,series(1,:)==17));mean(series(2,series(1,:)==29))]);
MeanVertical = rescale([mean(series(3,series(1,:)==6));mean(series(3,series(1,:)==10));...
    mean(series(3,series(1,:)==17));mean(series(3,series(1,:)==29))]);
STDHorizental = [std(series(2,series(1,:)==6));std(series(2,series(1,:)==10));...
    std(series(2,series(1,:)==17));std(series(2,series(1,:)==29))];
STDVertical = [std(series(3,series(1,:)==6));std(series(3,series(1,:)==10));...
    std(series(3,series(1,:)==17));std(series(3,series(1,:)==29))];

result.withoutOutlier(p,:) = {MeanHorizontal,MeanVertical,STDHorizontal,STDVertical};

end

%==========================================================================
%-Mean values of positions across subjects
%==========================================================================

a.withOutlier = table2array(result.withOutlier(1:numel(subjs),:));
a.withoutOutlier = table2array(result.withoutOutlier(1:numel(subjs),:));

result.withOutlier(p+1,:) = ...
    {[mean(reshape(cell2mat(a.withOutlier(1:p,1)),[4,numel(subjs)])')]',...
    [mean(reshape(cell2mat(a.withOutlier(1:p,2)),[4,numel(subjs)])')]',...
    [mean(reshape(cell2mat(a.withOutlier(1:p,3)),[4,numel(subjs)])')]',...
    [mean(reshape(cell2mat(a.withOutlier(1:p,4)),[4,numel(subjs)])')]'};

result.withoutOutlier(p+1,:) = ...
    {[mean(reshape(cell2mat(a.withoutOutlier(1:p,1)),[4,numel(subjs)])')]',...
    [mean(reshape(cell2mat(a.withoutOutlier(1:p,2)),[4,numel(subjs)])')]',...
    [mean(reshape(cell2mat(a.withoutOutlier(1:p,3)),[4,numel(subjs)])')]',...
    [mean(reshape(cell2mat(a.withoutOutlier(1:p,4)),[4,numel(subjs)])')]'};

%==========================================================================
%-T-Test
%==========================================================================

data = reshape(cell2mat(a.withoutOutlier(:,1)),[4,numel(subjs)])';
%[h,p] = ttest2(data(:,1),data(:,4),'Tail','right');
[h,p] = ttest(data(:,1),0.5,"Tail","right");
[h,p] = ttest(data(:,4),0.5,"Tail","left");

%==========================================================================
%-Linear Fit
%==========================================================================

eq = 'a*x+b';
options = fitoptions(eq,'Lower',[-Inf -Inf],'Upper',[Inf Inf],'StartPoint',[0 0]);
Y = [data(:,1);data(:,2);data(:,3);data(:,4)];
X = [repmat(1,30,1);repmat(2,30,1);repmat(3,30,1);repmat(4,30,1)];

boxplot(data,'Labels',params.num.label);
lines = findobj(gcf,'type','line');
boxes = findobj(gca,'Tag','Box');
for j=1:length(boxes)
    patch(get(boxes(j),'XData'),get(boxes(j),'YData'),'blue','FaceAlpha',.2);
end
set(lines,'Color','black');
set(findobj(gca,'type','line'),'linew',5)
hold on

plot(X,Y,'.','MarkerSize',80,'Color',[.5 .5 .5]);
hold on
plot(mean(data),'.','MarkerSize',80,'Color','red');
hold on
yline(0.5,'LineStyle','--','LineWidth',5);
hold on

[fitobj,goodness] = fit(X,Y,eq,options);
mdl = fitlm(X,Y);

plot(X,fitobj(X),'LineWidth',5,'Color',[0 0 0]);
hold on

yticks([]);
ylabel(['Rightward               ' '               Leftward']);

%==========================================================================
%-Create RDM
%==========================================================================
for p = 1:numel(subjs)
    for i = 1:numel(params.num.uval)
        for j = 1:numel(params.num.uval)
            % with outlier
            v_h = result.withOutlier{p,1}{1};
            v_v = result.withOutlier{p,2}{1};
            rdm.withOutlier.horizontal(p,i,j) = abs(v_h(i)-v_h(j));
            rdm.withOutlier.vertical(p,i,j) = abs(v_v(i)-v_v(j));
            % without outlier
            v_h = result.withoutOutlier{p,1}{1};
            v_v = result.withoutOutlier{p,2}{1};
            rdm.withoutOutlier.horizontal(p,i,j) = abs(v_h(i)-v_h(j));
            rdm.withoutOutlier.vertical(p,i,j) = abs(v_v(i)-v_v(j));
        end
    end
end

rdm.withOutlier.horizontalMean = squeeze(mean(rdm.withOutlier.horizontal(:,:,:),1));
rdm.withOutlier.verticalMean = squeeze(mean(rdm.withOutlier.vertical(:,:,:),1));
rdm.withoutOutlier.horizontalMean = squeeze(mean(rdm.withoutOutlier.horizontal(:,:,:),1));
rdm.withoutOutlier.verticalMean = squeeze(mean(rdm.withoutOutlier.vertical(:,:,:),1));

xy_ds = cmdscale(rdm.withOutlier.horizontalMean);
plot(xy_ds(1:end,1),xy_ds(1:end,2),'.','MarkerSize',10);
text(xy_ds(1,1),xy_ds(1,2),'Six');
text(xy_ds(2,1),xy_ds(2,2),'Ten');
text(xy_ds(3,1),xy_ds(3,2),'Seventeen');
text(xy_ds(4,1),xy_ds(4,2),'Twenty-nine');

figure;
xy_ds = cmdscale(rdm.withoutOutlier.horizontalMean);
plot(xy_ds(1:end,1),xy_ds(1:end,3),'.','MarkerSize',10);
text(xy_ds(1,1),xy_ds(1,3),'Six');
text(xy_ds(2,1),xy_ds(2,3),'Ten');
text(xy_ds(3,1),xy_ds(3,3),'Seventeen');
text(xy_ds(4,1),xy_ds(4,3),'Twenty-nine');

figure;
xy_ds = cmdscale(rdm.withOutlier.verticalMean);
plot(xy_ds(1:end,1),xy_ds(1:end,3),'.','MarkerSize',10);
text(xy_ds(1,1),xy_ds(1,3),'Six');
text(xy_ds(2,1),xy_ds(2,3),'Ten');
text(xy_ds(3,1),xy_ds(3,3),'Seventeen');
text(xy_ds(4,1),xy_ds(4,3),'Twenty-nine');

figure;
xy_ds = cmdscale(rdm.withoutOutlier.verticalMean);
plot(xy_ds(1:end,1),xy_ds(1:end,3),'.','MarkerSize',10);
text(xy_ds(1,1),xy_ds(1,3),'Six');
text(xy_ds(2,1),xy_ds(2,3),'Ten');
text(xy_ds(3,1),xy_ds(3,3),'Seventeen');
text(xy_ds(4,1),xy_ds(4,3),'Twenty-nine');