%% This script generates the information of all stimuli for a subject and save it a mat file.

function [FMRI,MEG] = MakeExp

% https://www.mathworks.com/matlabcentral/answers/140748-explain-this-one-line-of-code
rand('state',sum(clock*100));

%% Load parameters
% load experiment parameters
% create stimuli based on my notebook's monitor specification
[SubjDets,ScreenDets,ExpDets,StimDets] = params(0,0,0);

%% Make dot details
% generate dot details for both fMRI and MEG experiments
% the same stimuli have been used in both fMRI and MEG experiments
% if number of experimental blocks are different in the fMRI and MEG experiment,
% two different loops must be used to generate dots detail

for i=1:ExpDets.numExpBlock.all
    DotDets.Sample=MakeDotSample(ScreenDets,StimDets,ExpDets);
    [DotDets.Match.Larger,DotDets.Match.Smaller]=MakeDotMatch(ScreenDets,StimDets,ExpDets);
    % convert dots' information from vis deg to pixel for MEG and fMRI experiemtn
    DotDetsFMRI.Sample=convVdPx('fmri',DotDets.Sample);
    DotDetsFMRI.Match.Larger=convVdPx('fmri',DotDets.Match.Larger);
    DotDetsFMRI.Match.Smaller=convVdPx('fmri',DotDets.Match.Smaller);
    DotDetsMEG.Sample=convVdPx('meg',DotDets.Sample);
    DotDetsMEG.Match.Larger=convVdPx('meg',DotDets.Match.Larger);
    DotDetsMEG.Match.Smaller=convVdPx('meg',DotDets.Match.Smaller);
    % save dots' details for fMRI and MEG experiment
    FMRI.ExpBlock(i).DotDets=DotDetsFMRI;
    MEG.ExpBlock(i).DotDets=DotDetsMEG;
end

%% Make sequence detail
% although the same stimuli have been used in both fMRI and MEG experiment,
% the order of stimuli is different in them

% generate the sequence for the entire session of the fMRI experiment
sequenceFMRI=MakeSeqFMRI(ExpDets);

% generate the sequence for the entire sessions of the MEG experiment
for sess=1:ExpDets.numSession.meg
    sequenceMEG{sess}=MakeSeqMEG(ExpDets);
end

% generate sequence detail for each experimental block (fMRI experiment)
for i=1:ExpDets.numExpBlock.fmri
    FMRI.ExpBlock(i).SeqDets=sequenceFMRI(i);
end

% generate sequence detail for each experimental block (MEG experiment)
for sess=1:ExpDets.numSession.meg
    for i=1:ExpDets.numExpBlock.meg
        MEG.ExpBlock(i).SeqDets{sess}=sequenceMEG{sess}(i);
    end
end
