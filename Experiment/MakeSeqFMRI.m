%__________________________________________________________________________
function SeqDets = MakeSeqFMRI(ExpDets)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Make stimuli and ISI sequence (sample and match) for a session

%-SeqStim1D: stimulus number (1D)
%-SeqStim3D: stimulus number (3D)
%-SeqStimISI: ISI, SeqStimType: Sample, Match-Larger, Match-Smaller
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-Create SeqStim1D, SeqStim3D, SeqStimISI, and SeqStimType
%--------------------------------------------------------------------------
%-generate list of index for smaller/larger match stimuli for a run
%-the smaller/larger match stimuli balanced across a run
for i=1:ExpDets.numRun.fmri
    SeqMatchType.Run{i}=randPermPick(2,ExpDets.numBlock.fmri*ExpDets.numMatch);
end

%-generate list of index for smaller/larger match stimuli for each
%experimental block
expBlockIter=1;
for i=1:ExpDets.numRun.fmri
   for ii=1:ExpDets.numBlock.fmri
       SeqMatchType.Block{expBlockIter}=SeqMatchType.Run{i}...
           (((ii-1)*ExpDets.numMatch)+1:ii*ExpDets.numMatch);
       expBlockIter=expBlockIter+1;
   end
end

for i=1:ExpDets.numExpBlock.fmri
    %-create sequence of index for sample stimlui (1D) in an experimental block
    SeqSample1D = randperm(ExpDets.numSample);
    
    %-create sequence of index for sample ISI in an experimental block
    SeqSampleISI = ExpDets.ISI.Sample.fmri(randPermPick(length(ExpDets.ISI.Sample.fmri),ExpDets.numSample));
    
    %-create sequence of index for match ISI in an experimental block
    SeqMatchISI = ExpDets.ISI.Match.fmri(randPermPick(length(ExpDets.ISI.Match.fmri),ExpDets.numMatch));
    
    %-specifies where to place match stimuli in an experimental block (SeqMatchPlace)
    matchDist=true;
    while(matchDist)
        matchDist=false;
        SeqMatchPlace=datasample(2:ExpDets.numSample,ExpDets.numMatch,'Replace',false);
        SeqMatchPlace=sort(SeqMatchPlace);
        for j=1:(ExpDets.numMatch-1)
            if SeqMatchPlace(j+1)-SeqMatchPlace(j)<=ExpDets.matchDist.fmri
                matchDist=true;
            end
        end
    end

    %-create sequence of stimuli (sample+match, 1D) for a session
    %-create sequence of ISI (sample+match) for a session
    numMatchIter=ExpDets.numMatch;
    ExpBlock(i).SeqStim1D=SeqSample1D;
    ExpBlock(i).SeqStimISI=SeqSampleISI;
    ExpBlock(i).SeqStimType=zeros(1,ExpDets.numSample);
    while numMatchIter>0
       ExpBlock(i).SeqStim1D = [ExpBlock(i).SeqStim1D(1:SeqMatchPlace(numMatchIter)) ...
           ExpBlock(i).SeqStim1D(SeqMatchPlace(numMatchIter)) ...
           ExpBlock(i).SeqStim1D(SeqMatchPlace(numMatchIter)+1:length(ExpBlock(i).SeqStim1D))];
       ExpBlock(i).SeqStimISI = [ExpBlock(i).SeqStimISI(1:SeqMatchPlace(numMatchIter)) ...
           SeqMatchISI(numMatchIter) ...
           ExpBlock(i).SeqStimISI(SeqMatchPlace(numMatchIter)+1:length(ExpBlock(i).SeqStimISI))];
       ExpBlock(i).SeqStimType = [ExpBlock(i).SeqStimType(1:SeqMatchPlace(numMatchIter)) ...
           SeqMatchType.Block{i}(numMatchIter) ...
           ExpBlock(i).SeqStimType(SeqMatchPlace(numMatchIter)+1:length(ExpBlock(i).SeqStimType))];
       numMatchIter = numMatchIter-1;
    end
    
    %-create sequence of stimuli (sample+match, 3D) for a session
    for j=1:ExpDets.numStim
        ExpBlock(i).SeqStim3D(j) = convInx1d3d(ExpBlock(i).SeqStim1D(j),ExpDets);
    end
end

SeqDets=ExpBlock;
