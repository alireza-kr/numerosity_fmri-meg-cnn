%__________________________________________________________________________
function [picks] = randDiamPick(MeanDiam,StdDiam,NumOfDots)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Create a sequence of dots' diameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

picks = [MeanDiam-StdDiam:((MeanDiam+StdDiam)-(MeanDiam-StdDiam))/(NumOfDots-1):MeanDiam+StdDiam];

end
