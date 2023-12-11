%__________________________________________________________________________
function [LDotDets,SDotDets] = MakeDotMatch(ScreenDets,StimDets,ExpDets)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Make match stimuli for a run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-Create list of dots diameters for larger match within each sample stimulus 
%--------------------------------------------------------------------------
DotDets = {};

for i=1:ExpDets.numNum              % 3 numbers
    for ii=1:ExpDets.numSize        % 3 sizes
        NumDots=StimDets.Num.LMatch(i);
        DotDiamAVG=StimDets.DotDiam.AVG.vd(ii);
        DotDiamSTD=StimDets.DotDiam.STD.vd(ii);
        DiamList{i,ii}=randDiamPick(DotDiamAVG,DotDiamSTD,NumDots);
    end
end

%-Generate larger match stimuli specification and save it
%--------------------------------------------------------------------------
for i=1:ExpDets.numTFA                      % iterate on TFAs
    for ii=1:ExpDets.numNum                 % iterate on the number of dots
        for iii=1:ExpDets.numSize           % iterate on the sizes
            
            [DotPosition,DotDiam] = MakeDot(100*StimDets.TFDiam.vd(i),StimDets.Num.LMatch(ii),...
                100*DiamList{ii,iii},100*StimDets.BuffMod.vd,StimDets.ImShape);              % determines dots' position
            
            DotPosition = DotPosition/100;
            DotDiam = DotDiam/100;
            
            %-i: TFA, ii: Number of dots, iii: Average diamater of dots
            index.tfa=i; index.num=ii; index.size=iii;
            
            %-values in vis deg
            DotDets.DotPosition.vd{index.num,index.size,index.tfa} = DotPosition;
            DotDets.DiamList.vd{index.num,index.size,index.tfa} = DotDiam;
            DotDets.TFDiam.vd{index.num,index.size,index.tfa} = StimDets.TFDiam.vd(i);
            DotDets.Density.vd{index.num,index.size,index.tfa} = ...
                StimDets.Num.LMatch(ii)/StimDets.TFA.vd(i);
            DotDets.TSA.vd{index.num,index.size,index.tfa} = ...
                StimDets.Num.LMatch(ii)*StimDets.DotArea.AVG.vd(iii);
            
        end
    end
end

%-larger match dots details
LDotDets=DotDets;

%-Create list of dots diameters for smaller match within each sample stimulus
%--------------------------------------------------------------------------
DotDets = {};

for i=1:ExpDets.numNum              % 3 numbers
    for ii=1:ExpDets.numSize        % 3 sizes
        NumDots=StimDets.Num.SMatch(i);
        DotDiamAVG=StimDets.DotDiam.AVG.vd(ii);
        DotDiamSTD=StimDets.DotDiam.STD.vd(ii);
        DiamList{i,ii}=randDiamPick(DotDiamAVG,DotDiamSTD,NumDots);
    end
end

%-Generate smaller match stimuli specification and save it
%--------------------------------------------------------------------------
for i=1:ExpDets.numTFA                      % iterate on TFAs
    for ii=1:ExpDets.numNum                 % iterate on the number of dots
        for iii=1:ExpDets.numSize           % iterate on the sizes
            
            %-determines dots' position
            [DotPosition,DotDiam] = MakeDot(100*StimDets.TFDiam.vd(i),StimDets.Num.SMatch(ii),...
                100*DiamList{ii,iii},100*StimDets.BuffMod.vd,StimDets.ImShape);
            
            DotPosition = DotPosition/100;
            DotDiam = DotDiam/100;
            
            %-i: TFA, ii: Number of dots, iii: Average diamater of dots
            index.tfa=i; index.num=ii; index.size=iii;
            
            %-values in vis deg
            DotDets.DotPosition.vd{index.num,index.size,index.tfa} = DotPosition;
            DotDets.DiamList.vd{index.num,index.size,index.tfa} = DotDiam;
            DotDets.TFDiam.vd{index.num,index.size,index.tfa} = StimDets.TFDiam.vd(i);
            DotDets.Density.vd{index.num,index.size,index.tfa} = ...
                StimDets.Num.SMatch(ii)/StimDets.TFA.vd(i);
            DotDets.TSA.vd{index.num,index.size,index.tfa} = ...
                StimDets.Num.SMatch(ii)*StimDets.DotArea.AVG.vd(iii);
            
        end
    end
end

%-smaller match dots details
SDotDets=DotDets;

end
