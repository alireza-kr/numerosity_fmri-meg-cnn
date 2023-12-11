%__________________________________________________________________________
function DotDets = MakeDotSample(ScreenDets,StimDets,ExpDets)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Make sample stimuli for a run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-Create list of dots diameters within each sample stimulus
%--------------------------------------------------------------------------
DotDets = {};

for i=1:ExpDets.numNum                  % 4 numbers
    for ii=1:ExpDets.numSize            % 4 sizes
        NumDots=StimDets.Num.Sample(i);
        DotDiamAVG=StimDets.DotDiam.AVG.vd(ii);
        DotDiamSTD=StimDets.DotDiam.STD.vd(ii);
        DiamList{i,ii}=randDiamPick(DotDiamAVG,DotDiamSTD,NumDots);
    end
end

%-Generate sample stimuli specification and save it
%--------------------------------------------------------------------------
for i=1:ExpDets.numTFA                  % iterate on TFAs
    for ii=1:ExpDets.numNum             % iterate on the number of dots
        for iii=1:ExpDets.numSize       % iterate on the sizes
            
            %-determines dots' position
            [DotPosition,DotDiam] = MakeDot(100*StimDets.TFDiam.vd(i),StimDets.Num.Sample(ii),...
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
                StimDets.Num.Sample(ii)/StimDets.TFA.vd(i);
            DotDets.TSA.vd{index.num,index.size,index.tfa} = ...
                StimDets.Num.Sample(ii)*StimDets.DotArea.AVG.vd(iii);
            
        end
    end
end

end
