%__________________________________________________________________________
function DotDetsOut = convVdPx(modality,DotDetsIn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Convert dot details from vis deg to pixel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(modality,'meg')
[SubjDets,ScreenDets,ExpDets,StimDets] = params(39,[20 15],[1440 1080]);
    
for i=1:ExpDets.numTFA                  % iterate on TFAs
    for ii=1:ExpDets.numNum             % iterate on the number of dots
        for iii=1:ExpDets.numSize       % iterate on the sizes
            %-i: TFA, ii: Number of dots, iii: Average diamater of dots
            index.tfa=i; index.num=ii; index.size=iii;
            
            %-values in pixel
            DotDetsIn.DotPosition.px{index.num,index.size,index.tfa} = ...
                DotDetsIn.DotPosition.vd{index.num,index.size,index.tfa}/ScreenDets.degperpix;
            DotDetsIn.DiamList.px{index.num,index.size,index.tfa} = ...
                DotDetsIn.DiamList.vd{index.num,index.size,index.tfa}/ScreenDets.degperpix;
            DotDetsIn.TFDiam.px{index.num,index.size,index.tfa} = StimDets.TFDiam.px(i);
            DotDetsIn.Density.px{index.num,index.size,index.tfa} = ...
                StimDets.Num.Sample(ii)/StimDets.TFA.px(i);
            DotDetsIn.TSA.px{index.num,index.size,index.tfa} = ...
                StimDets.Num.Sample(ii)*StimDets.DotArea.AVG.px(iii);
        end
    end
end

elseif strcmp(modality,'fmri')
[SubjDets,ScreenDets,ExpDets,StimDets] = params(65,[34.9 19.6],[1920 1080]);

for i=1:ExpDets.numTFA                  % iterate on TFAs
    for ii=1:ExpDets.numNum             % iterate on the number of dots
        for iii=1:ExpDets.numSize       % iterate on the sizes
            %-i: TFA, ii: Number of dots, iii: Average diamater of dots
            index.tfa=i; index.num=ii; index.size=iii;
            
            %-values in pixel
            DotDetsIn.DotPosition.px{index.num,index.size,index.tfa} = ...
                DotDetsIn.DotPosition.vd{index.num,index.size,index.tfa}/ScreenDets.degperpix;
            DotDetsIn.DiamList.px{index.num,index.size,index.tfa} = ...
                DotDetsIn.DiamList.vd{index.num,index.size,index.tfa}/ScreenDets.degperpix;
            DotDetsIn.TFDiam.px{index.num,index.size,index.tfa} = StimDets.TFDiam.px(i);
            DotDetsIn.Density.px{index.num,index.size,index.tfa} = ...
                StimDets.Num.Sample(ii)/StimDets.TFA.px(i);
            DotDetsIn.TSA.px{index.num,index.size,index.tfa} = ...
                StimDets.Num.Sample(ii)*StimDets.DotArea.AVG.px(iii);
        end
    end
end

else
    error('Invalid input!');
end

DotDetsOut = DotDetsIn;

end
