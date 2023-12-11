%__________________________________________________________________________
function inx3d = convInx1d3d(inx1d,ExpDets)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Convert a number to three index (Number, Size, TFA)
%-Convert a 1D index to 3D index. The firs index is number index,
%the second index is size index and the third index is TFA index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inx=0;
for i=1:ExpDets.numNum
   for ii=1:ExpDets.numSize
      for iii=1:ExpDets.numTFA
          inx=inx+1;
          if inx1d==inx
              %i: number index, ii: size index, iii: tfa index
              inx3d={[i,ii,iii]};
              break;
          end
      end
   end
end

end
