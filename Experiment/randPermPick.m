%__________________________________________________________________________
function [picks] = randPermPick(rng,len)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-Create a random sequence of length=len and range=rng in which almost
%all numbers distributed evenly
%-This function is used for determining shuffling list of dots size,
%and generation of sample ISI sequence and match sequence

%-rng: range of number
%-len: length of the list of number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if len <= rng
        pick = randperm(rng);
        picks = pick(1:len);
    else
        picks_temp = [];
        bin = len/rng;
        for i = 1:floor(bin)
            pick = randperm(rng);
            picks_temp = [picks_temp ; pick'];
        end
    
        pick = randperm(rng);
        pick_temp = pick(1:(len-((floor(bin))*rng)));
        picks_temp = [picks_temp ; (pick_temp)'];
    
        picks = picks_temp';
    end
    
    %-This extra line was added to the code by me to shuffle the picks
    %inorder to mix each bin with each other
    %https://www.mathworks.com/matlabcentral/answers/317259-how-to-randomise-numbers-in-a-vector
    picks = picks(randperm(numel(picks)));
    
end
