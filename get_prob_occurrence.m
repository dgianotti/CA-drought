function p_occ = get_prob_occurrence(occ_history,trans_probs)
% The function get_prob_occurrence gives the probability of occurrence on a
% vector of observations (all from the same day of year). 
%
% The input occ_states is the "history" leading up to that day, an N x 5 
% matrix of occurrence. [0 0 0 1 1] means that it rained the last two days,
% but not the three before that.
%
% The input trans_probs is a vector of transition probabilities
% (unknown length) directly from the AICOccurModel (see AICOccurModel.m)
% for the station. I'll have to see how those were stored...
%
% The output, p_occ, is the probability of occurrence for each observation
% in the vector, given its history.

% Okay, so here's the key to the trans_probs:
% If the length of the trans_probs is 1, then it's a 0-order model, and the
% trans_probs is the probability of rain.
% If the length is 2, then it's a 1st order model, and the first value is
% the probability of rain given NO RAIN (0) yesterday, and the second is
% probability of rain given RAIN (1) yesterday.
% If the length is 4, then 2nd order, {p001, p011, p101, p111}, where the
% last index value (always 1) means probability of rain TODAY.

% We're going to look at this as seperate cases for the length of
% trans_probs vector (ie the chain order).

N = size(occ_history,1);
p_occ = nan(N,1);

switch length(trans_probs)
    case 1 % zero-order chain
        % history doesn't matter, p_occ is trans_probs!
        p_occ(:) = trans_probs;
    case 2 % 1st order chain, [p01; p11]
        p_occ( occ_history(:,5) == 0 ) = trans_probs(1);
        p_occ( occ_history(:,5) == 1 ) = trans_probs(2);
    case 4 % 2nd order chain, [p001; p011; p101; p111]
        % make matching vectors:
        match_xxx00 = all( (occ_history(:,4:5) == repmat([0 0], [N,1])), 2 );
        p_occ(match_xxx00) = trans_probs(1);
        
        match_xxx01 = all( (occ_history(:,4:5) == repmat([0 1], [N,1])), 2 );
        p_occ(match_xxx01) = trans_probs(2);
        
        match_xxx10 = all( (occ_history(:,4:5) == repmat([1 0], [N,1])), 2 );
        p_occ(match_xxx10) = trans_probs(3);
        
        match_xxx11 = all( (occ_history(:,4:5) == repmat([1 1], [N,1])), 2 );
        p_occ(match_xxx11) = trans_probs(4);
        
    case 8 % third order chain, [p0001; p0011; p0101; p0111; p1001; p1011; p1101; p1111]
        % make matching vectors:
        match_xx000 = all( (occ_history(:,3:5) == repmat([0 0 0], [N,1])), 2 );
        p_occ(match_xx000) = trans_probs(1);
        
        match_xx001 = all( (occ_history(:,3:5) == repmat([0 0 1], [N,1])), 2 );
        p_occ(match_xx001) = trans_probs(2);

        match_xx010 = all( (occ_history(:,3:5) == repmat([0 1 0], [N,1])), 2 );
        p_occ(match_xx010) = trans_probs(3);
        
        match_xx011 = all( (occ_history(:,3:5) == repmat([0 1 1], [N,1])), 2 );
        p_occ(match_xx011) = trans_probs(4);

        match_xx100 = all( (occ_history(:,3:5) == repmat([1 0 0], [N,1])), 2 );
        p_occ(match_xx100) = trans_probs(5);

        match_xx101 = all( (occ_history(:,3:5) == repmat([1 0 1], [N,1])), 2 );
        p_occ(match_xx101) = trans_probs(6);

        match_xx110 = all( (occ_history(:,3:5) == repmat([1 1 0], [N,1])), 2 );
        p_occ(match_xx110) = trans_probs(7);

        match_xx111 = all( (occ_history(:,3:5) == repmat([1 1 1], [N,1])), 2 );
        p_occ(match_xx111) = trans_probs(8);

    case 16 % fourth order chain, 
        % [p00001; p00011; p00101; p00111; p01001; p01011; p01101; p01111;
        %  p10001; p10011; p10101; p10111; p11001; p11011; p11101; p11111]
        
        % make matching vectors:
        match_x0000 = all( (occ_history(:,2:5) == repmat([0 0 0 0], [N,1])), 2 );
        p_occ(match_x0000) = trans_probs(1);
   
        match_x0001 = all( (occ_history(:,2:5) == repmat([0 0 0 1], [N,1])), 2 );
        p_occ(match_x0001) = trans_probs(2);

        match_x0010 = all( (occ_history(:,2:5) == repmat([0 0 1 0], [N,1])), 2 );
        p_occ(match_x0010) = trans_probs(3);
        
        match_x0011 = all( (occ_history(:,2:5) == repmat([0 0 1 1], [N,1])), 2 );
        p_occ(match_x0011) = trans_probs(4);
        
        match_x0100 = all( (occ_history(:,2:5) == repmat([0 1 0 0], [N,1])), 2 );
        p_occ(match_x0100) = trans_probs(5);
        
        match_x0101 = all( (occ_history(:,2:5) == repmat([0 1 0 1], [N,1])), 2 );
        p_occ(match_x0101) = trans_probs(6);
        
        match_x0110 = all( (occ_history(:,2:5) == repmat([0 1 1 0], [N,1])), 2 );
        p_occ(match_x0110) = trans_probs(7);
        
        match_x0111 = all( (occ_history(:,2:5) == repmat([0 1 1 1], [N,1])), 2 );
        p_occ(match_x0111) = trans_probs(8);
        
        match_x1000 = all( (occ_history(:,2:5) == repmat([1 0 0 0], [N,1])), 2 );
        p_occ(match_x1000) = trans_probs(9);
        
        match_x1001 = all( (occ_history(:,2:5) == repmat([1 0 0 1], [N,1])), 2 );
        p_occ(match_x1001) = trans_probs(10);
        
        match_x1010 = all( (occ_history(:,2:5) == repmat([1 0 1 0], [N,1])), 2 );
        p_occ(match_x1010) = trans_probs(11);
        
        match_x1011 = all( (occ_history(:,2:5) == repmat([1 0 1 1], [N,1])), 2 );
        p_occ(match_x1011) = trans_probs(12);
        
        match_x1100 = all( (occ_history(:,2:5) == repmat([1 1 0 0], [N,1])), 2 );
        p_occ(match_x1100) = trans_probs(13);
        
        match_x1101 = all( (occ_history(:,2:5) == repmat([1 1 0 1], [N,1])), 2 );
        p_occ(match_x1101) = trans_probs(14);
        
        match_x1110 = all( (occ_history(:,2:5) == repmat([1 1 1 0], [N,1])), 2 );
        p_occ(match_x1110) = trans_probs(15);
        
        match_x1111 = all( (occ_history(:,2:5) == repmat([1 1 1 1], [N,1])), 2 );
        p_occ(match_x1111) = trans_probs(16);

    case 32 % fifth order chain! Yikes! 
        % [p000001; p000011; p000101; p000111; p001001; p001011; p001101; p001111;
        %  p010001; p010011; p010101; p010111; p011001; p011011; p011101; p011111;
        %  p100001; p100011; p100101; p100111; p101001; p101011; p101101; p101111;
        %  p110001; p110011; p110101; p110111; p111001; p111011; p111101; p111111]
        
        % make matching vectors:
        match_00000 = all( (occ_history(:,1:5) == repmat([0 0 0 0 0], [N,1])), 2 );
        p_occ(match_00000) = trans_probs(1);
  
        match_00001 = all( (occ_history(:,1:5) == repmat([0 0 0 0 1], [N,1])), 2 );
        p_occ(match_00001) = trans_probs(2);
  
        match_00010 = all( (occ_history(:,1:5) == repmat([0 0 0 1 0], [N,1])), 2 );
        p_occ(match_00010) = trans_probs(3);
  
        match_00011 = all( (occ_history(:,1:5) == repmat([0 0 0 1 1], [N,1])), 2 );
        p_occ(match_00011) = trans_probs(4);
  
        match_00100 = all( (occ_history(:,1:5) == repmat([0 0 1 0 0], [N,1])), 2 );
        p_occ(match_00100) = trans_probs(5);
  
        match_00101 = all( (occ_history(:,1:5) == repmat([0 0 1 0 1], [N,1])), 2 );
        p_occ(match_00101) = trans_probs(6);
  
        match_00110 = all( (occ_history(:,1:5) == repmat([0 0 1 1 0], [N,1])), 2 );
        p_occ(match_00110) = trans_probs(7);
  
        match_00111 = all( (occ_history(:,1:5) == repmat([0 0 1 1 1], [N,1])), 2 );
        p_occ(match_00111) = trans_probs(8);
  
        match_01000 = all( (occ_history(:,1:5) == repmat([0 1 0 0 0], [N,1])), 2 );
        p_occ(match_01000) = trans_probs(9);
  
        match_01001 = all( (occ_history(:,1:5) == repmat([0 1 0 0 1], [N,1])), 2 );
        p_occ(match_01001) = trans_probs(10);
  
        match_01010 = all( (occ_history(:,1:5) == repmat([0 1 0 1 0], [N,1])), 2 );
        p_occ(match_01010) = trans_probs(11);
  
        match_01011 = all( (occ_history(:,1:5) == repmat([0 1 0 1 1], [N,1])), 2 );
        p_occ(match_01011) = trans_probs(12);
  
        match_01100 = all( (occ_history(:,1:5) == repmat([0 1 1 0 0], [N,1])), 2 );
        p_occ(match_01100) = trans_probs(13);
  
        match_01101 = all( (occ_history(:,1:5) == repmat([0 1 1 0 1], [N,1])), 2 );
        p_occ(match_01101) = trans_probs(14);
  
        match_01110 = all( (occ_history(:,1:5) == repmat([0 1 1 1 0], [N,1])), 2 );
        p_occ(match_01110) = trans_probs(15);
  
        match_01111 = all( (occ_history(:,1:5) == repmat([0 1 1 1 1], [N,1])), 2 );
        p_occ(match_01111) = trans_probs(16);
  
        match_10000 = all( (occ_history(:,1:5) == repmat([1 0 0 0 0], [N,1])), 2 );
        p_occ(match_10000) = trans_probs(17);
  
        match_10001 = all( (occ_history(:,1:5) == repmat([1 0 0 0 1], [N,1])), 2 );
        p_occ(match_10001) = trans_probs(18);
  
        match_10010 = all( (occ_history(:,1:5) == repmat([1 0 0 1 0], [N,1])), 2 );
        p_occ(match_10010) = trans_probs(19);
  
        match_10011 = all( (occ_history(:,1:5) == repmat([1 0 0 1 1], [N,1])), 2 );
        p_occ(match_10011) = trans_probs(20);
  
        match_10100 = all( (occ_history(:,1:5) == repmat([1 0 1 0 0], [N,1])), 2 );
        p_occ(match_10100) = trans_probs(21);
  
        match_10101 = all( (occ_history(:,1:5) == repmat([1 0 1 0 1], [N,1])), 2 );
        p_occ(match_10101) = trans_probs(22);
  
        match_10110 = all( (occ_history(:,1:5) == repmat([1 0 1 1 0], [N,1])), 2 );
        p_occ(match_10110) = trans_probs(23);
  
        match_10111 = all( (occ_history(:,1:5) == repmat([1 0 1 1 1], [N,1])), 2 );
        p_occ(match_10111) = trans_probs(24);
  
        match_11000 = all( (occ_history(:,1:5) == repmat([1 1 0 0 0], [N,1])), 2 );
        p_occ(match_11000) = trans_probs(25);
  
        match_11001 = all( (occ_history(:,1:5) == repmat([1 1 0 0 1], [N,1])), 2 );
        p_occ(match_11001) = trans_probs(26);
  
        match_11010 = all( (occ_history(:,1:5) == repmat([1 1 0 1 0], [N,1])), 2 );
        p_occ(match_11010) = trans_probs(27);
  
        match_11011 = all( (occ_history(:,1:5) == repmat([1 1 0 1 1], [N,1])), 2 );
        p_occ(match_11011) = trans_probs(28);
  
        match_11100 = all( (occ_history(:,1:5) == repmat([1 1 1 0 0], [N,1])), 2 );
        p_occ(match_11100) = trans_probs(29);
  
        match_11101 = all( (occ_history(:,1:5) == repmat([1 1 1 0 1], [N,1])), 2 );
        p_occ(match_11101) = trans_probs(30);
  
        match_11110 = all( (occ_history(:,1:5) == repmat([1 1 1 1 0], [N,1])), 2 );
        p_occ(match_11110) = trans_probs(31);
  
        match_11111 = all( (occ_history(:,1:5) == repmat([1 1 1 1 1], [N,1])), 2 );
        p_occ(match_11111) = trans_probs(32);

end % switch
        
        if any(isnan(p_occ))
            error('The p_occ values whould all be probabilities, but at least one is nan!  Error!\n');
        elseif any( p_occ > 1 )
            error('The p_occ values should all be probabilities, but at least one is greater than 1! Error! \n');
        elseif any(p_occ < 0 )
            error('The p_occ values should all be probabilities, but at least one is less than 0! Error! \n');
        end