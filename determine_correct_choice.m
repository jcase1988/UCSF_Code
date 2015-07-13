function correct_choice = determine_correct_choice(offered,choice,correct_bool)
    %outputs the stim ID of the correct sound
    
    for iTrial = 1:length(correct_bool)
        if correct_bool(iTrial)
            correct_choice(iTrial,1) = offered(iTrial,choice(iTrial));
        else
            correct_choice(iTrial,1) = offered(iTrial,~(choice(iTrial)-1)+1);
        end
    end
    
end
