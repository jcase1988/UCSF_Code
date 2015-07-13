function correct_bool = determine_correct_bool(choice,offered)
for i=1:length(choice)
    c = choice(i);
    d = (~(choice(i)-1))+1;
    if values(offered(i,c)) > values(offered(i,d))
        correct_bool(i,1) = 1;
    else
        correct_bool(i,1) = 0;
    end
end
end