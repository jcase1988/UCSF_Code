function data = hgf(data)

    %notch filters
    data.data = notchfilter(data.data,data.sampFreq,120);
   
    for iElec = 1:size(data.data,1)

        % High gamma filter
        data.data(iElec,:) = abs(my_hilbert(data.data(iElec,:), data.sampFreq, 70, 150)).^2;
    end
end

