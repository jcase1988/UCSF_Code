function data = notf(data)
    
    data.data = notchfilter(data.data,data.sampFreq,60);
    data.data = notchfilter(data.data,data.sampFreq,120);
    data.data = notchfilter(data.data,data.sampFreq,180);
end