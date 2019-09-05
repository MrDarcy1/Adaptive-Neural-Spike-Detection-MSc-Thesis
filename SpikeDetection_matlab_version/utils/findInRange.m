function out = findInRange(data,range)
    data=data(data<range(2));
    out=data(data>=range(1));
end

