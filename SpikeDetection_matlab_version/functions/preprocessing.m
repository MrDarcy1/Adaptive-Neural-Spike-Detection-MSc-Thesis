function out = preprocessing(data,Mode,par)
    switch Mode
        case 'aso'
            out=ASO(data,par{1},par{2});
        case 'neo'
            out=NEO(data,par{1},par{2});
    end
end

