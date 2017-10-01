%% Matrix fractional norm implementation 
function f = fNorm(X, Frac)
    elemFrac = X.^Frac;
    upperLimit = 100;
    if Frac < 0
        elemFrac(find(elemFrac>upperLimit)) = upperLimit;
    end
    f = sum(sum(elemFrac));
end