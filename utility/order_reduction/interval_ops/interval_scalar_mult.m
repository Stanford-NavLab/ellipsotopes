function I = interval_scalar_mult(s,I)
    I = s.*I ;
    if s < 0
        I = [I(:,2), I(:,1)] ;
    end
end