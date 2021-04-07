function coef = vec_nchoosek(n)
    k = 1:n ;
    coef = [1 cumprod((n-k+1)./k)] ;
end
    