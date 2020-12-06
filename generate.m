function y = generate(x,h,n)
    i = 1;
    y = 0;
    N = length(x);
    while i <= N/2
        y = y + quant(h(i)*(x(i)+x(N-i)),1/2^n);
        i = i + 1;
    end
end

