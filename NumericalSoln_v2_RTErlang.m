function RTErlang = NumericalSoln_v2_RTErlang(x,lambda,k,U)
syms w n m

CDF = double((1 - symsum((1/factorial(n))*exp(-lambda*U)*(lambda*U)^n,n,0,k-1)));
logErlang = k*log(lambda) + (k-1)*log(x) - lambda*x - double(symsum(log(w),w,2,k-1));
RTErlang = exp(logErlang)./CDF;
end
