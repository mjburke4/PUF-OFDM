K = 200;
tp1 = zeros(1,K); ti1 = zeros(1,K);
tp0 = zeros(1,K); ti0 = zeros(1,K);

for i=1:K
    o1 = test_one_symbol_smoke(10,true);
    o0 = test_one_symbol_smoke(10,false);
    tp1(i)=o1.Tperm_abs; ti1(i)=o1.Tim.norm;
    tp0(i)=o0.Tperm_abs; ti0(i)=o0.Tim.norm;
end

mean(tp1), mean(tp0)
mean(ti1), mean(ti0)
