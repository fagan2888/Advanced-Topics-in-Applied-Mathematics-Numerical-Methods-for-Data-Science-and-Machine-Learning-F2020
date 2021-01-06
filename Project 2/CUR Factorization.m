%Exercise 3 CUR factorization
clear all
[M,y] = readdata();
[nrow , ncol] = size(M);
fro_norm = zeros(9,8);
rel_err = zeros(9,8);
%Set alpha
for alpha = 1:8
for k = 2:10
    fr = zeros(100,1);
    rl = zeros(100,1);
    for iter = 1:100
    %Eckart-Young best approximation-Truncated SVD
    [U, S, V] = svd(M,'econ');
    Uk = U(:,1:k);
    Vtr = V';
    Vk = Vtr(1:k,:);
    Sk = S(1:k,1:k);
    X0 = Uk*Sk;
    Y0 = Vk;
    Mk = X0*Y0;
    %Frobenius norm
    %fro_norm(k-1,alpha) = norm(M-Mk,'fro');
    fr(iter) = norm(M-Mk,'fro');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    c = alpha*k;
    r = c;
    %column select
    C = ex3_columnselect(M,k,c);
    %row select
    R = ex3_columnselect(M',k,r);
    %matrix U
    U = pinv(C)*M*pinv(R)';
    %Relative error
    %rel_err(k-1,alpha) = norm(M-C*U*R','fro')/norm(M-Mk,'fro');
    rl(iter) = norm(M-C*U*R','fro')/norm(M-Mk,'fro');
    end
    fro_norm(k-1,alpha) = sum(fr)/100;
    rel_err(k-1,alpha) = sum(rl)/100;
end

end

hold on
for alpha = 1:8
    
    %plot(2:10,fro_norm(1:9,alpha))
    plot(2:10,rel_err(1:9,alpha))
end
%xlabel("k")
%ylabel("Frobenius norm")
xlabel("k")
ylabel("Relative error")
legend('alpha = 1','alpha = 2', 'alpha = 3', 'alpha = 4','alpha = 6','alpha = 6', 'alpha = 7', 'alpha = 8');

%For k = 10 and alpha = 6 we get the smallest relative error
[m n] = min(rel_err(:));
[x y z] = ind2sub(size(rel_err),n);
%After inspection a reasonable pair of values is k = 7 and alpha = 6 