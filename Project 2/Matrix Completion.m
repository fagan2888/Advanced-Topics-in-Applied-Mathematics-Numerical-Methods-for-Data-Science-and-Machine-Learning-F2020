%Exercise 2 Matrix Completion
clear all
MovieRankings36 = readtable('MovieRankings36.csv');
MR36 = MovieRankings36{:,:};
[nrow , ncol] = size(MR36);
MRproj = MR36;
MRproj(isnan(MRproj)) = 0;
%(a) Low-rank factorization A=XY'
[U,S,V] = svd(MRproj,'econ');
plot(1:20,diag(S))
legend('eigenvalues of MR36')
%Choose k by looking the scree plot
k = 3;
%Initialization step
lambda = 0;
Uk = U(:,1:k);
Vk = V(1:k,:);
Sk = S(1:k,1:k);
sl = functions_ex2(Sk,lambda,k);
X0 = Uk*sqrt(sl);
Y0 = Vk'*sqrt(sl);
Y0 = Y0';
%Alternating iteration
X = zeros(nrow,k);
Y = zeros(k,ncol);
N = 50;
obj_fnc = zeros(N+1,1);
%fro_norm = zeros(N+1,1);
obj_fnc(1) = functions1_ex2(X0,Y0', MRproj,lambda);
%fro_norm(1) = norm(MR36-X0*Y0,'fro');

%find(MR36(1,:)>0);
for iter = 1:N
    %LLS for matrix X
 alpha = zeros(nrow,ncol);

    for i = 1:nrow
        [nn, nm] = size(find(MRproj(i,:)>0));
        alpha(i,1:nm) = MR36(i,MRproj(i,:)>0);
     Y0_omega = zeros(k,nm);
    for j = 1:k
        Y0_omega(j,:) = Y0(j,MRproj(i,:)>0);
    end
    X(i,:) = pinv(Y0_omega*(Y0_omega') + lambda*eye(k))*Y0_omega*(alpha(i,1:nm))';
    end



%update
%obj_fnc = functions1_ex2(X,Y0, MR36,lambda);


%LLS for matrix Y
    beta = zeros(nrow,ncol);
    
for i = 1:ncol
    [nn, nm] = size(find(MRproj(:,i)>0));
    beta(1:nn,i) = MRproj(MRproj(:,i)>0,i);
X0_omega = zeros(nn,k);
    for j = 1:k
        X0_omega(:,j) = X(MRproj(:,i)>0,j);
    end
    Y(:,i) = pinv((X0_omega')*(X0_omega) + lambda*eye(k))*(X0_omega')*(beta(1:nn,i));
end

%update
obj_fnc(iter+1) = functions1_ex2(X,Y', MRproj,lambda);
Y0 = Y;
%fro_norm(iter+1) = norm(MR36-X*Y,'fro');
end

%set(gca,'Yscale','log')
plot(1:(N+1),obj_fnc)
title('objective function')
%hold on
 %hold off
legend('? = 0', ' ? = 1','? = 10','? = 50')

%plot(1:(N+1),fro_norm)
%title('Frobenius norm')


%(b) Nuclear norm trick
clear all
MovieRankings36 = readtable('MovieRankings36.csv');
MR36 = MovieRankings36{:,:};
[nrow , ncol] = size(MR36);
MRproj = MR36;
MRproj(isnan(MRproj)) = 0;
%MR36(isnan(MR36))= 0;
[U,S,V] = svd(MRproj,'econ');
%k = 7;
k = 3;
%Initialization step
lambda = 1;
Uk = U(:,1:k);
Vk = V(1:k,:);
Sk = S(1:k,1:k);
sl = functions_ex2(Sk,lambda,k);
X0 = Uk*sqrt(sl);
Y0 = Vk'*sqrt(sl);
%sl = functions_ex2(S,lambda,ncol);
%initialization step
M0 = X0*Y0';
M = M0;
R = MRproj-M;
%R(isnan(MR36)) = 0;
%objective function
N = 100;
obj_fnc2 = zeros(N+1,1);
obj_fnc2(1) = functions2_ex2(R,MRproj,M,lambda);

for iter = 1:N

R(isnan(MR36)) = 0;
%Update step
[Unew, Snew, Vnew] = svd(M + R,'econ');
dim = size(nonzeros(Snew));
Snew = Snew(1:dim,1:dim);
M = Unew*functions_ex2(Snew,lambda,dim)*Vnew';
R = MRproj-M;
%Objective function
obj_fnc2(iter+1) = functions2_ex2(R,MR36,M,lambda);
end

plot(1:(N+1),obj_fnc2)
%hold on




