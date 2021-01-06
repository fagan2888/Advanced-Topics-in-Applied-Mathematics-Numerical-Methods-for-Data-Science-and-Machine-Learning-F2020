%Exercise 4 text categorization
clear all
[M,y] = readdata();
[~ , ncol] = size(M);

k = 7;
alpha = 6;
%c = alpha*k;
%r = c;
Mmean = mean(M,2);        % find mean
M = M - Mmean*ones(1,ncol); 
 %Eckart-Young best approximation-Truncated SVD
    [U, S, V] = svd(M,'econ');
    Uk = U(:,1:k);
    Vtr = V';
    Vk = Vtr(1:k,:);
    Sk = S(1:k,1:k);
    %X0 = Uk*Sk;
    %Y0 = Vk;
    %Mk = X0*Y0;
    
%maximal leverage
pi = zeros(ncol,1);

for i = 1:ncol
    sum = sumfun(Vk,k,i);
    pi(i) = (1/k)*sum;
end
[pivals, idx] = sort(pi, 'descend');
% 8, 18, 46, 52, 69, 147, 219
maxidx = idx(1:k);
%features
% 8 -->  the
% 18 --> and
% 46 --> for
% 52 --> gif
% 36 --> this
% 147 --> contact
% 219 --> with

%No help 

%Criterion: Pick the columns with the highest leverage scores
idx10000 = idx(1:10000);
Mnew = M(:,idx10000);
[nrow , ncol] = size(Mnew);
M_mean = mean(Mnew,2);        % find mean
Mnew = Mnew - M_mean*ones(1,ncol); 
%(a)
[Ua, Sa, Va] = svd(Mnew,'econ');
%U2 = Ua(:,1:2);
%S2 = Sa(1:2,1:2);
%V2 = Vtr(1:2,:);
%V2 = V2';
C = Ua(:,1:2)*Sa(1:2,1:2); 
%label = zeros(nrow, 1);
Ind = find(y == -1);
Fl = find(y == 1);
figure(1);
hold on
% scatter(p1,p2,7,y,'filled')
plot(C(Ind,1),C(Ind,2),'.','color','b','Markersize',20);
plot(C(Fl,1),C(Fl,2),'.','color','r','Markersize',20);
xlabel('PC1'); ylabel('PC2')

%(b)
idx5 = idx(1:5);
M5 = M(:,idx5);
[nrow , ncol] = size(M5);
M5mean = mean(M5,2);        % find mean
M5 = M5 - M5mean*ones(1,ncol); 
[Ub, Sb, Vb] =svd(M5,'econ');%principal components
C = M5(:,1:2)*Vb(1:2,1:2);   % first 3 coefficients for each point, same as U(:,1:3)'*A;
label = zeros(nrow, 1);
Ind = find(y == -1);
Fl = find(y == 1);
figure(1);
hold on
% scatter(p1,p2,7,y,'filled')
plot(C(Fl,1),C(Fl,2),'.','color','r','Markersize',20);
plot(C(Ind,1),C(Ind,2),'.','color','b','Markersize',20);
xlabel('PC1'); ylabel('PC2')