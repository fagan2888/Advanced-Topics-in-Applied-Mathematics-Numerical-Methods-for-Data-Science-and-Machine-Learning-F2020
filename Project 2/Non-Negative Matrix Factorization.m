%Exercise 1 NMF
%a = fopen('MovieRankings36.csv');
%b = textscan(a, '%s %s %s %s', 'delimiter', ',', 'CollectOutput',true);
MovieRankings36 = MovieRankings36{:,:};
%[row, col] = find(isnan(MR));

%Complete submatrix (11x13)
MR=MovieRankings36([1,2,3,4,5,8,11,13,16,25,29],[1,2,3,5,6,9,10,11,13,16,18,19,20]);
[nrow , ncol] = size(MR);
%Computing optimal k
clust = zeros(nrow,nrow);
for i=1:nrow
 clust(:,i) = kmeans(MR,i);
end
evalclusters(MR,clust,'CalinskiHarabasz')
%optimal K is 10, second best is 9
%K=2 is good enough
indMR = kmeans(MR,2);
%ends up giving dramatically different cluster indices every time you run
%it and so it shows that the number of data points are not enough to give
%more solid results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute NMF (A~WH) with k=2
k = 2;
%(a)Projected GD
%initialization step
[~,H,~,W] = kmeans(MR,k);
%W = 5*rand(nrow,2);
%H = 5*rand(2,ncol);
R = MR-W*H;
%step size
alpha = .1;
%Number of steps
N = 100;
%tolerance
tolPGD = zeros(N+1,1);
tolPGD(1) = norm(MR-W*H, 'fro');
%number of steps
%iter = 0;
for i=1:N
    %PGD update step
    Wnew = max(W+alpha*R*H',0);
    Hnew = max(H+alpha*W'*R,0);
    tolPGD(i+1) = norm(MR-Wnew*Hnew, 'fro');
    
    %update
    W = Wnew;
    H = Hnew;
    R = MR-W*H;
end
fprintf("The best Frobenius norm difference between MR and WH is %f\n", tolPGD(N+1))
%plot(1:(N+1),tolPGD)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(b) Lee-Seung scheme

%initialization step
[~,H,~,W] = kmeans(MR,k);
%R = MR-W*H;
%tolerance
tolLS = zeros(N+1,1);
tolLS(1) = norm(MR-W*H, 'fro');
%number of steps
%iter = 0;
for i=1:N
    %LS update step
    Wnew = (W.*(MR*(H')))./(W*H*(H'));
    Hnew = (H.*((Wnew')*MR))./((Wnew')*Wnew*H);
    tolLS(i+1) = norm(MR-Wnew*Hnew, 'fro');
    
    %update
    W = Wnew;
    H = Hnew;
    %R = MR-W*H;
end
fprintf("The best Frobenius norm difference between MR and WH for Lee-Seung is %f\n", min(tolLS))
figure
plot(1:(N+1),tolPGD,1:(N+1),tolLS)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(c) PGD and then LS scheme
%PGD initialization
%W = 5*rand(nrow,2);
%H = 5*rand(2,ncol);
[~,H,~,W] = kmeans(MR,k);
R = MR-W*H;
%step size
alpha = .001;
%tolerance
tolPGD = zeros(N+1,1);
tolPGD(1) = norm(MR-W*H, 'fro');
%number of steps
%iter = 0;
for i=1:N
    %PGD update step
    Wnew = max(W+alpha*R*H',0);
    Hnew = max(H+alpha*W'*R,0);
    tolPGD(i+1) = norm(MR-Wnew*Hnew, 'fro');
    
    %update
    W = Wnew;
    H = Hnew;
    R = MR-W*H;
end


Wls = W;
Hls = H;
R = MR-Wls*Hls;
%tolerance
tolPGD_LS = zeros(N+1,1);
tolPGD_LS(1) = norm(MR-Wls*Hls, 'fro');
%number of steps
%iter = 0;
for i=1:N
    %LS update step
    Wnew = (Wls.*(MR*(Hls')))./(Wls*Hls*(Hls'));
    Hnew = (Hls.*((Wnew')*MR))./((Wnew')*Wnew*Hls);
    tolPGD_LS(i+1) = norm(MR-Wnew*Hnew, 'fro');
    
    %update
    Wls = Wnew;
    Hls = Hnew;
    R = MR-Wls*Hls;
end

figure
plot(10:(N+1),tolPGD(10:(N+1)),10:(N+1),tolLS(10:(N+1)),10:(N+1),tolPGD_LS(10:(N+1)))
legend('PGD','LS','PGD+LS')
%set(gca,'Yscale','log')