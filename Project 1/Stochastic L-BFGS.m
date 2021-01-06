function Ex4
%% read data
A2012 = readmatrix('A2012.csv');
A2016 = readmatrix('A2016.csv');
% Format for A2012 and A2016:
% FIPS, County, #DEM, #GOP, then <str> up to Unemployment Rate
str = ["Median Income", "Migration Rate", "Birth Rate",...
"Death Rate", "Bachelor Rate", "Unemployment Rate","log(#Votes)"];
%
% remove column county that is read by matlab as NaN
A2012(:,2) = [];
A2016(:,2) = [];
%% Remove rows with missing data
A = A2016;
% remove all rows with missing data
ind = find(~isfinite(A(:,2)) |  ~isfinite(A(:,3)) | ~isfinite(A(:,4)) ...
    | ~isfinite(A(:,5)) | ~isfinite(A(:,6)) | ~isfinite(A(:,7)) ...
    | ~isfinite(A(:,8)) | ~isfinite(A(:,9)));
A(ind,:) = [];
%% select CA, OR, WA, NJ, NY counties
% ind = find((A(:,1)>=6000 & A(:,1)<=6999))%  ...  %CA
 % | (A(:,1)>=53000 & A(:,1)<=53999) ...        %WA
 % | (A(:,1)>=34000 & A(:,1)<=34999)); %...        %NJ  
  %| (A(:,1)>=36000 & A(:,1)<=36999) ...        %NY
  %| (A(:,1)>=41000 & A(:,1)<=41999));          %OR
% A = A(ind,:);

[n,dim] = size(A);

%% assign labels: -1 = dem, 1 = GOP
idem = find(A(:,2) >= A(:,3));
igop = find(A(:,2) < A(:,3));
num = A(:,2)+A(:,3);
label = zeros(n,1);
label(idem) = -1;
label(igop) = 1;

% %% select max subset of data with equal numbers of dem and gop counties
ngop = length(igop);
ndem = length(idem);
if ngop > ndem
    rgop = randperm(ngop,ndem);
    Adem = A(idem,:);
    Agop = A(igop(rgop),:);
    A = [Adem;Agop];
else
    rdem = randperm(ndem,ngop);
    Agop = A(igop,:);
    Adem = A(idem(rdem),:);
    A = [Adem;Agop];
end  
[n,dim] = size(A);
idem = find(A(:,2) >= A(:,3));
igop = find(A(:,2) < A(:,3));
num = A(:,2)+A(:,3);
label = zeros(n,1);
label(idem) = -1;
label(igop) = 1;

%% set up data matrix and visualize

X = [A(:,4:9),log(num)];
X(:,1) = X(:,1)/1e4;
% select three data types that distinguish dem and gop counties the most
i1 = 1; % Median Income
i2 = 7; % log(# votes)
i3 = 5; % Bachelor Rate

%% rescale data to [0,1] and visualize
XX = X(:,[i1,i2,i3]); % data matrix
% rescale all data to [0,1]
xmin = min(XX(:,1)); xmax = max(XX(:,1));
ymin = min(XX(:,2)); ymax = max(XX(:,2));
zmin = min(XX(:,3)); zmax = max(XX(:,3));
X1 = (XX(:,1)-xmin)/(xmax-xmin);
X2 = (XX(:,2)-ymin)/(ymax-ymin);
X3 = (XX(:,3)-zmin)/(zmax-zmin);
XX = [X1,X2,X3];

fsz = 16;
[n, dim] =size(XX);
lam = 0.01;
Y = (label*ones(1,dim + 1)).*[XX,ones(n,1)];
w = [-1;-1;1;1];
func = @(I,Y,w)fun0(I,Y,w,lam);
gfun = @(I,Y,w)gfun0(I,Y,w,lam);
Hvec = @(I,Y,w,v)Hvec0(I,Y,w,v,lam);

Ng = 10;
Nh = 20;
a = 5;
gam = 0.9; % line search step factor
jmax = ceil(log(1e-14)/log(gam)); % max # of iterations in line search
eta = 0.5; % backtracking stopping criterion factor
tol = 1e-10;
m = 5; % the number of steps to keep in memory
%% 
x0 = [1,1,-1,-1];  %initial guess
%
s = zeros(4,m);
y = zeros(4,m);
rho = zeros(1,m);
gnorm = zeros(1,1000);
%
x = x0';
f = zeros(1000,1);
%Stochastic
Ig = randperm(n,Ng);
g = gfun(Ig,Y,x);
f(1) = func(Ig,Y,x);
gnorm(1) = norm(g);
a = 1;
xnew = x - a*g;

gnew = gfun(Ig,Y,xnew);
s(:,1) = xnew - x;
%Stochastic 
Ih = randperm(n,Nh);
y(:,1) = gfun(Ih,Y,xnew) - gfun(Ih,Y,x);
rho(1) = 1/(s(:,1)'*y(:,1));
x = xnew;
g = gnew;
nor = norm(g);
gnorm(2) = nor;
iter = 1;
%Frequency of the pair updates
M = 10;
for iter = 1:1000
    if iter < m
        I = 1 : iter;
        p = finddirection(g,s(:,I),y(:,I),rho(I));
    else
        p = finddirection(g,s,y,rho);
    end
%     [a,j] = linesearch(x,p,g,func,eta,gam,jmax);
%     if j == jmax
%         p = -g;
%         [a,j] = linesearch(x,p,g,func,eta,gam,jmax);
%     end

    a = 1/iter;
    step = a*p;
    xnew = x + step; 
        if mod(iter,M) == 0
    Ig = randperm(n,Ng);
    gnew = gfun(Ig,Y,xnew);

    s = circshift(s,[0,1]); 
    y = circshift(y,[0,1]);
    rho = circshift(rho,[0,1]);
    s(:,1) = step;
    Ih = randperm(n,Nh);
    y(:,1) = gfun(Ih,Y,xnew) - gfun(Ih,Y,x);
    disp(iter)

    if(iter>m)
        l = double(mod(iter,m-1));
        rho(l+1) = 1/(step'*y(:,1));
        
    end
    end
    x = xnew;
    g = gnew;



    nor = norm(g);

    iter = iter + 1;
    gnorm(iter) = nor;
    f(iter) = func(I,Y,x);
end
f = f(f>0);
fprintf('L-BFGS: %d iterations,f = %d ,norm(g) = %d\n',iter,f(end),nor);
fsz = 16;
%%
figure;
hold on;
grid;
niter = length(f);
plot((0:niter-1)',f,'Linewidth',2);
set(gca,'Fontsize',fsz);
xlabel('k','Fontsize',fsz);
ylabel('f','Fontsize',fsz);
%%
figure;
hold on;
grid;
niter = length(gnorm);
plot((0:niter-1)',gnorm,'Linewidth',2);
set(gca,'Fontsize',fsz);
set(gca,'YScale','log');
xlabel('k','Fontsize',fsz);
ylabel('|| stoch grad f||','Fontsize',fsz);

end

