
%% set up optimization problem
[n,dim] = size(XX);
lam = 0.01;
Y = (label*ones(1,dim + 1)).*[XX,ones(n,1)];
w = [-1;-1;1;1];
fun = @(I,Y,w)fun0(I,Y,w,lam);
gfun = @(I,Y,w)gfun0(I,Y,w,lam);
Hvec = @(I,Y,w,v)Hvec0(I,Y,w,v,lam);

[w,f,gnorm] = SINewton(fun,gfun,Hvec,Y,w);
%[w,f,gnorm] = SINewton(fun,gfun,Hvec,M1,x);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n,dim] = size(XX);

%soft-margin matrices in objective function
M1 = [diag(ones(dim,1)),zeros(dim,n+1);zeros(n+1,dim),zeros(n+1,n+1)];
M2 = [zeros(1,dim+1),ones(1,n)];

fun = @(M1,x)fun_init(M1,x,M2);
gfun = @(M1,x)gfun_init(M1,x,M2);
Hvec = @(M1,p)Hvec_init(M1,p);
Hfun = @(M1)Hfun_init(M1);

%Constraints
Y0 = (label*ones(1,dim + 1)).*[XX,ones(n,1)];

A = [Y0, diag(ones(n,1));zeros(n,dim+1),  diag(ones(n,1))];

b = [ones(n,1);zeros(n,1)];

x = [w;zeros(n,1)];
%Active-constraints
W = find(A*x-b<=0);


[xiter,lm] = ASM_ex1(x,gfun,Hfun,A,b,W);

[niter, dimiter] = size(xiter);
x_idem = find(xiter(:,dimiter)<0);
x_idem(x_idem>58) = [];
x_igop = find(xiter(:,dimiter)>0);
x_igop(x_igop>58) = [];


%x_zero = find(xiter(:,59)==0);
%x_zero = find(x_zero<=58);

%SVN figure
close all
figure;
hold on; grid;
plot3(X(x_idem,i1),X(x_idem,i2),X(x_idem,i3),'.','color','b','Markersize',20);
plot3(X(x_igop,i1),X(x_igop,i2),X(x_igop,i3),'.','color','r','Markersize',20);
view(3)
fsz = 16;
set(gca,'Fontsize',fsz);
xlabel(str(i1),'Fontsize',fsz);
ylabel(str(i2),'Fontsize',fsz);
zlabel(str(i3),'Fontsize',fsz);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = fun_init(M1,x,M2)
C = 10^3;
f = 0.5*x'*M1*x + C*M2*x;
end

function g = gfun_init(M1,x,M2)
C = 10^3;
g = M1*x + C*M2';
end

function H = Hvec_init(M1,p)
H = M1*p;
end

function H = Hfun_init(M1)
H = M1;
end
