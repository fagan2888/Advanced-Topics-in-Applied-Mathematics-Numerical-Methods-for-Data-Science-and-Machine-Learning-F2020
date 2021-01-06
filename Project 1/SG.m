%SG
function [a, normgrad, x, f] = SG(batch,Y, w, fun, gfun)%rho)
k = 0;
kmax =1000;
grad = zeros(4,1);
n = size(Y,1);
B = ceil(n/batch);
a = 1; %decreasing step size
%a = .1; %initial stepsize;
x = w;
normgrad = zeros(kmax,1);
f = zeros(kmax + 1,1);
I = randperm(n,batch);
f(1) = fun(I,Y,x);
X = zeros( kmax, 4);
%rerr = 1;
%normb = norm(b);
   for l = 1:kmax
            I = randperm(n,batch);
            grad = gfun(I,Y,x);
            normgrad(l) = norm(grad);
            x = x - a*grad;
            a = max(0, a - 1/(kmax+1));
            l = l + 1;
            f(l) = fun(I,Y,x);
    
    %rerr = norm(r)/normb;
   end
end

function f = fun0(I,Y,w,lam)
f = sum(log(1 + exp(-Y(I,:)*w)))/length(I) + 0.5*lam*w'*w;
end
%%