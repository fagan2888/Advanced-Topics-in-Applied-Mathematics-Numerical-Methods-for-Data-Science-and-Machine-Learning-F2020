function p = stochastic_finddirection(g,s,y,rho,x,xnew)
% input: g = gradient dim-by-1
% s = matrix dim-by-m, s(:,i) = x_{k-i+1}-x_{k-i}
% y = matrix dim-by-m, y(:,i) = g_{k-i+1}-g_{k-i}
% rho is 1-by-m, rho(i) = 1/(s(:,i)'*y(:,i))
m = size(s,2);
a = zeros(m,1);  
for i = 1 : m
    Ig = randperm(n,Ng);
    g = gfun(Ig,Y,x);
    a(i) = rho(i)*s(:,i)'*g;
    Ih = randperm(n,Nh);
    y(:,i) = gfun(Ih,Y,xnew) - gfun(Ih,Y,x);
    g = g - a(i)*y(:,i);
end
gam = s(:,1)'*y(:,1)/(y(:,1)'*y(:,1)); % H0 = gam*eye(dim)
g = g*gam;
for i = m :-1 : 1
    aux = rho(i)*y(:,i)'*g;
    g = g + (a(i) - aux)*s(:,i);
end
p = -g;
end