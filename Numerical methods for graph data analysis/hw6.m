%Exercise 1
%Find roots of f(S)=1-S-exp(-zS) and <s> numerically (via Newton Method)
f11 = @(S,z) 1-S-exp(-z.*S);
df11 = @(S,z) -1+exp(-z.*S);
f12 = @(S,z) 1/(1-z+z*S);
%df12 = @(S,z) (-z)/((1-z+z*S)^2);
gss = 2/3;
ite = 100;
tol = 1e-4;
%i=1;
approx11 = zeros(10,1);
approx12 = zeros(10,1);
%z=0;
%z=1;
x = linspace(.01,4,10);
    %[apprx,iter] = newt(f,df,gss,ite,tol,z);
[~,length] = size(x);
      %f = f(gss,z(i));
      for i = 1:length
      z = x(i);
      apprx11 = newt0(f11,df11,gss,ite,tol,z);
      approx11(i) = apprx11;
      if apprx11 <= 0
          approx11(i) = 0;
      end
      end
      %gss = 100;
      %tol = 1;
      for j = 1:length
      z = x(j);
      %apprx12 = newt0(f12,df12,gss,ite,tol,z);
      %approx12(j) = apprx12;
      approx12(j) = f12(approx11(j),z);
      if approx12 <= 0
          approx12(j) = 0;
      end
      end

%legend('Numerically','Experimentally')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Find S(z) and <s> experimentally
%Create ER random graphs
N = 1000;
x = linspace(.01,4,10);
[~,length] = size(x);
nongsmean = zeros(length,1);
meansize = zeros(length,1);
      for i = 1:length
      z = x(i);
      p = z/(N-1);
      nongiant2 = zeros(100,1);
      nongiant3 = zeros(100,1);
      giantsize = zeros(100,1);
      nongs = zeros(100,1);
      countCOMP = zeros(100,1);
        for r = 1:100
            G = rand(N,N) < p;
            G = triu(G,1);
            G = G + G';
            G = digraph(G);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            c = conncomp(G,'Type','weak');
                for l=1:N
                    if c(l) == 1
                        giantsize(r) = giantsize(r) + 1;
                    end
                    if c(l) == 2
                        nongiant2(r) = nongiant2(r) + 1;
                    end
                    if c(l) == 3
                        nongiant3(r) = nongiant3(r) + 1;
                    end
                end
            %nongs(r) = (N-giantsize(r))/max(c);
            nongs(r) = (nongiant2(r)+nongiant3(r))/2;
            v = dfsearch(G,1,'allevents');
            v = table2array(v(:,2:4));
            [rowdim, ~] = size(v);
            countCOMP(r) = rowdim;
            for j = 1:rowdim
                if num2str(v(j,1)) == "NaN"
                    countCOMP(r) = countCOMP(r) - 1;
                end
            end            
        end
        meansize(i) = sum(countCOMP)/100;
        nongsmean(i) = sum(nongs)/100;
      end
meansize = meansize./(2*N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot the fraction S of the vertices of a giant component
figure
plot(x,approx11,'-o','LineWidth', 2)
hold on
plot(x,meansize,'g','LineWidth', 2)
title('fraction S of vertices')
xlabel('z')
ylabel('f(S)')
legend('Numerically','Experimentally')

%title('numerical vs experimental')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot the average size <s> of a non-giant component
figure
plot(x,approx12,'-o','LineWidth', 2)
hold on
plot(x,nongsmean/10,'g','LineWidth', 2)
title('Mean size of a non-giant component')
xlabel('z')
ylabel('<s>')
legend('Numerically','Experimentally')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Exercise 2
z = 4;
%average shortest path
ll = zeros(4,1);
tl = zeros(4,1);
NN = zeros(4,1);
for i = 1:4
    N = 2^(i+9);
    p = z/(N-1);
    NN(i) = N;
    %Theoritical average shortest path
    tl(i) = log(N)/log(z);
    %Generate ER(N,p)
    G = rand(N,N) < p;
    G = triu(G,1);
    G = G + G';
    G = digraph(G);
    %Choose randomly 100 vertices
    rows = randperm(N, 100);
    %shortest path for each randomly selected vertex
    ddd = zeros(100,1);
    for r = 1:100
        dd = 0;
        %shortest path for each vertex
        for l = 1:N
        [P,d] = shortestpath(G,rows(r),l);
        if d == Inf
            d = 0;
        end
        dd = dd + d;
        end
        ddd(r) = dd/N;
    %v = bfsearch(G,rows(r),'allevents');
    %c = conncomp(G,'Type','weak');
    end
    ll(i) = sum(ddd)/100;
end

%Plot the average shortest path l(n)
figure
plot(NN,tl,'-o','LineWidth', 2)
hold on
plot(NN,ll,'g','LineWidth', 2)
title('Average shortest path')
xlabel('N')
ylabel('l(N)')
legend('Theoritically','Experimentally')
