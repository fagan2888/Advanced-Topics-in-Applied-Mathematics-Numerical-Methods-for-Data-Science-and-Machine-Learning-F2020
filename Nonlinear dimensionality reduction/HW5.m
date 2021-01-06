%HW5

%Scurve Data
 data3 = MakeScurveData();

%Perturbed Scurve Data
[n , dim] = size(data3);
perturbed_data3 = data3 + sqrt(.02)*randn(n,dim);
figure(2);
plot3(perturbed_data3(:,1),perturbed_data3(:,2),perturbed_data3(:,3),'.','Markersize',20);
daspect([1,1,1]);
set(gca,'fontsize',16);
view(3);
grid
save('PerturbedScurveData.mat','perturbed_data3');

%Make Emoji Data
dataemoji = MakeEmojiData();
%PCA
[coeff1,score1,latent1] = pca(data3);
[coeff2,score2,latent2] = pca(perturbed_data3);
[coeff3,score3,latent3] = pca(dataemoji);
figure(3);
%Scurvy data visualization
plot3(score1(:,1),score1(:,2),score1(:,3),'.','Markersize',20);
title("PCA on Scurvy Data")
%Perturbed Scurvy data visualization
figure(4);
plot3(score2(:,1),score2(:,2),score2(:,3),'.','Markersize',20);
title("PCA on Perturbed Scurvy Data (variance = 0.02")
%Emoji data visualization
figure(5);
plot3(score3(:,1),score3(:,2),score3(:,3),'.','Markersize',20);
title("PCA on emoji data")

%ISOMAP
dat = load('PerturbedScurveData.mat');
X = dat.perturbed_data3;
[n,~] = size(X);
%% compute pairwise distances
d = zeros(n);
e = ones(n,1);
for i = 1 : n
d(i,:) = sqrt(sum((X - e*X(i,:)).^2,2));
end
iso1 = IsoMap(d,'k',20);

%LLE
lle1 = lle(data3',10,3);
lle1 = lle1';
figure(6);
plot3(lle1(:,1),lle1(:,2),lle1(:,3),'.','Markersize',20);
daspect([1,1,1]);
set(gca,'fontsize',16);
view(3);
grid
title("LLE on Scurvy data")

lle2 = lle(perturbed_data3',10,3);
lle2 = lle2';
figure(6);
plot3(lle2(:,1),lle2(:,2),lle2(:,3),'.','Markersize',20);
daspect([1,1,1]);
set(gca,'fontsize',16);
view(3);
grid
title("LLE on Perturbed Scurvy data")

dataemoji((4*32+1):(5*32),:) = [];
dataemoji((12*32+1):(13*32),:) = [];
dataemoji((20*32+1):(21*32),:) = [];
dataemoji((28*32+1):(29*32),:) = [];
lle3 = lle(dataemoji',10,1600);
lle3 = lle3';
figure(6);
plot3(lle3(:,1),lle3(:,2),lle3(:,3),'.','Markersize',20);
daspect([1,1,1]);
set(gca,'fontsize',16);
view(3);
grid
title("LLE on Emoji data")

%t-SNE
%[Y,loss] = tsne(X,'Algorithm','exact','Perplexity',30,'Exaggeration',4);
[Y,loss] = tsne(dataemoji,'Algorithm','exact','NumDimensions',3);
gscatter(Y(:,1),Y(:,2),Y(:,3))
title("t-SNE on Emoji Data")

%Diffusion Maps

%% visualize the data
[N ,~] = size(dataemoji);
cmap = jet(N);
scatter3(perturbed_data3(:,1),perturbed_data3(:,2),perturbed_data3(:,3),20,cmap);
title('Scurve data');

%% Changing these values will lead to different nonlinear embeddings
knn    = ceil(0.03*N); % each patch will only look at its knn nearest neighbors in R^d
sigma2 = 100; % determines strength of connection in graph... see below

%% now let's get pairwise distance info and create graph 
m                = size(dataemoji,1);
dt               = squareform(pdist(dataemoji));
[srtdDt,srtdIdx] = sort(dt,'ascend');
dt               = srtdDt(1:knn+1,:);
nidx             = srtdIdx(1:knn+1,:);

% nz   = dt(:) > 0;
% mind = min(dt(nz));
% maxd = max(dt(nz));

% compute weights
tempW  = exp(-dt.^2/sigma2); 

% build weight matrix
i = repmat(1:m,knn+1,1);
W = sparse(i(:),double(nidx(:)),tempW(:),m,m); 
W = max(W,W'); % for undirected graph.

% The original normalized graph Laplacian, non-corrected for density
ld = diag(sum(W,2).^(-1/2));
DO = ld*W*ld;
DO = max(DO,DO');%(DO + DO')/2;

% get eigenvectors
[v,d] = eigs(DO,10,'la');

eigVecIdx = nchoosek(2:4,2);
for i = 1:size(eigVecIdx,1)
    figure,scatter(v(:,eigVecIdx(i,1)),v(:,eigVecIdx(i,2)),20,cmap)
    title('Nonlinear embedding');
    xlabel(['\phi_',num2str(eigVecIdx(i,1))]);
    ylabel(['\phi_',num2str(eigVecIdx(i,2))]);
end

figure,subplot(1,2,1)
scatter3(dataemoji(:,1),dataemoji(:,2),dataemoji(:,3),20,cmap);
title('Emoji data');
subplot(1,2,2)
scatter3(v(:,2),v(:,3),v(:,4),20,cmap)
title('Nonlinear embedding')

