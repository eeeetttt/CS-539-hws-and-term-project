%% Rejection Sampling
%Sample from a step function over [0,1]:
X = sampleDist(@(x)1.3*(x>=0&x<0.7)+0.3*(x>=0.7&x<=1),...
1.3,1e6,[0,1],true);
%Sample from a normal distribution over [-5,5]:
X = sampleDist(@(x) 1/sqrt(2*pi) *exp(-x.^2/2),...
1/sqrt(2*pi),1e6,[-5,5],true);


%% MCMC
%https://www.mathworks.com/matlabcentral/fileexchange/47912-markov-chain-monte-carlo-sampling-of-posterior-distribution
data=randn(100,1)*2+3;
logmodelprior=@(m)0; %log(normpdf(m(2),0,2)); %use a flat prior.
loglike=@(m)sum(log(normpdf(data,m(1),m(2))));
minit=[0 1];
m=mcmc(minit,loglike,logmodelprior,[.2 .5],10000);
m(1:100,:)=[]; %crop drift
plotmatrix(m);

close all
%% MCMC
% https://www.mathworks.com/matlabcentral/fileexchange/49820-ensemble-mcmc-sampler
%define problem:
mu = [5;-3;6];
C = [.5 -.4 0;-.4 .5 0; 0 0 1];
iC=pinv(C);
logPfuns={@(m)-0.5*sum((m-mu)'*iC*(m-mu))}

%make a set of starting points for the entire ensemble of walkers
minit=randn(length(mu),length(mu)*2);

%Apply the MCMC hammer
[models,logP]=gwmcmc(minit,logPfuns,100000);
models(:,:,1:floor(size(models,3)*.2))=[]; %remove 20% as burn-in
models=models(:,:)'; %reshape matrix to collapse the ensemble member dimension
scatter(models(:,1),models(:,2))
prctile(models,[5 50 95])

scatter(models(:,1),models(:,3))
prctile(models,[5 50 95])
