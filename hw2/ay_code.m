clear all
% noise term
mu_ek     = 0;
sigma_ek  = 0.1;
num_samps = 100;
ek        = normrnd(mu_ek,sigma_ek, num_samps,1);
% x

mu_x      = 0;
sigma_x   = 1;
xi        = normrnd(mu_x,sigma_x, num_samps,1);
% y
ti = 1 + 0.01 * xi -2 * xi.*xi + ek;
% by now, we created ti, xi

% now, I work on parameter estimation
lambda = linspace(0,10,100);
eps    = 0.0001;
phi_x  = [ones(100,1) xi xi.*xi];
alpha  = 1e-4;
WS     = zeros(100,3);
for l = 1:100
    t_lambda = lambda(l);
    w_r      = randn(3,1);
    for iter=1:10000
        err  = ti- phi_x * w_r;
        dw_0 = sum(err);
        d1   = 1/max(eps,abs(w_r(2)));
        dw_1 = sum(phi_x(:,2).*err) - t_lambda * w_r(2) * d1;
        d2   = 1/max(eps,abs(w_r(3)));
        dw_2 = sum(phi_x(:,3).*err) - t_lambda * w_r(3) * d2;
        
        w_r(1) = w_r(1) + alpha * dw_0;
        w_r(2) = w_r(2) + alpha * dw_1;
        w_r(3) = w_r(3) + alpha * dw_2;
    end
    WS(l,:)=w_r;
end

figure(1)
subplot(3,1,1)
plot(lambda,WS(:,1),'linewidth',2)
ylabel('W0');
xlabel('Lambda');
subplot(3,1,2)
plot(lambda,WS(:,2),'linewidth',2)
ylabel('W1');
xlabel('Lambda');
subplot(3,1,3)
plot(lambda,WS(:,3),'linewidth',2)
ylabel('W2');
xlabel('Lambda');
