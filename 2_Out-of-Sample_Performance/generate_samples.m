function [ xi ] = generate_samples(N,d,C,mu)
sigma = 0.25*ones(d,1);
xi = exp(mvnrnd(mu, diag(sigma)*C*diag(sigma), N)');
end

