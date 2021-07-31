function [ output_args ] = run_d_N_exact(d_iter,N_iter,i)
N_ = [5; 10; 20; 40; 80; 160; 320; 640; 1280];
d_ = [1; 2; 4; 8; 16; 32; 64; 128; 256; 512; 1024];
solver = 'mosek';
options = sdpsettings('verbose', 0, 'dualize', 1, 'solver', solver);
foldername = sprintf('results_new/%d_%d', d_iter, N_iter);
filename = sprintf('results_new/%d_%d/d_%d_N_%d_i_%d.mat', d_iter, N_iter, d_iter, N_iter, i);
if exist(filename, 'file') ~= 2
    return;
end
load(filename);
N = N_(N_iter);
d = d_(d_iter);
load('inputs');
epsilon = 1/sqrt(N);
tau = 0;
N_1 = 0;
N_2 = 1;
data.Q = zeros(N_2,d);
W = 1;
data.r = 1;
data.W = [1; 1];
data.S = eye(d);
data.t = ones(d,1);
data.N = N;
a = a{d_iter,i};
b = b{d_iter,i};
data.xi = xi;
data.T = [a'; zeros(1,d)];
data.h = [-b; 0];

yalmip('clear');
data_exact = data;
[F,obj] = exact_constraints_SOCP(epsilon,data_exact,0);
diagnostics_exact = solvesdp([F], obj, options);
obj_exact = double(obj);
obj_exact


filename = sprintf('results_new/%d_%d/d_%d_N_%d_i_%d_exact', d_iter, N_iter, d_iter, N_iter, i);
save(filename,'obj_exact','diagnostics_exact');

end

