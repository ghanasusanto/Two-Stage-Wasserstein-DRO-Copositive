function [ output_args ] = run_d_N(d_iter,N_iter,i)
rng('shuffle');
solver = 'mosek';
options = sdpsettings('verbose', 0, 'dualize', 0, 'solver', solver);
N_ = [5; 10; 20; 40; 80; 160; 320; 640; 1280];
d_ = [1; 2; 4; 8; 16; 32; 64; 128; 256; 512; 1024];

N = N_(N_iter);
d = d_(d_iter);
load('inputs');
epsilon = 1/sqrt(N);
tau = 0;

A = A{d_iter,i};
b = b{d_iter,i};
L =length(b);
N_1 = 0;
N_2 = L;
data.Q = zeros(N_2,d);
data.r = ones(L,1);
data.W = [eye(L); eye(L)];
data.S = eye(d);
data.t = ones(d,1);
data.N = N;
data.xi = rand(d,N);
xi = data.xi;

data.T = [A; zeros(L,d)];
data.h = [-b; zeros(L,1)];
foldername = sprintf('results_new/%d_%d', d_iter, N_iter);
mkdir(foldername);
filename = sprintf('results_new/%d_%d/d_%d_N_%d_i_%d', d_iter, N_iter, d_iter, N_iter, i);

yalmip('clear');
data_robust = data;
[F,obj] = primal_constraints(epsilon,tau,data_robust,0);
diagnostics_robust = solvesdp([F], obj, options);
obj_robust = double(obj);

yalmip('clear');
data_dec_rules = data;
[F,obj] = dec_rules_constraints(epsilon,data_dec_rules,0);
diagnostics_dec_rules = solvesdp([F], obj, options);
obj_dec_rules = double(obj);

yalmip('clear');
data_exact = data;
[F,obj] = exact_constraints_SOCP(epsilon,data_exact,0);
diagnostics_exact = solvesdp([F], obj, options);
obj_exact = double(obj);

%fprintf('obj_robust: %f obj_dec_rules: %f analytical: %f improvement: %f \n', obj_robust, obj_dec_rules, analytical_value, (obj_dec_rules/obj_robust-1)*100);
save(filename,'obj_robust','obj_dec_rules','obj_exact','diagnostics_robust','diagnostics_dec_rules','diagnostics_exact','xi');
quit;
end

