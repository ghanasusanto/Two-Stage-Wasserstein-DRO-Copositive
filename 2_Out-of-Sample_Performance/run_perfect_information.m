function [ output_args ] = run_perfect_information(beta,cost,budget,holding_cost,shortage_cost)
clc
solver = 'mosek';
options = sdpsettings('verbose', 0, 'dualize', 0, 'solver', solver);
rng('default');
nIter = 100;

nProd = 3;
d = nProd;
N_1 = nProd;
N_2 = nProd+1;
h_ = holding_cost * ones(nProd,1);
s_ = shortage_cost * ones(nProd,1);
data.Q = zeros(N_2,d);
W = [eye(nProd); eye(nProd)];
e_ = [zeros(nProd,1); 1];
data.r = e_;
data.W = [e_'; (e_-[ones(nProd,1);0])'; W zeros(size(W,1),1)];
data.S = zeros(1,d);
data.t = 0;
data.N = 1;
c=cost*ones(N_1,1);
load('inputs');

obj_perfect_ = zeros(nIter,1);

for iter=1:nIter
    data_OOS = data;
    data_OOS.N = 20000;
    data_OOS.xi = OOS_xi{iter};

    yalmip('clear');
    alpha = sdpvar;
    x = sdpvar(nProd,1);
    F = [x>=0, sum(x)<=budget];
    obj1 = alpha;
    data_OOS.T = [zeros(1,d); zeros(1,d); diag(s_); -diag(h_)];
    data_OOS.h = [0; -alpha; -diag(s_)*x; diag(h_)*x];
    [G,obj2] = SAA_constraints_new(data_OOS);
    obj = obj1 + 1/beta*obj2;
    solvesdp(G, obj, options);
    obj_perfect = double(obj);
    fprintf ('iter %d perfect: %f \n', iter, obj_perfect);

    obj_perfect_(iter) = obj_perfect;
end
clear xi OOS_xi mu C N_
filename = sprintf('results/perfect_information');
save(filename);

end

