function [ output_args ] = run_N_i(N,i)
clc
addpath /homes/gah11:/homes/gah11/Projects/SeDuMi_1_3:/homes/gah11/Projects/SDPT3-4.0/Solver:/homes/gah11/Projects/mosek/7/toolbox/r2013a:/vol/bitbucket/gah11/ibm/ILOG/CPLEX_Studio126/cplex/matlab:/homes/gah11/Projects/tools/gurobi562/linux64/matlab:/homes/gah11/Projects/tools/yalmip/operators:/homes/gah11/Projects/tools/yalmip/modules/sos:/homes/gah11/Projects/tools/yalmip/modules/global:/homes/gah11/Projects/tools/yalmip/modules/moment:/homes/gah11/Projects/tools/yalmip/modules/parametric:/homes/gah11/Projects/tools/yalmip/solvers:/homes/gah11/Projects/tools/yalmip:/homes/gah11/Projects/tools/yalmip/modules:/homes/gah11/Projects/tools/yalmip/demos:/homes/gah11/Projects/tools/yalmip/extras
N_ = [5; 10; 20; 40; 80; 160; 320; 640; 1280];
[beta,cost,budget,holding_cost,shortage_cost] = generate_data();
solver = 'mosek';
options = sdpsettings('verbose', 0, 'dualize', 0, 'solver', solver);
rng('default');
Eps = [0; 2^(-8); 2^(-7); 2^(-6); 2^(-5); 2^(-4); 2^(-3); 2^(-2); 2^(-1); 0.75];

nProd = 3;
tau = 0;
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
data.N = N;
c=cost*ones(N_1,1);
load('inputs');

data.xi = xi{find(N_==N),i};
data_OOS = data;
data_OOS.N = 20000;
data_OOS.xi = OOS_xi{i};
data_cv = data;
data_cv_OOS = data;

K = 3;
index = crossvalind('Kfold',N,K);
ave_obj = zeros(length(Eps),1);
for epsIter = 1:length(Eps)
    eps = Eps(epsIter);
    objs_ = zeros(K,1);
    for k=1:K
        %% bootstrap
%         while(true)
%             idx = unidrnd(N,N,1);
%             OOS_idx = setdiff([1:N],idx);
%             if (OOS_idx)
%                 break
%             end
%         end
%         cv_xi = data.xi(:,idx);
%         cv_xi_OOS = data.xi(:,OOS_idx);
%         data_cv.N = size(cv_xi,2);
%         data_cv.xi = cv_xi;
%         data_cv_OOS.N = size(cv_xi_OOS,2);
%         data_cv_OOS.xi = cv_xi_OOS;
        
        
        %% crossvalidation
        cv_xi = data.xi(:,index~=k);
        cv_xi_OOS = data.xi(:,index==k);
        data_cv.N = size(cv_xi,2);
        data_cv.xi = cv_xi;
        data_cv_OOS.N = size(cv_xi_OOS,2);
        data_cv_OOS.xi = cv_xi_OOS;
% 
        %% in sample CV
        yalmip('clear');
        alpha = sdpvar;
        x = sdpvar(nProd,1);
        F = [x>=0, sum(x)<=budget];
        obj1 = alpha;
        data_robust = data_cv;
        data_robust.Q = zeros(nProd,d);
        data_robust.W = W;
        data_robust.r = ones(nProd,1);
        data_robust.T = [diag(s_); -diag(h_)];
        data_robust.h = [-diag(s_)*x; diag(h_)*x];
        [G,obj2] = primal_constraints_new(eps,tau,data_robust,alpha);
        obj = obj1 + 1/beta*obj2;
        solvesdp([F, G], obj, options);
        x_cv = double(x);
        
        %% out sample CV
        yalmip('clear');
        alpha = sdpvar;
        x = x_cv;
        obj1 = alpha;
        data_cv_OOS.T = [zeros(1,d); zeros(1,d); diag(s_); -diag(h_)];
        data_cv_OOS.h = [0; -alpha; -diag(s_)*x; diag(h_)*x];
        [G,obj2] = SAA_constraints_new(data_cv_OOS);
        obj = obj1 + 1/beta*obj2;
        solvesdp(G, obj, options);
        objs_(k) = double(obj);
        %fprintf ('crossvalidation-> N %f N_OOS %f, k %d epsilon %f: obj %f \n', data_cv.N, data_cv_OOS.N, k, eps, objs_(k));
        
    end
    ave_obj(epsIter) = mean(objs_);
end
[tmp,best_id] = min(ave_obj);
epsilon = Eps(best_id);
fprintf ('Wassertein-> N %d iter %d epsilon %f: ', N, i, epsilon);

data.xi = xi{find(N_==N),i};
data_OOS = data;
data_OOS.N = 20000;
data_OOS.xi = OOS_xi{i};

%% SAA problem
yalmip('clear');
alpha = sdpvar;
x = sdpvar(nProd,1);
F = [x>=0, sum(x)<=budget];
obj1 = alpha;
data.T = [zeros(1,d); zeros(1,d); diag(s_); -diag(h_)];
data.h = [0; -alpha; -diag(s_)*x; diag(h_)*x];
[G,obj2] = SAA_constraints_new(data);
obj = obj1 + 1/beta*obj2;
diagnostics_SAA = solvesdp([F, G], obj, options);
obj_SAA = double(obj);
x_SAA = double(x);

%% OOS SAA
yalmip('clear');
alpha = sdpvar;
x = x_SAA;
obj1 = alpha;
data_OOS.T = [zeros(1,d); zeros(1,d); diag(s_); -diag(h_)];
data_OOS.h = [0; -alpha; -diag(s_)*x; diag(h_)*x];
[G,obj2] = SAA_constraints_new(data_OOS);
obj = obj1 + 1/beta*obj2;
diagnostics_SAA_OOS = solvesdp(G, obj, options);
obj_SAA_OOS = double(obj);

%% robust problem
yalmip('clear');
alpha = sdpvar;
x = sdpvar(nProd,1);
F = [x>=0, sum(x)<=budget];
obj1 = alpha;
data_robust = data;
data_robust.Q = zeros(nProd,d);
data_robust.W = W;
data_robust.r = ones(nProd,1);
data_robust.T = [diag(s_); -diag(h_)];
data_robust.h = [-diag(s_)*x; diag(h_)*x];
[G,obj2] = primal_constraints_new(epsilon,tau,data_robust,alpha);
obj = obj1 + 1/beta*obj2;
diagnostics_robust = solvesdp([F, G], obj, options);
obj_robust = double(obj);
x_robust = double(x);
fprintf ('robust: %f SAA: %f ', obj_robust, obj_SAA);

%% OOS robust
yalmip('clear');
alpha = sdpvar;
x = x_robust;
obj1 = alpha;
data_OOS.T = [zeros(1,d); zeros(1,d); diag(s_); -diag(h_)];
data_OOS.h = [0; -alpha; -diag(s_)*x; diag(h_)*x];
[G,obj2] = SAA_constraints_new(data_OOS);
obj = obj1 + 1/beta*obj2;
diagnostics_robust_OOS = solvesdp(G, obj, options);
obj_robust_OOS = double(obj);

fprintf ('OOS robust: %f SAA: %f improvement: %f%% \n', obj_robust_OOS, obj_SAA_OOS, (obj_SAA_OOS/obj_robust_OOS-1)*100);
filename = sprintf('results/%d/N_%d_i_%d_Wasserstein', N, N, i);
clear xi OOS_xi mu C N_
save(filename,'obj_robust','obj_SAA','obj_robust_OOS','obj_SAA_OOS','x_robust','x_SAA','diagnostics_SAA','diagnostics_robust','epsilon');
end

