function [F,obj] = exact_constraints_SOCP(epsilon,data,obj)
%PRIMAL_CONSTRAINTS Generate primal constraints
N = data.N;
xi = data.xi;
Q = data.Q;
r = data.r;
h = data.h;
L = length(data.r);
T = data.T;
W = data.W;
S = data.S;
t = data.t;
M = size(T,1);
J = size(S,1);
d = size(S,2);
N_2 = size(Q,1);

T_ = data.T(1:L,:);
h_ = data.h(1:L,:);
lambda = sdpvar;
constraints = {};
constraints{end+1} = [lambda >= 0];
obj = obj + epsilon^2*lambda;
s = sdpvar(N,1);
constraints{end+1} = [s >= 0];

pad_template = strcat('%0',num2str(L),'s');
nPiece = 2^L;
T = cell(2^L,1);
h = cell(2^L,1);
for k=1:nPiece
    curr = char(dec2base(k-1,2));
    curr = sprintf(pad_template, curr);
    curr_index = (double(curr)-47)';
    T{k} = sum(T_((curr_index==1),:),1)'+zeros(d,1);
    h{k} = sum(h_(curr_index==1))'+0;
end
constraints_ = cell(N,1);
for i=1:N
    curr_xi = xi(:,i);
    obj = obj + 1/N * s(i);
    constraints_1 = cell(nPiece,1);
    for k=1:nPiece
        theta = sdpvar(d,1);eta = sdpvar(d,1);
        constraints_1{k} = [theta>=0, eta>=0, ...
            s(i)-h{k}+lambda*norm(curr_xi)^2-eta'*ones(d,1)>=0, ...
            norm([(-T{k}-2*lambda*curr_xi-theta+eta); s(i)-h{k}+lambda*norm(curr_xi)^2-eta'*ones(d,1)-lambda])<=lambda+s(i)-h{k}+lambda*norm(curr_xi)^2-eta'*ones(d,1)];
    end
    constraints_{i} = [constraints_1{:}];
end
F = [constraints{:}, constraints_{:}];
end

