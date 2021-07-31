function [F,obj] = dec_rules_constraints(epsilon,data,obj)
%PRIMAL_CONSTRAINTS Generate primal constraints
N = data.N;
xi = data.xi;
Q = data.Q;
r = data.r;
L = length(data.r);
h = data.h;
T = data.T;
W = data.W;
S = data.S;
t = data.t;
M = size(T,1);
J = size(S,1);
d = size(S,2);
N_2 = size(Q,1);

lambda = sdpvar;
constraints = {};
constraints{end+1} = [lambda >= 0];
obj = obj + epsilon^2*lambda;
s = sdpvar(N,1);
sum_Y = 0;
sum_y = 0;
sum_y_0 = 0;
for k=1:L
    Y{k} = sdpvar(d);
    y{k} = sdpvar(d,1);
    y_0{k} = sdpvar;
    sum_Y = sum_Y + Y{k};
    sum_y = sum_y + y{k};
    sum_y_0 = sum_y_0 + y_0{k};
end
for k=1:L
    theta = sdpvar(d,1); eta = sdpvar(d,1);
    constraints{end+1} = [theta >= 0, eta >= 0];
    constraints{end+1} = [[Y{k} 0.5*(y{k}-T(k,:)'-theta+eta); 0.5*(y{k}-T(k,:)'-theta+eta)' y_0{k}-h(k)-eta'*ones(d,1)] >= 0];
    theta = sdpvar(d,1); eta = sdpvar(d,1);
    constraints{end+1} = [theta >= 0, eta >= 0];
    constraints{end+1} = [[Y{k} 0.5*(y{k}-theta+eta); 0.5*(y{k}-theta+eta)' y_0{k}-eta'*ones(d,1)] >= 0];
end

for i=1:N
    curr_xi = xi(:,i);
    psi = sdpvar(d,1);
    phi = sdpvar(d,1);
    constraints{end+1} = [psi >= 0, phi >= 0];
    obj = obj + 1/N * s(i);
    constraints{end+1} = [[lambda*eye(d)-sum_Y -0.5*(sum_y+2*lambda*curr_xi+phi-psi); -0.5*(sum_y+2*lambda*curr_xi+phi-psi)' s(i)-sum_y_0+lambda*norm(curr_xi)^2-psi'*ones(d,1)] >= 0];  
end
F = [constraints{:}];
end