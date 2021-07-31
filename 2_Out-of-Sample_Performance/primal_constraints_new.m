function [F,obj] = primal_constraints(epsilon,tau,data,alpha)
%PRIMAL_CONSTRAINTS Generate primal constraints
N = data.N;
xi = data.xi;
Q = data.Q;
r = data.r;
h = data.h;

T = data.T;
W = data.W;
S = data.S;
t = data.t;
M = size(T,1);
J = size(S,1);
d = size(S,2);
N_2 = size(Q,1);

Q_ = [Q];
r_ = [r];
T_ = [T];
h_ = [h];
W_ = [W];


lambda = sdpvar;
constraints = {};
constraints{end+1} = [lambda >= 0];
obj = epsilon^2*lambda;
s = sdpvar(N,1);
r_square = r_ .* r_;
psi = sdpvar(N_2,N,'full');
phi = sdpvar(N_2,N,'full');
constraints{end+1} = [s>=0];
for i=1:N
    curr_xi = xi(:,i);
    obj = obj + 1/N * s(i);
    A = [lambda*eye(d)+Q_'*diag(phi(:,i))*Q_ -0.5*T_'-Q_'*diag(phi(:,i))*W_' -lambda*curr_xi-0.5*Q_'*psi(:,i); 
        (-0.5*T_'-Q_'*diag(phi(:,i))*W_')' W_*diag(phi(:,i))*W_'+tau*eye(size(W_,1)) 0.5*(W_*psi(:,i)-h_);
        (-lambda*curr_xi-0.5*Q_'*psi(:,i))' 0.5*(W_*psi(:,i)-h_)' s(i)+alpha-(r_'*psi(:,i) - lambda*norm(curr_xi)^2 + phi(:,i)'*r_square)]; 
    NN = sdpvar(length(A));
    constraints{end+1} = [NN(:) >= 0, A - NN >= 0];
end
F = [constraints{:}];
end

