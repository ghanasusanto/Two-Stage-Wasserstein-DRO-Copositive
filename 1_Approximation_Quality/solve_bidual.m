function [obj] = solve_bidual(options,epsilon,tau,data,obj)
%SOLVE_BIDUAL Generate primal constraints
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

Q_ = [Q; S];
r_ = [r; -t];
T_ = [T; zeros(J,d)];
h_ = [h; zeros(J,1)];
W_ = [W zeros(M,J); zeros(J,size(W,2)) -eye(J)];

cons = 0;
constraints = {};

for i=1:N
    curr_xi = xi(:,i);
    Y = sdpvar(d,M+J);
    gamma = sdpvar(M+J,1);
    Gamma = sdpvar(M+J);
    mu = sdpvar(d,1);
    Omega = sdpvar(d);
    obj = obj + 1/N*(trace(T_*Y) + h_'*gamma - tau*trace(Gamma));
    cons = cons + 1/N*(trace(Omega)-2*curr_xi'*mu+curr_xi'*curr_xi);
    constraints{end+1} = Q_*mu+r_==W_'*gamma;
    for j=1:N_2+J
        curr_Q = Q_(j,:)';
        curr_W = W_(:,j);
        constraints{end+1} = curr_Q'*Omega*curr_Q-2*curr_Q'*Y*curr_W+curr_W'*Gamma*curr_W==r_(j)^2;
    end
    B = [Omega Y mu; Y' Gamma gamma; mu' gamma' 1];
    constraints{end+1} = [B >= 0, B(:) >= 0];
end
constraints{end+1} = cons <= epsilon^2;
diagnostics = solvesdp([constraints{:}], -obj, options);
obj = double(obj);
diagnostics
end

