function [F,obj] = primal_constraints_moments_bounds_new(data,alpha,gamma1,gamma2)
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

mu = mean(data.xi')';
Sigma = (data.xi-repmat(mu,1,N))*(data.xi-repmat(mu,1,N))'/N;
constraints = {};
s = sdpvar;
t = sdpvar;
m = sdpvar(d,1);
M = sdpvar(d);
r_square = r_ .* r_;
psi = sdpvar(N_2,1,'full');
phi = sdpvar(N_2,1,'full');
obj = s + t;
A = [M -0.5*T_' 0.5*m;
    (-0.5*T_')' W_*diag(phi)*W_' 0.5*(W_*psi-h_);
    (0.5*m)' 0.5*(W_*psi-h_)' s+alpha-(r_'*psi + phi'*r_square)];
NN = sdpvar(length(A));
constraints{end+1} = [NN(:) >= 0, A - NN >= 0];
B = [M 0.5*m; 0.5*m' s];
NN = sdpvar(length(B));
constraints{end+1} = [NN(:) >= 0, B - NN >= 0];
constraints{end+1} = [t>=trace((gamma2*Sigma+mu*mu')*M)+m'*mu+sqrt(gamma1)*norm(Sigma^0.5*(m+2*M*mu))];
F = [constraints{:}];
end

