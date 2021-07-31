function [F,obj] = SAA_constraints(data)
%SAA_CONSTRAINTS Generate SAA constraints
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

obj = 0;
F = [];
constraints = {};
T_ = kron(speye(N),T);
xi_ = reshape(xi,N*d,1);
h_ = repmat(h,N,1);
W_ = kron(speye(N),W);
r_ = repmat(r,N,1);
Y = sdpvar(N_2,N,'full');
Y_ = reshape(Y,N*N_2,1);
obj = 1/N*r_'*Y_;
%for i=1:N
%    curr_xi = xi(:,i);
%    obj = obj + 1/N*((Q*curr_xi+r)'*Y(:,i));
    %constraints{end+1} = T*curr_xi+h <= W*Y(:,i);
%end
constraints{end+1} = T_*xi_+h_ <= W_*Y_;
F = [constraints{:}];
end

