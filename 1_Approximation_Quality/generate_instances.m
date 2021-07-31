clc
clear all
rng('default');

N_ = [5; 10; 20; 40; 80; 160; 320; 640];
d_ = [1; 2; 4; 8; 16; 32; 64];
for d_iter = 1:length(d_)
    d = d_(d_iter);
    for i=1:100
        L = unidrnd(ceil(log(d+1)));
        A{d_iter,i} = rand(L,d);
        b{d_iter,i} = rand(L,1).*sum(A{d_iter,i},2);
    end
end
save('inputs', 'A', 'b');
return
