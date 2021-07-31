function [ output_args ] = run_from_to(d_iter,N_iter,from,to)
N_ = [5; 10; 20; 40; 80; 160; 320; 640];
d_ = [1; 2; 4; 8; 16; 32; 64];

for i=from:to
    filename = sprintf('results_new/%d_%d/d_%d_N_%d_i_%d_exact.mat', d_iter, N_iter, d_iter, N_iter, i);
    if exist(filename, 'file') == 2
        continue;
    end
    fprintf ('d %d N %d i %d \n', d_iter, N_iter, i);
    run_d_N_exact(d_iter,N_iter,i);
end
end

