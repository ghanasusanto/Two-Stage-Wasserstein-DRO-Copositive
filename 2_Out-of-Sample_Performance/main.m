clc
clear all
N_ = [5; 10; 20; 40; 80; 160; 320; 640; 1280];
[beta,cost,budget,holding_cost,shortage_cost] = generate_data();
run_perfect_information(beta,cost,budget,holding_cost,shortage_cost);
return
for N=N_'
    run_vs_moment_model_N(N,N_,beta,cost,budget,holding_cost,shortage_cost);
    run_N(N,N_,beta,cost,budget,holding_cost,shortage_cost);
end

return

