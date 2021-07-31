clear all

N_ = [5; 10; 20; 40; 80; 160; 320; 640];
d_ = [1; 2; 4; 8; 16; 32; 64];
folder_name = 'results_new';
load('inputs');
for d_iter=1:length(d_)
    OBJ_Copositive{d_iter} = zeros(100,length(N_));
    OBJ_Exact{d_iter} = zeros(100,length(N_));
    SOLVETIME_Copositive{d_iter} = zeros(100,length(N_));
    OBJ_Dec_Rules{d_iter} = zeros(100,length(N_));
    SOLVETIME_Dec_Rules{d_iter} = zeros(100,length(N_));
    Analytical_Value{d_iter} = zeros(100,length(N_));
    for N_iter=1:length(N_)
        for i=1:100
%             if (0.90*sum(a{d_iter,i}) <= b{d_iter,i})
%                 continue;
%             end
            curr = sprintf('%d_%d/d_%d_N_%d_i_%d.mat', d_iter, N_iter, d_iter, N_iter, i);
            filename = strcat(folder_name,'/',curr);
            if exist(filename, 'file') == 2
                load(filename);
                OBJ_Copositive{d_iter}(i,N_iter) = obj_robust;
                SOLVETIME_Copositive{d_iter}(i,N_iter) = diagnostics_robust.solvertime;
                if (diagnostics_robust.problem == 4)
                    OBJ_Copositive{d_iter}(i,N_iter) = NaN;
                end
                OBJ_Dec_Rules{d_iter}(i,N_iter) = obj_dec_rules;
                SOLVETIME_Dec_Rules{d_iter}(i,N_iter) = diagnostics_dec_rules.solvertime;
%                 if (diagnostics_dec_rules.problem == 4)
%                     OBJ_Dec_Rules{d_iter}(i,N_iter) = NaN;
%                 end

%                Analytical_Value{d_iter}(i,N_iter) = analytical_value;
                
%             end
%             curr = sprintf('%d_%d/d_%d_N_%d_i_%d_exact.mat', d_iter, N_iter, d_iter, N_iter, i);
%             filename = strcat(folder_name,'/',curr);
%             if exist(filename, 'file') == 2
%                 load(filename);
                OBJ_Exact{d_iter}(i,N_iter) = obj_exact;
            end
        end
    end
    Improvement{d_iter} = (OBJ_Dec_Rules{d_iter}-OBJ_Copositive{d_iter})./(OBJ_Copositive{d_iter})*100;
    Gap_Copositive{d_iter} = (OBJ_Copositive{d_iter}-OBJ_Exact{d_iter})./(OBJ_Exact{d_iter})*100;
    Gap_Dec_Rules{d_iter} = (OBJ_Dec_Rules{d_iter}-OBJ_Exact{d_iter})./(OBJ_Exact{d_iter})*100; 
end
for d_iter=1:length(d_)
    solved{d_iter} = 1-isnan(Improvement{d_iter});
    solved_Dec_Rules{d_iter} =  1-isnan(OBJ_Dec_Rules{d_iter});
    percent_solved{d_iter} = (mean(solved{d_iter}))*100;
    Improvement{d_iter}(isnan(Improvement{d_iter})) = 0;
    solved_gap_copositive{d_iter} = 1-isnan(Gap_Copositive{d_iter});
    Gap_Copositive{d_iter}(isnan(Gap_Copositive{d_iter})|isinf(Gap_Copositive{d_iter})) = 0;
    solved_gap_dec_rules{d_iter} = 1-isnan(Gap_Dec_Rules{d_iter});
    Gap_Dec_Rules{d_iter}(isnan(Gap_Dec_Rules{d_iter})|isinf(Gap_Dec_Rules{d_iter})) = 0;
    SOLVETIME_Copositive{d_iter}(solved{d_iter}==0) = 0;
    SOLVETIME_Dec_Rules{d_iter}(isnan(SOLVETIME_Dec_Rules{d_iter})) = 0;
    
    Mean_Improvement{d_iter} = sum(Improvement{d_iter})./percent_solved{d_iter};
    Mean_Solvetime_Copositive{d_iter} = sum(SOLVETIME_Copositive{d_iter})./percent_solved{d_iter};
    Mean_Solvetime_Dec_Rules{d_iter} = sum(SOLVETIME_Dec_Rules{d_iter})./sum(solved_Dec_Rules{d_iter});
    Mean_Gap_Copositive{d_iter} = sum(Gap_Copositive{d_iter})./sum(solved_gap_copositive{d_iter});
    Mean_Gap_Dec_Rules{d_iter} = sum(Gap_Dec_Rules{d_iter})./sum(solved_gap_dec_rules{d_iter});
end

return

d_ = 1;
q = 5;
figure(1)
set(gca,'fontsize',20)
hold on
line([N_(1) N_(end)], [0 0],'Color','k','LineStyle',':');
p1 = plot(N_,Mean_Improvement{d_},'Color','b','LineWidth',5,'LineStyle','--');
plot(N_,prctile(Improvement{d_},100-q),'Color','b','LineWidth',2,'LineStyle','--');
plot(N_,prctile(Improvement{d_},q),'Color','b','LineWidth',2,'LineStyle','--');
xlabel('I');
ylabel('Improvement over Decision Rules Approximation (%)');
set(gca,'XScale','log');
set(gca, 'XTick', N_);
set(gca, 'XLim', [N_(1) N_(end)]);
grid on
box on


return
