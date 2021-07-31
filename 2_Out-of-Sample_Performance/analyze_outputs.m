clear all
N_ = [10; 20; 40; 80; 160; 320; 640; 1280];
folder_name = 'results/crossvalidation';
OBJ_Wasserstein = zeros(100,length(N_));
OBJ_Moments = zeros(100,length(N_));
OBJ_SAA = zeros(100,length(N_));
OBJ_Perfect = zeros(100,length(N_));
epsilons = zeros(100,length(N_));

curr = sprintf('perfect_information');
filename = strcat(folder_name,'/',curr);
load(filename);
OBJ_Perfect = repmat(obj_perfect_,1,length(N_));

for iterN=1:length(N_)
    N = N_(iterN);
    for i=1:100
    curr = sprintf('%d/N_%d_i_%d_Wasserstein.mat',N,N,i);
    filename = strcat(folder_name,'/',curr);
    if exist(filename, 'file') == 2
        load(filename);
        OBJ_Wasserstein(i,iterN) = obj_robust_OOS;
        epsilons(i,iterN) = epsilon;
    end
    
    curr = sprintf('%d/N_%d_i_%d_Moments.mat',N,N,i);
    filename = strcat(folder_name,'/',curr);
    if exist(filename, 'file') == 2
        load(filename);
        OBJ_Moments(i,iterN) = obj_robust_OOS;
        OBJ_SAA(i,iterN) = obj_SAA_OOS;
    end

    end
end

Improvement_Wasserstein = (OBJ_SAA./OBJ_Wasserstein-1)*100;
Improvement_Moments = (OBJ_SAA./OBJ_Moments-1)*100;

Suboptimality_Wasserstein = (OBJ_Wasserstein./OBJ_Perfect-1)*100;
Suboptimality_Moments = (OBJ_Moments./OBJ_Perfect-1)*100;
Suboptimality_SAA = (OBJ_SAA./OBJ_Perfect-1)*100;

Mean_Wasserstein = mean(Improvement_Wasserstein);
Mean_Moments =  mean(Improvement_Moments);
Mean_Suboptimality_Wasserstein = mean(Suboptimality_Wasserstein);
Mean_Suboptimality_Moments = mean(Suboptimality_Moments);
Mean_Suboptimality_SAA = mean(Suboptimality_SAA);
mean_epsilon = mean(epsilons);

q = 20;
figure(1)
set(gca,'fontsize',40)
hold on
line([0.75 8.25], [0 0],'Color','k','LineStyle',':');
upper = prctile(Improvement_Moments,100-q);
lower = prctile(Improvement_Moments,q);
for i=1:8
    line([i i], [lower(i) upper(i)],'Color','k','Linewidth',2);  
    line([i-0.2 i+0.2], [lower(i) lower(i)],'Color','k','Linewidth',2);  
    line([i-0.2 i+0.2], [upper(i) upper(i)],'Color','k','Linewidth',2);  
end
p1 = plot([1:8],Mean_Moments,'Color','b','LineWidth',7);

% plot(N_,prctile(Improvement_Moments,100-q),'Color','b','LineWidth',2,'LineStyle','--');
% plot(N_,prctile(Improvement_Moments,q),'Color','b','LineWidth',2,'LineStyle','--');
xlabel('I');
ylabel('Improvement over SAA (%)');
%set(gca,'XScale','log');
set(gca, 'XTickLabel', N_);
set(gca, 'XLim', [0.75 8.25]);
%set(gca,'xtick',1:8, 'xticklabel',{'10','20','40','80','160','320','640','1280'}) 
set(gca, 'YLim', [-30 70]);
grid on
box on 

figure(2)
set(gca,'fontsize',40)
hold on
line([0.75 8.25], [0 0],'Color','k','LineStyle',':');
upper = prctile(Improvement_Wasserstein,100-q);
lower = prctile(Improvement_Wasserstein,q);
for i=1:8
    line([i i], [lower(i) upper(i)],'Color','k','Linewidth',2);  
    line([i-0.2 i+0.2], [lower(i) lower(i)],'Color','k','Linewidth',2);  
    line([i-0.2 i+0.2], [upper(i) upper(i)],'Color','k','Linewidth',2);  
end
p2 = plot([1:8],Mean_Wasserstein,'Color','b','LineWidth',7);
%plot(N_,prctile(Improvement_Wasserstein,100-q),'Color','b','LineWidth',2,'LineStyle','--');
%plot(N_,prctile(Improvement_Wasserstein,q),'Color','b','LineWidth',2,'LineStyle','--');
xlabel('I');
ylabel('Improvement over SAA (%)');
%legend([p1 p2], 'Chebyshev DRO', 'Wasserstein DRO');
%set(gca,'XScale','log');
set(gca, 'XTickLabel', N_);
%set(gca, 'XLim', [N_(1) N_(end)]);
set(gca, 'XLim', [0.75 8.25]);
set(gca, 'YLim', [-30 70]);
grid on
box on 

figure(3)
set(gca,'fontsize',40)
hold on
line([0.75 8.25], [0 0],'Color','k','LineStyle',':');
upper = prctile(Suboptimality_SAA,100-q);
lower = prctile(Suboptimality_SAA,q);
for i=1:8
    line([i i], [lower(i) upper(i)],'Color','k','Linewidth',2);  
    line([i-0.2 i+0.2], [lower(i) lower(i)],'Color','k','Linewidth',2);  
    line([i-0.2 i+0.2], [upper(i) upper(i)],'Color','k','Linewidth',2);  
end
p1 = plot([1:8],Mean_Suboptimality_SAA,'Color','b','LineWidth',7);
%plot(N_,prctile(Suboptimality_SAA,100-q),'Color','b','LineStyle','--','LineWidth',2);
%plot(N_,prctile(Suboptimality_SAA,q),'Color','b','LineStyle','--','LineWidth',2);
xlabel('I');
ylabel('Optimality Gap (%)');
%set(gca,'XScale','log');
set(gca, 'XTickLabel', N_);
set(gca, 'XLim', [0.75 8.25]);
set(gca, 'YLim', [0 100]);
grid on
box on

figure(4)
set(gca,'fontsize',40)
hold on
line([0.75 8.25], [0 0],'Color','k','LineStyle',':');
upper = prctile(Suboptimality_Moments,100-q);
lower = prctile(Suboptimality_Moments,q);
for i=1:8
    line([i i], [lower(i) upper(i)],'Color','k','Linewidth',2);  
    line([i-0.2 i+0.2], [lower(i) lower(i)],'Color','k','Linewidth',2);  
    line([i-0.2 i+0.2], [upper(i) upper(i)],'Color','k','Linewidth',2);  
end
p2 = plot([1:8],Mean_Suboptimality_Moments,'Color','b','LineWidth',7);
%plot(N_,prctile(Suboptimality_Moments,100-q),'Color','b','LineStyle','--','LineWidth',2);
%plot(N_,prctile(Suboptimality_Moments,q),'Color','b','LineStyle','--','LineWidth',2);
xlabel('I');
ylabel('Optimality Gap (%)');
%set(gca,'XScale','log');
set(gca, 'XTickLabel', N_);
set(gca, 'XLim', [0.75 8.25]);
set(gca, 'YLim', [0 100]);
grid on
box on

figure(5)
set(gca,'fontsize',40)
hold on
line([0.75 8.25], [0 0],'Color','k','LineStyle',':');
upper = prctile(Suboptimality_Wasserstein,100-q);
lower = prctile(Suboptimality_Wasserstein,q);
for i=1:8
    line([i i], [lower(i) upper(i)],'Color','k','Linewidth',2);  
    line([i-0.2 i+0.2], [lower(i) lower(i)],'Color','k','Linewidth',2);  
    line([i-0.2 i+0.2], [upper(i) upper(i)],'Color','k','Linewidth',2);  
end
p3 = plot([1:8],Mean_Suboptimality_Wasserstein,'Color','b','LineWidth',7);
%plot(N_,prctile(Suboptimality_Wasserstein,100-q),'Color','b','LineStyle','--','LineWidth',2);
%plot(N_,prctile(Suboptimality_Wasserstein,q),'Color','b','LineStyle','--','LineWidth',2);
xlabel('I');
ylabel('Optimality Gap (%)');
%set(gca,'XScale','log');
set(gca, 'XTickLabel', N_);
set(gca, 'XLim', [0.75 8.25]);
set(gca, 'YLim', [0 100]);
grid on
box on
