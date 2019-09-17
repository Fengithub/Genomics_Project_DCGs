dat_normal = csvread('tcga_normal_refined_imputed_subdata4000.csv',1,1);
dat_cancer = csvread('tcga_cancer_refined_imputed_subdata4000.csv',1,1);
dat_normal = dat_normal';
dat_cancer = dat_cancer';
cov_normal = cov(dat_normal);
cov_cancer = cov(dat_cancer);
rho = 0.1;
maxIt = 2;
tol = 10e6;
glasso_normal = graphicalLasso(cov_normal, rho, maxIt, tol);
% glasso_cancer = graphicalLasso(cov_cancer, rho, maxIt, tol);
save W_normal W
save Theta_normal Theta

