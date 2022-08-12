function [P2, W] = clus_fin_update(rho, lam, lam2, eta, c, X, true_labs)



[n,p]=size(X);  CCC=max(true_labs);

%% Step 1: contruct multiple doubly stochastic similarity matrices using Gaussian kernels.
[Kernels, id]= func_doubly(X);

%% Step 2: Obtain the intermediae target matrix
[~, W, P, ~]=clus_sim_update2_2(CCC, c, rho, lam, id, eta, single(Kernels)) ;

%% Step 3: Obtain the final target matrix
[V, temp, evs]=eig1(P, CCC);
[~, P2, ~, ~]=clus_sim_update0_3(CCC, lam2, eta, V*V');

end





