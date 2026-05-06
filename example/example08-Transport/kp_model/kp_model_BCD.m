%% useful_tools
clear
useful_matrices(["sigma","tau"])
%%
syms k_x k_y k_z real
syms v t m alpha eta real 
%%
H_kp = HK(2,2);
H_kp = H_kp...
     + Term(t*k_x, sigma_0)...
     + Term(v*k_y, sigma_x)...
     + Term(v*eta*k_x, sigma_y)...
     + Term(m/2 - alpha*(k_x^2+k_y^2), sigma_z);
%%
v = 1;
t = 0.5;
alpha = 1;
eta = -1;
m = 0.2;
%%
H_kp_n = H_kp.Subsall();
H_kp_n.bandplot([-0.2 0.2]);
%% BCD
BCDCAR = BCD_k(H_kp_n, [1 2 2], [0 0.1 0], [-0.5:0.1:0.5],'eps',1e-6);