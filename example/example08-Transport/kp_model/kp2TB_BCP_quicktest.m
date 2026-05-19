%%
syms k_x k_y k_z real
syms v1 v2 v3 w real 
%%
H_kp = HK(2,2);
H_kp = H_kp...
     + Term(v1*k_x, sigma_x)...
     + Term(v2*k_y, sigma_y)...
     + Term(v3*k_z, sigma_z)...
     + Term(w*k_x^2, sigma_0)...
     + Term(w*k_y^2, sigma_0)...
     + Term(w*k_z^2, sigma_0);
%%
v1 = 1;
v2 = 1;
v3 = 1;
delta = 40e-3;
alpha = 1.2;%m* = 1.2 m_e, alpha = 1.2
w = 3.81/alpha;% % $E(k)[\mathrm{eV}] \approx \frac{3.81}{\alpha} k^2$
%%
H_kp_n = H_kp.Subsall();
H_TB = H_kp.kp2TB();
H_TB_n = H_TB.Subsall();
H_TB_n.Gen_hr('wannier90_hr.dat');
H_TB_n = HR.from_wannier90('wannier90_hr.dat');
H_TB_n.Rm = eye(3);

H_TB_n.bandplot([-0.2 0.2]);
%%
krange = 1; % just around the Gamma point
kcube_bulk = krange .* [-0.5 -0.5 0; 1 0 0; 0 1 0; 0 0 0];
NK1 = 500;%
NK2 = NK1;
% 
% [klist_cart, klist_frac, klist_k1, klist_k2, klist_k3, Grid] = kmeshgen(H_kp_n.Rm, kcube_bulk, ...
%     "Nk1", NK1, "Nk2", NK2,"full_edge",true);
    [klist_cart,~,~,~,~,Grid] = kmeshgen(H_TB_n.Rm, ...
        "Nk1",NK1, ...
        "Nk2",NK2, ...,
        "kstart",kcube_bulk(1,:), ...
        "kdir1", kcube_bulk(2,:), ...
        "kdir2", kcube_bulk(3,:), ...
        "kdir3", kcube_bulk(4,:), ...
        "full_edge",true);


%[klist_cart, klist_frac] = kcubegen3D('Rm', H_kp_n.Rm, 'KCUBE_BULK', kcube_bulk, 'nk', [NK1 NK2 1]);
% $G_{b c, n}=2 \operatorname{Re} \sum_{\ell \neq n} \frac{\left\langle v_b\right\rangle_{n \ell}\left\langle v_c\right\rangle_{\ell n}}{\left(\varepsilon_n-\varepsilon_{\ell}\right)^3}$
%% BCP
BCPCAR = BCP(H_TB_n, [1 1], klist_cart, 'ncore',1,'eps',1e-6);
%
BCPCAR_xx_1 = reshape(BCPCAR(:,1),[NK1,NK2]);
%%
X = Grid(:,:,1);
Y = Grid(:,:,2);
Z = Grid(:,:,3);

syms k_x k_y k_z real;
k = sqrt(k_x.^2+k_y.^2+k_z.^2);
Gxx = -1/4 * (k.^2.*k_y.^2 + k_x.^2.*k_z.^2)./(k.^5.*(k_x.^2+k_y.^2));
Gxx_0 = subs(Gxx,k_z,0);
Gxx_0_function = matlabFunction(Gxx_0, 'Vars', [k_x, k_y, k_z]);
G_0_xx_th = Gxx_0_function(X,Y,Z);

% 
% xlabel("k_x");
% ylabel("k_y");
% title(['G_{xx}, N:',num2str(NK1),'*',num2str(NK2)]);
% xlim(krange*[-0.5,0.5]);
% ylim(krange*[-0.5,0.5]);
% colormap("cool")
%%
colorcut = 0.000001;
minValue = min(G_0_xx_th,[],"all");
% clim(axs(1),[minValue*colorcut 0]);
G_0_xx_th(G_0_xx_th<minValue*colorcut) = minValue*colorcut;
BCPCAR_xx_1(BCPCAR_xx_1<minValue*colorcut) = minValue*colorcut;
%%
[fig,axs] = Figs(1,2);
set(fig ,'Position', [0, 0, 1440, 1080]);
% 添加整体标题
sgtitle('Berry Connection Polarizability (BCP) Comparison', 'FontSize', 14, 'FontWeight', 'bold');
% 创建子图并绘制等高线图

contourf(axs(1),X, Y, G_0_xx_th,'LineStyle', 'none');
title(axs(1),['Theoretical $G_0^{xx}$',', N:',num2str(NK1),'*',num2str(NK2)], 'Interpreter', 'latex');
colorbar(axs(1));
axis(axs(1),'equal');



contourf(axs(2),X, Y, BCPCAR_xx_1,'LineStyle', 'none'); 
title(axs(2),['Numerical $G_0^{xx}$',', N:',num2str(NK1),'*',num2str(NK2)], 'Interpreter', 'latex');
colorbar(axs(2));
axis(axs(2),'equal');
% clim(axs(2),[minValue*colorcut 0]);