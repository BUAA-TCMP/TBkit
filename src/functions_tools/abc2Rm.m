function Rm = abc2Rm(a, b, c, alpha, beta, gamma)
% LATTICE_VECTORS 将晶格参数 (a,b,c,alpha,beta,gamma) 转换为直角坐标系下的 3x3 晶格矩阵
% 输入：
%   a, b, c : 晶格常数（长度单位：埃）
%   alpha, beta, gamma : 晶轴夹角（单位：度）
% 输出：
%   lattice       : 3x3 矩阵，每一行为一个晶格矢量（与 VASP POSCAR 格式一致）

    % 第一个基矢沿 x 轴
    v1 = [a, 0, 0];

    % 第二个基矢在 xy 平面内
    v2_x = b * cosd(gamma);
    v2_y = b * sind(gamma);
    v2 = [v2_x, v2_y, 0];

    % 第三个基矢的 z 分量由体积和几何关系确定
    % 使用三斜晶系标准公式
    % 参考：https://en.wikipedia.org/wiki/Fractional_coordinates#Conversion_to_Cartesian_coordinates

    v3_x = c * cosd(beta);
    v3_y = c * (cosd(alpha) - cosd(beta)*cosd(gamma)) / sind(gamma);
    v3_z = c * sqrt(1 - cosd(beta)^2 - ((cosd(alpha) - cosd(beta)*cosd(gamma))/sind(gamma))^2 );

    v3 = [v3_x, v3_y, v3_z];

    % 组合成 3x3 矩阵（每行是一个晶格矢量，与 VASP 一致）
    Rm = [v1; v2; v3];
end