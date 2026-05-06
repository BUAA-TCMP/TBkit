function [Rm_new, R] = standardize_Rm(Rm_old)
% ALIGN_LATTICE_TO_STANDARD
% 输入 L: 3x3 矩阵，每行为晶格矢量 a, b, c（如 VASP POSCAR）
% 输出 L_new: 旋转后的 3x3 矩阵，满足 a 沿 x 轴，b 在 xy 平面

    % 提取列向量（转置以便按列处理）
    a = Rm_old(1,:)';  % 列向量
    b = Rm_old(2,:)';
    c = Rm_old(3,:)';

    % Step 1: e1 = a / |a|
    e1 = a / norm(a);

    % Step 2: 从 b 中减去 e1 方向分量，得到垂直于 e1 的部分
    b_perp = b - dot(b, e1) * e1;
    if norm(b_perp) < 1e-12
        error('Lattice vector b is collinear with a. Cannot define unique xy plane.');
    end
    e2 = b_perp / norm(b_perp);

    % Step 3: e3 = e1 × e2 （右手系）
    e3 = cross(e1, e2);
    e3 = e3 / norm(e3);  % 数值稳定性

    % 构造旋转矩阵 R（将原坐标系向量变换到新坐标系）
    R = [e1'; e2'; e3'];  % 每行是一个新基矢在原坐标系中的表示
    % 注意：R 是正交矩阵，R * v 就是 v 在新坐标系下的坐标

    % 应用旋转到所有晶格矢量（作为列）
    M_new = R * Rm_old';   % 旋转后的列向量

    % 转换回 VASP 行格式
    Rm_new = M_new';       % 每行是一个晶格矢量

    tol = 1e-6;
    Rm_new(abs(Rm_new) < tol) = 0;
    R(abs(R) < tol) = 0;
end