%[text] # Symmetry Adapted Tensor
%[text] ## Step1: Getting the symmetry operators
%[text] ### Specifying the operators manually
% P = Oper.inversion();
% Mx = Oper.mirror([1,0,0]);
% My = Oper.mirror([0,1,0]);
% Mz = Oper.mirror([0,0,1]);
% C6x = Oper.rotation(1/6,[1,0,0]);
% C6z = Oper.rotation(1/6,[0,0,1]);
% C3x = Oper.rotation(1/3,[1,0,0]);
% C3z = Oper.rotation(1/3,[0,0,1]);
% S6x = Mx*C6x;
% S6z = Mz*C6z;
%%
%[text] ### Using our inner database
Gen_list_xyz1 = MSG_read([63, 457]);
%%
%[text] ### Reading from a mcif format file which has non-standard BNS magnetic unit cell
cif_struct = cif_read("MnTe_non_std.mcif");
Gen_list_xyz2 = cif_struct.symop_xyz;
easy_axis = cif_struct.easy_axis_xyz; % we have assumed it is a collinear magnet

% rotate the easy axis and the Gen_list_xyz to the x-direction to check the results
Raxis_inplane = sym([easy_axis(1:2)'; -easy_axis(2), easy_axis(1)]);
Raxis = Oper(blkdiag(Raxis_inplane, 1));
for i = 1:length(Gen_list_xyz2)
    Gen_list_xyz2(i) = Raxis * Gen_list_xyz2(i) * Raxis^(-1);
end
%%
%[text] ## Step2: Defining the Jahn symbol of the tensor
% AHC_1st     = 'a{V2}';
% AHC_2nd_QMD = 'a{V2}V';
% AHC_2nd_BCD = '{V2}V';
% SHC         = 'eV3';

jahn_symbol_Str = 'eV3';
Tensor = jahn_symbol(jahn_symbol_Str);

for i = 1:length(Gen_list_xyz1)
    Tensor = group_transformation(Tensor, Gen_list_xyz1(i));
end

[~, ~, SymMatDisplay] = pretty(Tensor,"Table"); %[output:7ea1fcb1]
SymMatDisplay %[output:16c11ae5]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:7ea1fcb1]
%   data: {"dataType":"text","outputData":{"text":"Independent Elements: 6\n","truncated":false}}
%---
%[output:16c11ae5]
%   data: {"dataType":"symbolic","outputData":{"name":"SymMatDisplay","value":"\\left(\\begin{array}{cccccccccc}\n\\mathrm{Tensor} & 11 & 21 & 31 & 12 & 22 & 32 & 13 & 23 & 33\\\\\n1 & 0 & 0 & 0 & 0 & 0 & \\chi_{1,3,2}  & 0 & \\chi_{1,2,3}  & 0\\\\\n2 & 0 & 0 & \\chi_{2,3,1}  & 0 & 0 & 0 & \\chi_{2,1,3}  & 0 & 0\\\\\n3 & 0 & \\chi_{3,2,1}  & 0 & \\chi_{3,1,2}  & 0 & 0 & 0 & 0 & 0\n\\end{array}\\right)"}}
%---
