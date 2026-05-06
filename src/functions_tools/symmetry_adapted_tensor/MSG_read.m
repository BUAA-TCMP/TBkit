function gen_list_xyz = MSG_read(BNS_number, options)
% read generators of a magnetic space group from the inner database by specifying the BNS number
% the output generators are all in Cartesian coordinate system 
arguments
    BNS_number (1,2) double
    options.read_translation logical = false
end
%% inner database
load("Table_of_Magnetic_Sym.mat", "Table_of_Magnetic_Sym");
id = Table_of_Magnetic_Sym.mag_space_group == BNS_number(1) & Table_of_Magnetic_Sym.number2 == BNS_number(2);
Table_of_Magnetic_Sym_check = Table_of_Magnetic_Sym(id,:);
oper_num = Table_of_Magnetic_Sym_check.Number_of_Symmetry_operations;

R_abc2xyz = BCbasis2xyz(BNS_number(1));
%%
gen_list_xyz = repmat(Oper(), 1, oper_num);
for i = 1:oper_num

    oper_i_list_form = Table_of_Magnetic_Sym_check.symmetry_operation_R{1}(i,:);
    oper_i_isUnitary = Table_of_Magnetic_Sym_check.unitary_antiunitary{1}(i)==1;
    
    oper_i_R = reshape(oper_i_list_form(1:9),[3,3])';
    oper_i_R = R_abc2xyz * oper_i_R * R_abc2xyz^(-1);

    if options.read_translation
        oper_i_t = Table_of_Magnetic_Sym_check.symmetry_operation_t{1}(i,:);
    else
        oper_i_t = [0 0 0];
    end
    oper_i = Oper(double(oper_i_R), NaN, oper_i_t, "conjugate",~oper_i_isUnitary);
    gen_list_xyz(i) = oper_i;
end
end