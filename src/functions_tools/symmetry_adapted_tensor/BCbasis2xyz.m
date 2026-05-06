function R_abc2xyz = BCbasis2xyz(MSG_number1)
%% 14 Bravais lattices, defined by The Mathematical Theory of Symmetry in Solids Representation Theory for Point Groups and Space Groups (Christopher Bradley, Arthur Cracknell)
if     ismember(MSG_number1, 1:2) % triclinic
    R_abc2xyz = eye(3);
elseif ismember(MSG_number1, [3,4,6,7,10,11,13,14]) % monoclinic-P
    R_abc2xyz = eye(3);
elseif ismember(MSG_number1, [5,8,9,12,15]) % monoclinic-C
    R_abc2xyz = [1/sqrt(2),  1/sqrt(2), 0; 1/sqrt(2), -1/sqrt(2), 0;  0, 0, 1];
elseif ismember(MSG_number1, [16:19,25:34,47:62]) % orthorhombic-P
    R_abc2xyz = [[0, 1, 0]; [-1, 0, 0]; [0, 0, 1]];
elseif ismember(MSG_number1, [20:21,35:41,63:68]) % orthorhombic-C
    R_abc2xyz = [[1/2, 1/2, 0]; [-1/2, 1/2, 0]; [0, 0, 1]];
elseif ismember(MSG_number1, [23:24,44:46,71:74]) % orthorhombic-I
    R_abc2xyz = [[1/2, -1/2, 1/2]; [1/2, -1/2, -1/2]; [1/2, 1/2, -1/2]];
elseif ismember(MSG_number1, [22,42:43,69:70]) % orthorhombic-F
    R_abc2xyz = [[1/2, 0, 1/2]; [0, -1/2, -1/2]; [1/2, 1/2, 0]];
elseif ismember(MSG_number1, [75:78,81,83:86,89:96,99:106,111:118,123:138]) % Tetragonal-P
    R_abc2xyz = eye(3);
elseif ismember(MSG_number1, [79:80,82,87:88,97:98,107:110,119:122,139:142]) % Tetragonal-I  
    R_abc2xyz = [[-1/2, 1/2, 1/2]; [1/2, -1/2, 1/2]; [1/2, 1/2, -1/2]];
elseif ismember(MSG_number1, [143:145,147,149:154,156:159,162:165]) % trigonal-P, same as hexagonal
    R_abc2xyz = [[0, sym(sqrt(3))/2, 0]; [-1, 1/2, 0]; [0, 0, 1]];
elseif ismember(MSG_number1, [146,148,155,160,161,166,167]) % trigonal-R
    R_abc2xyz = [[0, sym(sqrt(3))/2, -sym(sqrt(3))/2]; [-1, 1/2, 1/2]; [1, 1, 1]];
elseif ismember(MSG_number1, 168:194) % hexagonal
    R_abc2xyz = [[0, sym(sqrt(3))/2, 0]; [-1, 1/2, 0]; [0, 0, 1]];
elseif ismember(MSG_number1, [195,198,200,201,205,207,208,212,213,215,218,221:224]) % Cubic-P
    R_abc2xyz = eye(3);
elseif ismember(MSG_number1, [196,202,203,209,210,216,219,225:228]) % Cubic-F
    R_abc2xyz = [[0, 1/2, 1/2]; [1/2, 0, 1/2]; [1/2, 1/2, 0]];
elseif ismember(MSG_number1, [197,199,204,206,211,214,217,220,229,230]) % Cubic-I  
    R_abc2xyz = [[-1/2, 1/2, 1/2]; [1/2, -1/2, 1/2]; [1/2, 1/2, -1/2]];
end