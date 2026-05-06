function [Lx, Ly, Lz] = Gen_LxLyLz(hr, options)
arguments
    hr = []
    options.fast logical = true
    options.write_to_QtransF logical = true
end
if isempty(hr)
    [orbL, elementL, quantumL, SpinfulFlag] = wout_read();
else
    orbL = hr.orbL;
    elementL = hr.elementL;
    quantumL = hr.quantumL;
    if sum(quantumL(:,4)) == 0 && sum(abs(quantumL(:,4))) ~= 0
        SpinfulFlag = true;
    else
        SpinfulFlag = false;
    end
end
%%
if options.fast
    L_p(1,:) = Y_lm(1,  0);
    L_p(2,:) = Y_lm(1,  1);
    L_p(3,:) = Y_lm(1, -1);

    Lx_p = (L_p.Lplus+L_p.Lminus)/sym(2);
    Ly_p = (L_p.Lplus-L_p.Lminus)/sym(2i);
    Lz_p =  L_p.Lz;

    L_d(1,:) = Y_lm(2,  0);
    L_d(2,:) = Y_lm(2,  1);
    L_d(3,:) = Y_lm(2, -1);
    L_d(4,:) = Y_lm(2,  2);
    L_d(5,:) = Y_lm(2, -2);

    Lx_d = (L_d.Lplus+L_d.Lminus)/sym(2);
    Ly_d = (L_d.Lplus-L_d.Lminus)/sym(2i);
    Lz_d =  L_d.Lz;

    change_site = ~(diff(orbL(:,1))==0 & diff(orbL(:,2))==0 & diff(orbL(:,3))==0);
    change_L    = diff(quantumL(:,1))~=0 | diff(quantumL(:,2))~=0;
    change_points = find(change_site | change_L) + 1;
    starts = [1; change_points];
    
    Lx_all = [];
    Ly_all = [];
    Lz_all = [];
    for i = 1:length(starts)
        switch quantumL(starts(i), 2)
            case 0
                Lx = 0;
                Ly = 0;
                Lz = 0;
            case 1
                Lx = Lx_p;
                Ly = Ly_p;
                Lz = Lz_p;
            case 2
                Lx = Lx_d;
                Ly = Ly_d;
                Lz = Lz_d;
        end
    
        Lx_all = blkdiag(Lx_all, Lx);
        Ly_all = blkdiag(Ly_all, Ly);
        Lz_all = blkdiag(Lz_all, Lz);
    end
    Lx = Lx_all;
    Ly = Ly_all;
    Lz = Lz_all;

    if SpinfulFlag
        I2 = sym(eye(2));
    
        Lx = kron(Lx, I2);
        Ly = kron(Ly, I2);
        Lz = kron(Lz, I2);
    end
else
    ZeroMat = zeros(size(orbL,1));
    for i = 1:size(orbL,1)
        for j = 1:size(orbL,1)
            ZeroMat(i,j) = ~all(orbL(i,:)==orbL(j,:));
        end
    end
    
    EqnumL = [elementL,quantumL];
    [EqnumLsort,sort_label] = sortrows(EqnumL,5,"descend");
    M  = TBkit.P2M(sort_label); 
    EQup = EqnumLsort(EqnumLsort(:,end) > 0,:);
    
    clear L
    for i = 1:size(EQup,1)
        ielement = EQup(i,1);
        %assume(symvar,'real');
        L(i,:) = Y_lm(EQup(i,3), EQup(i,4) ,'An',ielement);
    end
    
    if SpinfulFlag
        I2 = sym(eye(2));
    
        Lx = kron(I2, (L.Lplus+L.Lminus)/sym(2));
        Ly = kron(I2, (L.Lplus-L.Lminus)/sym(2i));
        Lz = kron(I2, (L.Lz));
    else
        Lx = (L.Lplus+L.Lminus)/sym(2);
        Ly = (L.Lplus-L.Lminus)/sym(2i);
        Lz = (L.Lz);
    end
    Lx = M'*Lx/M';
    Ly = M'*Ly/M';
    Lz = M'*Lz/M';
    % 
    Lx(logical(ZeroMat)) = 0;
    Ly(logical(ZeroMat)) = 0;
    Lz(logical(ZeroMat)) = 0;
end
%% write to QtransF, which uses inner uudd basis
if options.write_to_QtransF
    WAN_NUM = length(Lx);
    
    udud2uudd = [1:2:(WAN_NUM-1),2:2:(WAN_NUM)];
    Lx_write = Lx(udud2uudd,udud2uudd);
    Ly_write = Ly(udud2uudd,udud2uudd);
    Lz_write = Lz(udud2uudd,udud2uudd);

    Lx_line = reshape(double(Lx_write), numel(Lx_write), 1);
    Ly_line = reshape(double(Ly_write), numel(Ly_write), 1);
    Lz_line = reshape(double(Lz_write), numel(Lz_write), 1);
    
    writematrix([real(Lx_line), imag(Lx_line)], "Lx.dat", 'Delimiter',' ')
    writematrix([real(Ly_line), imag(Ly_line)], "Ly.dat", 'Delimiter',' ')
    writematrix([real(Lz_line), imag(Lz_line)], "Lz.dat", 'Delimiter',' ')
end
end
