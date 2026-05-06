function [Out,SymMat,SymMatDisplay,n_indp] = pretty(Tensor,mode,options)
arguments
    Tensor
    mode  {mustBeMember(mode,{'EQ','latex','Table'})} =  'EQ';
    options.silent logical = false;
end
n_ele = Tensor.n_ele;
tensor_rank = Tensor.tensor_rank;

tensor_symmed_mat = reshape(Tensor.tensor0, n_ele, n_ele);
[tensor_symmed_reduce, independent_ele] = rref(tensor_symmed_mat');
tensor_names = sym('chi_', ones(1,tensor_rank)*3, 'real');
tensor_names_vector = reshape(tensor_names, n_ele, 1);
n_indp = length(independent_ele);
switch mode
    case 'EQ'
        if ~options.silent
            fprintf('Independent Elements: %d\n',n_indp);
        end
        for i = 1:n_ele
            if ~all(tensor_symmed_reduce(1:n_indp,i)==0)
                eq = sum(tensor_symmed_reduce(1:n_indp,i) .* tensor_names_vector(independent_ele));
                disp(string(tensor_names_vector(i))+ ...
                    " = "+ ...
                    string(eq))
            end
        end
        Out = '';
        SymMat = '';
    case 'latex'
         Table = pretty(Tensor,'Table');
         %
         Out = table2latex(Table);
         
    case 'Table'
        if ~options.silent
            fprintf('Independent Elements: %d\n',n_indp);
        end
        tensor_names_vector_reduce = tensor_names_vector;
         for i = 1:n_ele
            if ~all(tensor_symmed_reduce(1:n_indp,i)==0)
                eq = sum(tensor_symmed_reduce(1:n_indp,i) .* tensor_names_vector(independent_ele));
                tensor_names_vector_reduce(i) = simplify(eq);
            else
                tensor_names_vector_reduce(i) = sym(0);

            end
        end
         tensor_names_vector_Str = string(tensor_names_vector);
         tensor_names_vector_StrL  = split(tensor_names_vector_Str,{'chi_','_'});
         tensor_names_vector_StrL(:,1) =[];
         Num_Mat = double(tensor_names_vector_StrL);
         % Choose How to sort
         RankList = 1:tensor_rank;
         SelectCols = RankList(1:floor(tensor_rank/2));
         Nselect = numel(SelectCols);
         [Num_Mat_sort,seqL] = sortrows(Num_Mat,SelectCols);
         tensor_names_vector_reduce_reseq = tensor_names_vector_reduce(seqL);

         NRows = 3^(tensor_rank-Nselect);
         NCols = 3^Nselect;
         tensor_names_vector_reduce_mat = reshape(tensor_names_vector_reduce_reseq,[ NRows NCols]).';

         Num_Mat_sort_str_Rows = num2str(Num_Mat_sort(:,Nselect+1:end));
         Num_Mat_sort_str_Cols = num2str(Num_Mat_sort(1:NRows:n_ele,1:Nselect));

         % 转换为 cell 数组（每行一个 cell）
         cellArray = cellstr(Num_Mat_sort_str_Rows(1:NRows,:));
         % 逐行去除空格（逻辑索引方法）
         RowsCellNoSpaces = cellfun(@(x) x(x ~= ' '), cellArray, 'UniformOutput', false);

         VariableNames = [{'Tensor'};RowsCellNoSpaces];

         % 转换为 cell 数组（每行一个 cell）
         cellArray = cellstr(Num_Mat_sort_str_Cols(1:NCols,:));
         % 逐行去除空格（逻辑索引方法）
         ColsCellNoSpaces = cellfun(@(x) x(x ~= ' '), cellArray, 'UniformOutput', false);
         
         A = sym(ColsCellNoSpaces);
         Out  = array2table( ...
             [A,tensor_names_vector_reduce_mat],'VariableNames',VariableNames);
         SymMat = tensor_names_vector_reduce_mat;
         symTitle = [str2sym('Tensor'), sym(RowsCellNoSpaces).'];
         symCol = sym(ColsCellNoSpaces);
         SymMatDisplay = [symTitle;[symCol,SymMat]];
    otherwise
         return;
end
end