function [Energy_list, node_index_list, klist_s_list] = find_node(occupi_band_index, EIGENCAR, node_cut, klist_s)
    % FIND_NODE - Identifies the nodes between consecutive bands where the energy difference is below a given threshold.
    %
    % Syntax:
    %   [Energy_list, node_index_list, klist_s_list] = find_node(occupi_band_index, EIGENCAR, node_cut, klist_s)
    %
    % Inputs:
    %   occupi_band_index  - The index of the occupied band. Default: norb/2
    %   EIGENCAR           - Eigenvalues of the system (matrix). Default: EIGENVAl_read()
    %   node_cut           - Energy difference threshold to identify nodes. Default: 0.0001
    %   klist_s            - List of k-points. Default: generated by kpathgen3D()
    %
    % Outputs:
    %   Energy_list        - List of energy values at the nodes
    %   node_index_list    - Indices of the nodes
    %   klist_s_list       - K-points corresponding to the nodes
    
    % Default arguments if not provided
    [norb, ~] = size(EIGENCAR);
    if nargin < 4
        POSCAR_read;
        [~, ~, klist_s, ~, ~] = kpathgen3D(Rm);  % Generate k-path if klist_s is not provided
    end
    if nargin < 3
        node_cut = 0.0001;  % Default node_cut
    end
    if nargin < 2
        EIGENCAR = EIGENVAl_read();  % Read eigenvalues if not provided
    end
    if nargin < 1
        occupi_band_index = norb / 2;  % Default band index if not provided
    end
    
    % Compute the energy differences between consecutive bands
    E_dif_list = EIGENCAR(occupi_band_index + 1, :) - EIGENCAR(occupi_band_index, :);
    
    % Initialize lists for the nodes
    node_index_list = [];
    Energy_list = [];
    klist_s_list = [];
    
    % Loop over all k-points and check for nodes
    for i = 1:length(E_dif_list)
        if E_dif_list(i) < node_cut  % If the energy difference is below the threshold
            node_index_list = [node_index_list, i];
            Energy_list = [Energy_list, (EIGENCAR(occupi_band_index + 1, i) + EIGENCAR(occupi_band_index, i)) / 2];  % Average energy at the node
            klist_s_list = [klist_s_list; klist_s(i, :)];  % Store the corresponding k-point
        end
    end
end
