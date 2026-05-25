%% Example: Constructing and plotting surface states from Wannier Hamiltonian
% This script demonstrates how to use Wannier functions to build a
% tight-binding Hamiltonian, compare with DFT results, and visualize
% both bulk and surface band structures and Fermi Surface.

%% Load Wannier Hamiltonian
% Copy the symmetrized Wannier Hamiltonian file to the working file

HR_TB = HR.from_wannier90('wannier90_hr.dat');   % Build Hamiltonian from Wannier90
HR_TB = HR_TB < 'POSCAR-rhobuhedrical-12';       % Attach lattice and atomic positions
HR_TB = HR_TB < 'KPOINTS_NP';                    % Attach k-point sampling

%% Resequence orbitals
HR_TB_bk = HR_TB;
HR_TB = HR_TB_bk.reseq([1 3 5 7 9 11 2 4 6 8 10 12]); % Reorder Wannier orbitals

%% Plot Wannier band structure
Efermi = 0;                                     % Fermi level (can shift as needed)
EIGENCAR = HR_TB.EIGENCAR_gen() - Efermi;       % Generate eigenvalues
[fig,ax] = bandplot(EIGENCAR,[-2,2], ...
    'title','wannier-dft','Color','b', ...
    'POSCAR','POSCAR - rhobuhedrical','KPOINTS','KPOINTS_NP');

%% Compare with DFT band structure
!cp EIGENVAL-primitive EIGENVAL;                % Load primitive cell EIGENVAL
EIGENCAR_DFT = EIGENVAL_read('vasp','EIGENVAL'); 
Efermi = 6.26717880;                            % DFT Fermi energy
bandplot(EIGENCAR_DFT-Efermi,[-5,10], ...
    'title',"Bismuth-DFT-vs-Model",'Color','r', ...
    'POSCAR','POSCAR - rhobuhedrical','KPOINTS','KPOINTS_NP','ax',ax);

%% Build supercell Hamiltonian
tran1 = [-1 1 0; 0 -1 1; 1 1 1];                % Supercell transformation matrix
HR_TBsup = HR_TB.rewrite.supercell_hr(tran1);   % Build supercell Hamiltonian

%% Resequence supercell orbitals
HR_TB_bks = HR_TBsup;
res_mat = [1:HR_TB.WAN_NUM/2 ...
           (2*HR_TB.WAN_NUM/2+1):(2*HR_TB.WAN_NUM/2+HR_TB.WAN_NUM/2) ...
           (4*HR_TB.WAN_NUM/2+1):(4*HR_TB.WAN_NUM/2+HR_TB.WAN_NUM/2) ...
           HR_TB.WAN_NUM/2+1:2*HR_TB.WAN_NUM/2 ...
           3*HR_TB.WAN_NUM/2+1:4*HR_TB.WAN_NUM/2 ...
           5*HR_TB.WAN_NUM/2+1:6*HR_TB.WAN_NUM/2];
HR_TBsup = HR_TB_bks.reseq(res_mat);

%% Supercell Wannier band structure
HR_TBsup = HR_TBsup < 'KPOINTS_hex_3D';
Efermi = 0;
EIGENCARs = HR_TBsup.EIGENCAR_gen() - Efermi;
[~, ax] = bandplot(EIGENCARs,[-5,10],'title','wannier-dft','Color','b', ...
    'POSCAR','POSCAR_6','KPOINTS','KPOINTS_hex_3D');

%% Supercell DFT band structure
!cp EIGENVAL-super EIGENVAL;
Efermis = 6.22237991;                           % DFT Fermi energy for supercell
EIGENCAR_DFTs = EIGENVAL_read('vasp','EIGENVAL',Efermis); 
bandplot(EIGENCAR_DFTs,[-5,10], ...
    'title',"Bismuth-super-DFT-vs-Model",'Color','r', ...
    'POSCAR','POSCAR_6','KPOINTS','KPOINTS_hex_3D','ax',ax);

%% Construct slab for surface states
!cp POSCAR_6 POSCAR;
repeatnum   = 6;        % Number of layers in slab
fin_dir     = 3;        % Surface normal direction
glue_edges  = false;    
vacuum_mode = true;     % Add vacuum separation

Kane_Mele_tot_n_slab = HR_TBsup.cut_piece(repeatnum,fin_dir,glue_edges,vacuum_mode);

%% Slab band structure with surface projection
ProjectionStruct.discrimination = 0.0001;
ProjectionStruct.center = [0.5,0.5,0.5];
ProjectionStruct.orientation = fin_dir;
ProjectionStruct.sign = 0;

Efermi = 0;
Kane_Mele_tot_n_slab = Kane_Mele_tot_n_slab < 'KPOINTS_hex_2D';
[EIGENCAR_slab,WAVECAR_disk,WEIGHTCAR_slab] = ...
    Kane_Mele_tot_n_slab.EIGENCAR_gen('LWAVE',true,'WEIGHTCAR',true,'ProjectionMethod','slab');

coolwarm = ColorMap.Matplotlib('coolwarm');
pbandplot(abs(WEIGHTCAR_slab),EIGENCAR_slab, ...
    'POSCAR','POSCAR_6','KPOINTS','KPOINTS_hex_2D', ...
    'cmap',@(varargin) coolwarm(size(coolwarm,1)/2+1:end,:), ... 
    'Ecut',[-1,3],'Linewidth',2);

%% Fermi surface plot
[EIGENCAR_edge,klist1,klist2] = ...
    Kane_Mele_tot_n_slab.EIGENCAR_gen_3D([131,131],[-1.,-1.,-0.;2 0 0;0 2 0],'fermi',Efermi);

Efermi = 0.1400;
[~,Ax]= Figs(1);
bandplot3d(EIGENCAR_edge(:,:,107:110)-Efermi,klist1,klist2, ...
    'Ecut',[0,0.01],'ax',Ax,'xlabel','k_x','ylabel','k_y','title',Efermi);

%% Moiré coupling estimation at Km with changing moiré lattice
HR_TB_bk1 = Kane_Mele_tot_n_slab < 'KPOINTS_Mcoupling';
klist = HR_TB_bk1.klist_frac();
theta = 2*pi/3;
Rz = [cos(theta), -sin(theta), 0;
      sin(theta),  cos(theta), 0;
      0,           0,          1];
klist_rot = (Rz * (klist*HR_TB_bk1.Gk)')'/HR_TB_bk1.Gk;

% Compute slab eigenstates on original and rotated k-paths
[EIGENCAR_slab,WAVECAR_disk]  = Kane_Mele_tot_n_slab.EIGENCAR_gen('LWAVE',true,'WEIGHTCAR',false,'ProjectionMethod','slab','klist',klist);
[EIGENCAR_slab1,WAVECAR_disk1]= Kane_Mele_tot_n_slab.EIGENCAR_gen('LWAVE',true,'WEIGHTCAR',false,'ProjectionMethod','slab','klist',klist_rot);

%% Evaluate the moire coupling factor <C3z>
clear SYMCAR1
kn = length(klist(:,1));
for i = 1:kn
    SYMCAR1(:,:,i) = ( double(WAVECAR_disk(:,107:108,i))' * double(WAVECAR_disk1(:,107:108,i)) ); % select the target surface band
end

% Average coupling strength
for k = 1:length(SYMCAR1)
    % spec_norm(k) = abs(mean(diag(SYMCAR1(:,:,k))));
    spec_norm(k) = norm((SYMCAR1(:,:,k)));
end
figure;
x_vec = linspace(0, 0.33, length(spec_norm));
plot(x_vec,spec_norm') 
ylim([0,1])
xlabel('moireBZ/BZ');
ylabel('Moire coupling factor at Km');