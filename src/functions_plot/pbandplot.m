function varargout = pbandplot(WEIGHTCAR_struct, EIGENCAR, klist_l, kpoints_l, kpoints_name, options, optionsplot)
%PBANDPLOT Projected Band Structure Plotting Tool
%   PBANDPLOT visualizes the projected band structure using provided weight data and eigenvalue information.
%   This function supports multiple plotting modes including patch, bubble, and refined bubble plots.
%
%   Syntax:
%       varargout = pbandplot(WEIGHTCAR_struct, EIGENCAR, klist_l, kpoints_l, kpoints_name, options, optionsplot)
%
%   Parameters:
%       WEIGHTCAR_struct - Projection weight data (struct/cell/array)
%       EIGENCAR         - Eigenvalue data (default: EIGENVAL_read())
%       klist_l          - k-path coordinates
%       kpoints_l        - Special k-points positions
%       kpoints_name     - Labels for special k-points
%       options          - Plot configuration options (see below)
%       optionsplot      - Advanced plotting parameters
%
%   Options Structure Fields:
%       ColorCut/ColorCutMinus - Scaling factors for positive/negative color mapping
%       Ecut            - Energy range [min, max]
%       fontname        - Font family for labels
%       KPOINTS/POSCAR  - Input files for k-path generation
%       cmap            - Colormap function/handle
%       title           - Plot title
%       LineSpec        - Line style specification
%       silent          - Suppress projection selection prompt
%
%   Outputs:
%       Returns figure and axis handles when requested
%
%   Example:
%       ax = pbandplot(WEIGHTCAR_data, EIGENCAR, 'cmap', @jet, 'Ecut', [-2 2]);
%
%   See also: EIGENVAL_read, POSCAR_read, kpathgen3D

%-------- Initialization & Argument Processing --------
arguments
    WEIGHTCAR_struct = [];
    EIGENCAR = EIGENVAL_read();
    klist_l double = [];
    kpoints_l double = [];
    kpoints_name string = [];
    options.ColorCut double = 1;
    options.ColorCutMinus double = -1;
    options.Ecut = [-3,3];
    options.fontname = 'Helvetica';
    options.KPOINTS = 'KPOINTS';
    options.POSCAR = 'POSCAR';
    options.cmap = @hsv;
    options.title = '';
    options.klist_l = [];
    options.kpoints_l = [];
    options.kpoints_name = [];
    options.xlabel = '';
    options.ylabel = 'E(eV)';
    options.LineSpec = '-';
    options.LineWidth = 1;
    options.MarkerEdgeColor = [];
    options.MarkerFaceColor = [];
    options.silent = false;
    options.plotMode = ''
    options.legend = []
    optionsplot.density = 1;
    optionsplot.WEIGHTCAR_factor = 1;
    optionsplot.ax = handle([]);
    optionsplot.filled = true;   
end

%-------- Kpath Validation --------
if isempty(klist_l)
    Rm = POSCAR_read(options.POSCAR);
    [~,klist_l,~,kpoints_l,kpoints_name] = kpathgen3D(Rm, options.KPOINTS);
end
[klist, kpoints_l, kpoints_name] = validateKPathInputs(options, klist_l, kpoints_l, kpoints_name);

%-------- Title Processing --------
if strcmp(options.title,'')
    [~,titlestring] = fileparts(pwd);
    titlestring = strrep(titlestring,'_','-');
else
    titlestring = options.title;
end

%-------- Color Scaling Setup --------
if options.ColorCutMinus == -1
    options.ColorCutMinus = options.ColorCut;
end


%-------- Multi-Projection Handling --------
if isempty(WEIGHTCAR_struct)
    [varargout{1:nargout}]  = handleMultiProjections(EIGENCAR,klist_l, kpoints_l,kpoints_name,klist,options, titlestring,optionsplot);
    return;
end

%-------- Plot Mode Determination --------
if strcmp(options.plotMode, 'RGB')
    if ~iscell(WEIGHTCAR_struct) || numel(WEIGHTCAR_struct)~=3
        error("For RGB mode, WEIGHTCAR must be a cell with 3 kinds of projs.")
    end
    plotMode = options.plotMode;
elseif strcmp(options.plotMode, 'RWB')
    if ~iscell(WEIGHTCAR_struct) || numel(WEIGHTCAR_struct)>2
        error("For RWB mode, WEIGHTCAR must be a cell with 1 or 2 kinds of projs.")
    end
    plotMode = options.plotMode;
else
    [plotMode, cmap, nProj, projNames, WEIGHTCAR] = resolvePlotMode(WEIGHTCAR_struct, options);
    %-------- Projection Selection --------
    selectedProjs = handleProjectionSelection(plotMode, nProj, options);
end

%-------- Axis Preparation --------
ax = prepareAxis(optionsplot);
Nbands=size(EIGENCAR,1);
xmin=klist_l(1);
xmax=klist_l(length(klist_l));
Xcut = [xmin,xmax];
%-------- Core Plotting Logic --------
switch plotMode
    case 'patch' % line
        ax = pband_plot_one_patch(klist,EIGENCAR,WEIGHTCAR,'cmap',cmap,'ax',ax,'LineWidth',options.LineWidth);
    case 'bubble_only'
        optionsplot.ax =ax;
        optionsplotcell = namedargs2cell(optionsplot);
        ax = pband_plot_one_bubble(klist,EIGENCAR,WEIGHTCAR,cmap,projNames,optionsplotcell{:});
        legend(ax);
    case {'bubble_rough','bubble_refined'}
        for Ei=1:Nbands
            plot(ax,klist,EIGENCAR(Ei,:),...
                'LineWidth',options.LineWidth,'Color',[0.1 0.1 0.1],'DisplayName',num2str(Ei),...
                'HandleVisibility','off');
        end
        optionsplot.ax =ax;
        optionsplotcell = namedargs2cell(optionsplot);
        ax = pband_plot_set(klist,EIGENCAR,WEIGHTCAR,projNames,selectedProjs,...
            'cmap',cmap,optionsplotcell{:});% waiting
        legend(ax);
    case 'RGB'
        X = klist;
        
        WEIGHT_R = WEIGHTCAR_struct{1};
        WEIGHT_G = WEIGHTCAR_struct{2};
        WEIGHT_B = WEIGHTCAR_struct{3};
        
        for i = 1:size(EIGENCAR,1)
            Y = [EIGENCAR(i,:); EIGENCAR(i,:)];
            Z = zeros(size(Y));
        
            WEIGHT_RGB_i = [WEIGHT_R(i,:)',WEIGHT_G(i,:)',WEIGHT_B(i,:)'];
            for j = 1:size(WEIGHT_RGB_i,1)
                WEIGHT_RGB_i(j,:) = WEIGHT_RGB_i(j,:)/norm(WEIGHT_RGB_i(j,:));
            end
            WEIGHT_RGB_i = reshape(WEIGHT_RGB_i, [1 size(WEIGHT_RGB_i)]);
            WEIGHT_RGB_i = repmat(WEIGHT_RGB_i,2,1,1);
            surf(ax, X,Y,Z, WEIGHT_RGB_i, ...
                'EdgeColor', 'interp', 'FaceColor', 'none', 'LineWidth', options.LineWidth, ...
                'HandleVisibility','off')
        end
        
        if ~isempty(options.legend)
            plot(ax, -[xmin xmax], [0 0], 'Color','r', 'LineWidth',2)
            plot(ax, -[xmin xmax], [0 0], 'Color','g', 'LineWidth',2)
            plot(ax, -[xmin xmax], [0 0], 'Color','b', 'LineWidth',2)
            legend(ax, options.legend)
        end
    case 'RWB'
        % red-white-blue
        X = klist;
        switch numel(WEIGHTCAR_struct)
            case 2 
                % red and blue
                WEIGHT = WEIGHTCAR_struct{1} - WEIGHTCAR_struct{2};
            case 1          
                % positive value to red, negative value to blue
                WEIGHT = WEIGHTCAR_struct{1};
        end
        
        colormap(ax, ColorMap.redblue(length(WEIGHT)));
        cmax = max(abs(WEIGHT),[],"all");
        clim(ax, [-cmax, cmax])
        
        for i = 1:size(EIGENCAR,1)
            Y = [EIGENCAR(i,:); EIGENCAR(i,:)];
            Z = zeros(size(Y)); 
        
            WEIGHT_i = WEIGHT(i,:);
            WEIGHT_i = repmat(WEIGHT_i,2,1);
            surf(ax, X,Y,Z,WEIGHT_i, ...
                'EdgeColor', 'interp', 'FaceColor', 'none', 'LineWidth', options.LineWidth, ...
                'HandleVisibility','off')
        end
        if ~isempty(options.legend)
            plot(ax, -[xmin xmax], [0 0], 'Color','r', 'LineWidth',2)
            plot(ax, -[xmin xmax], [0 0], 'Color','b', 'LineWidth',2)
            legend(ax, options.legend)
        end
end


%-------- Reference Lines & Labels --------

ax = set_reference(kpoints_l,kpoints_name,Xcut,options.Ecut,...
    'ax',ax,...
    'fontname',options.fontname ...
    );

%--------  fbug  --------
if iscell(titlestring )
    titlestring = cell2mat(titlestring);
end

title(ax,titlestring);
%-------- Output Handling --------
if nargout > 2
    plot_data_coll.EIGENCAR = EIGENCAR;
    plot_data_coll.Ecut = Ecut;
    plot_data_coll.titlestring = titlestring;
    plot_data_coll.color = color;
    plot_data_coll.klist_l = klist_l;
    plot_data_coll.kpoints_l = kpoints_l;
    plot_data_coll.kpoints_name = kpoints_name;
    plot_data_coll.fontname = options.fontname;
    plot_data_coll.WEIGHTCAR_struct = WEIGHTCAR_struct;
    varargout{3} = plot_data_coll;
end
[varargout{1:2}] = deal(ax.Parent, ax);
end

%% Helper Functions
function [plotMode, cmap, nProj, projNames, WEIGHTCAR] = resolvePlotMode(WEIGHTCAR_struct, options)
% RESOLVEPLOTMODE Determine plotting mode and process weights
% Process colormap

if isstruct(WEIGHTCAR_struct)
    % Structure array handling (refined bubble mode)
    plotMode = 'bubble_refined';
    nProj = length(WEIGHTCAR_struct);
    projNames = arrayfun(@(x) x.displayname, WEIGHTCAR_struct, 'UniformOutput', false);
    WEIGHTCAR = {WEIGHTCAR_struct.WEIGHTCAR};
elseif iscell(WEIGHTCAR_struct)
    % Cell array handling (basic bubble mode)
    nProj = length(WEIGHTCAR_struct);
    WEIGHTCAR = WEIGHTCAR_struct;
    if length(WEIGHTCAR_struct) > 1
        plotMode = 'bubble_rough' ;
        projNames(nProj ,:) = "";
    else
        plotMode = 'bubble_only' ;
        projNames(nProj ,:) = "";
    end
    for i = 1:nProj
        projNames(i,:) = strcat('A',string(i));
    end
elseif isnumeric(WEIGHTCAR_struct)
    % Matrix input handling (patch mode)
    plotMode = 'patch';
    nProj = 1;
    projNames = {'Total Projection'};
    WEIGHTCAR = WEIGHTCAR_struct;
    maxBCplus = options.ColorCut*max((WEIGHTCAR),[],'all');
    maxBCMinus = options.ColorCutMinus*min((WEIGHTCAR),[],'all');
    maxBC = max(abs(maxBCplus),abs(maxBCMinus));
    WEIGHTCAR(WEIGHTCAR>maxBCplus) = maxBCplus;
    WEIGHTCAR(WEIGHTCAR<maxBCMinus) = maxBCMinus;
    if isa(options.cmap,'function_handle')
        try
            cmap = options.cmap(64,[-maxBC,maxBC]);
        catch
            cmap = options.cmap(64);
        end
    end
    return;
else
    error('Invalid WEIGHTCAR_struct format');
end

if isa(options.cmap, 'function_handle')
    cmap = options.cmap(nProj);
else
    cmap = options.cmap;
end

end

function selectedProjs = handleProjectionSelection(plotMode, nProj, options)
if nProj ~= 9
    % warning("It seems like that the PBAND info does not come from vaspkit!")
end
if ~options.silent && (strcmp(plotMode,"bubble_only") || strcmp(plotMode,"bubble_rough") ||strcmp(plotMode,"bubble_refined"))
    prompt = "0) all(or just press Enter)\n1) s    2) py    3) pz    4) px\n"...
        +"5) dxy  6) dyz   7) dz2   8) dxz   9) dx2-y2:\n";
    selectedProjs = input(prompt);
    if isempty(selectedProjs) 
        selectedProjs = 1:nProj;
    elseif isscalar(selectedProjs)
        if selectedProjs==0 
            selectedProjs = 1:nProj;
        end
    end
else
    selectedProjs = 1:nProj;
end
end

function [klist, kpoints_l, kpoints_name] = validateKPathInputs(options, klist_l, kpoints_l, kpoints_name)
% Validate and override k-path parameters
klist = options.klist_l;
if isempty(klist), klist = klist_l; end
if ~isempty(options.kpoints_l), kpoints_l = options.kpoints_l; end
if ~isempty(options.kpoints_name), kpoints_name = options.kpoints_name; end
end

function varargout = handleMultiProjections(EIGENCAR, klist_l, kpoints_l,kpoints_name,klist,options, titlestring,optionsplot)
    % HANDLEMULTIPROJECTIONS Handles multiple projections for plotting.
    % Generate WEIGHTCAR structures and name list
    [WEIGHTCAR_struct_cell, Name_list, ~] = WEIGHTCAR_read_dir('PBAND');
    [~, nWEIGHTCAR_struct_cell] = size(WEIGHTCAR_struct_cell);
    
    % Loop through each WEIGHTCAR structure for processing
    WEIGHTCAR_BK = {};
    projNames_BK = {};
    for i = 1:nWEIGHTCAR_struct_cell
        % Split the name to extract the element name
        element_name = split(Name_list{i},'_');
        titlestring_new = titlestring + "-" + element_name(1);
        
        % Prompt user for projection selection if not in silent mode
        if ~options.silent
            disp("Please choose which projections of " + element_name(1) + " will be plotted:");
            disp("(array format like 1:4, [1,4,5]... is supported)");
        end

        % Set the title for the current plot
        options.title = titlestring_new;
        % propertyCell = namedargs2cell(options); % Convert named arguments to cell array

        WEIGHTCAR_struct = WEIGHTCAR_struct_cell{i}; % Get the current WEIGHTCAR structure
        WEIGHTCAR_struct(end) = []; % Remove the last entry from the structure
        
        % Generate the plot based on the current WEIGHTCAR structure
        %ax_temp = pbandplot(WEIGHTCAR_struct, EIGENCAR, klist_l, kpoints_l, kpoints_name, propertyCell{:});
        [plotMode, cmap, nProj, projNames, WEIGHTCAR] = resolvePlotMode(WEIGHTCAR_struct, options);
        selectedProjs = handleProjectionSelection(plotMode, nProj, options);
        %-------- Axis Preparation --------
        [~,ax_temp] =  Figs(1,1);
        Nbands=size(EIGENCAR,1);
        xmin=klist_l(1);
        xmax=klist_l(length(klist_l));
        Xcut = [xmin,xmax];
        for Ei=1:Nbands
            plot(ax_temp,klist,EIGENCAR(Ei,:),...
                'LineWidth',1.0,'Color',[0.1 0.1 0.1],'DisplayName',num2str(Ei),...
                'HandleVisibility','off');
        end
        optionsplot.ax =ax_temp;
        optionsplotcell = namedargs2cell(optionsplot);
        ax_temp = pband_plot_set(klist,EIGENCAR,WEIGHTCAR,projNames,selectedProjs,...
            'cmap',cmap,optionsplotcell{:});% waiting
        %-------- Reference Lines & Labels --------
        ax_temp = set_reference(kpoints_l,kpoints_name,Xcut,options.Ecut,...
            'ax',ax_temp,...
            'fontname',options.fontname ...
            );
        legend(ax_temp);
        %title(ax_temp,titlestring_new);
        % Store the figure and axis handles
        fig(i) = ax_temp.Parent; % Get the parent figure of the axis
        ax(i) = ax_temp; % Store the axis handle
        WEIGHTCAR_BK = [WEIGHTCAR_BK,WEIGHTCAR(selectedProjs)];
        projNames_BK = [projNames_BK,projNames(selectedProjs)];
    end
    if nWEIGHTCAR_struct_cell>1
        nProj = length(WEIGHTCAR_BK);
        if isa(options.cmap, 'function_handle')
            cmap = options.cmap(nProj);
        else
            cmap = options.cmap;
        end
        [~,ax(i+1)] =  Figs(1,1);
        Nbands=size(EIGENCAR,1);
        xmin=klist_l(1);
        xmax=klist_l(length(klist_l));
        Xcut = [xmin,xmax];
        for Ei=1:Nbands
            plot(ax(i+1),klist,EIGENCAR(Ei,:),...
                'LineWidth',1.0,'Color',[0.1 0.1 0.1],'DisplayName',num2str(Ei),...
                'HandleVisibility','off');
        end
        optionsplot.ax = ax(i+1);
        optionsplotcell = namedargs2cell(optionsplot);
        ax(i+1) = pband_plot_set(klist,EIGENCAR,WEIGHTCAR_BK,projNames_BK,1:nProj,...
            'cmap',cmap,optionsplotcell{:});% waiting
        %-------- Reference Lines & Labels --------
        ax(i+1) = set_reference(kpoints_l,kpoints_name,Xcut,options.Ecut,...
            'ax',ax(i+1),...
            'fontname',options.fontname ...
            );
        legend(ax(i+1));
        %--------  fbug  --------
        if iscell(titlestring )
            titlestring = cell2mat(titlestring);
        end
        title(ax(i+1),titlestring);
    end
    % Return handles based on the number of requested output arguments
    if nargout == 2
        varargout{1} = fig; % Return figure handles
        varargout{2} = ax;  % Return axis handles
    elseif nargout == 1
        varargout{1} = ax;  % Return only axis handles
    end
end

function ax = prepareAxis(optionsplot)
if isempty(optionsplot.ax)
    [~,ax] =  Figs(1,1);
else
    ax = optionsplot.ax;
end
end

