function varargout = bandplot3d(EIGENCAR_3D, klistX, klistY, options)
    % bandplot3d - Visualizes a 3D band structure with optional weight data.
    %
    % Usage:
    %   [ax] = bandplot3d(EIGENCAR_3D, klistX, klistY, options)
    %   [fig, ax] = bandplot3d(EIGENCAR_3D, klistX, klistY, options)
    %
    % Inputs:
    %   EIGENCAR_3D - 3D matrix of energy values for plotting
    %   klistX, klistY - k-point grids for the x and y axes
    %   options - (optional) struct of options with the following fields:
    %     - 'refined' (default: true)
    %     - 'Ecut' (default: [-3, 3])
    %     - 'WEIGHTCAR' (default: [])
    %     - 'ax' (default: empty, creates a new figure if not provided)
    %     - 'fontname' (default: 'Helvetica')
    %     - 'cmap' (default: turbo)
    %     - 'density' (default: 1)
    %     - 'title' (default: '')
    %     - 'xlabel', 'ylabel', 'zlabel' (default: 'k_1', 'k_2', 'E(eV)')
    %     - 'LineSpec', 'LineWidth', 'MarkerSize', 'EdgeColor', 'FaceAlpha'
    %     - 'view' (default: [0, 90])

    % Argument validation and default settings
    arguments
        EIGENCAR_3D double = [];
        klistX double = [];
        klistY double = [];
        options.refined logical = true;
        options.Ecut double = [-3, 3];
        options.WEIGHTCAR double = [];
        options.ax handle = handle([]);
        options.fontname string = 'Helvetica';
        options.KPOINTS string = 'KPOINTS';
        options.POSCAR string = 'POSCAR';
        options.cmap  = turbo;
        options.density double = 1;
        options.title string = '';
        options.klist_l double = [];
        options.kpoints_l double = [];
        options.kpoints_name string = [];
        options.xlabel string = 'k_1';
        options.ylabel string = 'k_2';
        options.zlabel string = 'E(eV)';
        options.LineSpec string = '-';
        options.LineWidth double = 1;
        options.MarkerSize double = 3;
        options.EdgeColor string = 'none';
        options.FaceAlpha double = 0.8;
        options.view double = [0, 90];
    end

    % Handle title
    if isempty(options.title)
        [~, titlestring] = fileparts(pwd);
        titlestring = strrep(titlestring, '_', '-');
    else
        titlestring = options.title;
    end

    % Prepare axes for plotting
    if isempty(options.ax)
        [Fig,ax]= Figs(1, 1); % Create new figure
    else
        if isvalid(options.ax)
            ax = options.ax; % Use provided axis handle
        else
            [Fig,ax]= Figs(1, 1); % Create new figure
        end
    end

    % If WEIGHTCAR is not provided, use absolute value of EIGENCAR_3D
    if isempty(options.WEIGHTCAR)
        WEIGHTCAR = abs(EIGENCAR_3D);
    else
        WEIGHTCAR = options.WEIGHTCAR;
    end

    hold(ax, 'on');
    
    % Refined plotting
    if options.refined
        for i = 1:size(EIGENCAR_3D, 3)
            h = surf(ax, klistX, klistY, EIGENCAR_3D(:,:,i), WEIGHTCAR(:,:,i));
            set(h, 'EdgeColor', options.EdgeColor, 'FaceAlpha', options.FaceAlpha);
        end
    end

    % Set plot limits and labels
    zlim(ax, options.Ecut);
    xlim(ax, [min(min(klistX)), max(max(klistX))]);
    ylim(ax, [min(min(klistY)), max(max(klistY))]);
    xlabel(ax, options.xlabel);
    ylabel(ax, options.ylabel);
    zlabel(ax, options.zlabel);
    colormap(ax, options.cmap);
    view(ax, options.view(1), options.view(2));
    axis equal   

    % Return output handles
    if nargout == 2
        varargout{1} = ax.Parent;
        varargout{2} = ax;
    elseif nargout == 1
        varargout{1} = ax;
    end
end
