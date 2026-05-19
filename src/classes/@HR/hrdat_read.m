function [dataArray, NRPT_list, NRPTS, NUM_WAN] = hrdat_read(filename)
%HRDAT_READ Reads Wannier90 Hamiltonian data from _hr.dat file

    if nargin < 1 || isempty(filename)
        filename = 'wannier90_hr.dat';
    end

    %% Read header
    fileID = fopen(filename, 'r');
    if fileID == -1
        error('Failed to open file: %s', filename);
    end

    headerData = textscan(fileID, '%d%*s', 2, ...
        'HeaderLines', 1, ...
        'Delimiter', ' ', ...
        'MultipleDelimsAsOne', true);

    fclose(fileID);

    if isempty(headerData{1}) || numel(headerData{1}) < 2
        error('Failed to read NUM_WAN and NRPTS from file: %s', filename);
    end

    NUM_WAN = headerData{1}(1);
    NRPTS   = double(headerData{1}(2));

    NRPT_lines = ceil(NRPTS / 15);

    %% Read NRPT_list
    fileID = fopen(filename, 'r');
    if fileID == -1
        error('Failed to open file: %s', filename);
    end

    fgetl(fileID);
    fgetl(fileID);
    fgetl(fileID);

    NRPT_list = zeros(NRPTS, 1);
    count = 0;

    for iline = 1:NRPT_lines
        thisLine = fgetl(fileID);

        if ~ischar(thisLine)
            fclose(fileID);
            error('Unexpected end of file while reading NRPT_list.');
        end

        vals = sscanf(thisLine, '%f');
        nval = numel(vals);

        if count + nval > NRPTS
            nval = NRPTS - count;
        end

        NRPT_list(count+1:count+nval) = vals(1:nval);
        count = count + nval;
    end

    fclose(fileID);

    if count ~= NRPTS
        error('Failed to read NRPT_list: expected %d entries, got %d.', ...
            NRPTS, count);
    end

    %% Read hopping parameters
    fileID = fopen(filename, 'r');
    if fileID == -1
        error('Failed to open file: %s', filename);
    end

    startRow = 4 + NRPT_lines;

    hoppingFormat = '%f %f %f %f %f %f %f';

    dataArray = textscan(fileID, hoppingFormat, ...
        'HeaderLines', startRow - 1, ...
        'Delimiter', ' ', ...
        'MultipleDelimsAsOne', true);

    fclose(fileID);

end