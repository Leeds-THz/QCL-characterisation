function process_spectra(fileprefix, filesuffix, varargin)

%% Handle input arguments
parser = inputParser;

% Define REQUIRED function arguments
addRequired(parser, 'fileprefix', @ischar);
addRequired(parser, 'filesuffix', @ischar);

% Define OPTIONAL function arguments 
default_controlvartype = 'current'; % Either 'current' or 'temperature'
default_fmin = 1; % Minimum frequency to plot [THz]
default_fmax = 5; % Maximum frequency to plot [THz]

addParameter(parser, 'controlvartype', default_controlvartype, @isstring);
addParameter(parser, 'fmin', default_fmin, @isnumeric);
addParameter(parser, 'fmax', default_fmax, @isnumeric);

% Parse all function arguments
parse(parser, fileprefix, filesuffix, varargin{:});
controlvartype = parser.Results.controlvartype;
fmin = parser.Results.fmin;
fmax = parser.Results.fmax;

%% Find data files
% Find a list of all files matching the filename pattern
nameformat = [fileprefix '*' filesuffix];
allfiles = ls(nameformat);
nfiles = size(allfiles,1);

if(nfiles <= 1)
    error('No input files found with name format: %s', nameformat);
end

%% Generate figures and set axis labels
fig_spectra = figure('Name', 'Spectral plot');
hold on;
ax_spectra = gca;
ax_spectra.XLabel.String = 'Frequency (THz)';
ax_spectra.YLabel.String = 'Spectral intensity (a.u.)';

%% Create storage for output variables

% The QCL control variable (either temperature [K] or current [mA])
controlvar = zeros(1,nfiles);

%% Loop through files and add to plots
for ifile = 1:nfiles
    %% Read data from input file
    filename = allfiles(ifile,:);
    [frequency, intensity] = read_dpt(filename,     ...
                                      'fmin', fmin, ...
                                      'fmax', fmax);
    
    %% Find the current or temperature from the filename
    % Strip prefix and suffix off the filename
    endpart   = strsplit(filename, fileprefix);
    frontpart = strsplit(endpart{2}, filesuffix);

    controlvar(ifile) = str2double(frontpart{1});

    yoffset = ifile - 1; % Vertical offset between plots

    % Add plot to the vertical scale
    plot(ax_spectra, frequency, 0.9 * intensity + yoffset);

    % Add label to plot
    unit = ' mA';
    if (strcmp(controlvartype, 'temperature'))
        unit = ' K';
    end

    text(fmin + 0.05 * (fmax - fmin), yoffset + 0.3, ...
        [num2str(controlvar(ifile)) unit]);
end

ax_spectra.XLim = [fmin fmax];
ax_spectra.YLim = [0 nfiles];

%% Print plots to file, in various formats
print(fig_spectra, 'spectra', '-dpdf', '-r600');
print(fig_spectra, 'spectra', '-deps', '-r600');
print(fig_spectra, 'spectra', '-dpng', '-r600');