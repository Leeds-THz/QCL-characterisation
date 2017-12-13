function [frequency, intensity] = read_dpt(filename, varargin)
% READ_DPT  Read THz spectral data from a Bruker DPT file
%
%    [frequency, intensity] = READ_DPT(filename)
%
%    The function reads the raw data (intensity vs wavenumber) from file
%    and returns the spectral data on a THz scale.
%
%    INPUT FILE FORMAT:
%       Input files should be Bruker DPT files containing 2 columns of
%       text:
%       1. Wavenumber (1/cm) in numerically decreasing order
%       2. Spectral intensity (arb. units)
%
%    PARAMETERS:
%       filename  - The name of the DPT file to read
%
%    RETURN VALUES:
%       frequency - The frequency in THz (increasing order)
%       intensity - The spectral intensity (a.u.) at each frequency

%% Handle input arguments
parser = inputParser;

% Define REQUIRED function arguments
addRequired(parser, 'filename', @ischar);

% Define OPTIONAL function arguments 
default_fmin = 1; % Minimum frequency to plot [THz]
default_fmax = 5; % Maximum frequency to plot [THz]

addParameter(parser, 'fmin', default_fmin, @isnumeric);
addParameter(parser, 'fmax', default_fmax, @isnumeric);

parse(parser, filename, varargin{:});
fmin = parser.Results.fmin;
fmax = parser.Results.fmax;

%% Read raw data from file
filedata = load(filename);

wavenumber = filedata(:,1); % [1/cm]
intensity  = filedata(:,2); % [a.u.]

% Convert wavenumber to frequency
frequency = wavenumber * physconst('LightSpeed') / 1e10; % [THz]

% Reverse data direction as Bruker system saves from high-to-low freq
frequency = fliplr(frequency);
intensity = fliplr(intensity);

% Find range of frequencies within the specified range
indices = frequency > fmin & frequency < fmax;

% Truncate frequency and intensity arrays to desired range
frequency = frequency(indices);
intensity = intensity(indices);

% Normalise intensity to 1
intensity = intensity / max(intensity);