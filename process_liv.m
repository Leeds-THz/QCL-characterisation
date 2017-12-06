function process_liv(fileprefix, filesuffix, varargin)
% PROCESS_LIV  Extract data from QCL LIV plots
%   PROCESS_LIV(fileprefix, filesuffix) generates LIV plots
%   PROCESS_LIV(fileprefix, filesuffix, 'ParamName', ParamValue ...)
%
%   The function reads all data files whose name matches the specified
%   format and generates plots of the QCL performance characteristics.
%
%   INPUT FILE FORMAT:
%     The input data files must each contain text in 3 columns:
%     1. The QCL drive current (A)
%     2. The QCL terminal voltage (V)
%     3. The THz detector output signal (assumed mV, but can be any unit)
%
%   OUTPUT PLOTS:
%     Figure 1 - the detector signal as a function of current
%     Figure 2 - the QCL voltage as a function of current
%
%   REQUIRED PARAMETERS:
%     fileprefix - The string that appears at the start of all file names
%     filesuffix - The string that appears at the end of all file names
%
%   OPTIONAL PARAMETERS:
%     'responsivity' - The responsivity of the detector [V/W].  This is
%                      used to scale the power-vs-current plot.  If
%                      unspecified, then a default value of 14400 is used.
%
%     'attenuation'  - The total (linear) optical attenuation factor used
%                      in the experimental measurements.  This is a number
%                      <= 1.  If unspecified, a value of '1' is assumed
%                      (i.e., no attenuation).
%
%     'noisefloor'   - The noise floor (mV) for the measurement.  This is
%                      used to determine the threshold current.  If not
%                      specified, a value of 1 mV is assumed.
%
%    EXAMPLES:
%      process_liv('LIV-', 'K.dat');
%        Processes all files in the current folder, with names in the
%        format 'LIV-*K.dat', where the '*' is taken as the heat-sink
%        temperature.
%
%      process_liv('L1071-', 'Kelvin.txt', 'responsivity', '100');
%        Processes all files in the current folder, with names in the
%        format 'L1071-*Kelvin.txt', where the '*' is taken as the
%        heat-sink temperature.  The responsivity is taken as 100 W/V.
%     
% (c) Alexander Valavanis
%     University of Leeds, 2017

%% Handle input arguments
parser = inputParser;

% Define REQUIRED function arguments
addRequired( parser, 'fileprefix', @ischar);
addRequired( parser, 'filesuffix', @ischar);

% Define OPTIONAL function arguments 
default_responsivity = 14400; % [V/W] (Default value for B2)
default_attenuation  = 1;     % Linear attenuation factor (a fraction < 1)
default_noise_floor  = 1;     % [mV]

addParameter(parser, 'responsivity', default_responsivity, @isnumeric);
addParameter(parser, 'attenuation',  default_attenuation,  @isnumeric);
addParameter(parser, 'noisefloor',   default_noise_floor,  @isnumeric);

% Parse all function arguments
parse(parser, fileprefix, filesuffix, varargin{:});
responsivity = parser.Results.responsivity;
attenuation  = parser.Results.attenuation;
noise_floor  = parser.Results.noisefloor;

%% Find data files
% Find a list of all files matching the filename pattern
nameformat = [fileprefix '*' filesuffix];
allfiles = ls(nameformat);
nfiles = size(allfiles,1);

if(nfiles <= 1)
    error('No input files found with name format: %s', nameformat);
end

%% Generate figures and set axis labels
fig_iv = figure('Name', 'I-V plot');
hold on;
ax_iv = gca;
ax_iv.XLabel.String = 'Current (A)';
ax_iv.YLabel.String = 'Voltage (V)';

fig_li = figure('Name', 'L-I plot');
hold on;
ax_li = gca;
ax_li.XLabel.String = 'Current (A)';
ax_li.YLabel.String = 'THz power (mW)';

fig_eff = figure('Name', 'Efficiency plot');
hold on;
ax_eff = gca;
ax_eff.XLabel.String = 'Current (A)';
ax_eff.YLabel.String = 'Efficiency (\%)';

fig_T_Ppk = figure('Name', 'Power-temperature plot');
hold on;
ax_T_Ppk = gca;
ax_T_Ppk.XLabel.String = 'Heat-sink temperature (K)';
ax_T_Ppk.YLabel.String = 'Peak THz power (mW)';

fig_T_Ith = figure('Name', 'Threshold current vs temperature');
hold on;
ax_T_Ith = gca;
ax_T_Ith.XLabel.String = 'Heat-sink temperature (K)';
ax_T_Ith.YLabel.String = 'Threshold current (A)';

%% Create storage for output variables
temperature       = zeros(1,nfiles); % [K]
peak_THz_power    = zeros(1,nfiles); % [mW]
threshold_current = zeros(1,nfiles); % [A]

%% Current range
I_range = [1000 -1000];

%% Loop through files and add to plots
for ifile = 1:nfiles
    %% Read data from input file
    filename = allfiles(ifile,:);
    filedata = load(filename);

    current = filedata(:,1); % [A]
    voltage = filedata(:,2); % [V]
    det_v   = filedata(:,3); % [mV]

    %% Update current range
    if(max(current) > I_range(2))
        I_range(2) = max(current);
    end
    
    if(min(current) < I_range(1))
        I_range(1) = min(current);
    end
    
    %% Find the temperature from the filename
    % Strip prefix and suffix off the filename
    endpart   = strsplit(filename, fileprefix);
    frontpart = strsplit(endpart{2}, filesuffix);
    
    temperature(ifile) = str2double(frontpart{1});

    %% Calculate input and output powers
    elec_power = current .* voltage; % [W]
    THz_power  = det_v / (responsivity * attenuation); % [mW]
    
    %% Find threshold current
    % Find first point at which the detector signal exceeds the noise
    % floor and define this as the threshold current
    threshold_current(ifile) = current(find(det_v > noise_floor, 1));
    efficiency = THz_power ./ (elec_power*1000);
    
    %% Find peak output power
    peak_THz_power(ifile) = max(THz_power);
    
    %% Update the LI, IV and efficiency plots
    plot(ax_iv,  current, voltage);
    plot(ax_li,  current, THz_power);
    plot(ax_eff, current, efficiency * 100);
end

%% Update current ranges on relevant axes
ax_iv.XLim  = I_range;
ax_li.XLim  = I_range;
ax_eff.XLim = I_range;

%% Plot the peak power as a function of heat-sink temperature
plot(ax_T_Ppk, temperature, peak_THz_power, 'ko');

%% Fit an exponential function to the threshold current - heat-sink data
fit_model  = 'I0 + I1*exp(x/T0)';
startpoint = [0 threshold_current(1) 10];
T_I_fit = fit(temperature', threshold_current',fit_model,...
              'StartPoint', startpoint);
I0_fitted = T_I_fit(temperature');
plot(ax_T_Ith, temperature', I0_fitted, 'k-');
plot(ax_T_Ith, temperature', threshold_current', 'ko');

%% Print plots to file, in various formats
print(fig_iv,    'I-V',        '-dpdf', '-r600');
print(fig_li,    'L-I',        '-dpdf', '-r600');
print(fig_eff,   'efficiency', '-dpdf', '-r600');
print(fig_T_Ppk, 'Ppk-T',      '-dpdf', '-r600');
print(fig_T_Ith, 'Ith-T',      '-dpdf', '-r600');

print(fig_iv,    'I-V',        '-dpng', '-r600');
print(fig_li,    'L-I',        '-dpng', '-r600');
print(fig_eff,   'efficiency', '-dpng', '-r600');
print(fig_T_Ppk, 'Ppk-T',      '-dpng', '-r600');
print(fig_T_Ith, 'Ith-T',      '-dpng', '-r600');

print(fig_iv,    'I-V',        '-deps', '-r600');
print(fig_li,    'L-I',        '-deps', '-r600');
print(fig_eff,   'efficiency', '-deps', '-r600');
print(fig_T_Ppk, 'Ppk-T',      '-deps', '-r600');
print(fig_T_Ith, 'Ith-T',      '-deps', '-r600');