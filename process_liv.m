function process_liv(fileprefix, filesuffix, varargin)
% PROCESS_LIV  Extract data from QCL LIV plots
%   PROCESS_LIV(fileprefix, filesuffix) generates LIV plots
%   PROCESS_LIV(fileprefix, filesuffix, 'ParamName', ParamValue ...)
%
%   The function reads all data files whose name matches a specified format
%   and plots:
%     Figure 1 - the detector signal as a function of current
%     Figure 2 - the QCL voltage as a function of current
%
%   Required parameters:
%     fileprefix - The string that appears at the start of all file names
%     filesuffix - The string that appears at the end of all file names
%
%   Optional parameters:
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
%                      specified, a value of 1 uV is assumed.
%
% (c) Alexander Valavanis
%     University of Leeds, 2017

parser = inputParser;

% Define REQUIRED function arguments
addRequired( parser, 'fileprefix', @ischar);
addRequired( parser, 'filesuffix', @ischar);

% Define OPTIONAL function arguments 
default_responsivity = 14400; % [V/W] (Default value for B2)
default_attenuation  = 1;     % Linear attenuation factor (a fraction < 1)
default_noise_floor  = 1e-3;  % [mW] (Default value = 10 uW)

addParameter(parser, 'responsivity', default_responsivity, @isnumeric);
addParameter(parser, 'attenuation',  default_attenuation,  @isnumeric);
addParameter(parser, 'noisefloor',   default_noise_floor,  @isnumeric);

% Parse all function arguments
parse(parser, fileprefix, filesuffix, varargin{:});
responsivity = parser.Results.responsivity;
attenuation  = parser.Results.attenuation;
noise_floor  = parser.Results.noisefloor;

%% Generate figures
% Plot of all I-V curves
figure;
iv_fig = gcf;
hold on;

% Plot of all L-I curves
figure
li_fig = gcf;
hold on;

% Plot of efficiency
figure
eff_fig = gcf;
hold on;

%% Find data files
% Find a list of all files matching the filename pattern
allfiles = ls([fileprefix '*' filesuffix]);
nfiles = size(allfiles,1);

%% Create storage for output variables
temperature    = zeros(1,nfiles); % [K]
peak_THz_power = zeros(1,nfiles); % [mW]

%% Loop through files and add to plots
for ifile = 1:nfiles
    filename = allfiles(ifile,:);

    % Strip prefix off the filename
    endpart  = strsplit(filename, fileprefix);
    
    % Strip suffix off the filename
    frontpart = strsplit(endpart{2}, filesuffix);

    % Now get the temperature from the file name
    temperature(ifile) = str2double(frontpart{1});

    filedata = load(filename);
   
    current = filedata(:,1); % [A]
    voltage = filedata(:,2); % [V]
    det_v   = filedata(:,3); % [mV]

    elec_power = current .* voltage; % [W]
    THz_power  = det_v / (responsivity * attenuation); % [mW]

    threshold_current(ifile) = current(find(THz_power > noise_floor, 1));
    efficiency = THz_power ./ (elec_power*1000);
    
    peak_THz_power(ifile) = max(THz_power);
    
    figure(iv_fig);
    plot(current, voltage);
      
    figure(li_fig);
    plot(current, THz_power);

    figure(eff_fig);
    plot(current, efficiency * 100);
end

figure;
plot(temperature, peak_THz_power);

figure;
T_I_fit = fit(temperature', threshold_current', 'exp1');
plot(T_I_fit, temperature', threshold_current');

