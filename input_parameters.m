% This file contains all the variables that should be modified at beginning of every cruise
% This is for onboard processing of AMT 31 in Nov/Dec 2024 on RRS James Cook
%-----------------------------

struct_levels_to_print(0);
warning("off");

# set octave"s graphics_toolkit to the one that works...
graphics_toolkit("gnuplot");

%-----------------------------
#CRUISE = "amt31"; # use for jday 333
CRUISE = "AMT31";  # use for jday 334 onwards
WAP_ROOT = CRUISE;

%-----------------------------
% Variables to be changed during cruise according to specific setups and user needs
%
% Dates
% jday=333 (20241128) is first day of data collection on JC-272 cruise
% Set subdirectories, ini/end dates and WAPhours - code is run by commenting/uncommenting each case

# Step 1: `ACS_both' -
#UWAY_WAP_SUBDIR = "ACS_both/"; (167 in port 3 122 in port 4)
#inidate = "20241128"; # Jday 333 is 28th Nov 2024
#enddate = "20241128";
#WAPhour = "15";

#
#UWAY_WAP_SUBDIR = "ACS_122/"; (122 in port 4)
#inidate = "20241128"; # Jday 334-337
#enddate = "20241202";
#WAPhour = "008";

#UWAY_WAP_SUBDIR = "ACS_167/";
#inidate = "20241202"; # Jday 337 onwards (day we switched out 167 and 122 in port 4)
#enddate = "20241203";
#WAPhour = "03";

#Step 1: `dark_counts_BB3' - run 1
#UWAY_WAP_SUBDIR = "dark_counts_BB3/";
#inidate = "20220602";
#enddate = "20220616";
#WAPhour = "007";
#WAPhour = "010";

UWAY_WAP_SUBDIR = "ACS_167/";
inidate = "20241128"; # Jday 337 onwards (day we switched out 167 and 122 in port 4)
enddate = "20241203";
WAPhour = "03";


% Parameters specific for Underway plotting/processing
% (this will change depending on specific section fo the cruise)
% Setup to automatically change based on UWAY_WAP_SUBDIR

% Implemented instruments to select from are
% {"ctd","acs","bb3","cstar","acs2","ac9","clam"}
if strcmp (UWAY_WAP_SUBDIR, "ACS_both/") == 1 % # acs2 == 167, acs1 ==122
      dh8_instruments = {"bb3", "ctd", "acs2", "acs" }; # -
      dh8_ports = {1,2,3,4};
      dh8_serialnumber = {1173, [], 122};
elseif strcmp(UWAY_WAP_SUBDIR, "ACS_122/") == 1 % tjor: this is the `default" config (i.e. without subdirectories inside WAP_extracted)
        # dh8_instruments = {"bb3", "ctd", "acs", "acs2"};  - full instrument list
      dh8_instruments = {"bb3", "ctd", "acs", }; # - neglecting 167, as it was not used for rest of cruiset
      dh8_ports = {1,2,4};
      dh8_serialnumber = {1173, [], 122};
elseif strcmp(UWAY_WAP_SUBDIR, "ACS_167/") == 1 % tjor: this is the `default" config (i.e. without subdirectories inside WAP_extracted)
        # dh8_instruments = {"bb3", "ctd", "acs", "acs2"};  - full instrument list
      dh8_instruments = {"bb3", "ctd", "acs2", }; # - neglecting 167, as it was not used for rest of cruiset
      dh8_ports = {1,2,4};
      dh8_serialnumber = {1173, [], 122};
elseif strcmp(UWAY_WAP_SUBDIR, "dark_counts_BB3/") == 1 %
     dh8_instruments = {"bb3", "ctd", "acs"};
     dh8_serialnumber = {1173, [], 122};

endif
%-----------------------------

%-----------------------------
% Paths
# MAIN_PATH = "/data/datasets/cruise_data/active/JC272/";
MAIN_PATH = "D:/AMT31/Optics_all/" # local path for laptop on cook
 fflush(stdout);
% MAIN_PATH = [MAIN_PATH, "/Data/", CRUISE,"/"];     % Root directory for current AMT cruise
PATH_DATA = 'Y:/JC272/Data/PML_optics/IOP_Underway/' % Data on ships network onboard James Cook
PATH_SOURCE = [MAIN_PATH, "Source/"];% Directory with all source code
OUT_PROC = [MAIN_PATH, "Processed/"];    % Output directory for processed oct and mat files
OUT_FIGS = [MAIN_PATH, "Figures/"];      % Output directory for figures

addpath([PATH_SOURCE]);

%-----------------------------
% Each directory will contain a series of subdirectories for each instrument
% (e.g. Underway, Optics_rig, BB3_ctd etc. etc.)
UWAY_DIR = '/'
# BB3_DIR = "BB3_ctd/"; -
# CTD_DIR = "Ship_CTD/"; - not present for DY151
% Specific data subdirectories with Underway
DATA_WAPPED = "WAP_Extracted/";
DATA_RAW = "Raw_underway/";
DATA_FLOW = "Flow/";

%-----------------------------
% calibration file dir
D_CAL_FILES = [PATH_DATA , "/WAP_fileconversions/cals/"];
% ACS calibration file
ACS_CAL_FILE_NAME = "acs122.dev"; %
ACS_CAL_FILE_NAME_2 = "acs167.dev";
%-----------------------------

%-----------------------------
% Ship"s system directories - these are specfic to a Discovery cruise metadata format (see AMT28 for JCR).
#PATH_SHIP = [PATH_DATA, "Underway/ship_uway/"]; % note: some cruises store as: [PATH_DATA, "ship_uway/"]
PATH_SHIP = ["Z:/Ship_Systems/Data/TechSAS/NetCDF/"] # netcdf directory files on James cook
PATH_GPS = [PATH_SHIP,'GPS/'];  %
PATH_ATT = [PATH_SHIP,'ATT/'];  %
PATH_DEPTH = [PATH_SHIP,'EA600/'];  %
PATH_TS = [PATH_SHIP,'SURFMETV3/']; %
PATH_TSG = [PATH_SHIP,'TSG/']; %

#----------------------------
# Input parameters for ship"s underway data
# file paths for GGA function
FUNC_GGA = []; # note: handle not used explictly for cookprocessing
# - however, this was done to make similar format in input paramters to the JCR crusies.

DIR_GPS = PATH_GPS;
FN_GPS =  '*position-POSMV_GPS*';

DIR_ATT = PATH_ATT;
FN_ATT = '*shipattitude-POSMV_ATT.att*';

DIR_DEPTH = PATH_DEPTH;
FN_DEPTH = '*EA640_DEPTH*';


# file paths for Oceanlogger function
FUNC_OL = []; # note: handle not used explictly for cook processing

DIR_TS = PATH_TS;
FN_SURF = '*Surf-SURFMET*';
FN_METDATA = '*MET-SURFMET*';
FN_LIGHT = '*Light-SURFMET*';

DIR_TSG = PATH_TSG;
FN_TSG = '*SBE45*';

# Path-related variables for step2
global DIR_STEP1 = [OUT_PROC UWAY_DIR "Step1/"];
global DIR_STEP2 = [OUT_PROC UWAY_DIR "Step2/"];
global DIR_STEP3 = [OUT_PROC UWAY_DIR "Step3/"];
global FN_ROOT_STEP2 = [DIR_STEP2 "proc_optics_" lower(CRUISE) "_"];


% Create path for saving figures
#   global fig_dir = [OUT_FIGS, UWAY_DIR];
global DIR_FIGS = [OUT_FIGS];

%DIR_FIGS
% Create directories if they do not exists
 %%  if ~exist(DIR_FIGS, "dir")
   %   mkdir(DIR_FIGS);
%   endif

 %  if ~exist(DIR_STEP2, "dir")
  %    mkdir(DIR_STEP2);
   % endif



%-----------------------------
% Parameters specific for Optics rig plotting/processing
%
% Whether cdt is saved as ASCII format (false for AMT26; true for AMT27)
ctdASCII = true;
% Limits for temperature and salinity profiles
Tlim = [0 20];
Slim = [33 35];
% Limits for absorption and attenuation profiles
alim = [0.1 0.5];
clim = [0.05 0.6];
chlac9lim = [0 5];
%-----------------------------

% Processors to be used by parcellfun in run_step1par.m
NProc = nproc() - 1;

% Name of processed file to be saved
fproc_name = ["optics_" lower(CRUISE) "_"];

% Limits for time-series plots
acs_raw_lim = [-0.03 0.1]; % acs
flow_lim = [20 45];        % flow rate
bb3_lim = [50 140];       % backscattering
SST_lim = [15 20];         % CTD temperature
SSS_lim = [35 36.5];        % CTD salinity
% Limits for step2 time series
acs_lim = [-0.05 0.3];
ac9_lim = acs_lim;
bb_opt_lim = [70 150];
cstar_lim = [0.75 0.85];
spectra_alim = [0.03];
spectra_clim = [1];
chl_lim = [0.01 5];

%-----------------------------
% Parameters specific for BB3 plotting/processing
%
% Limits for bb3 profiles
bb3lim = [50 300];
%-----------------------------

%-----------------------------
% Parameters specific for underway transect
%
latlim = 54;
trans_SST = [01 30];
trans_SSS = [33 38];
trans_chl = [0.01 5];
trans_cp650 = [0.01 0.3];
%-----------------------------

%-----------------------------
% useful functions
movavg = inline("filter(1/mavgwd*ones(1, mavgwd), 1, x)", "x", "mavgwd"); % this is equivalent to a standard moving average with a window size of mavgwd
%-----------------------------





