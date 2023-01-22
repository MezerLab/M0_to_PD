function [logname]=M0_ToPD(outDir, Qmap_file, M0_file, BM_file, seg_file,Qmap_factor)
%
% This function performs the M0-PD fit. The local fits are joined together into
% one PD image.
%
% Intputs: 
%       - outDir - path of the directory in which the final PD image will
%                   be formed in.
%       - Qmap_file - path of the directory of a map to which the PD image
%                   will be fitted to.
%       - M0_file - path of the directory of the M0 file to transform to
%                   PD (in weighted images, PDw image can be used here).
%       - seg_file - path to a segmentation file that contains typically 3
%                   different segment types (CSF, GM, WM).
%       - Qmap_factor - a scalar, to be determined with respect to the
%                   ratio between the order of the Qmap_file and the 
%                   M0_file values. For T1 and M0, Qmap_factor=1.
%
% Outputs:
%       - logname - path of the directory the opt variable is saved at (opt
%                   is a structure moving between functions, containing the
%                   parameters for the fit).
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel
%   2022
%
%

%%

degrees = 2;                    % the polynom degree used to fit the gain
outMm = [1 1 1];                % the resolution (in mm) on which to perform the fit
boxSize = 14;                   % the box size on which to perform the fit
Inclusion_Criteria = [0.7 100]; % which boxes to exclude?
percent_overlap = 0.5;          % percent of overlap between boxes

[logname] = PD_Fit_saveParams(outDir,degrees,M0_file,Qmap_file,BM_file,...
                            seg_file,outMm,boxSize,percent_overlap,...
                            Inclusion_Criteria, Qmap_factor);


%% III. Perform the fit

load(logname);
fprintf('\n fitting gain ...              \n');

Fit_PD_Qmap(opt);               % use the linear relationship between the map and M0 to fing the gain


%% IV. Build the local fits
%Join the local overlap area to one PD image.
fprintf('\n Building the PD map...              \n');

RepErrThreshold=[];
ErrorThresh=0.01;
PrcCutOff=[];
                              
logname=buildPD_SM(logname,RepErrThreshold,PrcCutOff,ErrorThresh);
%

