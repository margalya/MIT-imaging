
function appData = createConsts(appData)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% physical constants + atomic constants
appData.consts.AtomType.str = {'Sodium', 'Lithium 6'}; %, 'Lithium 6 |1>'};
appData.consts.AtomType.Sodium = 1;
appData.consts.AtomType.Lithium6 = 2;
% appData.consts.AtomType.Lithium6S1 = 3;

appData.consts.wavelength{1} = 589.1583264e-9; %Sodium D2 wavelength [m]
appData.consts.wavelength{2} = 670.977338e-9; %Lithium 6 D2 line wavelength [m]

appData.consts.scatcross0{1} = 3*(589.1583264e-9)^2/2/pi; %Sodium D2 line resonante scattering cross section [m^2]
appData.consts.scatcross0{2} = 3*(670.977338e-9)^2/2/pi; %Lithium 6 D2 line resonant scattering cross section [m^2]
appData.consts.hbar = 1.054571628e-34; %Planck's reduced constant, [J*sec]

appData.consts.MNa = 3.817e-26; %sodium mass [Kg]
appData.consts.MLi6 = 9.9883414e-27; %Lithium 6 mass [Kg]
appData.consts.Kb = 1.38064852e-23; % Boltzmann constant [J/K]
appData.consts.linew{1} = 9.7946e6; % Sodium D2 linewidth [Hz]
appData.consts.linew{2} = 5.8724e6; % Lithium 6 D2 linewidth [Hz]
appData.consts.Isat{1} = 6.2600*10; % Sodium saturation intensity for sigma+ light I_sat(mF =±2 -> mF'=±3) [W/m^2]
appData.consts.Isat{2} = 2.54*10; % Lithium 6 saturation intensity for sigma+ light I_sat(mF =±2 -> mF'=±3) [W/m^2]
                                                    % line width of the Rb87 cooling transition: Lambda = 2*pi*linew
% appData.consts.scatcross0 = 3*appData.consts.wavelength^2/2/pi;   % scattering cross section for Rb87 in m^2
% appData.consts.Isat = 16.6933; %35.7713 [W/m^2] Saturation Intensity, for isotropic light field. For sigma+/- light is it 16.6933 [W/m^2]
% appData.consts.scatcross = appData.consts.scatcross0;      %  resonant imaging 
% appData.consts.defaultDetuning = 0;
% consts.scatcross = scatcross0 * 1/(1+(detun*2/linew)^2);   % off resonance scattering cross section


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% defaults

appData.consts.plotTypes.default = 1;
appData.consts.ROIUnits.default = 2;
appData.consts.ROIUnits.defaultSizeX = 0.5;
appData.consts.ROIUnits.defaultSizeY = 0.5;
appData.consts.calcAtomsNo.default = 1;
appData.consts.plotSetting.default = 1;
appData.consts.fitTypes.default = 2;%12;
appData.consts.saveParams.default = 1;
appData.consts.saveParamValDefault = 0;

appData.consts.defaultAvgWidth = 2; %on each side
appData.consts.maxAvgWidth = 8;
appData.save.defaultDir = ['F:\My Documents\Experimental\' datestr(now, 29)];%'C:\Documents and Settings\broot\Desktop\shimi';
appData.consts.defaultStrLVFile_Save = 'F:\My Documents\Experimental\LVData_Save.txt'; %'D:\My Documents\MATLAB\lv data.txt';
appData.consts.defaultStrLVFile_Load = 'F:\My Documents\Experimental\LVData_Load.txt'; %'D:\My Documents\MATLAB\lv data.txt';
appData.consts.defaultDetuning = 0;

appData.consts.winXPos = 1;
appData.consts.winYPos = 1;
% appData.consts.ftpAddress = '132.72.8.2';
% appData.consts.ftpDir = 'LabViewProjects/Ramp_Integration/Interface39.5/TXT/';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loops
appData.consts.availableLoops.str = {'Gen. Loop' 'TOF', 'LT' 'SG'};%, 'Heating', 'SG'};
appData.consts.availableLoops.generalLoop = 1;
appData.consts.availableLoops.TOF = 2;
appData.consts.availableLoops.LT = 3;
% appData.consts.availableLoops.heating = 4;
appData.consts.availableLoops.SG = 4;
appData.consts.availableLoops.createFncs = {@GeneralLoop.create, @TOF.create, @LT.create, @SternGerlach.create};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loops defaults
appData.consts.loops.iterationsStr = {'Iterate Measurement' 'Iterate Loop' 'Random Iterations'};
appData.consts.loops.options = struct('Interpreter', 'tex', 'WindowStyle', 'normal', 'Resize', 'off');

appData.consts.loops.TOF.noIterations = '1';
appData.consts.loops.TOF.TOFTimes = '5:2:13';

appData.consts.loops.LT.noIterations = '1';
appData.consts.loops.LT.LTTimes = '[0:1:20]*1e3';

appData.consts.loops.SG.noIterations = '1';
appData.consts.loops.SG.SGTimes = '[0:1:20]*1e3';
appData.consts.loops.SG.RFRampNo = 15;
% appData.consts.loops.SG.RFRampNo = 18;

appData.consts.loops.GenLoop.saveFolder = 'F:\My Documents\Experimental\Loops';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cameras

appData.consts.cameraTypes.str = {'PixelFly', 'PixelFly X axis', 'Andor', 'BlacfFly'};
appData.consts.cameraTypes.pixelfly = 1; %plug axis, 1:1 magnification
appData.consts.cameraTypes.pixelflyXaxis = 2; % x axis - side bucket window, magnification = 244.6/43.4
appData.consts.cameraTypes.andor = 3;
appData.consts.cameraTypes.blackfly = 4;
appData.consts.cameraTypes.default = 1;

% check pixelfly parameters
appData.consts.cameras{appData.consts.cameraTypes.pixelfly}.string = 'PixelFly';
appData.consts.cameras{appData.consts.cameraTypes.pixelfly}.dir = 'C:\Documents and Settings\broot\Desktop\shimi';
appData.consts.cameras{appData.consts.cameraTypes.pixelfly}.fileName = 'tmpPic';
appData.consts.cameras{appData.consts.cameraTypes.pixelfly}.fileFormat = 'TIF';
appData.consts.cameras{appData.consts.cameraTypes.pixelfly}.fileReadFunction = @pixelflyReadFunction;
appData.consts.cameras{appData.consts.cameraTypes.pixelfly}.darkPicStr = '';
appData.consts.cameras{appData.consts.cameraTypes.pixelfly}.magnification = 1;
appData.consts.cameras{appData.consts.cameraTypes.pixelfly}.xPixSz = ...
    6.45e-6 / appData.consts.cameras{appData.consts.cameraTypes.pixelfly}.magnification;
appData.consts.cameras{appData.consts.cameraTypes.pixelfly}.yPixSz = ...
    6.45e-6 / appData.consts.cameras{appData.consts.cameraTypes.pixelfly}.magnification;
appData.consts.cameras{appData.consts.cameraTypes.pixelfly}.width = 1024;
appData.consts.cameras{appData.consts.cameraTypes.pixelfly}.height = 1392;
appData.consts.cameras{appData.consts.cameraTypes.pixelfly}.rotate = 0;
appData.consts.cameras{appData.consts.cameraTypes.pixelfly}.bits = 14;
appData.consts.cameras{appData.consts.cameraTypes.pixelfly}.chipStart = 1;

appData.consts.cameras{appData.consts.cameraTypes.pixelflyXaxis}.string = 'PixelFly X axis';
appData.consts.cameras{appData.consts.cameraTypes.pixelflyXaxis}.dir = 'C:\Documents and Settings\broot\Desktop\shimi';
appData.consts.cameras{appData.consts.cameraTypes.pixelflyXaxis}.fileName = 'tmpPic';
appData.consts.cameras{appData.consts.cameraTypes.pixelflyXaxis}.fileFormat = 'TIF';
appData.consts.cameras{appData.consts.cameraTypes.pixelflyXaxis}.fileReadFunction = @pixelflyReadFunction;
appData.consts.cameras{appData.consts.cameraTypes.pixelflyXaxis}.darkPicStr = '';
appData.consts.cameras{appData.consts.cameraTypes.pixelflyXaxis}.magnification = 250/49.4; %ration of achromat lenses' back focal lengths
appData.consts.cameras{appData.consts.cameraTypes.pixelflyXaxis}.xPixSz = ...
    6.45e-6 / appData.consts.cameras{appData.consts.cameraTypes.pixelflyXaxis}.magnification;
appData.consts.cameras{appData.consts.cameraTypes.pixelflyXaxis}.yPixSz = ...
    6.45e-6 / appData.consts.cameras{appData.consts.cameraTypes.pixelflyXaxis}.magnification;
appData.consts.cameras{appData.consts.cameraTypes.pixelflyXaxis}.width = 1024;
appData.consts.cameras{appData.consts.cameraTypes.pixelflyXaxis}.height = 1392;
appData.consts.cameras{appData.consts.cameraTypes.pixelflyXaxis}.rotate = 0;
appData.consts.cameras{appData.consts.cameraTypes.pixelflyXaxis}.bits = 14;
appData.consts.cameras{appData.consts.cameraTypes.pixelflyXaxis}.chipStart = 1;


appData.consts.cameras{appData.consts.cameraTypes.andor}.string = 'Andor';
appData.consts.cameras{appData.consts.cameraTypes.andor}.dir = 'C:\Documents and Settings\broot\Desktop\shimi';
appData.consts.cameras{appData.consts.cameraTypes.andor}.fileName = 'tmpPic';
appData.consts.cameras{appData.consts.cameraTypes.andor}.fileFormat = 'FITS';
appData.consts.cameras{appData.consts.cameraTypes.andor}.fileReadFunction = @andorReadFunction;
appData.consts.cameras{appData.consts.cameraTypes.andor}.darkPicStr = '';
appData.consts.cameras{appData.consts.cameraTypes.andor}.magnification = 41/49.4; %f2/f1 = 41/49.4 First lens is OptoSigma Visible Spectrum Achromat 49.4mm Focal Length, second lens is Special optics 54-36-41@671&1064nm, EFL = 41 mm 
appData.consts.cameras{appData.consts.cameraTypes.andor}.xPixSz = ...
    6.5e-6 / appData.consts.cameras{appData.consts.cameraTypes.andor}.magnification;
appData.consts.cameras{appData.consts.cameraTypes.andor}.yPixSz = ...
    6.5e-6 / appData.consts.cameras{appData.consts.cameraTypes.andor}.magnification;
appData.consts.cameras{appData.consts.cameraTypes.andor}.width = 2560;
appData.consts.cameras{appData.consts.cameraTypes.andor}.height = 2160;
appData.consts.cameras{appData.consts.cameraTypes.andor}.rotate = 0;
appData.consts.cameras{appData.consts.cameraTypes.andor}.bits = 16;
appData.consts.cameras{appData.consts.cameraTypes.andor}.chipStart = 1;
appData.consts.cameras{appData.consts.cameraTypes.andor}.firstImageNo = 1;
appData.consts.cameras{appData.consts.cameraTypes.andor}.secondImageNo = 2;
appData.consts.cameras{appData.consts.cameraTypes.andor}.photonPerADU = 1.39/0.97 ;% 0.28 / 0.55; % Andor photon per ADU, calcualted as Gain/QE. Gain has units of e-/ADU, and QE has units of e-/photon

% check pixelfly parameters
appData.consts.cameras{appData.consts.cameraTypes.blackfly}.string = 'Blackfly';
appData.consts.cameras{appData.consts.cameraTypes.blackfly}.dir = 'C:\Documents and Settings\broot\Desktop\shimi';
appData.consts.cameras{appData.consts.cameraTypes.blackfly}.fileName = 'tmpPic';
appData.consts.cameras{appData.consts.cameraTypes.blackfly}.fileFormat = 'TIF';
appData.consts.cameras{appData.consts.cameraTypes.blackfly}.fileReadFunction = @blackflyReadFunction;
appData.consts.cameras{appData.consts.cameraTypes.blackfly}.darkPicStr = '';
appData.consts.cameras{appData.consts.cameraTypes.blackfly}.magnification = 41/49.4;
appData.consts.cameras{appData.consts.cameraTypes.blackfly}.xPixSz = ...
    6.9e-6 / appData.consts.cameras{appData.consts.cameraTypes.blackfly}.magnification;
appData.consts.cameras{appData.consts.cameraTypes.blackfly}.yPixSz = ...
    6.9e-6 / appData.consts.cameras{appData.consts.cameraTypes.blackfly}.magnification;
appData.consts.cameras{appData.consts.cameraTypes.blackfly}.width = 540;
appData.consts.cameras{appData.consts.cameraTypes.blackfly}.height = 720;
appData.consts.cameras{appData.consts.cameraTypes.blackfly}.rotate = 0;
appData.consts.cameras{appData.consts.cameraTypes.blackfly}.bits = 12;
appData.consts.cameras{appData.consts.cameraTypes.blackfly}.chipStart = 1;
appData.consts.cameras{appData.consts.cameraTypes.blackfly}.photonPerADU = 0.35/0.50./100; % Blackfly photon per ADU, calcualted as Gain/QE. Gain has units of e-/ADU, and QE has units of e-/photon. Etra factor 100 is for the 40dB Gain mode of the Blackfly

% appData.consts.cameras{appData.consts.cameraTypes.idsMain}.string = 'IDS - main';
% appData.consts.cameras{appData.consts.cameraTypes.idsMain}.dir =  'C:\Program Files\IDS\uEye\Samples'; %'D:\My Documents\MATLAB\imaging v4.3';%
% appData.consts.cameras{appData.consts.cameraTypes.idsMain}.fileName = 'uEye_Image_';
% appData.consts.cameras{appData.consts.cameraTypes.idsMain}.fileFormat = 'bmp';
% appData.consts.cameras{appData.consts.cameraTypes.idsMain}.fileReadFunction = @idsReadFunction;
% appData.consts.cameras{appData.consts.cameraTypes.idsMain}.darkPicStr = 'dark.bmp';
% appData.consts.cameras{appData.consts.cameraTypes.idsMain}.magnification = 300/200;%100/300;
% appData.consts.cameras{appData.consts.cameraTypes.idsMain}.xPixSz = ...
%     6e-6 / appData.consts.cameras{appData.consts.cameraTypes.idsMain}.magnification;
% appData.consts.cameras{appData.consts.cameraTypes.idsMain}.yPixSz = ...
%     6e-6 / appData.consts.cameras{appData.consts.cameraTypes.idsMain}.magnification;
% appData.consts.cameras{appData.consts.cameraTypes.idsMain}.width = 480;
% appData.consts.cameras{appData.consts.cameraTypes.idsMain}.height = 752;
% appData.consts.cameras{appData.consts.cameraTypes.idsMain}.rotate = 90;
% appData.consts.cameras{appData.consts.cameraTypes.idsMain}.bits = 8;
% appData.consts.cameras{appData.consts.cameraTypes.idsMain}.chipStart = 40;
% appData.consts.cameras{appData.consts.cameraTypes.idsMain}.firstImageNo = 1;
% appData.consts.cameras{appData.consts.cameraTypes.idsMain}.secondImageNo = 2;
% % consts.cameras{consts.cameraTypes.ids} = ids;
% 
% appData.consts.cameras{appData.consts.cameraTypes.idsSecond}.string = 'IDS - second';
% appData.consts.cameras{appData.consts.cameraTypes.idsSecond}.dir = 'C:\Program Files\IDS\uEye\Samples';%'F:\Jonathan\Documents\MATLAB\';
% appData.consts.cameras{appData.consts.cameraTypes.idsSecond}.fileName = 'uEye_Image_';
% appData.consts.cameras{appData.consts.cameraTypes.idsSecond}.fileFormat = 'bmp';
% appData.consts.cameras{appData.consts.cameraTypes.idsSecond}.fileReadFunction = @idsReadFunction;
% appData.consts.cameras{appData.consts.cameraTypes.idsSecond}.darkPicStr = 'dark.bmp';
% appData.consts.cameras{appData.consts.cameraTypes.idsSecond}.magnification = 1/2.55;
% appData.consts.cameras{appData.consts.cameraTypes.idsSecond}.xPixSz = ...
%     6e-6 / appData.consts.cameras{appData.consts.cameraTypes.idsSecond}.magnification;
% appData.consts.cameras{appData.consts.cameraTypes.idsSecond}.yPixSz = ...
%     6e-6 / appData.consts.cameras{appData.consts.cameraTypes.idsSecond}.magnification;
% appData.consts.cameras{appData.consts.cameraTypes.idsSecond}.width = 480;
% appData.consts.cameras{appData.consts.cameraTypes.idsSecond}.height = 752;
% appData.consts.cameras{appData.consts.cameraTypes.idsSecond}.rotate = 90;
% appData.consts.cameras{appData.consts.cameraTypes.idsSecond}.bits = 8;
% appData.consts.cameras{appData.consts.cameraTypes.idsSecond}.chipStart = 1;
% appData.consts.cameras{appData.consts.cameraTypes.idsSecond}.firstImageNo = 1;
% appData.consts.cameras{appData.consts.cameraTypes.idsSecond}.secondImageNo = 2;
% % consts.cameras{consts.cameraTypes.ids} = ids;
% 
% appData.consts.cameras{appData.consts.cameraTypes.prosilica}.string = 'prosilica';
% appData.consts.cameras{appData.consts.cameraTypes.prosilica}.dir = 'C:\Program Files\Prosilica';
% appData.consts.cameras{appData.consts.cameraTypes.prosilica}.fileName = 'snap';
% appData.consts.cameras{appData.consts.cameraTypes.prosilica}.fileFormat = 'tiff';
% appData.consts.cameras{appData.consts.cameraTypes.prosilica}.fileReadFunction = @prosilicaReadFunction;
% appData.consts.cameras{appData.consts.cameraTypes.prosilica}.darkPicStr = 'dark.bmp';
% appData.consts.cameras{appData.consts.cameraTypes.prosilica}.magnification = 300/200;%100/300;
% appData.consts.cameras{appData.consts.cameraTypes.prosilica}.xPixSz = ...
%     3.45e-6 / appData.consts.cameras{appData.consts.cameraTypes.prosilica}.magnification;
% appData.consts.cameras{appData.consts.cameraTypes.prosilica}.yPixSz = ...
%     3.45e-6 / appData.consts.cameras{appData.consts.cameraTypes.prosilica}.magnification;
% appData.consts.cameras{appData.consts.cameraTypes.prosilica}.width = 2448;
% appData.consts.cameras{appData.consts.cameraTypes.prosilica}.height = 1050;
% appData.consts.cameras{appData.consts.cameraTypes.prosilica}.rotate = 0;
% appData.consts.cameras{appData.consts.cameraTypes.prosilica}.bits = 16;
% appData.consts.cameras{appData.consts.cameraTypes.prosilica}.chipStart =102;%82 ;%was 131 before hight calibration
% appData.consts.cameras{appData.consts.cameraTypes.prosilica}.firstImageNo = 0;
% appData.consts.cameras{appData.consts.cameraTypes.prosilica}.secondImageNo = 1;
% 
% appData.consts.cameras{appData.consts.cameraTypes.prosilicaC}.string = 'prosilica-c';
% appData.consts.cameras{appData.consts.cameraTypes.prosilicaC}.dir = 'C:\Documents and Settings\broot\Desktop\shimi\prosilica';
% appData.consts.cameras{appData.consts.cameraTypes.prosilicaC}.fileName = 'snap';
% appData.consts.cameras{appData.consts.cameraTypes.prosilicaC}.fileFormat = 'tiff';
% appData.consts.cameras{appData.consts.cameraTypes.prosilicaC}.fileReadFunction = @prosilicaCReadFunction;
% appData.consts.cameras{appData.consts.cameraTypes.prosilicaC}.darkPicStr = 'dark.bmp';
% appData.consts.cameras{appData.consts.cameraTypes.prosilicaC}.magnification = 3.27;%100/300;
% appData.consts.cameras{appData.consts.cameraTypes.prosilicaC}.xPixSz = ...
%     3.45e-6 / appData.consts.cameras{appData.consts.cameraTypes.prosilicaC}.magnification;
% appData.consts.cameras{appData.consts.cameraTypes.prosilicaC}.yPixSz = ...
%     3.45e-6 / appData.consts.cameras{appData.consts.cameraTypes.prosilicaC}.magnification;
% % appData.consts.cameras{appData.consts.cameraTypes.prosilicaC}.width = 480;
% % appData.consts.cameras{appData.consts.cameraTypes.prosilicaC}.height = 752;
% appData.consts.cameras{appData.consts.cameraTypes.prosilicaC}.rotate = 90;
% appData.consts.cameras{appData.consts.cameraTypes.prosilicaC}.bits = 16;
% appData.consts.cameras{appData.consts.cameraTypes.prosilicaC}.chipStart = 1;
% appData.consts.cameras{appData.consts.cameraTypes.prosilicaC}.firstImageNo = 1;
% appData.consts.cameras{appData.consts.cameraTypes.prosilicaC}.secondImageNo = 2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

appData.consts.calcAtomsNo.str = {'Real', 'Theoretical', 'Theoretical - Full'};
appData.consts.calcAtomsNo.real = 1;
appData.consts.calcAtomsNo.theoretical = 2;
appData.consts.calcAtomsNo.theoreticalFull = 3;
appData.consts.calcAtomsNo.calcs = { CalcReal CalcTheoretical CalcTheoreticalFull};%, Theoretical};
appData.consts.plotSetting.str = {'Default Setting', 'Last Setting'};
appData.consts.plotSetting.defaults = 1;
appData.consts.plotSetting.last = 2;

appData.consts.fitTypes.str = {'XY Cut', 'Only Maximum', '1D Gaussian', '1D Gaussian (Only Y)', '2-1D Gaussian (Only Y)', '2D Gaussian', ...
    '1D TF', '1D TF (Only Y)', '2D TF', '1D BiModal', '2D BiModal', 'Stern-Gerlach', 'Fermi Radial', '2D Fermi', 'Fringes (only Y)', 'Custom Fit'};
appData.consts.fitTypes.XYCut = 1;
appData.consts.fitTypes.onlyMaximum = 2;
appData.consts.fitTypes.oneDGaussian = 3;
appData.consts.fitTypes.oneDGaussianOnlyY = 4;
appData.consts.fitTypes.twoOneDGaussianOnlyY = 5;
appData.consts.fitTypes.twoDGaussian = 6;
appData.consts.fitTypes.oneDTF = 7;
appData.consts.fitTypes.oneDTFOnlyY = 8;
appData.consts.fitTypes.twoDTF = 9;
appData.consts.fitTypes.oneDBiModal = 10;
appData.consts.fitTypes.twoDBiModal = 11;
appData.consts.fitTypes.SG = 12;
appData.consts.fitTypes.FermiRadial = 13;
appData.consts.fitTypes.Fermi2D = 14;
appData.consts.fitTypes.fringesY = 15;
appData.consts.fitTypes.customFit = 16;
appData.consts.fitTypes.fits = {FitXYCut, FitOnlyMax, Fit1DGaussian, Fit1DGaussianOnlyY, ...
    Fit21DGaussianOnlyY, Fit2DGaussian, FitThomasFermi1D, FitThomasFermi1DOnlyY, FitThomasFermi2D, ...
    FitBiModal1D, FitBiModal2D, FitSternGerlach, FitFermiRadial, Fit2DFermi, FitFringesY, FitCustom};

appData.consts.plotTypes.str = {'Absorption', 'ROI', 'With Atoms', 'Without Atoms', 'Dark'};
appData.consts.plotTypes.absorption = 1;
appData.consts.plotTypes.ROI = 2;
appData.consts.plotTypes.withAtoms = 3;
appData.consts.plotTypes.withoutAtoms = 4;
appData.consts.plotTypes.dark = 5;
appData.consts.plotTypes.plots = {Absorption, ROI, WithAtoms, WithoutAtoms, Dark};

appData.consts.ROIUnits.str = {'Sigma (Sx Sy)', 'mm (Sx Sy)', 'Size [mm] (Sx Sy Cx Cy)'};
appData.consts.ROIUnits.sigma = 1;
appData.consts.ROIUnits.mm = 2;
appData.consts.ROIUnits.size = 3;
appData.consts.ROIUnits.ROITypes = {Sigma, MM, Size};

% appData.consts.availableMonitoring.str = {'Atoms Num.', 'X Position', 'Y Position'}; % when adding, change appData.monitoring.monitoringData size
% appData.consts.availableMonitoring.atomNum = 1;
% appData.consts.availableMonitoring.xPos = 2;
% appData.consts.availableMonitoring.yPos = 3;

appData.consts.saveParams.str = {'Other...', 'TOF [ms]', 'Dark Time [ms]', 'X-Bias', 'RF freq [MHz]','Osc time [ms]', 'Other Param'};
appData.consts.saveParams.other = 1;
appData.consts.saveParams.TOF = 2;
appData.consts.saveParams.darkTime = 3;
appData.consts.saveParams.XBias = 4;
appData.consts.saveParams.RFfreq = 5;
appData.consts.saveParams.OscTime = 6;
appData.consts.saveParams.otherParam = 7; %ALWAYS the last
appData.consts.saveOtherParamStr = '';
appData.consts.commentStr = '';

appData.consts.availableAnalyzing.str = {'Temperature', 'Temperature 1 Shot', 'Gravity', 'Life Time (1 exp)', 'Life Time (2 exp)',...
    'Two-body loss', 'Three-body loss', 'Atom No', 'Photon Count', 'Photons Per Atom', 'OD', 'X Position', 'Y Position', 'Size X', 'Size Y', 'Aspect Ratio', 'Aspect Ratio + fit', 'T/T_F','BEC fraction',...
    'Delta_y', 'Pic Mean','SG', 'mF1', 'mF1 Rabi Freq', 'SG Y Position', 'lambda', 'phi', ...
    'Visibility', 'Norm. Vis.', 'Chirp', 'Save Param Val', 'Analyze Folders', 'Principal Axes', 'Image parameters', 'Class property'};
appData.consts.availableAnalyzing.temperature = 1;
appData.consts.availableAnalyzing.temperatureSingleShot = 2;
appData.consts.availableAnalyzing.gravity = 3;
appData.consts.availableAnalyzing.lifeTime1 = 4;
appData.consts.availableAnalyzing.lifeTime2 = 5;
appData.consts.availableAnalyzing.twoBody = 6;
appData.consts.availableAnalyzing.threeBody = 7;
appData.consts.availableAnalyzing.atomNo = 8;
appData.consts.availableAnalyzing.photonCount = 9;
appData.consts.availableAnalyzing.photonsPerAtom = 10;
appData.consts.availableAnalyzing.OD = 11;
appData.consts.availableAnalyzing.xPos = 12;
appData.consts.availableAnalyzing.yPos = 13;
appData.consts.availableAnalyzing.sizeX = 14;
appData.consts.availableAnalyzing.sizeY = 15;
appData.consts.availableAnalyzing.aspectRatio = 16;
appData.consts.availableAnalyzing.aspectRatioExpFit = 17;
appData.consts.availableAnalyzing.ToverTF = 18;
appData.consts.availableAnalyzing.BECfraction = 19;
appData.consts.availableAnalyzing.deltaY_2 = 20;
appData.consts.availableAnalyzing.picMean = 21;
appData.consts.availableAnalyzing.SG = 22;
appData.consts.availableAnalyzing.mF1 = 23;
appData.consts.availableAnalyzing.mF1RabiFreq = 24;
appData.consts.availableAnalyzing.SGyPos = 25;
appData.consts.availableAnalyzing.lambda = 26;
appData.consts.availableAnalyzing.phi = 27;
appData.consts.availableAnalyzing.visibility = 28;
appData.consts.availableAnalyzing.NormVis = 29;
appData.consts.availableAnalyzing.Chirp = 30;
appData.consts.availableAnalyzing.SaveParamVal = 31;
appData.consts.availableAnalyzing.analyzeFolders = 32;
appData.consts.availableAnalyzing.principalAxes = 33;
appData.consts.availableAnalyzing.imageParams = 34;
appData.consts.availableAnalyzing.classProperty = 35;

appData.consts.pmLineSpec.str = {'o', 'o-', '-', 'o, mean+/-std', 'o, mean+/-sem', 'Average Param'};
appData.consts.pmLineSpec.o = 1;
appData.consts.pmLineSpec.oHyphen = 2;
appData.consts.pmLineSpec.Hyphen = 3;
appData.consts.pmLineSpec.oMeanStd = 4;
appData.consts.pmLineSpec.oMeanSem = 5;
appData.consts.pmLineSpec.AverageParam = 6;
% appData.consts.pmLineSpec.SaveToWS = 7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function file = idsReadFunction(appData, num)
if ( strcmp(appData.data.camera.string, appData.consts.cameraTypes.str{appData.consts.cameraTypes.idsMain}) ~= 1 && ...
        strcmp(appData.data.camera.string, appData.consts.cameraTypes.str{appData.consts.cameraTypes.idsSecond}) ~= 1 )
    warndlg('Trying to read IDS file.', 'Warning', 'nonmodal');
end
fileName = [appData.data.camera.dir '\' appData.data.camera.fileName  num2str(num) '.' appData.data.camera.fileFormat];
file = imread(fileName, appData.data.camera.fileFormat);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function file = prosilicaReadFunction(appData, num)
if ( strcmp(appData.data.camera.string, appData.consts.cameraTypes.str{appData.consts.cameraTypes.prosilica}) ~= 1 )
    warndlg('Trying to read Prosilica file.', 'Warning', 'nonmodal');
end
fileName = [appData.data.camera.dir '\' appData.data.camera.fileName  num2str(num) '.' appData.data.camera.fileFormat];
file = imread(fileName, appData.data.camera.fileFormat)/2^4;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function file = prosilicaCReadFunction(appData, num)
if ( strcmp(appData.data.camera.string, appData.consts.cameraTypes.str{appData.consts.cameraTypes.prosilicaC}) ~= 1 )
    warndlg('Trying to read Prosilica-c file.', 'Warning', 'nonmodal');
end
fileName = [appData.data.camera.dir '\' appData.data.camera.fileName  num2str(num) '.' appData.data.camera.fileFormat];
file = imread(fileName, appData.data.camera.fileFormat);
file=file(:,:,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [atoms, back, dark] = andorReadFunction(fileName)
% if ( strcmp(appData.data.camera.string, appData.consts.cameraTypes.str{appData.consts.cameraTypes.andor}) ~= 1)
%     warndlg('Trying to read Andor file.', 'Warning', 'nonmodal');
% end
% fileName = [appData.data.camera.dir '\' appData.data.camera.fileName  num2str(num) '.' appData.data.camera.fileFormat];
% file = imread(fileName);
% fid = fopen(fileName);
% file = fread(fid, [512 512], 'uint32');
% fclose(fid);
% fitsread
data = fitsread(fileName);
% clean = data(:,:,1); % sensor cleaning image
atoms = data(:,:,2);
back = [];
dark = data(:,:,3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [atoms, back, dark] = pixelflyReadFunction(fileName)
% if ( strcmp(appData.data.camera.string, appData.consts.cameraTypes.str{appData.consts.cameraTypes.pixelfly}) ~= 1)
%     warndlg('Trying to read pixelfly file.', 'Warning', 'nonmodal');
% end

numimgs = size(imfinfo(fileName),1); %number of images to read. 2 = fluorescence; 3 = absoprtion
% numimgs=2;
% imdata = imdata[:,:,0] + 256*imdata[:,:,1] + 65536*imdata[:,:,2]
% inverting xCamera software data compression(?) for tiff, taken from Aviv's Odysseus code
atoms = uint16(imread(fileName, 'tiff' , 1)); atoms = atoms(:,:,1) + 256.*atoms(:,:,2) + 65536.*atoms(:,:,3);
if numimgs == 2 %fluorescence
    back =[]; %send empty output
    dark = uint16(imread(fileName, 'tiff' , 2)); dark = dark(:,:,1) + 256.*dark(:,:,2) + 65536.*dark(:,:,3);
elseif numimgs ==3 %absoprtion
    back = uint16(imread(fileName, 'tiff' , 2)); back = back(:,:,1) + 256.*back(:,:,2) + 65536.*back(:,:,3);
    dark = uint16(imread(fileName, 'tiff' , 3)); dark = dark(:,:,1) + 256.*dark(:,:,2) + 65536.*dark(:,:,3);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [atoms, back, dark] = blackflyReadFunction(fileName)
% assumes that camera is used for flurescence imaging.
% Shot #1 is for sensor cleaning, and is omited

% numimgs = size(imfinfo(fileName),1); 
% numimgs=2;
atoms = uint16(imread(fileName, 'tiff' , 2)); 
dark = uint16(imread(fileName, 'tiff' , 3));
back = [];
% if numimgs == 2 %fluorescence
%     back =[]; %send empty output
%     dark = uint16(imread(fileName, 'tiff' , 2));
% elseif numimgs ==3 %absoprtion
%     back = uint16(imread(fileName, 'tiff' , 2));
%     dark = uint16(imread(fileName, 'tiff' , 3));
% end

end


end


