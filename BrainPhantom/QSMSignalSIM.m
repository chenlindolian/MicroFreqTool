
%% simulate brain phantom signal and fit deltaf,p1,p2,Cf using TE-dependent frequency fitting method
clc,clear all
addpath('~/Desktop/WM/code/SupportFunction/nifti/')
addpath('~/Desktop/WM/code/SupportFunction/')
addpath('~/Desktop/WM/code/Exampledata/HCFM/')
addpath('~/Desktop/WM/code/Exampledata/ExampleImage/')
addpath('~/Desktop/WM/code/SimulateSignal/')

load('BrainSimMacroSus.mat') %
load('BrainSimMicroFreq.mat')
load('PhantomMask.mat');

%% simulate signal voxel by voxel

Params.B0 = 7;              % Tesla
Params.gamma = 42.577e6;      % Hz/T

Params.fov = [224, 224, 126];        % mm     
Params.sizeVol = [224, 224, 126];  

Params.voxSize = Params.fov./Params.sizeVol;    % mm

chi11 = Simchitotal(:,:,:,1,1);
chi12 = Simchitotal(:,:,:,1,2);
chi13 = Simchitotal(:,:,:,1,3);
chi22 = Simchitotal(:,:,:,2,2);
chi23 = Simchitotal(:,:,:,2,3);
chi33 = Simchitotal(:,:,:,3,3);

chi11_tissue = Simchi(:,:,:,1,1);
chi12_tissue = Simchi(:,:,:,1,2);
chi13_tissue = Simchi(:,:,:,1,3);
chi22_tissue = Simchi(:,:,:,2,2);
chi23_tissue = Simchi(:,:,:,2,3);
chi33_tissue = Simchi(:,:,:,3,3);

[chiavg, chiani, chiparl, chiperp] = chitensor2chiaapp(chi11, chi12, chi13, chi22, chi23, chi33);

[chiavg_tissue, ~, ~, ~] = chitensor2chiaapp(chi11_tissue, chi12_tissue, chi13_tissue, chi22_tissue, chi23_tissue, chi33_tissue);

%% Forward Simulation

OriAngle = '20-40';
OriAngleRLarray = [0, -20, -20, -20, 40, 40, 40];
OriAngleAParray = [0, 0, 0, 0, 0, 0, 0];
OriAngleFHarray = [0, 0, 120, 240, 0, 120, 240];

H0Lab = [0, 0, 1]';
H0Sub = zeros(3, length(OriAngleRLarray));

OrientInd = 1;    
Params.AngRL = OriAngleRLarray(OrientInd);
Params.AngAP = OriAngleAParray(OrientInd);
Params.AngFH = OriAngleFHarray(OrientInd);

Params.TAng = Rmatrix_arb(Params.AngRL, Params.AngAP, Params.AngFH);     % rotation matrix

%% QSM forward method
% D = conv_kernel_rot_c0(Params, Params.TAng);
% 
% D = fftshift(D);
% 
% deltaB_total_QSM = real(ifftn(fftn(chiavg).*D)).*Params.gamma.*Params.B0*1e-6;  % total chi
% deltaB_tissue_QSM = real(ifftn(fftn(chiavg_tissue).*D)).*Params.gamma.*Params.B0*1e-6;  % tissue chi
%% STI forward method
H0Sub(:,OrientInd) = Params.TAng'*H0Lab;
H0 = H0Sub(:,OrientInd);
deltaB_total = SimulateFieldFrom3DTensor(Simchitotal, Params, H0); % unit Hz
deltaB_tissue = SimulateFieldFrom3DTensor(Simchi, Params, H0); % unit Hz

Params.voxSize = [1 1 1];
Params.TEs = [1:1:30]*1e-3;
Params.SNR = 100;
PhaseOffset = 1; %multiplier term of a quadratic phase over the brain

% PhaseOffset
[c , w ] = centerofmass(Eve_GREMag0);
[y,x,z] = meshgrid([1:Params.sizeVol(2)]-c(2), [1:Params.sizeVol(1)]-c(1), [1:Params.sizeVol(3)]-c(3));
temp = (x/w(1)).^2 + (y/w(2)).^2 + (z/w(3)).^2 ;
PhaseOffset = - temp/(max(temp(BrainMask==1))-min(temp(BrainMask==1)))*pi*PhaseOffset;

 % simulate MRI signal

 Params.TR = 0.075; %75ms
 Params.theta = 15; %flip angle

 Signal = zeros([Params.sizeVol, length(Params.TEs)]);
 deltaB_micro = zeros([Params.sizeVol, length(Params.TEs)]);

 testMicro = 0;
 testMacro = 0;
 testtotal = 1;

 for iecho = 1:length(Params.TEs)
     deltaB_micro(:,:,:,iecho) = BrainMicroFreq(:,:,:,iecho)-Cfmap;
     Eve_GREMag1 = Eve_GREMag0.*exp(-Params.TEs(1).*EveR2star).*(1-exp(-Params.TR.*Eve_R1)).*sind(Params.theta)./(1-cosd(Params.theta).*exp(-Params.TR.*Eve_R1)); % add R1 and R2* effect
     minGREMag = min(Eve_GREMag1(:));
     maxGREMag = max(Eve_GREMag1(:));
     if testMicro == 1
        Eve_GREMag = Eve_GREMag0.*exp(-Params.TEs(iecho).*EveR2star).*(1-exp(-Params.TR.*Eve_R1)).*sind(Params.theta)./(1-cosd(Params.theta).*exp(-Params.TR.*Eve_R1)); % add R1 and R2* effect
        Eve_GREMag = (Eve_GREMag - minGREMag)/(maxGREMag - minGREMag); %normalize for SNR
        Signal(:,:,:,iecho) = Eve_GREMag.*exp(1i*(2*pi*Params.TEs(iecho)*deltaB_micro(:,:,:,iecho)+PhaseOffset))...
                             +1/Params.SNR*((randn(Params.sizeVol))+1i*randn(Params.sizeVol))*max(Eve_GREMag(:)); % add noise in real and imaginary part
     end
     if testMacro == 1
        Eve_GREMag = Eve_GREMag0.*exp(-Params.TEs(iecho).*EveR2star).*(1-exp(-Params.TR.*Eve_R1)).*sind(Params.theta)./(1-cosd(Params.theta).*exp(-Params.TR.*Eve_R1));% add R1 and R2* effect
        Eve_GREMag = (Eve_GREMag - minGREMag)/(maxGREMag - minGREMag);  %normalize for SNR
        Signal(:,:,:,iecho) = Eve_GREMag.*exp(1i*(2*pi*Params.TEs(iecho)*deltaB_total+PhaseOffset))...
                              +1/Params.SNR*((randn(Params.sizeVol))+1i*randn(Params.sizeVol))*max(Eve_GREMag(:)); % add noise in real and imaginary part
     end
     if testtotal == 1
        Eve_GREMag = Eve_GREMag0.*exp(-Params.TEs(iecho).*EveR2star).*(1-exp(-Params.TR.*Eve_R1)).*sind(Params.theta)./(1-cosd(Params.theta).*exp(-Params.TR.*Eve_R1));% add R1 and R2* effect
        Eve_GREMag = (Eve_GREMag - minGREMag)/(maxGREMag - minGREMag); %normalize for SNR
        Signal(:,:,:,iecho) = Eve_GREMag.*exp(1i*(2*pi*Params.TEs(iecho)*(deltaB_total+deltaB_micro(:,:,:,iecho))+PhaseOffset))...
                             +1/Params.SNR*((randn(Params.sizeVol))+1i*randn(Params.sizeVol))*max(Eve_GREMag(:)); % add noise in real and imaginary part
     end   
 end

 GREMag = abs(Signal); 
 GREPhase = angle(Signal); 

 %% phase unwrapping
GREPhase = GREPhase(:,:,:,2:2:30); %acquire 15 echoes data
GREMag = GREMag(:,:,:,2:2:30);
Params.TEs = [2:2:30]*1e-3;
Params.nEchoes = 15;
Params.nDynamics = 1;
PhaseUnwrp = zeros(size(GREPhase));
GREPhase = double(GREPhase);
TV = TVOP;
PhaseUnwrapMethods = 'Romeo';

% remove phase offset
PhaseOffsetIndex = 1;
if PhaseOffsetIndex == 1
    DeltaPhase = angle(GREMag(:,:,:,2).*exp(1i*GREPhase(:,:,:,2)).*...
        (GREMag(:,:,:,1).*exp(-1i*GREPhase(:,:,:,1))));  
    dTE = Params.TEs(2) - Params.TEs(1);    % deltaTE
   if mod(Params.TEs(1), dTE) ~= 0
    % do unwrapping on DeltaPhase
        DeltaPhase = phase_unwrap_path_mex(DeltaPhase);  
   end
   PhaseOffsetFit= angle(exp(1i*GREPhase(:,:,:,1)).*exp(-1i*Params.TEs(1)/dTE*DeltaPhase));
   
   GREPhase = angle(GREMag.*exp(1i.*GREPhase).*GREMag.*exp(-1i.*PhaseOffsetFit)); 
end
% 
%  % remove unrelymask and phase unwrapping
switch PhaseUnwrapMethods
    case 'Path'
      for echo_ind = 1:length(Params.TEs)
         PhaseUnwrp(:,:,:,echo_ind) = phase_unwrap_path_mex(GREPhase(:,:,:,echo_ind));       
     end
   case 'Laplacian'
     for echo_ind = 1:length(Params.TEs)
         PhaseUnwrp(:,:,:,echo_ind) = phase_unwrap_laplacian(GREPhase(:,:,:,echo_ind), Params, [], 2);      
     end 
    case 'Romeo'
       Params.mag = GREMag;
       Params.mask = double(BrainMask);

       Params.calculate_B0 = false;
       Params.phase_offset_correction = 'off';

       Params.additional_flags = {'-v'};

       [PhaseUnwrp, ~] = ROMEO(GREPhase, Params);

    otherwise
      PhaseUnwrp = GREPhase;
end

MaskErode = maskErode;

% unit to Hz
 Freq = zeros(size(PhaseUnwrp));
 for iEcho = 1:length(Params.TEs)   
     TEsSE = Params.TEs(iEcho);
     Freq(:,:,:,iEcho) = PhaseUnwrp(:,:,:,iEcho)./(2*pi*TEsSE);  % in Hz
 end 
 
%background remove
Params.echoNums = 1:15;

Params.lbv.tol = 0.0001;     % default 0.1
Params.lbv.depth = -1;      % default
Params.lbv.peel = 0;        % similar to mask erosion
FreqLocal = zeros(size(Freq));
Background = zeros(size(Freq));
BgRemovalMethods = 'LBV';
switch BgRemovalMethods
      case 'LBV'
          for echoii = 1:length(Params.echoNums)
           FreqLocal(:,:,:,Params.echoNums(echoii)) = LBV(Freq(:,:,:,Params.echoNums(echoii)), MaskErode, Params.sizeVol, Params.voxSize, ...
                                                                Params.lbv.tol, Params.lbv.depth, Params.lbv.peel);
           Background(:,:,:,Params.echoNums(echoii)) = Freq(:,:,:,Params.echoNums(echoii)) - FreqLocal(:,:,:,Params.echoNums(echoii));
          end
       case 'VSHARP'
          radiusArray = 1:15;
          Params.thresh_tsvd = 0.05;
          [FreqLocal, Background, mask_eval] = VSHARP_k(Freq, BrainMask, radiusArray, Params.thresh_tsvd, Params, []);
          multiWaitbar('CloseAll'); 
       case 'PDF'
           for echoii = 1:length(Params.echoNums)
               DPWeight = Eve_GREMag0;
               FreqLocal(:,:,:,Params.echoNums(echoii)) = dipolefit_v2(Freq(:,:,:,Params.echoNums(echoii)), BrainMask, DPWeight, Params);
               Background(:,:,:,Params.echoNums(echoii)) = Freq(:,:,:,Params.echoNums(echoii)) - PhaseLocal(:,:,:,Params.echoNums(echoii));
           end
       case 'LBV+VSHARP'
           for echoii = 1:length(Params.echoNums)
           FreqLocal(:,:,:,Params.echoNums(echoii)) = LBV(Freq(:,:,:,Params.echoNums(echoii)), BrainMask, Params.sizeVol, Params.voxSize, ...
                                                                Params.lbv.tol, Params.lbv.depth, Params.lbv.peel);
           end
           radiusArray = 1:15;
           Params.thresh_tsvd = 0.05;
           [FreqLocal, Background, mask_eval] = VSHARP_k(FreqLocal, BrainMask, radiusArray, Params.thresh_tsvd, Params, []);
           multiWaitbar('CloseAll');
       otherwise
           FreqLocal = Freq;
end


 weightmap = AverageEchoWeight(GREMag,Params.TEs);
 %% calculate curve for frequency difference map
[Ny, Nx, Nz, NEcho] = size(GREMag);
x0 = [0 100 0 0];
ub = [10 1000 1 +Inf];
lb = [-10 80 0 -Inf];
TEs = Params.TEs;

% options.Algorithm = 'levenberg-marquardt';
options.Algorithm = 'trust-region-reflective';
options.MaxFunEvals = 8000;
options.MaxIter = 1000;

BrainDifferenceFreq = zeros(size(MaskErode));
curvep1 = zeros(size(MaskErode));
curvep2 = zeros(size(MaskErode));
Cfmap = zeros(size(MaskErode));
resmap = zeros(size(MaskErode));

parfor kk = 1:Nz
    for ii = 1:Ny
       for jj = 1:Nx

          if MaskErode(ii,jj,kk)

             FreqSig = squeeze(FreqLocal(ii,jj,kk,:))';
             weight= squeeze(weightmap(ii,jj,kk,:))';
             weightedFun = ...  
                        @(x) (Freqmodel(x, TEs) - FreqSig).*weight;
            [x,~,res] = lsqnonlin(weightedFun, x0, lb, ub,options);

            BrainDifferenceFreq(ii,jj,kk) = x(1);
            curvep1(ii,jj,kk) =  x(2);
            curvep2(ii,jj,kk) =  x(3);
            Cfmap(ii,jj,kk) =  x(4);
            resmap(ii,jj,kk) =  sqrt(mean(res.^2));
          end
          
       end
    end
end
save('TotalFreqFitSNR100BrainR1R2RomeoMCPC3DsScaleLBVweightEcho1-15.mat','chiavg_tissue','FreqLocal','deltaB_tissue','GREMag','maskWM','maskCSF','MaskErode','BrainDifferenceFreq','weightmap','GREPhase','curvep1','curvep2','Cfmap','resmap','Params');

function ydata = Freqmodel(x,iTE)
         ydata = x(1).*(exp(-x(2)*iTE)./(x(3)*(exp(-x(2)*iTE)-1)+1))+x(4);
 end
