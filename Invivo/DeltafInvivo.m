%% TE-dependent frequency fitting with regularization for in vivo data 

clc,clear all
addpath('~/Desktop/WM/code/SupportFunction/nifti/')
addpath('~/Desktop/WM/code/SupportFunction/')

load('data.mat'); %load your own unwrapphase and magnitude image
load('brainmask.mat') % load brain mask

% remove phase offset
PhaseOffsetFit = ((PhaseUnwrp(:,:,:,2)*Params.TEs(1) - PhaseUnwrp(:,:,:,1)*Params.TEs(2))/(Params.TEs(1)-Params.TEs(2)));
   for iecho = 1:length(Params.TEs)
       PhaseUnwrp(:,:,:,iecho) = PhaseUnwrp(:,:,:,iecho) - PhaseOffsetFit;
   end

% do LBV
Params.lbv.tol = 0.001;     % default 0.1
Params.lbv.depth = -1;      % default
Params.lbv.peel = 0;        % similar to mask erosion
PhaseLocal = zeros(size(PhaseUnwrp));
Background = zeros(size(PhaseUnwrp));
BgRemovalMethods = 'LBV';
Params.echoNums = 1:Params.nEchoes;

% unit to Hz
FreqLocal = zeros(size(PhaseLocal));
FreqBg = zeros(size(PhaseLocal));
Freq = zeros(size(PhaseLocal));
 for selectedEcho = 1:length(Params.TEs)   
     TEsSE = Params.TEs(selectedEcho);
     Freq(:,:,:,selectedEcho) = PhaseUnwrp(:,:,:,selectedEcho)./(2*pi*TEsSE);  % in Hz = 
 end 

for echoii = 1:length(Params.echoNums)
    FreqLocal(:,:,:,Params.echoNums(echoii)) = LBV(Freq(:,:,:,Params.echoNums(echoii)),maskErode, Params.sizeVol, Params.voxSize, ...
                                                                Params.lbv.tol, Params.lbv.depth, Params.lbv.peel);
    Background(:,:,:,Params.echoNums(echoii)) = Freq(:,:,:,Params.echoNums(echoii)) - FreqLocal(:,:,:,Params.echoNums(echoii));
end

GREMag = GREMag.*maskErode;
weightmap = AverageEchoWeight(GREMag,Params.TEs);
% %% load gibbs ringing removal local freq maps
% nii = load_untouch_nii('LocalFreq_unring_echo1.nii.gz');
% FreqLocal(:,:,:,1) = permute(nii.img,[2,1,3]);
% 
% nii = load_untouch_nii('LocalFreq_unring_echo2.nii.gz');
% FreqLocal(:,:,:,2) = permute(nii.img,[2,1,3]);
% 
% nii = load_untouch_nii('LocalFreq_unring_echo3.nii.gz');
% FreqLocal(:,:,:,3) = permute(nii.img,[2,1,3]);
% 
% nii = load_untouch_nii('LocalFreq_unring_echo4.nii.gz');
% FreqLocal(:,:,:,4) = permute(nii.img,[2,1,3]);

maxSNRecho = 14;
weightsum = AverageEchoWeight(GREMag(:,:,:,6:15),Params.TEs(6:15));
FreqLocal_average = sum(FreqLocal(:,:,:,6:15).*weightsum,4);

% mimage(FreqLocal(:,:,55,:),-15,15,3,7)

%% calculate curve for frequency difference map
x0 = [1 100 0 0];
ub = [10 1000 1 +Inf];
lb = [-10 80 0 -Inf];
[Ny, Nx, Nz, NEcho] = size(FreqLocal);
TEs = Params.TEs;

BrainDifferenceFreq = zeros(size(maskErode));
curvep1 = zeros(size(maskErode));
curvep2 = zeros(size(maskErode));
Cfmap = zeros(size(maskErode));
resmap = zeros(size(maskErode));

% options.Algorithm = 'levenberg-marquardt';
options.Algorithm = 'trust-region-reflective';
options.MaxFunEvals = 8000;
options.MaxIter = 1000;
lambda = 0.01;
parfor kk = 1:Nz
    for ii = 1:Ny
       for jj = 1:Nx
            if maskErode(ii,jj,kk)
            FreqSig = squeeze(FreqLocal(ii,jj,kk,:))';
            FreqSig_avg = FreqLocal_average(ii,jj,kk);
            weight= squeeze(weightmap(ii,jj,kk,:))';
            weightedFun = ...  
                        @(x) [reshape((Freqmodel(x, TEs) - FreqSig).*weight,[],1); lambda*(x(4)-FreqSig_avg)];   
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

save('BrainMicroFreqRomeoPhaseoffsetScaleLBVReg0p005.mat','BrainDifferenceFreq','FreqLocal','maskErode','curvep1','curvep2','Cfmap','resmap','FreqLocal_average');

function ydata = Freqmodel(x,iTE)
         ydata = x(1).*(exp(-x(2)*iTE)./(x(3)*(exp(-x(2)*iTE)-1)+1))+x(4);
end