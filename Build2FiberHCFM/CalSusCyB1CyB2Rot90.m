% Calculate susceptibility tensor map for two fiber bundle crossing at 50-90° in a step of 5°.
% In the fiber_1 frame
% Fiber 2 can be set in the y-z plane 
clc,clear all
pathname_main = '~/Desktop/WM/code/Exampledata/HCFM/HCFM2FiberSus';
t = 0:0.1:2*pi+0.1;

theta = 90;
phi = 0;
Params.B0 = 7;              % Tesla
Params.gamma = 42.577e6;      % Hz/T
H0 = [sind(theta)*cosd(phi) sind(theta)*sind(phi) cosd(theta)]';

% myelin (required: T2, xi, xa)
Params.myelin.Mag= 0.7; 
Params.myelin.T2 = 8*1e-3;
Params.myelin.T1 = 500*1e-3;
Params.myelin.proton_density = 1; 

Params.myelin.chii = -0.10;  % myelin anisotropic susceptibility (ppm)
Params.myelin.chia = -0.10;  % myelin isotropic susceptibility (ppm)

% intra axonal (required: T2)
Params.intra_axonal.Mag = 1;
Params.intra_axonal.T2 = 36*1e-3;
Params.intra_axonal.T1 = 1.5;
Params.intra_axonal.proton_density= 1; 
Params.intra_axonal.chii= 0; 

% extra axonal (required: T2)
Params.extra_axonal.Mag = 1;
Params.extra_axonal.T2 = 36*1e-3;
Params.extra_axonal.T1 = 1.5;
Params.extra_axonal.proton_density = 1; 
Params.extra_axonal.chii = 0; 
Params.E = 0.02;

%% rotate bundle 2
 Rot_Center = 1;             % 0: Rotation around the Origin, 1: Rotation around the Center
 
 RotAngArray = [0, 0, 0;
                -5, 0, 0;
                -10, 0, 0;
                -15, 0, 0;
                -20, 0, 0;
                -25, 0, 0;
                -30, 0, 0;
                -35, 0, 0;
                -40, 0, 0;]; % RL, AP, FH
 OriNum = size(RotAngArray,1);
 % creat cubic mask
 [xgrid,ygrid,zgrid] = meshgrid(1:200,1:200,1:200);

 vob.Length = 160;
 vob.CenterX = 100;
 vob.CenterY = 100;
 vob.CenterZ = 100;
 Mask = VObjCube(vob,xgrid,ygrid,zgrid);


for SimInd = 1:OriNum
    load('~/Desktop/WM/code/Exampledata/HCFM/3DCylinderBundel1.mat')
    PhantomCy1.model = model;
    PhantomCy1.FOV = FOV;
    PhantomCy1.axonlist = axonlist;

    [Cy1_chi11,Cy1_chi12,Cy1_chi13,Cy1_chi22,Cy1_chi23,Cy1_chi33,~] = createTensorFrom3DModel(PhantomCy1,Params); %Bundel1 Chi

    load('~/Desktop/WM/code/Exampledata/HCFM/3DCylinderBundel2Rot90.mat')
    PhantomCy2.model = model;
    PhantomCy2.FOV = FOV;
    PhantomCy2.axonlist = axonlist;

   [Cy2_chi11,Cy2_chi12,Cy2_chi13,Cy2_chi22,Cy2_chi23,Cy2_chi33,~] = createTensorFrom3DModel(PhantomCy2,Params); %Bundel2 Chi
    in_ref = imref3d(size(PhantomCy2.model));
    
    if Rot_Center == 1
  % ------------  Set original to the center, for in_ref
      XWorldLimits = in_ref.XWorldLimits - in_ref.ImageExtentInWorldX/2;
      YWorldLimits = in_ref.YWorldLimits - in_ref.ImageExtentInWorldY/2;
      ZWorldLimits = in_ref.ZWorldLimits - in_ref.ImageExtentInWorldZ/2;
      in_ref = imref3d(size(PhantomCy2.model), XWorldLimits, YWorldLimits, ZWorldLimits);
    end
      t3d = eye(4);
      t3dTAng = Rmatrix_ang(-RotAngArray(SimInd, 1), -RotAngArray(SimInd, 2), -RotAngArray(SimInd, 3)); % to match affine3d
      t3d(1:3,1:3) = t3dTAng; 
      tform3d = affine3d(t3d);        % affine transformation

      [Cy2_chi11_tform, out_ref11] = imwarp(Cy2_chi11, in_ref, tform3d, 'OutputView', in_ref);
      [Cy2_chi12_tform, out_ref12] = imwarp(Cy2_chi12, in_ref, tform3d, 'OutputView', in_ref);
      [Cy2_chi13_tform, out_ref13] = imwarp(Cy2_chi13, in_ref, tform3d, 'OutputView', in_ref);
      [Cy2_chi22_tform, out_ref22] = imwarp(Cy2_chi22, in_ref, tform3d, 'OutputView', in_ref);
      [Cy2_chi23_tform, out_ref23] = imwarp(Cy2_chi23, in_ref, tform3d, 'OutputView', in_ref);
      [Cy2_chi33_tform, out_ref33] = imwarp(Cy2_chi33, in_ref, tform3d, 'OutputView', in_ref);
      [PhantomCy2.model_tform, out_ref] = imwarp(PhantomCy2.model-1, in_ref, tform3d, 'nearest','OutputView', in_ref);
      PhantomCy2.model_tform = PhantomCy2.model_tform + 1;
      Cy1_chi11 = Cy1_chi11(:,:,51:250);
      Cy1_chi12 = Cy1_chi12(:,:,51:250);
      Cy1_chi13 = Cy1_chi13(:,:,51:250);
      Cy1_chi22 = Cy1_chi22(:,:,51:250);
      Cy1_chi23 = Cy1_chi23(:,:,51:250);
      Cy1_chi33 = Cy1_chi33(:,:,51:250);


      Cy2_chi11_tform = Cy2_chi11_tform(51:250,:,:);
      Cy2_chi12_tform = Cy2_chi12_tform(51:250,:,:);
      Cy2_chi13_tform = Cy2_chi13_tform(51:250,:,:);
      Cy2_chi22_tform = Cy2_chi22_tform(51:250,:,:);
      Cy2_chi23_tform = Cy2_chi23_tform(51:250,:,:);
      Cy2_chi33_tform = Cy2_chi33_tform(51:250,:,:);

      total_chi = zeros(200,200,200,3,3);

      total_chi(:,1:100,:,1,1) = Cy1_chi11;
      total_chi(:,101:200,:,1,1) = Cy2_chi11_tform;

      total_chi(:,1:100,:,1,2) = Cy1_chi12;
      total_chi(:,101:200,:,1,2) = Cy2_chi12_tform;

      total_chi(:,1:100,:,1,3) = Cy1_chi13;
      total_chi(:,101:200,:,1,3) = Cy2_chi13_tform;

      total_chi(:,1:100,:,2,2) = Cy1_chi22;
      total_chi(:,101:200,:,2,2) = Cy2_chi22_tform;

      total_chi(:,1:100,:,2,3) = Cy1_chi23;
      total_chi(:,101:200,:,2,3) = Cy2_chi23_tform;

      total_chi(:,1:100,:,3,3) = Cy1_chi33;
      total_chi(:,101:200,:,3,3) = Cy2_chi33_tform;

      total_model = zeros(200,200,200);
      total_model(:,1:100,:) = PhantomCy1.model(:,:,51:250);
      total_model(:,101:200,:) = PhantomCy2.model_tform(51:250,:,:);

       Phantom.model = total_model;
       Phantom.mask = Mask;
       Phantom.FOV = FOV;

      SaveName = [pathname_main,'/317Fiber2BundleCrossAngle',num2str(90 - (SimInd-1)*5),'.mat'];

      save(SaveName,'total_chi','Phantom','Mask');
end

     