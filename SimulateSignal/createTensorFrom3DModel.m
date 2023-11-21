function [chi11,chi12,chi13,chi22,chi23,chi33,phimap]= createTensorFrom3DModel(model,Params)
 
 % for cylinder along z axis, also test effective for other orientation
 
 axonlist = model.axonlist;
 dim = size(model.model);
 phimap = zeros(dim);
 chi = zeros([dim,3,3]);

 myelin_chii = [Params.myelin.chii,0,0;0,Params.myelin.chii,0;0,0,Params.myelin.chii];
 myelin_chia = [Params.myelin.chia,0,0;0,-Params.myelin.chia/2,0;0,0,-Params.myelin.chia/2];
 
 if (isfield(Params, 'extra_axonal') && isfield(Params.extra_axonal, 'chii')) 
    extra_axonal_chii = [Params.extra_axonal.chii 0 0; 0 Params.extra_axonal.chii 0; 0 0 Params.extra_axonal.chii];
 end
 
 if (isfield(Params, 'intra_axonal') && isfield(Params.intra_axonal, 'chii')) 
    intra_axonal_chii = [Params.intra_axonal.chii 0 0; 0 Params.intra_axonal.chii 0; 0 0 Params.intra_axonal.chii];
 end

 sigma= 2;
 smooth_mask = imgaussfilt3(model.model, sigma, 'FilterSize', 3, 'Padding', 'replicate','FilterDomain', 'spatial');
 [~, gradient_phi, gradient_elevation] = imgradient3(smooth_mask);
 phi = (pi/180)*(gradient_phi + 90);   
 gradient_elevation_degree = (pi/180)*gradient_elevation;
 theta = gradient_elevation_degree; 

 %% Cylinder Z model
  for Num = 1:length(axonlist)
     sub_myelin = axonlist(Num).myelin;
     sub_axon = axonlist(Num).axon;

     for k = 1:size(sub_myelin,1)


         phi_rot = mod(phi(sub_myelin(k,1), sub_myelin(k,2),sub_myelin(k,3))+2*pi,2*pi)-pi;
         theta_rot = theta(sub_myelin(k,1), sub_myelin(k,2),sub_myelin(k,3));

         Rz = [sin(phi_rot) cos(phi_rot) 0; cos(phi_rot) -sin(phi_rot) 0; 0 0 1];
         Ry = [cos(theta_rot) 0 sin(theta_rot);0 1 0; -sin(theta_rot) 0 cos(theta_rot)];
         R = Rz*Ry;

         chi(sub_myelin(k,1), sub_myelin(k,2),sub_myelin(k,3),:,:) =  myelin_chii + R*myelin_chia*R';
         phimap(sub_myelin(k,1), sub_myelin(k,2),sub_myelin(k,3)) = phi_rot;

     end
     if (isfield(Params, 'intra_axonal') && isfield(Params.intra_axonal, 'chii'))
        for k = 1:size(sub_axon,1)
            chi(sub_axon(k,1),sub_axon(k,2),sub_axon(k,3),:,:) = intra_axonal_chii;
        end
     end

  end

   chi11 = chi(:,:,:,1,1);
   chi12 = chi(:,:,:,1,2);
   chi13 = chi(:,:,:,1,3);
   chi22 = chi(:,:,:,2,2);
   chi23 = chi(:,:,:,2,3);
   chi33 = chi(:,:,:,3,3);
