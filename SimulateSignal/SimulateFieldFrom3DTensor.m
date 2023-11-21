 function deltaB = SimulateFieldFrom3DTensor(chi, model, Params, H0)

 dims = size(model.model);
 Ny = dims(1);
 Nx = dims(2);
 Nz = dims(3);

 chi11 = chi(:,:,:,1,1);
 chi12 = chi(:,:,:,1,2);
 chi13 = chi(:,:,:,1,3);
 chi22 = chi(:,:,:,2,2);
 chi23 = chi(:,:,:,2,3);
 chi33 = chi(:,:,:,3,3);

 [kx,ky,kz] = meshgrid(-Nx/2:Nx/2-1, -Ny/2:Ny/2-1, -Nz/2:Nz/2-1);

 kx = (kx / model.FOV(2));
 ky = (ky / model.FOV(1));
 kz = (kz / model.FOV(3));

 KX_Grid = single(fftshift(kx));
 KY_Grid = single(fftshift(ky));
 KZ_Grid = single(fftshift(kz));
 
 KSq = KX_Grid.^2 + KY_Grid.^2 + KZ_Grid.^2; 
 
 h1 = H0(1);                                        
 h2 = H0(2);
 h3 = H0(3);

 KHdKSq = (KX_Grid*h1 + KY_Grid*h2 + KZ_Grid*h3)./(KSq+eps);

 a11 = h1.^2/3 - KHdKSq.*KX_Grid*h1;
 a22 = h2.^2/3 - KHdKSq.*KY_Grid*h2;
 a33 = h3.^2/3 - KHdKSq.*KZ_Grid*h3;

 a12 = 2*h1.*h2/3 - KHdKSq.*(KX_Grid*h2 + KY_Grid*h1);
 a13 = 2*h1.*h3/3 - KHdKSq.*(KX_Grid*h3 + KZ_Grid*h1);
 a23 = 2*h2.*h3/3 - KHdKSq.*(KY_Grid*h3 + KZ_Grid*h2); 
 
 chi11k = fftn(chi11);                     % tensor components in k space
 chi12k = fftn(chi12);
 chi13k = fftn(chi13);
 chi22k = fftn(chi22);
 chi23k = fftn(chi23);
 chi33k = fftn(chi33);
 
 deltaBk = a11.*chi11k + a12.*chi12k + a13.*chi13k + ...
            a22.*chi22k + a23.*chi23k + a33.*chi33k ;
        
 deltaB = real(ifftn(deltaBk)).*Params.gamma.*Params.B0*1e-6; 

 signal_model = model.model;
   
 myelin_index = (signal_model == 2); 

 deltaB = deltaB + Params.E*myelin_index*Params.gamma.*Params.B0*1e-6;
end