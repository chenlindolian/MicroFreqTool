%function [chiavg, chiani, chiparl, chiperp, eigvs] = chitensor2chiaapp(chi11, chi12, chi13, chi22, chi23, chi33, pathname)
function [chiavg, chiani, chiparl, chiperp] = chitensor2chiaapp(chi11, chi12, chi13, chi22, chi23, chi33)
% [chiavg, chiani, eigvs] = chitensor2chiaa(chi11, chi12, chi13, chi22, chi23, chi33, pathname)
% 
%% Author: Xu Li
% Affiliation: Radiology @ JHU
% Email address: xuli@mri.jhu.edu
% 
% calculate the the chi_average and
% chi_anisotropy from full susceptibility tensor . Assuming cylindrical symmetric susceptibility tensor.
%
%
% OUTPUT:
% 
% chiavg:  average susceptibility (1/3*trace)
% chiani:  susceptibility anisotropy (chi_parallel - chi_perpendicular)
% eigvs:    structure defining the eigen vectors 
%           should have eigvs.eig0, eigvs.eig1, eigvs.eig2, which are
%           filenames referenced to eigenvector data files exported from DTI
%
% INPUT:
% chi_ij:   full susceptibility tensor components (symmetric tensor)

[Ny, Nx, Nz] = size(chi11);

chiavg = zeros(size(chi11));
chiani = zeros(size(chi11));
chiparl = zeros(size(chi11));
chiperp = zeros(size(chi11));

% eigvs = struct('eig0', [pathname, 'EigenVec_0.dat'], ...
%                 'eig1', [pathname, 'EigenVec_1.dat'], ...
%                   'eig2', [pathname, 'EigenVec_2.dat']);
% 
% eig0 = zeros(3, Ny, Nx, Nz);
% eig1 = eig0;
% eig2 = eig0;
              
% loop through each voxel
h = waitbar(0, 'calculating average and anisotropy...');
for iz = 1:Nz
    for ix = 1:Nx
        for iy = 1:Ny
            
            chiT = [chi11(iy, ix, iz), chi12(iy, ix, iz), chi13(iy, ix, iz); ...
                    chi12(iy, ix, iz), chi22(iy, ix, iz), chi23(iy, ix, iz); ...
                    chi13(iy, ix, iz), chi23(iy, ix, iz), chi33(iy, ix, iz)];
            [V, D] = eig(chiT);
            
            chidiag = diag(D);
            [B, Ind] = sort(chidiag);
            
            chiavg(iy, ix, iz) = 1/3*sum(B);        % average
                
            chiani(iy, ix, iz) = B(end) - B(1);     % anisotropy
            
            chiparl(iy, ix, iz) = B(end);           % parallel
            
            chiperp(iy, ix, iz) = (B(1) + B(2))./2; % perpendicular
            
            eig0(:, iy, ix, iz) = V(:, Ind(end));   % priciple eigen vector
            eig1(:, iy, ix, iz) = V(:, Ind(end-1));
            eig2(:, iy, ix, iz) = V(:, Ind(end-2));
            
        end
    end
    
    waitbar(iz/Nz);
end
close(h)

%% save eigenvectors
% fid = fopen(eigvs.eig0, 'w', 'l');
% fwrite(fid, eig0(:), 'float');
% fclose(fid);
% 
% fid = fopen(eigvs.eig1, 'w', 'l');
% fwrite(fid, eig1(:), 'float');
% fclose(fid);
% 
% fid = fopen(eigvs.eig2, 'w', 'l');
% fwrite(fid, eig2(:), 'float');
% fclose(fid);
