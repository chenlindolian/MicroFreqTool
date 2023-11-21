%% create one model including 156/161 random distributed hollow cylinders along Z axis

% no overlapping regions, one point belongs to one hollow cylinder
% simulate 3d parallel fiber regions 
% build your own cylinder bundle: https://github.com/neuropoly/axonpacking

addpath('~/Desktop/WM/code/SupportFunction/')
load('~/Desktop/WM/code/Exampledata/HCFM/cylinderbundle.mat');
Xsize = 100;
Ysize = 200;
Zsize = 300;
model = ones(Ysize,Xsize,Zsize);
FOV = [0.2 0.2 0.2];
[xgrid,ygrid,zgrid] = meshgrid(1:Xsize,1:Ysize,1:Zsize);

FaceNum = 100;
Length = 300;
centerZ(1:numberofCylinder) = 150;
% Center X, Center Y, Center Z

%% for Z cylinder
for id = 1:numberofCylinder
    [Xm,Ym,Zm] = cylinder(r_outer(id),FaceNum); 
    [Xa,Ya,Za] = cylinder(r_outer(id)*gratio(id),FaceNum); 
    CXm(:,:,id)=Xm+center(1,id);
    CYm(:,:,id)=Ym+center(2,id);
    CZm(:,:,id)=(Zm-0.5)*Length+centerZ(id);
    CXa(:,:,id)=Xa+center(1,id);
    CYa(:,:,id)=Ya+center(2,id);
    CZa(:,:,id)=(Za-0.5)*Length+centerZ(id);
%     CXm(:,:,id) = Xm+center(1,id);
%     CYm(:,:,id) = (Ym - 0.5)*Length + centerZ(id);
%     CZm(:,:,id) = Zm + center(2,id);
%     CXa(:,:,id) = Xa + center(1,id);
%     CYa(:,:,id) = (Ya - 0.5)*Length + centerZ(id);
%     CZa(:,:,id) = Za + center(2,id);
    fvc1m = surf2patch(CXm(:,:,id),CYm(:,:,id),CZm(:,:,id));
    fvc1a = surf2patch(CXa(:,:,id),CYa(:,:,id),CZa(:,:,id));
    mask1m = vert2mask(fvc1m.vertices,xgrid,ygrid,zgrid);
    mask1a = vert2mask(fvc1a.vertices,xgrid,ygrid,zgrid);
    mask = mask1m+mask1a;
    [row,col,slice] = ind2sub(size(mask),find(mask == 1)); %% find y,x,z
     axonlist(id).myelin(:,1) = row;
     axonlist(id).myelin(:,2) = col;
     axonlist(id).myelin(:,3) = slice;
    [row,col,slice] = ind2sub(size(mask),find(mask == 2)); %% find y,x,z
     axonlist(id).axon(:,1) = row;
     axonlist(id).axon(:,2) = col;
     axonlist(id).axon(:,3) = slice;
     axonlist(id).r_mean = r_outer(id); 
     axonlist(id).g_ratio = gratio(id);
     model = model+mask;
end

for n = 1:numberofCylinder
    figure(1)
    hold on
    grid on;
    axis equal;
    axis([0 100 1 200 0 300]);
    view(-50,30);
    s1 = surf(CXm(:,:,n),CYm(:,:,n),CZm(:,:,n),'FaceColor','#000000');
    s1.EdgeColor = 'none';
    % s1.LineStyle = 'none';
    hold on
    s2 = surf(CXa(:,:,n),CYa(:,:,n),CZa(:,:,n),'FaceColor','#80B3FF');
    s2.EdgeColor = 'none';
    set(gca,'XDir','Normal','YDir','Normal','ZDir','Normal'),

    xlabel('{\itx}','FontSize',25,'LineWidth',2,'fontname','Times New Roman')
    ylabel('{\ity}','FontSize',25,'LineWidth',2,'fontname','Times New Roman')
    zlabel('{\itz}','FontSize',25,'LineWidth',2,'fontname','Times New Roman')
    set(gca,'LineWidth',2,'FontSize',25);
end
