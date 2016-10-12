% a three dimensional phase field modeling ferroelectric polarization switching
% under external electric fields and mechanical stresses 
% this programme is updated on Dec. 2 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discretion of real space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx = 128;
ny = 2;
nz = 64;
hx = nx;
hy = ny;
hz = nz;
e=(2*pi/nx)^2;
x = linspace(1,hx,nx);  % the discretion of axis X in real space
y = linspace(1,hy,ny);  % the discretion of axis Y in real space
z = linspace(1,hz,nz);  % the discretion of axis Z in real space
fx = pi/nx;    %the size of grid in Fourier space along axis kx
fy = pi/ny;    %the size of grid in Fourier space along axis ky
fz = pi/nz;    %the size of grid in Fourier space along axis kz
kx = fx*[(0:nx/2-1) (-nx/2:-1)];  %the discretion of axis kx in Fourier space
ky = fy*[(0:ny/2-1) (-ny/2:-1)];  %the discretion of axis ky in Fourier space
kz = fz*[(0:nz/2-1) (-nz/2:-1)];  %the discretion of axis kz in Fourier space
% spectral differentiation in FFT algorithum 
vect = ones(nx,ny,nz);
for i=1:1:nx; 
   k2x(i,:,:) =2*(kx(i)*vect(i,:,:));
   xx(i,:,:) = x(i)*vect(i,:,:);
end
for i=1:1:ny;
   k2y(:,i,:) =2*(ky(i)*vect(:,i,:));
   yy(:,i,:) = y(i)*vect(:,i,:);
end
for i=1:1:nz;
   k2z(:,:,i) =2*(kz(i)*vect(:,:,i));
   %zz(:,:,i) = z(i)*vect(:,:,i);
end
%  k2x = gsingle(k2x);
%  k2y = gsingle(k2y);
%  k2z = gsingle(k2z);

kk = k2x.^2 + k2y.^2 + k2z.^2;
kk(1,1,1) = e;
aa = (k2x.^2./kk);
bb = (k2y.^2./kk);
cc = (k2z.^2./kk);
ab = (k2x.*k2y./kk);
ac = (k2x.*k2z./kk);
bc = (k2y.*k2z./kk);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial distribution of composition

comp = 0.50*ones(nx,ny,nz);
comp(:,:,nz/4+1:nz/2)=0.40;
comp(:,:,nz/2+1:nz*3/4)=0.50;
comp(:,:,nz*3/4+1:nz)=0.40;

% some parameteres
% evolving controlling parameters
timetotal  = 10000;
step = 0.02;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%normalized coefficients used in Landau functions for PZT

temperature = 298.0/748.0;
% a20 = 0.61+1.52*comp;
% tc = 0.65+0.35*comp;
a40 = 3.0;
a60 = 1.99;
yita = 1.5;
c_1 = 0.38;
c_m = 0.45;

% normalized coefficients used in gradient energy functions
g11 = 1.6;
g12 = 0.0;
g44 = 0.8;
g44b = 0.8;
% normalized coefficients used in external mechanical stress fields
v = 0.3; % poison ratio
appl = 0; % under the plane strain assumption
applxx = appl ; % applied stress 
applyy = -appl*v/(1-v);    % applied stress
applxy = 0;     % applied stress
misfit_strain11 = -0.01;
misfit_strain22 = -0.01;
misfit_strain12 = 0;
scale_elas = 0.3;
% normalized coefficients in dipole-like electrostatic interaction
e_para = 1;
e_field = 2.213e8;
% normalized electrostrictive coefficients
q11 = 0.047;
q12 = -0.017;
q44 = 0.013;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % It is used to determine the heterogeneous elasticity of multiple
% % inclusions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Green function tensor for infinite cubic anisotropic elasticity at a
% given global coordinate system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x3_star = [0 0 1]; % the given epitaxial orientation (0 0 1)
x1_star = [1 0 0];
x2_star = [0 1 0];
% corresponding components of rotation_matrix 
rotation = [1 0 0;0 1 0;0 0 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%x3_star = [1 1 0]; % the given epitaxial orientation (1 1 0)
%x1_star = [-1 1 0];
%x2_star = [0 0 1];
% corresponding components of rotation_matrix 
%rotation = [-1/sqrt(2) 1/sqrt(2) 0;0 0 1;1/sqrt(2) 1/sqrt(2) 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%x3_star = [1 1 1]; % the given epitaxial orientation (1 1 1)
%x1_star = [-1 1 0];
%x2_star = [-1 -1 2];
% corresponding components of rotation_matrix 
%rotation = [-1/sqrt(2) 1/sqrt(2) 0;-1/sqrt(6) -1/sqrt(6) 2/sqrt(6);...
%        1/sqrt(3) 1/sqrt(3) 1/sqrt(3)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%x3_star = [3 1 1]; % the given epitaxial orientation (3 1 1)
%x1_star = [-2 3 3];
%x2_star = [0 1 -1];
% corresponding components of rotation_matrix 
%rotation = [-2/sqrt(22) 3/sqrt(22) 3/sqrt(22);0 1/sqrt(2) -1/sqrt(2);...
%        3/sqrt(11) 1/sqrt(11) 1/sqrt(11)];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % K matrix
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K1_matrix = rotation.^2;
K2_matrix = [rotation(1,2)*rotation(1,3) rotation(1,3)*rotation(1,1) rotation(1,1)*rotation(1,2);...
        rotation(2,2)*rotation(2,3) rotation(2,3)*rotation(2,1) rotation(2,1)*rotation(2,2);...
        rotation(3,2)*rotation(3,3) rotation(3,3)*rotation(3,1) rotation(3,1)*rotation(3,2)];
K3_matrix = [rotation(2,1)*rotation(3,1) rotation(2,2)*rotation(3,2) rotation(2,3)*rotation(3,3);...
        rotation(3,1)*rotation(1,1) rotation(3,2)*rotation(1,2) rotation(3,3)*rotation(1,3);...
        rotation(1,1)*rotation(2,1) rotation(1,2)*rotation(2,2) rotation(1,3)*rotation(2,3)];
K4_matrix = [rotation(2,2)*rotation(3,3)+rotation(2,3)*rotation(3,2) rotation(2,3)*rotation(3,1)+rotation(2,1)*rotation(3,3) rotation(2,1)*rotation(3,2)+rotation(2,2)*rotation(3,1);...
        rotation(3,2)*rotation(1,3)+rotation(3,3)*rotation(1,2) rotation(3,3)*rotation(1,1)+rotation(3,1)*rotation(1,3) rotation(3,1)*rotation(1,2)+rotation(3,2)*rotation(1,1);...
        rotation(1,2)*rotation(2,3)+rotation(1,3)*rotation(2,2) rotation(1,3)*rotation(2,1)+rotation(1,1)*rotation(2,3) rotation(1,1)*rotation(2,2)+rotation(1,2)*rotation(2,1)];
K_matrix = [K1_matrix 2*K2_matrix;K3_matrix K4_matrix]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scaled elastic moduli for cubic crystal in the original coordinate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C44 = 669;
C11 = 1054;
C12 = 482;
cubic_b =C11-C12-2*C44; 
moduli_star = [C11 C12 C12 0 0 0;C12 C11 C12 0 0 0;C12 C12 C11 0 0 0;...
        0 0 0 C44 0 0;0 0 0 0 C44 0;0 0 0 0 0 C44];
moduli = K_matrix*moduli_star*K_matrix';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% in the original coordinate system
k2x_star = rotation(1,1)*k2x + rotation(2,1)*k2y + rotation(3,1)*k2z;
k2y_star = rotation(1,2)*k2x + rotation(2,2)*k2y + rotation(3,2)*k2z;
k2z_star = rotation(1,3)*k2x + rotation(2,3)*k2y + rotation(3,3)*k2z;
kk_star = k2x_star.^2 + k2y_star.^2 + k2z_star.^2;
kk_star(1,1,1) = e;
func_f = 1./(1+(C12+C44)*(k2x_star.^2./kk_star./(C44 + cubic_b*k2x_star.^2./kk_star) + k2y_star.^2./kk_star./(C44 + cubic_b*k2y_star.^2./kk_star) +...
    k2z_star.^2./kk_star./(C44 + cubic_b*k2z_star.^2./kk_star)));
% if cubic ==1
%     func_f(1,1,1) = C44/(C12 +2*C44);
% end
omega11_star = (1./(C44 + cubic_b*k2x_star.^2./kk_star)-...
    (C12+C44)*(k2x_star.^2./kk_star)./(C44 + cubic_b*k2x_star.^2./kk_star)./(C44 + cubic_b*k2x_star.^2./kk_star).*func_f);
omega22_star = (1./(C44 + cubic_b*k2y_star.^2./kk_star)-...
    (C12+C44)*(k2y_star.^2./kk_star)./(C44 + cubic_b*k2y_star.^2./kk_star)./(C44 + cubic_b*k2y_star.^2./kk_star).*func_f);
omega33_star = (1./(C44 + cubic_b*k2z_star.^2./kk_star)-...
    (C12+C44)*(k2z_star.^2./kk_star)./(C44 + cubic_b*k2z_star.^2./kk_star)./(C44 + cubic_b*k2z_star.^2./kk_star).*func_f);
omega12_star = -(C12+C44)*(k2x_star.*k2y_star./kk_star)./(C44 + cubic_b*k2x_star.^2./kk_star)./(C44 + cubic_b*k2y_star.^2./kk_star).*func_f;
omega21_star = omega12_star;
omega13_star = -(C12+C44)*(k2x_star.*k2z_star./kk_star)./(C44 + cubic_b*k2x_star.^2./kk_star)./(C44 + cubic_b*k2z_star.^2./kk_star).*func_f;
omega31_star = omega13_star;
omega23_star = -(C12+C44)*(k2y_star.*k2z_star./kk_star)./(C44 + cubic_b*k2y_star.^2./kk_star)./(C44 + cubic_b*k2z_star.^2./kk_star).*func_f;
omega32_star = omega23_star;
% in the present coordinate system
% a temporary matrix
matrix11 = rotation(1,1)*omega11_star + rotation(1,2)*omega21_star + rotation(1,3)*omega31_star;
matrix12 = rotation(1,1)*omega12_star + rotation(1,2)*omega22_star + rotation(1,3)*omega32_star;
matrix13 = rotation(1,1)*omega13_star + rotation(1,2)*omega23_star + rotation(1,3)*omega33_star;
matrix21 = rotation(2,1)*omega11_star + rotation(2,2)*omega21_star + rotation(2,3)*omega31_star;
matrix22 = rotation(2,1)*omega12_star + rotation(2,2)*omega22_star + rotation(2,3)*omega32_star;
matrix23 = rotation(2,1)*omega13_star + rotation(2,2)*omega23_star + rotation(2,3)*omega33_star;
matrix31 = rotation(3,1)*omega11_star + rotation(3,2)*omega21_star + rotation(3,3)*omega31_star;
matrix32 = rotation(3,1)*omega12_star + rotation(3,2)*omega22_star + rotation(3,3)*omega32_star;
matrix33 = rotation(3,1)*omega13_star + rotation(3,2)*omega23_star + rotation(3,3)*omega33_star;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omega11 = matrix11*rotation(1,1) + matrix12*rotation(1,2) + matrix13*rotation(1,3);
omega12 = matrix11*rotation(2,1) + matrix12*rotation(2,2) + matrix13*rotation(2,3);
omega13 = matrix11*rotation(3,1) + matrix12*rotation(3,2) + matrix13*rotation(3,3);
%omega21 = matrix21*rotation(1,1) + matrix22*rotation(1,2) + matrix23*rotation(1,3);
omega22 = matrix21*rotation(2,1) + matrix22*rotation(2,2) + matrix23*rotation(2,3);
omega23 = matrix21*rotation(3,1) + matrix22*rotation(3,2) + matrix23*rotation(3,3);
%omega31 = matrix31*rotation(1,1) + matrix32*rotation(1,2) + matrix33*rotation(1,3);
%omega32 = matrix31*rotation(2,1) + matrix32*rotation(2,2) + matrix33*rotation(2,3);
omega33 = matrix31*rotation(3,1) + matrix32*rotation(3,2) + matrix33*rotation(3,3);
% isotropic Green function
% omega11 = ((1-1/2/(1-v)*aa));
% omega22 = ((1-1/2/(1-v)*bb));
% omega33 = ((1-1/2/(1-v)*cc));
% omega23 = -1/2/(1-v)*bc;
% omega32 = omega23;
% omega13 = -1/2/(1-v)*ac;
% omega31 = omega13;
% omega12 = -1/2/(1-v)*ab;
% omega21 = omega12;
clear omega11_star omega12_star omega13_star omega22_star omega21_star omega23_star omega13_star omega32_star omega33_star
clear matrix11 matrix12 matrix13 matrix21 matrix22 matrix23 matrix31 matrix32 matrix33
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% output files defined
nameout1 = 'lro1.txt';
namein1 = 'lro1';
nameout2 = 'lro2.txt';
namein2 = 'lro2';
nameout3 = 'lro3.txt';
namein3 = 'lro3';
for j = 1:1:10;
     ge(j) = num2str(j-1);% Number to string conversion
 end
for j = 1:1:10;
     shi(j) = num2str(j-1);
 end
for k = 1:1:10;
    bai(k) = num2str(k-1);
end
period = 2000;

% initial distribution of two order parameters
change = 0.001;
lro1 = 0 + change*randn(nx,ny,nz);
lro2 = 0 + change*randn(nx,ny,nz);
lro3 = 0 + change*randn(nx,ny,nz);

e_field_x = zeros(nx,ny,nz);
e_field_y = zeros(nx,ny,nz);
e_field_z = zeros(nx,ny,nz);

direction_x = 1/sqrt(2);
direction_y = 0;
direction_z = 1/sqrt(2);
D(1) = 0;

for t = 1:1:15
    for i = 1:1:timetotal;
    
    k_lro1 = fftn(lro1);
    k_lro2 = fftn(lro2);
    k_lro3 = fftn(lro3);
    
    express11 = (0.61+1.52*comp).*(temperature-(0.65+0.35*comp)).*lro1 + a40*(comp-c_1).*(lro1.^2+lro2.^2+lro3.^2).*lro1- a40*yita*(comp-c_m).*lro1.^3 +a60*(lro1.^2+lro2.^2+lro3.^2).^2.*lro1;
    express21 = (0.61+1.52*comp).*(temperature-(0.65+0.35*comp)).*lro2 + a40*(comp-c_1).*(lro1.^2+lro2.^2+lro3.^2).*lro2- a40*yita*(comp-c_m).*lro2.^3 +a60*(lro1.^2+lro2.^2+lro3.^2).^2.*lro2;
    express31 = (0.61+1.52*comp).*(temperature-(0.65+0.35*comp)).*lro3 + a40*(comp-c_1).*(lro1.^2+lro2.^2+lro3.^2).*lro3- a40*yita*(comp-c_m).*lro3.^3 +a60*(lro1.^2+lro2.^2+lro3.^2).^2.*lro3;
    
    k_express11 = fftn(express11);
    k_express21 = fftn(express21);
    k_express31 = fftn(express31);
    
    k_express12 = g11*k2x.^2.*k_lro1 + g12*(k2x.*k2y.*k_lro2 + k2z.*k2x.*k_lro3) + g44*(k2y.^2.*k_lro1 + k2z.^2.*k_lro1 + k2x.*k2y.*k_lro2 + k2z.*k2x.*k_lro3) + g44b*(k2y.^2.*k_lro1 + k2z.^2.*k_lro1 - k2x.*k2y.*k_lro2 - k2z.*k2x.*k_lro3);
    k_express22 = g11*k2y.^2.*k_lro2 + g12*(k2y.*k2z.*k_lro3 + k2x.*k2y.*k_lro1) + g44*(k2z.^2.*k_lro2 + k2x.^2.*k_lro2 + k2y.*k2z.*k_lro3 + k2x.*k2y.*k_lro1) + g44b*(k2z.^2.*k_lro2 + k2x.^2.*k_lro2 - k2y.*k2z.*k_lro3 - k2x.*k2y.*k_lro1);
    k_express32 = g11*k2z.^2.*k_lro3 + g12*(k2z.*k2x.*k_lro1 + k2y.*k2z.*k_lro2) + g44*(k2x.^2.*k_lro3 + k2y.^2.*k_lro3 + k2z.*k2x.*k_lro1 + k2y.*k2z.*k_lro2) + g44b*(k2x.^2.*k_lro3 + k2y.^2.*k_lro3 - k2z.*k2x.*k_lro1 - k2y.*k2z.*k_lro2);
   
    k_express13 = e_para*(k_lro1.*k2x + k_lro2.*k2y + k_lro3.*k2z).*k2x./kk - fftn(e_field_x);
    k_express23 = e_para*(k_lro1.*k2x + k_lro2.*k2y + k_lro3.*k2z).*k2y./kk - fftn(e_field_y);
    k_express33 = e_para*(k_lro1.*k2x + k_lro2.*k2y + k_lro3.*k2z).*k2z./kk - fftn(e_field_z);
    % the driven force of the contribution of elastic energy release
    % derivation of eigen strain to polarisation
    virt11_1 = 2*q11*lro1; 
    virt22_1 = 2*q12*lro1; 
    virt33_1 = 2*q12*lro1;
    virt12_1 = q44*lro2;   
    virt23_1 = 0;          
    virt31_1 = q44*lro3;
    
    virt11_2 = 2*q12*lro2; 
    virt22_2 = 2*q11*lro2; 
    virt33_2 = 2*q12*lro2;
    virt12_2 = q44*lro1;   
    virt23_2 = q44*lro3;   
    virt31_2 = 0;
    
    virt11_3 = 2*q12*lro3; 
    virt22_3 = 2*q12*lro3; 
    virt33_3 = 2*q11*lro3;
    virt12_3 = 0;          
    virt23_3 = q44*lro2;   
    virt31_3 = q44*lro1;
    
    % eigenstrain due to polarization
    virt11 = q11*lro1.^2 + q12*(lro2.^2 + lro3.^2); 
    virt22 = q11*lro2.^2 + q12*(lro3.^2 + lro1.^2); 
    virt33 = q11*lro3.^2 + q12*(lro1.^2 + lro2.^2); 
    virt23 = q44*lro2.*lro3; 
    virt13 = q44*lro3.*lro1; 
    virt12 = q44*lro1.*lro2;
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % general Hook's law 
 
    sigma11 = moduli(1,1)*fftn(virt11) + moduli(1,2)*fftn(virt22) + moduli(1,3)*fftn(virt33) +...
        moduli(1,4)*fftn(2*virt23) + moduli(1,5)*fftn(2*virt13) +moduli(1,6)*fftn(2*virt12);
    sigma22 = moduli(2,1)*fftn(virt11) + moduli(2,2)*fftn(virt22) + moduli(2,3)*fftn(virt33) +...
        moduli(2,4)*fftn(2*virt23) + moduli(2,5)*fftn(2*virt13) +moduli(2,6)*fftn(2*virt12);
    sigma33 = moduli(3,1)*fftn(virt11) + moduli(3,2)*fftn(virt22) + moduli(3,3)*fftn(virt33) +...
        moduli(3,4)*fftn(2*virt23) + moduli(3,5)*fftn(2*virt13) +moduli(3,6)*fftn(2*virt12);
    sigma23 = moduli(4,1)*fftn(virt11) + moduli(4,2)*fftn(virt22) + moduli(4,3)*fftn(virt33) +...
        moduli(4,4)*fftn(2*virt23) + moduli(4,5)*fftn(2*virt13) +moduli(4,6)*fftn(2*virt12);
    sigma13 = moduli(5,1)*fftn(virt11) + moduli(5,2)*fftn(virt22) + moduli(5,3)*fftn(virt33) +...
        moduli(5,4)*fftn(2*virt23) + moduli(5,5)*fftn(2*virt13) +moduli(5,6)*fftn(2*virt12); 
    sigma12 = moduli(6,1)*fftn(virt11) + moduli(6,2)*fftn(virt22) + moduli(6,3)*fftn(virt33) +...
        moduli(6,4)*fftn(2*virt23) + moduli(6,5)*fftn(2*virt13) +moduli(6,6)*fftn(2*virt12);
    
    
    heterstrain11 = real(ifftn(aa.*omega11.*sigma11 + ab.*omega11.*sigma12 + ac.*omega11.*sigma13 +...
        aa.*omega12.*sigma12 + ab.*omega12.*sigma22 + ac.*omega12.*sigma23 +...
        aa.*omega13.*sigma13 + ab.*omega13.*sigma23 + ac.*omega13.*sigma33));
    heterstrain22 = real(ifftn(ab.*omega12.*sigma11 + bb.*omega12.*sigma12 + bc.*omega12.*sigma13 +...
        ab.*omega22.*sigma12 + bb.*omega22.*sigma22 + bc.*omega22.*sigma23 +...
        ab.*omega23.*sigma13 + bb.*omega23.*sigma23 + bc.*omega23.*sigma33));
    heterstrain33 = real(ifftn(ac.*omega13.*sigma11 + bc.*omega13.*sigma12 + cc.*omega13.*sigma13 +...
        ac.*omega23.*sigma12 + bc.*omega23.*sigma22 + cc.*omega23.*sigma23 +...
        ac.*omega33.*sigma13 + bc.*omega33.*sigma23 + cc.*omega33.*sigma33));
    heterstrain12 = 1/2*real(ifftn(ab.*omega11.*sigma11 + bb.*omega11.*sigma12 + bc.*omega11.*sigma13 +...
        ab.*omega12.*sigma12 + bb.*omega12.*sigma22 + bc.*omega12.*sigma23 +...
        ab.*omega13.*sigma13 + bb.*omega13.*sigma23 + bc.*omega13.*sigma33 +...
        aa.*omega12.*sigma11 + ab.*omega12.*sigma12 + ac.*omega12.*sigma13 +...
        aa.*omega22.*sigma12 + ab.*omega22.*sigma22 + ac.*omega22.*sigma23 +...
        aa.*omega23.*sigma13 + ab.*omega23.*sigma23 + ac.*omega23.*sigma33));
    heterstrain13 = 1/2*real(ifftn(ac.*omega11.*sigma11 + bc.*omega11.*sigma12 + cc.*omega11.*sigma13 +...
        ac.*omega12.*sigma12 + bc.*omega12.*sigma22 + cc.*omega12.*sigma23 +...
        ac.*omega13.*sigma13 + bc.*omega13.*sigma23 + cc.*omega13.*sigma33 +...
        aa.*omega13.*sigma11 + ab.*omega13.*sigma12 + ac.*omega13.*sigma13 +...
        aa.*omega23.*sigma12 + ab.*omega23.*sigma22 + ac.*omega23.*sigma23 +...
        aa.*omega33.*sigma13 + ab.*omega33.*sigma23 + ac.*omega33.*sigma33));
    heterstrain23 = 1/2*real(ifftn(ac.*omega12.*sigma11 + bc.*omega12.*sigma12 + cc.*omega12.*sigma13 +...
        ac.*omega22.*sigma12 + bc.*omega22.*sigma22 + cc.*omega22.*sigma23 +...
        ac.*omega23.*sigma13 + bc.*omega23.*sigma23 + cc.*omega23.*sigma33 +...
        ab.*omega13.*sigma11 + bb.*omega13.*sigma12 + bc.*omega13.*sigma13 +...
        ab.*omega23.*sigma12 + bb.*omega23.*sigma22 + bc.*omega23.*sigma23 +...
        ab.*omega33.*sigma13 + bb.*omega33.*sigma23 + bc.*omega33.*sigma33));
        
%     strain11 = (mean(mean(mean(virt11))) + heterstrain11 - virt11);
%     strain22 = (mean(mean(mean(virt22))) + heterstrain22 - virt22);
%     strain33 = (mean(mean(mean(virt33))) + heterstrain33 - virt33);
%     strain12 = (mean(mean(mean(virt12))) + heterstrain12 - virt12);
%     strain23 = (mean(mean(mean(virt23))) + heterstrain23 - virt23);
%     strain13 = (mean(mean(mean(virt13))) + heterstrain13 - virt13);
    strain11 = (misfit_strain11 + heterstrain11 - virt11);
    strain22 = (misfit_strain22 + heterstrain22 - virt22);
    strain33 = (0 + heterstrain33 - virt33);
    strain12 = (misfit_strain12 + heterstrain12 - virt12);
    strain23 = (0 + heterstrain23 - virt23);
    strain13 = (0 + heterstrain13 - virt13);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    stress11 = moduli(1,1)*(strain11) + moduli(1,2)*(strain22) + moduli(1,3)*(strain33) +...
        moduli(1,4)*(2*strain23) + moduli(1,5)*(2*strain13) +moduli(1,6)*(2*strain12);
    stress22 = moduli(2,1)*(strain11) + moduli(2,2)*(strain22) + moduli(2,3)*(strain33) +...
        moduli(2,4)*(2*strain23) + moduli(2,5)*(2*strain13) +moduli(2,6)*(2*strain12);
    stress33 = moduli(3,1)*(strain11) + moduli(3,2)*(strain22) + moduli(3,3)*(strain33) +...
        moduli(3,4)*(2*strain23) + moduli(3,5)*(2*strain13) +moduli(3,6)*(2*strain12);
    stress23 = moduli(4,1)*(strain11) + moduli(4,2)*(strain22) + moduli(4,3)*(strain33) +...
        moduli(4,4)*(2*strain23) + moduli(4,5)*(2*strain13) +moduli(4,6)*(2*strain12);
    stress13 = moduli(5,1)*(strain11) + moduli(5,2)*(strain22) + moduli(5,3)*(strain33) +...
        moduli(5,4)*(2*strain23) + moduli(5,5)*(2*strain13) +moduli(5,6)*(2*strain12); 
    stress12 = moduli(6,1)*(strain11) + moduli(6,2)*(strain22) + moduli(6,3)*(strain33) +...
        moduli(6,4)*(2*strain23) + moduli(6,5)*(2*strain13) +moduli(6,6)*(2*strain12);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    k_express14 = -scale_elas*fftn(stress11.*virt11_1 + stress22.*virt22_1 + stress33.*virt33_1 + 2*stress12.*virt12_1 + 2*stress23.*virt23_1 + 2*stress13.*virt31_1); 
    k_express24 = -scale_elas*fftn(stress11.*virt11_2 + stress22.*virt22_2 + stress33.*virt33_2 + 2*stress12.*virt12_2 + 2*stress23.*virt23_2 + 2*stress13.*virt31_2);
    k_express34 = -scale_elas*fftn(stress11.*virt11_3 + stress22.*virt22_3 + stress33.*virt33_3 + 2*stress12.*virt12_3 + 2*stress23.*virt23_3 + 2*stress13.*virt31_3); 
    
    
    mid_k_lro1 = k_lro1 - step*(k_express11 + k_express12 + k_express13 + k_express14);
    lro1 = real(ifftn(mid_k_lro1))+ change*randn(nx,ny,nz);
    mid_k_lro2 = k_lro2 - step*(k_express21 + k_express22 + k_express23 + k_express24);
    lro2 = real(ifftn(mid_k_lro2))+ change*randn(nx,ny,nz);
    mid_k_lro3 = k_lro3 - step*(k_express31 + k_express32 + k_express33 + k_express34);
    lro3 = real(ifftn(mid_k_lro3))+ change*randn(nx,ny,nz);   
 end
    
    % record the lro1  profile
       name = [bai(round((t-mod(t,100))/100)+1) shi(round((mod(t,100)-mod(t,10))/10)+1) ge(mod(t,10)+1) nameout1];% 
       nameb = [bai(round((t-mod(t,100))/100)+1) shi(round((mod(t,100)-mod(t,10))/10)+1) ge(mod(t,10)+1) namein1];% 
       fid = fopen(name,'w');
       fprintf(fid,'TITLE = "Phase field variable for lro1" \n');
       fprintf(fid,strcat('VARIABLES=','"',nameb,'"'));
       fprintf(fid, '\n');
       string1 = strcat('zone i=',num2str(nx),',', 'j=',num2str(ny),',', 'k=', num2str(nz));
       fprintf(fid, string1);
       fprintf(fid, '\n');
        for ii=1:1:nz;
            for jj=1:1:ny;           
                    fprintf(fid, '%10.4f',double(lro1(:,jj,ii)));  
                    fprintf(fid, '\n');
            end
         end
        fclose(fid);
        
         % record the lro2  profile
       name = [bai(round((t-mod(t,100))/100)+1) shi(round((mod(t,100)-mod(t,10))/10)+1) ge(mod(t,10)+1) nameout2];% 
       nameb = [bai(round((t-mod(t,100))/100)+1) shi(round((mod(t,100)-mod(t,10))/10)+1) ge(mod(t,10)+1) namein2];% 
       fid = fopen(name,'w');
       fprintf(fid,'TITLE = "Phase field variable for lro2" \n');
       fprintf(fid,strcat('VARIABLES=','"',nameb,'"'));
       fprintf(fid, '\n');
       string1 = strcat('zone i=',num2str(nx),',', 'j=',num2str(ny),',', 'k=', num2str(nz));
       fprintf(fid, string1);
       fprintf(fid, '\n');
         for ii=1:1:nz;
            for jj=1:1:ny;           
                    fprintf(fid, '%10.4f',double(lro2(:,jj,ii)));  
                    fprintf(fid, '\n');
            end
         end
        fclose(fid);
        
         % record the lro3 profile
       name = [bai(round((t-mod(t,100))/100)+1) shi(round((mod(t,100)-mod(t,10))/10)+1) ge(mod(t,10)+1) nameout3];% 
       nameb = [bai(round((t-mod(t,100))/100)+1) shi(round((mod(t,100)-mod(t,10))/10)+1) ge(mod(t,10)+1) namein3];% 
       fid = fopen(name,'w');
       fprintf(fid,'TITLE = "Phase field variable for lro3" \n');
       fprintf(fid,strcat('VARIABLES=','"',nameb,'"'));
       fprintf(fid, '\n');
       string1 = strcat('zone i=',num2str(nx),',', 'j=',num2str(ny),',', 'k=', num2str(nz));
       fprintf(fid, string1);
       fprintf(fid, '\n');
        for ii=1:1:nz;
            for jj=1:1:ny;           
                    fprintf(fid, '%10.4f',double(lro3(:,jj,ii)));  
                    fprintf(fid, '\n');
            end
         end
        fclose(fid);
        
e_field_x = direction_x*t*vect;
e_field_y = direction_y*t*vect;
e_field_z = direction_z*t*vect;

D(t+1) = mean(mean(mean(sqrt(lro1.^2 + lro2.^2 + lro3.^2)))) + 8.85e-12*e_field*t;

dielectric = 1/8.85e-12*(D(t+1) - D(t))/e_field;
end
    
    
    
    
    
    
  

