
%复现文章结果
% Weak Measurement of Elliptical Dipole Moments by C-Point Splitting
% C:\Users\Administrator\Desktop\程序示例及修改\斯托克斯参量\iscat

clear;
clc;
nm = 1e-9;
lamda = 633*nm;
%高度要注意调节，越高则倏逝波越小
d = 40*nm;
k0 = 2*pi/lamda;
n = 1;
NA = n*sin(pi/2);

kx = linspace(-n*k0,n*k0,2001);
ky = linspace(-n*k0,n*k0,2001);
[KX,KY] = meshgrid(kx,ky);
K0 = ones(2001)*k0;
KZ = sqrt(K0.^2-KX.^2-KY.^2);
KZ2 = sqrt(n^2*K0.^2-KX.^2-KY.^2);
C = exp(1j*KZ*d).*sqrt(n^2*K0.^2-KX.^2-KY.^2)./KZ;
TP = 2*n*KZ./(KZ2 + n^2*KZ);
TS = 2*KZ./(KZ + KZ2);
p = [0,1,1j]; 
%p = [1,0,0.2i];
    
Efp = (p(1)*C.*KX.*KZ)./(sqrt(KX.^2+KY.^2).*K0) + (p(2)*C.*KY.*KZ)./(sqrt(KX.^2+KY.^2).*K0) - p(3)*C.*sqrt(KX.^2+KY.^2)./K0;
Efs = (-1*p(1)*C.*KY./sqrt(KX.^2+KY.^2)) + p(2)*C.*KX./sqrt(KX.^2+KY.^2);

I = abs(Efs).^2 + abs(Efp).^2;

Elp = Efs + 1j*Efp;
Erp = Efs - 1j*Efp;

subplot(2,2,1)
imagesc(kx/k0,ky/k0,I);title('I')
colormap("jet")
colorbar
axis xy

subplot(2,2,2)
imagesc(kx/k0,ky/k0,abs(Erp).^2);title('Irp')
colorbar
axis xy

subplot(2,2,3)
imagesc(kx/k0,ky/k0,abs(Elp).^2);title('Ilp')
colorbar
axis xy

subplot(2,2,4)
imagesc(kx/k0,ky/k0,abs(Elp).^2+abs(Erp).^2);
colorbar
axis xy