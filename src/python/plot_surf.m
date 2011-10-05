clear all; close all; clc;

cdim = 153;
sdim = 81;

[bs] = load('blade_surf.dat');
[bsm] = load('blade_surf_mod.dat');

x = bs(:,1);
y = bs(:,2);
z = bs(:,3);

xm = bsm(:,1);
ym = bsm(:,2);
zm = bsm(:,3);

x = reshape(x,sdim,cdim);
y = reshape(y,sdim,cdim);
z = reshape(z,sdim,cdim);

xm = reshape(xm,sdim,cdim);
ym = reshape(ym,sdim,cdim);
zm = reshape(zm,sdim,cdim);

% d = sqrt((xm-x).^2 + (ym-y).^2 + (zm-z).^2);

% surf(x,y,z,'edgecolor','none')
% colorbar

surf(xm, ym, zm)

return

f = load('normal_field.dat');
x = f(:,1);
y = f(:,2);
z = f(:,3);
p = f(:,4);

x = reshape(x,sdim,cdim);
y = reshape(y,sdim,cdim);
z = reshape(z,sdim,cdim);
p = reshape(p,sdim,cdim);

surf(x(:,1:67),y(:,1:67),z(:,1:67),p(:,1:67),'edgecolor','none')
shading interp
view(0,-90)
title('Suction Surface')
caxis([min(min(p)), max(max(p))])
colorbar
axis equal

figure;
surf(x(:,67:end),y(:,67:end),z(:,67:end),p(:,67:end),'edgecolor','none')
shading interp
view(0,-90)
title('Pressure Surface')
caxis([min(min(p)), max(max(p))])
colorbar
axis equal
