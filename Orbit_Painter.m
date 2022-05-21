R0 = 6371; %radius of the Earth
AA = readmatrix('C:\Qt\ModellingaSpacecraft\build-ModellingaSpacecraft-Desktop_Qt_6_2_4_MinGW_64_bit-Debug\Orbit_Data.txt');
BB = readmatrix('C:\Qt\ModellingaSpacecraft\build-ModellingaSpacecraft-Desktop_Qt_6_2_4_MinGW_64_bit-Debug\Orbit_Data_new.txt');
X = AA(1:end, 1);
Y = AA(1:end, 2);
Z = AA(1:end, 3);
Vx = AA(1:end, 4);
Vy = AA(1:end, 5);
Vz = AA(1:end, 6);
x0 = X(1);
y0 = Y(1);
z0 = Z(1);
X2 = BB(1:end, 1);
Y2 = BB(1:end, 2);
Z2 = BB(1:end, 3);
Vx2 = BB(1:end, 4);
Vy2 = BB(1:end, 5);
Vz2 = BB(1:end, 6);


%modelling the Earth
[x_e, y_e, z_e] = sphere;
x_e = x_e*R0;
y_e = y_e*R0;
z_e = z_e*R0;

%Drawing the orbit of the SpaceCraft
figure
%plot3(0, 0, 0, 'o');
surfl(x_e, y_e, z_e);
xlabel('X, km');
ylabel('Y, km');
zlabel('Z, km');
hold on
plot3(X, Y, Z, 'r');
plot3(x0, y0, z0, '*g');
plot3(X2, Y2, Z2, '--r');
legend('Earth', 'Orbit', 'Spacecraft', 'New orbit');
