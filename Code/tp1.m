%RAFFARA-BRUNENGO : TP CALIBRAGE

clear all;
close all;
clc;


point_mire = dlmread('ImageMire.txt','',4,1);

point_rep_objet = point_mire(:,1:3);

point_rep_image = point_mire(:,4:5);
coord_point_princ = 256*ones(10,2);

point_rep_image_tilda = point_rep_image - coord_point_princ;

A = ones(10,7);
for k=1:10
    for i=1:3
        A(k,i) = (point_rep_image_tilda(k,2))*point_rep_objet(k,i);
    end
    i=1;
    for j=5:7
        A(k,j) = -1 * (point_rep_image_tilda(k,1))*point_rep_objet(k,i);
        i=i+1;
    end
    A(k,4) = point_rep_image_tilda(k,2);
end

%On peut afficher la matrice A ici
%disp(A);

G = transpose(A) * A;
V = transpose(A) * point_rep_image_tilda(:,1);
L = inv(G) * V;

o2_c = 1/sqrt(L(5)^2+L(6)^2+L(7)^2);
beta = o2_c * sqrt(L(1)^2+L(2)^2+L(3)^2);
o1_c = (L(4) * o2_c) / beta;
r11 = (L(1) * o2_c) / beta;
r12 = (L(2) * o2_c) / beta;
r13 = (L(3) * o2_c) / beta;
r21 = (L(5) * o2_c) / beta;
r22 = (L(6) * o2_c) / beta;
r23 = (L(7) * o2_c) / beta;
r1=[r11,r12,r13];
r2=[r21,r22,r23];
r3=cross(r1,r2);

phy = -atan(r2(3)/r3(3));
gamma = - atan(r12/r11);
omega = atan(r13/(-r23*sin(phy)+r3(3)*cos(phy)));

%On a ici tous les parametres extrinseques

B = ones(10,2);
for k=1:10
    B(k,1) = point_rep_image_tilda(k,2);
    B(k,2) = -1*(r21*point_rep_objet(k,1) +r22*point_rep_objet(k,2) +r23*point_rep_objet(k,3) +o2_c);
end

R = ones(10,1);
for k=1:10
    R(k,1) = -1 * point_rep_image_tilda(k,2) * (r3(1)*point_rep_objet(k,1)+ r3(2)*point_rep_objet(k,2) +r3(3)*point_rep_objet(k,3));
end

G1 = transpose(B)*B;
V1 = transpose(B)*R;
L1 = inv(G1)*V1;

o3_c = L1(1);
f2 = L1(2);

%Commentaire emploi de 10 billes: Suffisant pour avoir un résultat exacte,
%6 étant le minimum.
%L'incertitude résdide plus sur le centre des billes 
%Il faudrait confirmer l'emploi des billes sur des échantillons bruités

