clear all
clc

% MATLAB Code for ECE830 Final Project - Robert Mulla
% Feedback System for 3-Generator Power System

%State Matricies

%Generator Voltages
V1 = 1.0566;
V2 = 1.0502;
V3 = 1.0170;

%Detla Values
e1 = degtorad(2.272);
e2 = degtorad(19.732);
e3 = degtorad(13.175);

%Inertia Constants H

H1 = 35.46;
H2 = 9.6;
H3 = 4.515;

%Damping Coefficient
D1 = 2.36;
D2 = 0.64;
D3 = 0.30;

%Y Bus Matrix can be converted into real (G) and reactive parts (B)
G = [0.846 0.287 0.210; 0.287 0.420 0.213; 0.210 0.213 0.277];
B = [2.988 1.513 1.226; 1.513 2.724 1.088; 1.226 1.088 2.368];

%syncronmous speeds Ws (give) = 377rad/sec

ws = 377;

% Partial derivatives used to create A matrix

df1e1 = 0;
df1e2 = 0;
df1e3 = 0;
df1w1 = 1;
df1w2 = 0;
df1w3 = 0;

df2e1 = 0;
df2e2 = 0;
df2e3 = 0;
df2w1 = 0;
df2w2 = 1;
df2w3 = 0;

df3e1 = 0;
df3e2 = 0;
df3e3 = 0;
df3w1 = 0;
df3w2 = 0;
df3w3 = 1;

df4e1 = (-ws/(2*H1))*(V1*(V2*(-G(1,2)*sin(e1-e2)+B(1,2)*cos(e1-e2))+V3*(-G(1,3)*sin(e1-e3)+B(1,3)*cos(e1-e3))));
df4e2 = (-ws/(2*H1))*(V1*(V2*(G(1,2)*sin(e1-e2)-B(1,2)*cos(e1-e2))));
df4e3 = (-ws/(2*H1))*(V1*(V3*(G(1,3)*sin(e1-e3)-B(1,3)*cos(e1-e3))));
df4w1 = (ws/(2*H1))*(-D1/ws);
df4w2 = 0;
df4w3 = 0;

df5e1 = (-ws/(2*H2))*(V2*(V1*(G(2,1)*sin(e2-e1)-B(2,1)*cos(e2-e1))));
df5e2 = (-ws/(2*H2))*(V2*(V1*(-G(2,1)*sin(e2-e1)+B(2,1)*cos(e2-e1))+V3*(-G(2,3)*sin(e2-e3)+B(2,3)*cos(e2-e3))));
df5e3 = (-ws/(2*H2))*(V2*(V3*(G(2,3)*sin(e2-e3)-B(2,3)*cos(e2-e3))));
df5w1 = 0;
df5w2 = (ws/(2*H2))*(-D2/ws);
df5w3 = 0;

df6e1 = (-ws/(2*H3))*(V3*(V1*(G(3,1)*sin(e3-e1)-B(3,1)*cos(e3-e1))));
df6e2 = (-ws/(2*H3))*(V3*(V2*(G(3,2)*sin(e3-e2)-B(3,2)*cos(e3-e2))));
df6e3 = (-ws/(2*H3))*(V3*(V2*(-G(3,2)*sin(e3-e2)+B(3,2)*cos(e3-e2))+V1*(-G(3,1)*sin(e3-e1)+B(3,1)*cos(e3-e1))));
df6w1 = 0;
df6w2 = 0;
df6w3 = (ws/(2*H3))*(-D3/ws);

% Calculating the A Matrix

A = [df1e1 df1e2 df1e3 df1w1 df1w2 df1w3;
    df2e1 df2e2 df2e3 df2w1 df2w2 df2w3;
    df3e1 df3e2 df3e3 df3w1 df3w2 df3w3;
    df4e1 df4e2 df4e3 df4w1 df4w2 df4w3;
    df5e1 df5e2 df5e3 df5w1 df5w2 df5w3;
    df6e1 df6e2 df6e3 df6w1 df6w2 df6w3]

df4pm1 = ws/(2*H1);
df5pm2 = ws/(2*H2);
df6pm3 = ws/(2*H3);

% B Matrix

B = [0 0 0; 0 0 0; 0 0 0; df4pm1 0 0; 0 df5pm2 0; 0 0 df6pm3]

% C Matrix

C = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0]

% D Matrix is just zero

D = [0]

% Find the eigenvalues of the A matrix

eig(A)

% Checking controlability by crating controlability matrix P

P = [B A*B A^2*B A^3*B A^4*B A^5*B];

% Check rank of controlability matrix, should be 6 to be controlable
fprintf('Rank of the controlability matrix');
fprintf('');
rank(P)

% Check observability by creating observability matrix Q

Q = [C' A'*C' (A^2)'*C' (A^3)'*C' (A^4)'*C' (A^5)'*C'];

% Check rank of Q
fprintf('Rank of the observability matrix');
fprintf('');
rank(Q)

% For output feedback first calculate the desired eigenvalues

Desiredeye = eig(A)-[0.7 ; 0.7 ; 0.7 ; 0.7 ; 0.7 ; 0.7]

% Calculate the phi and psi matrix for a closed loop

phiev1 = inv(Desiredeye(1) * eye(6) - A);
phiev2 = inv(Desiredeye(2) * eye(6) - A);
phiev3 = inv(Desiredeye(3) * eye(6) - A);
phiev4 = inv(Desiredeye(4) * eye(6) - A);
phiev5 = inv(Desiredeye(5) * eye(6) - A);
phiev6 = inv(Desiredeye(6) * eye(6) - A);

psiev1 = C*phiev1*B;
psiev2 = C*phiev2*B;
psiev3 = C*phiev3*B;
psiev4 = C*phiev4*B;
psiev5 = C*phiev5*B;
psiev6 = C*phiev6*B;

% calculate G' E' and K' matrix for output feedback

Gprime = [psiev1(:,1) psiev2(:,2) psiev3(:,3)]

Eprime = [1 0 0; 0 1 0; 0 0 1]

Kprime = -Eprime*inv(Gprime)

% Checking closed loop eigenvalues

Anil = A-B*Kprime*C

eig(Anil)

% State Feedback Calculations--------------------------
% State feedback psi values

Spsiev1 = phiev1*B;
Spsiev2 = phiev2*B;
Spsiev3 = phiev3*B;
Spsiev4 = phiev4*B;
Spsiev5 = phiev5*B;
Spsiev6 = phiev6*B;


G = [Spsiev1(:,1) Spsiev2(:,2) Spsiev3(:,3) Spsiev4(:,1) Spsiev5(:,2) Spsiev6(:,3)]

E = [1 0 0 1 0 0; 0 1 0 0 1 0; 0 0 1 0 0 1]

K = -E*inv(G)

Ao = A-B*K

fprintf('State Feedback eigenvalues:');
fprintf('');
eig(Ao)

% Check State equation for controlability and observability

P = [B Ao*B Ao^2*B Ao^3*B Ao^4*B Ao^5*B];

% Check rank of controlability matrix, should be 6 to be controlable
fprintf('rank of controlability matrix for state feedback');
fprintf('');

rank(P)

% Check observability by creating observability matrix Q

Q = [C' Ao'*C' (Ao^2)'*C' (Ao^3)'*C' (Ao^4)'*C' (Ao^5)'*C'];

% Check rank of Q
fprintf('rank of observability matrix for state feedback');
fprintf('');


rank(Q)

% Observer ---------------------

% The desired eigenvalues for the observer
fprintf('desired observer eigenvalues');
fprintf('');
DesiredeyeObs = eig(A)-[0.7 ; 0.7 ; 0.7 ; 0.7 ; 0.7 ; 0.7]*10

% Calculate the phi and Cphi matrix for the observer

phiev1Obs = inv(DesiredeyeObs(1) * eye(6) - A);
phiev2Obs = inv(DesiredeyeObs(2) * eye(6) - A);
phiev3Obs = inv(DesiredeyeObs(3) * eye(6) - A);
phiev4Obs = inv(DesiredeyeObs(4) * eye(6) - A);
phiev5Obs = inv(DesiredeyeObs(5) * eye(6) - A);
phiev6Obs = inv(DesiredeyeObs(6) * eye(6) - A);

Cphiev1 = C*phiev1Obs;
Cphiev2 = C*phiev2Obs;
Cphiev3 = C*phiev3Obs;
Cphiev4 = C*phiev4Obs;
Cphiev5 = C*phiev5Obs;
Cphiev6 = C*phiev6Obs;

Gc = [Cphiev1(1,:) ; Cphiev2(2,:) ; Cphiev3(3,:) ; Cphiev4(1,:) ; Cphiev5(2,:) ; Cphiev6(3,:)];

Ec = [1 0 0 ; 0 1 0; 0 0 1 ; 1 0 0; 0 1 0 ; 0 0 1];

% observer matrix

L = -1*inv(Gc)*Ec

% Check observer eigenvalues
fprintf('checked observer eigenvalues');

Ac=A-L*C

eig(A-L*C)

% Check observer for controlability and observability

Pc = [B Ac*B Ac^2*B Ac^3*B Ac^4*B Ac^5*B];

% Check rank of controlability matrix, should be 6 to be controlable
fprintf('rank of controlability matrix for observer');
fprintf('');

rank(Pc)

% Check observability of observer by creating observability matrix Q

Qc = [C' Ac'*C' (Ac^2)'*C' (Ac^3)'*C' (Ac^4)'*C' (Ac^5)'*C'];

% Check rank of Q
fprintf('rank of observability matrix for observer');
fprintf('');

rank(Qc)

% Reducing the controllers used:

% New B matrix

Breduced = [0 ; 0 ; 0 ; df4pm1 ; 0 ; 0 ]

Spsiev1 = phiev1*Breduced;
Spsiev2 = phiev2*Breduced;
Spsiev3 = phiev3*Breduced;
Spsiev4 = phiev4*Breduced;
Spsiev5 = phiev5*Breduced;
Spsiev6 = phiev6*Breduced;


Greduced = [Spsiev1(:,1) Spsiev2(:,1) Spsiev3(:,1) Spsiev4(:,1) Spsiev5(:,1) Spsiev6(:,1)]

Ereduced = [1 1 1 1 1 1]

Kreduced = -Ereduced*inv(Greduced)

Aoreduced = A-Breduced*Kreduced

fprintf('State Feedback eigenvalues:');
fprintf('');
eig(Ao)

% reduced outputs for the observer

Creduced = [1 0 0 0 0 0];

Cphiev1 = Creduced*phiev1Obs;
Cphiev2 = Creduced*phiev2Obs;
Cphiev3 = Creduced*phiev3Obs;
Cphiev4 = Creduced*phiev4Obs;
Cphiev5 = Creduced*phiev5Obs;
Cphiev6 = Creduced*phiev6Obs;

Gcred = [Cphiev1(1,:) ; Cphiev2(1,:) ; Cphiev3(1,:) ; Cphiev4(1,:) ; Cphiev5(1,:) ; Cphiev6(1,:)];

Ecred = [1 ; 1; 1 ; 1 ; 1 ; 1];

% observer matrix

Lred = -1*inv(Gcred)*Ecred

% Check reduced observer eigenvalues
fprintf('checked reduced observer eigenvalues');

AcReduced=A-Lred*Creduced

eig(AcReduced)
