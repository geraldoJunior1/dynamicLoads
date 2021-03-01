%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCRIPT PARA A DETERMINAÇÃO DE CARGAS DINÂMICAS  %
%      OCASIONADAS POR MANOBRAS E RAJADAS         %
%                versao beta-1.0                  %
%                                                 %
%       Desenvolvido por: Gerlado Junior          %
%               Harpia Aerodesign                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; 
clear all;

%caso de voo
velocidade = 150; %velocidade do aviao
rho0 = 1.225; %densidade
rsigma = 0.8; %sweep medio da asa
V = velocidade / rsigma; %velocidade aparente
gambiarra = 2;
%setagem da análise
t0 = 0; %tempo inicial
tfin = 5.0;  %tempo final
dt = 0.005; %incremento

%n mexer aqui
t = [t0: dt: tfin]'; 

[N, dummy] = size(t); 
[Nsim, dummy] = size(t);

%profundor
t_deflex = 1.0; %intervalo de tempo de ir de 0 ate o angulo de deflexao
ang_deflex = 2.0; %angulo de deflexao

%rajada
vel_rajada = 5.0; % velocidade da rajada
tamanho_rajada = 250.0; % comprimento da rajada

% fase de voo, normalmalmente usa o que ta aqui mesmo
n = 1.0; 
q_pr = 0;

%frequencia do modo de flexão da fuselagem
f_e = 4.0; 

% coeficientes aerodinâmicos
aw = 4.5; %coeficiente aw
at = 3.2;  %coeficiente aw do tail
ae = 1.5;  %coeficiete aw do profundor
alpha_0 = -0.03; %incidencia
cm0 = - 0.03; %coeficiente de momento com alpha 0
cd = 0.1; %coeficiente de arrasto
k_epsilon = 0.35; %d_epslon/d_alpha, normalmente entre 0.35-0.4 (downwash)


%parâmetros do avião
m = 10000; %TOW

area_asa = 30.0; %area de asa
area_tail = 7.5; %area do tail
enverga = 7.5; %envergadura
corda = 2.0; %corda

W = m * 9.81;
mf = 0.15 * m;
mw = 0.3 * m; 
mc = 0.4 *m; 
mt = 0.15 * m;

%estimativas, mas se tiver esses valores, é bom colocar
lw = 0.3*corda;
lt = 3.5 * corda;
la = 0.125 * corda; 
le = 0.125 * corda;

%aqui n mexe
lmw = lw - la - le; 
lf = (mt * lt- mw * lmw) / mf;
lwt = lw + lt; 
%ln = lf; 
%lm = 0.375 * corda; 
%mu = mw / 2 / enverga;

% Momentos de inercia
Iyfuse = mf * lf^2 + mt * lt^2; 
Iyasa = mw * (corda / 3)^2; 
Iy = Iyfuse + Iyasa + mw * lmw^2; 
%lyasa = sqrt(Iyasa / mw);
%ly = sqrt(Iy / m);

%profundor set
eta = zeros(N, 1);
npulse = t_deflex / dt + 1;
ang_deflex = ang_deflex*pi / 180;
eta(1:npulse) = ang_deflex*ones(npulse, 1);

%
xg = V*t; 
dx = V*dt; 
Nb = round(lwt / dx+1); 
N = Nsim + Nb - 1;

%rajadinhas
Nd = round(tamanho_rajada / dx + 1); 
Nbd = Nb + Nd - 1; 
wg = zeros(N, 1);
wgw = zeros(N, 1);
wgt = zeros(N, 1);
wg(Nb:Nbd) = (vel_rajada / 2)*(1 - cos(2*pi*xg(1:Nd) / tamanho_rajada)); % vetor de velocidade no centro de massa
wgw(1:Nsim) = wg(Nb:N);
wgt(1:Nsim) = wg(1:Nsim); % vetor de velocidade na asa e no tail

%flexão de fus
kappa_e0 = 1;
gamma_e0 = 0; 
A = 0;
B = 0;

% norm1
kappa_tip = kappa_e0*(1 + A) + gamma_e0*(1 + B)*(corda - corda/4 - la); 
kappa_e0 =kappa_e0/kappa_tip;

% flexao de fus solucao
X = [mf mt; -mf*lf mt*lt]; 
Y = [- (mw + mc)*kappa_e0;

mw*lmw*kappa_e0];
Z = X\Y; 
kappa_ec = 1; 
kappa_ef = Z(1);
kappa_et = Z(2); 
gamma_et = 2*(kappa_et - kappa_ec) / lt - gamma_e0;

% baguncinha da fuselagem

me = mf*kappa_ef^2 + mw*kappa_e0^2 + mc*kappa_ec^2 + mt* kappa_et^2;
omega_e = 2*pi*f_e; 

ke = omega_e^2*me;

% derivadas aerodinamicas
J1 = gamma_e0*(1 + B / 2);
J2 = kappa_e0*(1 + A / 3) - la*gamma_e0*(1 + B / 2);
J3 = gamma_e0*kappa_e0*(1 + A / 3 + B / 2 + A*B / 4) - la* gamma_e0^2*(1 + B + B^2/ 3);

% derivadas aerodinamicas no eixo inercial
z0 = -0.5*rho0*velocidade^2*[- area_asa*aw + area_tail*at*k_epsilon]*alpha_0;
zalpha = -0.5*rho0*velocidade^2*[area_asa*aw + area_tail*at*(1 - k_epsilon)];
zq = -0.5*rho0*velocidade*area_tail*at*lt;
zeta = -0.5*rho0*velocidade^2*area_tail*ae;
dot_Zz = -0.5*rho0*velocidade*(area_asa*aw + area_tail*at*(1 - k_epsilon));
zgw = -0.5*rho0*velocidade*area_asa*aw;
zgt = -0.5*rho0*velocidade*area_tail*at*(1 - k_epsilon);

% derivadas aerodinamicas no eixo inercial em pitch
m0w = 0.5*rho0*velocidade^2*area_asa*corda*cm0 -0.5*rho0*velocidade^2*area_asa*aw*lw*alpha_0; 
m0t= -0.5*rho0*velocidade^2*area_tail*at*k_epsilon*lt*alpha_0;
m0 = m0w + m0t;
malpha = 0.5*rho0*velocidade^2*[area_asa*aw*lw - area_tail*at*(1 - k_epsilon)*lt]; 
mq =-0.5*rho0*velocidade*area_tail*at*lt^2;
meta = -0.5*rho0*velocidade^2*area_tail*ae*lt;
dor_Mz = 0.5*rho0*velocidade*(area_asa*aw*lw - area_tail*at*lt*(1 - k_epsilon));
mgw = 0.5*rho0*velocidade*area_asa*aw*lw;
mgt = -0.5*rho0*velocidade*area_tail*at*lt*(1 - k_epsilon);

% derivadas aerodinamicas no eixo inercial em rajada
zw = -0.5*rho0*velocidade*(area_asa*aw + area_tail*at*(1 - k_epsilon)+ area_asa*cd);
mw = 0.5*rho0*velocidade*(area_asa*aw*lw - area_tail*at*lt*(1 - k_epsilon));

% derivadas aerodinamicas considerando como flexivel
ze = 0.5*rho0*velocidade^2*(-area_asa*aw*J1 - area_tail*at*gamma_et);
dot_ze = -0.5*rho0*velocidade*area_tail*at*kappa_et;
me = 0.5*rho0*velocidade^2*(area_asa*aw*lw*J1 - area_tail*at*lt*gamma_et); 
dot_me = -0.5*rho0*velocidade*area_tail*at*lt*kappa_et;
q0 = 0.5*rho0*velocidade^2*(area_asa*aw*J2 - area_tail*at*k_epsilon*kappa_et) *alpha_0;
qalpha = 0.5*rho0*velocidade^2*(-area_asa*aw*J2 - area_tail*at*(1 - k_epsilon) *kappa_et);
qq = -0.5*rho0*velocidade*area_tail*at*lt*kappa_et;
qeta = -0.5*rho0*velocidade^2*area_tail*ae*kappa_et;
qe = 0.5*rho0*velocidade^2*(-area_asa*aw*J3 - area_tail*at*gamma_et*kappa_et); 
dot_qz =0.5*rho0*velocidade*(-area_asa*aw*J2 - area_tail*at*(1 - k_epsilon) *kappa_et);
dot_qe = -0.5*rho0*velocidade*area_tail*at*kappa_et^2; 
qgw = -0.5*rho0*velocidade*area_asa*aw*J2;
qgt = -0.5*rho0*velocidade*area_tail*at*kappa_et;

% derivadas aerodinamicas considerando como flexivel em rajada
qw = 0.5*rho0*velocidade*(- area_asa*aw*J2 - area_tail*at*(1 - k_epsilon) *kappa_et);


%%%%% CORPO RIGIDO %%%%%%

% Equacao de movimento pra o caso rigido
Aaa = -[zeta zalpha; meta malpha];
Ccc = [1 ; 0]; 
Ddd = [zq; mq]; 
Eee = [z0; m0]; 
Bbb = inv(Aaa) *(Ccc*(n*W) + Ddd*q_pr + Eee); 
Bdeg = Bbb*180 / pi;
prof_trim = Bdeg(1);
inciaio = Bdeg(2);

%%%%% CORPO ELASTICO %%%%%%

% Equacao de movimento pra o caso elastico
Aelastico = - [zeta zalpha ze; meta malpha me; qeta qalpha qe - ke];
Celastico = [1 ; 0 ; 0]; 
Delastico = [zq ; mq ; qq];
Eelastico = [z0 ; m0 ; q0];
Belastico = inv(Aelastico)*(Celastico*(n*W) + Delastico*q_pr + Eelastico); 
Bdeg = Belastico*180 / pi;
prof_trim_elastico = Bdeg(1);
inciaio_elastico = Bdeg(2);

%%%%%%%%%%%%% RESPOSTA COM O PROFUNDOR %%%%%%%%%%%%%

%condicao inicial
U_e = velocidade; 
W_e = 0;

% equacoes de movimento linearizaras na mecanica de voo
M = [m 0;0 Iy]; 
C = [-zw -(m*U_e + zq); -mw -mq];
F = [zeta; meta];
MC = inv(M)*C; 
MF = inv(M)*F;

%Resposta no referencial corpo
input = [t, eta];
sim('resp');

w = rate(:,1); % velocidade downward
q = rate(:,2); % pitch rate
wdot = acc(:,1); % aceleracoes downward
theta = delocs(:,2);
alpha = w / velocidade; % incidencia
gamma = theta - alpha; % pertubacao - path angle
az = wdot - q*U_e; % fator de carga
 
%gambiarra
t2 = [t0: dt/gambiarra: tfin]';
topa=t2;

%plot da resposta em relacao ao corpo
figure;
plot(t2, q*180/pi,'k-', t2, alpha *180/pi, 'k:')
title('resposta da incidencia e Pitch Rate para o input do profundor') 
xlabel('t(s)'); 
ylabel('Pitch Rate (deg/s), incidencia (deg)'); 
legend('Pitch Rate','incidencia')

figure;
plot(t2, az/9.81, 'k-')
title('Fator de carga com o input do profundor');
xlabel ('t(s)');
ylabel('fator de carga(g)')

figure;
plot(t2, theta*180/pi, 'k-', t2, gamma*180/pi, 'k:')
title('Pitch e Flight Path Angles com o input do profundor')
xlabel('t (s)');


%%%%%%%%%%%%% RESPOSTA EM RAJADA %%%%%%%%%%%%%

% Equações de movimento
M = [m 0; 0 Iy]; 
CC = - [dot_Zz zq; dor_Mz mq]; 
KK = - [0 zalpha; 0 malpha];
FW = [zgw; mgw]; 
FT = [zgt; mgt];
MI = inv(M); 
MC = MI*CC; 
MK = MI*KK; 
MFW = MI*FW; 
MFT = MI*FT;

x=size(t)
inputW = [t, wgw(1:Nsim)];
inputT = [t, wgt(1:Nsim)]; 
sim('rajada') 

% deslocamentos, velocidade e acelerações no centro de massa, fus e tail

zc = delocs(:,1); 
dot_zc = rate(:,1); 
ddot_zc = acc(:,1); 
theta = delocs(:,2); 
dot_theta = rate(:,2); 
ddot_theta = acc(:,2);
zf = zc - lf*theta; 
ddot_zf = ddot_zc - lf*ddot_theta; 
zt = zc + lt*theta; 
ddot_zt = ddot_zc + lt*ddot_theta;

%%plots com o eixo no corpo

t = topa;
figure;
plot(t, theta*180/pi, 'k-', t, dot_theta *180/pi, 'k:') 
title('resposta em pitch para a rajada')
xlabel('t(s)'); 
ylabel('Pitch (deg) e Pitch Rate (deg/s)'); 
legend('Angulo de pitch', 'Pitch rate')

figure;
plot(t, ddot_zc/9.81, 'k-', t, ddot_zt/9.81, 'k:') 
title('fator de carga com a rajada (canolli assoprando)') 
xlabel('t(s)'); 
ylabel('fator de carga(g)'); 
legend('Centro de massa', 'Tail')

%centro de massa
U_E = U_e + W_e*theta;
W_E = -U_e*theta + W_e + w;

% tail
U_E_tp = U_e + W_e*theta + w.* theta; 
W_E_tp = -U_e*theta + W_e + w + lt*q;

% Plots em relação ao eixo de cordenadas no solo

figure
plot(t, w, 'k-', t, W_E, 'k:')
title('Velocidade vertical com o input do profundor (referencial vento/solo)')
xlabel('t(s)'); 
ylabel('Velocidade (m/s)') 
legend('velocidade - vento', 'velocidade - solo')

figure;
plot(t, W_E, t, W_E_tp, 'k-')
title('Velocidade vertical com o input do profundor (referencial solo)')
xlabel('t(s)'); 
ylabel('Velocidade(m/s)'); 
legend ('Centro de massa', 'Tail')

sim('convert');

% resposta dop centro de massa no refericial solo
figure;
plot(X_E, -Z_E, 'k-')
title('perfil de voo com o input do profundor (referencial solo)') 
xlabel('deslocamento horizontal(m)'); 
ylabel('deslocamento vertical(m)')

