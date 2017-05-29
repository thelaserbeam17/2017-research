%% Time Domain Reflectometry With Water Trees
clear;
clc;

%% Local Water Tree Model
e_0 = 8.85e-12;
f = 60;                                     % Hz
w = 2*pi*f;                                 % rad/s
cond_water = 1e-7;                          % conductivity of water
e_water = 81-i*cond_water/(2*pi*f*e_0);     % complex permiativity of water
e_xlpe = 2.3-i*.001;                        % complex permiativity of XLPE

% water tree
kw = .75;                                   % water tree concentration
hw = .4;                                    % moisture content
q_w = kw*hw;                                % water content in water tree
D = 1/4;                                    % geometry of water tree
e_wt = e_xlpe*(1+q_w*(e_water-e_xlpe)/(e_xlpe+D*(1-q_w)*(e_water-e_xlpe))); % water tree permeativity

w = .2;                                     % length parameter of water treed region
e_total = w*e_wt + (1-w)*e_xlpe;            % total permiativity of water treed region

%% Single Core Transmission Line Model
cond_cu = 5.96e7;                           % conductivity of copper
u_0 = (4*pi)*1e-7 ;                         % permeability of free space
r_in = 5e-3;                                   % radius of conductor
r_out = 100e-3;                                  % radius of insulation
C_0 = 2*pi*e_0/log(r_out/r_in);             % geometric capcitance

% Degraded Cable Section
C_deg_cmplx = 2*pi*e_0*e_total/log(r_out/r_in);   % complex capacitance
C_ins_deg = real(C_deg_cmplx);                     % capacitance of insulation for single core cable (F/m)
L_deg = u_0*e_0/C_0;                                 % Inductance of conductor (H/m)
R_deg = 1/(pi*r_in^2*cond_cu);                  % Resistance of conductor (ohms/m)
G_deg = real(2*pi*f*C_ins_deg);                     % Conductance of insulation (S/m)

% Good condition Section
C_cmplx = 2*pi*e_0*e_xlpe/log(r_out/r_in);   % complex capacitance
C_ins = real(C_cmplx);                     % capacitance of insulation for single core cable (F/m)
L = u_0*e_0/C_0;                                 % Inductance of conductor (H/m)
R = 1/(pi*r_in^2*cond_cu);                  % Resistance of conductor (ohms/m)
G = real(2*pi*f*C_ins);                     % Conductance of insulation (S/m)

% Model Length 
length = 100;                               % 100 meter cable
deg_coef = .2;                              % 20 meters of cable is degraded due to water tree

Z_good = (1-.2)*((R + i*w*L)/(G + i*w*C_ins))^.5;   % good condition cable chracteristic impedance (ohms)
Z_deg = .2*((R + i*w*L)/(G + i*w*C_ins))^.5;        % degraded cables characteristic impedance     (ohms)

Z_L = 50;                                           % terminating load (ohms)

