%% Time Domain Reflectometry With Water Trees
clear;
clc;

%% Local Water Tree Model
e_0 = 8.85e-12;
c = 3e8;                                    % speed of light m/s
f = 60;                                     % Hz
w = 2*pi*f;                             % rad/s
cond_water = 1e-7;                          % conductivity of water
e_water = 81-i*cond_water/(2*pi*f*e_0);     % complex permiativity of water
e_xlpe = 2.3-i*.001;                        % complex permiativity of XLPE

% water tree
kw = .75;                                   % water tree concentration
hw = .4;                                    % moisture content
q_w = kw*hw;                                % water content in water tree
D = 1/4;                                    % geometry of water tree
e_wt = e_xlpe*(1+q_w*(e_water-e_xlpe)/(e_xlpe+D*(1-q_w)*(e_water-e_xlpe))); % water tree permeativity

x = .22;                                     % length parameter of water treed region
%e_total = w*e_wt + (1-w)*e_xlpe;            % total permiativity of water treed region
e_total = 1/(x/e_wt + (1-x)/e_xlpe);

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

%% Model  
length = .01;                               % 1 centimeter PUL set up
deg_coef = .9;                              % Governs what percent of cable is degrade

Z_good = ((R + i*w*L)/(G + i*w*C_ins))^.5;   % good condition cable chracteristic impedance (ohms)
Z_deg = ((R_deg + i*w*L_deg)/(G_deg + i*w*C_ins_deg))^.5;        % degraded cables characteristic impedance     (ohms)

Z_L = 1e100;                                           % terminating load (ohms)

v_p = 1/(L*C_ins)^.5;                           % propagation velocity of good condition
v_p_deg = 1/(L_deg*C_ins_deg)^.5;               % propagation velocity of water treed region

gamma_good = ((R + i*w*L)*(G + i*w*C_ins))^.5;                  % proagation constant good
gamma_deg = ((R_deg + i*w*L_deg)*(G_deg + i*w*C_ins_deg))^.5;   % propagation constant wter treed

% Good condition
a11 = cosh(gamma_good*length); 
a12 = Z_good*sinh(gamma_good*length);
a21 = 1/Z_good*sinh(gamma_good*length);
a22 = cosh(gamma_good*length);

T1 = [a11 a12;                                          % ABCD parameters for good section of cable
      a21 a22];   

% Degraded
a11_d = cosh(gamma_deg*length); 
a12_d = Z_deg*sinh(gamma_deg*length);
a21_d = 1/Z_deg*sinh(gamma_deg*length);
a22_d = cosh(gamma_deg*length);

T2 = [a11_d a12_d;
      a21_d a22_d];  

% Overall transmission Line properties
N = 10000;           % sample size for PUL sections - will yeild 100 m line for .01 m samples
x = zeros(1, N);
 
for k = 1:N
    if(k == 1)
        T = T1*T1;
    end
    if(k < (1-deg_coef)*N)
        T = T*T1;
    end
    if(k >= (1-deg_coef)*N)
        T = T*T2;
        x(k) = 1;
    end           
end

Z_trans = (T(1,2)/T(2,1))^.5;                                 % overall characteristic impedance of line
v_trans = length * 1/(imag(T(1,2))*imag(T(2,1))/w^2)^.5 * N;      % overall propagation velocity of line
gamma_trans = 1/length * (T(1,2)*T(2,1))^.5 / N;                  % overall propagation constant of line


% Trasnmission Line (Sensor at begining)
trans_length = N*length;                                        % length of line in meters

ref_1 = (Z_deg - Z_good)/(Z_deg + Z_good);                        % Reflection coefitient at first boandary
ref_end = (Z_L - Z_good)/(Z_L + Z_good);                          % Reflection coefitient at end of cable

t_p = trans_length/v_trans;                                     % time of pusle propegation along one lenght of cable 

t_ref_1 = trans_length*(1-deg_coef)/v_trans;                    % time to reflection

time = linspace(0, t_p, N);                                       % time vector
V_p = zeros(1, N);
V_in = zeros(1, N);
V_tr = zeros(1, N);
V_ref = zeros(1, N);

noise = wgn(N, 1, 0);

% input pulse
for k = 1:N
    if(time(k) <= 5e-9)
        V_in(k) = 240 + noise(k);
    end
    if(time(k) > 5e-9)
        V_in(k) = 0 + noise(k);
    end 
end
figure;
plot(time, V_in);

distance = linspace(0, trans_length, N);
V_p(1) = V_in(1);
ref_occured = 0;

for k = 1:N
    if(x(k) == 1 && ref_occured ~= 2)
        ref_occured = 1;
    end
    if(x(k) ~= 1)
        ref_occured = 0;
    end
    
    if(ref_occured ~= 1 && k ~= N)
        V_p(k+1) = V_p(k)*exp(-gamma_trans*distance(k));
        V_ref(k) = noise(k);
    end
    if(ref_occured == 1)
        V_ref(k) = ref_1*V_p(k) + noise(k);
        V_p(k+1) = V_p(k)*(1+ref_1);
        ref_occured = 2;
    end
    if(k == N)
        V_ref(k) = ref_end*V_p(k) + noise(k);
    end
end

figure;
plot(v_trans*time, V_ref);
figure;
plot(v_trans*time, V_p);

