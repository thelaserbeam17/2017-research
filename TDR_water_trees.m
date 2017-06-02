%% Time Domain Reflectometry With Water Trees
clear;
clc;
n = 100000;
freq = linspace(1e4, 1e8, n);
e = zeros(1, n);
e_t = zeros(1, n);
e_w = zeros(1, n);

%% Local Water Tree Model
e_0 = 1e-9/(36*pi);
c = 3e8;                                    % speed of light m/s
f = 60;                                     % Hz
w = .0001;                             % rad/s
cond_water = .0052;                          % conductivity of water
e_water = 81-i*cond_water/(2*pi*f*e_0);     % complex permiativity of water
e_xlpe = 2.3-i*.001;                        % complex permiativity of XLPE

% water tree
kw = .75;                                   % water tree concentration
hw = .4;                                    % moisture content
q_w = .09;                                % water content in water tree
D = 1/5;                                    % geometry of water tree
e_wt = e_xlpe*(1+q_w*(e_water-e_xlpe)/(e_xlpe+D*(1-q_w)*(e_water-e_xlpe))); % water tree permeativity

y = .25;                                     % length parameter of water treed region
%e_total = x*e_wt + (1-x)*e_xlpe;            % total permiativity of water treed region
e_total = 1/(y/e_wt + (1-y)/e_xlpe);

for j = 1:n
    e(j) = 81-i*cond_water/(2*pi*freq(j)*e_0);
    e_w(j) = e_xlpe*(1+q_w*(e(j)-e_xlpe)/(e_xlpe+D*(1-q_w)*(e(j)-e_xlpe)));
    e_t(j) = 1/(y/e_w(j) + (1-y)/e_xlpe);
end

figure;
subplot(2, 1, 1);
semilogx(freq, real(e_t));
title('Real Permiativity')
ylim([2.3 2.7]);

subplot(2, 1, 2);
semilogx(freq, abs(imag(e_t)));
title('imaginary permiativity')
%% Single Core Transmission Line Model
cond_cu = 5.96e7;                           % conductivity of copper
u_0 = (4*pi)*1e-7 ;                         % permeability of free space
r_in = 10e-3;                                   % radius of conductor
r_out = 30e-3;                                  % radius of insulation
C_0 = 2*pi*e_0/log(r_out/r_in);             % geometric capcitance

% Degraded Cable Section
C_deg_cmplx = 2*pi*e_0*e_total/log(r_out/r_in);   % complex capacitance
C_ins_deg = real(C_deg_cmplx);                     % capacitance of insulation for single core cable (F/m)
L_deg = u_0*e_0/C_0;                                 % Inductance of conductor (H/m)
R_deg = 1/(pi*r_in^2*cond_cu);                  % Resistance of conductor (ohms/m)
G_deg = -imag(2*pi*f*C_ins_deg);                     % Conductance of insulation (S/m)

% Good condition Section
C_cmplx = 2*pi*e_0*e_xlpe/log(r_out/r_in);   % complex capacitance
C_ins = real(C_cmplx);                     % capacitance of insulation for single core cable (F/m)
L = u_0*e_0/C_0;                                 % Inductance of conductor (H/m)
R = 1/(pi*r_in^2*cond_cu);                  % Resistance of conductor (ohms/m)
G = -imag(2*pi*f*C_ins);                     % Conductance of insulation (S/m)

%% Model  
length = .01;                               % 1 centimeter PUL set up
deg_coef = .00005;                              % Governs what percent of cable is degrade
deg_pos = .5;                               % Governs position alon cable degredation is

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
deg_flag = 0;
 
for k = 1:N
    if(k == 1)
        T = T1*T1;
    end
    if(k < deg_pos*N && deg_flag == 0)
        T = T*T1;
    end
    if(k >= deg_pos*N)
        if(k <= (deg_pos + deg_coef)*N)
            T = T*T2;
            x(k) = 1;
            deg_flag = 1;
        end
        if(k > (deg_pos + deg_coef)*N)
            T = T*T1;
            deg_flag = 0;
        end
    end           
end

Z_trans = (T(1,2)/T(2,1))^.5;                                 % overall characteristic impedance of line
v_trans = length * 1/(imag(T(1,2))*imag(T(2,1))/w^2)^.5 * N;      % overall propagation velocity of line
gamma_trans = 1/length * (T(1,2)*T(2,1))^.5 / N;                  % overall propagation constant of line

% Trasnmission Line (Sensor at begining)
trans_length = N*length;                                        % length of line in meters

ref_1 = (Z_deg - Z_good)/(Z_deg + Z_good);                        % Reflection coefitient at start of first boandary
ref_2 = (Z_good - Z_deg)/(Z_deg + Z_good);                        % Reflection coefitient at end of first boandary
ref_end = (Z_L - Z_good)/(Z_L + Z_good);                          % Reflection coefitient at end of cable

t_p = trans_length/v_trans;                                     % time of pusle propegation along one lenght of cable 

t_ref_1 = trans_length*(1-deg_coef)/v_trans;                    % time to reflection

time = linspace(0, t_p, N);                                       % time vector
V_p = zeros(1, N);
V_in = zeros(1, N);
V_tr = zeros(1, N);
V_ref = zeros(1, N);

noise = zeros(1, N); %wgn(N, 1, 0);

% input pulse
for k = 1:N
    if(time(k) <= 5e-9)
        V_in(k) = 240 + noise(k);
    end
    if(time(k) > 5e-9)
        V_in(k) = 0 + noise(k);
    end 
end

distance = linspace(0, trans_length, N);
V_p = V_in;
ref_occured = 0;
% ref_occured chart
%       0         1           2            3
%    no ref   ref nw-yw   ref yw-nw    

for k = 1:N
    if(k > 1)
        if(x(k) == 1 && x(k-1) == 0 && ref_occured ~= 3)
            ref_occured = 1;
        end
        if(x(k) == 0 && x(k-1) == 1 && ref_occured ~= 3)
            ref_occured = 2;
        end
    
        if(x(k) ~= 1 && x(k-1) ~= 1)
            ref_occured = 0;
        end
    end
    if(ref_occured ~= 1 && k ~= N)
        V_p(k+1) = V_p(k)*exp(-gamma_trans*distance(k));
        V_ref(k) = noise(k);
    end
    if(ref_occured == 1)
        V_ref(k) = ref_1*V_p(k) + noise(k);
        V_p(k+1) = V_p(k)*(1+ref_1);
        ref_occured = 0;
    end
    if(ref_occured == 2)
        V_ref(k) = ref_2*V_p(k) + noise(k);
        V_p(k+1) = V_p(k)*(1+ref_2);
        ref_occured = 0;        
    end
    if(k == N)
        V_ref(k) = ref_end*V_p(k) + noise(k);
    end
end

figure;
subplot(2, 2, 1);
plot(time, V_in);
title('Input Pulse');
xlabel('time (s)');

subplot(2, 2, 2);
plot(v_trans*time, V_p);
title('Propagation Amplitude');
xlabel('Length across line (m)');

subplot(2, 2, [3 4]);
plot(v_trans*time, V_ref);
title('Reflection Pattern')
xlabel('Length across line (m)');

%% LIRA Good Condition Cable
n = 10000;
freq = linspace(1, 10e7, n);  
t = linspace(0, 10e-7, n);
ch1 = zeros(1, n);

Z_dut = zeros(1, n);
theta = zeros(1, n);
v_phase = zeros(1, n);

Z_dut_d = zeros(1, n);
theta_d = zeros(1, n);
v_phase_d = zeros(1, n);

l_trans = 100;
y = .25;

for k = 1:n
    % Good condition cable
    C_cmplx = 2*pi*e_0*e_xlpe/log(r_out/r_in);   % complex capacitance
    C_ins = real(C_cmplx);                     % capacitance of insulation for single core cable (F/m)
    L = u_0*e_0/C_0;                                 % Inductance of conductor (H/m)
    R = 1/(pi*r_in^2*cond_cu);                  % Resistance of conductor (ohms/m)
    G = - freq(k) * imag(C_cmplx);                     % Conductance of insulation (S/m)
    
    % degraded cable
    e_water = 81-i*cond_water/(freq(k)*e_0);
    e_wt = e_xlpe*(1+q_w*(e_water-e_xlpe)/(e_xlpe+D*(1-q_w)*(e_water-e_xlpe)));
    e_total = 1/(y/e_wt + (1-y)/e_xlpe);    
    
    C_deg_cmplx = 2*pi*e_0*e_total/log(r_out/r_in);   % complex capacitance
    C_ins_deg = real(C_deg_cmplx);                     % capacitance of insulation for single core cable (F/m)
    L_deg = u_0*e_0/C_0;                                 % Inductance of conductor (H/m)
    R_deg = 1/(pi*r_in^2*cond_cu);                  % Resistance of conductor (ohms/m)
    G_deg = -imag(2*pi*f*C_ins_deg);                     % Conductance of insulation (S/m)
    
    % LIRA math
    % ////////////////////
    % Good coniditon
    Z1 = 50;
    Z0 = ( (R + i*freq(k)*L)/(G + i*freq(k)*C_ins) )^.5;
    ref_L = ( Z1 - Z0 ) / ( Z1 + Z0 );
    g = ( (R + i*freq(k)*L) * (G + i*freq(k)*C_ins) )^.5; 
    ref_d = ref_L * exp( -2 * g * l_trans);
    
    v_phase(k) = freq(k) / imag(g);
    
    Z_dut(k) = Z0 * ( (1 + ref_d) / (1 - ref_d) );
    
    theta(k) = angle( Z_dut(k) );
    
    % degraded
    Z0_deg = ( (R_deg + i*freq(k)*L_deg)/(G_deg + i*freq(k)*C_ins_deg) )^.5;
    ref_L = ( Z1 - Z0_deg ) / ( Z1 + Z0_deg );
    g_deg = ( (R_deg + i*freq(k)*L_deg) * (G_deg + i*freq(k)*C_ins_deg) )^.5; 
    ref_d = ref_L * exp( -2 * g_deg * l_trans);
    
    v_phase_d(k) = freq(k) / imag(g_deg);
    
    Z_dut_d(k) = Z0_deg * ( (1 + ref_d) / (1 - ref_d) );
    
    theta_d(k) = angle( Z_dut_d(k) );    
end
h = ifft(theta);
h_d = ifft(theta_d);
f_vec = imag(fft(theta));
f_prime = 0;

for k = 1:n
    if(f_vec(k) > f_prime)
        f_prime = k/(2*pi);
    end
end

f_prime = 2*l_trans/max(v_phase);

figure;
subplot(3, 1, 1);
plot(freq, Z_dut);
hold on;
plot(freq, Z_dut_d, 'red');
title('Magnitude');
xlabel('Frequency (Rad/s)');
legend('Blue = New Cable', 'Red = Degredaed Cable');

subplot(3, 1, 2);
plot(freq, 180*theta/pi);
hold on;
plot(freq, 180*theta_d/pi, 'red');
title('Phase');
xlabel('Frequency (Rad/s)');

subplot(3, 1, 3)
plot(max(v_phase)*t/2, imag(h));
hold on;
plot(max(v_phase_d)*t/2, imag(h_d), 'red');
title(' Time Domain ');
xlabel('Cable Line Position (m)');

figure;
plot(freq, theta_d-theta);
title('Phase Differrence Between New & Degraded Cables');
xlabel('Frequency (Rad/s)');
ylabel('Phase new - Phase degraded');
ylim([-.5 .5]);