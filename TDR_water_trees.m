%% Time Domain Reflectometry With Water Trees
clear;
clc;

%% Local Water Tree Model
e_0 = 8.85e-12;
f = 60;                                     % Hz
cond_water = 1e-7;                         % conductivity of water
e_water = 81-i*cond_water/(2*pi*f*e_0);     % complex permiativity of water
e_xlpe = 2.3-i*.001;                        % complex permiativity of XLPE

% water tree
kw = .75;                                   % water tree concentration
hw = .4;                                    % moisture content
q_w = kw*hw;                                % water content in water tree
D = 1/4;                                    % geometry of water tree
e_wt = e_xlpe*(1+q_w*(e_water-e_xlpe)/(e_xlpe+D*(1-q_w)*(e_water-e_xlpe))); % water tree permeativity

w = .2;                                     % length parameter of water treed region
e_total = w*e_wt + (1-w)*e_xlpe;             % total permiativity of water treed region

%% Single Core Transmission Line Model
