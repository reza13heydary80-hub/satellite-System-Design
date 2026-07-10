clc; clear; close all;

%% =========================
%  CONSTANTS
% =========================
mu  = 398600.4418;       
Re  = 6378.137;          
J2  = 1.08263e-3;        
omega_earth = 7.2921159e-5; 

%% =========================
%  ORBIT DESIGN (SSO FOR IRAN)
% =========================
altitude = 650;              
a = Re + altitude;           
e = 0.01;                   

%True Sun-Synchronous inclination
i = deg2rad(97.8);           

%  LTAN selection (typical imaging orbit)
LTAN = 10.5; % 10:30 AM

%  Align orbit with Iran (~50°E)
RAAN0 = deg2rad(50 + 15*(12 - LTAN));  

omega = 0;      
M0 = 0;         

%% =========================
%  DERIVED PARAMETERS
% =========================
n = sqrt(mu/a^3);        
p = a*(1 - e^2);

%  J2 nodal precession (SSO condition)
RAAN_dot = - (3/2)*J2*(Re^2/p^2)*n*cos(i);

%% =========================
%  TIME
% =========================
T = 2*pi/n;
tspan = linspace(0, 10*T, 8000); % more orbits for dense ground track

%% =========================
%  STORAGE
% =========================
r_eci = zeros(3, length(tspan));
lat = zeros(1, length(tspan));
lon = zeros(1, length(tspan));

%% =========================
%  PROPAGATION
% =========================
for k = 1:length(tspan)

    t = tspan(k);

    % Mean anomaly
    M = M0 + n*t;

    % Solve Kepler
    E = M;
    for iter = 1:6
        E = E - (E - e*sin(E) - M)/(1 - e*cos(E));
    end

    % True anomaly
    nu = 2*atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));

    % Radius
    r = a*(1 - e*cos(E));

    % Orbital plane position
    r_pf = [r*cos(nu); r*sin(nu); 0];

    %  Time-varying RAAN (SSO physics)
    RAAN = RAAN0 + RAAN_dot * t;

    % Rotation matrices
    R3_W = [cos(RAAN) -sin(RAAN) 0;
            sin(RAAN)  cos(RAAN) 0;
            0          0         1];

    R1_i = [1 0 0;
            0 cos(i) -sin(i);
            0 sin(i)  cos(i)];

    R3_w = [cos(omega) -sin(omega) 0;
            sin(omega)  cos(omega) 0;
            0           0          1];

    Q = R3_W * R1_i * R3_w;

    % ECI position
    r_eci(:,k) = Q * r_pf;

    %% =========================
    %  ECI → ECEF
    % =========================
    theta = omega_earth * t;

    R3_theta = [ cos(theta) sin(theta) 0;
                -sin(theta) cos(theta) 0;
                 0          0          1];

    r_ecef = R3_theta * r_eci(:,k);

    %% =========================
    %  LAT / LON
    % =========================
    x = r_ecef(1); y = r_ecef(2); z = r_ecef(3);

    lon(k) = atan2(y, x);
    lat(k) = atan2(z, sqrt(x^2 + y^2));
end

lat = rad2deg(lat);
lon = rad2deg(lon);

%% =========================
%  FIX LONGITUDE JUMPS
% =========================
lon = mod(lon + 180, 360) - 180;

%% =========================
%  REALISTIC EARTH (TEXTURE)
% =========================
figure;
[X,Y,Z] = sphere(200);


% https://visibleearth.nasa.gov/images/57752/blue-marble-land-surface-shallow-water-and-shaded-topography
cdata = imread('land_ocean_ice_2048.jpg');

surf(Re*X, Re*Y, Re*Z, ...
    'FaceColor','texturemap', ...
    'CData',cdata, ...
    'EdgeColor','none');

hold on;
plot3(r_eci(1,:), r_eci(2,:), r_eci(3,:), 'r', 'LineWidth',1.5);

axis equal;
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('SSO Orbit (650 km) - Aligned with Iran');
grid on;
light; lighting phong;

%% =========================
%  GROUND TRACK (WITH IRAN)
% =========================
figure;
load coastlines;

plot(coastlon, coastlat, 'k'); hold on;
plot(lon, lat, '.r');

% Highlight Iran region (approx box)
lon_iran = [44 63 63 44 44];
lat_iran = [25 25 40 40 25];
plot(lon_iran, lat_iran, 'b', 'LineWidth',2);

xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
title('SSO Ground Track over Iran');
grid on;
xlim([-180 180]);
ylim([-90 90]);