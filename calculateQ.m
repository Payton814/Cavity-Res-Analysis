


function Q = calculateQ(n, l, a, b, delta)

%% This code performs a monte carlo integral for the magnetic field in the cavity resonance

%% n is the mode number
%% l is the length of the cavity
%% a is the width of the cavity
%% b is the height of the cavity
%% delta is the cavity walls skin depth

%% I'm assuming l > a > b

rng('shuffle')
% Monte carlo to determine effective volume of cylindrical sample in C-band cavity
%b = 29.08;  %mm  WR229
%a = 2*b;
%l = 12*25.4;

%% speed of light is 300 mm/ns
f0 = 300/2*sqrt(1/a^2 + (n/l)^2);


% origin at wall edge and end of WG
E0 = 1.0;

NumSamples = 10000000;
Hx2 = zeros(NumSamples, 1);
Hz2 = zeros(NumSamples, 1);

%% Cavity Specific Parameters
k = 2*pi*f0/300; % Wavenumber [1/mm]
Z0 = 376.730313412; % Impedance of free space [Ohms]
Beta = sqrt(k^2 - (pi/a)^2); %% This used to have a (n/l)^2 term and the pi wasnt there. This was changed 051825. Look at Pozar for details
Z_TE = k*Z0/Beta;
%fprintf('1 over Z_TE: %d\n', 1/Z_TE^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% H surface integral
%% There are 3 unique surface integrals

%% Perform integral on yz wall
z = unifrnd(0,l,NumSamples,1);
Hz2 = sin(n*pi*z/l).^2;
I1 = 2*mean(Hz2)*(pi/(k*Z0*a))^2*l*b;

%% Perform integral on xy 
x = unifrnd(0,a,NumSamples,1);
Hx2 = sin(pi*x/a).^2;
I2 = 2*mean(Hx2)*(1/Z_TE)^2*a*b;

%% Perform integral on xz
x = unifrnd(0,a,NumSamples,1);
z = unifrnd(0,l,NumSamples,1);
Hx2 = sin(pi*x/a).^2.*cos(n*pi*z/l).^2;
Hz2 = cos(pi*x/a).^2.*sin(n*pi*z/l).^2;
I3 = 2*(mean(Hx2)*(1/Z_TE)^2 + mean(Hz2)*(pi/(k*Z0*a))^2)*a*l;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Volume integral
I4 = (1/Z0^2)*a*b*l/4;




Q = (2/delta)*I4/(I1 + I2 + I3);


end