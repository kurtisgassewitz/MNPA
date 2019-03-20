
%Initialize variables and matricies

clear;

G = zeros(8,8); 
C = zeros(8,8); 
F = zeros(8,1); 
R1 =1;
R2 = 2;
R3 = 10;
R4 = 0.1;
R0 = 1000;
cap = 0.25;
L = 0.2;
alpha = 100;

% Create G matrix 
G(1,1) = 1/R1;
G(1,2) = -1/R1;
G(2,1) = -1/R1;
G(2,2) = 1/R1 + 1/R2;
G(3,3) = 1/R3;
G(4,4) = 1/R4;
G(4,5) = -1/R4;
G(5,4) = -1/R4;
G(5,5) = 1/R0 + 1/R4;
G(1,6) = 1;
G(6,1) = 1;
G(2,7) = 1;
G(7,2) = 1;
G(3,7) = -1;
G(7,3) = -1;
G(4,8)=1;
G(8,3) = -alpha/R3;
G(8,4) = 1;

% Create C matrix
C(1,1)= cap;
C(2,2)= cap;
C(1,2)= -cap;
C(2,1)= -cap;
C(7,7)= -L;

% Solve F matrix
Vin = linspace(-10,10,100);
V3 = zeros(length(Vin),1);
V0 = zeros(length(Vin),1);
for i = 1:length(Vin)
    F(6) = Vin(i);
    V = G\F;
    V3(i) = V(3);
    %V0 is voltage at N5
    V0(i) = V(5);
end

figure(1)
plot(Vin,V3);
xlabel('Vin')
ylabel('V3')
title('V3 vs Vin sweep')

figure(2)
plot(Vin,V0);
xlabel('Vin')
ylabel('V0')
title('V0 vs Vin sweep')

w = 2*pi*linspace(0,80,100);
V0 = zeros(length(w),1);
gain = zeros(length(w),1);

for i = 1:length(w)
    s = 1i*w(i);
    M = inv((G +((s).*C)))*F; 
    V0(i) = abs(M(5));
    gain(i) = 20*log10(abs(V0(i))/abs(M(1)));
end

figure(3)
plot(w,V0);
xlabel('w (rads/sec)')
ylabel('V0')
title('AC plot for V0')

figure(4)
semilogx(w,gain);
xlabel('w (rads/sec)')
ylabel('V0/V1 (dB)')
title('Gain');

V0 = zeros(length(w),1);
gain = zeros(length(w),1);

for i = 1:length(gain)
    pert = randn()*0.05;
    C(1,1) = cap*pert;
    C(1,2) = -cap*pert;
    C(2,1) = -cap*pert;
    C(2,2) = cap*pert;

    s = 1i*2*pi*pi;
    M = inv((G +((s).*C)))*F; 
    V0(i) = abs(M(5));
    gain(i) = 20*log10((V0(i))/abs(M(1)));
end

figure(5);
hist(gain,100);
xlabel('Gain')
ylabel('Counts')
title('Hist C')