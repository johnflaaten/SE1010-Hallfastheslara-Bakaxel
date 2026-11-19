%%%% Lastfall b

% Parametrar
% Anger deras belopp och är snittet av parameterintervallet
% När bara belopp anges behöver man vara konsekvent med tecken i ekvationerna

run konstanter.m

%% Normalkraft bakhjulen - Punkt I

FL = 0.5 * rho_luft * c * A * v_b^2;                                % Luftmotstånd

Nb = (FL*(h1+h) + m*g*df - m*a2*h) / (db + df);                     % Normalkraft för båda bakhjul 
Nb1 = Nb / 2;                                                       % Normalkraft för 1 hjul

%% Tvär- och Lagerkrafter

% y-led
FB = m*a2 - FL;                                                     % Bromskraft vid hjulens kontakt med underlag (Stort B)
FB1 = FB / 2;                                                       % För ett hjul

Fb = dh*FB / (2*rb);                                                % Kraft från bromsskivor (Kommer från att hjul och skiva har samma men motsatt vridmoment)
Fb1 = Fb / 2;                                                       % För en bromsskiva

Ry = Fb - FB;                                                       % Total lagerkraft
Ry1 = Ry / 2;                                                       % För ett lager

% z-led
Rz = Nb;                                                            % Lagerkraft vertikalt
Rz1 = Nb1;

%% Normalspänningsdiagram

% x = [0 L];
% normalStress = [0 0];                                               % Finns ingen normalkraft i lastfall b
% 
% figure('Name', 'Lastfall b');                                       
% tiledlayout(3, 2);                                                  % Skapar en plotmatris med form 3x2
% 
% nexttile; plot(x, normalStress, "b-"); grid on; title("Normalkrafter");
% ylabel('Kraft [N]'); xlabel('x [m]');

%% Tvärkraftsdiagram

% y-led
xPositions = linspace(0, L, L*1000);                                % Antal steg = L*1000 så att steglängd är 1 mm
shearForceY = -FB1 + Fb1 *(xPositions>=bb) - Ry1 *(xPositions>=b1) - Ry1 *(xPositions>=L-b1) + Fb1 *(xPositions>=L-bb) - FB1 * (xPositions==L);
                          % ^^^^^^^^^^^^^ Gör att termen endast appliceras när villkoret uppfylls

figure('Name', 'Lastfall b');                                       
tiledlayout(3, 2); 
nexttile;  hold on; 
plot(xPositions, shearForceY, "b"); 
title("Tvärkrafter"); ylabel('Tvärkraft [N]'); 
xlabel('x [m]'); grid on;

% z-led
shearForceZ = Nb1 - Rz1 *(xPositions>=b1) - Rz1 *(xPositions>=L-b1) + Nb1 *(xPositions==L);

plot(xPositions, shearForceZ, "r");


% Resultanttvärkraft
shearForceYZ = hypot(shearForceY, shearForceZ);                     % Total tvärkraft

plot(xPositions, shearForceYZ, 'g--');
legend("V_{y}", "V_{z}", 'V_{res}');

%% Böjmomentsdiagram

bendingmomentY = cumtrapz(xPositions, shearForceZ);                 % Använder kumulativa integralen av tvärkraftsdiagrammet
bendingmomentZ = cumtrapz(xPositions, shearForceY);
bendingTotal = hypot(bendingmomentY, bendingmomentZ);               % Resulterande böjmoment     

nexttile; hold on;
plot(xPositions, bendingmomentY, "b");
plot(xPositions, bendingmomentZ, "r");
plot(xPositions, bendingTotal, "g--");
title("Böjmoment"); ylabel("Böjmoment [Nm]"); xlabel("x [m]");
legend("M_{y}", "M_{z}", 'M_{b}'); grid on;

%% Vridmomentsdiagram

torque = -FB1*dh/2 + Fb1*rb * (xPositions>=bb) + Fb1*rb * (xPositions>=L-bb) - FB1*dh/2 * (xPositions>=L);  
% Vridmoment från däcken är negativ m.a.p x-axeln och från bromsskivorna motsatt

nexttile; plot(xPositions, torque, 'b'); grid on; title('Vridmoment');
xlabel('x [m]'); ylabel('T_x [Nm]');


%% Nominella spänningar - Punkt II
D_guess = 0.050;
i_1 = round(b1*1000); i_2 = round((L-b1)*1000); i_3 = round(L*1000);    % Indexgränserna för olika sektioner
xPos = [xPositions(1:i_1), xPositions(i_1:i_2), xPositions(i_2:i_3)];   % I punkterna b1 och L-b1 vill jag ha två värden fårn lilla och stora axeln

% Böjspänning
sigma_b = 32/pi/D_guess^3*[bendingTotal(1:i_1)/0.6^3, bendingTotal(i_1:i_2), bendingTotal(i_2:i_3)/0.6^3];  % Spänning_b = My/I (y är radien på axeln)

% Vridskjuvspänning
tau_T = 16/pi/D_guess^3*[torque(1:i_1)/0.6^3, torque(i_1:i_2), torque(i_2:i_3)/0.6^3]; % Spänning_T = Tr/J (r är radien på axeln)

% Böjskjuvspänning
tau_V = 16/3/pi/D_guess^2*[shearForceYZ(1:i_1)/0.6^2, shearForceYZ(i_1:i_2), shearForceYZ(i_2:i_3)/0.6^2]; % Spänning_V =  

% OBS! Böjspänning och Vridskjuvspänning är för axelns ytteryta och Böjskjuvspänning för axelns centrum

% Diagram
nexttile; hold on;
plot(xPos, abs(sigma_b), 'b'); plot(xPos, abs(tau_T), 'r'); plot(xPos, abs(tau_V), 'g');
title('Nominella spänningar, Antagen diameter: 50 mm'); grid on; ylim([-1e7, 5e7]);
xlabel('x [m]'); ylabel('Spänning [Pa]'); legend('Böjspänning','Vridskjuvspänning','Böjskjuvspänning');

%% Nominell Effektivspänning - Beräkna diameter då spänningskoncentrationer ej beaktas

% Beräkna diameter för nominell effektivspänning m.a.p axelns ytteryta
sigma_allowed = sigma_s / ns;                                       % Spänningsgränsen med säkerhetsfaktor mot plastisk deformation
 
% Faktor nominell effektivspänning på lilla axeln (von Mises) (Första sektionen)                           
S1 = sqrt(bendingTotal(1:i_1).^2 + 3/4*torque(1:i_1).^2);           % Har brytit ut 32/pi/d^3 och tagit bort den, det som blir kvar är en konstant faktor
[sigma_max_small, idx_sigma_small] = max(S1);                       % Beräkna största spänning och position

% Faktor nominell effektivspänning på stora axeln (Andra sektionen)
S2 = sqrt(bendingTotal(i_1:i_2).^2 + 3/4*torque(i_1:i_2).^2);
[sigma_max_large, idx_sigma_large] = max(S2);   
idx_sigma_large = i_1 + idx_sigma_large - 1;

% Faktor nominell effektivspänning på lilla axeln (Tredje sektionen)
S3 = sqrt(bendingTotal(i_2:i_3).^2 + 3/4*torque(i_2:i_3).^2);

% Beräkna diameter
D_small_nom = (32/(pi*sigma_allowed) * sigma_max_small)^(1/3);    % Omskrivning av von Mises utan att beakta spänningskoncentrationer 
D_large_nom = (32/(pi*sigma_allowed) * sigma_max_large)^(1/3);
 
D_nom = max([D_large_nom D_small_nom/0.6]);                   % Väljer den som uppfyller båda

% Nominell effektivspänningsdiagram
sigma_nom = [32/(pi*(D_nom*0.6)^3)*S1, 32/(pi*D_nom^3)*S2, 32/(pi*(D_nom*0.6)^3)*S3];

nexttile; hold on; grid on;
plot(xPos, sigma_nom, 'g-');


%% Effektivspänning - Beräkna diameter med spänningskoncentrationer - Punkt III o IV

% Beräkna spänningskoncentrationsfaktor Ktb och KtT
% d = (32/pi/sigma_allowed*sqrt((K_tb_b*bendingTotal(i_1+1))^2 + 3/4*(K_tT_b * torque(i_1+1))^2))^(1/3);
% disp(d/0.6*1000)
% disp(d*1000)

% Beräkna Diameter
D_small = (32/pi/sigma_allowed*sqrt((K_tb_b*bendingTotal(i_1+1))^2 + 3/4*(K_tT_b * torque(i_1+1))^2))^(1/3); % Omskrivning av lokal von Mises - Kollar vid lager

if D_small < D_small_nom                                            % Ifall koncentrerade spänningen < max(spänning_nom) kommer dom ha samma diameter
    D_small = D_small_nom; end
D = D_small / 0.6;

conc_stress = 32/pi/(0.6*D)^3*sqrt((K_tb_b*bendingTotal(i_1))^2 + 3/4*(K_tT_b * torque(i_1))^2);

% Effektivspänningsdiagram
S1 = 32/pi/(0.6*D)^3*sqrt((bendingTotal(1:i_1)).^2 + 3/4*torque(1:i_1).^2);       % Nominell effektivspänning på första sektionen
S1(end) = conc_stress;                                                            % Vid lager är spänningskoncentrationen

S2 = 32/pi/D^3*sqrt(bendingTotal(i_1:i_2).^2 + 3/4*torque(i_1:i_2).^2);

S3 = 32/pi/(0.6*D)^3*sqrt(bendingTotal(i_2:i_3).^2 + 3/4*torque(i_2:i_3).^2);
S3(1) = conc_stress;

plot(xPos, [S1, S2, S3], 'm--');
yline(sigma_allowed, 'r--'); ylim([0 sigma_allowed*1.15]);
title(sprintf('Nominell Effektivspänning D_{nom} = %.2f mm vs Effektivspänning D_{eff} = %.2f mm',D_nom*1000, D*1000))
xlabel('x [m]'); ylabel('\sigma_{eff} [Pa]');
legend('\sigma_{eff,nom}', '\sigma_{eff}', '\sigma_{tillåten}');

%% Utmattningsanalys - Roterande Böjning - Punkt V
sigma_fatigue_allowed = sigma_u / nu;                       % Gräns för roterande böjmomentsspänning innan utmattning

% Beräkna Anvisningsfaktorn Kf
d = (K_f_b*32/(sigma_fatigue_allowed*pi)*bendingTotal(i_1+1))^(1/3);
disp(d/0.6*1000)
disp(d*1000)

% Nominell roterande böjspänning
d_fat_nom = max((32/(sigma_fatigue_allowed*pi)*bendingTotal(1:i_2)).^(1/3));
D_fat_nom = d_fat_nom / 0.6;

% Roterande böjningspänning med anvisningseffekt
d_fat = (K_f_b*32/(sigma_fatigue_allowed*pi)*bendingTotal(i_1+1))^(1/3);

if d_fat < d_fat_nom                                        % Ifall koncentrerade spänningen < max(spänning_nom) kommer dom ha samma diameter
    d_fat = d_fat_nom; end


D_fat = d_fat / 0.6;
sigma_fatigue_max = K_f_b*32/pi/(d_fat)^3*bendingTotal(i_1);

% Roterande Böjningsspänningsdiagram

sigma_rotation = [32/pi/(d_fat)^3*bendingTotal(1:i_1-1), sigma_fatigue_max...
                  32/pi/(D_fat)^3*bendingTotal(i_1:i_2), sigma_fatigue_max...
                  32/pi/(d_fat)^3*bendingTotal(i_2+1:i_3)];

sigma_rotation_nom = [32/pi/(d_fat_nom)^3*bendingTotal(1:i_1),...
                      32/pi/(D_fat_nom)^3*bendingTotal(i_1:i_2),...
                      32/pi/(d_fat_nom)^3*bendingTotal(i_2:i_3)];

nexttile; hold on; grid on;
plot(xPos, sigma_rotation_nom, 'g-')
plot(xPos, sigma_rotation, 'm--');
yline(sigma_fatigue_allowed, 'r--'); ylim([0, 1.15*sigma_fatigue_allowed]);
title(sprintf('Nominell Roterande böjspänning, D_{utm,nom} = %.2f mm vs Roterande böjspänning, D_{utm} = %.2f mm', D_fat_nom*1000, D_fat*1000));

xlabel('x [m]'); ylabel('\sigma_a [Pa]');
legend('\sigma_a_{,nom}', '\sigma_a', '\sigma_{tillåten}');