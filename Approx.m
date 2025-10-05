C = Constants;
radialDistance = linspace(0, 3);

function d = D(r, z)
    d = sqrt(r^2 + z^2);
end

function dip = Dipole(r, z_r, z_v, sigma_tr, alpha)
    d_r = D(r, z_r);
    d_v = D(r, z_v); 
    term1 = (alpha*z_r*(1+sigma_tr*d_r)*exp(-1*sigma_tr*d_r))/(4*pi*(d_r^3));
    term2 = (alpha*z_v*(1+sigma_tr*d_v)*exp(-1*sigma_tr*d_r))/(4*pi*(d_v^3));
    dip = term1 - term2; 
end

function Multipole_Red = MultipoleRed(r, C)
    sigma_t = C.sigma_a_epi_red + C.sigma_s_epi;
    sigma_tr = sqrt(3*C.sigma_a_epi_red*sigma_t);
    Multipole_Red = Dipole(r, C.z_r_epi_red, C.z_v_epi_red, sigma_tr, C.alpha_epi_red) + Dipole(r, C.z_r_derm_red, C.z_v_derm_red, sigma_tr, C.alpha_derm_red)...
                    + Dipole(r, C.z_r_epi_red, C.z_v_epi_red, sigma_tr, C.alpha_epi_red) + Dipole(r, C.z_r_derm_red, C.z_v_derm_red, sigma_tr, C.alpha_derm_red);
end

function Multipole_Green = MultipoleGreen(r, C)
    sigma_t = C.sigma_a_epi_green + C.sigma_s_epi;
    sigma_tr = sqrt(3*C.sigma_a_epi_green*sigma_t);
    Multipole_Green = Dipole(r, C.z_r_epi_green, C.z_v_epi_green, sigma_tr, C.alpha_epi_green) + Dipole(r, C.z_r_derm_green, C.z_v_derm_green, sigma_tr, C.alpha_derm_green)...
                    + Dipole(r, C.z_r_epi_green, C.z_v_epi_green, sigma_tr, C.alpha_epi_green) + Dipole(r, C.z_r_derm_green, C.z_v_derm_green, sigma_tr, C.alpha_derm_green);
end

function Multipole_Blue = MultipoleBlue(r, C)
    sigma_t = C.sigma_a_epi_blue + C.sigma_s_epi;
    sigma_tr = sqrt(3*C.sigma_a_epi_blue*sigma_t);
    Multipole_Blue = Dipole(r, C.z_r_epi_blue, C.z_v_epi_blue, sigma_tr, C.alpha_epi_blue) + Dipole(r, C.z_r_derm_blue, C.z_v_derm_blue, sigma_tr, C.alpha_derm_blue)...
                    + Dipole(r, C.z_r_epi_blue, C.z_v_epi_blue, sigma_tr, C.alpha_epi_blue) + Dipole(r, C.z_r_derm_blue, C.z_v_derm_blue, sigma_tr, C.alpha_derm_blue);
end


fplot(@(r) MultipoleRed(r, C), [0, 5], 'Color', 'r', 'LineWidth', 2);
hold on
fplot(@(r) MultipoleGreen(r, C), [0, 5], 'Color', 'g', 'LineWidth', 2);
fplot(@(r) MultipoleBlue(r, C), [0, 5], 'Color', 'b', 'LineWidth', 2);
hold off
legend('Crveni profil', 'Zeleni profil', 'Plavi profil');
title('Difuzni profili')
xlabel('Radialna udaljenost (mm)')
ylabel('R(r)') 
