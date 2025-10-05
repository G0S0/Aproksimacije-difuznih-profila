C = Constants;
g = @(r,v) (exp(-r.^2/2*v))./(2*pi*v);

sigma_t_green = C.sigma_a_epi_green + C.sigma_s_epi;
sigma_tr_green = sqrt(3*C.sigma_a_epi_green*sigma_t_green);

function d = D(r, z)
    d = sqrt(r.^2 + z^2);
end

function dip = Dipole(r, z_r, z_v, sigma_tr, alpha)
    d_r = D(r, z_r);
    d_v = D(r, z_v); 
    term1 = (alpha*z_r*(1+sigma_tr*d_r).*exp(-1*sigma_tr.*d_r))./(4*pi*(d_r.^3));
    term2 = (alpha*z_v*(1+sigma_tr*d_v).*exp(-1*sigma_tr.*d_r))./(4*pi*(d_v.^3));
    dip = term1 - term2; 
end

function Multipole_Green = MultipoleGreenFunc(r, C)
    sigma_t = C.sigma_a_epi_green + C.sigma_s_epi;
    sigma_tr = sqrt(3*C.sigma_a_epi_green*sigma_t);
    Multipole_Green = Dipole(r, C.z_r_epi_green, C.z_v_epi_green, sigma_tr, C.alpha_epi_green) + Dipole(r, C.z_r_derm_green, C.z_v_derm_green, sigma_tr, C.alpha_derm_green)...
                    + Dipole(r, C.z_r_epi_green, C.z_v_epi_green, sigma_tr, C.alpha_epi_green) + Dipole(r, C.z_r_derm_green, C.z_v_derm_green, sigma_tr, C.alpha_derm_green);
end


MultipoleGreen = @(r, C) ...
    Dipole(r, C.z_r_epi_green, C.z_v_epi_green, sigma_tr_green, C.alpha_epi_green) + Dipole(r, C.z_r_derm_green, C.z_v_derm_green, sigma_tr_green, C.alpha_derm_green)...
                    + Dipole(r, C.z_r_epi_green, C.z_v_epi_green, sigma_tr_green, C.alpha_epi_green) + Dipole(r, C.z_r_derm_green, C.z_v_derm_green, sigma_tr_green, C.alpha_derm_green);

rVals = linspace(0,4,501);

Single = @(params) sum( rVals .* ...
   (params(1) * g(rVals, abs(params(2))) - MultipoleGreen(rVals, C)).^2 );

two = @(params) sum( rVals .* ...
   (params(1) * g(rVals, abs(params(2))) + params(3) * g(rVals, abs(params(4))) - MultipoleGreen(rVals, C)).^2 );

four = @(params) sum( rVals .* ...
    (params(1) * g(rVals, abs(params(2))) + params(3) * g(rVals, abs(params(4))) + params(5) * g(rVals, abs(params(6))) +...
  + params(7) * g(rVals, abs(params(8))) - MultipoleGreen(rVals, C)).^2 );

six = @(params) sum( rVals .* ...
  (params(1) * g(rVals, abs(params(2))) + params(3) * g(rVals, abs(params(4))) + params(5) * g(rVals, abs(params(6))) +...
  + params(7) * g(rVals, abs(params(8))) + params(9) * g(rVals, abs(params(10))) + params(11) * g(rVals, abs(params(12))) - MultipoleGreen(rVals, C)).^2 );

eight = @(params) sum( rVals .* ...
   (params(1) * g(rVals, abs(params(2))) + params(3) * g(rVals, abs(params(4))) + params(5) * g(rVals, abs(params(6))) +...
  + params(7) * g(rVals, abs(params(8))) + params(9) * g(rVals, abs(params(10))) + params(11) * g(rVals, abs(params(12)))...
  + params(13) * g(rVals, abs(params(14))) + params(15) * g(rVals, abs(params(16))) - MultipoleGreen(rVals, C)).^2 );


options = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e5, 'Display', 'iter');

x02 = [0.15, 0.08, 0.15,  0.4];
%twoFit = fminsearch(two, x02, options);
%w12 = twoFit(1);
%v12 = twoFit(2);
%w22 = twoFit(3);
%v22 = twoFit(4);

x04 = [0.6, 0.04, 0.6, 0.08, 0.6, 0.16, 0.6, 0.32];
fourFit = fminsearch(four, x04, options);
w14 = fourFit(1)
v14 = fourFit(2)
w24 = fourFit(3)
v24 = fourFit(4)
w34 = fourFit(5)
v34 = fourFit(6)
w44 = fourFit(7)
v44 = fourFit(8)

%x06 = [0.5, 0.04, 0.5, 0.08, 0.6, 0.16, 0.6, 0.32, 0.6, 0.64, 0.4, 0.80];
%sixFit = fminsearch(six, x06, options);
%w16 = sixFit(1);
%v16 = sixFit(2);
%w26 = sixFit(3);
%v26 = sixFit(4);
%w36 = sixFit(5);
%v36 = sixFit(6);
%w46 = sixFit(7);
%v46 = sixFit(8);
%w56 = sixFit(9);
%v56 = sixFit(10);
%w66 = sixFit(11);
%v66 = sixFit(12);


%x08 = [0.8, 0.04, 0.4, 0.08, 0.3, 0.16, 0.6, 0.32, 0.5, 0.64, 0.7, 0.80, 0.7, 0.96, 0.7, 0.8];
%eightFit = fminsearch(eight, x08, options);
%w18 = eightFit(1);
%v18 = eightFit(2);
%w28 = eightFit(3);
%v28 = eightFit(4);
%w38 = eightFit(5);
%v38 = eightFit(6);
%w48 = eightFit(7);
%v48 = eightFit(8);
%%w58 = eightFit(9);
%v58 = eightFit(10);
%w68 = eightFit(11);
%v68 = eightFit(12);
%w78 = eightFit(13);
%%v78 = eightFit(14);
%w88 = eightFit(15);
%v88 = eightFit(16);

g8sum = @(r)...
    w18 *g(r, v18) + w28 *g(r, v28) + w38 *g(r, v38) + w44 *g(r, v48) + ...
    w58 *g(r, v58) + w68 *g(r, v68) + w78 *g(r, v78) + w88 *g(r, v88);

g6sum = @(r)...
    w16 *g(r, v16) + w26 *g(r, v26) + w36 *g(r, v36) + w46 *g(r, v46) + ...
    w56 *g(r, v56) + w66 *g(r, v66);

g4sum = @(r)...
    w14 *g(r, v14) + w24 *g(r, v24) + w34 *g(r, v34) + w44 *g(r, v44);

g2sum = @(r)...
    w12 *g(r, v12) + w22 *g(r, v22);

fplot(@(r) MultipoleGreenFunc(r, C), [0, 4], 'Color', 'g', 'LineWidth', 2);
hold on
%fplot(@(rs) g8sum(rs), [0, 4], 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2) %white
fplot(@(rs) g6sum(rs), [0, 4], 'Color', [1, 0, 1], 'LineWidth', 2) %magenta
fplot(@(rs) g4sum(rs), [0, 4], 'Color', [1, 1, 0], 'LineWidth', 2) %yellow
fplot(@(rs) g2sum(rs), [0, 4], 'Color', [0, 1, 1], 'LineWidth', 2) %cyan
hold off
legend('Zeleni profil', 'Suma 6 Gaussa', 'Suma 4 Guassa', 'Suma 2 Gaussa');
title('Aproksimacija zelenog difuznog profila')
xlabel('Radialna udaljenost (mm)')
ylabel('R(r)') 

weightSumNegative = w14 + w24 + w34 + w44

w14NormNegative = w14/weightSumNegative
w24NormNegative = w24/weightSumNegative
w34NormNegative = w34/weightSumNegative
w44NormNegative = w44/weightSumNegative

weightSumNormNegative = w14NormNegative+ w24NormNegative+ w34NormNegative+ w44NormNegative

normalizedGaussGreen = @(r)...
    w14NormNegative*g(r, v14) + w24NormNegative*g(r, v24) + w34NormNegative*g(r, v34) + ...
    w44NormNegative*g(r, v44);

fplot(@(rs) normalizedGaussGreen(rs), [0, 4], 'Color', [0, 1, 1], 'LineWidth', 2) %cyan