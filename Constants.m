 classdef Constants
    properties( Constant = true )
        z_r_epi_red = 0.8;  
        z_v_epi_red = -6.1023226346791475;
        z_r_derm_red = 1.05;
        z_v_derm_red = -6.588309429932237;

        z_r_epi_green = 0.8;  
        z_v_epi_green = -4.624960702437489;
        z_r_derm_green = 1.05;
        z_v_derm_green = -5.596793751087997;

        z_r_epi_blue = 0.8;  
        z_v_epi_blue = -3.638943642497341;
        z_r_derm_blue = 1.05;
        z_v_derm_blue = -5.144688561514153;


        sigma_s_epi = 0.5;
        sigma_s_derm = 0.75;

        sigma_a_epi_red = 0.127942724679755;
        sigma_a_derm_red = 0.03259680804656661;

        sigma_a_epi_green = 0.29894832784535047;
        sigma_a_derm_green = 0.20325706841266522;

        sigma_a_epi_blue = 0.4764179515058391;
        sigma_a_derm_blue = 0.30850865498704155;


        alpha_derm_red = 0.9583478903677984;
        alpha_epi_red = 0.7962509642181067;

        alpha_derm_green = 0.7867762273705217;
        alpha_epi_green = 0.6258227003846777;

        alpha_derm_blue = 0.7085440411530519;
        alpha_epi_blue = 0.5120757962600916;
    end
 end