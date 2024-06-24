To reproduce the results of Figures 8 to 13 in the supplement, please run

- Figure8_Figure11_order_15.m, Figure8_Figure11_order_20.m, Figure8_Figure11_order_25.m, Figure8_Figure11_order_28.m, and Figure11_Figure12_Figure13_compare_BIC_cor.m to store results of Figures 8 and 11 in the supplement
- Figure9_Figure12_wine_red.m, Figure9_Figure12_wine_white.m, and Figure11_Figure12_Figure13_compare_BIC_cor.m to store results of Figures 9 and 12 in the supplement
- Figure10_Figure13_piston_40.m, Figure10_Figure13_piston_60.m, Figure10_Figure13_piston_80.m, Figure10_Figure13_piston_90.m, and Figure11_Figure12_Figure13_compare_BIC_cor.m to store results of Figures 10 and 13 in the supplement

Other listed files are used in the main codes mentioned above:

- plma_kcv_linearM_approximate_lambda_CV_all_Annals.m, CV_plM_sim2_linear_Annals.m, and CV_plM_sim2_linear_real.m are used for the proposed method
- pairwise_loss.m, pairwise_loss_real.m, pairwise_loss_approximate_Annals_1.m, pairwise_loss_approximate_Annals_2.m, pairwise_loss_approximate_Annals_3.m, and pairwise_loss_approximate_Annals_4.m are used to calculate ranking loss
- parfor_prograss.m is used for parallel program
- kcv_nlsM_linear.m is used for the competing method proposed by Zhang and Liu (2023)
- order prority data.xlsx, group score.xlsx, and Piston100.m are the data that needs to be used in the main code
