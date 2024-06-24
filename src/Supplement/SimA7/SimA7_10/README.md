To reproduce the results of Figure 14 in the supplement, please run

- Figure14_10_50.m, Figure14_10_100.m, and Figure14_10_200.m to store results of Figure 14 in the supplement
- Note: When running the above codes to store the results, we need to change the path ("C:/Users/dell/Desktop") of the data  to the current desktop path of the reader. The rest of the path ("2023.0257/scr/...") remains unchanged.

Other listed files are used in the main codes mentioned above:

- plma_kcv_linear_approximate_lambda_CV_all_Annals_survival.m, CV_pl_sim2_linear_Annals_survival.m, and CV_pl_sim2_linear_real_survival.m are used for the proposed method
- pairwise_loss.m, pairwise_loss_real.m, pairwise_loss_approximate_Annals_1.m, pairwise_loss_approximate_Annals_2.m, pairwise_loss_approximate_Annals_3.m, and pairwise_loss_approximate_Annals_4.m are used to calculate ranking loss
- parfor_prograss.m is used for parallel program
- order.R and screen.R are used for the proposed method to select appropriate genes 
- kcv_nls_linear_survival.m is used for the competing method proposed by Zhang and Liu (2023)
- NKI_cleaned.xlsx is NKI data that needs to be used in the main code
