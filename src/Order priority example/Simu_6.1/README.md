To reproduce the results of Figure 3, please run
- Order15.m, Order20.m, Order25.m, and Order28.m to store results of Figure 3
- Note: When running the above codes to store the results, we need to change the path ("C:/Users/dell/Desktop") to the output data at the end of the code to the current desktop path of the reader. The rest of the path ("2023.0257/scr/output/...") remains unchanged

Other listed files are used in the main codes mentioned above:
- plma_kcv_linear_approximate_lambda_CV_all_Annals.m, CV_pl_sim2_linear_Annals.m, and CV_pl_sim2_linear_real.m are used for the proposed method
- pairwise_loss.m, pairwise_loss_real.m, pairwise_loss_approximate_Annals_1.m, pairwise_loss_approximate_Annals_2.m, pairwise_loss_approximate_Annals_3.m, and pairwise_loss_approximate_Annals_4.m are used to calculate ranking loss
- parfor_prograss.m is used for parallel program
- kcv_nls_linear.m is used for the competing method proposed by Zhang and Liu (2023)
- order prority data.xlsx is order priority data that needs to be used in the main code
