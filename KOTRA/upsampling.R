getwd()

#case1
num_up = table(final_case1$y)[2] - table(final_case1$y)[1]
samp_idx =1:table(final_case1$y)[1]
samp_idx = sample(samp_idx,num_up,replace=TRUE)
up_sampled = final_case1[final_case1$y == 0,][samp_idx,]
dim(up_sampled)

table(final_case1$y)
final_case1_up = rbind(final_case1,up_sampled)
dim(final_case1_up)

table(final_case1$y)
table(final_case1_up$y)

#case2
num_up = table(final_case2_imp$y)[2] - table(final_case2_imp$y)[1]
samp_idx =1:table(final_case2_imp$y)[1]
samp_idx = sample(samp_idx,num_up,replace=TRUE)
up_sampled = final_case2_imp[final_case2_imp$y == 0,][samp_idx,]
dim(up_sampled)

table(final_case2_imp$y)
final_case2_up = rbind(final_case2_imp,up_sampled)
dim(final_case2_up)

table(final_case2_imp$y)
table(final_case2_up$y)


#case4
num_up = table(final_case4_imp$y)[2] - table(final_case4_imp$y)[1]
samp_idx =1:table(final_case4_imp$y)[1]
samp_idx = sample(samp_idx,num_up,replace=TRUE)
up_sampled = final_case4_imp[final_case4_imp$y == 0,][samp_idx,]
dim(up_sampled)

table(final_case4_imp$y)
final_case4_up = rbind(final_case4_imp,up_sampled)
dim(final_case4_up)

table(final_case4_imp$y)
table(final_case4_up$y)

save(final_case1_up,final_case2_up,final_case4_up,file="case_data_upsampled.Rdata")

