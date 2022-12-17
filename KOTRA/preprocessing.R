rm(list=ls())
getwd()
# oldloc = Sys.getlocale("LC_ALL")
# Sys.setlocale("LC_ALL", "korean")

library(readxl)
library(tidyverse)
library(missForest)

#--------------------------------------------------------------------#
# data reading
#--------------------------------------------------------------------#

dat20 = read_excel("data/최종데이터_예측용_20_전달용.xlsx")[-1] # index 제거 후 불러오기
matrix(colnames(dat20),ncol=1)
head(dat20)

# column명 전부 영어로 변경

colnames(dat20) = c(read.table("data/colnames.txt")[[1]])
colnames(dat20)
head(dat20)
dim(dat20) # 73303    37
#View(dat20)
str(dat20)
##
#--------------------------------------------------------------------#
# 함수
#--------------------------------------------------------------------#



# 연속형 범주형 나누기

divide = function(data){
  cont = c()
  cat = c()
  for(idx in 1:ncol(data)){
    char = is.character(data[,idx,drop=T])
    if(char == TRUE){
      cat = c(cat,colnames(data)[idx])
    }
    else{
      cont = c(cont,colnames(data)[idx])
    }
  }
  
  # cont
  # cat
  # length(cont)
  # length(cat)

  return(list("cont"=cont,"cat"=cat))
}

show_cat_tbl = function(names,data){
  
  for(i in 1:length(names)){
    name = names[i]
    cat("*************\n")
    cat(name,"\n")
    print(table(data[name]))
  }
}

#--------------------------------------------------------------------#
# case 분류
#--------------------------------------------------------------------#
# case 1 kotra
dat20_case1 = dat20 %>% select(y,
                               tot_exp,
                               nn_exp,
                               cnt_exp,
                               itm6,
                               itm4, 
                               itm2,
                               itm1,
                               memb)
dim(dat20_case1) # 9개

case1_lst = divide(data=dat20_case1)
case1_cont = case1_lst$cont
case1_cat = case1_lst$cat
case1_cont
case1_cat

show_cat_tbl(case1_cat,dat20_case1) # 범주형변수 테이블


# case1_mean = apply(case1_cont,2,mean)
# case1_sd = apply(case1_cont,2,sd)
# tmp= sweep(case1_cont, 2, case1_mean, "-")
# colSums(tmp)
# tmp = sweep(tmp,2, case1_sd, "/")
# apply(tmp, 2, sd)
# scaled_case1 = tmp
# summary(scaled_case1)
# 

# continuous variables(item1 ~ item4 are removed)
case1_cont[1:4]
final_case1 = dat20_case1[case1_cont[1:4]]

# categorical variables

final_case1$y = factor(dat20_case1[case1_cat][,1,drop=TRUE],labels=c(0,1)) # 지속: 1
final_case1$memb = factor(dat20_case1[case1_cat][,2,drop=TRUE],labels=c(0,1)) # 지속: 1

head(final_case1)

#--------------------------------------------------------------------#
## case 2 kotra
idx_kotra = (dat20$memb == "회원")
sum(idx_kotra)

dat20_case2 = dat20[idx_kotra, ] %>% select(y, 
                                            tot_exp,
                                            nn_exp,
                                            cnt_exp,
                                            itm6,
                                            itm4, 
                                            itm2,
                                            itm1,
                                            rev,
                                            op,
                                            np,
                                            age,
                                            empl,
                                            webp,
                                            entp_last,
                                            entp_curr,
                                            entp_spec,
                                            voc,
                                            bk_num,
                                            bk_num_unq)
dim(dat20_case2) # 20개
case2_lst = divide(dat20_case2)
case2_cont = case2_lst$cont
case2_cat = case2_lst$cat
case2_cont
case2_cat

show_cat_tbl(case2_cat,dat20_case2)


# continuous variables(item1 ~ item4 are removed)
case2_cont[! case2_cont %in% c("itm1","itm2","itm4")]
final_case2 = dat20_case2[case2_cont[! case2_cont %in% c("itm1","itm2","itm4")]]

# categorical variables

final_case2$y = factor(dat20_case2[case2_cat][,1,drop=TRUE],labels=c(0,1)) # 지속: 1
final_case2$webp = factor(dat20_case2[case2_cat][,2,drop=TRUE],labels=c(0,1)) # 지속: 1


# age,empl을 채우고, rev NA는 제거
temp_case2 = final_case2[!is.na(final_case2$rev),]
summary(temp_case2) # rev NA 제거 : age(1802), empl(451) NA 존재

##### MICE IMPUTATION(age,empl) #####
for_imp = temp_case2 %>% select(-y)
for_imp = for_imp[!is.na(for_imp$rev),]



case2_imp = missForest(data.matrix(for_imp))
summary(case2_imp$ximp)

hist(temp_case2$age)
hist(case2_imp$ximp[,"age"])
hist(temp_case2$empl)
hist(case2_imp$ximp[,"empl"])
summary(temp_case2$age)
summary(case2_imp$ximp[,"age"])
summary(temp_case2$empl)
summary(case2_imp$ximp[,"empl"])

plot(temp_case2$age,log(temp_case2$rev+1),main="before")
plot(case2_imp$ximp[,"age"],log(case2_imp$ximp[,"rev"]+1),main="after")

plot(log(temp_case2$empl+1),log(temp_case2$rev+1),main="before")
plot(log(case2_imp$ximp[,"empl"]+1),log(case2_imp$ximp[,"rev"]+1),main="after")

final_case2_imp= tibble(data.frame(case2_imp$ximp))
final_case2_imp$y = temp_case2$y
summary(final_case2_imp)

#--------------------------------------------------------------------#
## case 4 kotra_entp
idx_entp = (dat20$entp_last>0) | (dat20$entp_curr>0)
sum(idx_entp)

dat20_case4 = dat20[idx_entp, ] %>% select(y, 
                                           tot_exp,
                                           nn_exp,
                                           cnt_exp,
                                           itm6,
                                           itm4, 
                                           itm2,
                                           itm1,
                                           rev,
                                           op,
                                           np,
                                           age,
                                           empl,
                                           webp,
                                           entp_last,
                                           entp_curr,
                                           entp_spec,
                                           voc,
                                           bk_num,
                                           bk_num_unq,
                                           ptcp_exp,
                                           ptcp_vch,
                                           ptcp_cnslt,
                                           ptcp_wc,
                                           ptcp_inq,
                                           ptcp_expo,
                                           ptcp_rgof,
                                           ptcp_vdc,
                                           buy_tot,
                                           buy_spc)
dim(dat20_case4)
case4_lst = divide(dat20_case4)
case4_cont = case4_lst$cont
case4_cat = case4_lst$cat
case4_cont
case4_cat


show_cat_tbl(case4_cat,dat20_case4)


# continuous variables(item1 ~ item4 are removed)
case4_cont[! case4_cont %in% c("itm1","itm2","itm4")]
final_case4 = dat20_case4[case4_cont[! case4_cont %in% c("itm1","itm2","itm4")]]

# categorical variables
case4_cat # 사실상 전부 binary라서 factor로 변경안해도 됨. y빼고
dat20_case4[case4_cat]
final_case4 = cbind(final_case4,dat20_case4[case4_cat])
summary(final_case4)
final_case4$y = factor(final_case4$y,labels=c(0,1))

# age,empl을 채우고, rev NA는 제거
temp_case4 = final_case4[!is.na(final_case4$rev),]
summary(temp_case4) # rev NA 제거 : age(1802), empl(451) NA 존재
dim(temp_case4)
idx_entp = (final_case2_imp$entp_last>0) | (final_case2_imp$entp_curr>0)
for_case4 = final_case2_imp[idx_entp,]%>% select(age,empl)

dim(temp_case4)
summary(temp_case4)

temp_case4$age = for_case4$age
temp_case4$empl = for_case4$empl
summary(temp_case4)

final_case4 = temp_case4
final_case4_imp = tibble(final_case4)

save(final_case1,final_case2_imp,final_case4_imp,file="case_data.Rdata")



