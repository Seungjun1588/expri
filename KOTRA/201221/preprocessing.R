#rm(list=ls())
getwd()
# 기존 데이터에 index와 작년수출금액 추가
# index가 겹치는 경우가 있음에 유의
# 결측치만 채우고  upsampling은 하지 말고
# 결과 : 인덱스 | 확률 | 예측결과
#       confusion matrix 
library(readxl)
library(tidyverse)
library(missForest)
library(mice)
library(naniar)
library(foreach)
library(doParallel)
#--------------------------------------------------------------------#
# data reading
#--------------------------------------------------------------------#


ex_cont20 = read_excel("연대_전달용2_1221/최종데이터_수출대륙_20_전달용.xlsx")[-1] # index 제거 후 불러오기
ex_item20 = read_excel("연대_전달용2_1221/최종데이터_수출품목_20_전달용.xlsx")[-1] # index 제거 후 불러오기
dat20 = read_excel("연대_전달용2_1221/최종데이터_예측용_20_전달용_인덱스추가.xlsx")[-1] # index 제거 후 불러오기

matrix(colnames(dat20),ncol=1) # ind와 작년수출금액 추가된 것을 확인.


# column명 전부 영어로 변경 ; ind와 ly_exp 추가

colnames(dat20) = c(read.table("colnames.txt")[[1]])
colnames(dat20)
head(dat20)
dim(dat20) # 73303    39
#View(dat20)
str(dat20)


#--------------------------------------------------------------------#
# 최종데이터_예측용_20_전달용_인덱스추가(dat20) 추가정보 확인
#--------------------------------------------------------------------#
# ind in dat20
length(unique(dat20$ind)) # 3개 정도가 인덱스가 중복되는 것으로 관찰된다. 
# 근데 어떻게 겹칠 수가 있는거지?
# 그럼 완전히 똑같은 데이터가 2번씩 들어가 있는건가?

# 확인
# unique_idx = rep(NA,nrow(dat20))
# for(i in 1:nrow(dat20)){
#   if(dat20$ind[i] %in% unique_idx){
#     print(dat20$ind[i])
#     print(dat20[dat20$ind == dat20$ind[i] ,])
#     print("********************")
#   }
#   
#   unique_idx[i] = dat20$ind[i]
# }
# ind 30318 하나만 4번 관찰된다. 
rm(unique_idx)

View(dat20[dat20$ind == 30318,])
# 차이점은 mn_td webp만 다르다. 



#### 일단 하나만 남기고 뺀다.
# 30319 ~ 30321 were removed
dat20[c(30319:30321),]
udat20 = dat20[-c(30319:30321),]
dim(udat20)
udat20[udat20$ind == 30318,]

#--------------------------------------------------------------------#
# 최종데이터_수출대륙_20_전달용(ex_cont20) 추가정보 확인
#--------------------------------------------------------------------#
head(ex_cont20)
length(unique(ex_cont20$ind))

# ex_cont20이 ind가 훨씬 많다..?
# NA도 있다
min(ex_cont20$ind)
max(ex_cont20$ind)

min(udat20$ind)
max(udat20$ind)

#  udat20 기준으로 ind를 골라서 사용해야할 듯?
# arrange
arncont20 = ex_cont20 %>% arrange(ind)
# dummy
start = min(arncont20$ind)
end = max(arncont20$ind)
dumm_len = length(unique(arncont20$cntt_name))
dummy = matrix(0,ncol=dumm_len,nrow=end)
colnames(dummy) = unique(arncont20$cntt_name)


# preprocessing
colnames(dummy)[8] ="NA"
colnames(dummy)
arncont20$cntt_name[is.na(arncont20$cntt_name)] = "NA"
sum(is.na(arncont20$cntt_name))

# making dummies
make_dummies = function(dummy_mat,data){
  i =1
  start = min(data$ind)
  end = max(data$ind)
  
  while(TRUE){# nrow(arncont20)
    if(i > nrow(data)) break
    if(start == data$ind[i]){
      dummy_mat[start,colnames(dummy_mat) == data$cntt_name[i]] = data$amt[i]
      i = i + 1
    }
    else{
      start = start + 1
    }
  }
  return(dummy_mat)
}
dummies = make_dummies(dummy,arncont20)
dim(dummies)
head(dummies)
head(arncont20)
tail(dummies)
tail(arncont20)

contdummies = cbind("ind"=1:nrow(dummies),dummies)
contdummies = data.frame(contdummies)

#--------------------------------------------------------------------#
# 최종데이터_수출품목_20_전달용(ex_cont20) 추가정보 확인
#--------------------------------------------------------------------#
head(ex_item20)
length(unique(ex_item20$ind))
length(unique(ex_item20$hscd)) # 1~97인데 code 77은 없네?

# arrange
arnitem20 = ex_item20 %>% arrange(ind)
# dummy
start = min(arnitem20$ind)
end = max(arnitem20$ind)
dumm_len = length(unique(arnitem20$hscd))
dummy = matrix(0,ncol=dumm_len,nrow=end)
colnames(dummy) = sort(unique(arnitem20$hscd))
sum(is.na(arncont20)) # 결측치는 전혀 없다 !

make_dummies = function(dummy_mat,data){
  i =1
  start = min(data$ind)
  end = max(data$ind)
  
  while(TRUE){# nrow(arncont20)
    if(i > nrow(data)) break
    if(start == data$ind[i]){
      dummy_mat[start,colnames(dummy_mat) == data$hscd[i]] = data$amt[i]
      i = i + 1
    }
    else{
      start = start + 1
    }
  }
  return(dummy_mat)
}

dummies = make_dummies(dummy,arnitem20)
dim(dummies)
head(dummies)
head(arnitem20)
tail(dummies)
tail(arnitem20)

itemdummies = cbind("ind"=1:nrow(dummies),dummies)
itemdummies = data.frame(itemdummies)
colnames(itemdummies) = c("ind",paste0("hscd",c(1:76,78:97)))

#--------------------------------------------------------------------#
# merge
#--------------------------------------------------------------------#
dim(udat20)
length(unique(udat20$ind))

dim(contdummies)
length(unique(contdummies$ind))

dim(itemdummies)
length(unique(itemdummies$ind))

tdat20 = udat20 %>% left_join(contdummies,by="ind") %>% left_join((itemdummies),by="ind")
dim(tdat20)
str(tdat20)

#--------------------------------------------------------------------#
# making binary
#--------------------------------------------------------------------#
cha_names = c()
for(i in 1:ncol(tdat20)){
  if(is.character(tdat20[,i,drop=T])){
    cha_names =c(cha_names,colnames(tdat20)[i])
  }
  
}
cha_names

# fr_inv
table(tdat20$fr_inv)
tdat20$fr_inv[tdat20$fr_inv == "모름"] = -1
tdat20$fr_inv[tdat20$fr_inv == "아니오"] = 0
tdat20$fr_inv[tdat20$fr_inv == "예"] = 1
table(tdat20$fr_inv)

# ino
table(tdat20$ino)
tdat20$ino[tdat20$ino == "모름"] = -1
tdat20$ino[tdat20$ino == "N"] = 0
tdat20$ino[tdat20$ino == "Y"] = 1
table(tdat20$ino)

# mn_td
table(tdat20$mn_td)
tdat20$mn_td[tdat20$mn_td == "모름"] = -1
tdat20$mn_td[tdat20$mn_td == "무역"] = 0
tdat20$mn_td[tdat20$mn_td == "제조"] = 1
tdat20$mn_td[tdat20$mn_td == "제조 및 무역"] = 2
tdat20$mn_td[tdat20$mn_td == "기타"] = 3
table(tdat20$mn_td)


# memb
table(tdat20$memb)
tdat20$memb[tdat20$memb == "비회원"] = 0
tdat20$memb[tdat20$memb == "회원"] = 1
table(tdat20$memb)

# y
table(tdat20$y)
tdat20$y[tdat20$y == "중단"] = 0
tdat20$y[tdat20$y == "지속"] = 1
table(tdat20$y)
str(tdat20)

for(i in 1:ncol(tdat20)){
  if(is.character(tdat20[,i,drop=T])){
    tdat20[,i,drop=T] = as.numeric(tdat20[,i,drop=T]) 
  }
  
}
str(tdat20)


# Q1 : 인덱스 골라내서 사용하면 되는건가?
# Q2 : hscd의 경우 범주가 너무 많은데 대충 비슷한거끼리 묶어달라고 할까?
# Q3: 두 데이터셋 모두 더미로 만들 수는 있는데, 각각의 amount는 어떻게 반영해야할까?
# Q4 : confustion matrix 구성방법
# Q5 : index 어떻게 해야 쉽게 처리할까
# Q6 : testset을 만들어야할 것 같다. 


#--------------------------------------------------------------------#
# case 분류
#--------------------------------------------------------------------#
# case 1 kotra
dat20_case1 = tdat20 %>% select(y,
                               tot_exp,
                               nn_exp,
                               cnt_exp,
                               itm6,
                               itm4, 
                               itm2,
                               itm1,
                               memb,colnames(contdummies),colnames(itemdummies))

dim(dat20_case1) # 114

# NA 전부 그냥 제거한 데이터 (즉, 기존)
final_case1 = dat20_case1
final_case1_org = dat20_case1

#--------------------------------------------------------------------#
## case 2 kotra
idx_kotra = (tdat20$memb == 1)
sum(idx_kotra)

dat20_case2 = tdat20[idx_kotra, ] %>% select(y, 
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
                                            bk_num_unq,colnames(contdummies),colnames(itemdummies))
dim(dat20_case2) # 125

summary(dat20_case2)
# age,empl을 채우고, rev NA는 제거
temp_case2 = dat20_case2[!is.na(dat20_case2$rev),]
summary(temp_case2) # rev NA 제거 : age(1802), empl(448) NA 존재

##### MICE IMPUTATION(age,empl) #####
for_imp = temp_case2 %>% select(-y)

registerDoParallel(cores=10)
case2_imp = missForest(data.matrix(for_imp),parallelize="forests") # 
summary(case2_imp$ximp)

final_case2_imp= tibble(data.frame(case2_imp$ximp))
final_case2_imp$y = temp_case2$y
summary(final_case2_imp)

# Imputation diagnostics
ggplot(temp_case2,
       aes(x = age,
           y = tot_exp)) +
  geom_miss_point() + 
  theme_dark()
ggplot(temp_case2,
       aes(x = empl,
           y = nn_exp)) +
  geom_miss_point() + 
  theme_dark()



col1 = adjustcolor("blue",alpha=0.3)
col2 = adjustcolor("red",alpha=0.3)
pdf("image/ipm_result_age.pdf") 
hist(temp_case2$age,ylim=c(0,3000),breaks=50,col=col1,main="Imputation result;Age",xlab="Age")
hist(final_case2_imp$age,add=TRUE,breaks=50,col=col2)
legend("topright",c("Before","After"),fill=c(col1,col2))
dev.off()

summary(final_case2_imp$age[is.na(for_imp$age)]) # 채워진 값들에 대한 summary
summary(temp_case2$age) # before
summary(final_case2_imp$age) #  after

#---------------
pdf("image/ipm_result_empl.pdf") 
hist(temp_case2$empl[temp_case2$empl<50],breaks=50,col=col1,main="Imputation result;empl",xlab="empl",freq=F)
hist(final_case2_imp$empl[final_case2_imp$empl<50],add=TRUE,breaks=50,col=col2,freq=F)
legend("topright",c("Before","After"),fill=c(col1,col2))
dev.off()

summary(final_case2_imp$empl[is.na(for_imp$empl)]) # 채워진 값들에 대한 summary
summary(temp_case2$empl)
summary(final_case2_imp$empl)


# NA 전부 그냥 제거한 데이터 (즉, 기존)
final_case2_org = na.omit(dat20_case2)
dim(final_case2_imp)
dim(final_case2_org)
#--------------------------------------------------------------------#
## case 4 kotra_entp
idx_entp = (tdat20$entp_last>0) | (tdat20$entp_curr>0)
sum(idx_entp)

dat20_case4 = tdat20[idx_entp, ] %>% select(y, 
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
                                           buy_spc,colnames(contdummies),colnames(itemdummies))
dim(dat20_case4)

# age,empl을 채우고, rev NA는 제거
temp_case4 = dat20_case4[!is.na(dat20_case4$rev),]
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
# NA 전부 그냥 제거한 데이터 (즉, 기존)
final_case4_org = na.omit(dat20_case4)

#---------------------------------------------------------------------------#
# ind 따로 저장.
final_case1_ind = final_case1 %>% select(ind)
final_case1 = final_case1 %>% select(-ind)
final_case1_org_ind = final_case1_org %>% select(ind)
final_case1_org = final_case1_org %>% select(-ind)

final_case2_imp_ind = final_case2_imp %>% select(ind)
final_case2_imp = final_case2_imp %>% select(-ind)
final_case2_org_ind = final_case2_org %>% select(ind)
final_case2_org = final_case2_org %>% select(-ind)

final_case4_imp_ind = final_case4_imp %>% select(ind)
final_case4_imp = final_case4_imp %>% select(-ind)
final_case4_org_ind = final_case4_org %>% select(ind)
final_case4_org = final_case4_org %>% select(-ind)

#---------------------------------------------------------------------------#
save(final_case1,final_case2_imp,final_case4_imp,file="201221/data/case_data_imp.Rdata")
save(final_case1_ind,final_case2_imp_ind,final_case4_imp_ind,file="201221/data/case_data_imp_ind.Rdata")
save(final_case1_org,final_case2_org,final_case4_org,file="201221/data/case_data_org.Rdata")






