#前置準備動作(倉庫大小，物品種類.數量.長寬高)
warehouse <- c(40,30,30)#倉庫長寬高
commodity_inf <- matrix(c(10,20,25,5,10,10,15,5,20,5,20,25,5,15,10), nrow = 3, ncol = 5)#物品種類資訊(五種)
commodity_qua <- c(2,3,1,1,3)#對應五種物品種類的物品數量，總共10個物品
correspond_table <- c(1,1,2,2,3,3,3,4,4,4,5,6,7,8,9,9,9,10,10,10)#對應表，高不動，物品旋轉看作兩種物品，例如1跟2是同一個物品
mutation_probability <- 5#突變機率
frequency <- 1000#執行次數

#最大跑過31750

#初始化(亂數決定五個染色體分別擺入物品順序)
initialization <- function(){
  
  qua_1 <- commodity_qua
  qua_2 <- commodity_qua
  qua_3 <- commodity_qua
  qua_4 <- commodity_qua
  qua_5 <- commodity_qua
  for(i in c(1:10)){
    
    a1 <- sample(20,1)
    a2 <- sample(20,1)
    a3 <- sample(20,1)
    a4 <- sample(20,1)
    a5 <- sample(20,1)
    while(qua_1[ceiling(correspond_table[a1]/2)]<1){a1 <- sample(20,1)}
    while(qua_2[ceiling(correspond_table[a2]/2)]<1){a2 <- sample(20,1)}
    while(qua_3[ceiling(correspond_table[a3]/2)]<1){a3 <- sample(20,1)}
    while(qua_4[ceiling(correspond_table[a4]/2)]<1){a4 <- sample(20,1)}
    while(qua_5[ceiling(correspond_table[a5]/2)]<1){a5 <- sample(20,1)}
    qua_1[ceiling(correspond_table[a1]/2)] <- qua_1[ceiling(correspond_table[a1]/2)]-1
    qua_2[ceiling(correspond_table[a2]/2)] <- qua_2[ceiling(correspond_table[a2]/2)]-1
    qua_3[ceiling(correspond_table[a3]/2)] <- qua_3[ceiling(correspond_table[a3]/2)]-1
    qua_4[ceiling(correspond_table[a4]/2)] <- qua_4[ceiling(correspond_table[a4]/2)]-1
    qua_5[ceiling(correspond_table[a5]/2)] <- qua_5[ceiling(correspond_table[a5]/2)]-1
    chromosomes_1[i] <<- correspond_table[a1]
    chromosomes_2[i] <<- correspond_table[a2]
    chromosomes_3[i] <<- correspond_table[a3]
    chromosomes_4[i] <<- correspond_table[a4]
    chromosomes_5[i] <<- correspond_table[a5]
    
  }
  
}



#交配(使用Partially-mapped crossover，一共交配兩次，因此會產出四個新chromosome)
crossover <- function(){
  
  range_1 <- sort(sample(9,2))
  temp_1 <- c()
  temp_2 <- c()
  sel_1_cou <- commodity_qua
  sel_2_cou <- commodity_qua
  

  
  #第一次交配
  for(i in c((range_1[1]+1):range_1[2])){
    sel_1_cou[ceiling(select_2[i]/2)] <- sel_1_cou[ceiling(select_2[i]/2)]-1
    sel_2_cou[ceiling(select_1[i]/2)] <- sel_2_cou[ceiling(select_1[i]/2)]-1
    chromosomes_1[i] <<- select_2[i]
    chromosomes_2[i] <<- select_1[i]
  }
  for(i in c(1:range_1[1])){
    if(sel_1_cou[ceiling(select_1[i]/2)]>0){
      chromosomes_1[i]<<-select_1[i]
      sel_1_cou[ceiling(select_1[i]/2)] <- sel_1_cou[ceiling(select_1[i]/2)]-1
    }else{
      temp_1 <- c(temp_1,i)
    }
    if(sel_2_cou[ceiling(select_2[i]/2)]>0){
      chromosomes_2[i]<<-select_2[i]
      sel_2_cou[ceiling(select_2[i]/2)] <- sel_2_cou[ceiling(select_2[i]/2)]-1
    }else{
      temp_2 <- c(temp_2,i)
    }
  }
  for(i in c((range_1[2]+1):10)){
    if(sel_1_cou[ceiling(select_1[i]/2)]>0){
      chromosomes_1[i]<<-select_1[i]
      sel_1_cou[ceiling(select_1[i]/2)] <- sel_1_cou[ceiling(select_1[i]/2)]-1
    }else{
      temp_1 <- c(temp_1,i)
    }
    if(sel_2_cou[ceiling(select_2[i]/2)]>0){
      chromosomes_2[i]<<-select_2[i]
      sel_2_cou[ceiling(select_2[i]/2)] <- sel_2_cou[ceiling(select_2[i]/2)]-1
    }else{
      temp_2 <- c(temp_2,i)
    }
  }
  a<-1
  for(i in c((range_1[1]+1):range_1[2])){
    if(sel_1_cou[ceiling(select_1[i]/2)]>0){
      chromosomes_1[temp_1[a]]<<-select_1[i]
      sel_1_cou[ceiling(select_1[i]/2)] <- sel_1_cou[ceiling(select_1[i]/2)]-1
      a<-a+1
    }
  }
  a<-1
  for(i in c((range_1[1]+1):range_1[2])){
    if(sel_2_cou[ceiling(select_2[i]/2)]>0){
      chromosomes_2[temp_2[a]]<<-select_2[i]
      sel_2_cou[ceiling(select_2[i]/2)] <- sel_2_cou[ceiling(select_2[i]/2)]-1
      a<-a+1
    }
  }
  
  
  range_2 <- sort(sample(9,2))
  temp_1 <- c()
  temp_2 <- c()
  sel_1_cou <- commodity_qua
  sel_2_cou <- commodity_qua
  
  #第二次交配
  for(i in c((range_2[1]+1):range_2[2])){
    sel_1_cou[ceiling(select_2[i]/2)] <- sel_1_cou[ceiling(select_2[i]/2)]-1
    sel_2_cou[ceiling(select_1[i]/2)] <- sel_2_cou[ceiling(select_1[i]/2)]-1
    chromosomes_3[i] <<- select_2[i]
    chromosomes_4[i] <<- select_1[i]
  }
  for(i in c(1:range_2[1])){
    if(sel_1_cou[ceiling(select_1[i]/2)]>0){
      chromosomes_3[i]<<-select_1[i]
      sel_1_cou[ceiling(select_1[i]/2)] <- sel_1_cou[ceiling(select_1[i]/2)]-1
    }else{
      temp_1 <- c(temp_1,i)
    }
    if(sel_2_cou[ceiling(select_2[i]/2)]>0){
      chromosomes_4[i]<<-select_2[i]
      sel_2_cou[ceiling(select_2[i]/2)] <- sel_2_cou[ceiling(select_2[i]/2)]-1
    }else{
      temp_2 <- c(temp_2,i)
    }
  }
  for(i in c((range_2[2]+1):10)){
    if(sel_1_cou[ceiling(select_1[i]/2)]>0){
      chromosomes_3[i]<<-select_1[i]
      sel_1_cou[ceiling(select_1[i]/2)] <- sel_1_cou[ceiling(select_1[i]/2)]-1
    }else{
      temp_1 <- c(temp_1,i)
    }
    if(sel_2_cou[ceiling(select_2[i]/2)]>0){
      chromosomes_4[i]<<-select_2[i]
      sel_2_cou[ceiling(select_2[i]/2)] <- sel_2_cou[ceiling(select_2[i]/2)]-1
    }else{
      temp_2 <- c(temp_2,i)
    }
  }
  a<-1
  for(i in c((range_2[1]+1):range_2[2])){
    if(sel_1_cou[ceiling(select_1[i]/2)]>0){
      chromosomes_3[temp_1[a]]<<-select_1[i]
      sel_1_cou[ceiling(select_1[i]/2)] <- sel_1_cou[ceiling(select_1[i]/2)]-1
      a<-a+1
    }
  }
  a<-1
  for(i in c((range_2[1]+1):range_2[2])){
    if(sel_2_cou[ceiling(select_2[i]/2)]>0){
      chromosomes_4[temp_2[a]]<<-select_2[i]
      sel_2_cou[ceiling(select_2[i]/2)] <- sel_2_cou[ceiling(select_2[i]/2)]-1
      a<-a+1
    }
  }
  
}



#突變(使用Frame-Shift mutation移碼突變 && Whole-gene mutation全基因突變，全基因突變的突變結果只影響突變物品的旋轉擺放，例如:1突變只會變2)
mutation <- function(){
  
  
  a<-sample(1:100,size=4,replace=T)
  if(a[1]<=mutation_probability){
    temp <- chromosomes_1[10]
    temp2 <- chromosomes_1[c(1:9)]
    chromosomes_1 <<- c(temp,temp2)
  }
  if(a[2]<=mutation_probability){
    temp <- chromosomes_2[10]
    temp2 <- chromosomes_2[c(1:9)]
    chromosomes_2 <<- c(temp,temp2)
  }
  if(a[3]<=mutation_probability){
    temp <- chromosomes_3[10]
    temp2 <- chromosomes_3[c(1:9)]
    chromosomes_3 <<- c(temp,temp2)
  }
  if(a[4]<=mutation_probability){
    temp <- chromosomes_4[10]
    temp2 <- chromosomes_4[c(1:9)]
    chromosomes_4 <<- c(temp,temp2)
  }
  
  
  
  for(i in c(1:10)){
    a<-sample(100,1)
    if(a<=mutation_probability){
      if((chromosomes_1[i]%%2)==1){
        chromosomes_1[i] <<- chromosomes_1[i]+1
      }else{chromosomes_1[i] <<- chromosomes_1[i]-1}
    }
    a<-sample(100,1)
    if(a<=mutation_probability){
      if((chromosomes_2[i]%%2)==1){
        chromosomes_2[i] <<- chromosomes_2[i]+1
      }else{chromosomes_2[i] <<- chromosomes_2[i]-1}
    }
    a<-sample(100,1)
    if(a<=mutation_probability){
      if((chromosomes_3[i]%%2)==1){
        chromosomes_3[i] <<- chromosomes_3[i]+1
      }else{chromosomes_3[i] <<- chromosomes_3[i]-1}
    }
    a<-sample(100,1)
    if(a<=mutation_probability){
      if((chromosomes_4[i]%%2)==1){
        chromosomes_4[i] <<- chromosomes_4[i]+1
      }else{chromosomes_4[i] <<- chromosomes_4[i]-1}
    }
  }
  
  
}


#計算適應值(先利用"下後左角優先堆疊定則"確認每個染色體的排序可以依序擺入幾個物品，再去計算每個染色體物品的整體大小)
fitness <- function(){
  
  L <- warehouse[1]
  W <- warehouse[2]
  L_Max <- 0
  i<-1
  V<-0
  repeat {
    if(L>fitness_commodity[1,chromosomes_1[i]]){
      if(W>fitness_commodity[2,chromosomes_1[i]]){
        W <- W - fitness_commodity[2,chromosomes_1[i]]
        if(L_Max<fitness_commodity[1,chromosomes_1[i]]){
          L_Max <- fitness_commodity[1,chromosomes_1[i]]
        }
      }else{
        L <- L-L_Max
        L_Max <- 0
        W <- warehouse[2]
      }
    }else{
      break
    }
    i <- i+1
    if(i>10){
      break
    }
  }
  for(j in c(1:(i-1))){
    V<-V+fitness_commodity[1,chromosomes_1[j]]*fitness_commodity[2,chromosomes_1[j]]*fitness_commodity[3,chromosomes_1[j]]
  }
  fitness_value[1] <<- V
  
  
  L <- warehouse[1]
  W <- warehouse[2]
  L_Max <- 0
  i<-1
  V<-0
  repeat {
    if(L>fitness_commodity[1,chromosomes_2[i]]){
      if(W>fitness_commodity[2,chromosomes_2[i]]){
        W <- W - fitness_commodity[2,chromosomes_2[i]]
        if(L_Max<fitness_commodity[1,chromosomes_2[i]]){
          L_Max <- fitness_commodity[1,chromosomes_2[i]]
        }
      }else{
        L <- L-L_Max
        L_Max <- 0
        W <- warehouse[2]
      }
    }else{
      break
    }
    i <- i+1
    if(i>10){
      break
    }
  }
  for(j in c(1:(i-1))){
    V<-V+fitness_commodity[1,chromosomes_2[j]]*fitness_commodity[2,chromosomes_2[j]]*fitness_commodity[3,chromosomes_2[j]]
  }
  fitness_value[2] <<- V
  
  
  L <- warehouse[1]
  W <- warehouse[2]
  L_Max <- 0
  i<-1
  V<-0
  repeat {
    if(L>fitness_commodity[1,chromosomes_3[i]]){
      if(W>fitness_commodity[2,chromosomes_3[i]]){
        W <- W - fitness_commodity[2,chromosomes_3[i]]
        if(L_Max<fitness_commodity[1,chromosomes_3[i]]){
          L_Max <- fitness_commodity[1,chromosomes_3[i]]
        }
      }else{
        L <- L-L_Max
        L_Max <- 0
        W <- warehouse[2]
      }
    }else{
      break
    }
    i <- i+1
    if(i>10){
      break
    }
  }
  for(j in c(1:(i-1))){
    V<-V+fitness_commodity[1,chromosomes_3[j]]*fitness_commodity[2,chromosomes_3[j]]*fitness_commodity[3,chromosomes_3[j]]
  }
  fitness_value[3] <<- V
  
  
  
  L <- warehouse[1]
  W <- warehouse[2]
  L_Max <- 0
  i<-1
  V<-0
  repeat {
    if(L>fitness_commodity[1,chromosomes_4[i]]){
      if(W>fitness_commodity[2,chromosomes_4[i]]){
        W <- W - fitness_commodity[2,chromosomes_4[i]]
        if(L_Max<fitness_commodity[1,chromosomes_4[i]]){
          L_Max <- fitness_commodity[1,chromosomes_4[i]]
        }
      }else{
        L <- L-L_Max
        L_Max <- 0
        W <- warehouse[2]
      }
    }else{
      break
    }
    i <- i+1
    if(i>10){
      break
    }
  }
  for(j in c(1:(i-1))){
    V<-V+fitness_commodity[1,chromosomes_4[j]]*fitness_commodity[2,chromosomes_4[j]]*fitness_commodity[3,chromosomes_4[j]]
  }
  fitness_value[4] <<- V
  
  
  L <- warehouse[1]
  W <- warehouse[2]
  L_Max <- 0
  i<-1
  V<-0
  repeat {
    if(L>fitness_commodity[1,chromosomes_5[i]]){
      if(W>fitness_commodity[2,chromosomes_5[i]]){
        W <- W - fitness_commodity[2,chromosomes_5[i]]
        if(L_Max<fitness_commodity[1,chromosomes_5[i]]){
          L_Max <- fitness_commodity[1,chromosomes_5[i]]
        }
      }else{
        L <- L-L_Max
        L_Max <- 0
        W <- warehouse[2]
      }
    }else{
      break
    }
    i <- i+1
    if(i>10){
      break
    }
  }
  for(j in c(1:(i-1))){
    V<-V+fitness_commodity[1,chromosomes_5[j]]*fitness_commodity[2,chromosomes_5[j]]*fitness_commodity[3,chromosomes_5[j]]
  }
  fitness_value[5] <<- V
  
}


#菁英政策
elite <- function(){
  chromosomes_5 <<- select_1
}


#選種(選擇適應值最佳的兩個染色體，不選擇兩個同時最高的)
select_c <- function(){
  
  r <- order(fitness_value)
  
  if(r[5]==1){select_1<<-chromosomes_1}
  if(r[5]==2){select_1<<-chromosomes_2}
  if(r[5]==3){select_1<<-chromosomes_3}
  if(r[5]==4){select_1<<-chromosomes_4}
  if(r[5]==5){select_1<<-chromosomes_5}
  if(r[4]==1){select_2<<-chromosomes_1}
  if(r[4]==2){select_2<<-chromosomes_2}
  if(r[4]==3){select_2<<-chromosomes_3}
  if(r[4]==4){select_2<<-chromosomes_4}
  if(r[4]==5){select_2<<-chromosomes_5}
  if(fitness_value[r[5]]==fitness_value[r[4]]){
    if(r[1]==1){select_2<<-chromosomes_1}
    if(r[1]==2){select_2<<-chromosomes_2}
    if(r[1]==3){select_2<<-chromosomes_3}
    if(r[1]==4){select_2<<-chromosomes_4}
    if(r[1]==5){select_2<<-chromosomes_5}
  }
}


#主程式
main <- function(){
  
  initialization()#初始化
  for(i in c(1:frequency)){#執行n次
    
    fitness()#計算fitness_value
    select_c()#選種
    crossover()#交配
    mutation()#突變
    elite()#菁英政策
    
  }
  fitness()#計算最後一次fitness_value(適應值)
  print(fitness_value[5])#輸出這次計算時最好結果的適應值
  print(chromosomes_5)#輸出這次計算時最好結果的放入順序
  
}