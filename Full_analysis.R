# Section0: Libraries ----------------------------------------

Packages = c("dplyr", "stringr", "readxl", "tidyr", "ggplot2",
             "reshape2", "corrplot", "viridis", "grid", "gridExtra", "gridGraphics",
             "janitor", "RColorBrewer","ggbeeswarm", "ggpubr", "vegan", "ggcorrplot", "randomcoloR",
             "scales", "networkD3", "gridGraphics", "grid", "plyr", "tibble","ggthemes", "circlize",
             "writexl", "mdthemes", "progress", "ggrepel", "ComplexHeatmap")

lapply(Packages, library, character.only = TRUE)

# Section1: Functions -----------------------------------------------------

# --- Merge PCR to merge replicate together ---
merge_PCR = function(data,min.rep=2,CPM=TRUE,show.corr=FALSE, cor.filter = 0.6, Exclude.cfDNA = TRUE){
  #Code to generate filtered data by replicate
  #min.rep: is the minimal number of replicates
  #CPM: normalize to counts per million (CPM) prior and after merging
  #note:remove samples with zero counts before running this
  samples <- substr(names(data),1,nchar(names(data))-1) 
  merged = data.frame(row.names = row.names(data))
  #Progress bar for processing samples:
  pb <- progress_bar$new(
    format = " Merging PCR replicate: [:bar] :percent",
    total = length(unique(samples)), clear = FALSE, width= 60)
  for(s in unique(samples)){
    pb$tick()
    d = data[,grep(s,names(data))]
    d = as.data.frame(d)
    if(CPM==TRUE & ncol(d)>=2 & nrow(d)>0) d=sweep(d,2,colSums(d)/1e6,`/`)
    if(ncol(d)>=2 & show.corr==TRUE & nrow(d)>1) {
      corr=cor(d,use="pair",method="p")
      #print(corr[1,2])
      diag(corr)=NA #replace diago
      if(Exclude.cfDNA==TRUE & sum(str_detect(names(d), "cfDNA"))<1){#bypass correlation for cfDNA sample
        if(mean(corr,na.rm = TRUE)<cor.filter){merged=cbind(merged,0);next}}
      if(Exclude.cfDNA==FALSE){
        if(mean(corr,na.rm = TRUE)<cor.filter){merged=cbind(merged,0);next}}
    }
    
    if(ncol(d)<=1){merged=cbind(merged,0)}
    else {merged=cbind(merged,(rowSums(d>0,na.rm = TRUE)>=min.rep)*rowSums(d,na.rm = TRUE))}
  }
  if(CPM==TRUE) merged=sweep(merged,2,colSums(merged)/1e6,`/`)
  names(merged)=unique(samples)
  merged = merged[,!is.na(colSums(merged))]
  return(merged)
}

# --- Merge Tumour samples ---
merge_tumor = function(data,min.rep=1,CPM=TRUE, Remove_NA_col=FALSE){
  tum.pieces = names(select(data,contains("Tum")))
  tum.samples <- paste0(apply(str_split_fixed(tum.pieces,"_",3)[,1:2],1, paste,collapse="_"),"_Tum")
  tum.merged = data.frame(row.names = row.names(data))
  pb <- progress_bar$new(
    format = " Merging tumour pieces as a single tumour: [:bar] :percent",
    total = length(unique(tum.samples)), clear = FALSE, width= 60)
  
  for(s in unique(tum.samples)){
    pb$tick()
    d = data[,grep(s,names(data))]
    if(Remove_NA_col==TRUE){d=d[!is.na(colSums(d))]}
    if(CPM==TRUE) d=sweep(d,2,colSums(d)/1e6,`/`)
    if(ncol(d)<=1){tum.merged=cbind(tum.merged,0)}
    else {tum.merged=cbind(tum.merged,(rowSums(d>0)>=min.rep)*rowSums(d))}
  }
  names(tum.merged)=unique(tum.samples)
  #tum.merged=tum.merged[,colSums(tum.merged)>0]
  if(CPM==TRUE) tum.merged=sweep(tum.merged,2,colSums(tum.merged)/1e6,`/`)
  
  ##now merge this with the non.tumor samples
  return(merge(tum.merged,select(data,!contains("Tum")),by="row.names",all=TRUE)%>%tibble::column_to_rownames(var="Row.names") )
}

#--- Merge experiment or mouse as single sample by row or col
merge_exp = function(data,grp=c("exp","mouse"),by = c("col", "row"), min.rep=1,CPM=TRUE){
  exp = str_split_fixed(names(data),"_",3)[,1]
  mice = str_split_fixed(names(data),"_",3)[,2]
  mice = mice[nchar(mice)>0]
  mice = paste(exp, mice, sep = "_")
  
  df.merged = data.frame(ID= 1:nrow(data), row.names = row.names(data))
  exp.merged = data.frame()
  if(by == "row") {df1 = data.frame()}
  if(by == "col") {df1 <- as.data.frame(matrix(NA, nrow = nrow(data)))}
  
  if(grp == "exp"){
    for(s in unique(exp)){
      print(s)
      d = as.data.frame(data[,grepl(paste0(s,"_"),names(data))])
      
      if( by == "row"){
        if(ncol(d)<=1){next}
        else{
          if(CPM==TRUE) d=sweep(d,2,colSums(d)/1e6,`/`);
          d = d[,colSums(is.na(d))<nrow(d)] #remove NA column
          df.merged$n = rowSums(d>0)
          df.merged$exp = s
          df1 = rbind(df1, df.merged)
          #df.merged=cbind(df.merged,(rowSums(d>0)>=min.rep)*rowSums(d))
        }
      }
      if( by == "col"){
        if(ncol(d)<=1){df1=cbind(df1,0)}
        else{
          if(CPM==TRUE) d=sweep(d,2,colSums(d)/1e6,`/`);
          d = d[,colSums(is.na(d))<nrow(d)] #remove NA column
          df.merged$n = rowSums(d>0)
          df.merged$exp = s
          df1 = cbind(df1, df.merged$n)
        }
      }
    }
    if (by == "col"){
      df1 = df1[,-1]
      names(df1)= unique(exp)}
  }
  if(grp == "mouse"){
    for(s in unique(mice)){
      print(s)
      d = as.data.frame(data[,grepl(paste0(s,"_"),names(data))])
      
      if( by == "row"){
        if(ncol(d)<=1){next}
        else{
          if(CPM==TRUE) d=sweep(d,2,colSums(d)/1e6,`/`);
          d = d[,colSums(is.na(d))<nrow(d)] #remove NA column
          df.merged$n = rowSums(d>0)
          df.merged$exp = s
          df1 = rbind(df1, df.merged)
          #df.merged=cbind(df.merged,(rowSums(d>0)>=min.rep)*rowSums(d))
        }
      }
      if( by == "col"){
        if(ncol(d)<=1){df1=cbind(df1,0)}
        else{
          if(CPM==TRUE) d=sweep(d,2,colSums(d)/1e6,`/`);
          d = d[,colSums(is.na(d))<nrow(d)] #remove NA column
          df.merged$n = rowSums(d>0)
          df.merged$exp = s
          df1 = cbind(df1, df.merged$n)
        }
      }
    }
    if (by == "col"){ df1 = df1[,-1]
    names(df1)= unique(mice)}
  }
  if(by == "col") {df1=df1[,colSums(df1)>0]}
  return(df1)
}

#--- Function to filter data per mouse by presence of barcode in primary tumour:
Filter_by_Tum = function(data, CPM=TRUE, filter=TRUE){
  exp = str_split_fixed(names(data),"_",3)[,1]
  mice = str_split_fixed(names(data),"_",3)[,2]
  mice = mice[nchar(mice)>0]
  mice = paste(exp, mice, sep = "_")
  pb <- progress_bar$new(
    format = " [:bar] :percent",
    total = length(unique(mice)), clear = FALSE, width= 60)
  df.merged = data.frame(row.names = row.names(data))
  for(s in unique(mice)){
    #print(s)
    pb$tick()
    d = as.data.frame(data[,grepl(paste0(s,"_"),names(data)),drop=F])
    DA = which(d[,grep("Tum", names(d)), drop=F]>0)
    if(filter==TRUE) d[-DA,]=0
    if(CPM==TRUE) d=sweep(d,2,colSums(d)/1e6,`/`)
    df.merged = cbind(df.merged,d)
  }
  return(df.merged)
}

#--- Merge samples I (ex:Lung LungI) as a single sample
merge_sampleI = function(data, CPM=TRUE){
  #data=df.PCR.Tum.Filt
  #CPM=TRUE
  Mouse.nb = str_split(names(data), "_", simplify = T)[,2]
  Mega.merged = data.frame(row.names = row.names(data))
  pb <- progress_bar$new(
    format = " [:bar] :percent",
    total = length(unique(Mouse.nb)), clear = FALSE, width= 60)
  for (i in unique(Mouse.nb)) {
    pb$tick()
    TEST.name = paste("_", i, "_", sep="")
    Ms.x = data[,grep(TEST.name, names(data)), drop=F]
    id.name <- paste(str_split(names(Ms.x[1]), "_", simplify = T)[,c(1,2)], collapse = "_")
    colnames(Ms.x) <- str_split(names(Ms.x), "_", simplify = T)[,3]
    names(Ms.x)[which(endsWith(names(Ms.x), 'I'))] = str_sub(names(Ms.x)[which(endsWith(names(Ms.x), 'I'))],1,nchar(names(Ms.x)[which(endsWith(names(Ms.x), 'I'))])-1)
    
    samples <- names(Ms.x)
    merged = data.frame(row.names = row.names(Ms.x))
    for(s in unique(samples)){
      d = Ms.x[,grep(s,names(Ms.x)), drop=F]
      if(ncol(d)==1){merged=cbind(merged,d);next}
      d = d[,!is.na(colSums(d)), drop=F]
      if(ncol(d)==1){merged=cbind(merged,d);next}
      if(ncol(d)==2){d=sweep(d,2,colSums(d)/1e6,`/`)}
      merged=cbind(merged,(rowSums(d>0)*rowSums(d)))
    }
    if(CPM==TRUE) merged=sweep(merged,2,colSums(merged)/1e6,`/`)
    names(merged) = paste(id.name,unique(samples), sep = "_")
    Mega.merged = cbind(Mega.merged,merged)
  }
  return(Mega.merged)
}

#--- Stacked histogram function: needs a dataframe with one column ID
Stacked_histo = function(data, PCT=TRUE, angle.x=0, my_col = "Default"){
  if(PCT==TRUE) {data= data/10000}
  if(sum(names(data)=="ID")==0){data$ID = 1:nrow(data)}
  
  combo_color = c(brewer.pal(11, "Spectral"), brewer.pal(11, "RdYlGn"), brewer.pal(11, "RdBu"), brewer.pal(11, "PuOr"), brewer.pal(11, "PiYG"))
  getPalette = colorRampPalette(combo_color)
  if(my_col=="Default"){COLOR = sample(getPalette(2608), 2608)} else {
    COLOR = my_col
  }
  
  data$COLOR = COLOR
  COLOR.to.save = COLOR
  data = data[rowSums(data[1:(ncol(data)-2)])>0,] #remove row w/o bc (ignoring last 2 col (ID, colours))
  colnames(data) =gsub("^[^_]*_","",names(data))
  
  COLOR= data$COLOR
  data = data[,-ncol(data)]
  
  mA = melt(data, id=c("ID")) 
  mA$ID = as.factor(mA$ID) #transform barcode number in factor
  
  mA = mA %>% group_by(variable) %>% mutate(pos = sum(value)-(cumsum(value)-0.5*value))
  #Add new colunm in mA table for position of each barcode (to later on be labelled on graph)
  
  #colourCount = length(unique(mA$ID))
  
  
  #write.table(COLOR, file=paste("COLOR_exp51.txt"), sep="\t", row.names = TRUE, col.names = TRUE)
  
  insert_minor <- function(major_labs, n_minor) {labs <- 
    c( sapply( major_labs, function(x) c(x, rep("", n_minor) ) ) )
  labs[1:(length(labs)-n_minor)]}
  
  #PLot:
  p123 = ggplot(data =mA, aes(y = value, x = variable, fill = ID, label= ID)) + 
    geom_bar(stat="identity",alpha=0.9) + #pos="fill"
    #geom_text(data = subset(mA, value>5), aes(y = pos, label = ID), size = 3) +
    theme_bw()+
    theme(legend.position = "none",
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line.y = element_line(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(color = "grey10", face="bold",size=10,angle = 0),
          axis.text.x = element_text(color = "black", face="bold",size=10,angle = angle.x,
                                     vjust = 0.3),
          #hjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_text(color = "grey10", size = 15))+
    xlab("")+ylab("Number of reads(%)")+ ggtitle("")+
    scale_fill_manual(values = COLOR )+
    scale_y_continuous(breaks = seq(0,100, 10),
                       labels = insert_minor (seq(0, 100, by=20),1),
                       limits = c(0,100.001),
                       expand = c(0,0))
  my_list = list(Stack.plot = p123, Colour= COLOR.to.save)
  return(my_list)
}
Bubble_plot = function(data, PCT=TRUE, Y="avg", data2, angle.x=0, my_col=FALSE){
  
  if(PCT==TRUE){data = data/10000} #CPM to percent
  if(Y=="avg"){data$Y = rowSums(data)/ncol(data)} else if (Y=="Random"){#Average all sample for Y order or Random order
    data$Y= sample(length(data[,1]), length(data[,1]))
  } else if(Y=="input") {data$Y = data2 #manual input of previously saved Y
  } else{data$Y = data[,1]}#Take first column as Y ref = Tum collapsed
  colnames(data) =gsub("^[^_]*_","",names(data))
  data$ID = as.numeric(rownames(data))
  
  data1 = data[rowSums(data[1:(ncol(data)-2)])>0,1:(ncol(data)-2)] #remove empty rows
  data1$Y = data[rownames(data1),"Y"]
  data1$ID = data[rownames(data1),"ID"]
  
  mG =melt(data1, id=c("ID","Y"))
  mG$ID = as.factor(mG$ID) 
  mG[mG==0] = NA 
  
  Pallett = unname(distinctColorPalette(500))
  Max.color = sample(Pallett, length(unique(data$ID)),replace = TRUE)
  
  if(my_col==FALSE){col.bbp = Max.color[data1$ID]} else {
    col.bbp = my_col[data1$ID]
  }
  
  Plot.bbp = ggplot()+
    #geom_point(aes(x=variable,y=Y,color=ID,size=value), data=subset(mG,variable=!"X25k_P0_1M", stat="identity"),alpha=0.5)+
    geom_quasirandom(data=mG,aes(x=variable,y=Y,color=ID,size=value),alpha=0.8, width = 0.2)+
    scale_size_continuous(range = c(1.5,15))+
    #geom_hline(aes(yintercept = as.vector(1:10 %o% 10^(1:-7))),color="gray", size=0.2,alpha=0.3)+
    scale_color_manual(values = col.bbp,guide="none")+ #breaks=length(unique(mG$ID)) removed
    xlab("")+ylab("Frequency in Tumours")+ ggtitle("")+labs(size="Frq %")+
    #scale_y_continuous(trans = 'log10',limits=c(1e-5,100),
    #                  breaks=c(100,10,1,1e-1,1e-2,1e-3,1e-4,1e-5),
    #                 labels=c("100","10","1","0.1","0.01","0.001","0.0001","0.00001"))+
    theme_minimal()+
    theme( panel.grid.major.x = element_blank() ,
           panel.grid.minor.y = element_blank(),
           axis.text.x = element_text(angle = angle.x))
  
  if(Y=="avg"){Plot.bbp = Plot.bbp +
    scale_y_continuous(trans = 'log10',limits=c(1e-5,100),
                       breaks=c(100,10,1,1e-1,1e-2,1e-3,1e-4,1e-5),
                       labels=c("100","10","1","0.1","0.01","0.001","0.0001","0.00001"))+
    geom_hline(aes(yintercept = as.vector(1:10 %o% 10^(1:-7))),color="gray", size=0.2,alpha=0.3)} else if(Y=="Random"){
      
      Plot.bbp = Plot.bbp +
        scale_y_continuous(limits = c(0,length(data[,1])))}
  Plot.bbp
  
  
  my_list = list(BBP = Plot.bbp, Colour = Max.color, Y.order = data$Y)
  return(my_list)
}

biomass_fun = function(data){
  #data=Ms.x
  x = matrix(data=NA, nrow = ncol(data), ncol = ncol(data))
  colnames(x) = gsub("^[^_]*_","",names(data))
  rownames(x) = gsub("^[^_]*_","",names(data))
  
  for (i in 1:ncol(data)) {
    for (j in 1:ncol(data)) {
      #print(c(i,j))
      if(sum(data[i])==0){x[i,j] =0} else {x[i,j] = sum(data[data[i]>0,][j])}
    }
  }
  return(x)
}
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

look.samples = function(data, sample.to.look){
  data[,grep(paste(sample.to.look, collapse = "|"), names(data))]
}

# Section2: Pre - Figures analysis -----------------------------------------------------

info = read_xlsx("Info_all_mice.xlsx")
weight = read_xlsx("Info_all_mice.xlsx", sheet=2)

df = read.csv("data_barcodes.csv")

##   Tresholding: Barcodes with 10 Reads or less set to 0
df[df<=10]=0

##   Tresholding: Sample with more than 10k
df= df[,colSums(df)>10000]

##  Filtering based on mean correlation in replicate > 0.6
df.PCR = merge_PCR(df, show.corr = TRUE, cor.filter = 0.6, Exclude.cfDNA = TRUE)

df.PCR.Tum = merge_tumor(df.PCR)

#reorder per barcode id aka row number last function order row by 1,10,100,1000,...)
df.PCR.Tum$id <- as.integer(row.names(df.PCR.Tum))#Reorder by id name
df.PCR.Tum = df.PCR.Tum[order(df.PCR.Tum$id), ]#Reorder by id name
df.PCR.Tum = df.PCR.Tum %>% select(-"id") #Remove id column

## Filtering of barcode present uniquely in the matched primary tumour
df.PCR.Tum.Filt = Filter_by_Tum(df.PCR.Tum)

## Merge samples labelled "I" (ex:Lung, LungI, same sample run twice, as a single sample)
df.PCR.Tum.Filt.I = merge_sampleI(df.PCR.Tum.Filt)



# Section3: Figures -----------------------------------------------------
#----------------   FIGURE 1 ==========================  ----


Figure1 = function(data.info, data.df.Tum, data.df, data.weight, tresh = 0){
  info.1 = data.info %>% filter(Clever_cutting==TRUE)
  NA_tum = vector()
  No_Datum = vector()
  df.Fig1 = data.frame()
  df.Fig1.F = data.frame()
  df.Fig1.Ratio = data.frame()

  pb <- progress_bar$new(
    format = " [:bar] :percent",
    total = length(unique(info.1 %>% .$'Mice #')), clear = FALSE, width= 60)
  
  for (i in info.1 %>% .$'Mice #'){
    pb$tick()
    Mouse.ID = paste("_", i, "_", sep="")
    
    Ms.x = data.df.Tum[,grep(Mouse.ID, names(data.df.Tum)), drop = F]
    if(ncol(Ms.x)==0){No_Datum = c(No_Datum,i);next}
    
    id.name <- names(Ms.x[1])
    
    colnames(Ms.x) <- str_split(names(Ms.x), "_", simplify = T)[,3]
    if(is.na(sum(Ms.x[,"Tum"]))==T){NA_tum = c(NA_tum,i);next}
    
    Ms.x.Tum = cbind(data.df.Tum[,grep(paste(Mouse.ID,"Tum",sep=""), names(data.df.Tum)), drop=F],
                     data.df[,grep(paste(Mouse.ID,"Tum",sep=""), names(data.df)), drop=F])
    colnames(Ms.x.Tum) <- str_split(names(Ms.x.Tum), "_", simplify = T)[,3]
    
    col.order = c(paste("Tum",sort(as.numeric(as.character(str_remove(colnames(Ms.x.Tum), "Tum")))), sep=""),"Tum")
    
    Ms.x.Tum2 = Ms.x.Tum[,col.order]
    
    #Info and creation df with centre edge info:
    Ms.info = info.1 %>% filter(`Mice #`==str_split(id.name, "_", simplify = T)[2])
    
    Ms.x.edges = paste("Tum", c(1:(as.numeric(as.character(unlist(str_split(Ms.info$Centre,","))))[1]-1)), sep="")
    Ms.x.centre = paste("Tum", c(as.numeric(as.character(unlist(str_split(Ms.info$Centre,","))))), sep="")
    
    Ms.x.Tum.temp = Ms.x.Tum2[rowSums(Ms.x.Tum2)>0,]
    Ms.x.Tum.temp$centre = rowSums(Ms.x.Tum.temp[,grep(paste(Ms.x.centre,collapse = "|"), names(Ms.x.Tum.temp)), drop=F]/length(Ms.x.centre))
    Ms.x.Tum.temp$edges = rowSums(Ms.x.Tum.temp[,grep(paste(Ms.x.edges,collapse = "|"), names(Ms.x.Tum.temp)), drop=F]/length(Ms.x.edges))
    
    e = biomass_fun(Ms.x.Tum.temp)
    
    weight.ms = data.weight %>% filter(Mouse_ID == i)
    weight.ms = weight.ms[,colSums(is.na(weight.ms)) == 0]
    
    if(nrow(weight.ms)==0){weight_match = rep(NA, ncol(Ms.x.Tum.temp))} else{
      weight_match = as.character(as.vector(weight.ms[,match(names(Ms.x.Tum.temp),names(weight.ms))[!is.na(match(names(Ms.x.Tum.temp),names(weight.ms)))]]))
      weight_add_NA = rep(NA, ncol(Ms.x.Tum.temp) - length(match(names(Ms.x.Tum.temp),names(weight.ms))[!is.na(match(names(Ms.x.Tum.temp),names(weight.ms)))]))
      weight_match = c(weight_match,weight_add_NA)
    }
    
    Ms.x.Tum.df = data.frame(nb_bc = colSums(Ms.x.Tum.temp>0),
                             pct_bc = colSums(Ms.x.Tum.temp>0)*100/colSums(Ms.x.Tum.temp[,"Tum", drop=F]>0),
                             Names= colnames(Ms.x.Tum.temp),
                             Shanon=vegan::diversity(t(Ms.x.Tum.temp)),
                             Biom_Tum = e[,"Tum"],
                             weight = weight_match,
                             Mouse = i,
                             Model = Ms.info$Cells,
                             Exp = Ms.info$Exp,
                             MI=Ms.info$MI,
                             Nb_cells = Ms.info$`Number of cells injected`)
    
    Ms.x.Tum.df = Ms.x.Tum.df %>% mutate(Loca = case_when(rownames(Ms.x.Tum.df) %in% Ms.x.centre ~"Centre",
                                                          rownames(Ms.x.Tum.df) %in% Ms.x.edges ~"Edges",
                                                          rownames(Ms.x.Tum.df) %in% "centre" ~ "P_centre",
                                                          rownames(Ms.x.Tum.df) %in% "edges" ~ "P_edges",
                                                          TRUE ~"TUM"  ))
    df.Fig1 = rbind(df.Fig1,Ms.x.Tum.df)
    
    Ms.x.Tum.temp = Ms.x.Tum.temp %>% mutate(Unique.loc = case_when(Ms.x.Tum.temp[,"centre"] >tresh & Ms.x.Tum.temp[,"edges"]>tresh ~"Both",
                                                                    Ms.x.Tum.temp[,"centre"] >tresh & Ms.x.Tum.temp[,"edges"]<=tresh ~"Core",
                                                                    Ms.x.Tum.temp[,"centre"] <=tresh & Ms.x.Tum.temp[,"edges"]>tresh ~ "Peri",
                                                                    TRUE ~"NA"  ))
    Ms.x.Tum.temp = Ms.x.Tum.temp %>% mutate(Ratio.Core = centre/edges)
    
    a = Ms.x.Tum.temp %>% filter(Unique.loc=="Core") %>% nrow()
    b = Ms.x.Tum.temp %>% filter(Unique.loc=="Peri") %>% nrow()
    c = Ms.x.Tum.temp %>% filter(Unique.loc=="Both") %>% nrow()
    
    d = mean(Ms.x.Tum.temp$Ratio.Core[is.finite(Ms.x.Tum.temp$Ratio.Core)])
    
    e = median(Ms.x.Tum.temp$Ratio.Core[is.finite(Ms.x.Tum.temp$Ratio.Core)])
    
    Ms.x.Tum.df.Loca = data.frame(Loca = c("Core","Peri","Both"),
                                  bc_nb = c(a,b,c),
                                  Ratio.Core.mean = d,
                                  Ratio.Core.median =e,
                                  Mouse = i,
                                  Model = Ms.info$Cells,
                                  Exp = Ms.info$Exp,
                                  MI=Ms.info$MI)
    Ms.x.Tum.df.Loca = Ms.x.Tum.df.Loca %>% mutate(bc_pct = bc_nb*100/sum(bc_nb))
    
    df.Fig1.F = rbind(df.Fig1.F, Ms.x.Tum.df.Loca)
    
    Ms.x.Tum.df.Ratio = data.frame(centre = Ms.x.Tum.temp$centre,
                                   edges = Ms.x.Tum.temp$edges,
                                   Unique.loc = Ms.x.Tum.temp$Unique.loc,
                                   Mouse = i,
                                   Model = Ms.info$Cells,
                                   Exp = Ms.info$Exp,
                                   MI=Ms.info$MI)
    
    df.Fig1.Ratio = rbind(df.Fig1.Ratio, Ms.x.Tum.df.Ratio)
  }
  
  df.Fig1$weight = str_replace_all(df.Fig1$weight, "\\*", "")
  
  return(list(df.Fig1,df.Fig1.F,df.Fig1.Ratio))
}

## --- [FIGURE 1D] Plot %bc core vs peri no cln splitting by model ----
df.Fig1 = Figure1(info, df.PCR.Tum, df.PCR, weight)[[1]]

## Data frame for figure, clone splitting experiment excluded (Figure supp5e-f):
df.Fig1.clean = df.Fig1 %>% filter(Loca == "Edges" | Loca=="Centre") %>% filter(!Exp %in% c("80A","80B"))
df.Fig1.clean = df.Fig1.clean %>% group_by(Model) %>% dplyr::mutate(count = n_distinct(Mouse)) %>% ungroup() %>% dplyr::mutate(Model2 = paste(.$Model, " (n=", .$count, ")", sep=""))

df.Fig1.clean = df.Fig1.clean %>% group_by(Model, Loca) %>% 
  dplyr::mutate(mean = mean(pct_bc),sd= sd(pct_bc),n = n(),SE = sd(pct_bc)/sqrt(n())) %>% ungroup()

anno_df4 = compare_means(pct_bc ~ Loca, group.by = "Model2", data = df.Fig1.clean, p.adjust.method = "holm", method = "t.test") %>%
  mutate(y_pos = 105, p.adj = format.pval(p.adj, digits = 2))

n.Ms.list = list()

n.Fig1.d = df.Fig1.clean %>% select(Model, Exp, Mouse) %>% group_by(Model) %>% 
  dplyr::mutate(count.Ms = n_distinct(Mouse), count.Exp = n_distinct(Exp)) %>%
  ungroup() %>% select(Model, count.Ms, count.Exp) %>% unique()

n.Ms.list[[1]] = n.Fig1.d
names(n.Ms.list)[1] = "Figure_1d"

df.Fig1.clean %>% group_by(Model,Loca) %>% dplyr::summarise( mean = mean(pct_bc))

ggplot(data=df.Fig1.clean,aes(x=Loca, y=pct_bc,))+
  geom_bar(data=unique(df.Fig1.clean[,c("Loca","pct_bc","Model2")]),aes(fill=Loca),position = "dodge", stat = "summary", fun = "mean")+
  geom_errorbar(data=unique(df.Fig1.clean[,c("Loca","pct_bc","Model2","mean", "sd")]),aes(ymin = mean - sd, ymax = mean + sd), width=0.2)+
  geom_quasirandom(colour="blue",width = 0.3)+
  facet_wrap(~Model2, nrow = 1)+
  geom_signif(data=anno_df4,aes(xmin = group1, xmax = group2, annotations = paste(p.signif, " (", p.format,")", sep=""), y_position = y_pos),textsize=3,manual= TRUE)+
  scale_y_continuous(limits = c(0,106), breaks = c(0,20,40,60,80,100),minor_breaks = seq(0 , 100, 10))+
  labs(title = "Primary tumour Core vs Periphery (Clone splt Exp removed AND Core pieces NOT pooled) MI removed")+
  scale_color_brewer(palette = "Dark2",direction = -1)+
  scale_fill_manual(values = c("#f26522", "#662d91"))+
  theme_classic()+
  theme(strip.background = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none")+
  xlab("")+ylab("Percentage of barcodes in Tumour (%)")+
  scale_x_discrete(labels=c("Periphery", "Core"), limit=c("Edges","Centre"))#5x13

## --- [FIGURE 1C] Plot %bc per pieces Ms82 example + bubble ----

Example = merge_PCR(look.samples(df,"_27_Tum"), show.corr = TRUE, cor.filter = 0.1, Exclude.cfDNA = TRUE)

Figure1(info %>% filter(`Mice #` == 27), df.PCR.Tum, Example, weight)[[1]] %>% filter(Mouse=="27" & Loca%in% c("Edges", "Centre")) %>% ggplot()+
  geom_col(aes(x=Names,y=pct_bc, fill=Loca))+
  scale_y_continuous(limits = c(0,100))+
  scale_fill_manual(values = c("#f26522", "#662d91"))+
  theme_classic()+
  xlab("")+ylab("Percentage of barcodes in tumour pieces (%)") #4x6

Example = merge_PCR(look.samples(df,"_27_Tum"), show.corr = TRUE, cor.filter = 0.1, Exclude.cfDNA = TRUE)
Example = Example[-which(rowSums(Example>0)==0),]
colnames(Example) = str_split(colnames(Example), "_", simplify = T)[,3]
labl = data.frame(bc=colSums(Example>0), y=102)

Example$Y = 1:nrow(Example)
Example$ID = as.numeric(rownames(Example))

mG =melt(Example, id=c("ID","Y"))
mG$ID = as.factor(mG$ID)
mG[mG==0] = NA

ggplot()+
  geom_quasirandom(data=mG,aes(x=variable,y=Y,color=ID,size=value),alpha=0.8, width = 0.07)+
  scale_size_continuous(range = c(1,15))+
  geom_text(data=labl,aes(x=rownames(labl), y=y, label=bc))+
  scale_color_manual(values = unname(distinctColorPalette(length(unique(mG$ID)))),guide="none")+
  xlab("")+ylab("Barcode ID")+ ggtitle("")+labs(size="Frq %")+
  theme_classic()#8x6


## --- [FIGURE 1E] Plot Barcode location in tum:   ----
df.Fig1.E = Figure1(info, df.PCR.Tum, df.PCR, weight)[[2]]

df.Fig1.E = df.Fig1.E %>% filter(!Exp %in% c("80A", "80B")) %>% # Clone splitting experiment removed
  group_by(Model) %>%
  dplyr::mutate(count = n_distinct(Mouse)) %>% ungroup() %>% dplyr::mutate(Model2 = paste(.$Model, " (n=", .$count, ")", sep=""))

df.Fig1.E = df.Fig1.E %>% group_by(Model, Loca) %>% 
  dplyr::mutate(mean = mean(bc_pct),sd= sd(bc_pct),n = n(),SE = sd(bc_pct)/sqrt(n())) %>% ungroup()

n.Fig1.e = df.Fig1.E %>% select(Model, Exp, Mouse) %>% group_by(Model) %>% 
  dplyr::mutate(count.Ms = n_distinct(Mouse), count.Exp = n_distinct(Exp)) %>%
  ungroup() %>% select(Model, count.Ms, count.Exp) %>% unique()

n.Ms.list[[length(n.Ms.list)+1]] = n.Fig1.e
names(n.Ms.list)[length(n.Ms.list)] = "Figure_1e"

# Anova with follow up Tukey
compare_means(bc_pct ~ Loca,
              group.by = "Model",
              data = df.Fig1.E,
              p.adjust.method = "holm",
              method = "anova") %>%
  mutate(p.adj = format.pval(p.adj, digits = 2))
# Follow up Tukey test
TukeyHSD(aov(bc_pct~Loca, data = df.Fig1.E %>% filter(Model == "MDA-231")))
TukeyHSD(aov(bc_pct~Loca, data = df.Fig1.E %>% filter(Model == "MDA-468")))
TukeyHSD(aov(bc_pct~Loca, data = df.Fig1.E %>% filter(Model == "PDX-1432C")))
TukeyHSD(aov(bc_pct~Loca, data = df.Fig1.E %>% filter(Model == "PDX-434")))
TukeyHSD(aov(bc_pct~Loca, data = df.Fig1.E %>% filter(Model == "PDX-T412")))

anno_df4 = compare_means(bc_pct ~ Loca, group.by = "Model2", data = df.Fig1.E, p.adjust.method = "holm", method = "t.test") %>%
  mutate(p.adj = format.pval(p.adj, digits = 2))
anno_df4$y_pos =rep(c(75,100,110),6)

df.Fig1.E %>%
  ggplot(aes(x=Loca, y=bc_pct))+
  geom_bar(aes(fill=Loca),position = "dodge", stat = "summary", fun = "mean")+
  geom_errorbar(aes(ymax = mean + sd, ymin = ifelse(mean - sd<0, 0,mean - sd)), width=0.2)+
  geom_quasirandom(colour="darkblue",width = 0.1)+
  geom_signif(data=anno_df4,aes(xmin = group1,
                                xmax = group2,
                                annotations = paste(p.signif, " (", p.format,")", sep=""),
                                y_position = y_pos),
              textsize=3,manual= TRUE, )+
  facet_wrap(~Model2, nrow = 2, ncol=3)+
  labs(title = "Barcode spatial location in tumour by compartment, Core, Periphery or Both")+
  scale_fill_manual(values = c("grey", "#f26522", "#662d91"))+
  theme_classic()+
  theme(strip.background = element_blank(), panel.grid.major.x = element_blank())+
  xlab("")+ ylab("Percentage of barcode")+
  scale_y_continuous(limits = c(0,120), breaks = c(0,20,40,60,80,100),minor_breaks = seq(0 , 100, 10))+
  scale_x_discrete(limit= c("Peri", "Core", "Both"), labels=c("Periphery", "Core", "Both"),) #6x7

## --- [FIGURE SUPP 1A] -----

df.Fig1 %>% filter(Loca == "TUM") %>% filter( Mouse %in% unique(df.Fig1.clean$Mouse)) %>%
  group_by(Model) %>% 
  dplyr::mutate(mean_sha = mean(nb_bc),sd_sha= sd(nb_bc),n_sha = n(),SE_sha = sd(nb_bc)/sqrt(n())) %>% ungroup() %>%
  ggplot(aes(x=Model, y=nb_bc))+
  geom_bar(position = "dodge", stat = "summary", fun = "mean")+
  geom_quasirandom(colour="blue",width = 0.2)+
  geom_errorbar(aes(ymin = mean_sha - sd_sha, ymax = mean_sha + sd_sha), width=0.2)+
  theme_classic()+
  ylab("Number of barcodes") #4.5x5

df.Fig1 %>% filter(Loca == "TUM") %>% filter( Mouse %in% unique(df.Fig1.clean$Mouse)) %>%
  group_by(Model) %>% dplyr::summarise(mean = mean(nb_bc), median = median(nb_bc), n = dplyr::n())

## --- [FIGURE SUPP 1E] Plot weight and % barcode: ----

df.Fig.weight = df.Fig1.clean[!is.na(df.Fig1.clean$weight),]
df.Fig.weight$weight = as.numeric(as.character(df.Fig.weight$weight))


r2<-ddply(df.Fig.weight,.(Model2),function(x) summary(lm(x$pct_bc ~ x$weight))$r.squared)
names(r2)<-c("Model2","r2")

n.Supp.Fig1.d = df.Fig.weight %>% select(Model, Exp, Mouse) %>% group_by(Model) %>% 
  dplyr::mutate(count.Ms = n_distinct(Mouse), count.Exp = n_distinct(Exp)) %>%
  ungroup() %>% select(Model, count.Ms, count.Exp) %>% unique()

n.Ms.list[[length(n.Ms.list)+1]] = n.Supp.Fig1.d
names(n.Ms.list)[length(n.Ms.list)] = "Supp_Figure_1d"

ggplot(data=df.Fig.weight,aes(x=weight, y=pct_bc))+
  geom_point(aes(colour= Loca))+
  geom_smooth(method = "lm", se = FALSE)+
  geom_text(data=r2,aes(x=75,y=100, label = paste("R^2: ", round(r2,digit=4),sep="")),parse=T)+
  facet_wrap(~Model2,ncol=1, scale="free")+
  theme_classic()+
  theme(strip.background = element_blank(), panel.grid.major.x = element_blank())+
  labs(title = "Percentage of barcodes per pieces by weight")+
  xlab("Weight (mg)")+
  ylab("Percentage of barcodes in Tumour (%)") #12x4


## --- [FIGURE SUPP 1C] Example specific mouse bubble plot and histogram: ----

# Examples = c(72,54,62, 67, 50, 78)

Exp= "_50_" # see above for other mouse number in examples
str_split(Exp, "_", simplify=T)[,2]
Ms.Tum.exp = df.PCR[,grep(paste(Exp,"Tum",sep=""), names(df.PCR)), drop=F]
colnames(Ms.Tum.exp) <- str_split(names(Ms.Tum.exp), "_", simplify = T)[,3]
col.order = c(paste("Tum",sort(as.numeric(as.character(str_remove(colnames(Ms.Tum.exp), "Tum")))), sep=""))
Ms.Tum.exp = Ms.Tum.exp[,col.order]

p1 = Bubble_plot(Ms.Tum.exp)$BBP+labs(title = paste("Mouse",str_split(Exp, "_", simplify=T)[,2]))

p1

## --- [FIGURE SUPP 1D] Heatmap tumour ----


examples.mice = c(72,54,62,67,50,78)
pdf(file= "Heatmaps_Examples_Mice.pdf", height = 10, width = 10)
for (i in examples.mice){
  print(i)
  TEST.name = paste("_", i, "_", sep="")
  
  Ms.x = df.PCR.Tum[,grep(TEST.name, names(df.PCR.Tum))]
  id.name <- names(Ms.x[1])
  colnames(Ms.x) <- str_split(names(Ms.x), "_", simplify = T)[,3]
  Ms.x.Tum = cbind(df.PCR.Tum[,grep(paste(TEST.name,"Tum",sep=""), names(df.PCR.Tum)), drop=F],
                   df.PCR[,grep(paste(TEST.name,"Tum",sep=""), names(df.PCR)), drop=F])
  colnames(Ms.x.Tum) <- str_split(names(Ms.x.Tum), "_", simplify = T)[,3]
  col.order = c(paste("Tum",sort(as.numeric(as.character(str_remove(colnames(Ms.x.Tum), "Tum")))), sep=""),"Tum")
  Ms.x.Tum2 = Ms.x.Tum[,col.order]
  Ms.x.Tum2 = Ms.x.Tum2[rowSums(Ms.x.Tum2[1:ncol(Ms.x.Tum2)-1])>0,]
  
  Ms.info = df.Fig1%>% filter(Mouse==i) %>% filter(Names=="Tum")
  Ms.info[1,"Model2", drop=T]
  M1 = log(as.matrix(Ms.x.Tum2)+1)
  p1 = heatmap(t(M1),Rowv=T, Colv = T, col= viridis(256, direction=1), scale = "none",
               main=paste("Ms",i, Ms.info[1,"Model", drop=T], "bc:", Ms.info[1,"nb_bc", drop=T]))
  p1
}
dev.off()

#  --- [FIGURE SUPP 2C] Scatter plot with lung mets: -----
df.Fig1.Lung=data.frame()

for (i in unique(df.Fig1.clean$Mouse)) {
  print(i)
  #i=41
  if(i %in% c(41,42)){next} # Skip mice with "LungNodes" present outside lungs
  info.1 = info %>% filter(Clever_cutting==TRUE)
  Ms.info = info.1 %>% filter(`Mice #`==i)
  TEST.name = paste("_", i, "_", sep="")
  
  Ms.x = df.PCR.Tum.Filt.I[,grep(TEST.name, names(df.PCR.Tum.Filt.I)), drop=FALSE]
  colnames(Ms.x) <- str_split(names(Ms.x), "_", simplify = T)[,3]
  names(Ms.x)
  if(sum(names(Ms.x)%in% "SNLung")>0){Ms.x = subset(Ms.x, select=-SNLung)}
  if(ncol(Ms.x)==0){No_Datum = c(No_Datum,i);next}
  id.name <- names(Ms.x[1])
  
  Ms.x.Core = paste("Tum", c(as.numeric(as.character(unlist(str_split(Ms.info$Core,","))))), sep="")
  Ms.x.Tum = df.PCR[,grep(TEST.name, names(df.PCR))]
  colnames(Ms.x.Tum) <- str_split(names(Ms.x.Tum), "_", simplify = T)[,3]
  
  # Check to merge core of the tumour when centre of the tumour was 2 or more samples
  if(length(Ms.x.Core)==1 & sum(Ms.x.Core %in% "TumNA")>0){next
  } else if(length(Ms.x.Core)==1 & sum(Ms.x.Core %in% names(Ms.x.Tum)==TRUE)>0){Ms.x.Tum$Core = Ms.x.Tum[,Ms.x.Core]
  } else if(length(Ms.x.Core)==1 & sum(Ms.x.Core %in% names(Ms.x.Tum)==FALSE)>0){next
  }else {Ms.x.Tum$Core = rowSums(Ms.x.Tum[,Ms.x.Core])
  }
  
  Ms.x = Ms.x[,!is.na(colSums(Ms.x)), drop=F]
  Ms.x = cbind(Ms.x, Core=Ms.x.Tum$Core)
  Ms.x = Ms.x[rowSums(Ms.x)>0,]
  
  aa=Ms.x
  Needles = c("Tum", "Core", "N1a", "N1b", "N2a", "N2b","N3a", "N3b", "N4a", "N4b")
  a=Ms.x[,!names(Ms.x) %in% Needles, drop=F]
  cond.lung = sum(str_detect(names(a), "Lung"))
  aa$Mets=(aa$Tum>0)+(as.vector(rowSums(a)>0))==2
  cond.lung
  if(cond.lung==0){aa$is.Lung= 0} else{
    aa$is.Lung = (aa$Tum>0)+(a[,str_detect(names(a), "Lung")]>0)==2}
  
  if(sum(str_detect(names(Ms.x), "Lung"))==0){next
  }else{
    Ms.x.temp= aa[,c("Tum", "Core","Mets", "is.Lung")]
    Ms.x.temp$LUNG = Ms.x[,str_detect(names(Ms.x), "Lung")]
    Ms.x.temp$mouse = i
    Ms.x.temp$Model = Ms.info$Cells
    Ms.x.temp$Exp = Ms.info$Exp
  }
  df.Fig1.Lung = rbind(df.Fig1.Lung, Ms.x.temp)
}

df.Fig1.Lung = df.Fig1.Lung %>% group_by(Model) %>% 
  dplyr::mutate(count = n_distinct(mouse)) %>% ungroup() %>% 
  dplyr::mutate(Model2 = paste(.$Model, " (n=", .$count, ")", sep=""))

df.Fig1.Lung.n= subset(df.Fig1.Lung,LUNG>0&Tum>0)

df.Fig1.Lung.n %>% group_by(Model) %>% 
  dplyr::summarise(n = n_distinct(mouse), n.exp =n_distinct(Exp)) 

ggplot(df.Fig1.Lung)+
  geom_point(data= subset(df.Fig1.Lung,LUNG>0&Tum>0),aes(x=LUNG, y=Tum, color=Core>0))+
  geom_quasirandom(data=subset(df.Fig1.Lung,LUNG==0 & Tum>0), aes(x=10, y=Tum, color=Core>0), width=0.4, groupOnX=TRUE,alpha=0.7)+
  scale_color_manual(values = c("dimgray", "darkorange"), labels=c("No", "Yes"))+scale_x_log10()+scale_y_log10()+xlab("Frequency in LUNG")+ylab("Frequency in Tumour")+
  labs(title = "Scatter plot Tumour vs Lung", color= "in Core")+
  facet_wrap(~Model2, nrow = 1)+theme_classic()#3.5x15

#  --- [FIGURE SUPP 3A + 3B] Barcode frq in Tum by Pct pieces:----
df.fig1.Excl = data.frame()

#Loop to calculate per barcode their appearance in number of piece and % /core
for (i in unique(df.Fig1.clean$Mouse)) {
  print(i)
  test_84 = df.PCR[,grep(paste("_",i,"_", sep=""), names(df.PCR))]
  test_84_Tum = select(test_84,contains("Tum"))
  test_84_TUM = df.PCR.Tum[,grep(paste("_",i,"_", sep=""), names(df.PCR.Tum)), drop=F]
  test_84_TUM = select(test_84_TUM,contains("Tum"))
  test_84_Tum = cbind(test_84_Tum,test_84_TUM[1])
  test_84_Tum = test_84_Tum[-which(rowSums(test_84_Tum>0)==0),] #remove bc with 0 reads in all samples
  colnames(test_84_Tum) =  str_split(colnames(test_84_Tum),"_", simplify = T)[,3]
  test_84_Tum$Excl = (rowSums(test_84_Tum>0)-1)*100/(ncol(test_84_Tum)-1)
  
  #INfo and creation df with centre edge info:
  info.1 = info %>% filter(Clever_cutting==TRUE)
  Ms.info = info.1 %>% filter(`Mice #`==i)
  
  Ms.x.outer.P = paste("Tum", c(as.numeric(as.character(unlist(str_split(Ms.info$Outer.P,","))))), sep="")
  Ms.x.Peri = paste("Tum", c(as.numeric(as.character(unlist(str_split(Ms.info$Peri,","))))), sep="")
  Ms.x.outer.C = paste("Tum", c(as.numeric(as.character(unlist(str_split(Ms.info$Outer.C,","))))), sep="")
  Ms.x.Core = paste("Tum", c(as.numeric(as.character(unlist(str_split(Ms.info$Core,","))))), sep="")
  
  # Check to merge core of the tumour when centre of the tumour was 2 or more samples
  if(length(Ms.x.Core)==1 & sum(Ms.x.Core %in% "TumNA")>0){next
  } else if(length(Ms.x.Core)==1 & sum(Ms.x.Core %in% names(test_84_Tum)==TRUE)>0){test_84_Tum$Core = test_84_Tum[,Ms.x.Core]
  } else if(length(Ms.x.Core)==1 & sum(Ms.x.Core %in% names(test_84_Tum)==FALSE)>0){next
  }else {test_84_Tum$Core = rowSums(test_84_Tum[,Ms.x.Core])
  } 
  
  Ms.x.Core %in% names(test_84_Tum) ==FALSE
  
  Peri.pieces = c(Ms.x.outer.C,Ms.x.Peri,Ms.x.outer.P)
  Peri.pieces = Peri.pieces[!Peri.pieces == "TumNA"]
  Peri.pieces = Peri.pieces[which(Peri.pieces %in% names(test_84_Tum))]
  
  test_84_Tum$n.Peri = rowSums(test_84_Tum[,Peri.pieces]>0)
  test_84_Tum$pct.Peri = test_84_Tum$n.Peri*100/length(Peri.pieces)
  
  test_84_df = data.frame(Exp = Ms.info$Exp, Mouse = i, Model = Ms.info$Cells,
                          ID.bc = rownames(test_84_Tum),
                          Biom.Tum = test_84_Tum$Tum,
                          Excl.Tum = test_84_Tum$Excl,
                          Core = test_84_Tum$Core,
                          n.Peri = test_84_Tum$n.Peri,
                          pct.Peri = test_84_Tum$pct.Peri)
  
  df.fig1.Excl = rbind(df.fig1.Excl,test_84_df)
}

n.Supp.Fig3.a = df.fig1.Excl %>% select(Model, Exp, Mouse) %>% group_by(Model) %>% 
  dplyr::mutate(count.Ms = n_distinct(Mouse), count.Exp = n_distinct(Exp)) %>%
  ungroup() %>% select(Model, count.Ms, count.Exp) %>% unique()

n.Ms.list[[length(n.Ms.list)+1]] = n.Supp.Fig3.a
names(n.Ms.list)[length(n.Ms.list)] = "Supp_Figure_3a"

ggplot()+
  geom_point(data= subset(df.fig1.Excl, n.Peri>=1 & Core>0), aes(x=Excl.Tum, y=Biom.Tum/10000), colour= "gray")+
  geom_point(data= subset(df.fig1.Excl, Core==0), aes(x=Excl.Tum, y=Biom.Tum/10000), colour= "darkblue")+
  geom_point(data= subset(df.fig1.Excl, n.Peri==0), aes(x=Excl.Tum, y=Biom.Tum/10000), colour= "darkorange")+
  scale_y_log10()+
  facet_wrap(~Model, nrow=1)+
  scale_color_viridis()+
  theme_classic()+
  xlab("Percentage of presence across tumour pieces (%)")+
  ylab("Barcode frequency in Tumour (%)")+
  labs(title = "Barcode frq in Tum by Pct pieces")#4x14

df.fig1.Excl$Only.Peri = !df.fig1.Excl$Core>0

p1 = df.fig1.Excl %>% filter(Only.Peri==FALSE) %>% ggplot()+
  geom_density(aes(x=pct.Peri),fill="#f26522", alpha =0.7)+
  facet_wrap(~Model, scales = "free_y", nrow=1)+
  theme_classic()+
  xlab("Percentage of presence across tumour pieces (%)")+
  ylab("Density")

p2 =df.fig1.Excl %>% filter(Only.Peri==TRUE) %>% ggplot()+
  geom_density(aes(x=pct.Peri),fill="#662d91", alpha =0.8)+
  facet_wrap(~Model, scales = "free_y", nrow = 1)+
  theme_classic()+
  xlab("Percentage of presence across tumour pieces (%)")+
  ylab("Density")

ggarrange(p1,p2, nrow = 2) #8x14


#----------------   FIGURE 2 ==========================    ----
## --- Processing Pre-figure -----

Needles = c("Tum", "N1a", "N1b", "N2a", "N2b","N3a", "N3b", "N4a", "N4b")

No_Needles = vector()
No_corr_Needles = vector()
df.Fig2 = data.frame()
df.corr = data.frame()

for (i in unique(df.Fig1.clean$Mouse)) {
  
  print(i)
  TEST.name = paste("_", i, "_", sep="")
  
  Ms.x = df.PCR.Tum.Filt[,grep(TEST.name, names(df.PCR.Tum.Filt)), drop=FALSE]
  if(ncol(Ms.x)==0){No_Datum = c(No_Datum,i);next}
  
  id.name <- names(Ms.x[1])
  
  colnames(Ms.x) <- str_split(names(Ms.x), "_", simplify = T)[,3]
  if(is.na(sum(Ms.x[,"Tum"]))==T){NA_tum = c(NA_tum,i);next}
  
  Ms.x.Needles = df.PCR.Tum.Filt[,grep(paste(paste(TEST.name,Needles, sep=""),collapse = "|"), names(df.PCR.Tum.Filt)), drop=F]
  if(ncol(Ms.x.Needles)<=1){No_Needles = c(No_Needles,i);next}
  colnames(Ms.x.Needles) <- str_split(names(Ms.x.Needles), "_", simplify = T)[,3]
  Ms.x.Needles = Ms.x.Needles[,!is.na(colSums(Ms.x.Needles))]
  
  Ms.x.Needles.temp = Ms.x.Needles[rowSums(Ms.x.Needles)>0,]
  e = biomass_fun(Ms.x.Needles.temp)
  
  Ms.x.Needles.filt = Ms.x.Needles.temp[,2:ncol(Ms.x.Needles.temp)]
  Ms.x.Needles.filt[Ms.x.Needles.filt<10000]=0 # FILTER BARCODE UNDER 1% in needles
  Ms.x.Needles.temp.1 = cbind(Ms.x.Needles.temp[,1, drop=F],Ms.x.Needles.filt)
  
  e1 = biomass_fun(Ms.x.Needles.temp.1)
  
  
  Ms.x.Needles.filt2 = Ms.x.Needles.temp[,2:ncol(Ms.x.Needles.temp)]
  Ms.x.Needles.filt2[Ms.x.Needles.filt2<100000]=0 # FILTER BARCODE UNDER 10% in needles
  Ms.x.Needles.temp.2 = cbind(Ms.x.Needles.temp[,1, drop=F],Ms.x.Needles.filt2)
  
  e2 = biomass_fun(Ms.x.Needles.temp.2)
  
  Needle.corr = Ms.x.Needles.filt %>% select(contains("a"))
  if(ncol(Needle.corr)<=1){No_corr_Needles = c(No_corr_Needles,i);next}

  r = as.data.frame(cor(Needle.corr)) %>% rownames_to_column() %>% pivot_longer(cols = -rowname) %>% unite(name, rowname, name, sep = "-")
  
  Ms.info = df.Fig1.clean %>% filter(Mouse==i) %>% select(c("Exp", "Model", "Model2", "Nb_cells")) %>% unique() %>% mutate(Mouse =i)
  
  
  d1 = cbind(r,Ms.info)
  
  df.corr = rbind(df.corr, d1)
  
  Ms.x.Needles.df = data.frame(nb_bc = colSums(Ms.x.Needles.temp>0),
                               pct_bc = colSums(Ms.x.Needles.temp>0)*100/colSums(Ms.x.Needles.temp[,"Tum", drop=F]>0),
                               Names= colnames(Ms.x.Needles.temp),
                               Shanon=diversity(t(Ms.x.Needles.temp)),
                               Biom_Tum = e[,"Tum"],
                               Biom_Tum1 = e1[,"Tum"],
                               Biom_Tum10 = e2[,"Tum"],
                               Ms.info)
  
  df.Fig2 = rbind(df.Fig2,Ms.x.Needles.df)
}

df.corr$x = str_split(df.corr$name,"-", simplify = T)[,1]
df.corr$y = str_split(df.corr$name,"-", simplify = T)[,2]


df.Fig2D=data.frame()
for (i in unique(df.Fig1.clean$Mouse)) {
  print(i)
  info.1 = info %>% filter(Clever_cutting==TRUE)
  Ms.info = info.1 %>% filter(`Mice #`==i)
  
  TEST.name = paste("_", i, "_", sep="")
  
  Ms.x.outer.P = paste("Tum", c(as.numeric(as.character(unlist(str_split(Ms.info$Outer.P,","))))), sep="")
  Ms.x.Peri = paste("Tum", c(as.numeric(as.character(unlist(str_split(Ms.info$Peri,","))))), sep="")
  Ms.x.outer.C = paste("Tum", c(as.numeric(as.character(unlist(str_split(Ms.info$Outer.C,","))))), sep="")
  Ms.x.Core = paste("Tum", c(as.numeric(as.character(unlist(str_split(Ms.info$Core,","))))), sep="")
  
  Ms.x.Tum = df.PCR[,grep(TEST.name, names(df.PCR)), drop=FALSE]
  colnames(Ms.x.Tum) <- str_split(names(Ms.x.Tum), "_", simplify = T)[,3]
  if(length(Ms.x.Core)==1 & sum(Ms.x.Core %in% "TumNA")>0){next
  } else if(length(Ms.x.Core)==1 & sum(Ms.x.Core %in% names(Ms.x.Tum)==FALSE)>0){next 
  } else if(length(Ms.x.Core)==1 & sum(Ms.x.Core %in% names(Ms.x.Tum)==TRUE)>0){Ms.x.Tum$Core = Ms.x.Tum[,Ms.x.Core]
  }else {Ms.x.Tum$Core = rowSums(Ms.x.Tum[,Ms.x.Core])
  }
  
  Ms.x = df.PCR.Tum.Filt[,grep(TEST.name, names(df.PCR.Tum.Filt)), drop=FALSE]
  if(ncol(Ms.x)==0){No_Datum = c(No_Datum,i);next}
  
  id.name <- names(Ms.x[1])
  
  colnames(Ms.x) <- str_split(names(Ms.x), "_", simplify = T)[,3]
  if(is.na(sum(Ms.x[,"Tum"]))==T){NA_tum = c(NA_tum,i);next}
  
  Ms.x.Needles = df.PCR.Tum.Filt[,grep(paste(paste(TEST.name,Needles, sep=""),collapse = "|"), names(df.PCR.Tum.Filt)), drop=F]
  if(ncol(Ms.x.Needles)<=1){No_Needles = c(No_Needles,i);next}
  colnames(Ms.x.Needles) <- str_split(names(Ms.x.Needles), "_", simplify = T)[,3]
  Ms.x.Needles = Ms.x.Needles[,!is.na(colSums(Ms.x.Needles))]
  Ms.x.Needles = cbind(Ms.x.Needles, Core =Ms.x.Tum$Core>0)
  Ms.x.Needles.temp = Ms.x.Needles[rowSums(Ms.x.Needles)>0,]
  if(sum(colnames(Ms.x.Needles) %in% c("N1a", "N2a"))<2){next} 
  
  Ms.x.Needles.N1N2 = Ms.x.Needles.temp[,c("N1a", "N2a", "Core", "Tum")]
  Ms.x.Needles.N1N2$mouse = i
  Ms.x.Needles.N1N2$Model = Ms.info$Cells
  Ms.x.Needles.N1N2$Exp = Ms.info$Exp
  
  df.Fig2D = rbind(df.Fig2D, Ms.x.Needles.N1N2)
}

## --- [FIGURE 2B] Scatter Plot Needle1 vs Tumour  ----

n.Fig2.b = df.Fig2D %>% select(Model, Exp, mouse) %>% group_by(Model) %>% 
  dplyr::mutate(count.Ms = n_distinct(mouse), count.Exp = n_distinct(Exp)) %>%
  ungroup() %>% select(Model, count.Ms, count.Exp) %>% unique()

n.Ms.list[[length(n.Ms.list)+1]] = n.Fig2.b
names(n.Ms.list)[length(n.Ms.list)] = "Figure_2b"

ggplot()+
  geom_hline(yintercept=1)+ geom_vline(xintercept = 1)+
  geom_point(data= subset(df.Fig2D,N1a>0&Tum>0),aes(x=N1a/10000, y=Tum/10000, color=Core), alpha=0.9)+
  #geom_point(data= subset(df.Fig2D,N2a>0&Tum>0),aes(x=N2a/10000, y=Tum/10000, color=Core),shape=2, alpha=0.9)+
  geom_quasirandom(data= subset(df.Fig2D,N1a==0&Tum>0),aes(x=0.0001, y=Tum/10000, color=Core ),width=0.5,groupOnX=TRUE,alpha=0.7)+
  scale_x_log10()+
  scale_y_log10()+
  facet_wrap(~Model, nrow = 2, ncol=3)+
  theme_classic()+
  scale_color_manual(values = c("grey30", "brown2"))+
  labs(title = "Scatter plot Needle1 vs Tumour")+
  xlab("Barcode frequencies in Needle1 (%)")+ylab("Barcode frequencies in Tumour (%)")#5x9


## --- [FIGURE 2C] Example histogram figure2  -----

Examples =c(50,70,86,61)

Ne = c("Tum", "N1a", "N2a", "N3a", "N4a")
Exp= "_61_"
Ms.Tum.exp = df.PCR.Tum.Filt[,grep(paste(Exp, Ne, sep="", collapse = "|"), names(df.PCR.Tum.Filt)), drop=F]
colnames(Ms.Tum.exp) <- str_split(names(Ms.Tum.exp), "_", simplify = T)[,3]

colSums(Ms.Tum.exp>0) # Number of barcodes
Stacked_histo(Ms.Tum.exp) #4x4

## --- [FIGURE 2D] Scatter Plot Needle1 vs Needle2  ----

ggplot()+
  geom_point(data=df.Fig2D,aes(x=N1a/10000+0.0001, y=N2a/10000+0.001, fill=Core, size=Tum),colour="black",pch=21, alpha=0.5)+
  scale_size_continuous(range = c(0.8,10))+
  scale_x_log10()+
  scale_y_log10()+
  facet_wrap(~Model, nrow = 1)+
  theme_classic()

## --- [FIGURE SUPP 4A] Pearson correlation all needles ----

n.Supp.Fig4.a = df.corr %>% select(Model, Exp, Mouse) %>% group_by(Model) %>% 
  dplyr::mutate(count.Ms = n_distinct(Mouse), count.Exp = n_distinct(Exp)) %>%
  ungroup() %>% select(Model, count.Ms, count.Exp) %>% unique()

n.Ms.list[[length(n.Ms.list)+1]] = n.Supp.Fig4.a
names(n.Ms.list)[length(n.Ms.list)] = "Supp_Figure_4a"

df.corr %>% filter(value<1) %>% ggplot()+
  geom_violin(aes(x=x, y=value))+
  geom_quasirandom(aes(x=x, y=value), color="grey50", alpha=0.5)+
  facet_wrap(~Model, nrow=1)+
  scale_color_brewer(palette = "Paired")+
  theme_classic()+
  labs(x="", y= "Pearson correlation")

## --- [FIGURE SUPP 4B] Mean Pearson correlation matrix all needles ----

data = df.corr %>% group_by(x,y,name,Model2) %>% dplyr::summarise(mean_pearson=mean(value))

ggplot(data, aes(x=x,y=y,fill=mean_pearson))+
  geom_tile()+
  scale_fill_viridis(option = "D")+
  theme_classic2()+
  theme(strip.background = element_blank(),panel.grid.major.x = element_blank())+
  facet_wrap(~Model2, nrow = 1)+
  labs(x="", y= "")

## --- [FIGURE SUPP 4C] Biomass captured with threshold  with needles ----

n.Supp.Fig4.c = df.Fig2 %>% filter(Names!="Tum") %>% filter(Names %in% c("N1a", "N2a", "N3a", "N4a"))%>%
  group_by(Model, Names) %>% select(Model, Exp, Mouse) %>% group_by(Model) %>% 
  dplyr::mutate(count.Ms = n_distinct(Mouse), count.Exp = n_distinct(Exp)) %>%
  ungroup() %>% select(Model, count.Ms, count.Exp) %>% unique()

n.Ms.list[[length(n.Ms.list)+1]] = n.Supp.Fig4.c
names(n.Ms.list)[length(n.Ms.list)] = "Supp_Figure_4c"

df.Supp.Fig2_C = df.Fig2 %>% filter(Names!="Tum") %>% filter(Names %in% c("N1a", "N2a", "N3a", "N4a"))%>% 
  select("Names", "Biom_Tum", "Biom_Tum1", "Biom_Tum10", "Model") %>%
  melt(id.vars=c("Names", "Model"))

df.Supp.Fig2_C.clean = df.Supp.Fig2_C %>% group_by(Model, Names, variable) %>% 
  dplyr::mutate(mean_biom = mean(value),sd_sha= sd(value),n_sha = n(),SE_sha = sd(value)/sqrt(n())) %>% ungroup()

# Stats on Supp figure 2C
anno_df4 = compare_means(value ~ variable, 
                         group.by = c("Names", "Model"),
                         data = df.Supp.Fig2_C.clean,
                         p.adjust.method = "holm",
                         method = "t.test") %>%
  mutate(p.adj = format.pval(p.adj, digits = 2))

anno_df4 = anno_df4 %>% arrange(Model, Names)
anno_df4$y_pos = c(rep(c(1300000,1200000, 1100000), 24))
anno_df4$x = c(rep(c(rep(1,3),rep(2,3),rep(3,3),rep(4,3)),6))

group_levels <- levels(as.factor(df.Supp.Fig2_C.clean$variable))
condition_levels <- levels(as.factor(df.Supp.Fig2_C.clean$Names))

group_dodge <- seq(-0.25, 0.25, length.out = length(group_levels))

get_dodge_pos <- function(group, condition_num) {
  group_idx <- match(group, group_levels)
  condition_num + group_dodge[group_idx]
}

anno_df4 <- anno_df4 %>%
  dplyr::mutate(
    xmin = get_dodge_pos(group1, x),
    xmax = get_dodge_pos(group2, x)
  )

ggplot(df.Supp.Fig2_C.clean, aes(x=Names, y=value))+
  stat_summary(aes(fill=variable),geom = "bar", position = "dodge", fun="mean", colour="black")+
  geom_errorbar(aes(group=variable, 
                    ymin = ifelse(mean_biom - sd_sha<0, 0,mean_biom - sd_sha),
                    ymax = mean_biom + sd_sha), position = "dodge")+
  stat_pvalue_manual(anno_df4,label = "p.signif", y.position = "y_pos",
                     tip.length = 0.01, bracket.size = 0.6 )+
  facet_wrap(~Model, nrow = 1)+
  scale_fill_brewer(palette = 1, direction = -1)+
  theme_classic()+
  labs(x="", y="Tumour biomass captured")


#----   FIGURE 3    ==========================  ----
## --- Processing Pre-figure ----

v= colSums(df.PCR.Tum.Filt.I[,grep("cfDNA", names(df.PCR.Tum.Filt.I))]) #CPM
m= colSums(df.PCR.Tum.Filt.I[,grep("cfDNA", names(df.PCR.Tum.Filt.I))]>0) #number of barcode
b= names(df.PCR.Tum.Filt.I[,grep("cfDNA", names(df.PCR.Tum.Filt.I))]) #sample name
cfDNA.all = data.frame(b, v, m)

## --- Processing Second version dataframe with cfDNA ----

cfDNA.filt = cfDNA.all
cfDNA.filt = cbind(cfDNA.filt, str_split(cfDNA.filt$b, "_", simplify = T))

tr= info %>% select(`Mice #`, Cells, Exp) #Info of mice number, experiment number and cells injected
tr=as.data.frame(tr)
ft = vector()
vec.exp = vector()
# Loop to create vector with tumour model and true experiment number for each sample
for (i in 1:nrow(cfDNA.filt)) {
  print(i)
  tiu = as.vector(tr[which(tr$`Mice #`==cfDNA.filt[i,5]),2])
  ft=c(ft,tiu)
  temp.exp = as.vector(tr[which(tr$`Mice #`==cfDNA.filt[i,5]),3])
  vec.exp=c(vec.exp,temp.exp)
}
cfDNA.filt$Model = ft
cfDNA.filt$Exp2 = vec.exp
names(cfDNA.filt) = c("b","v","m","Exp","Ms","cfDNA","Model", "Exp2")

cfDNA.filt.cs = cfDNA.filt %>% filter(Exp %in% c("Exp80")) # Data frame for Figure SUPP 6E clone splitting experiment

cfDNA.filt = cfDNA.filt %>% filter(!Exp %in% c("Exp80")) #clone splitting experiment
cfDNA.filt = cfDNA.filt %>% filter(!cfDNA %in% c("cfDNAr", "cfDNAr1", "cfDNAr2","cfDNAf"))#cfDNA on tumour resection

## --- [FIGURE 3A] Success cfDNA recovery ---- 

cf.info = info.1 %>% select("Exp", 'Mice #', "Project", 'Injected cells', "Cells", "MI", 
                            "Cln.Spl", "Barcode", "Tumour removed",
                            '1st Bleeding Blood  Volume','2nd Bleeding  Volume', 'Terminal Endbleed volume', "Protocol")

names(cf.info) = c("Exp", 'Mice#', "Project", 'Injected cells', "Cells", "MI", "Cln.Spl", "Barcode", "Tum.removed",
                   'cf1.Volume','cf2.Volume', 'cf3.volume', "meth")

cf.info1 = cf.info %>% filter(str_detect(cf.info$meth, "cfDNA") & 
                                str_detect(cf.info$Barcode, "Genetic") &
                                MI =="IMFP" &
                                Cln.Spl == "No"&
                                is.na(Tum.removed))

cf.info2 = cf.info1 %>% select('Mice#', "Cells","cf1.Volume", "cf2.Volume", "cf3.volume")
names(cf.info2) = c("Ms", "Model", "cfDNA1", "cfDNA2", "cfDNA3")

cf.info3 = melt(cf.info2, id=c("Ms", "Model"))

cf.info3$value2 = as.numeric(as.character(str_replace_all(cf.info3$value, "[^[:alnum:]]", ""))) #remove volume with '*' at the end

ggplot()+
  geom_tile(data=cf.info3, aes(x=variable, y=as.factor(Ms), fill=value2>0), colour = "black")+
  geom_tile(data=subset(cfDNA.filt,m>0), aes(x=cfDNA, y=Ms),colour = "black", fill="darkgreen", alpha=0.5)+
  labs(title = "CfDNA samples by Model")+
  facet_wrap(~Model, scale="free_y", nrow = 1)+
  scale_fill_manual(na.value = NA, values="gray")+
  theme_void()+xlab("")+ylab("")+
  theme(legend.position = "none")# 4x15


## --- [FIGURE 3B] Example graph cfDNA  -----

Ms.x = df.PCR.Tum.Filt[,grep("_77", names(df.PCR.Tum.Filt))]
colnames(Ms.x) <- str_split(names(Ms.x), "_", simplify = T)[,3]
Ms.x = Ms.x[,!is.na(colSums(Ms.x)), drop=F]
Ms.x = Ms.x[,c("cfDNA1", "cfDNA2", "cfDNA3", "Tum")]
colSums(Ms.x>0)
Stacked_histo(Ms.x)[1]


df.Fig3 = data.frame()

for (i in unique(df.Fig1.clean$Mouse)) {
  print(i)
  Ms.name = paste("_", i, "_", sep="")
  cf = c("Tum", "cfDNA", "Lung")
  
  info.1 = info %>% filter(Clever_cutting==TRUE)
  Ms.info = info.1 %>% filter(`Mice #`==i)
  
  Ms.x = df.PCR.Tum.Filt.I[,grep(Ms.name, names(df.PCR.Tum.Filt.I)), drop=F]
  if(ncol(Ms.x)==0){next}
  
  id.name <- names(Ms.x[1])
  
  colnames(Ms.x) <- str_split(names(Ms.x), "_", simplify = T)[,3]
  if(is.na(sum(Ms.x[,"Tum"]))==T){next}
  
  Ms.x.cf = df.PCR.Tum.Filt.I[,grep(paste(paste(Ms.name,cf, sep=""), collapse = "|"), names(df.PCR.Tum.Filt.I)), drop=F]
  if(ncol(Ms.x.cf)<=1){next}
  colnames(Ms.x.cf) <- str_split(names(Ms.x.cf), "_", simplify = T)[,3]
  Ms.x.cf = Ms.x.cf[,!is.na(colSums(Ms.x.cf)), drop=F]
  Ms.x.cf.temp = Ms.x.cf[rowSums(Ms.x.cf)>0,]
  e = biomass_fun(Ms.x.cf.temp)
  #if(sum(str_detect(names(Ms.x.cf.temp), "cfDNA3"))==0){next}
  
  Ms.x.cf.df = data.frame(nb_bc = colSums(Ms.x.cf.temp>0),
                          pct_bc = colSums(Ms.x.cf.temp>0)*100/colSums(Ms.x.cf.temp[,"Tum", drop=F]>0),
                          Names= colnames(Ms.x.cf.temp),
                          Shanon=diversity(t(Ms.x.cf.temp)),
                          Biom_Tum = e[,"Tum"],
                          Mouse = i,
                          Exp =Ms.info$Exp,
                          Model =Ms.info$Cells)
  
  df.Fig3 = rbind(df.Fig3,Ms.x.cf.df)
}

#Function to compute error bar by selecting column name and group
error.df = function(df, filt.col, grp.b, uniq, fil.value){
  # fil.col = Name of the column do to filtering on
  # fil.val = Name of the value to filter
  # grp.b   = Name of the col to group_by
  # uniq    = Name col to do error bar calculation on
  
  df %>%
    filter(!!filt.col == fil.value) %>%
    group_by(!!grp.b) %>%
    dplyr::mutate(mean_val = mean(!!uniq),
                  sd_val= sd(!!uniq),
                  n_val = n(),
                  SE_val = sd(!!uniq)/sqrt(n())) %>%
    ungroup()
  
}

df.Fig3.errorbar.Biom = error.df(df.Fig3,
                                 filt.col=sym("Names"),
                                 grp.b=sym("Model"),
                                 uniq=sym("Biom_Tum"),
                                 fil.value = "cfDNA3")
## --- [FIGURE 3E] Tumour biomass captured from cfDNA  -----

df.Fig3.errorbar.Biom %>% 
  ggplot(aes(x=Model, y=Biom_Tum))+
  geom_bar(stat = "summary", fun="mean")+
  geom_point()+
  geom_errorbar(aes(ymin = ifelse(mean_val - sd_val<0, 0,mean_val - sd_val),ymax = mean_val + sd_val),width=0.2)+
  labs(title = "Tum Biomass captured", x="", y="Tumour biomass captured in cfDNA3")+
  theme_classic()#4.5x4.5

#number of mice plotted in graph:
n.Fig3.d = df.Fig3.errorbar.Biom %>% select(Model, Exp, Mouse) %>% group_by(Model) %>% 
  dplyr::mutate(count.Ms = n_distinct(Mouse), count.Exp = n_distinct(Exp)) %>%
  ungroup() %>% select(Model, count.Ms, count.Exp) %>% unique()

n.Ms.list[[length(n.Ms.list)+1]] = n.Fig3.d
names(n.Ms.list)[length(n.Ms.list)] = "Figure_3d.e"

## --- [FIGURE 3D] Percentage of barcode from tumour captured in cfDNA  -----

df.Fig3.errorbar.pct.bc = error.df(df.Fig3, sym("Names"), sym("Model"), sym("pct_bc"), "cfDNA3")

df.Fig3.errorbar.pct.bc %>% 
  ggplot(aes(x=Model, y=pct_bc))+
  geom_bar(stat = "summary", fun="mean")+
  geom_point()+labs(title = "Percentage barcode")+
  geom_errorbar(aes(ymin = ifelse(mean_val - sd_val<0, 0,mean_val - sd_val),ymax = mean_val + sd_val),width=0.2)+
  theme_classic()#4.5x4.5

## --- [FIGURE SUPP 5B] Number of barcode from tumour captured in cfDNA  -----

df.Fig3.errorbar.nb_bc = error.df(df.Fig3, sym("Names"), sym("Model"), sym("nb_bc"), "cfDNA3")
df.Fig3.errorbar.nb_bc %>% filter(Names == "cfDNA3") %>%
  ggplot(aes(x=Model, y=nb_bc))+
  geom_bar(stat = "summary", fun="mean")+
  geom_point()+labs(title = "Number of bc")+
  geom_errorbar(aes(ymin = ifelse(mean_val - sd_val<0, 0,mean_val - sd_val),ymax = mean_val + sd_val),width=0.2)+
  theme_classic()#4.5x4.5

df.Fig3.Supp = df.Fig3 %>% filter(Names=="cfDNA1"|Names=="cfDNA2"|Names=="cfDNA3") %>%
  filter(Model=="PDX-T412"|Model=="PDX-4295")

n.Supp.Fig5.a = df.Fig3.Supp %>% select(Model, Exp, Mouse) %>% group_by(Model) %>% 
  dplyr::mutate(count.Ms = n_distinct(Mouse), count.Exp = n_distinct(Exp)) %>%
  ungroup() %>% select(Model, count.Ms, count.Exp) %>% unique()

n.Ms.list[[length(n.Ms.list)+1]] = n.Supp.Fig5.a
names(n.Ms.list)[length(n.Ms.list)] = "Supp.Figure_5b"


## --- [FIGURE SUPP 5A] Tumour biomass capture at different cfDNA timepoint -----

df.Fig3.Supp.1B = df.Fig3.Supp %>% filter(Names=="cfDNA1"|Names=="cfDNA2"|Names=="cfDNA3")
df.Fig3.Supp.1B.er = error.df(df.Fig3.Supp.1B, sym("Names"), sym("Model"), sym("Biom_Tum"), c("cfDN3", "cfDNA2","CfDNA1"))


df.Fig3.Supp.1B.er = df.Fig3.Supp.1B %>% group_by(Model, Names) %>% 
  dplyr::mutate(mean_biom = mean(Biom_Tum),sd_sha= sd(Biom_Tum),n_sha = n(),SE_sha = sd(Biom_Tum)/sqrt(n())) %>% ungroup()

# ANOVA on cfDNA sample per models
compare_means(Biom_Tum ~ Names,
                         group.by = "Model",
                         data = df.Fig3.Supp.1B.er,
                         p.adjust.method = "holm",
                         method = "anova") %>%
  mutate(p.adj = format.pval(p.adj, digits = 2))
# Follow up Tukey test
TukeyHSD(aov(Biom_Tum~Names, data = df.Fig3.Supp.1B.er%>% filter(Model == "PDX-T412") %>% select(Biom_Tum, Names)))

ggplot(data=df.Fig3.Supp.1B.er, aes(x=Names, y=Biom_Tum))+
  geom_bar(stat = "summary", fun="mean")+
  geom_point()+
  geom_errorbar(aes(ymin = ifelse(mean_biom - sd_sha<0, 0,mean_biom - sd_sha),
                   ymax = mean_biom + sd_sha),width=0.2, position = "dodge")+
  labs(title = "Tum Biomass captured")+
  facet_wrap(~Model, nrow = 1)+
  theme_classic()#4.5x4.5

## --- [FIGURE 3C & SUPP 5C] Scatter plots cfDNA ---- 

df.Fig3D=data.frame()
for (i in unique(df.Fig3$Mouse)) {
  print(i)
  info.1 = info %>% filter(Clever_cutting==TRUE)
  Ms.info = info.1 %>% filter(`Mice #`==i)
  
  TEST.name = paste("_", i, "_", sep="")
  Ms.x = df.PCR.Tum.Filt.I[,grep(TEST.name, names(df.PCR.Tum.Filt.I))]
  colnames(Ms.x) <- str_split(names(Ms.x), "_", simplify = T)[,3]
  if(ncol(Ms.x)==0){No_Datum = c(No_Datum,i);next}
  id.name <- names(Ms.x[1])
  
  Ms.x.Core = paste("Tum", c(as.numeric(as.character(unlist(str_split(Ms.info$Core,","))))), sep="")
  Ms.x.Tum = df.PCR[,grep(TEST.name, names(df.PCR))]
  colnames(Ms.x.Tum) <- str_split(names(Ms.x.Tum), "_", simplify = T)[,3]
  
  if(sum(Ms.x.Core %in% "TumNA")>0){next
  } else if(length(Ms.x.Core)==1& sum(Ms.x.Core %in% names(Ms.x.Tum)==TRUE)>0){Ms.x.Tum$Core = Ms.x.Tum[,Ms.x.Core]
  }else {Ms.x.Tum$Core = rowSums(Ms.x.Tum[,Ms.x.Core])
  }
  
  Ms.x = Ms.x[,!is.na(colSums(Ms.x))]
  Ms.x = cbind(Ms.x, Core=Ms.x.Tum$Core)
  Ms.x = Ms.x[rowSums(Ms.x)>0,]
  
  aa=Ms.x
  Needles = c("Tum", "Core", "N1a", "N1b", "N2a", "N2b","N3a", "N3b", "N4a", "N4b")
  a=Ms.x[,!names(Ms.x) %in% Needles, drop=F]
  cond.lung = sum(str_detect(names(a), "Lung"))
  a=a[,!str_detect(names(a), "cfDNA"), drop=F]
  aa$Mets=(aa$Tum>0)+(as.vector(rowSums(a)>0))==2
  if(cond.lung==0){aa$is.Lung= 0} else{
    aa$is.Lung = (aa$Tum>0)+(a[,str_detect(names(a), "Lung")]>0)==2}
  
  if(sum(c("Tum", "Core", "cfDNA3", "Mets", "is.Lung") %in% names(aa))<5){next
  }else{
    Ms.x.temp= aa[,c("Tum", "Core", "cfDNA3", "Mets", "is.Lung")]
    Ms.x.temp$mouse = i
    Ms.x.temp$Model = Ms.info$Cells
    Ms.x.temp$Exp = Ms.info$Exp
  }
  
  df.Fig3D = rbind(df.Fig3D, Ms.x.temp)
}

df.Fig3D = df.Fig3D %>% group_by(Model) %>% 
  dplyr::mutate(count = n_distinct(mouse)) %>% ungroup() %>% 
  dplyr::mutate(Model2 = paste(.$Model, " (n=", .$count, ")", sep=""))

df.Fig3D.nb= subset(df.Fig3D,cfDNA3>0&Tum>0)
n.Supp.Fig5.a = df.Fig3D.nb %>% group_by(Model) %>% 
  dplyr::summarise(n = n_distinct(mouse), exp.n = n_distinct(Exp))

n.Ms.list[[length(n.Ms.list)+1]] = n.Supp.Fig5.a
names(n.Ms.list)[length(n.Ms.list)] = "Figure_3c"


# [FIGURE 3C]
ggplot(df.Fig3D)+
  geom_point(data= subset(df.Fig3D,cfDNA3>0&Tum>0),aes(x=cfDNA3, y=Tum, color=is.Lung>0))+
  geom_quasirandom(data=subset(df.Fig3D,cfDNA3==0 & Tum>0), aes(x=10, y=Tum, color=is.Lung>0), width=0.9, groupOnX=TRUE,alpha=0.7)+
  #geom_quasirandom(data=subset(df.Fig3D,cfDNA3>0 & Tum==0), aes(x=cfDNA3, y=1e-5), groupOnX=FALSE,alpha=0.5)+
  scale_color_manual(values = c("dimgray", "darkorange"), labels=c("No", "Yes"))+scale_x_log10()+scale_y_log10()+xlab("Frequency in cfDNA3")+ylab("Frequency in Tumour")+
  labs(title = "Scatter plot cfDNA3 vs Tumour in Lung", color= "in Lung")+
  facet_wrap(~Model2, nrow = 1)+theme_classic()#3.5x15

# [FIGURE SUPP 5c]
ggplot(df.Fig3D)+
  geom_point(data= subset(df.Fig3D,cfDNA3>0&Tum>0),aes(x=cfDNA3, y=Tum, color=Core>0))+
  geom_quasirandom(data=subset(df.Fig3D,cfDNA3==0 & Tum>0), aes(x=10, y=Tum, color=Core>0), width=0.9, groupOnX=TRUE,alpha=0.7)+
  #geom_quasirandom(data=subset(df.Fig3D,cfDNA3>0 & Tum==0), aes(x=cfDNA3, y=1e-5), groupOnX=FALSE,alpha=0.5)+
  scale_color_hue(direction = -1, labels=c("Peri", "Core"))+scale_x_log10()+scale_y_log10()+xlab("Frequency in cfDNA3")+ylab("Frequency in Tumour")+
  labs(title = "Scatter plot cfDNA3 vs Tumour", color= "Position")+
  facet_wrap(~Model2, nrow = 1)+theme_classic()#3.5x15

#----------------   FIGURE 4 ==========================    ----
## --- Processing datum for Fig4 ----
set.seed(112233)
df.Fig4 = data.frame()
for (i in unique(df.Fig1.clean$Mouse)) {
  print(i)
  TEST.name = paste("_", i, "_", sep="")
  
  info.1 = info %>% filter(Clever_cutting==TRUE)
  Ms.info = info.1 %>% filter(`Mice #`==i)
  
  Ms.x = df.PCR.Tum.Filt.I[,grep(TEST.name, names(df.PCR.Tum.Filt.I)), drop=FALSE]
  if(ncol(Ms.x)==0){No_Datum = c(No_Datum,i);next}
  
  id.name <- names(Ms.x[1])
  
  colnames(Ms.x) <- str_split(names(Ms.x), "_", simplify = T)[,3]
  if(is.na(sum(Ms.x[,"Tum"]))==T){NA_tum = c(NA_tum,i);next}
  
  Ms.x = Ms.x[,!is.na(colSums(Ms.x)), drop=F]
  
  if(sum(names(Ms.x) %in% "cfDNA3")==0 & sum(names(Ms.x) %in% c("N1a", "N1b", "N2a", "N2b","N3a", "N3b", "N4a", "N4b"))==0){next}
  if(sum(names(Ms.x) %in% "cfDNA3")==0 | sum(names(Ms.x) %in% c("N1a", "N1b", "N2a", "N2b","N3a", "N3b", "N4a", "N4b"))==0){Sampl="one";Ms.x$cfDNA3 = 0} else {Sampl="both"}
  
  Ms.x.1 = subset( Ms.x, select = -Tum )
  Ms.x.1[Ms.x.1<10000]=0 # FILTER BARCODE UNDER 1% in needles
  Ms.x.2 = cbind(Ms.x[,"Tum", drop=F],Ms.x.1)
  
  if(sum(names(Ms.x) %in% c("N1a","N2a","N3a","N4a"))>=2){
    combo = sample(names(Ms.x)[names(Ms.x) %in% c("N1a","N2a","N3a","N4a")],2)
    nee.combo = as.data.frame(Ms.x[,combo[1]]+Ms.x[,combo[2]])
    names(nee.combo) = "Na.combo"
    Ms.x = cbind(Ms.x,nee.combo)
  }
  
  if(sum(names(Ms.x) %in% c("N1b","N2b","N3b","N4b"))>=2){
    combo = sample(names(Ms.x)[names(Ms.x) %in% c("N1b","N2b","N3b","N4b")],2)
    nee.combo.b = as.data.frame(Ms.x[,combo[1]]+Ms.x[,combo[2]])
    names(nee.combo.b) = "Nb.combo"
    Ms.x = cbind(Ms.x,nee.combo.b)
  }
  
  if(sum(names(Ms.x.2) %in% c("N1a","N2a","N3a","N4a"))>=2){
    combo = sample(names(Ms.x.2)[names(Ms.x.2) %in% c("N1a","N2a","N3a","N4a")],2)
    nee.combo = as.data.frame(Ms.x.2[,combo[1]]+Ms.x.2[,combo[2]])
    names(nee.combo) = "Na.combo"
    Ms.x.2 = cbind(Ms.x.2,nee.combo)
  }
  
  if(sum(names(Ms.x.2) %in% c("N1b","N2b","N3b","N4b"))>=2){
    combo = sample(names(Ms.x.2)[names(Ms.x.2) %in% c("N1b","N2b","N3b","N4b")],2)
    nee.combo.b = as.data.frame(Ms.x.2[,combo[1]]+Ms.x.2[,combo[2]])
    names(nee.combo.b) = "Nb.combo"
    Ms.x.2 = cbind(Ms.x.2,nee.combo.b)
  }
  
  if(sum(names(Ms.x) %in% "cfDNA3")==1){
    
    test= Ms.x[,(colnames(Ms.x) %in% c("N1a", "N1b", "N2a", "N2b","N3a", "N3b", "N4a", "N4b"))]+Ms.x[,"cfDNA3"]
    names(test) = paste(names(test),"cf", sep = ".")
    Ms.x = cbind(Ms.x, test)
    if(sum(names(Ms.x) %in% c("CTC","Blood"))==1){
      test.2= test[,(colnames(test) %in% c("N1a.cf", "N1b.cf", "N2a.cf", "N2b.cf","N3a.cf", "N3b.cf", "N4a.cf", "N4b.cf"))]+Ms.x[,(colnames(Ms.x) %in% c("CTC","Blood"))]
      names(test.2) = paste(names(test.2),"CTC", sep = ".")
      Ms.x = cbind(Ms.x, test.2)
    }
  }
  
  if(sum(names(Ms.x.2) %in% "cfDNA3")==1){
    
    test= Ms.x.2[,(colnames(Ms.x.2) %in% c("N1a", "N1b", "N2a", "N2b","N3a", "N3b", "N4a", "N4b"))]+Ms.x.2[,"cfDNA3"]
    names(test) = paste(names(test),"cf", sep = ".")
    Ms.x.2 = cbind(Ms.x.2, test)
    if(sum(names(Ms.x) %in% c("CTC","Blood"))==1){
      test.2= test[,(colnames(test) %in% c("N1a.cf", "N1b.cf", "N2a.cf", "N2b.cf","N3a.cf", "N3b.cf", "N4a.cf", "N4b.cf"))]+Ms.x[,(colnames(Ms.x) %in% c("CTC","Blood"))]
      names(test.2) = paste(names(test.2),"CTC", sep = ".")
      Ms.x.2 = cbind(Ms.x.2, test.2)
    }
  }
  
  e = biomass_fun(Ms.x)
  e1 = biomass_fun(Ms.x.2)
  
  Ms.x.df = data.frame(nb_bc = colSums(Ms.x>0),
                       pct_bc = colSums(Ms.x>0)*100/colSums(Ms.x[,"Tum", drop=F]>0),
                       Names= colnames(Ms.x),
                       Shanon=diversity(t(Ms.x)),
                       Biom_Tum = e[,"Tum"],
                       Corr_Tum = cor(log10(Ms.x+1))[,1],
                       Corr_Tum2 = cor(log2(Ms.x+1))[,1],
                       Corr_Tum3 = cor(Ms.x)[,1],
                       Corr_Tum4 = cor(log10(Ms.x+1), use = "pairwise.complete.obs")[,1],
                       Biom_Tum1 = e1[,"Tum"],
                       #Biom_Tum10 = e2[,"Tum"],
                       Mouse = i,
                       Exp =Ms.info$Exp,
                       Model =Ms.info$Cells,
                       Sampling = Sampl)
  df.Fig4 = rbind(df.Fig4,Ms.x.df)
}

df.Fig4$Raw.Names = df.Fig4$Names
df.Fig4$Names = str_replace(df.Fig4$Names, c("N1a|N2a|N3a|N4a"), "Needles1")
df.Fig4$Names = str_replace(df.Fig4$Names, c("N1b|N2b|N3b|N4b"), "Needles2")

df.Fig4$Names = str_replace(df.Fig4$Names, c("N1a.cf|N2a.cf|N3a.cf|N4a.cf"), "Needles1.cf")
df.Fig4$Names = str_replace(df.Fig4$Names, c("N1b.cf|N2b.cf|N3b.cf|N4b.cf"), "Needles2.cf")

df.Fig4$Names = str_replace(df.Fig4$Names, c("N1a.cf.CTC|N2a.cf.CTC|N3a.cf.CTC|N4a.cf.CTC"), "Needles1.cf.CTC")
df.Fig4$Names = str_replace(df.Fig4$Names, c("N1b.cf.CTC|N2b.cf.CTC|N3b.cf.CTC|N4b.cf.CTC"), "Needles2.cf.CTC")

## --- [FIGURE SUPP 6A]  Correlation with Tumour ----

error.df.no.filt = function(df, uniq, ...){
  # uniq name col to do error bar calculation on
  # ... multiple col name to group
  groupCol <- quos(...)
  df %>%
    group_by(!!!groupCol) %>%
    dplyr::mutate(mean_val = mean(!!uniq),
                  sd_val= sd(!!uniq),
                  n_val = n(),
                  SE_val = sd(!!uniq)/sqrt(n())) %>%
    ungroup()
  
}

df.Fig4.B = df.Fig4 %>% filter(Sampling=="both") %>% filter(Names=="Needles1"|Names=="cfDNA3"|Names=="Needles2")
df.Fig4.B.er = error.df.no.filt(df.Fig4.B, uniq=sym("Corr_Tum"), Model, Names)

n.Fig4.b = df.Fig4.B.er %>% select(Model, Exp, Mouse) %>% group_by(Model) %>% 
  dplyr::mutate(count.Ms = n_distinct(Mouse), count.Exp = n_distinct(Exp)) %>%
  ungroup() %>% select(Model, count.Ms, count.Exp) %>% unique()

n.Ms.list[[length(n.Ms.list)+1]] = n.Fig4.b
names(n.Ms.list)[length(n.Ms.list)] = "Figure_4b"

# ANOVA on cfDNA sample per models
compare_means(Corr_Tum ~ Names,
              group.by = "Model",
              data = df.Fig4.B.er,
              p.adjust.method = "holm",
              method = "anova") %>%
  mutate(p.adj = format.pval(p.adj, digits = 2))
# Follow up Tukey test
TukeyHSD(aov(Corr_Tum~Names, data = df.Fig4.B.er%>% filter(Model == "MDA-231") %>% select(Corr_Tum, Names)))
TukeyHSD(aov(Corr_Tum~Names, data = df.Fig4.B.er%>% filter(Model == "PDX-1432C") %>% select(Corr_Tum, Names)))
TukeyHSD(aov(Corr_Tum~Names, data = df.Fig4.B.er%>% filter(Model == "PDX-434") %>% select(Corr_Tum, Names)))
TukeyHSD(aov(Corr_Tum~Names, data = df.Fig4.B.er%>% filter(Model == "PDX-T412") %>% select(Corr_Tum, Names)))

ggplot(data=df.Fig4.B.er, aes(x=Names, y=Corr_Tum))+
  geom_bar(aes(fill=Names), stat = "summary", fun=mean)+
  geom_point(colour="grey50")+
  geom_errorbar(aes(ymin = ifelse(mean_val - sd_val<0, 0,mean_val - sd_val),ymax = mean_val + sd_val),width=0.2)+
  facet_wrap(~Model, nrow=1)+
  theme_classic()+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_colorblind()+
  xlab("")+ylab("Pearson correlation")+labs(title = "Correlation with Primary tumour (Log10)")#5x13

## --- [FIGURE 4B]  Percentage barcodes captured ----
df.Fig4.A = df.Fig4 %>% filter(Sampling=="both"|Sampling=="one") %>% filter(Names=="Tum"|Names=="Needles1"|Names=="cfDNA3"|Names=="Needles2")

df.Fig4.A.er = error.df.no.filt(df.Fig4.A, uniq=sym("pct_bc"), Model, Names)

n.Fig4.a = df.Fig4.A.er %>% select(Model, Exp, Mouse) %>% group_by(Model) %>% 
  dplyr::mutate(count.Ms = n_distinct(Mouse), count.Exp = n_distinct(Exp)) %>%
  ungroup() %>% select(Model, count.Ms, count.Exp) %>% unique()

n.Ms.list[[length(n.Ms.list)+1]] = n.Fig4.a
names(n.Ms.list)[length(n.Ms.list)] = "Figure_4a"

# ANOVA on cfDNA sample per models
compare_means(pct_bc ~ Names,
              group.by = "Model",
              data = df.Fig4.A.er %>% filter(Names!="Tum"),
              p.adjust.method = "holm",
              method = "anova") %>%
  mutate(p.adj = format.pval(p.adj, digits = 2))
# Follow up Tukey test
TukeyHSD(aov(pct_bc~Names, data = df.Fig4.A.er%>% filter(Model == "MDA-231" &Names!="Tum") %>% select(pct_bc, Names)))
TukeyHSD(aov(pct_bc~Names, data = df.Fig4.A.er%>% filter(Model == "MDA-468"&Names!="Tum") %>% select(pct_bc, Names)))
TukeyHSD(aov(pct_bc~Names, data = df.Fig4.A.er%>% filter(Model == "PDX-1432C"&Names!="Tum") %>% select(pct_bc, Names)))
TukeyHSD(aov(pct_bc~Names, data = df.Fig4.A.er%>% filter(Model == "PDX-4295"&Names!="Tum") %>% select(pct_bc, Names)))
TukeyHSD(aov(pct_bc~Names, data = df.Fig4.A.er%>% filter(Model == "PDX-434"&Names!="Tum") %>% select(pct_bc, Names)))
TukeyHSD(aov(pct_bc~Names, data = df.Fig4.A.er%>% filter(Model == "PDX-T412"&Names!="Tum") %>% select(pct_bc, Names)))


ggplot(data=df.Fig4.A.er)+
  geom_bar(aes(x=Names, y=pct_bc, fill=Names), stat = "summary", fun=mean)+
  geom_quasirandom(aes(x=Names, y=pct_bc), colour="grey50")+
  facet_wrap(~Model, nrow=1)+
  theme_classic()+
  scale_fill_colorblind()+scale_x_discrete(limits = c("Tum", "cfDNA3", "Needles1", "Needles2"))+
  geom_errorbar(aes(x=Names, y=pct_bc, ymin = ifelse(mean_val - sd_val<0, 0,mean_val - sd_val),ymax = mean_val + sd_val),width=0.2)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("")+ylab("% of barcodes captured")+labs(title = "Percentages of barcodes captured") #5x13

## --- [FIGURE 4F]  Biomass captured no treshold  ----

df.Fig4.F = df.Fig4 %>% filter(Sampling=="both") %>% filter(Names=="Needles1"|Names=="cfDNA3"|Names=="Needles2")
df.Fig4.F.er = error.df.no.filt(df.Fig4.F, uniq=sym("Biom_Tum"), Model, Names)

# ANOVA on cfDNA sample per models
compare_means(Biom_Tum ~ Names,
              group.by = "Model",
              data = df.Fig4.F.er,
              p.adjust.method = "holm",
              method = "anova") %>%
  mutate(p.adj = format.pval(p.adj, digits = 2))
# Follow up Tukey test
TukeyHSD(aov(Biom_Tum~Names, data = df.Fig4.F.er%>% filter(Model == "PDX-1432C") %>% select(Biom_Tum, Names)))
TukeyHSD(aov(Biom_Tum~Names, data = df.Fig4.F.er%>% filter(Model == "PDX-434") %>% select(Biom_Tum, Names)))
TukeyHSD(aov(Biom_Tum~Names, data = df.Fig4.F.er%>% filter(Model == "PDX-T412") %>% select(Biom_Tum, Names)))

ggplot(data=df.Fig4.F.er)+
  geom_bar(aes(x=Names, y=Biom_Tum, fill=Names), stat = "summary", fun=mean)+
  geom_point(aes(x=Names, y=Biom_Tum), colour="grey50")+
  facet_wrap(~Model, nrow=1)+
  theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_colorblind()+
  geom_errorbar(aes(x=Names, y=Biom_Tum, ymin = ifelse(mean_val - sd_val<0, 0,mean_val - sd_val),ymax = mean_val + sd_val),width=0.2)+
  xlab("")+ylab("Tumour biomass")+labs(title = "Primary tumour Biomass captured")#5x13

## --- [FIGURE 4H]  Biomass captured combination N+cf / N+cf+CTC / N+N ----

df.Fig4.H = df.Fig4 %>% filter(Sampling=="both") %>% filter(Names=="Needles1"|Names=="cfDNA3"|Names=="Needles2"|Names=="Needles1.cf"|Names=="Needles2.cf")

df.Fig4.H.er = error.df.no.filt(df.Fig4.H, uniq=sym("Biom_Tum"), Model, Names)

n.Fig4.h = df.Fig4.H.er %>% select(Model, Exp, Mouse) %>% group_by(Model) %>% 
  dplyr::mutate(count.Ms = n_distinct(Mouse), count.Exp = n_distinct(Exp)) %>%
  ungroup() %>% select(Model, count.Ms, count.Exp) %>% unique()

n.Ms.list[[length(n.Ms.list)+1]] = n.Fig4.h
names(n.Ms.list)[length(n.Ms.list)] = "Figure_4h"

ID = df.Fig4.H.er %>% select(Model, Exp, Mouse) %>% distinct()
names(ID) = c("Model", "Exp", "Mouse")
ID$Figure = "Figure_4h"
Mouse_ID= rbind(Mouse_ID, ID)

# ANOVA on cfDNA sample per models

compare_means(Biom_Tum ~ Names,
              group.by = "Model",
              data = df.Fig4.H.er %>% filter(Sampling=="both") %>%filter(Names=="Needles1"|Names=="cfDNA3"|Names=="Needles1.cf"),
              p.adjust.method = "holm",
              method = "anova") %>%
  mutate(p.adj = format.pval(p.adj, digits = 2))
# Follow up Tukey test
TukeyHSD(aov(Biom_Tum~Names, data = df.Fig4.H.er%>% filter(Model == "MDA-468") %>% filter(Sampling=="both") %>% filter(Names=="Needles1"|Names=="cfDNA3"|Names=="Needles1.cf")%>% select(Biom_Tum, Names)))
TukeyHSD(aov(Biom_Tum~Names, data = df.Fig4.H.er%>% filter(Model == "PDX-1432C") %>% filter(Sampling=="both") %>%filter(Names=="Needles1"|Names=="cfDNA3"|Names=="Needles1.cf")%>% select(Biom_Tum, Names)))
TukeyHSD(aov(Biom_Tum~Names, data = df.Fig4.H.er%>% filter(Model == "PDX-434") %>% filter(Sampling=="both") %>%filter(Names=="Needles1"|Names=="cfDNA3"|Names=="Needles1.cf")%>% select(Biom_Tum, Names)))
TukeyHSD(aov(Biom_Tum~Names, data = df.Fig4.H.er%>% filter(Model == "PDX-T412") %>% filter(Sampling=="both") %>%filter(Names=="Needles1"|Names=="cfDNA3"|Names=="Needles1.cf")%>% select(Biom_Tum, Names)))

df.Fig4.H.er %>% filter(Sampling=="both") %>%filter(Names=="Needles1"|Names=="cfDNA3"|Names=="Needles1.cf") %>%
  ggplot(aes(x=Names, y=Biom_Tum))+
  geom_bar(aes(fill=Names), stat = "summary", fun=mean)+
  geom_point(colour="grey50")+
  geom_errorbar(aes(ymin = ifelse(mean_val - sd_val<0, 0,mean_val - sd_val),ymax = mean_val + sd_val),width=0.2)+
  facet_wrap(~Model, nrow=1)+
  theme_classic()+
  scale_fill_colorblind()+
  scale_y_continuous(breaks=c(0,250000,500000,750000,1000000))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("")+ylab("Tumour biomass")+labs(title = "Primary tumour Biomass captured with combination")#5x13

## --- [FIGURE SUPP 7B] Biomass captured combination N+cf+CTC ----

ms.cf.CTC = df.Fig4 %>% filter(Names == "Needles1.cf.CTC"|Names=="Needles2.cf.CTC") %>% select(Mouse) %>% unique()

df.Fig4.supp4.2.B = df.Fig4 %>% filter(Mouse %in% ms.cf.CTC$Mouse) %>%
  filter(Names=="Needles1"|Names=="cfDNA3"|Names=="Needles2"|Names=="Needles1.cf"|Names=="Needles2.cf"|Names=="Needles1.cf.CTC"|Names=="Needles2.cf.CTC")

df.Fig4.supp4.2.B.er = error.df.no.filt(df.Fig4.supp4.2.B, uniq=sym("Biom_Tum1"), Model, Names)


# comparison of tumour biomass capture by combination of cfDNA+Needles(Needles.cf) vs cfDNA+Needles+CTCs (Needles.cf.CTC)
  # For needles1 (deep) and needles2 (shallow)
compare_means(Biom_Tum1 ~ Names,
              group.by = "Model",
              data = df.Fig4.supp4.2.B %>% filter(Model == "MDA-231" & Names %in% c("Needles1.cf", "Needles1.cf.CTC", "Needles2.cf", "Needles2.cf.CTC")),
              p.adjust.method = "holm",
              method = "t.test",
              paired = TRUE) %>%
  mutate(p.adj = format.pval(p.adj, digits = 2))

ggplot(data=df.Fig4.supp4.2.B.er %>% filter(Model == "MDA-231"), aes(x=Names, y=Biom_Tum1))+
  geom_bar(aes(fill=Names), stat = "summary", fun=mean)+
  geom_point(colour="grey50")+
  geom_errorbar(aes(ymin = ifelse(mean_val - sd_val<0, 0,mean_val - sd_val),ymax = mean_val + sd_val),width=0.2)+
  facet_wrap(~Model, nrow=1, scales = "free_x")+
  theme_classic()+
  scale_fill_colorblind()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("")+ylab("Tumour biomass")+labs(title = "Primary tumour Biomass captured with combination adding CTCs (1%)")#5x8


## --- [FIGURE SUPP 6C] Biomass captured combination N+N ----

df.Fig4.supp4.1D = df.Fig4 %>% filter(Sampling=="both") %>% filter(Names=="Needles1"|Names=="Needles2"|Names=="Na.combo"|Names=="Nb.combo")
df.Fig4.supp4.1D.er = error.df.no.filt(df.Fig4.supp4.1D, uniq=sym("Biom_Tum"), Model, Names)

n.Supp.Fig6.b = df.Fig4.supp4.1D.er %>% select(Model, Exp, Mouse) %>% group_by(Model) %>% 
  dplyr::mutate(count.Ms = n_distinct(Mouse), count.Exp = n_distinct(Exp)) %>%
  ungroup() %>% select(Model, count.Ms, count.Exp) %>% unique()

n.Ms.list[[length(n.Ms.list)+1]] = n.Supp.Fig6.b
names(n.Ms.list)[length(n.Ms.list)] = "Supp_Figure_6b"

anno_df4 = compare_means(Biom_Tum ~ Names,
              group.by = "Model",
              data = df.Fig4.supp4.1D.er %>% filter(Model!="MDA-468" & Sampling=="both") %>% filter(Names=="Needles1"|Names=="Na.combo"),
              p.adjust.method = "holm",
              method = "t.test") %>%
  mutate(y_pos = 1050000, p.adj = format.pval(p.adj, digits = 2))

df.Fig4.supp4.1D.er %>% filter(Sampling=="both") %>% filter(Names=="Needles1"|Names=="Na.combo") %>%
  ggplot(aes(x=Names, y=Biom_Tum))+
  geom_bar(aes(fill=Names), stat = "summary", fun=mean)+
  geom_point(colour="grey50")+
  geom_errorbar(aes(ymin = ifelse(mean_val - sd_val<0, 0,mean_val - sd_val),ymax = mean_val + sd_val),width=0.2)+
  geom_signif(data=anno_df4,aes(xmin = group1, xmax = group2, annotations = paste(p.signif, " (", p.format,")", sep=""), y_position = y_pos),textsize=3,manual= TRUE)+
  facet_wrap(~Model, nrow=1, scales = "free_x")+
  theme_classic()+
  scale_fill_colorblind()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(limits=c("Needles1","Na.combo"))+
  xlab("")+ylab("Tumour biomass")+labs(title = "Primary tumour Biomass captured with Needles combo")#5x13

## --- [FIGURE SUPP 5D] Lung Biomass, correlation ----
Ms.Lung.cf.N = df.Fig3 %>% filter(Names=="Lung") %>% .$Mouse

cf = c("Tum", "cfDNA", "Lung")

No_cfDNA = vector()
No_corr_cfDNA = vector()
df.Fig4.Lung = data.frame()

for (i in Ms.Lung.cf.N) {
  print(i)
  TEST.name = paste("_", i, "_", sep="")
  
  info.1 = info %>% filter(Clever_cutting==TRUE)
  Ms.info = info.1 %>% filter(`Mice #`==i)
  
  Ms.x = df.PCR.Tum.Filt.I[,grep(TEST.name, names(df.PCR.Tum.Filt.I))]
  if(ncol(Ms.x)==0){No_Datum = c(No_Datum,i);next}
  
  id.name <- names(Ms.x[1])
  
  colnames(Ms.x) <- str_split(names(Ms.x), "_", simplify = T)[,3]
  if(is.na(sum(Ms.x[,"Tum"]))==T){NA_tum = c(NA_tum,i);next}
  Ms.x = Ms.x[,!is.na(colSums(Ms.x)), drop=F]
  Ms.x = Ms.x[rowSums(Ms.x)>0,]
  
  # --- 1% FILTER On Needles and cfDNA ---
  test = Ms.x[,!names(Ms.x) %in% "Lung"]
  test[test<10000]=0
  Ms.x.2= cbind(Lung=Ms.x$Lung, test)
  
  Ms.x.2 = Ms.x.2[match(names(Ms.x), names(Ms.x.2))]
  
  e = biomass_fun(Ms.x)
  e.1 = biomass_fun(Ms.x.2)
  
  Ms.x.Needles.df = data.frame(nb_bc = colSums(Ms.x>0),
                               pct_bc = colSums(Ms.x>0)*100/colSums(Ms.x[,"Tum", drop=F]>0),
                               Names= colnames(Ms.x),
                               Shanon=diversity(t(Ms.x)),
                               Biom_Lung = e[,"Lung"],
                               Biom_Lung1 = e.1[,"Lung"],
                               Corr_Lung = cor(Ms.x)[,"Lung"],
                               #Biom_Tum10 = e2[,"Tum"],
                               Mouse = i,
                               Exp =Ms.info$Exp,
                               Model =Ms.info$Cells)
  
  df.Fig4.Lung = rbind(df.Fig4.Lung,Ms.x.Needles.df)
  
}

df.Fig4.Lung$Names = str_replace(df.Fig4.Lung$Names, c("N1a|N2a|N3a|N4a"), "Needles1")
df.Fig4.Lung$Names = str_replace(df.Fig4.Lung$Names, c("N1b|N2b|N3b|N4b"), "Needles2")

df.Fig4.supp2.c.d = df.Fig4.Lung %>% filter(Names=="Needles1"|Names=="Needles2"|Names=="cfDNA3")
df.Fig4.supp2.c.d.er = error.df.no.filt(df.Fig4.supp2.c.d, uniq=sym("Biom_Lung"), Model, Names)

# ANOVA on cfDNA sample per models
compare_means(Biom_Lung ~ Names,
              group.by = "Model",
              data = df.Fig4.supp2.c.d.er,
              p.adjust.method = "holm",
              method = "anova")
# Follow up Tukey test
TukeyHSD(aov(Biom_Lung~Names, data = df.Fig4.supp2.c.d.er%>% filter(Model == "PDX-T412")%>% select(Biom_Lung, Names)))

ggplot(data=df.Fig4.supp2.c.d.er, aes(x=Names, y=Biom_Lung))+
  geom_bar(aes(fill=Names), position="dodge", stat = "summary", fun="mean")+
  geom_quasirandom(groupOnX = T)+labs(title = "Lung Biomass captured ")+
  geom_errorbar(aes(ymin = ifelse(mean_val - sd_val<0, 0,mean_val - sd_val),ymax = mean_val + sd_val),width=0.2)+
  facet_wrap(~Model, nrow = 1)+
  theme_classic()#6x8

## --- [FIGURE 4 C, D, E] Example histogram ----

# Mice examples in manuscript : c("30", "47", "78")

i=78
TEST.name = paste("_", i, "_", sep="")
info.1 = info %>% filter(Clever_cutting==TRUE)
Ms.info = info.1 %>% filter(`Mice #`==i)
Ms.x = df.PCR.Tum.Filt.I[,grep(TEST.name, names(df.PCR.Tum.Filt.I))]
id.name <- names(Ms.x[1])
colnames(Ms.x) <- str_split(names(Ms.x), "_", simplify = T)[,3]
Ms.x = Ms.x[,!is.na(colSums(Ms.x)), drop=F]
Ms.x= Ms.x[,c("Tum", "cfDNA3", "N1a", "N2a")]
colSums(Ms.x>0)
Stacked_histo(Ms.x)[1] #6x4


## --- [FIGURE SUPP 7A] Visualisation Mice with cfDNA + Needles and combinasion Circle plots----

Fig4.select.expamples = df.Fig4 %>% filter(Names=="Tum"|Names=="Needles1"|Names=="cfDNA3"|Names=="Needles1.cf") %>%
  group_by(Model,Mouse,Names) %>% dplyr::summarise(n=mean(Biom_Tum/10000), n1=mean(Biom_Tum1/10000))

ms.cfDNA.circle = Fig4.select.expamples %>% filter(Names == "cfDNA3" & n > 0) %>% .$Mouse

layout(matrix(1:30, 6, 5))
for (i in unique(ms.cfDNA.circle)) {
  #i=363
  print(i)
  mini.df= Fig4.select.expamples %>% filter(Mouse==i)
  mini.df = mini.df %>% arrange(Names)
  category = mini.df$Names
  percent = round(mini.df$n1, digits = 1)
  color = viridis(4, direction = -1)
  
  circos.par("start.degree" = 90, cell.padding = c(0, 0, 0, 0))
  circos.initialize("a", xlim = c(0, 100)) # 'a` just means there is one sector
  circos.track(ylim = c(0.5, length(percent)+0.5), track.height = 0.8, 
               bg.border = NA, panel.fun = function(x, y) {
                 xlim = CELL_META$xlim
                 circos.segments(rep(xlim[1], 4), 1:4, rep(xlim[2], 4), 1:4, col = "#CCCCCC")
                 circos.rect(rep(0, 4), 1:4 - 0.45, percent, 1:4 + 0.45, col = color, border = "white")
                 circos.text(rep(xlim[1], 4), 1:4,paste(category, " - ", percent, "%"), 
                             facing = "downward", adj = c(1.05, 0.5), cex = 0.8) 
                 breaks = seq(0, 100, by = 5)
                 circos.axis(h = "top", major.at = breaks, labels = paste0(breaks, "%"),labels.cex = 0.6)
               })
  title(paste(mini.df$Model[1],"_", mini.df$Mouse[1], sep=""))
  circos.clear()
}#10x10


## --- [FIGURE SUPP 6D] Biomass captured with linked point in combination ----

df.Fig.supp6.d = data.frame()

for (i in unique(df.Fig1.clean$Mouse)) {
  print(i)
  TEST.name = paste("_", i, "_", sep="")
  
  info.1 = info %>% filter(Clever_cutting==TRUE)
  Ms.info = info.1 %>% filter(`Mice #`==i)
  
  Ms.x = df.PCR.Tum.Filt.I[,grep(TEST.name, names(df.PCR.Tum.Filt.I)), drop=FALSE]
  if(ncol(Ms.x)==0){No_Datum = c(No_Datum,i);next}
  
  id.name <- names(Ms.x[1])
  
  colnames(Ms.x) <- str_split(names(Ms.x), "_", simplify = T)[,3]
  if(is.na(sum(Ms.x[,"Tum"]))==T){NA_tum = c(NA_tum,i);next}
  
  Ms.x = Ms.x[,!is.na(colSums(Ms.x)), drop=F]
  
  if(sum(names(Ms.x) %in% "cfDNA3")==0 & sum(names(Ms.x) %in% c("N1a", "N1b", "N2a", "N2b","N3a", "N3b", "N4a", "N4b"))==0){next}
  if(sum(names(Ms.x) %in% "cfDNA3")==0 | sum(names(Ms.x) %in% c("N1a", "N1b", "N2a", "N2b","N3a", "N3b", "N4a", "N4b"))==0){Sampl="one";Ms.x$cfDNA3 = 0} else {Sampl="both"}
  
  Ms.x.1 = subset( Ms.x, select = -Tum )
  Ms.x.1[Ms.x.1<10000]=0 # FILTER BARCODE UNDER 1% in needles
  Ms.x.2 = cbind(Ms.x[,"Tum", drop=F],Ms.x.1)
  
  if(sum(names(Ms.x) %in% c("N1a","N2a","N3a","N4a"))>=2){
    combo = sample(names(Ms.x)[names(Ms.x) %in% c("N1a","N2a","N3a","N4a")],2)
    nee.combo = as.data.frame(Ms.x[,combo[1]]+Ms.x[,combo[2]])
    names(nee.combo) = "Na.combo"
    Ms.x = cbind(Ms.x,nee.combo)
  }
  
  if(sum(names(Ms.x) %in% c("N1b","N2b","N3b","N4b"))>=2){
    combo = sample(names(Ms.x)[names(Ms.x) %in% c("N1b","N2b","N3b","N4b")],2)
    nee.combo.b = as.data.frame(Ms.x[,combo[1]]+Ms.x[,combo[2]])
    names(nee.combo.b) = "Nb.combo"
    Ms.x = cbind(Ms.x,nee.combo.b)
  }
  
  if(sum(names(Ms.x.2) %in% c("N1a","N2a","N3a","N4a"))>=2){
    combo = sample(names(Ms.x.2)[names(Ms.x.2) %in% c("N1a","N2a","N3a","N4a")],2)
    nee.combo = as.data.frame(Ms.x.2[,combo[1]]+Ms.x.2[,combo[2]])
    names(nee.combo) = "Na.combo"
    Ms.x.2 = cbind(Ms.x.2,nee.combo)
  }
  
  if(sum(names(Ms.x.2) %in% c("N1b","N2b","N3b","N4b"))>=2){
    combo = sample(names(Ms.x.2)[names(Ms.x.2) %in% c("N1b","N2b","N3b","N4b")],2)
    nee.combo.b = as.data.frame(Ms.x.2[,combo[1]]+Ms.x.2[,combo[2]])
    names(nee.combo.b) = "Nb.combo"
    Ms.x.2 = cbind(Ms.x.2,nee.combo.b)
  }
  
  if(sum(names(Ms.x) %in% "cfDNA3")==1){
    
    test= Ms.x[,(colnames(Ms.x) %in% c("N1a", "N1b", "N2a", "N2b","N3a", "N3b", "N4a", "N4b"))]+Ms.x[,"cfDNA3"]
    names(test) = paste(names(test),"cf", sep = ".")
    Ms.x = cbind(Ms.x, test)
    if(sum(names(Ms.x) %in% c("CTC","Blood"))==1){
      test.2= test[,(colnames(test) %in% c("N1a.cf", "N1b.cf", "N2a.cf", "N2b.cf","N3a.cf", "N3b.cf", "N4a.cf", "N4b.cf"))]+Ms.x[,(colnames(Ms.x) %in% c("CTC","Blood"))]
      names(test.2) = paste(names(test.2),"CTC", sep = ".")
      Ms.x = cbind(Ms.x, test.2)
    }
  }
  
  if(sum(names(Ms.x.2) %in% "cfDNA3")==1){
    
    test= Ms.x.2[,(colnames(Ms.x.2) %in% c("N1a", "N1b", "N2a", "N2b","N3a", "N3b", "N4a", "N4b"))]+Ms.x.2[,"cfDNA3"]
    names(test) = paste(names(test),"cf", sep = ".")
    Ms.x.2 = cbind(Ms.x.2, test)
    if(sum(names(Ms.x) %in% c("CTC","Blood"))==1){
      test.2= test[,(colnames(test) %in% c("N1a.cf", "N1b.cf", "N2a.cf", "N2b.cf","N3a.cf", "N3b.cf", "N4a.cf", "N4b.cf"))]+Ms.x[,(colnames(Ms.x) %in% c("CTC","Blood"))]
      names(test.2) = paste(names(test.2),"CTC", sep = ".")
      Ms.x.2 = cbind(Ms.x.2, test.2)
    }
  }
  
  e = biomass_fun(Ms.x)
  e1 = biomass_fun(Ms.x.2)
  
  Ms.x.df = data.frame(nb_bc = colSums(Ms.x>0),
                       pct_bc = colSums(Ms.x>0)*100/colSums(Ms.x[,"Tum", drop=F]>0),
                       Names= colnames(Ms.x),
                       Shanon=diversity(t(Ms.x)),
                       Biom_Tum = e[,"Tum"],
                       Corr_Tum = cor(log10(Ms.x+1))[,1],
                       Corr_Tum2 = cor(log2(Ms.x+1))[,1],
                       Corr_Tum3 = cor(Ms.x)[,1],
                       Corr_Tum4 = cor(log10(Ms.x+1), use = "pairwise.complete.obs")[,1],
                       Biom_Tum1 = e1[,"Tum"],
                       #Biom_Tum10 = e2[,"Tum"],
                       Mouse = i,
                       Exp =Ms.info$Exp,
                       Model =Ms.info$Cells,
                       Sampling = Sampl)
  df.Fig.supp6.d  = rbind(df.Fig.supp6.d ,Ms.x.df)
}

df.Fig.supp6.d$Raw.Names = df.Fig.supp6.d$Names
df.Fig.supp6.d$Names = str_replace(df.Fig.supp6.d$Names, c("N1a|N2a|N3a|N4a"), "Needles1")
df.Fig.supp6.d$Names = str_replace(df.Fig.supp6.d$Names, c("N1b|N2b|N3b|N4b"), "Needles2")

df.Fig.supp6.d$Names = str_replace(df.Fig.supp6.d$Names, c("N1a.cf|N2a.cf|N3a.cf|N4a.cf"), "Needles1.cf")
df.Fig.supp6.d$Names = str_replace(df.Fig.supp6.d$Names, c("N1b.cf|N2b.cf|N3b.cf|N4b.cf"), "Needles2.cf")

df.Fig.supp6.d$Names = str_replace(df.Fig.supp6.d$Names, c("N1a.cf.CTC|N2a.cf.CTC|N3a.cf.CTC|N4a.cf.CTC"), "Needles1.cf.CTC")
df.Fig.supp6.d$Names = str_replace(df.Fig.supp6.d$Names, c("N1b.cf.CTC|N2b.cf.CTC|N3b.cf.CTC|N4b.cf.CTC"), "Needles2.cf.CTC")

df.Fig.supp6.d$Names3 = str_split(df.Fig.supp6.d$Raw.Names, "\\.", simplify = T)[,1]

df.Fig.supp6.d$Mouse =as.numeric(df.Fig.supp6.d$Mouse)
df.Fig.supp6.d$Exp = as.numeric(df.Fig.supp6.d$Exp)


df.Fig4.2 = df.Fig.supp6.d %>% 
  filter(Names %in% c("Needles1","cfDNA3","Needles2","Needles1.cf","Needles2.cf")) %>%
  arrange(Exp, Mouse, Raw.Names)


v= df.Fig4.2 %>% .$Names3

paired.N = 9999:(9998+length(v))

for (i in 1:length(v)) {#loop to pair multiple needle samples with matching cfDNA
  print(i)
  #i=3
  if(v[i] == v[i+1]){
    paired.N[i]  = i
    paired.N[i+1] = i
  }
}


ms.cfDNA.linked = df.Fig4.2 %>% filter(Names == "cfDNA3" & Biom_Tum1 > 0) %>% .$Mouse

df.Fig4.2$Paired.spl = paired.N

df.Fig4.2 = df.Fig4.2 %>% group_by(Model,Mouse,Names3) %>% 
  dplyr::mutate(dif = (max(Biom_Tum1, na.rm = T)-min(Biom_Tum1, na.rm = T))/10000)%>% ungroup()

# Stats needles 1
df.Fig4.2.stat = df.Fig4.2 %>% filter( Mouse %in% ms.cfDNA.linked) %>%
  filter(Model=="PDX-T412"|Model=="MDA-231"|Model=="PDX-1432C") %>%
  filter(Names %in% c("Needles1", "Needles1.cf"))

compare_means(Biom_Tum1 ~ Names,
              group.by = "Model",
              data = df.Fig4.2.stat,
              p.adjust.method = "holm",
              method = "t.test",
              paired = TRUE) %>%
  mutate(p.adj = format.pval(p.adj, digits = 2))
# Stats needles 2
df.Fig4.2.stat = df.Fig4.2 %>% filter( Mouse %in% ms.cfDNA.linked) %>%
  filter(Model=="PDX-T412"|Model=="MDA-231"|Model=="PDX-1432C") %>%
  filter(Names %in% c("Needles2", "Needles2.cf"))

compare_means(Biom_Tum1 ~ Names,
              group.by = "Model",
              data = df.Fig4.2.stat,
              p.adjust.method = "holm",
              method = "t.test",
              paired = TRUE) %>%
  mutate(p.adj = format.pval(p.adj, digits = 2))

# Plot in figure
df.Fig4.2 %>% filter( Mouse %in% ms.cfDNA.linked) %>%
  filter(Model=="PDX-T412"|Model=="MDA-231"|Model=="PDX-1432C") %>%
  ggplot(aes(x=Names, y=Biom_Tum1))+
  geom_bar(aes(fill=Names),stat = "summary", fun=mean)+
  geom_line(aes(group=Paired.spl, colour=dif), size=1.2)+
  geom_point(colour="grey20")+
  facet_wrap(~Model)+
  theme_classic()+
  scale_fill_grey()+
  scale_color_viridis()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("")+ylab("Tumour biomass")+labs(title = "Primary tumour Biomass captured with combination (1%)")#6x11


## --- [FIGURE 4G and SUPP 6B] Scatter plot Needle/cfDNA ----
df.Fig4.scatter=data.frame()
for (i in unique(df.Fig1.clean$Mouse)) {
  print(i)
  TEST.name = paste("_", i, "_", sep="")
  
  info.1 = info %>% filter(Clever_cutting==TRUE)
  Ms.info = info.1 %>% filter(`Mice #`==i)
  
  Ms.x = df.PCR.Tum.Filt.I[,grep(TEST.name, names(df.PCR.Tum.Filt.I)), drop= FALSE]
  if(ncol(Ms.x)==0){No_Datum = c(No_Datum,i);next}
  
  id.name <- names(Ms.x[1])
  
  colnames(Ms.x) <- str_split(names(Ms.x), "_", simplify = T)[,3]
  if(is.na(sum(Ms.x[,"Tum"]))==T){NA_tum = c(NA_tum,i);next}
  
  Ms.x = Ms.x[,!is.na(colSums(Ms.x)), drop=F]
  
  if(sum(names(Ms.x) %in% "cfDNA3")==0 & sum(names(Ms.x) %in% c("N1a", "N1b", "N2a", "N2b","N3a", "N3b", "N4a", "N4b"))==0){next}
  #if(sum(names(Ms.x) %in% "cfDNA3")==0 | sum(names(Ms.x) %in% c("N1a", "N1b", "N2a", "N2b","N3a", "N3b", "N4a", "N4b"))==0){Sampl="one";Ms.x$cfDNA3 = 0} else {Sampl="both"}
  
  Ms.x.temp = Ms.x[rowSums(Ms.x)>0,]
  if(sum(colnames(Ms.x.temp) %in% c("N1a", "N2a"))<2){next}
  if(sum(names(Ms.x.temp) %in% "cfDNA3")==0){next}
  
  Ms.x.scatter = Ms.x.temp[,c("N1a", "N2a", "cfDNA3", "Tum")]
  Ms.x.scatter$mouse = i
  Ms.x.scatter$Model = Ms.info$Cells
  Ms.x.scatter$Exp = Ms.info$Exp
  
  df.Fig4.scatter = rbind(df.Fig4.scatter, Ms.x.scatter)
  
}

n.Fig4.g = df.Fig4.scatter %>% select(Model, Exp, mouse) %>% group_by(Model) %>% 
  dplyr::mutate(count.Ms = n_distinct(mouse), count.Exp = n_distinct(Exp)) %>%
  ungroup() %>% select(Model, count.Ms, count.Exp) %>% unique()

n.Ms.list[[length(n.Ms.list)+1]] = n.Fig4.g
names(n.Ms.list)[length(n.Ms.list)] = "Figure_4g"

## --- [FIGURE 4G] Plot ----

ggplot()+
  geom_point(data= subset(df.Fig4.scatter,cfDNA3>0&N1a>0),aes(x=cfDNA3, y=N1a, size=Tum), colour="darkorange3", alpha=0.9)+
  geom_quasirandom(data= subset(df.Fig4.scatter,cfDNA3==0&N1a>0),aes(x=10, y=N1a, size=Tum),width=0.2,groupOnX=TRUE,alpha=0.7)+
  geom_quasirandom(data= subset(df.Fig4.scatter,cfDNA3>0&N1a==0),aes(x=cfDNA3, y=10, size=Tum),width=0.2,groupOnX=F,alpha=0.7)+
  scale_x_log10()+scale_y_log10()+
  facet_wrap(~Model, nrow = 1)+
  theme_classic()+labs(title = "Scatter plot cfDNA3 vs Needle1a (size = Tum)")#3.5x14

## --- [FIGURE SUPP 6B] Plot ----

mD = melt(df.Fig4.scatter, id=c("Tum", "mouse", "Model", "Exp"))
mD[mD==0] = NA
mD$value[mD$value>0]=0.1

mD %>% filter(Model!= "MDA-468") %>% ggplot()+
  geom_quasirandom(aes(x=variable, y=Tum, size=value, colour=variable), alpha=0.8)+
  scale_y_log10()+scale_size_continuous(range = c(1,1.1))+
  scale_color_brewer(palette ="Dark2", direction = -1)+
  facet_wrap(~Model, nrow = 1)+
  labs(title = "Barcode presence in Biopsies according to Tumour frq")+
  theme_classic()

## --- [FIGURE SUPP 6E]  Clone splitting: Tum growth and collection: ----

plot.6 = info %>% select(Exp,`Mice #`,`Days since injection`,`Volume of Primary (mm3)`)
plot.6 = plot.6 %>% filter(Exp=="80A"|Exp=="80B")
names(plot.6) = c("Exp", "Ms", "Days", "Vol")
plot.6$Vol = as.numeric(plot.6$Vol)

ggplot(plot.6)+
  geom_point(aes(x=Days, y=Vol, colour=Exp))+
  geom_text_repel(aes(x=Days, y=Vol, label=Ms),max.overlaps = 80)+
  theme_classic()

## --- [FIGURE SUPP 6E]  cfDNA in clone splitting experiment ----

cs.info = read_xlsx("Clone_splitting_info.xlsx")

cs.info = cs.info %>% select("Timepoint", "Exp", 'Mice #', "Project", 'Injected cells', "Cells", "MI", 
                             "Cln.Spl", "Barcode",'Terminal Endbleed volume', "Protocol")

names(cs.info) = c("Timepoint", "Exp", 'Mice#', "Project", 'Injected cells', "Cells", "MI", "Cln.Spl", "Barcode",
                   'blood_volume', "meth") 

cs.info = cs.info %>% select("Timepoint", "Exp", 'Mice#', "Cells","blood_volume")
names(cs.info) = c("Timepoint", "Exp2", "Ms", "Model", "blood_volume")

cs.info$Ms = as.character(cs.info$Ms)

cfDNA.filt.cs = cfDNA.filt.cs %>% filter(Exp2 %in% c("80A", "80B") & m>0) %>% left_join(cs.info, by = c("Ms", "Exp2", "Model"))

ggplot()+
  geom_tile(data=cs.info, aes(x=Timepoint, y=as.factor(Ms), colour=blood_volume>0, colour="grey"), fill = "grey90")+
  geom_tile(data=cfDNA.filt.cs, aes(x=Timepoint, y=Ms),colour = "black", fill="darkgreen", alpha=0.5)+
  facet_wrap(Exp2~Timepoint, scales = "free")

## --- [FIGURE SUPP 6F] clone splitting experiments, Histograms with ctrl and correlation heatmaps : ----
df.clone.split = df.PCR[,grep(paste(c("Exp80_"), collapse = "|"), names(df.PCR))]
df.clone.split.Tum = df.PCR.Tum[,grep(paste(c("Exp80_"), collapse = "|"), names(df.PCR.Tum))]

df.Fig1.Cln = Figure1(info, df.clone.split.Tum, df.clone.split, weight)[[1]]

Ms80A = df.Fig1.Cln %>% filter(Exp=="80A") %>% .$Mouse %>% unique()
Ms80A = paste0("_", Ms80A, "_")

Exp80A = df.clone.split.Tum[,grep(paste(Ms80A, collapse = "|"), names(df.clone.split.Tum))]
Exp80A = Exp80A[,grep(paste(c("_Tum", "cfDNA"), collapse = "|"), names(Exp80A))]

Control = df.clone.split.Tum[,grep("99_Ctrl", names(df.clone.split.Tum))]
Exp80A$Ctrl <- rowMeans(Control[,c(1,2)], na.rm=TRUE)
Exp80A = Exp80A[,!is.na(colSums(Exp80A))]
Exp80A = Exp80A %>% select(Ctrl,
                           Exp80_1_Tum, Exp80_2_Tum, Exp80_3_Tum, Exp80_4_Tum, #T1 Tumours
                           Exp80_5_Tum, Exp80_6_Tum, Exp80_7_Tum, Exp80_8_Tum, #T2 Tumours
                           Exp80_9_Tum, Exp80_10_Tum, Exp80_11_Tum, Exp80_12_Tum, #T3 Tumours
                           Exp80_5_cfDNA2_, Exp80_6_cfDNA2_, Exp80_8_cfDNA2_, #cfDNA2
                           Exp80_9_cfDNA3_, Exp80_10_cfDNA3_, Exp80_11_cfDNA3_, Exp80_12_cfDNA3_) # cfDNA3

Stacked_histo(Exp80A)[1]#6x15

Exp80A = Exp80A %>% select(Exp80_5_Tum,
                           Exp80_6_Tum,
                           Exp80_8_Tum,
                           Exp80_9_Tum,
                           Exp80_10_Tum,
                           Exp80_11_Tum,
                           Exp80_12_Tum,
                           Exp80_5_cfDNA2_,
                           Exp80_6_cfDNA2_,
                           Exp80_8_cfDNA2_,
                           Exp80_9_cfDNA3_,
                           Exp80_10_cfDNA3_,
                           Exp80_11_cfDNA3_,
                           Exp80_12_cfDNA3_)

names(Exp80A) = paste(str_split(names(Exp80A), "_", simplify = T)[,2], str_split(names(Exp80A), "_", simplify = T)[,3], sep = "_")
ComplexHeatmap::Heatmap(cor(Exp80A), col=viridis(256))#6x8



Ms80B = df.Fig1.Cln %>% filter(Exp=="80B") %>% .$Mouse %>% unique()
Ms80B = paste0("_", Ms80B, "_")
Exp80B = df.clone.split.Tum[,grep(paste(Ms80B, collapse = "|"), names(df.clone.split.Tum))]
Exp80B = Exp80B[,grep(paste(c("_Tum", "cfDNA"), collapse = "|"), names(Exp80B))]

Control = df.clone.split.Tum[,grep("100_Ctrl", names(df.clone.split.Tum))]
Exp80B$Ctrl <- rowMeans(Control[,c(1,2)], na.rm=TRUE)
Exp80B = Exp80B[,!is.na(colSums(Exp80B))]
Exp80B = Exp80B %>% select(Ctrl,
                           Exp80_13_Tum, Exp80_14_Tum, Exp80_15_Tum, Exp80_16_Tum,
                           Exp80_17_Tum, Exp80_18_Tum, Exp80_19_Tum, Exp80_20_Tum,
                           Exp80_21_Tum, Exp80_22_Tum, Exp80_23_Tum, Exp80_24_Tum,
                           Exp80_18_cfDNA2_,Exp80_20_cfDNA2_,
                           Exp80_21_cfDNA3_,Exp80_22_cfDNA3_,Exp80_23_cfDNA3_,Exp80_24_cfDNA3_)

Stacked_histo(Exp80B)[1]#6x15

Exp80B = Exp80B %>% select(Exp80_18_Tum,
                           Exp80_20_Tum,
                           Exp80_21_Tum,
                           Exp80_22_Tum,
                           Exp80_23_Tum,
                           Exp80_24_Tum,
                           Exp80_18_cfDNA2_,
                           Exp80_20_cfDNA2_,
                           Exp80_21_cfDNA3_,
                           Exp80_22_cfDNA3_,
                           Exp80_23_cfDNA3_,
                           Exp80_24_cfDNA3_)

names(Exp80B) = paste(str_split(names(Exp80B), "_", simplify = T)[,2], str_split(names(Exp80B), "_", simplify = T)[,3], sep = "_")
ComplexHeatmap::Heatmap(cor(Exp80B), col=viridis(256))

