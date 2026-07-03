###Analisi di network PHQ-9 (GGM) e RFQ as moderator (MNM)####

### Library ###
library(bootnet) 
library(qgraph)
library(networktools)
library(mgm)
library(tidyverse)
library(parallel)
library(EGAnet)       
library(igraph)   

?bootnet
?qgraph
?networktools
?mgm
?tidyverse
?parallel
?FactorGraph
?netSimulator
?EGA
?`EGAnet-package`

###Datasets x GGM e MNM####
dual_trauma <- read_csv("Data/dual_trauma.csv")
phq_items   <- paste0("phq0", 1:9)
data_ggm <- dual_trauma %>%
  select(all_of(phq_items)) %>%
  filter(complete.cases(.))  #969 part (female=481, male=479, prefer not to say=9)


#dataset per MNM con RFQ_tot 
data_rfqtot <- dual_trauma %>% select(rfq_tot, all_of(phq_items)) %>% filter(complete.cases(.))
#mean-centered data
m_rfqtot <- mean(data_rfqtot$rfq_tot)

#dataset per MNM con RFQ02 
data_rfq02 <- dual_trauma %>% select(rfq02, all_of(phq_items)) %>% filter(complete.cases(.))
data_rfq05 <- dual_trauma %>% select(rfq05, all_of(phq_items)) %>%  filter(complete.cases(.))
data_rfq07 <- dual_trauma %>% select(rfq07, all_of(phq_items)) %>%  filter(complete.cases(.))
data_rfq12 <- dual_trauma %>% select(rfq12, all_of(phq_items)) %>%  filter(complete.cases(.))
data_rfq13 <- dual_trauma %>% select(rfq13, all_of(phq_items)) %>%  filter(complete.cases(.))


#mean-centered data
m_rfq02 <- mean(data_rfq02$rfq02)

###Exploratory Graph Analysis of the phq-9####

ega_result <- EGA(
  data      = data_ggm,
  model     = "glasso",       # usa graphical LASSO (come EBICglasso)
  algorithm = "walktrap",     # algoritmo di community detection
  plot.EGA  = TRUE            # mostra il plot
)

boot_ega <- bootEGA(
  data      = data_ggm,
  model     = "glasso",
  algorithm = "walktrap",
  iter      = 1000,
  type      = "parametric",
  ncores    = max(1, parallel::detectCores() - 1)
)

#unidimensional--> usare solo strength e EI (non Bridge)




###GGM PHQ-9####
colours_phq <- rep("lightblue", 9)
node_labels <- paste0("PHQ0", 1:9)

network_phq <- estimateNetwork(data_ggm, default = "EBICglasso")

plot_network_phq <- qgraph(network_phq$graph, layout = "spring", color = colours_phq,
                           labels = node_labels,
                           edge.color = "#1565C0",
                           width = 7 * 1.4,
                           height = 7,
                           border.width = 0.5,
                           label.cex = 0.7,
                           label.color = "black",
                           label.prop = 0.9,
                           fade = T,
                           curveAll = T,
                           curveDefault = 0.5,
                           edge.labels = T,
                           font = 2)

qgraph(plot_network_phq, filetype='jpeg', dpi = 300)


###stabilità edges####
btnt_np <- bootnet(network_phq,
                   type = "nonparametric",
                   statistics = "edge",
                   nBoots = 1000,
                   verbose = FALSE,
                   nCores = 1)

plot(btnt_np,
     order = "sample")

plot(btnt_np, 
     plot = "difference",
     order = "sample")

###network stability CS-coefficent#### 
b1 <- bootnet(network_phq, boots = 1000, nCores = 4, statistics = c("strength", "expectedInfluence"), nBoots = 100, type = "case")

corStability(b1)

plot(b1, statistics = c("strength", "expectedInfluence"))
ggsave(filename = "stability_plot.jpeg", dpi = 300)

###centrality plots###
c1 <- centralityPlot(plot_network_phq, include = c("Strength", "ExpectedInfluence"), orderBy = "Strength")

Sim1 <- netSimulator(
  input = network_phq, 
  default = "EBICglasso",
  nCases = seq(200, 1000, 200),
  nReps = 100)

plot(Sim1)


myedges <- getWmat(network_phq)
write.csv(myedges, "results/mynetworkedges.csv")

#conta tutti gli edges
sum(myedges != 0)


###MODERATED NETWORK ANALYSIS####

#MOD with RFQ_TOT#

mgm_rfq_tot <- mgm(data = as.matrix(data_rfqtot),
                   type = rep("g", 10),
                   level = rep(1,10),
                   moderators = 1,
                   lambdaSel = "EBIC",
                   lambdaGam = 0.5,
                   ruleReg = "AND",
                   scale = TRUE,
                   verbatim = T)


mgm_rfq_tot$interactions$indicator
showInteraction(mgm_rfq_tot, int = c(1,7,10)) #0.103 positive

#visualizzazione factor Graph#
FactorGraph(mgm_rfq_tot, labels = c("rfq_tot", node_labels), 
            layout = "spring",
            shapeSizes = c(6,3),
            colors =   c("lightblue", "lightgreen", "red"))

FactorGraph(mgm_rfq_tot, labels = c("rfq_tot", node_labels), 
            layout = "spring",
            Nodewise = TRUE,
            shapeSizes = c(6,3),
            edge.labels = TRUE,
            colors =   c("lightblue", "lightgreen", "red"))




            
#MOD with singular items RFQ#

mgm_rfq02 <- mgm(data = as.matrix(data_rfq02),
                   type = rep("g", 10),
                   level = rep(1,10),
                   moderators = 1,
                   lambdaSel = "EBIC",
                   lambdaGam = 0.5,
                   ruleReg = "AND",
                   scale = TRUE,
                   verbatim = T)

mgm_rfq05 <- mgm(data = as.matrix(data_rfq05),
                   type = rep("g", 10),
                   level = rep(1,10),
                   moderators = 1,
                   lambdaSel = "EBIC",
                   lambdaGam = 0.5,
                   ruleReg = "AND",
                   scale = TRUE,
                   verbatim = T)
            
mgm_rfq07 <- mgm(data = as.matrix(data_rfq07),
                   type = rep("g", 10),
                   level = rep(1,10),
                   moderators = 1,
                   lambdaSel = "EBIC",
                   lambdaGam = 0.5,
                   ruleReg = "AND",
                   scale = TRUE,
                   verbatim = T)          
mgm_rfq12 <- mgm(data = as.matrix(data_rfq12),
                 type = rep("g", 10),
                 level = rep(1,10),
                 moderators = 1,
                 lambdaSel = "EBIC",
                 lambdaGam = 0.5,
                 ruleReg = "AND",
                 scale = TRUE,
                 verbatim = T)          
mgm_rfq13 <- mgm(data = as.matrix(data_rfq13),
                 type = rep("g", 10),
                 level = rep(1,10),
                 moderators = 1,
                 lambdaSel = "EBIC",
                 lambdaGam = 0.5,
                 ruleReg = "AND",
                 scale = TRUE,
                 verbatim = T)  


mgm_rfq02$interactions$indicator
mgm_rfq05$interactions$indicator
mgm_rfq07$interactions$indicator
mgm_rfq12$interactions$indicator
mgm_rfq13$interactions$indicator

#Unica tripla selezionata da RFQ02#
showInteraction(mgm_rfq02, int = c(1,7,10)) #0.118 positive
mgm_rfq02$pairwise$wadj[1, 7]  #peso diretto tra rfq02 e phq06 (0.05)
mgm_rfq02$pairwise$wadj[1, 10]  #nessun arco pairwise tra rfq02 e phq09
mgm_rfq02$rawlambda      # pesi grezzi per ogni regressione nodewise
mgm_rfq02$interactions$weights[[1]]  
mgm_rfq02$interactions$weights[[2]]
#factorgraph RFQ02 moderatore#
FactorGraph(mgm_rfq02, labels = c("rfq02", node_labels), 
            layout = "spring",
            shapeSizes = c(6,3),
            colors =   c("lightblue", "lightgreen", "red"))

FactorGraph(mgm_rfq02, labels = c("rfq02", node_labels), 
            layout = "spring",
            Nodewise = TRUE,
            shapeSizes = c(6,3),
            edge.labels = TRUE,
            colors =   c("lightblue", "lightgreen", "red"))

###Reti CONDIZIONATE con vari livelli di RFQ02####

m_rfq02 <- mean(data_rfq02$rfq02)
vals_orig <- c(1, 3, 4)
vals_cent <- vals_orig - m_rfq02 # m_rfq02 calcolato prima

l_mgm_cond <- list()
for (k in 1:3) {
  l_mgm_cond[[k]] <- condition(mgm_rfq02,
                               values = list("1" = vals_cent[k]))
}


# Step 2: scala comune (maximum)

v_max <- sapply(1:3, function(k) max(abs(l_mgm_cond[[k]]$pairwise$wadj)))
max_w <- max(v_max)


# Step 3: layout comune (stessi nodi, stessa posizione)

layout_comune <- averageLayout(
  l_mgm_cond[[1]]$pairwise$wadj,
  l_mgm_cond[[2]]$pairwise$wadj,
  l_mgm_cond[[3]]$pairwise$wadj
)


# Step 4: etichette e colori archi
# PHQ06 (col 7) - PHQ09 (col 10) in rosso

etichette <- c("Mai(nessun abuso)",
               "Qualche volta",
               "Spesso (abuso frequente)")

n <- 9
edge_col <- matrix("#1565C0", n, n)
edge_col[6, 9] <- "#D62728"
edge_col[9, 6] <- "#D62728"


# Step 5: plot 3 reti affiancate

par(mfrow = c(1, 3))

for (k in 1:3) {
  # estrai solo i nodi PHQ (rimuovi riga/col del moderatore)
  wadj_k <- l_mgm_cond[[k]]$pairwise$wadj[2:10, 2:10]
  lay_k  <- layout_comune[2:10, ]
  
  peso_target <- round(wadj_k[6, 9], 2)
  
  qgraph(wadj_k,
         layout         = lay_k,
         maximum        = max_w,
         color          = rep("lightblue", 9),
         labels         = node_labels,
         vsize          = 14,
         edge.color     = edge_col,
         width          = 7 * 1.4,
         height         = 7,
         border.width   = 1.5,
         label.cex      = 0.7,
         label.prop     = 0.9,
         edge.labels    = TRUE,
         edge.label.cex = 1.5,
         edge.label.bg  = TRUE,
         cut            = 0.05,
         fade           = TRUE,
         curveAll       = TRUE,
         curveDefault   = 0.5,
         mar            = c(3, 3, 5, 3))
         
}

par(mfrow = c(1, 1))




?resample
###RESAMPLE di rfq_tot e rfq02####
res_rfq_tot <- resample(object = mgm_rfq_tot,
                        data = as.matrix(data_rfqtot),
                        nB = 500)

hist(res_rfq_tot$bootParameters[7,10, ],
     main = "",
     xlab = "Parameter Estimate")


plotRes(object = res_rfq_tot,
        quantiles = c(0.05, .95))



res_rfq02 <- resample(object = mgm_rfq02,
                        data = as.matrix(data_rfq02),
                        nB = 500)

png("hist_PHQ06_PHQ09.png", width = 800, height = 600)
hist(res_rfq02$bootParameters[7, 10, ],
     main = "Bootstrap PHQ06-PHQ09",
     xlab = "Parameter Estimate")
dev.off()


plotRes(object = res_rfq02,
        quantiles = c(0.05, .95))



