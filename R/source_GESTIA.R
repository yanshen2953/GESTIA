######################souce code of GESTIA
#' @title calculate the relative positions of genes in a pathway
#' @param id.vertex A vector containing the names of the genes
#' @param tmp.g An igraph object which was constructed from the global network by deleting the unrelated genes of the pathway
#' @description Calculate the ratio of upstream/downstream genes in a pathway. For scattered pathway/gene set, the ratio is weighted by the proportion of the sub-graph.
#' @return A dataframe consisted of two variables, the first is the proportion of downstream genes, the second is the proportion of upstream genes.
#' @export
relativePosition<-function(id.vertex,tmp.g){
  tmp.d<-c()
  tmp.u<-c()
  for(i in id.vertex){
    dist.tmp.d <- distances(tmp.g, v=V(tmp.g)[name==i], to=V(tmp.g), weights=NA,mode="out") #mode=out means downstream
    dist.tmp.u <- distances(tmp.g, v=V(tmp.g)[name==i], to=V(tmp.g), weights=NA,mode="in") #upstream
    dist.tmp.d[is.infinite(dist.tmp.d)]<- -1
    dist.tmp.u[is.infinite(dist.tmp.u)]<- -1
    n.d<-length(dist.tmp.d[dist.tmp.d>0])
    n.u<-length(dist.tmp.u[dist.tmp.u>0])
    n.o<-length(dist.tmp.d)
    if(n.d + n.u == 0){ #this is when the gene is totally isolated
      tmp.d<-c(tmp.d,1/n.o)
      tmp.u<-c(tmp.u,1/n.o)
    }else{
      tmp.d<-c(tmp.d,n.d/(n.u+n.d)*(n.u+n.d+1)/n.o) #n.d/(n.u+n.d) is the relative position of this gene in the sub-network
      #(n.u+n.d+1)/n.o is the proportion of this sub-network to the whole network of the gene set
      #if the genes of this gene set are all connected, then (n.u+n.d+1)/n.o = 1
      #this will be used in the calculation of the impact score of genes outside this gene set
      tmp.u<-c(tmp.u,n.u/(n.u+n.d)*(n.u+n.d+1)/n.o) #this will be used as coefficient in the calculation of the genes inside this gene set
    }
  }
  names(tmp.d)<-id.vertex
  names(tmp.u)<-id.vertex
  return(data.frame(tmp.d,tmp.u))
}

#' @title Calculate impact of pathway A over pathway B
#' @param ig.all Global network constructed from KEGG
#' @param degree.a The degree of genes in pathway A
#' @param degree.b The degree of genes in pathway B
#' @param rp.a.d The proportion of downstream genes of the genes in pathway A
#' @param rp.b.d The proportion of downstream genes of the genes in pathway B
#' @param rp.a.u The proportion of upstream genes of the genes in pathway A
#' @param rp.b.u The proportion of upstream genes of the genes in pathway B
#' @param penal.isolated The penalization factor for the isolated genes. Default is 0.1
#' @param cutImpactSub The cutoff for impact score. Impact score smaller than this cutoff will be set to 0. Default 0.
#' @param plot Whether to plot the splitted network of the two pathways. TRUE or FALSE. Default FALSE.
#' @param n.permu Number of permutations for calculating the empirical null distribution. Default 30.
#' @param minGeneNum The cutoff for number of genes in pathways. Pathways with less genes (after id conversion) than this cutoff will be disgarded.
#' @param exportData Whether or not to export more details in the calculation. Default FALSE.
#' @description This function calls the core of the calculation function, and calculate the raw and permutated impact scores, then calculate p-value and GESTIA scores.
#' @return List of results including: p-value, Normalized score,ScoreA before normalization,ScoreB before normalization,Mean permutation A,Mean permutation B
#' @export
relativeImpact2<-function(ig.all,degree.a,degree.b,rp.a.d,rp.b.d,rp.a.u,rp.b.u,
                          penal.isolated=0.1,cutImpactSub=0,plot=FALSE,
                          n.permu=30,minGeneNum=5,exportData=F){ #relative position ->rp
  symbol.a<-names(degree.a)
  symbol.b<-names(degree.b)
  v.all<-V(ig.all)#get the vertice names of the total graph, here equal to gene names
  ids.a<-unlist(v.all$name[v.all$name %in% symbol.a]) #the gene names of pathway a included in the total graph
  ids.b<-unlist(v.all$name[v.all$name %in% symbol.b]) #the gene names of pathway b included in the total graph
  if(min(c(length(ids.a),length(ids.b))) <= minGeneNum){
    print("Not enough genes presenting in the map")
    return(NA)
  }else{
    ig.ab<-mergePw(ig.all,list(symbol.a,symbol.b))
    tmp.result<-core_rImp(ig.ab,degree.a,degree.b,rp.a.d,rp.b.d,rp.a.u,rp.b.u,penal.isolated, cutImpactSub,exportData)
    tmp.score<-tmp.result[[1]]
    tmp.a<-tmp.result[[2]]
    tmp.b<-tmp.result[[3]]
    if(length(intersect(symbol.a,symbol.b))>0){#if there are intersections between the two pathways
      result.split<-split_intersected_pathways(ig.all,symbol.a,symbol.b)
      ig.split<-result.split[[1]]
      genes1.spl<-result.split[[2]]
      genes2.spl<-result.split[[3]]
      degree.a.spl<-degree.a
      names(degree.a.spl)<-genes1.spl #
      rp.a.d.spl<-rp.a.d
      names(rp.a.d.spl)<-genes1.spl
      rp.a.u.spl<-rp.a.u
      names(rp.a.u.spl)<-genes1.spl
    }else{
      ig.split<-ig.ab
      genes1.spl<-symbol.a
      genes2.spl<-symbol.b
      degree.a.spl<-degree.a
      names(degree.a.spl)<-genes1.spl
      rp.a.d.spl<-rp.a.d
      names(rp.a.d.spl)<-genes1.spl
      rp.a.u.spl<-rp.a.u
      names(rp.a.u.spl)<-genes1.spl
    }
    a.permu<-c()
    b.permu<-c()
    for(i in 1:n.permu){
      g.new.spl<-perturb_pwPos(ig.split,genes1.spl,genes2.spl)
      splitResult<-core_rImp(g.new.spl,degree.a.spl,degree.b,rp.a.d.spl,rp.b.d,rp.a.u.spl,rp.b.u,penal.isolated,cutImpactSub)####
      a.permu<-c(a.permu,splitResult[[2]])
      b.permu<-c(b.permu,splitResult[[3]])
    }
    if(plot){
      ig.split<-set_splitPwColor(ig.split,genes1.spl,genes2.spl)
      l<-layout_with_fr(ig.split)
      plot(ig.split,
           edge.arrow.size=.8,vertex.label.dist=1,vertex.size=5,layout=l)
    }
    tmp.a.bk<-tmp.a
    tmp.b.bk<-tmp.b
    tmp.a<-adjustImpactScore(tmp.a,mean(a.permu))
    tmp.b<-adjustImpactScore(tmp.b,mean(b.permu))
    score.normalized<-tmp.a-tmp.b
    if(is.na(tmp.score)){
      tmp.score<-0
    }
    p.a<-pValueCal(tmp.a.bk,tmp.b.bk,a.permu,b.permu)
    if(exportData){
      result<-list(p.a,score.normalized,tmp.a.bk,tmp.b.bk,mean(a.permu),mean(b.permu),tmp.result[[4]],tmp.result[[5]])
      names(result)<-c("p-value","Normalized score","ScoreA before normalization",
                       "ScoreB before normalization","Mean permutation A","Mean permutation B",
                       "GeneA's impact on gene set B","GeneB's impact on gene set A")
      
    }else{
      result<-list(p.a,score.normalized,tmp.a.bk,tmp.b.bk,mean(a.permu),mean(b.permu))
      names(result)<-c("p-value","Normalized score","ScoreA before normalization",
                       "ScoreB before normalization","Mean permutation A","Mean permutation B")
    }
    return(result)
  }
}

#' @title Adjust the impact score
#' @param impact.score The raw impact score
#' @param meanImpScore the mean of the permutated impact score
#' @description This function first check whether both of the impact score and permutated scores are 0, then adjust the normalized GESTIA score
#' @return GESTIA score after adjustment
#' @export
adjustImpactScore<-function(impact.score,meanImpScore){
  if(impact.score == 0 & meanImpScore == 0){
    score.adj<-0
  }else{
    tmp<-impact.score - meanImpScore
    score.adj<-log(exp(tmp)+1,exp(1))
  }
  return(score.adj)
}

#' @title Calculate p value
#' @param score.a The raw impact score of pathway A on B
#' @param score.b The raw impact score of pathway B on A
#' @param a.permu The permutated impact scores of pathway A on B
#' @param b.permu The permutated impact scores of pathway B on A
#' @description Calculate the p-value based on the raw and permutated impact scores.
#' @return p-value
#' @export
pValueCal<-function(score.a,score.b,a.permu,b.permu){
  score.raw<-score.a - score.b
  permu.sub<-a.permu - b.permu
  if(score.raw > 0){
    tmp.p<-length(permu.sub[permu.sub >= score.raw])/length(permu.sub)
  }else{
    tmp.p<-length(permu.sub[permu.sub <= score.raw])/length(permu.sub)
  }
  if(score.a == 0 & score.b == 0){ #no interactions between A and B
    tmp.p<-1
  }
  return(tmp.p)
}

#' @title Calculate each gene's impact(from pathway A) over pathwayB: gene impact over pathway->giop
#' @param ig.tmp The combined network of the two pathways
#' @param ids.1.tmp The gene symbols of pathway A
#' @param ids.2.tmp The gene symbols of pathway B
#' @param rp.2.tmp The relative position of genes in pathway B, i.e. the proportion of downstream genes for each gene in pathway B
#' @description Calculate the impact of each gene in pathway A on pathway B
#' @return A vector containing the impact of each gene in pathway A
#' @export
giop<-function(ig.tmp,ids.1.tmp,ids.2.tmp,rp.2.tmp){
  impact.1.tmp<-c() #initialize result vector
  for(i in ids.1.tmp){
    dist.tmp <- distances(ig.tmp, v=V(ig.tmp)[name==i], to=V(ig.tmp), weights=NA,mode="out") #get the distance of each gene in the graph to gene i, limit to downstream distance
    dist.names<-colnames(dist.tmp)[dist.tmp==1] #only the direct downstream nodes are considered, because relative position(rp) has already considered downstream nodes
    dist.names<-dist.names[dist.names %in% ids.2.tmp] #limit the downstream nodes to be the ones in pathway 2
    value.single<-sum(rp.2.tmp[names(rp.2.tmp) %in% dist.names]) #sum of the downstream nodes' relative position; the nodes are limited to pathway 2
    impact.1.tmp<-c(impact.1.tmp,value.single) #put the score into the result vector
  }
  names(impact.1.tmp)<-ids.1.tmp #annotate gene names
  return(impact.1.tmp)
}

#' @title GESTIA calculation
#' @param g.tmp Global network constructed from KEGG
#' @param d1.symbol Gene names of pathway A
#' @param d2.symbol Gene names of pathway B
#' @param penal.isolated Penalization factor for isolated genes. Default 0.1
#' @param cutImpact Cutoff for impact score. Scores less than this cutoff will be set to 0. Default 0.
#' @param plot Whether to plot the splitted network. Default FALSE
#' @param n.permu Number of permutations for calculating the empirical null distribution
#' @param minGeneNum Cutoff for number of genes in pathways after id convertion
#' @param exportData Whether to export detail of calculation result. Default FALSE.
#' @description The main function of GESTIA.
#' @return List of results including: p-value, Normalized score,ScoreA before normalization,ScoreB before normalization,Mean permutation A,Mean permutation B
#' @export
gestia<-function(g.tmp,d1.symbol,d2.symbol,penal.isolated=0.1,cutImpact=0,
                        plot=FALSE,n.permu=30,minGeneNum=5,exportData=F){
  v.tmp<-V(g.tmp)
  index1.tmp<-which(!unlist(v.tmp$name) %in% d1.symbol)
  index2.tmp<-which(!unlist(v.tmp$name) %in% d2.symbol)
  g1.tmp<-delete_vertices(g.tmp,index1.tmp)
  g2.tmp<-delete_vertices(g.tmp,index2.tmp)
  ids1.tmp<-unlist(as.vector(V(g1.tmp)$name))
  ids2.tmp<-unlist(as.vector(V(g2.tmp)$name))
  tmp.num<-min(c(length(ids1.tmp),length(ids2.tmp)))
  print(paste("Min gene number in one of the two pathway is ",tmp.num,sep=""))
  if(min(c(length(ids1.tmp),length(ids2.tmp))) <= minGeneNum){
    tmp.num<-min(c(length(ids1.tmp),length(ids2.tmp)))
    return(NA)
  }else{
    rp1.tmp<-relativePosition(ids1.tmp,g1.tmp)
    rp2.tmp<-relativePosition(ids2.tmp,g2.tmp)
    rp1.d<-as.vector(rp1.tmp$tmp.d)
    rp1.u<-as.vector(rp1.tmp$tmp.u)
    names(rp1.d)<-rownames(rp1.tmp)
    names(rp1.u)<-rownames(rp1.tmp)
    rp2.d<-as.vector(rp2.tmp$tmp.d)
    rp2.u<-as.vector(rp2.tmp$tmp.u)
    names(rp2.d)<-rownames(rp2.tmp)
    names(rp2.u)<-rownames(rp2.tmp)
    d1.degree<-degree(g1.tmp)
    d2.degree<-degree(g2.tmp)
    result<-relativeImpact2(g.tmp,d1.degree,d2.degree,rp1.d,rp2.d,rp1.u,rp2.u,
                            penal.isolated,cutImpact,plot,n.permu,minGeneNum,exportData)
    return(result)
  }
}

#' @title Merge two pathways' network
#' @param ig.all The global network
#' @param list.pw A list of the pathways to merge. For example, list.pw[[1]] contains gene symbols of pathway A's genes.
#' @description Merge two pathways' network by deleting unrelated genes. This strategy can maintain the inter-pathway interactions.
#' @return An igraph object of the merged network
#' @export
mergePw<-function(ig.all,list.pw){ #list.pw: lists of pathways to be merged. each single list contains the gene symbols of the pathway
  if(length(list.pw) == 1){
    print("Only one list provided. No need to merge, creating igraph object according to the single list.")
    ig.result<-delete_vertices(ig.all,which(!unlist(V(ig.all)$name) %in% list.pw[[1]]))
  }else if(length(list.pw) >1){
    name.all<-unlist(V(ig.all)$name)
    name.final<-name.all
    for(i in 1:length(list.pw)){
      genes.tmp<-name.all[!name.all %in% list.pw[[i]]]
      name.final<-intersect(name.final,genes.tmp)
    }
    ig.result<-delete_vertices(ig.all,name.final)
  }
  return(ig.result)
}

#' @title Visualize two pathways' components
#' @param ig.all The global network
#' @param genesa The gene symbols of pathway A
#' @param genesb The gene symbols of pathway B
#' @param colorset The color set for the shared genes, genes of pathway A, and the genes of pathway B. Default yellowgreen, tomato, gold
#' @param isolated Whether to plot the isolated genes. Default FALSE
#' @param ly Layout of the igraph. Available options: kk, fr, random, grid
#' @description Visualization of two pathways.
#' @return Plot directly.
#' @export
visual2pw<-function(ig.all,genesa,genesb,colorset=c("yellowgreen","tomato","gold"),
                    isolated=FALSE,ly=c("kk","fr","random","grid")){
  ig.ab<-mergePw(ig.all,list(genesa,genesb))
  v.ab<-V(ig.ab)
  type.commuA<-ifelse(unlist(v.ab$name) %in% genesa,1,3)
  type.commuB<-ifelse(unlist(v.ab$name) %in% genesb,1,2)
  type.all<-pmax(type.commuA,type.commuB)
  V(ig.ab)$community<-type.all
  colrs <- adjustcolor( colorset, alpha=.6)
  if(length(ly)==4){
    ly="kk"
  }
  if(isolated){
    if(ly == "kk"){
      l <- layout_with_kk(ig.ab)
    }else if(ly == "fr"){
      l <- layout_with_fr(ig.ab)
    }else if(ly == "random"){
      l <- layout_randomly(ig.ab)
    }else if(ly == "grid"){
      l <- layout_on_grid(ig.ab)
    }
    plot(ig.ab, vertex.color=colrs[V(ig.ab)$community],
         edge.arrow.size=.8,vertex.label.dist=1,vertex.size=5,layout=l)
  }else{
    name.isolated<-names(degree(ig.ab))[degree(ig.ab) == 0]
    index.tmp<-which(unlist(V(ig.ab)$name) %in% name.isolated)
    ig.tmp<-delete_vertices(ig.ab,index.tmp)
    if(ly == "kk"){
      l <- layout_with_kk(ig.tmp)
    }else if(ly == "fr"){
      l <- layout_with_fr(ig.tmp)
    }else if(ly == "random"){
      l <- layout_randomly(ig.tmp)
    }else if(ly == "grid"){
      l <- layout_on_grid(ig.tmp)
    }
    plot(ig.tmp, vertex.color=colrs[V(ig.tmp)$community],
         edge.arrow.size=.8,vertex.label.dist=1,vertex.size=5,layout=l)
  }
}

#' @title Calculate matrix of pairwise pathway impact of row over column
#' @param list.gs The list of genes of each pathway/gene set.
#' @param ig.net The global network
#' @param imCut The cutoff of impact scores. Default 0
#' @param penal.isolated The penalization factor for isolated genes. Default 0.1
#' @param plot Whether to plot the splitted network. Default FALSE
#' @param n.permu The number of permutations for calculation of empirical null distribution
#' @description A wrapper to calculate the GESTIA matrix of a list of pathways. The core input is the list of pathways
#' @return A list of two matrices, the normalized GESTIA scores and the p-values
#' @export
calculateMa<-function(list.gs,ig.net,imCut = 0,
                      penal.isolated=0.1,
                      plot=FALSE,n.permu=30){
  genes.net<-V(ig.net)$name
  name.gs<-names(list.gs)
  num.gs<-length(list.gs)
  num.exist<-c()
  for(i in 1:num.gs){
    num.exist<-c(num.exist,length(intersect(genes.net,list.gs[[i]])))
  }
  index.fail<-which(num.exist <= 5)
  if(length(index.fail)>0){
    print(paste0("Excluded gene sets without sufficient number of existing genes in the network: ",
                 paste(index.fail,collapse = ",")))
    list.gs<-list.gs[!1:num.gs %in% index.fail]
    name.gs<-name.gs[!1:num.gs %in% index.fail]
    num.gs<-length(list.gs)
  }
  ma.norm<-matrix(rep(0,num.gs*num.gs),ncol=num.gs)
  ma.ori<-matrix(rep(0,num.gs*num.gs),ncol=num.gs)
  for(i in 1:num.gs){
    for(j in 1:num.gs){
      if(i == j){
        ma.norm[i,j]<-0
        ma.ori[i,j]<-0
      }else{
        result<-gestia(ig.net,list.gs[[i]],list.gs[[j]],penal.isolated = penal.isolated,
                              cutImpact = imCut, plot = plot, n.permu = n.permu)
        ma.norm[i,j]<-result[[2]]
        ma.ori[i,j]<-result[[1]]
      }
      print(paste0("Done for pw pair [",i,",",j,"] in ",num.gs," rows, ",
                   name.gs[i]," Vs ",name.gs[j],
                   "score: ",round(ma.norm[i,j],4)," ",Sys.time()))
    }
  }
  rownames(ma.norm)<-name.gs
  colnames(ma.norm)<-name.gs
  rownames(ma.ori)<-name.gs
  colnames(ma.ori)<-name.gs
  return(list(ma.norm,ma.ori))
}

#' @title Split two intersected pathways
#' @param ig.all The global network
#' @param genes1 The gene symbols of pathway A's genes
#' @param genes2 The gene symbols of pathway B's genes
#' @description Split two intersected pathways by duplicating the shared genes. Genes1 will be modified to contain fake vertices with "fa_"
#' @return A list containing: the igraph object of the splitted network, the new genes1, the vector genes2
#' @export
split_intersected_pathways<-function(ig.all,genes1,genes2){
  ig.ab<-mergePw(ig.all,list(genes1,genes2))
  genes1<-genes1[genes1 %in% unlist(V(ig.ab)$name)] #exclude genes not presenting in the graph
  genes2<-genes2[genes2 %in% unlist(V(ig.ab)$name)]
  vertices2add<-intersect(genes1,genes2) #the order of intersecting genes are based on genes1
  vertices.fake<-paste("fa_",vertices2add,sep="")
  ig.new<-ig.ab
  for(i in 1:length(vertices2add)){
    e.from<-ends(ig.ab,E(ig.ab)[from(match(vertices2add[i],V(ig.ab)$name))])
    e.to<-ends(ig.ab,E(ig.ab)[to(match(vertices2add[i],V(ig.ab)$name))])
    e.from[,1]<-vertices.fake[i]
    e.to[,2]<-vertices.fake[i]
    e2add<-c(as.character(t(e.from)),as.character(t(e.to)))
    ig.new<-add_vertices(ig.new,1,name = vertices.fake[i],
                         shape = "circle",
                         color = "blue",
                         pathway = "fake")
    ig.new<-add_edges(ig.new,e2add)
  }
  genes1[genes1 %in% vertices2add]<-vertices.fake #the order of genes1's intersected genes are the same with vertices.fake. For genes2, it's not always the same
  return(list(ig.new,genes1,genes2))
}

#' @title Permutation of the inter-pathway interactions
#' @param g.ab The network of the two pathways
#' @param gss1 The gene symbols of pathway A's genes
#' @param gss2 The gene symbols of pathway B's genes
#' @description Randomly assign the edges between two pathways, maintains the number of edges in/out the pathways, remains the structure of each pathway
#' @return An igraph object of the permutated network of A and B
#' @export
perturb_pwPos<-function(g.ab,gss1,gss2){
  e.ab<-ends(g.ab,E(g.ab))
  e.a2b<-ends(g.ab,E(g.ab)[from(match(gss1,unlist(V(g.ab)$name)))]) #all the edges coming out of a's nodes
  e.b2a<-ends(g.ab,E(g.ab)[from(match(gss2,unlist(V(g.ab)$name)))]) #all the edges coming out of b's nodes
  e.a2b<-matrix(e.a2b[e.a2b[,2] %in% gss2,],ncol=2)
  e.b2a<-matrix(e.b2a[e.b2a[,2] %in% gss1,],ncol=2)
  n.a2b<-dim(e.a2b)[1]
  n.b2a<-dim(e.b2a)[1]
  pool.a2b<-expand.grid(gss1,gss2)
  pool.b2a<-expand.grid(gss2,gss1)
  eRd.a2b<-as.matrix(pool.a2b[sample(1:dim(pool.a2b)[1],n.a2b),])
  eRd.b2a<-as.matrix(pool.b2a[sample(1:dim(pool.b2a)[1],n.b2a),])
  e.ab.bk<-e.ab
  e.ab[e.ab.bk[,1] %in% e.a2b[,1] & e.ab.bk[,2] %in% e.a2b[,2],] <-eRd.a2b #replace the edges (a->b) by random a->b edges
  e.ab[e.ab.bk[,1] %in% e.b2a[,1] & e.ab.bk[,2] %in% e.b2a[,2],] <-eRd.b2a
  e2add<-c(as.character(t(e.ab)))
  g.ab.new<-delete_edges(g.ab,E(g.ab))
  g.ab.new<-add_edges(g.ab.new,e2add)
  return(g.ab.new)
}

#' @title Calculation of the raw impact score
#' @param ig.ab The network of the two pathways
#' @param degree.a The degree of each genes in pathway A
#' @param degree.b The degree of each genes in pathway B
#' @param rp.a.d The proportion of downstream genes for each gene in pathway A
#' @param rp.b.d The proportion of downstream genes for each gene in pathway B
#' @param rp.a.u The proportion of upstream genes for each gene in pathway A
#' @param rp.b.u The proportion of upstream genes for each gene in pathway B
#' @param penal.isolated The penalization factor for isolated genes. Default 0.1
#' @param cutImpactSub The cutoff for impact score. Default 0
#' @param exportData Whether to output details of the calculation. Default F
#' @description Calculation of the raw impact score. The core of the GESTIA
#' @return A list containing: Raw GESTIA score, the impact score of A on B, and the impact score of B on A
core_rImp<-function(ig.ab,degree.a,degree.b,rp.a.d,rp.b.d,rp.a.u,rp.b.u,
                    penal.isolated=0.1,cutImpactSub=0,exportData=F){
  symbol.a<-names(degree.a)
  symbol.b<-names(degree.b)
  v.ab<-V(ig.ab)#get the vertice names of the total graph, here equal to gene names
  ids.a<-unlist(v.ab$name[v.ab$name %in% symbol.a]) #the gene names of pathway a included in the total graph
  ids.b<-unlist(v.ab$name[v.ab$name %in% symbol.b]) #the gene names of pathway b included in the total graph
  #consider gene's relative position in pathway a when calculating impact factor
  impact.a<-giop(ig.ab,ids.a,ids.b,rp.b.d) #calculate each gene's impact from pathway a over pathway b
  impact.b<-giop(ig.ab,ids.b,ids.a,rp.a.d) #calculate each gene's impact from pathway b over pathway a
  impact.a<-impact.a[order(names(impact.a))]
  impact.b<-impact.b[order(names(impact.b))]
  rp.a.u<-rp.a.u[order(names(rp.a.u))]
  rp.b.u<-rp.b.u[order(names(rp.b.u))]
  impactFinal.a<-impact.a*rp.a.u
  impactFinal.b<-impact.b*rp.b.u
  impactFinal.a[impactFinal.a <= cutImpactSub] = 0
  impactFinal.b[impactFinal.b <= cutImpactSub] = 0
  a<-sum(impactFinal.a) #take the sum
  b<-sum(impactFinal.b)
  tmp.score<-(a-b)/(a+b)*sqrt(a^2+b^2) #calculate final score
  if(is.na(tmp.score)){
    tmp.score<-0
  }
  if(exportData){
    result<-list(tmp.score,a,b,impactFinal.a,impactFinal.b)
  }else{
    result<-list(tmp.score,a,b)
  }
  return(result)
}

#' @title Set color for the splitted pathways for better visualization
#' @param ig.split The splitted network, an igraph object
#' @param gss1 The gene symbols of pathway A's genes
#' @param gss2 The gene symbols of pathway B's genes
#' @param verticesColor Whether to paint the vertices with color. Default TRUE
#' @description Set the color for the vertices and edges in the splitted network.
#' @return An igraph object with colors
#' @export
#set color for the split pathways for better visualization. Vertices and edges a->b are colored by "tomato", the other way around are colored by "gold"
set_splitPwColor<-function(ig.split,gss1,gss2,verticesColor=TRUE){
  if(verticesColor){
    V(ig.split)$color[V(ig.split)$name %in% gss1] <-"tomato"
    V(ig.split)$color[V(ig.split)$name %in% gss2] <-"gold"
  }
  e.a2b<-ends(ig.split,E(ig.split)[from(match(gss1,unlist(V(ig.split)$name)))]) #all the edges coming out of a's nodes
  e.b2a<-ends(ig.split,E(ig.split)[from(match(gss2,unlist(V(ig.split)$name)))]) #all the edges coming out of b's nodes
  e.a2b<-matrix(e.a2b[e.a2b[,2] %in% gss2,],ncol=2)
  e.b2a<-matrix(e.b2a[e.b2a[,2] %in% gss1,],ncol=2)
  e.ab<-ends(ig.split,E(ig.split))
  color.e<-rep("grey",dim(e.ab)[1])
  color.e[e.ab[,1] %in% e.a2b[,1] & e.ab[,2] %in% e.a2b[,2]]<-"tomato"
  color.e[e.ab[,1] %in% e.b2a[,1] & e.ab[,2] %in% e.b2a[,2]]<-"gold"
  E(ig.split)$color<-color.e
  return(ig.split)
}

#' @title Log transformation of the matrix
#' @param ma The matrix to be transformed
#' @param cut The cutoff for each value in the matrix after transformation
#' @description Log transformation of the GESTIA matrix, with sign maintained.
#' @return The matrix after transformation
log2Transform<-function(ma,cut=0.1){
  ma.sign<-sign(ma)
  ma<-log2(abs(ma)+1)*ma.sign
  for(i in 1:dim(ma)[1]){
    ma[i,ma[i,] < cut] = 0
  }
  return(ma)
}

#' @title Convert GESTIA score matrix to igraph object
#' @param ls.ma The list of two matrices, first is the GESTIA scores, second is the p-values
#' @param pCut The cutoff for p-value
#' @param scoreCut The cutoff for GESTIA score
#' @description Filter and then convert GESTIA score matrix to an igraph object, using the matrix as adjacency matrix
#' @return An igraph object of the GESTIA score matrix
#' @export
prepareMaPlot<-function(ls.ma,pCut=0.05,scoreCut=0){
  ma.sc.score<-ls.ma[[1]]
  ma.sc.p<-ls.ma[[2]]
  n.ma<-dim(ma.sc.score)[1]
  for(i in 1:n.ma){
    ma.sc.score[i,ma.sc.p[i,] > pCut]<-0
  }
  for(i in 1:n.ma){
    ma.sc.score[i,ma.sc.score[i,] < scoreCut]<-0
  }
  sf<-sum(ma.sc.score)/length(which(ma.sc.score>0))
  ig.pw<-graph_from_adjacency_matrix(ma.sc.score/sf,weighted = TRUE)
  E(ig.pw)$width<-E(ig.pw)$weight
  
  isolated<-names(degree(ig.pw))[degree(ig.pw) == 0]
  id.tmp<-which(unlist(V(ig.pw)$name) %in% isolated)
  ig.pw.tmp<-delete_vertices(ig.pw,id.tmp)
  return(ig.pw.tmp)
}