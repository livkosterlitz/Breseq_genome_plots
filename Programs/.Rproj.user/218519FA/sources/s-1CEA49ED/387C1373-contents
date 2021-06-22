# Written by Ben Kerr see Lindsey et. al. 2013 https://www.nature.com/articles/nature11879
# Updated by Olivia Kosterlitz 

#This program will plot the mutations in a set of strains from a set of treatments
#The arguments are described below

library("optparse")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--plotting_options"), type="character", default=NULL, 
              help="plotting options file name", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Command line arguments:
#-------------------
## raw.seq.df is the data frame with all the mutations (must have the following columns: 
# treatment, strain, base, mut.type and delta.bases (note: each mutation would be listed as its own row)
raw.seq.df <- read.csv(opt$file, header=TRUE)

#-------------------
## plotting.options is the data frame with all the plotting options (each row  
# passes a customization)
## Rows
# xmin is the minimum x value of the plot (good default = -0.15)
# ymin is the minimum y value of the plot (good default = -0.1)
# ymax is the maximum y value of the plot (good default = 1.05)
# Note: xmax is  is automatically calculated in the script
# N.bases is number of bases in the gene or the genome. This is the length of the reference sequence used to call variants for each strain.
# section.colors is a vector giving the colors of the highlighted sections. Pass atleast two colors for alternate highlighting. (Good default = snow2 & snow3)
# bp.fill is an option to add a bar to the figure where bp in the magnified sections are highlighted. To turn on pass TRUE. To turn off pass FALSE. 
# bp.colors is a vector giving the colors for the bp bar. Pass atleast two colors for alternate highlighting. (Good default = snow2 & snow3)
# treatments is a vector of all the treatments as they appear in the data frame. Automatically use all treatments in dataframe if FALSE. To plot a subset pass the treatment names. 
# treatment.colors is a vector of the colors for each treatment. If FALSE colors will be automatically set. 
# treatment.labels is a vector with the labels to place on the mutation table. If FALSE labels will be automatically set from the input data set. To customize pass desired treatment names.
# gene.color is the color of the gene in the figure (good default="ivory"). This is the top bar in the figure. 
# gene.width is the width of the gene in the figure (good default=0.03). This is the width of the top bar. 
# gene.section.border is the color of the border to outline the gene in the top gene bar. To turn off pass FALSE. 
# mag.border is the color of the border to outline the magnifying polygons. To turn off pass FALSE. 
# min.gap.length is minimum number of bases consecutive mutations need to spaced in order to be in different magnified sections (good default=30bp). Can also pass 'gene' to have variants within the same gene be plotted within the same 'magnified' section. 
# section.buffer is number of bases subtracted and added to the ends of each section to buffer the sections in the figure (good default = 1bp)
# mag.spacer is width of the spacing between the gene and the magnified section bar of the table of mutations (default = 0.05)
# tab.spacer is width of the spacing between the magnified section bar and the table of mutations (good default = 0.025). this is also the section where the bp bar is placed if turned on (see bp.fill)
# treatment.label.cex is the size of the text of the treatment labels (good default = 1.0)
# treatment.label.spacer is the spacing between the table and the treatment labels (good default = 0)
# mutation.colors are the colors of the different types of mutations in the table
# hist.spacer is the spacing between the table and the histogram (default = 0.02)
# hist.unit.height is the height of the histogram
# hist.axis.spacer is the spacing between the histogram and its axis (good default = 0.02) 
# hist.ticks.length is the length of the ticks on the histogram axis (default = 0.02)
# hist.numbers.cex is the size of the text of the numbers on the histogram axis (default = 1) 
# hist.numbers.spacer is the spacing between the ticks and the numbers (default = 0.05)
# hist.label is the label on the histogram axis
# hist.label.spacer is the spacing between the numbers and label on the histogram (default=0.1) 
# hist.label.cex is the size of the text of the label on the histogram axis (default=1)
# row.spec.width is the width for strains if you want to override the automatic calculation. To keep automatic calculation pass 0.  
# border.color is the color outlining the magnified sections within the strain rows. (good default = lightyellow4)
# treatment.border.color is the color outlining the treatments within the strain rows. (good default = black)
# output_filename is to customize the PDF output filename
#-------------------

plotting.options <- read.csv(opt$plotting_options, header = TRUE, na.strings=c(""))


# parse the plotting options file to customize the plotting
N.bases <- as.numeric(plotting.options$Value[plotting.options$Options == "N.bases"])
section.colors <- unname(unlist(plotting.options[plotting.options$Options == "section.colors", 3:length(plotting.options)]))
section.colors <- section.colors[complete.cases(section.colors)]
bp.fill <- plotting.options$Value[plotting.options$Options == "bp.fill"]
bp.colors <- unname(unlist(plotting.options[plotting.options$Options == "bp.colors", 3:length(plotting.options)]))
bp.colors <- bp.colors[complete.cases(bp.colors)]
treatments <- unname(unlist(plotting.options[plotting.options$Options == "treatments", 3:length(plotting.options)]))
treatments <- treatments[complete.cases(treatments)]
if (treatments[1] == FALSE) {
  treatments <- unique(raw.seq.df$treatment) # default is all of the treatments. User can pass a subset of the samples if desired
}
treatment.colors <- unname(unlist(plotting.options[plotting.options$Options == "treatment.colors", 3:length(plotting.options)]))
treatment.colors <- treatment.colors[complete.cases(treatment.colors)]
if (treatment.colors[1] == FALSE) {
  library('RColorBrewer')
  treatment.colors <- brewer.pal(length(treatments),"Set3")
}
treatment.labels <- unname(unlist(plotting.options[plotting.options$Options == "treatment.labels", 3:length(plotting.options)]))
treatment.labels <- treatment.labels[complete.cases(treatment.labels)]
if (treatment.labels[1] == FALSE) {
  treatment.labels <- treatments
}
gene.color = plotting.options$Value[plotting.options$Options == "gene.color"]
gene.width = as.numeric(plotting.options$Value[plotting.options$Options == "gene.width"])
gene.section.border = plotting.options$Value[plotting.options$Options == "gene.section.border"]
mag.border = plotting.options$Value[plotting.options$Options == "mag.border"]
min.gap.length = plotting.options$Value[plotting.options$Options == "min.gap.length"]
if (min.gap.length != 'gene') {
  min.gap.length = as.numeric(min.gap.length)
}
section.buffer = as.numeric(plotting.options$Value[plotting.options$Options == "section.buffer"])
mag.spacer= as.numeric(plotting.options$Value[plotting.options$Options == "mag.spacer"])
tab.spacer= as.numeric(plotting.options$Value[plotting.options$Options == "tab.spacer"])
treatment.label.cex= as.numeric(plotting.options$Value[plotting.options$Options == "treatment.label.cex"])
treatment.label.spacer= as.numeric(plotting.options$Value[plotting.options$Options == "treatment.label.spacer"])
snp.mutation.color = plotting.options$Value[plotting.options$Options == "snp.mutation.color"]
ins.mutation.color = plotting.options$Value[plotting.options$Options == "ins.mutation.color"]
del.mutation.color = plotting.options$Value[plotting.options$Options == "del.mutation.color"]
hist.spacer = as.numeric(plotting.options$Value[plotting.options$Options == "hist.spacer"])
hist.unit.height= as.numeric(plotting.options$Value[plotting.options$Options == "hist.unit.height"])
hist.axis.spacer= as.numeric(plotting.options$Value[plotting.options$Options == "hist.axis.spacer"])
hist.ticks.length= as.numeric(plotting.options$Value[plotting.options$Options == "hist.ticks.length"])
hist.numbers.cex= as.numeric(plotting.options$Value[plotting.options$Options == "hist.numbers.cex"])
hist.numbers.spacer= as.numeric(plotting.options$Value[plotting.options$Options == "hist.numbers.spacer"])
hist.label = plotting.options$Value[plotting.options$Options == "hist.label"]
hist.label.spacer= as.numeric(plotting.options$Value[plotting.options$Options == "hist.label.spacer"])
hist.label.cex= as.numeric(plotting.options$Value[plotting.options$Options == "hist.label.cex"])
xmin= as.numeric(plotting.options$Value[plotting.options$Options == "xmin"])
xmax = 1+hist.unit.height+hist.spacer+0.05
ymin= as.numeric(plotting.options$Value[plotting.options$Options == "ymin"])
ymax= as.numeric(plotting.options$Value[plotting.options$Options == "ymax"])
row.spec.width= as.numeric(plotting.options$Value[plotting.options$Options == "row.spec.width"])
border.color = plotting.options$Value[plotting.options$Options == "border.color"]
treatment.border.color = plotting.options$Value[plotting.options$Options == "treatment.border.color"]
output_filename = plotting.options$Value[plotting.options$Options == "output_filename"]

  
#start file
pdf(file = paste(output_filename, ".pdf", sep = ""),
    width = 15,
    height = 8)

#Reserve space for the plot
plot(c(xmin,xmax),c(ymin,ymax),col="white",axes=FALSE,ylab="",xlab="")

#Draw full gene
rect(0, 1-gene.width, 1, 1, col=gene.color) #xleft, ybottom, xright, ytop

#Pick the specified treatments out of the data frame
focal.raw.seq.df<-raw.seq.df[is.element(raw.seq.df$treatment,treatments),]	

#Construct site vectors/parameters
focal.raw.seq.df.order <- focal.raw.seq.df[order(focal.raw.seq.df$base),] # orders the mutations by position
sites <- sort(unique(focal.raw.seq.df$base)) #gets all the positions for sites that will be plotted that are unique
genes <- vector(mode="character", length=length(sites)) # a vector for the gene names for each mutation
for (i in c(1:length(sites))) {
  genes[i] <- toString(focal.raw.seq.df$gene[focal.raw.seq.df$base==sites[i]][1])
}
Nsit <- length(sites) #number of variant sites that will be plotted

#Determine the number and boundaries of genetic sections to plot
sec.init <- (sites[1]-section.buffer) #grabs the first site with a mutation and subtracts the amount of bp before the variant to start plotting
if(Nsit == 1) { #if the number of variant site is one
  sec.fin <- (sites[1]+section.buffer) 
  N.sec <- 1
}
if(Nsit > 1) { #if there is more than one variant site to plot
  N.sec <- 1 # counts the number of sections to be drawn
  j <- 1
  for(i in 2:Nsit) { #loop through the number of variants starting at 2. For example, i will be 2 then 3 then 4
    if (min.gap.length != "gene") {
      if((sites[i]-sites[i-1]) > min.gap.length) { #if the current mutation is more than the space specified between variants to be plotted in the same panel. Do this loop!
        sec.init <- c(sec.init, sites[i]-section.buffer) #add to the section initiation vector the starting position for the panel
        ifelse(j==1, sec.fin <- sites[i-1] + section.buffer, sec.fin <- c(sec.fin,sites[i-1]+section.buffer)) # if on the first section initiate the sec.fin vector else append to the vector
        j <- j+1
        N.sec <- N.sec+1
      }
    } else {
      if (genes[i] != genes[i-1]) { #if the current mutation is not in the same gene as the previous variant, Do this loop!
        sec.init<-c(sec.init,sites[i]-section.buffer) #add to the section initiation vector the starting position for the panel
        ifelse(j==1,sec.fin<-sites[i-1]+section.buffer,sec.fin<-c(sec.fin,sites[i-1]+section.buffer))
        j<-j+1
        N.sec<-N.sec+1
      }
    }
  }
  sec.fin<-c(sec.fin,sites[Nsit]+section.buffer)
}	

section.colors <- rep(section.colors, ceiling(N.sec/length(section.colors)))
#Draw sections on gene
for(i in 1:N.sec) {
  rect(sec.init[i]/N.bases, 1-gene.width, sec.fin[i]/N.bases, 1, col=section.colors[i], border=gene.section.border) #xleft, ybottom, xright, ytop
}	
lines(c(0,1),c(1,1))

#Compute the lengths of each section and total length of all sections
new.run <- numeric(N.sec) # length (bp) of each section
tot <- 0 # total length (bp) of all sections
for(i in 1:N.sec) {
  new.run[i] <- (sec.fin[i]-sec.init[i])+1 #this needed a plus one
  tot <- tot + new.run[i]
}

#Put together "new" sections that are confluent on the interval [0,1], preserving relative lengths
new.sec.init <- numeric(N.sec)
new.sec.fin <- numeric(N.sec)
tot2 <- 0
for(i in 1:N.sec) {
  new.sec.init[i] <- tot2
  new.sec.fin[i] <- tot2+(new.run[i]/tot)
  tot2 <- tot2 + (new.run[i]/tot)
}

#"Magnify" the sections
for(i in 1:N.sec) {
  polygon(c(sec.init[i]/N.bases, new.sec.init[i], new.sec.fin[i], sec.fin[i]/N.bases),
          c(1-gene.width, 1-gene.width-mag.spacer, 1-gene.width-mag.spacer, 1-gene.width),
          col = section.colors[i], border = mag.border) #the x vector c(top left, bottom left, bottom right, top right), the y vector (top left, bottom left, bottom right, top right)
}
lines(c(0,1),c(1-gene.width,1-gene.width))


#Function to organize strains in a list first by number of mutations
#then by location of the first mutation (this makes the mutation table
#a bit prettier)
OrderStrains<-function(X.df) {
  Xstrains<-sort(unique(X.df$strain))
  len<-length(Xstrains)
  x<-rep(0,len)
  y<-rep(0,len)
  z<-rep(0,len)
  SO.df<-data.frame(strain=x,base=y,N.mut=z)
  for(i in 1:len) {
    SO.df$strain[i]<-X.df[X.df$strain==Xstrains[i],]$strain[1]
    SO.df$base[i]<-min(X.df[X.df$strain==Xstrains[i],]$base)
    nm<-length(X.df[X.df$strain==Xstrains[i],]$base)
    count=0
    for(j in X.df[X.df$strain==Xstrains[i],]$mut.type) {
      if(j == 'del') {
        count <- count + 1
      }
    }
    SO.df$N.mut[i] <- (nm - (count/2))
  }
  SO.df<-SO.df[order(SO.df$N.mut,SO.df$base),]
  SO.df
}

#Organize strains in all treatments and collect information on the number
#of mutations in each.	
for(i in 1:length(treatments)) { # loop through the treatments
  T.df<-OrderStrains(raw.seq.df[raw.seq.df$treatment==treatments[i],]) # pass all of the mutations for a particular treatment to the OrderStrains function
  ifelse(i==1, strains<-T.df$strain, strains<-c(strains,T.df$strain))
  ifelse(i==1, mutHist<-T.df$N.mut, mutHist<-c(mutHist,T.df$N.mut))
}

Nstr<-length(strains) #the number of strains to be plotted

# calculate the row width to fill the y-axis
if (row.spec.width == 0) {
  row.width<-(1-gene.width-mag.spacer-tab.spacer)/(Nstr+1)  #the width of each row (strain)
} else {
  row.width <- row.spec.width
}

#Draw confluent sections at row height for the magnified bar
for(i in 1:N.sec) {
  rect(new.sec.init[i], 1-gene.width-mag.spacer-row.width, new.sec.fin[i], 1-gene.width-mag.spacer, col=section.colors[i]) #xleft, ybottom, xright, ytop
}


#Put together a list of all the bases in the table
focal.bases <- c((sec.init[1]):(sec.fin[1]))
if(N.sec>1) {
  for(i in 2:N.sec) {
    nex <- c((sec.init[i]):(sec.fin[i]))
    focal.bases <- c(focal.bases, nex)
  }
}

#Draw the proper background colors for each strain
current.treatment <- raw.seq.df[raw.seq.df$strain==strains[1],]$treatment[1]
for(j in 1:length(treatments)) {
  if(current.treatment==treatments[j]) {
    current.color <- treatment.colors[j]
  }
}
row.of.last.color <- 0
for(i in 1:Nstr) {
  if(raw.seq.df[raw.seq.df$strain==strains[i],]$treatment[1] != current.treatment) {	
    rect(0, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i-1)), 1, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(row.of.last.color)), col=current.color, border=FALSE) #xleft, ybottom, xright, ytop
    for(j in 1:length(treatments)) {
      if(raw.seq.df[raw.seq.df$strain==strains[i],]$treatment[1]==treatments[j]) {
        current.color<-treatment.colors[j]
      }
    }
    row.of.last.color<-(i-1)
    current.treatment<-raw.seq.df[raw.seq.df$strain==strains[i],]$treatment[1]
  }
  if(i==Nstr) {
    rect(0, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i)), 1, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(row.of.last.color)), col=current.color, border=FALSE)
  }	
}

#draw in the tab.spacer bp rectangles
if (bp.fill != FALSE){
  lfb<-length(focal.bases)
  bp.colors <- rep(bp.colors, ceiling(lfb/length(bp.colors)))
  bp_size <- 1/lfb
  for (pos in 1:lfb) {
    rect((pos-1)*bp_size, 1-gene.width-mag.spacer-row.width-tab.spacer, pos*bp_size, 1-gene.width-mag.spacer-row.width, col= bp.colors[pos], border=NA)
  }
  rect(0, 1-gene.width-mag.spacer-row.width-tab.spacer, 1, 1-gene.width-mag.spacer-row.width, col= NA, border='black')
}

#Draw the mutations for each strain
for(i in 1:Nstr) {
  mutations<-raw.seq.df[raw.seq.df$strain==strains[i],]$base
  type<-raw.seq.df[raw.seq.df$strain==strains[i],]$mut.type
  lfb<-length(focal.bases)
  for(j in 1:length(mutations)) {
    pos<-which(focal.bases==mutations[j])
    if(type[j]=='snp') {
      mut.col <- snp.mutation.color
    }  
    if(type[j]=='ins') {
      mut.col <- ins.mutation.color
    }
    if(type[j]=='del') {
      mut.col <- del.mutation.color
    }
    rect((pos-1)/lfb, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i)), pos/lfb,1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i-1)), col=mut.col, border=NA)
    if(type[j]=='del') {
      db <- raw.seq.df[raw.seq.df$strain==strains[i] & raw.seq.df$base==mutations[j],]$delta.bases
      if(db>0) {
        end.pos<-which(focal.bases==(mutations[j]+db))
        rect((pos-1)/lfb, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i)), end.pos/lfb, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i-1)), col=mut.col, border=NA)
        #        for(p in pos:end.pos) {
        #          rect((p-1)/lfb, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i)), p/lfb, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i-1)), col=mut.col, border=NA)
        #        }
        if(length(raw.seq.df[raw.seq.df$strain==strains[i] & 
                             raw.seq.df$base==(mutations[j]+db),]$delta.bases) != 1) {
          print(paste("There is an error; strain ", strains[i], 
                      " should have a deletion ending at ", mutations[j]+db, sep=""))
        }
      }
    }
  }
}

#Draw boundaries around all magnified sections in the strain rows
for(i in 1:N.sec) {
  rect(new.sec.init[i], 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*Nstr), new.sec.fin[i], 1-gene.width-mag.spacer-row.width-tab.spacer, border = border.color)
}

#Draw lines between the strains in the mutation table
row.of.last.color<-0
for(i in 1:Nstr) {
  if(i!=Nstr) {	
    rect(0, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i-1)),1, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i)), border = border.color)
    for(j in 1:length(treatments)) {
      if(raw.seq.df[raw.seq.df$strain==strains[i],]$treatment[1]==treatments[j]) {
        current.color<-treatment.colors[j]
      }
    }
    row.of.last.color<-(i-1)
    current.treatment<-raw.seq.df[raw.seq.df$strain==strains[i],]$treatment[1]
  }
  if(i==Nstr) {
    rect(0,1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i)),1,
         1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(row.of.last.color)),border=border.color)
  }	
}

#Draw lines around the treatments 
current.treatment <- raw.seq.df[raw.seq.df$strain==strains[1],]$treatment[1]
for(j in 1:length(treatments)) {
  if(current.treatment==treatments[j]) {
    current.color<-treatment.colors[j]
  }
}
row.of.last.color<-0

for(i in 1:Nstr) {
  if(raw.seq.df[raw.seq.df$strain==strains[i],]$treatment[1] != current.treatment) {	
    rect(0, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i-1)), 1, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(row.of.last.color)), border=treatment.border.color)
    for(j in 1:length(treatments)) {
      if(raw.seq.df[raw.seq.df$strain==strains[i],]$treatment[1]==treatments[j]) {
        current.color<-treatment.colors[j]
      }
    }
    row.of.last.color<-(i-1)
    current.treatment<-raw.seq.df[raw.seq.df$strain==strains[i],]$treatment[1]
  }
  if(i==Nstr) {
    rect(0,1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i)),1,
         1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(row.of.last.color)), border=treatment.border.color)
  }	
}


#Draw histogram to the right of the table
hist.unit.height <- hist.unit.height/max(mutHist)
for(i in 1:Nstr) {
  for(j in 1:length(treatments)) {
    if(raw.seq.df[raw.seq.df$strain==strains[i],]$treatment[1]==treatments[j]) {
      color<-treatment.colors[j]
    }
  }	
  rect(1+hist.spacer, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i)), 1+hist.spacer+mutHist[i]*hist.unit.height, 1-gene.width-mag.spacer-row.width-tab.spacer-(row.width*(i-1)),
       col=color)
}


#Plot the treatment labels to the left of the table
Tlen <- numeric(length(treatments))
Totlen <- 0
for(i in 1:length(treatments)) {
  Tlen[i] <- length(unique(raw.seq.df[raw.seq.df$treatment==treatments[i],]$strain))
  Totlen <- (Totlen+Tlen[i])
}
Ty<-numeric(length(treatments))
RunTot<-0
if (row.spec.width == 0) {
  for(i in 1:length(treatments)) {
    Ty[i]<-(1-gene.width-mag.spacer-row.width-tab.spacer)*(1-RunTot-0.5*(Tlen[i]/Totlen))
    text(-treatment.label.spacer, Ty[i], treatment.labels[i], col="black",cex=treatment.label.cex, pos = 2)
    RunTot<-(RunTot+(Tlen[i]/Totlen))
  }
} else {
  for(i in 1:length(treatments)) {
    mutant.space <- row.spec.width * Nstr
    Ty[i]<-(1-gene.width-mag.spacer-row.spec.width-tab.spacer)-(mutant.space-(mutant.space*(1-RunTot-0.5*(Tlen[i]/Totlen))))
    text(-treatment.label.spacer, Ty[i], treatment.labels[i],  
         col="black",cex=treatment.label.cex, pos = 2)
    RunTot<-(RunTot+(Tlen[i]/Totlen))
  }
}

#Plot the histogram axis, numbers, and label
if (row.spec.width == 0) {
  arrows(1+hist.spacer,-hist.axis.spacer,1+hist.spacer+(max(mutHist)+1)*hist.unit.height,-hist.axis.spacer,length=2*hist.unit.height)
  for(i in 1:max(mutHist)) {
    lines(c(1+hist.spacer+(i)*hist.unit.height,1+hist.spacer+(i)*hist.unit.height),
          c(-hist.axis.spacer,-hist.axis.spacer-hist.ticks.length))
    text(1+hist.spacer+(i)*hist.unit.height,-hist.numbers.spacer,i,cex=hist.numbers.cex,pos=1)
  }
  text(1+hist.spacer+(max(mutHist)+1)*.5*hist.unit.height,-hist.label.spacer,hist.label,cex=hist.label.cex,pos=1)
} else {
  Bottom.row <- 1-gene.width-mag.spacer-row.spec.width-tab.spacer-mutant.space
  arrows(1+hist.spacer,Bottom.row-hist.axis.spacer,1+hist.spacer+(max(mutHist)+1)*hist.unit.height,Bottom.row-hist.axis.spacer,length=2*hist.unit.height)
  for(i in 1:max(mutHist)) {
    lines(c(1+hist.spacer+(i)*hist.unit.height,1+hist.spacer+(i)*hist.unit.height),
          c(Bottom.row-hist.axis.spacer,Bottom.row-hist.axis.spacer-hist.ticks.length))
    text(1+hist.spacer+(i)*hist.unit.height,Bottom.row-hist.numbers.spacer,i,cex=hist.numbers.cex,pos=1)
  }
  text(1+hist.spacer+(max(mutHist)+1)*.5*hist.unit.height,Bottom.row-hist.label.spacer,hist.label,cex=hist.label.cex,pos=1)
}
#Add size labels 
text(0, 1, paste(0, "Mb"), cex=hist.label.cex, pos=3)
text(1, 1, paste(round(N.bases/1000000, digits = 1), "Mb"),cex=hist.label.cex, pos=3)

dev.off()



