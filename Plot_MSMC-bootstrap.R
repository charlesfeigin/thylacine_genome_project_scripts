
library(ggplot2)

# Parse MSMC output, with output prefix pfx, scaling for mutation rate mu
# and generations g. mu has units substitutions/base/gen.
parse_msmc <- function(pfx, g=3, ref.len=0.97*3027.64e6, chr.arms=12, round=50) {

	file.final <- sprintf("%s.final.txt", pfx)
	file.loop <- sprintf("%s.loop.txt", pfx)
	file.log <- sprintf("%s.log", pfx)

	if (!file.exists(file.final) || !file.exists(file.loop) || !file.exists(file.log)) {
		return (NA)
	}

	loops <- read.table(file.loop, header=F)

	loops <- read.table(file.loop, header=F, stringsAsFactors=F,
			  col.names=c("rho", "ll", "time_boundaries", "lambdas"))
	time_boundaries <- lapply(strsplit(loops$time_boundaries, ","), as.double)
	lambdas <- lapply(strsplit(loops$lambdas, ","), as.double)
	rounds <- length(lambdas)

	if (round > rounds) {
		print(sprintf("%s: only found %d rounds!", file.loop, rounds))
		return (NA)
	}

	n <- length(lambdas[[1]])
	time_boundary <- time_boundaries[[round]][1:n]
	lambda <- lambdas[[round]]
	rho <- loops$rho[round]

	f <- file(file.log, "r")
	while (length(line <- readLines(f, n=1)) > 0) {
		if (grepl("^mutationRate", line)) {
			fields <- strsplit(line, "[ ]+")[[1]]
			theta <- as.numeric(fields[2])
			break
		}
	}
	close(f)

	nu = theta / rho
	mu = nu * chr.arms / ref.len

	time <- g * time_boundary / mu
	N <- 1.0 / (2*mu*lambda)

	return (data.frame(time, N))
}

dir1 <- "/localscratch/AussieGenomes.msmc_stuff/Thylacine/paleomix_MEM/"
dir2 <- "/localscratch/AussieGenomes.msmc_stuff/TasmanianDevil/paleomix_MEM/"
dirlist <- c(dir1, dir2)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colours <- cbbPalette
labels <- c("Thylacine", "Devil")
plt <- ggplot()
nreps <- 100 # number of bootstrap replicates

rounds <- c(50)
#rounds <- c(50, 20)

if (length(rounds) > 1) {
	old_labels <- labels
	li <- 0
	for (r in rounds) {
		for (l in old_labels) {
			li <- li +1
			labels[li] <- sprintf("%s (i=%d)", l,r)
		}
	}
}

li <- 0
for (r in rounds) {
	for (d in 1:length(dirlist)) {
		li <- li +1
		label = labels[li]
		dir = dirlist[d]

		path <- sprintf("%s/msmc/msmc1.est-r.i50.out", dir)
		m <- parse_msmc(path, round=r)
		if (!inherits(m, "data.frame")) {
			print(sprintf("%s.*: file(s) not found", path))
			return
		}
		m$label = label
		plt <- plt + geom_step(data=m, mapping=aes(time,N,colour=label), size=1)

	########for (n in 0:(nreps-1)) {
	########	path <- sprintf("%s/msmc_bootstrap/bootstrap_%d/msmc1.out", dir, n)
	########	m <- parse_msmc(path, round=r)
	########	if (!inherits(m, "data.frame")) {
	########		print(sprintf("%s.*: file(s) not found", path))
	########		next
	########	}
	########	m$label = label
	########	plt <- plt + geom_step(data=m, mapping=aes(time,N,colour=label), alpha=0.05)
	########}
	}
}

breaks <- function(a) {x<-c(); y<-10; while(y<a[1]) y<-y*10; while(y<=a[2]) {x<-c(x,y); y<-y*10} ; return(x)}
mbreaks <- function(a) {x<-c(); y<-10; while(y<a[1]) y<-y*10; while(y<=a[2]) {x<-c(x,seq(y*2,y*10,y)); y<-y*10} ; return(x)}
lfunc <- function(a) parse(text=sprintf("10^%.0f",log10(a)))
plt <- plt +
	xlab("time (yBP)") +
	ylab(expression(N[e])) +
	scale_x_log10(minor_breaks=mbreaks, breaks=breaks, labels=lfunc) +
	scale_y_log10(minor_breaks=mbreaks, breaks=breaks, labels=lfunc, limits=c(NA, 10^8)) +
	annotation_logticks() +
#	annotation_logticks(short=unit(-1,"mm"), mid=unit(-2,"mm"), long=unit(-3,"mm")) + # outward pointing tick marks are broken for logticks
	scale_colour_manual("", breaks=labels, values=colours) +
	theme_bw() +
	theme(panel.grid.major=element_blank(),
	      panel.grid.minor=element_blank(),
	      legend.background=element_blank(),
	      #legend.background=element_rect(colour='black'),
	      legend.title=element_blank(),
	      legend.justification=c(0,1),
	      legend.position=c(0.02,1.00),
	      legend.key.width=unit(2, "line"),
	      text=element_text(size=20))

#fn <- sprintf("Thylacine-Devil.bootstrap%s.pdf", c("", sprintf("-%dand%diterations",rounds[1],rounds[2]))[length(rounds)])
fn <- "Thylacine-Devil.pdf"
ggsave(fn, width=12, height=8, units="in", paper="a4r")
