#Random

data=c()
for (i in 1:cells) {
	data=sample(1:positions, ribos, replace=TRUE)
	vector <- append(data,sample(1:positions, ribos, replace=FALSE))
}

hx=hist(vector,breaks=seq(1,positions,l=positions+1),plot=FALSE)
plot(hx$counts)
