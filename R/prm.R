prm <-
function(X,y,a,fairct=4,opt="l1m"){

# PRM Partial Robust M-regression estimator

n=nrow(X)
p=ncol(X)

if (p>n){
	dimensions <- 1
	dimension <- p-n
	ressvd <- svd(t(X))
	X <- ressvd$v%*%diag(ressvd$d)
	n=nrow(X)
	p=ncol(X)
}
else {dimensions <- 0}

if (opt=="l1m"){
	require(pcaPP)
	mx <- l1median(X)   
}
else { mx <- apply(X,2,median)}
my <- median(y)

Xmc <- scale(X,center=mx,scale=FALSE)
wx <- sqrt(apply(Xmc^2,1,sum))
wx <- wx/median(wx)
wx <- 1/((1+abs(wx/fairct))^2)

wy <- abs(y-my)
wy <- wy/median(wy)
wy <- 1/((1+abs(wy/fairct))^2)

w <- wx*wy
Xw <- X*sqrt(w)
yw <- y*sqrt(w)

loops <- 1
ngamma <- 10^5
difference <- 1

require(pls)
while ((difference>0.01) && loops<30){
	ngammaold <- ngamma
	spls <- mvr(yw~Xw,ncomp=a,method="simpls")
	b <- spls$coef[,,a]
	gamma <- t(t(y)%*%spls$sco)
	T <- spls$sco/sqrt(w)
	r <- y-T%*%gamma
	rc <- r-median(r)
	r <- rc/median(abs(rc))
	wy <- 1/((1+abs(r/fairct))^2)
	if (opt=="l1m"){mt <- l1median(T)}
	else {mt <- apply(T,2,median)}
	dt <- T-mt
	wt <- sqrt(apply(dt^2,1,sum))
	wt <- wt/median(wt)
	wt <- 1/((1+abs(wt/fairct))^2)
	ngamma <- sqrt(sum(gamma^2))
	difference <- abs(ngamma-ngammaold)/ngamma
	w <- drop(wy*wt)
	Xw <- X*sqrt(w)
	yw <- y*sqrt(w)
#print(difference)
#print(loops)
	loops <- loops+1
}

if (dimensions==1){
	b <- drop(ressvd$u%*%b)
	yfit <- as.matrix(X)%*%t(ressvd$u)%*%b
}
else {
	yfit <- as.matrix(X)%*%b
}
list(coef=b,wy=wy,wt=wt,scores=T,loadings=spls$loadings,
	fitted.values=yfit)
}

