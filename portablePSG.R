###PSG Analysis 簡易PSG用
###2015/07/23
###2015/12/13 revise
###All rights reserved by Taro Ueno

#PSGの時系列データ（睡眠stage）からWake, N1, N2, N3, REMの累積分布のグラフを描く
#IN:行が時間（30秒エポック）、1列の時系列データ
#aにデータをimportしておく ex) a <- matrix(scan("test.csv", what="character"), ncol=1, byrow=T)
#WやRなどの文字なので、scanの中に、what="character"を入れる


#PSGの時系列データを元に、それぞれのステージを1で表す行列を作成
#1列目：Wake、2列目：N1、3列目：N2、4列目：N3、5列目：Rの行列で出力
#それぞれのstageを1に、他を0とする。
stage <- function(x,y) {
	y <- matrix(0, ncol=5, nrow=length(x[,1]))

		for (i in 1:length(x[,1])) {
			if (x[i,1] == "W") {
			y[i,1] <- 1
			}
			if (x[i,1] == "1") {
			y[i,2] <- 1
			}
			if (x[i,1] == "2") {
			y[i,3] <- 1
			}
			if (x[i,1] == "3") {
			y[i,4] <- 1
			}
			if (x[i,1] == "R") {
			y[i,5] <- 1
			}
			}

			return(y)
			}


#1の持続する時間=bout lengthをチェック　持続的に
Rbout1 <- function(x) {
	for (m in 1:length(x[1,])) {
				for (n in 2:length(x[,m])) {
					if (x[n,m] == 1) {
						x[n,m] <- x[n-1,m] + 1
						}
							}
							}
							return(x)
							}


#持続時間のみに変換
Rbout2 <- function(x) {
	for (m in 1:length(x[1,])) {
		for (n in 1:length(x[,m])-1) {
			if (x[n+1,m] != 0) {
				x[n,m] <- 0
				}
				}
				}
				return(x)
				}


#0を除去するために使用
#Rcumのために使用
OnlyR <- function(x) {
	for (i in 1:length(x[1,])) {
		a <- matrix(x[,i], ncol = 1)
		b <- matrix(a[order(a[,1], decreasing = T), ], ncol = 1)
		x[,i] <- b[,1]
		}
		return(x)
		}


#0を削除し行列を圧縮
compress <- function(x) {
	for (i in 1:length(x[,1])) {
		if  ((max(x[i,]) >0) && (max(x[i+1,]) == 0)) {
			a <- i
			}
			}
			x <- x[1:a, ]

			return(x)
			}


#Rcumのver3
Rcum3 <- function(x) {
	a <- matrix(0, ncol = length(x[1,]), nrow = max(x)) #ゼロ行列を作成

	for (m in 1:length(x[1,])) {
				for (n in 1:length(x[,m])) {
					for (o in 1:length(x[,m])) {
					if (x[o,m] >= n) {
						a[n,m] <- o
						}
						}
						}
						}
						return(a)
						}


#transition数を計算　行にtransition元、列にtransition先 W, N1, N2, N3, Rの順
#stage関数で作ったものをxに、compress関数で作ったものをyに入れる
transition <- function(x,y) {
	a <- matrix(0, ncol=5, nrow=5)

	for (j in 1:5){
		for (i in 2:length(x[,1])){
			if(x[i,j]==0 && x[i-1,j]==1){
				if(x[i,1]==1){a[j,1] <- a[j,1]+1}
				if(x[i,2]==1){a[j,2] <- a[j,2]+1}
				if(x[i,3]==1){a[j,3] <- a[j,3]+1}
				if(x[i,4]==1){a[j,4] <- a[j,4]+1}
				if(x[i,5]==1){a[j,5] <- a[j,5]+1}
				}
				}
				}
				colnames(a) <- c("Wake", "N1", "N2", "N3", "REM")
				rownames(a) <- c("Wake->", "N1->", "N2->", "N3->", "REM->")
				return(a)
				}

#normalized transition rateを計算　行にtransition元、列にtransition先 W, N1, N2, N3, Rの順
#stage関数で作ったものをxに、compress関数で作ったものをyに入れる
transitionrate <- function(x,y){
		a <- matrix(0, ncol=5, nrow=5)

	for (j in 1:5){
		for (i in 2:length(x[,1])){
			if(x[i,j]==0 && x[i-1,j]==1){
				if(x[i,1]==1){a[j,1] <- a[j,1]+1}
				if(x[i,2]==1){a[j,2] <- a[j,2]+1}
				if(x[i,3]==1){a[j,3] <- a[j,3]+1}
				if(x[i,4]==1){a[j,4] <- a[j,4]+1}
				if(x[i,5]==1){a[j,5] <- a[j,5]+1}
				}
				}
				a[j,] <- a[j,]/sum(a[j,]) #normalize
				}
				colnames(a) <- c("Wake", "N1", "N2", "N3", "REM")
				rownames(a) <- c("Wake->", "N1->", "N2->", "N3->", "REM->")
				return(a)
				}


#Mignotの判定式（3N2/N3→2N1/W, 6epoch以上のN1/W bout, 5N1/W→2REM）
#stage関数で作ったものをxに、compress関数で作ったものをyに入れる
Mignot <- function(x,y) {
	a <- matrix(0, ncol=3)

	b <- matrix(0, ncol=3, nrow=length(x[,1]))#1列目がN2/N3, 2列目がN1/W, 3列目がREM
	b[,1] <- apply(x[,3:4],1,sum)
	b[,2] <- apply(x[,1:2],1,sum)
	b[,3] <- x[,5]
	for(i in 5:length(b[,1])){
		if(b[i,2]==1 && b[i-1,2]==1 && b[i-2,1]==1 && b[i-3,1]==1 && b[i-4,1]==1) {
			a[1,1] <- a[1,1] + 1
			}
			}

	a[1,2] <- length(which(y[,1:2] >=6))

	for(i in 7:length(b[,1])){
		if(b[i,3]==1 && b[i-1,3]==1 && b[i-2,2]==1 && b[i-3,2]==1 && b[i-4,2]==1 && b[i-5,2]==1 && b[i-6,2]==1) {
			a[1,3] <- a[1,3] + 1
			}
			}
	colnames(a) <- c("3N2/N3->2N1/W", "N1/W>=6", "5N1/W->2R")
	return(a)
	}

#以下、実際の解析
#フォルダ内のxlsxファイルからヘッダーを除外し、csvファイルを作成
#フォルダにまとめた全xlsxファイルに対して解析を実施
library(xlsx)

files <- list.files()
for(file.name in files){
	if(regexpr('¥¥.xlsx$', file.name) <0){

data <- read.xlsx(file.name,1)
data1 <- cbind(as.matrix(data[,2]), as.matrix(data[,3]))
data2 <- data1[5:length(data1[,1]),] #最初のヘッダー部分を除く

s <- matrix(which(data2[,2] == "W"), ncol=1) ##最初のWを除く
t <- rbind(s,0) - rbind(0,s)
data3 <- data2[min(which(t[,1] != 1)):length(data2[,1]),]

for(i in 1:(length(data3[,1])-1)){ ##重複しているエポックにNAを入れる　重複がない場合用にelse breakを追加した(151213)
	if(data3[i,1]==data3[i+1,1]){
		data3[i,1] <- NA
		}
		else{
			break
			}
		}
data4 <- na.omit(data3) #重複しているエポックを除く

file.name2  <- sub('\\.[^.]*', ".csv", file.name) # ピリオド以降を ".csv" に置換

write.csv(file=file.name2, data4,  row.names=F)##入眠前のWは残っている
}

a1 <- matrix(scan(file.name2, what="character", sep = ","), ncol=2, byrow=T)
a <- matrix(a1[2:length(a1[,1]),2], ncol=1)
b <- stage(a)
c <- Rbout1(b)
d <- Rbout2(c)
e <- OnlyR(d)
f <- compress(e)
g <- Rcum3(f)
h <- transition(b,f)
i <- Mignot(b,f)

#テキストで出力
head <- c("Wake", "N1", "N2", "N3", "REM")
file.name3  <- sub('\\.[^.]*', "PSGanalysis.tsv", file.name) # ピリオド以降を "PSGanalysis.tsv" に置換
write.table(g, file = file.name3, quote = F, row.names = F, sep = "\t", col.names=head)
file.name4  <- sub('\\.[^.]*', "TransitionRate.tsv", file.name) # ピリオド以降を "TransitionRate.tsv" に置換
write.table(h, file = file.name4, quote = F, row.names = T, sep = "\t", col.names=T)
file.name5  <- sub('\\.[^.]*', "Mignot.tsv", file.name) # ピリオド以降を "Mignot.tsv" に置換
write.table(i, file = file.name5, quote = F, row.names = F, sep = "\t", col.names=T)

#グラフ化
par(mfcol = c(2,3)) #8行1列の形でプロット出力
par(mai = c(0.5, 0.5, 0.5, 0.5)) #余白設定
file.name6  <- sub('\\.[^.]*', ".pdf", file.name) # ピリオド以降を ".pdf" に置換
pdf(file.name6)
plot(g[,1]/g[1,1], xlab="epoch", ylab="cum prob", main="Wake", ty="l")
plot(g[,2]/g[1,2], xlab="epoch", ylab="cum prob", main="N1", ty="l")
plot(g[,3]/g[1,3], xlab="epoch", ylab="cum prob", main="N2", ty="l")
plot(g[,4]/g[1,4], xlab="epoch", ylab="cum prob", main="N3", ty="l")
plot(g[,5]/g[1,5], xlab="epoch", ylab="cum prob", main="REM", ty="l")
dev.off()
dev.off()
}
