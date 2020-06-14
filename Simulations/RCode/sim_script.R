wd = '/home/ZijunKe/Research/MetaR/Simulation4/RCode'
setwd(wd)
fid <- 0
for(sdPhii in 1:2){
for(Ri in 1:3){
for(SRi in 1:1){
for(RRi in 1:2){
for(NSi in 1:3){
for(Nbari in 1:2){
	fid <- fid +1
	newname <- paste('MetaRR',fid,'.R',sep='')
	repfolder <- paste("sed 's/RiR/", Ri,"/g;
	s/sdPhiR/",sdPhii,"/g;
	s/SRiP/",SRi,"/g;
	s/RRiD/",RRi,"/g;
	s/NSiR/",NSi,"/g;
	s/NbariR/",Nbari,"/g' RunCRC.R > ", newname, "",sep='')
	system(repfolder)
	s.name <- paste('VL',sdPhii-1,'VPhi',sdPhii-1,sep='')
	cond.name <- paste('R',Ri,'SR',SRi,'RR',RRi,'NS',NSi,'Nbar',Nbari,sep='')
	write.table(t(c(fid,s.name,cond.name)),'FID.txt',append = T,row.names = F,col.names = F,quote = F)
}
}
}
}
}
}
