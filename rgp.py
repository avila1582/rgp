import sys
group=sys.argv[1]
def sgrp(gl,ori,dnaA,gi):
	strand={True:"D",False:"G"}
	g=sum(gi[:2])/2
	if g>ori:
		gori=g-ori
	else:
		gori=g-ori+gl
	if abs(gori)>gl/2:
		gori-=gl/2
	if abs(gori)>gl:
		gori+=gl
	dnaAstrand=dnaA[2]=="+"
	gstrand=gi[2]=="+"
	if dnaAstrand:
		if ori>g:
			gori=gl-ori+g
		else:
			gori=g-ori
	else:
		gstrand=not gstrand
		if ori>g:
			gori=ori-g
		else:
			gori=ori+gl-g
	if gori>gl/2:
		gori=gori-gl
		gstrand= not gstrand
	val=round(gori/gl*200,2)
	return(str(val)+strand[gstrand])
#group="AT21"
chroms=[]
lens={}
oris={}
dnaAs={}
for line in open("results/"+group+"/downloaded_ori.csv","r"):
	cols=[a.strip().strip("\"") for a in line.split("\t")]
	print(cols)
	gl=int(cols[3])
	oril=[int(a.replace(",","")) for a in cols[8].split(" ... ")]
	if oril[1]<oril[0]:
		ori=oril[0]+int(((gl-oril[0])+oril[1])/2)
		if ori>=gl:
			ori-=gl
	else:
		ori= int(sum(oril)/2)		
	lens[cols[1]]=gl
	print(cols[1])
	print(ori)
	oris[cols[1]]=ori
	dnaAs[cols[1]]=[int(a) for a in cols[4:6]]+[cols[6]]
	chroms.append(cols[1])
for ch  in chroms:
	w=open("rgp/"+group+"/"+ch+".rgp","w")
	for line in open("gff/"+group+"/"+ch+".gff","r"):
		if line[0]!="#":
			cols=[a.strip().strip("\"") for a in line.split("\t")]
			if cols[2]=="gene" or cols[2]=="pseudogene":
				cigar=cols[8].split(";")
				vals=[x.split("=")[1] for x in cigar if x.find("=")!=-1]
				keys=[x.split("=")[0] for x in cigar if x.find("=")!=-1]
				if "locus_tag" in keys:
					lctg= vals[keys.index("locus_tag")]
					rt=[int(a) for a in cols[3:5]]+[cols[6]]
					w.write(lctg+"\t"+sgrp(lens[ch],oris[ch],dnaAs[ch],rt)+"\n")

