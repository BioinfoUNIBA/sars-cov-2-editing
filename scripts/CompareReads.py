import pysam
import sys, math

try:
	rfile1=sys.argv[1]
	rfile2=sys.argv[2]
	badfile=sys.argv[3]
except: sys.exit('<BWA.reads> <GSNAP.reads> <output bad file base>')


def getReads(infile):
	d={}
	nr,pe=0,0
	f=open(infile)
	for i in f:
		l=(i.strip()).split()
		if d.has_key(l[0]):
			d[l[0]][l[1]]=[int(l[2]),int(l[3]),int(l[4]),int(l[5]),float(l[6]),[int(x) for x in l[7].split(',')],[int(x) for x in l[8].split(',')]]
			pe+=1
		else:
			d[l[0]]={'1':[],'2':[]}
			d[l[0]][l[1]]=[int(l[2]),int(l[3]),int(l[4]),int(l[5]),float(l[6]),[int(x) for x in l[7].split(',')],[int(x) for x in l[8].split(',')]]
			nr+=1
	f.close()
	return d,nr,pe

pe=1
bwa,bwa_r,bwa_pe=getReads(rfile1)
if bwa_pe==0: pe=0
sys.stderr.write('BWA: reads %i - PE: %i\n' %(bwa_r,bwa_pe))
gsnap,gsnap_r,gsnap_pe=getReads(rfile2)
sys.stderr.write('GSNAP: reads %i - PE: %i\n' %(gsnap_r,gsnap_pe))

onlybwa=set(bwa.keys()) - set(gsnap.keys())
onlygsnap=set(gsnap.keys()) - set(bwa.keys())
com=set(bwa.keys()) & set(gsnap.keys())

sys.stderr.write('BWA only: %i\n' %(len(onlybwa)))
sys.stderr.write('GSNAP only: %i\n' %(len(onlygsnap)))
sys.stderr.write('Common reads: %i\n' %(len(com)))

br,bp=0,0
o=open(badfile+'.badReads','w')
o1=open(badfile+'.badPos','w')
for i in onlybwa: o.write(i+'\n')
for i in onlygsnap: o.write(i+'\n')
br+=len(onlybwa)
br+=len(onlygsnap)
for i in com:
	if bwa.has_key(i) and gsnap.has_key(i):
		bw=bwa[i]
		gn=gsnap[i]
		if pe==1:
			if len(bw['1'])==0 or len(bw['2'])==0: 
				o.write(i+'\n')
				br+=1
				continue
			if len(gn['1'])==0 or len(gn['2'])==0: 
				o.write(i+'\n')
				br+=1
				continue
			#if bw['1'][5][0]==0: continue
			#print bw
			#print gn
			ov1,ov2=0,0
			r1,r2=bw['1'],gn['1']
			if r1[0]<=r2[0]<=r1[1] or r1[0]<=r2[1]<=r1[1]: ov1+=1
			r3,r4=bw['2'],gn['2']
			#print r3,r4
			if r3[0]<=r4[0]<=r3[1] or r3[0]<=r4[1]<=r3[1]: ov2+=1
			#print ov1,ov1	
			if ov1==1 and ov2==1:
				p1=set(r1[5])
				p2=set(r2[5])
				onlyp1=p1 - p2
				onlyp2=p2 - p1
				bpos=[x for x in onlyp1 if x!=0]+[x for x in onlyp2 if x!=0]
				for j in bpos: o1.write(str(j)+' '+i+' 1\n')
				bp+=len(bpos)
				p3=set(r3[5])
				p4=set(r4[5])
				onlyp3=p3 - p4
				onlyp4=p4 - p3
				bpos2=[x for x in onlyp3 if x!=0]+[x for x in onlyp4 if x!=0]
				for j in bpos2: o1.write(str(j)+' '+i+' 2\n')
				bp+=len(bpos2)
			else:
				o.write(i+'\n')
				br+=1
		elif pe==0: #single end
			rt='1'
			if len(bw['1'])>0: r1,rt=bw['1'],'1'
			elif len(bw['2'])>0: r1,rt=bw['2'],'2'
			if len(gn['1'])>0: r2,rt=gn['1'],'1'
			elif len(gn['2'])>0: r2,rt=gn['2'],'2'
			ov1=0
			if r1[0]<=r2[0]<=r1[1] or r1[0]<=r2[1]<=r1[1]: ov1+=1
			if ov1==1:
				p1=set(r1[5])
				p2=set(r2[5])
				onlyp1=p1 - p2
				onlyp2=p2 - p1
				bpos=[x for x in onlyp1 if x!=0]+[x for x in onlyp2 if x!=0]
				for j in bpos: o1.write(str(j)+' '+i+' %s\n' %(rt))
				bp+=len(bpos)				
			else:
				o.write(i+'\n')
				br+=1
	else:
		o.write(i+'\n')	
		br+=1

o.close()
sys.stderr.write('Bad reads: %i\n' %(br))
sys.stderr.write('Bad sites: %i\n' %(bp))

