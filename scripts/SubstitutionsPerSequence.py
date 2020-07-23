import pysam
import sys, math

try:
	bamfile=sys.argv[1]
	sval=sys.argv[2]
	#fastafile=sys.argv[2]
	estPercRead=sys.argv[3]
except: sys.exit('<BAM file> <strand orientation: 0, 1 or 2> <estimate read perc y/n>')


covidseq=''
f=open('/lustrehome/epicardi/home/covid/hyper_editing_scripts/covidgenome/nc045512.fa')
for i in f:
	if i.startswith('>'): continue
	if i.strip()=='': continue
	covidseq+=(i.strip()).upper()
f.close()

MAPQ=30
MQUAL=30
percRead=0.05 # 5% della read da condiserare
estpr=0
if estPercRead.lower() in ['y','ye','yes']: estpr=1

def baseCount(seq):
	s=seq.upper()
	sl=len(seq)*1.0
	for i in 'ACGT':
		if s.count(i)==0: return 1
		if s.count(i)/sl < 0.1: return 1
		if s.count(i)/sl > 0.6: return 1
	return 0

def mean(val):
	try: v=(sum(val)*1.0)/len(val)
	except: v=0.0
	return v

def makeCluster(allcoord):
	cluster=[]
	remaining=[]
	c1=allcoord[0][0]
	c2=allcoord[0][1]
	for i in range(len(allcoord)):
		if allcoord[i]!=(c1,c2):
			if c1<=allcoord[i][0]<=c2:
				cluster.append(allcoord[i])
				if allcoord[i][1]>c2:
					c2=allcoord[i][1]
			else:
				remaining.append(allcoord[i])
		else:
			cluster.append((c1,c2))
	return (c1,c2),remaining

def Clusters(interval):
	coords=[]
	interval.sort()
	while len(interval)!=0:
		coord,interval=makeCluster(interval)
		coords.append(coord)
	return coords
	
def getp(start,end,seq):
	al=(end-start)+1
	sl=len(seq)*1.0
	try: p=(al/sl)*100
	except: p=0.0
	return '%.1f' %(p) 

def comp(sub):
	d={'A':'T','T':'A', 'G':'C', 'C':'G', 'N':'N'}
	sub=sub.upper()
	return d[sub[0]]+d[sub[1]]

def comp2(sub):
	d={'A':'T','T':'A', 'G':'C', 'C':'G', 'N':'N'}
	sub=sub.upper()
	if d.has_key(sub): return d[sub]
	else: return 'N'

if sval not in ['0','1','2']: sys.exit('Strand orientation values: 0, 1 or 2')

s=[]
for i in 'ACGT':
	for j in 'ACGT':
		if i!=j: s.append(i+j)
s.sort()
d={}
reads={}
cls={}
for i in s: cls[i]=[0,[]]
outfile=bamfile[:-4]

good_clusters=[]
o=open(outfile+'.reads','w')
o2=open(outfile+'.cls','w')
o3=open(outfile+'.Badreads','w')
nreads=0
bam = pysam.AlignmentFile(bamfile, "rb" )
xx,pd=0,0
for read in bam.fetch("NC_045512.2", 1, 29903):
	if xx==1000: break
	if read.is_paired: pd+=1
	xx+=1
pe=1
if pd>0: pe=2
if pe==1: sys.stderr.write('BAM with Single End Reads.\n')
elif pe==2: sys.stderr.write('BAM with Paired End Reads.\n')

nucUP=[]
nucDW=[]
context=[]
for read in bam.fetch("NC_045512.2", 1, 29903):
	#if nreads==10: break
	#print read #SRR11550014.8368142
	if read.is_unmapped: continue
	if pe==2:
		if not read.is_paired: continue
		if not read.is_proper_pair: continue
	if read.mapping_quality < MAPQ: continue
	if not read.has_tag('NM'): continue
	nm=read.get_tag('NM')
	seq=read.query_sequence
	strand='+'
	isrev='0'
	if read.is_reverse: isrev='1'
	rt='2'
	if read.is_read1: rt='1'
	if sval=='2':
		if rt=='1':
			if read.is_reverse: strand='+'
			else: strand='-'
		elif rt=='2':
			if read.is_reverse: strand='-'
			else: strand='+'
	elif sval=='1':
		if rt=='1':
			if read.is_reverse: strand='-'
			else: strand='+'
		elif rt=='2':
			if read.is_reverse: strand='+'
			else: strand='-'
	else: strand='+'	
	name=read.query_name
	if not reads.has_key(name): reads[name]={'1':[],'2':[]}
	#print read.query_alignment_start+1, read.query_alignment_end
	#print read.reference_start+1, read.reference_end
	paln=getp(read.query_alignment_start+1, read.query_alignment_end,seq)
	if nm==0: 
		line=[name,	rt, read.reference_start+1, read.reference_end,read.query_alignment_start+1,read.query_alignment_end, paln, '0', '0', '0','0',strand,'0']
		o.write(' '.join([str(x) for x in line])+'\n')
		reads[name][rt].append((read.reference_start+1, read.reference_end,[],[]))
		continue
	pos,sub,rpos,q,underq,upn,dwn,subseq=[],[],[],[],0,[],[],[]
	qual=read.query_qualities
	rseq=read.query_sequence
	#filtro per la read
	if baseCount(rseq)==1: continue
	#sys.stderr.write(str(read)+'\n')
	aln=read.get_aligned_pairs(with_seq=True)
	dupdw={}
	#(1, 275, 'T')
	for j in aln:
		if j[2]==None: continue
		if j[0]==None: continue
		if j[1]==None: continue
		if j[2].islower()==True:
			ch=j[2].upper()+seq[j[0]] #base change
			if read.query_alignment_start+1 <=j[0]+1<= read.query_alignment_end:
				if qual[j[0]] < MQUAL:
					underq+=1
					q.append(qual[j[0]])
					continue
				pos.append(str(j[1]+1))
				if strand=='+': sub.append(ch)
				else: sub.append(comp(ch))
				rp=((j[0]+1)-(read.query_alignment_start+1))
				rpos.append(str(rp+1))
				q.append(qual[j[0]])
				#upn.append(rseq[(rp-1)])
				#try: dwn.append(rseq[(rp+1)])
				#except: dwn.append('N')
				#subseq.append(rseq[rp-5:rp+6])
				#base quality
		dupdw[j[1]+1]=j[2].upper()
			#print rt,j[1]+1,j[0]+1, ((j[0]+1)-(read.query_alignment_start+1))+1, ch
			#pileupread.alignment.query_sequence[pileupread.query_position].upper(),pileupread.alignment.query_qualities[pileupread.query_position]
	mq=mean(q)
	if len(pos)>0 and underq==0:
		line=[name,	rt, read.reference_start+1, read.reference_end , read.query_alignment_start+1,read.query_alignment_end, paln, ','.join(pos), ','.join(rpos), ','.join(sub),str(underq),strand,isrev,str(mq)]
		cigDic=dict(read.cigartuples)
		o.write(' '.join([str(x) for x in line])+'\n')
		#exclude reads with insertions and deletions
		if cigDic.has_key('1'): continue
		if cigDic.has_key('2'): continue
		if len(set(sub))==1 and len(sub)>=2:
			ralen=(read.query_alignment_end-(read.query_alignment_start+1))+1
			if estpr: percRead= 3.2/ralen
			else: percRead= 0.05
			if (len(sub)*1.0)/ralen < percRead: continue
			if float(paln) < 80: continue # 80% della parte allineata
			n_at_end=int(round(ralen*0.2)) # basi del 20% della read
			nrpos=[int(x) for x in rpos]
			nrpos.sort()
			#print line
			#print rseq
			if (read.query_alignment_start+1)<=nrpos[0]<=((read.query_alignment_start+1)+n_at_end)-1 and (read.query_alignment_start+1)<=nrpos[-1]<=((read.query_alignment_start+1)+n_at_end)-1: continue
			if (read.query_alignment_end - n_at_end)+1<=nrpos[0]<=read.query_alignment_end and (read.query_alignment_end - n_at_end)+1<=nrpos[-1]<=read.query_alignment_end: continue
			for k in pos:
				if strand=='-':
					try: dwn.append(comp2(dupdw[int(k)-1]))
					except: dwn.append(comp2(covidseq[int(k)-1]))
					try: upn.append(comp2(dupdw[int(k)+1]))
					except: upn.append(comp2(covidseq[int(k)+1]))
				else:
					try: upn.append(dupdw[int(k)-1])
					except: upn.append(covidseq[int(k)-1])
					try: dwn.append(dupdw[int(k)+1])
					except: dwn.append(covidseq[int(k)+1])
			if sub[0] in ['AG','TC']:
				nucUP+=upn
				nucDW+=dwn
				nch=sub[0]
				for k in range(len(pos)):
					context.append(pos[k]+' '+nch+' '+upn[k]+nch[0]+dwn[k])
			#print aln
			#rimozione di posizioni alle ends
			#newcc=((read.query_alignment_start+1)+5,read.query_alignment_end-5)
			#nrpos=[int(x) for x in rpos]
			#npos=[x for x in nrpos if newcc[0]<=x<=newcc[1]]
			#if len(npos)>=2:
			#print line
			#print [rseq.count('A'),rseq.count('C'),rseq.count('G'),rseq.count('T')]
			o2.write(' '.join([str(x) for x in line])+'\n')
			cls[sub[0]][0]+=1 
			cls[sub[0]][1].append(len(sub))
			good_clusters.append(line)
		elif len(set(sub)) > 2 and len(sub) > 2:
			o3.write(name+'\n')
		for k in range(len(rpos)):
			if d.has_key(int(rpos[k])): d[int(rpos[k])].append(sub[k])
			else: d[int(rpos[k])]=[sub[k]]
	elif len(pos)>0 and underq>0:
		line=[name,	rt, read.reference_start+1, read.reference_end , read.query_alignment_start+1,read.query_alignment_end, paln, ','.join(pos), ','.join(rpos), ','.join(sub),str(underq),strand,str(mq)]
		o.write(' '.join([str(x) for x in line])+'\n')
		if len(set(sub)) > 2 and len(sub) > 2:
			o3.write(name+'\n')
	elif len(pos)==0 and underq>0:
		line=[name,	rt, read.reference_start+1, read.reference_end , read.query_alignment_start+1,read.query_alignment_end, paln,  '0', '0', '0', str(underq),strand,str(mq)]
		o.write(' '.join([str(x) for x in line])+'\n')
	else:
		line=[name,	rt, read.reference_start+1, read.reference_end,read.query_alignment_start+1,read.query_alignment_end, paln, '0', '0', '0','0',strand,str(mq)]
		o.write(' '.join([str(x) for x in line])+'\n')
	reads[name][rt].append((read.reference_start+1, read.reference_end,[int(x) for x in pos],sub))

o.close()
o3.close()
o2.write('\nCluster Distro per base change:\n')
for i in s:
	o2.write('%s %i %s\n' %(i,cls[i][0],','.join([str(x) for x in sorted(cls[i][1])])))
o2.close()

o=open(outfile+'.dis','w')
o.write('Pos '+' '.join(s)+' All\n')
for i in sorted(d):
	v=[d[i].count(x) for x in s]
	all=sum(v)
	o.write(str(i)+' '+' '.join([str(x) for x in v])+' '+str(all)+'\n')

o.write('\n>overlap\n')
nreadov=0
allbov=0.0
poverr=[]
for i in reads:
	if reads[i]['1']==[] or reads[i]['2']==[]: continue 
	r1,r2=reads[i]['1'][0],reads[i]['2'][0]
	if r1[0]<=r2[0]<=r1[1] or r1[0]<=r2[1]<=r1[1]:
		nreadov+=1
		#print i,r1,r2
		if r1[0]<=r2[0]<=r1[1]: #and not r1[0]<=r2[1]<=r1[1]: 
			st=r2[0]
			if r2[0]<=r1[1]<=r2[1]: ed=r1[1]
			else: ed=r2[1]
			bov=(ed-st)+1
			allbov+=bov
			pov={}			
			for j in range(len(r1[2])):
				if st<=r1[2][j]<=ed:
					if pov.has_key(r1[2][j]): pov[r1[2][j]].append(r1[3][j])
					else: pov[r1[2][j]]=[r1[3][j]]
			for j in range(len(r2[2])):
				if st<=r2[2][j]<=ed:
					if pov.has_key(r2[2][j]): pov[r2[2][j]].append(r2[3][j])
					else: pov[r2[2][j]]=[r2[3][j]]
			for k in pov:
				if len(pov[k])==1: poverr.append(pov[k][0])
		elif r1[0]<=r2[1]<=r1[1]: #and not r1[0]<=r2[0]<=r1[1]:
			if r2[0]<=r1[0]<=r2[1]: st=r1[0]
			else: st=r2[0]
			ed=r2[1]
			bov=(ed-st)+1
			allbov+=bov
			pov={}
			for j in range(len(r1[2])):
				if st<=r1[2][j]<=ed:
					if pov.has_key(r1[2][j]): pov[r1[2][j]].append(r1[3][j])
					else: pov[r1[2][j]]=[r1[3][j]]
			for j in range(len(r2[2])):
				if st<=r2[2][j]<=ed:
					if pov.has_key(r2[2][j]): pov[r2[2][j]].append(r2[3][j])
					else: pov[r2[2][j]]=[r2[3][j]]
			for k in pov:
				if len(pov[k])==1: poverr.append(pov[k][0])
			#print i,r1,r2,st,ed

try: err=len(poverr)/allbov
except: err=0.0
o.write('Reads in overlap: %i\n' %(nreadov))
o.write('Bases in overlap: %i\n' %(allbov))
o.write('Errors in overlap: %i\n' %(len(poverr)))
o.write('Error Rate: %.6f\n' %(err))
o.write('Error types:\n')
for i in s:
	try: errval=(poverr.count(i)*1.0)/allbov
	except: errval=0.0
	o.write('%s %i %.6f\n' %(i,poverr.count(i),errval))
o.close()

o4=open(outfile+'.sum.cls','w')

up=[nucUP.count('A'), nucUP.count('C'), nucUP.count('G'), nucUP.count('T')]
dw=[nucDW.count('A'), nucDW.count('C'), nucDW.count('G'), nucDW.count('T')]

n='ACGT'
x=0
for i in up:
	try: co=(i*1.0)/sum(up)
	except: co=0.0
	o4.write('%s_up %.6f\n' %(n[x],co))
	x+=1
x=0
for i in dw:
	try: co=(i*1.0)/sum(dw)
	except: co=0.0
	o4.write('%s_up %.6f\n' %(n[x],co))
	x+=1

#good_clusters=[]
#f=open('SRR11550045.cls')
#for i in f:
#	if i.strip()=='': break
#	l=(i.strip()).split()
#	good_clusters.append(l)
#f.close()

nsubs={'+':{},'-':{}}

for i in 'ACGT':
	for j in 'ACGT':
		if i!=j: 
			nsubs['+'][i+j]=[[],[]]
			nsubs['-'][i+j]=[[],[]]
			
subs=[x for x in sorted(nsubs['+'].keys())]
#print nsubs
for i in good_clusters:
	cc=[int(x) for x in i[7].split(',') if x!='']
	cc.sort()
	coord=(cc[0],cc[-1])
	s=i[9].split(',')[0]
	strand=i[11]
	#print coord,s,strand
	nsubs[strand][s][0].append(coord)
	nsubs[strand][s][1]+=cc

for i in subs:
	if len(nsubs['+'][i][0])>0:
		cls=Clusters(nsubs['+'][i][0])
		pp=set(nsubs['+'][i][1])
		for j in cls:
			ne=[x for x in pp if j[0]<=x<=j[1]]
			r=['>',i,j[0],j[1],'+',(j[1]-j[0])+1,len(ne),','.join([str(x) for x in ne])]
			o4.write(' '.join([str(x) for x in r])+'\n')
	if len(nsubs['-'][i][0])>0:
		cls=Clusters(nsubs['-'][i][0])
		pp=set(nsubs['-'][i][1])
		for j in cls:
			ne=[x for x in pp if j[0]<=x<=j[1]]
			r=['>',i,j[0],j[1],'-',(j[1]-j[0])+1,len(ne),','.join([str(x) for x in ne])]
			o4.write(' '.join([str(x) for x in r])+'\n')
o4.close()

o5=open(outfile+'.context','w')
for i in context:
	o5.write(i+'\n')
o5.close()
