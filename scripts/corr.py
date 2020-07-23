import pysam
import sys

try:
	bamfile=sys.argv[1]
	redifile=sys.argv[2]
	fastafile=sys.argv[3]
	stranded=sys.argv[4]
	badposfile=sys.argv[5]
	badreadfile=sys.argv[6]
except: sys.exit('<BAM file> <reditools file> <reference fasta> <strand 1,2,0> <.badPos file> <.badReads files separated by commas>')

if stranded=='0':
	unchange1=0
	unchange2=0
	getstrand=0
	import scipy.stats as stats
elif stranded=='1':
	unchange1=1
	unchange2=0
	getstrand=1
elif stranded=='2':
	unchange1=0
	unchange2=1
	getstrand=1
else:
	unchange1=0
	unchange2=0
	getstrand=0
	import scipy.stats as stats

def comp(s):
	a={'A':'T','T':'A','C':'G','G':'C'}
	ss=''
	for i in s.upper():
		if a.has_key(i): ss+=a[i]
		elif i==' ': ss+=' '
		elif i=='-': ss+='-'
		else: ss+='N'
	return ss	
	
def infstr(strand):
	p=strand.count('1')
	m=strand.count('0')
	try: percp=float(p)/(p+m)
	except: percp=0.0
	percm=1.0-percp
	if percp>=0.7: return '1'
	elif percm>=0.7: return '0'
	else: return '2'

def getUpDw(seqlen):
	if 50 <= seqlen <= 75: return 3,5
	elif 76 <= seqlen <= 99: return 5,6
	elif 100 <= seqlen <= 101: return 5,10
	elif 102 <= seqlen <= 121: return 8,12
	elif 121 <= seqlen <= 151: return 10,15
	elif seqlen >= 152: return 15,15
	return 5,5

singleRep=(29871,29903)

br={}
for file in badreadfile.split(','):
	f=open(file)
	for i in f: br[i.strip()]=0
	f.close()

bp={}
f=open(badposfile)
for i in f:
	l=(i.strip()).split()
	pp=int(l[0])
	if bp.has_key(pp):
		bp[pp][(l[1],int(l[2]))]=0
		#bp[pp].append((l[1],int(l[2])))
	else: 
		bp[pp]={}
		bp[pp][(l[1],int(l[2]))]=0
f.close()

MAPQ=30
MQUAL=30
MINCOV=20
MAX_DEPTH=100000

kv={} #conosciute varianti
f=open('/lustrehome/epicardi/home/covid/knownVariants.txt')
for i in f:
	l=(i.strip()).split('\t')
	kv[int(l[3])]=0
f.close()

d={}
f=open(redifile)
for i in f:
	if i.startswith('Region'): continue
	l=(i.strip()).split('\t')
	d[int(l[1])]=l
f.close()

bam = pysam.AlignmentFile(bamfile, "rb" )
xx,pd=0,0
for read in bam.fetch("NC_045512.2", 1, 29903):
	if xx==1000: break
	if read.is_paired: pd+=1
	xx+=1

pe=0
if pd>0: pe=1
if pe==0: sys.stderr.write('BAM with Single End Reads.\n')
elif pe==1: sys.stderr.write('BAM with Paired End Reads.\n')

o=open(redifile[:-4]+'.corr.edi','w')
x=0
#bam = pysam.AlignmentFile(bamfile, "rb" )
fasta=pysam.Fastafile(fastafile)
for pileupcolumn in bam.pileup("NC_045512.2", 1, 29903, stepper='nofilter', ignore_overlaps=False, max_depth=MAX_DEPTH):
	pos=pileupcolumn.pos+1
	if not d.has_key(pos): continue
	if kv.has_key(pos): continue
	if d[pos][7]=='': continue
	#if pos!=5710: continue #26371
	if singleRep[0]<=pos<=singleRep[1]: continue #elimina pos in simple repeats
	ref=fasta.fetch("NC_045512.2",pileupcolumn.reference_pos,pileupcolumn.reference_pos+1).upper()
	badpos = {}
	if bp.has_key(pos): badpos=bp[pos]
	#if badpos!={}:
	#print badpos
	rt=1
	seq,names,subs,reads1mod,reads,strand='',{},{},[],{},''
	#names={read_name:[subs]} --> {'SRR11454613.2344642': ['TT', 'TT'], 'SRR11454613.525978': ['TC', 'TC'], 'SRR11454613.1449697': ['TT']}
	for pileupread in pileupcolumn.pileups:
		r1=pileupread.alignment.is_read1
		r2=pileupread.alignment.is_read2
		if r1: rt=1
		if r2: rt=2
		t='1'
		#if pileupread.alignment.is_reverse: st='-'
		if getstrand:
			#usa le info del mapping se strand oriented 
			if pileupread.alignment.is_read1:
				if unchange1:
					if pileupread.alignment.is_reverse: t='0'
					else: t='1'
				else:
					if pileupread.alignment.is_reverse: t='1'
					else: t='0'
			elif pileupread.alignment.is_read2:
				if unchange2:
					if pileupread.alignment.is_reverse: t='0'
					else: t='1'
				else:
					if pileupread.alignment.is_reverse: t='1'
					else: t='0'
			else: # for single ends
				if unchange1:
					if pileupread.alignment.is_reverse: t='0'
					else: t='1'
				else:
					if pileupread.alignment.is_reverse: t='1'
					else: t='0'					
		if br.has_key(pileupread.alignment.query_name): continue
		if badpos.has_key((pileupread.alignment.query_name,rt)): continue
		flag=pileupread.alignment.flag
		if pileupread.is_del: continue
		if pileupread.alignment.is_qcfail: continue
		if pileupread.alignment.is_supplementary: continue
		if pileupread.alignment.is_secondary: continue
		#if pileupread.alignment.has_tag('SA'): continue
		if pileupread.alignment.has_tag('NH'):
			if pileupread.alignment.get_tag('NH') > 1: continue
		if pileupread.alignment.is_duplicate: continue
		if pe:
			if not pileupread.alignment.is_paired: continue
		if pe:
			if not pileupread.alignment.is_proper_pair: continue
		if pileupread.alignment.is_duplicate: flag=flag-1024
		if pileupread.alignment.is_secondary: flag=flag-256
		if pe:
			if flag in [67,131,115,179]: continue
		if pileupread.alignment.mapping_quality < MAPQ: continue
		#if pileupread.alignment.has_tag('SA'): continue
		s,q,rname=pileupread.alignment.query_sequence[pileupread.query_position].upper(),pileupread.alignment.query_qualities[pileupread.query_position],pileupread.alignment.query_name
		#print ('\tbase in read %s = %s' %(pileupread.alignment.query_name,pileupread.alignment.query_sequence[pileupread.query_position]))
		startaln = pileupread.alignment.query_alignment_start+1
		endaln = pileupread.alignment.query_alignment_end
		# remove up and down positions
		rup,rdw=getUpDw(len(pileupread.alignment.query_sequence))
		nstartaln = startaln + rup
		nendaln = endaln - rdw
		#print s,q,rname,pileupread.query_position+1,startaln,endaln,nstartaln,nendaln
		if not nstartaln <= pileupread.query_position+1 <= nendaln : continue
		if q < MQUAL: continue
		#if ref==s: seq+=1
		fwRv=1 # forward
		if pileupread.alignment.is_reverse: fwRv=2 # reverse
		sub=ref+s
		#print sub
		if ref!=s:
			#print ref+s
			aln=pileupread.alignment.get_aligned_pairs(with_seq=True)
			if pileupread.query_position in [x[0] for x in aln[:5]]: continue
			if pileupread.query_position in [x[0] for x in aln[-5:]]: continue
			up=pileupread.query_position-5
			dw=pileupread.query_position+5
			#print up
			#print dw
			if up<0: up=0
			for k in aln[up:pileupread.query_position]:
				if k[0]==None: continue
				if k[1]==None: continue
				if k[2]==None: continue
			for k in aln[pileupread.query_position+1:dw+1]:
				if k[0]==None: continue
				if k[1]==None: continue
				if k[2]==None: continue			
		if not subs.has_key(sub): subs[sub]=[0,0,0,0] #n. subs,overalp,forward,reverse
		#print sub
		if names.has_key(rname):
			names[rname][0].append(sub)
			names[rname][1].append(fwRv)
		else: 
			names[rname]=[[sub],[fwRv]]
		if sub!=ref+ref:
			if reads.has_key(rname):
				reads[rname].append((rt,pileupread.alignment.query_sequence))
			else:
				reads[rname]=[(rt,pileupread.alignment.query_sequence)]
		strand+=t
		seq+=s
		#rnames.append((rname))
	xx=0
	#print names
	#print names #'SRR11454613.1080669': [['TG'], [1]]
	for i in names:
		if len(names[i][0])>1: #overlap
			if names[i][0][0]==names[i][0][1]: # sub confermata su overlap
				subs[names[i][0][0]][0]+=1
				subs[names[i][0][0]][1]+=1
		else:
			if names[i][0][0]!=ref+ref:
				reads1mod.append(reads[i])
				if names[i][1][0]==1: subs[names[i][0][0]][2]+=1
				else: subs[names[i][0][0]][3]+=1
				xx+=1
			subs[names[i][0][0]][0]+=1
	corrstrand='2'
	if getstrand:
		corrstrand=infstr(strand)
		errstrand=[]
		for k in range(len(seq)):
			#print seq[k],strand[k]
			if strand[k]!=corrstrand: errstrand.append(ref+seq[k])
		for k in errstrand:
			subs[k][0]-=1
	#print pos,subs,strand,len(strand),corrstrand
	#print d[pos]
	if len(subs)==1 and subs.keys()[0]==ref+ref: continue
	if len(subs)==0: continue
	#print [(subs[x][0],x) for x in subs if x!=ref+ref]
	#print subs
	msub=max([(subs[x][0],x) for x in subs if x!=ref+ref])[1]
	#obsSub = d[pos][7].split()[0]
	if subs[msub][0] < 4: continue
	if not getstrand:
		if subs[msub][2]==0 or subs[msub][3] ==0:
			if subs[msub][1]==0: continue
		tot=subs[msub][2]+subs[msub][3]
		v=tot*0.5
		#print tot,v
		try: rat,pval = stats.fisher_exact([[v, tot-v], [subs[msub][2], subs[msub][3]]])
		except: pval=1.0
		if pval < 0.05: continue
		#print pos, subs,tot,v,pval
		#if subs[msub][2]==0 or subs[msub][3] ==0: continue
	if subs[msub][0] >= 4: 
		#print subs, msub
		cov=sum([subs[x][0] for x in subs])
		if cov<MINCOV: continue
		a,c,g,t=0,0,0,0
		if subs.has_key(ref+'A'): a=subs[ref+'A'][0]
		if subs.has_key(ref+'C'): c=subs[ref+'C'][0]
		if subs.has_key(ref+'G'): g=subs[ref+'G'][0]
		if subs.has_key(ref+'T'): t=subs[ref+'T'][0]
		#print a,c,g,t, cov
		try: den=subs[msub[0]+ref][0]
		except: den=0
		try: freq=float(subs[msub][0])/(subs[msub][0]+den)
		except: freq=0.0
		line=d[pos]
		#print line
		if getstrand:
			if corrstrand=='2':
				line[2]=comp(ref)
				line[3]=corrstrand
				line[4]=str(cov)
				line[6]=str([t,g,c,a])
				line[7]=comp(msub)
				line[8]=str(freq)
			elif corrstrand=='1':
				line[2]=ref
				line[3]=corrstrand
				line[4]=str(cov)
				line[6]=str([a,c,g,t])
				line[7]=msub
				line[8]=str(freq)
		else:
			line[4]=str(cov)
			line[6]=str([a,c,g,t])
			line[7]=msub
			line[8]=str(freq)
		#print line
		o.write('\t'.join(line)+'\n')
o.close()

#	print reads1mod
	#print rnames, len(rnames)