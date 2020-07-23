import sys
s={}
for i in 'ACGT':
	for j in 'ACGT':
		if i!=j: s[i+j]=0
n={}
x=0
for i in 'ACGT':
	n[i]=x
	x+=1
all=0
f=open(sys.argv[1])
for i in f:
	if i.startswith('Reg'): continue
	l=(i.strip()).split('\t')
	if l[7]=='-': continue
	sub=l[7].split()[0]
	s[sub]+=1
	all+=1
f.close()


h,cc,all_,vals=[],[],[],[]
for i in sorted(s.keys()):
	try: v=(s[i]/float(all))*100
	except: v=0.0
	#print i,s[i],all,v
	h.append(i)
	cc.append(str(s[i]))
	all_.append(str(all))
	vals.append(str(v))

print '\t'.join(h)
print '\t'.join(cc)
print '\t'.join(all_)
print '\t'.join(vals)
