import sys
f=open(sys.argv[1])
n,x,xx=0.0,0.0,0.0
for i in f:
	l=(i.strip()).split('\t')
	v=int(l[2])
	if v>0:
		x+=1
		xx+=v
	n+=1
f.close()

try: dp=xx/n
except: dp=0.0
try: cov=x/n
except: cov=0.0

print 'Depth:',dp
print 'Coverage:',cov

