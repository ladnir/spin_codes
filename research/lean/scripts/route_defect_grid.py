exec(open('scripts/route_defect_probe.py',encoding='utf-8').read().split('for x in')[0])
xs=[.104,.15,.2,.25,.3,.35,.4,.45,.5]
for a,b in zip(xs,xs[1:]):
 def d(x,i):return pairs[i][0]*x+pairs[i][1]+x*math.log(x)+(1-x)*math.log(1-x)+math.log(2)/2
 i=min(range(len(pairs)), key=lambda i:max(d(a,i),d(b,i)))
 print(a,b,i,d(a,i),d(b,i))
