import re,math,fractions
s=open('SpinCodes/Majorant/RefinedData.lean',encoding='utf-8').read().split('def refined')[0]
pairs=[tuple(float(fractions.Fraction(t.replace(' ',''))) for t in p) for p in re.findall(r'\(\(([-\d /]+)\), \(([-\d /]+)\)\)',s)]
for x in [0.104,.105,.11,.12,.13,.15,.175,.2,.25,.3,.35,.4,.45,.5]:
 vals=[a*x+c for a,c in pairs]; i=min(range(len(vals)),key=vals.__getitem__); h=-x*math.log(x)-(1-x)*math.log(1-x)
 print(x,i, vals[i]-h+math.log(2)/2)

