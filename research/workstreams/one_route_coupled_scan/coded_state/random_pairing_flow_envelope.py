"""Analytical quadratic envelopes for the Gaussian flow integrand.

Floating feasibility screen, not directed certification. Standard library,
bounded serial execution, no output files or large arrays.
"""
import math
import shutil
import time
from itertools import permutations,product
from pathlib import Path
from random_pairing_flow_quadratic import (
    TYPES,LAW,P,OUTPUT,LOSS,LOG_WEIGHT_NUMERATORS,dot,variational,moments,
    tables,probability_colors,endpoint_allowance,
)


def norm(x):return math.sqrt(dot(x,x))


def local_data(h,u,theta,eta,kappa,center=None):
    n=sum(h);_,proposal,_=variational(h,u,theta,eta,kappa)
    if center is None:center,_=moments(h,tables(proposal))
    assert len(center)==5 and abs(sum(center)-n)<1e-9
    beta=math.sqrt(2*kappa/n)
    basis=[[0.]*4 for _ in range(5)]
    for j in range(4):
        for i in range(j+1):basis[i][j]=1/math.sqrt((j+1)*(j+2))
        basis[j+1][j]=-(j+1)/math.sqrt((j+1)*(j+2))
    projected={key:[beta*sum(basis[c][i]*out[c] for c in range(5)) for i in range(4)] for key,out in OUTPUT.items()}
    base=[u[c]+2*kappa*center[c]/n for c in range(5)]
    entries={};core=0.;core_gradient=[0.]*4
    for a,b in TYPES:
        logits=[math.log(float(LAW[j]))+theta*LOSS[a,b,j]+eta*j+dot(base,OUTPUT[a,b,j]) for j in range(3)]
        offset=max(logits);weights=[math.exp(x-offset) for x in logits];mass=sum(weights);weights=[x/mass for x in weights]
        f=offset+math.log(mass)
        gradient=[sum(weights[j]*projected[a,b,j][i] for j in range(3)) for i in range(4)]
        covariance=[[sum(weights[j]*projected[a,b,j][i]*projected[a,b,j][k] for j in range(3))-gradient[i]*gradient[k]
                     for k in range(4)] for i in range(4)]
        diameter=max(norm([x-y for x,y in zip(projected[a,b,j],projected[a,b,k])]) for j in range(3) for k in range(3))
        numerator,denominator=probability_colors(h,(a,b));count=(n//2)*numerator/denominator
        centered=[([projected[a,b,j][i]-gradient[i] for i in range(4)],weights[j]) for j in range(3)]
        entries[a,b]=(f,gradient,covariance,diameter,count,centered)
        core+=count*f
        for i in range(4):core_gradient[i]+=count*gradient[i]
    core-=eta*n*float(P)+dot(u,h)+kappa*(dot(h,h)+dot(center,center))/n
    for i in range(4):core_gradient[i]-=beta*sum(basis[c][i]*center[c] for c in range(5))
    return entries,core,core_gradient


def switches(entries,active):
    results=[];seen=set()
    for a,b,x,c in product(active,repeat=4):
        coefficients={}
        for pair,sign in (((x,a),1),((x,b),-1),((c,a),-1),((c,b),1)):
            pair=tuple(sorted(pair));coefficients[pair]=coefficients.get(pair,0)+sign
        key=tuple(sorted((pair,coefficient) for pair,coefficient in coefficients.items() if coefficient))
        if not key or key in seen:continue
        seen.add(key)
        value=sum(coef*entries[pair][0] for pair,coef in key)
        gradient=[sum(coef*entries[pair][1][i] for pair,coef in key) for i in range(4)]
        hessian=[[sum(coef*entries[pair][2][i][j] for pair,coef in key) for j in range(4)] for i in range(4)]
        frobenius=math.sqrt(sum(x*x for row in hessian for x in row))
        results.append((abs(value),norm(gradient),frobenius,key))
    return results


def switch_envelope(entries,subset,radius):
    options=switches(entries,subset)
    if not options:return 0.,0.,0.
    curvature=[]
    for _,_,hessian,key in options:
        residual=sum(abs(coef)*math.expm1(entries[pair][3]*radius)*sum(entries[pair][2][i][i] for i in range(4)) for pair,coef in key)
        curvature.append(hessian+residual)
    return max(x[0] for x in options),max(x[1] for x in options),max(curvature)


def logdet_positive(matrix):
    a=[row.copy() for row in matrix];result=0.
    for j in range(len(a)):
        pivot=a[j][j]
        if pivot<=0:return None
        result+=math.log(pivot)
        for i in range(j+1,len(a)):
            for k in range(j+1,len(a)):a[i][k]-=a[i][j]*a[j][k]/pivot
    return result


def eigensystem(matrix):
    size=len(matrix);a=[row.copy() for row in matrix]
    vectors=[[float(i==j) for j in range(size)] for i in range(size)]
    for _ in range(300):
        value,p,q=max((abs(a[i][j]),i,j) for i in range(size) for j in range(i+1,size))
        if value<1e-18:break
        angle=.5*math.atan2(2*a[p][q],a[q][q]-a[p][p]);c=math.cos(angle);s=math.sin(angle)
        app,aqq,apq=a[p][p],a[q][q],a[p][q]
        a[p][p]=c*c*app-2*s*c*apq+s*s*aqq;a[q][q]=s*s*app+2*s*c*apq+c*c*aqq;a[p][q]=a[q][p]=0.
        for i in range(size):
            if i not in (p,q):
                oldp,oldq=a[i][p],a[i][q]
                a[i][p]=a[p][i]=c*oldp-s*oldq;a[i][q]=a[q][i]=s*oldp+c*oldq
            oldp,oldq=vectors[i][p],vectors[i][q]
            vectors[i][p]=c*oldp-s*oldq;vectors[i][q]=s*oldp+c*oldq
    return [a[i][i] for i in range(size)],vectors


def affine_ellipsoid(entries,subset):
    """Find, then numerically check M >= vv^T for every affine switch v."""
    options=switches(entries,subset)
    vectors=[]
    for _,_,_,key in options:
        vectors.append([sum(coef*entries[pair][0] for pair,coef in key)]
                       +[sum(coef*entries[pair][1][i] for pair,coef in key) for i in range(4)])
    if not vectors:return [[0.]*5 for _ in range(5)]
    epsilon=1e-14
    covariance=[[sum(v[i]*v[j] for v in vectors)/len(vectors)+epsilon*int(i==j) for j in range(5)] for i in range(5)]
    best=(float('inf'),None)
    for iteration in range(80):
        eigen,basis=eigensystem(covariance);root=[math.sqrt(max(x,epsilon/10)) for x in eigen];trace=sum(root)
        leverage=[sum(sum(basis[i][j]*v[i] for i in range(5))**2/root[j] for j in range(5)) for v in vectors]
        scale=max(1.,max(leverage)/trace)*(1+1e-9)
        objective=scale*trace*trace
        if objective<best[0]:
            matrix=[[scale*trace*sum(basis[i][k]*root[k]*basis[j][k] for k in range(5)) for j in range(5)] for i in range(5)]
            best=(objective,matrix)
        vertex=vectors[max(range(len(vectors)),key=lambda i:leverage[i])];step=2/(iteration+5)
        covariance=[[(1-step)*covariance[i][j]+step*(vertex[i]*vertex[j]+epsilon*int(i==j)) for j in range(5)] for i in range(5)]
    matrix=best[1]
    for v in vectors:
        assert logdet_positive([[matrix[i][j]-v[i]*v[j] for j in range(5)] for i in range(5)]) is not None
    return matrix


def solve(matrix,vector):
    size=len(vector);a=[row.copy()+[value] for row,value in zip(matrix,vector)]
    for j in range(size):
        pivot=max(range(j,size),key=lambda i:abs(a[i][j]));a[j],a[pivot]=a[pivot],a[j]
        factor=a[j][j];assert abs(factor)>1e-15;a[j]=[x/factor for x in a[j]]
        for i in range(size):
            if i!=j:
                factor=a[i][j];a[i]=[x-factor*y for x,y in zip(a[i],a[j])]
    return [row[-1] for row in a]


def core_curvature(entries,radius):
    # Jensen bounds the centered tilt denominator below by one. Bound each
    # numerator separately, retaining its rank-one covariance contribution.
    return [[sum(value[4]*weight*math.exp(norm(vector)*radius)*vector[i]*vector[j]
                 for value in entries.values() for vector,weight in value[5])
             for j in range(4)] for i in range(4)]


def ellipsoid_shell(h,entries,core,gradient,radius,ellipsoids,lower=0.,truncate=False):
    active=tuple(a for a,count in enumerate(h) if count);cache={};best=(float('inf'),None)
    curvature=core_curvature(entries,radius)
    for priority in permutations(active):
        for slack in (.02,.05,.1,.2,.4):
            remaining=sum(h)//2-1;budget=[[0.]*5 for _ in range(5)]
            for index,a in enumerate(priority):
                if remaining<=0:break
                subset=tuple(sorted(priority[index:]))
                if subset not in cache:cache[subset]=switch_envelope(entries,subset,radius)[2]/2
                steps=min(h[a],remaining);remaining-=steps;matrix=ellipsoids[subset]
                for i in range(5):
                    for j in range(5):budget[i][j]+=steps*(1+slack)*matrix[i][j]/8
                remainder=steps*(1+1/slack)*cache[subset]**2*radius**2/8
                for i in range(1,5):budget[i][i]+=remainder
            precision=[[float(i==j)-curvature[i][j]-2*budget[i+1][j+1] for j in range(4)] for i in range(4)]
            determinant=logdet_positive(precision)
            if determinant is None:continue
            linear=[gradient[i]+2*budget[0][i+1] for i in range(4)];mean=solve(precision,linear)
            value=core+budget[0][0]-.5*determinant+.5*dot(linear,mean)
            if lower or truncate:
                eigen,_=eigensystem(precision);least=min(eigen);largest=max(eigen)
                if least<=0:continue
                x=least*max(0.,lower-norm(mean))**2/2
                mass=(1+x)*math.exp(-x)
                if truncate:
                    y=largest*(radius+norm(mean))**2/2
                    mass-=(1+y)*math.exp(-y)
                assert mass>0
                value+=math.log(mass)
            if value<best[0]:best=(value,dict(priority=priority,slack=slack,constant=budget[0][0]))
    return best


def shell_envelope(h,entries,core,gradient,radius,lower=0.):
    active=tuple(a for a,count in enumerate(h) if count);cache={};best=(float('inf'),None)
    curvature=[[sum(value[4]*math.exp(value[3]*radius)*value[2][i][j] for value in entries.values())
                for j in range(4)] for i in range(4)]
    eigen_upper=max(sum(abs(x) for x in row) for row in curvature)
    for priority in permutations(active):
        remaining=sum(h)//2-1;c0=c1=c2=0.
        for index,a in enumerate(priority):
            if remaining<=0:break
            subset=tuple(sorted(priority[index:]))
            if subset not in cache:cache[subset]=switch_envelope(entries,subset,radius)
            delta,slope,hessian=cache[subset];steps=min(h[a],remaining);remaining-=steps
            c0+=steps*delta*delta/8
            c1+=steps*2*delta*slope/8
            c2+=steps*(slope*slope+delta*hessian+slope*hessian*radius+hessian*hessian*radius*radius/4)/8
        linear=c1+norm(gradient)
        for epsilon in (0.00001,0.00003,0.0001,0.0003,0.001,0.003,0.01,0.03,0.1):
            precision=[[float(i==j)*(1-2*c2-epsilon)-curvature[i][j] for j in range(4)] for i in range(4)]
            determinant=logdet_positive(precision)
            least=1-2*c2-epsilon-eigen_upper
            if determinant is None or least<=0:continue
            value=core+c0+linear*linear/(2*epsilon)-determinant/2
            if lower:
                x=least*lower*lower/2;value+=math.log1p(x)-x
            if value<best[0]:best=(value,dict(priority=priority,epsilon=epsilon,c0=c0,c1=c1,c2=c2))
    return best,cache


def main():
    assert shutil.disk_usage(Path.cwd()).free>=20*(1<<30),'disk guard'
    began=time.monotonic();h=(185,188,83,56,0)
    u=[x/10**9 for x in LOG_WEIGHT_NUMERATORS];theta=.1;eta=.2;kappa=.08
    entries,core,gradient=local_data(h,u,theta,eta,kappa)
    allowance=endpoint_allowance(u,theta,eta,kappa)
    print(dict(core_at_zero=core,core_gradient_norm=norm(gradient),allowance=allowance),flush=True)
    active=tuple(a for a,count in enumerate(h) if count)
    subsets={tuple(sorted(priority[index:])) for priority in permutations(active) for index in range(len(active))}
    ellipsoids={subset:affine_ellipsoid(entries,subset) for subset in subsets}
    for radius in (3.,4.,5.,6.):
        best,cache=shell_envelope(h,entries,core,gradient,radius)
        print(dict(radius=radius,central_log_integral_envelope=best[0],parameters=best[1],
                   full_color_range=cache[tuple(a for a,count in enumerate(h) if count)],
                   floating_screen_not_certificate=True),flush=True)
        improved=ellipsoid_shell(h,entries,core,gradient,radius,ellipsoids)
        print(dict(radius=radius,ellipsoid_log_envelope=improved[0],parameters=improved[1]),flush=True)
    boundaries=(0.,4.,4.5,5.,6.,8.,12.,16.)
    total=0.
    for lower,radius in zip(boundaries,boundaries[1:]):
        bound=ellipsoid_shell(h,entries,core,gradient,radius,ellipsoids,lower,True)
        print(dict(shell=(lower,radius),log_bound=bound[0]),flush=True)
        total+=math.exp(bound[0])
    print(dict(shell_sum_log=math.log(total),tail_beyond_16_still_separate=True),flush=True)
    assert time.monotonic()-began<25
    print(dict(seconds=time.monotonic()-began))


if __name__=='__main__':main()
