"""Binary64 transfer envelope retaining the first nonzero feedback orbit.

Additional coordinates dominate a normalized distribution of Bx for one
known input shape. They persist through empty lazy steps; an active step
returns them to the old coarse representation. Not an outward certificate.
"""
from math import comb
import numpy as np
from scipy.sparse import csr_matrix


def memory_data(data,windows=8):
    spectrum,allowed,moments,atoms,cancel=data
    shapes=sorted(allowed)
    choices=np.array([windows*len(allowed[s]) for s in shapes],dtype=np.int64)
    nonzero=np.array([choices[i]-atoms[s][0] for i,s in enumerate(shapes)],dtype=np.int64)
    rows=[];columns=[];counts=[]
    for i,s in enumerate(shapes):
        for q,c in atoms[s].items():
            if q:
                rows.append(i);columns.append(q);counts.append(c)
    frequencies=csr_matrix((np.array(counts,dtype=np.int64),(rows,columns)),shape=(len(shapes),1<<19))
    collisions=(frequencies@frequencies.T).toarray()
    assert np.array_equal(collisions,collisions.T)
    for i,s in enumerate(shapes):
        assert collisions[i,i]==sum(c*c for q,c in atoms[s].items() if q)
    minimum=np.array([min(cancel[s]) for s in shapes])
    # [target shape, source memory]. The source draw is conditioned nonzero.
    pair=collisions.T.astype(float)/(choices[:,None]*nonzero[None,:])
    for target,s in enumerate(shapes):
        assert np.max(pair[target])<=max(c for q,c in atoms[s].items() if q)/choices[target]+1e-15
    return shapes,choices,nonzero,minimum,pair


def regions(data,extra,tilt,steps=256,check=False):
    spectrum,allowed,moments,atoms,cancel=data
    shapes,choices,nonzero,minimum,pair=extra
    levels=sorted(spectrum)
    size=len(shapes);m=(1<<19)-1
    powers=np.exp(-tilt*np.arange(145))
    fresh=np.array([0,0]+[spectrum[v]/(2*m) for v in levels])
    empty=np.zeros((7,7));empty[0,0]=1
    empty[1,1]=powers[48]/2;empty[1]+=powers[48]*fresh
    for i,v in enumerate(levels):
        empty[i+2,i+2]=powers[v]/2;empty[i+2]+=powers[v]*fresh
    mem_diagonal=powers[minimum]/2
    mem_empty=np.zeros((size,7))
    for i,s in enumerate(shapes):
        exact_moment=sum(sum(counts.values())*powers[v] for v,counts in cancel[s].items())/nonzero[i]
        assert exact_moment<=powers[minimum[i]]+1e-15
        mem_empty[i]=exact_moment*fresh
    active=np.zeros((size,7,7))
    mem_active=np.zeros((size,size,7))
    activation=np.zeros(size)
    for a,s in enumerate(shapes):
        weight=sum(s);q=atoms[s][0]/choices[a]
        active[a,0,0]=q*powers[weight]
        activation[a]=(1-q)*powers[weight]
        moment=powers[48-weight]
        maxatom=max(count for syndrome,count in atoms[s].items() if syndrome)
        least=min(w for row in cancel[s].values() for w in row)
        cancellation=min(moment,maxatom/choices[a]*powers[least])
        for i in range(1,7):
            if i>=2:
                v=levels[i-2]
                moment=sum(count*powers[w] for w,count in moments[s][v].items())/(choices[a]*spectrum[v])
                cancellation=sum(count*powers[w] for w,count in cancel[s][v].items())/(choices[a]*spectrum[v])
            active[a,i,0]=cancellation/2+moment/(2*m)
            active[a,i,1]=moment/2
            active[a,i]+=moment*fresh
        # Triangle inequality bounds the weighted mass of each memory orbit.
        # Its lazy return-to-zero coefficient uses the exact pair collision.
        f=powers[minimum-weight]
        c=powers[np.maximum(minimum-weight,least)]*pair[a]
        mem_active[a,:,0]=c/2+f/(2*m)
        mem_active[a,:,1]=f/2
        mem_active[a]+=f[:,None]*fresh[None,:]
    core_power=np.eye(7)
    memory_power=np.zeros_like(mem_empty)
    diagonal_power=np.ones(size)
    core=np.zeros_like(active)
    memory=np.zeros_like(mem_active)
    born=np.zeros(size)
    for _ in range(steps):
        core=core@empty+core_power@active
        core[:,0,:]+=born[:,None]*mem_empty
        memory=memory@empty+np.einsum('si,aij->asj',memory_power,active,optimize=True)+diagonal_power[None,:,None]*mem_active
        born=born*mem_diagonal+activation
        memory_power=memory_power@empty+diagonal_power[:,None]*mem_empty
        diagonal_power*=mem_diagonal
        core_power=core_power@empty
    if check:
        dense_empty=np.zeros((7+size,7+size))
        dense_empty[:7,:7]=empty
        dense_empty[7:,:7]=mem_empty
        dense_empty[7:,7:]=np.diag(mem_diagonal)
        expected=np.linalg.matrix_power(dense_empty,steps)
        assert np.allclose(expected[:7,:7],core_power,rtol=1e-12,atol=1e-18)
        assert np.allclose(expected[7:,:7],memory_power,rtol=1e-12,atol=1e-18)
        assert np.allclose(expected[7:,7:],np.diag(diagonal_power),rtol=1e-12,atol=1e-18)
        for a in sorted({0,min(10,size-1),min(30,size-1),min(50,size-1),size-1}):
            dense_active=np.zeros_like(dense_empty)
            dense_active[:7,:7]=active[a]
            dense_active[7:,:7]=mem_active[a]
            dense_active[0,7+a]=activation[a]
            direct=sum((np.linalg.matrix_power(dense_empty,j)@dense_active@np.linalg.matrix_power(dense_empty,steps-1-j)
                        for j in range(steps)),np.zeros_like(dense_empty))
            assembled=np.zeros_like(direct)
            assembled[:7,:7]=core[a];assembled[7:,:7]=memory[a];assembled[0,7+a]=born[a]
            assert np.allclose(direct,assembled,rtol=1e-12,atol=1e-18)
    return (shapes,core_power,memory_power,diagonal_power,core/steps,memory/steps,born/steps)


def self_test(data,extra):
    for tilt in (0,.0005,.002):
        for steps in (1,2,5):
            transfers=regions(data,extra,tilt,steps=steps,check=True)
            shapes,empty,mem_empty,diagonal,core,memory,born=transfers
            size=7+len(shapes)
            current=np.arange(1,1+3*size,dtype=float).reshape(3,size)/(3*size)
            actual=actions(current,transfers,.75)
            for b in range(1,len(shapes[0])+1):
                terms=[]
                for a,s in enumerate(shapes):
                    if sum(w!=0 for w in s)!=b:
                        continue
                    matrix=np.zeros((size,size))
                    matrix[:7,:7]=core[a];matrix[7:,:7]=memory[a];matrix[0,7+a]=born[a]
                    terms.append((current@matrix.T)*.75**sum(s))
                assert np.allclose(actual[b],np.maximum.reduce(terms),rtol=1e-12,atol=1e-18)
    print('Memory-region recurrence and backward actions match dense direct calculations',flush=True)


def actions(current,transfers,penalty):
    shapes,empty,empty_memory,diagonal,core,memory,born=transfers
    # Each active shape has only seven coarse output coordinates and its
    # own memory coordinate, reached only from zero. Avoid dense 76x76 maps.
    result=[np.concatenate((current[:,:7]@empty.T,
                            current[:,:7]@empty_memory.T+current[:,7:]*diagonal),axis=1)]
    scale=np.array([penalty**sum(s) for s in shapes])
    for b in range(1,len(shapes[0])+1):
        ids=[i for i,s in enumerate(shapes) if sum(w!=0 for w in s)==b]
        combined=np.concatenate((core[ids],memory[ids]),axis=1)
        values=(current[:,:7]@combined.reshape(-1,7).T).reshape(len(current),len(ids),7+len(shapes))
        values[:,:,0]+=current[:,7+np.array(ids)]*born[ids]
        values*=scale[ids][None,:,None]
        result.append(np.max(values,axis=1))
    return result
