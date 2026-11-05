"""Integer directed-cut lower bound for fixed target rows in a time window.

Drops incoming boundary terms, keeps half-full total activity as a
Lagrange constraint. Feasible flow is verified before returning a bound.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import maximum_flow,breadth_first_order


def cut_bound(target,partners,mu,scale=16,symmetric=False):
    L,n=target.shape;N=L*n
    assert N%2==0 and mu>=0 and scale>0 and partners.shape==target.shape
    assert not symmetric or scale%2==0
    selected=np.flatnonzero(target[1:].reshape(-1))+n
    aux=N+np.arange(len(selected),dtype=np.int32)
    source=N+len(aux);sink=source+1;vertices=sink+1
    finite=mu*N+scale*(int(target.sum())+(L-1)*n)
    infinity=finite+1
    assert infinity<2**31 and vertices<2**31
    rows=[];cols=[];caps=[]

    def add(a,b,c):
        a,b,c=np.broadcast_arrays(a,b,c)
        rows.append(np.asarray(a,dtype=np.int32).reshape(-1))
        cols.append(np.asarray(b,dtype=np.int32).reshape(-1))
        caps.append(np.asarray(c,dtype=np.int32).reshape(-1))

    add(source,np.arange(N,dtype=np.int32),mu)
    add(np.flatnonzero(target[0]),sink,scale)
    if len(selected):
        times=selected//n;positions=selected%n
        add(selected,aux,infinity)
        add(selected-n,aux,infinity)
        add((times-1)*n+partners[times,positions],aux,infinity)
        add(aux,sink,scale)
    previous=np.arange(N-n,dtype=np.int32);following=np.arange(n,N,dtype=np.int32)
    if symmetric:
        # Exactly the same death count when consecutive layers have equal weight.
        add(previous,following,scale//2);add(following,previous,scale//2)
    else:add(previous,following,scale)
    graph=coo_matrix((np.concatenate(caps),(np.concatenate(rows),np.concatenate(cols))),
                     shape=(vertices,vertices),dtype=np.int64).tocsr()
    graph.eliminate_zeros()
    result=maximum_flow(graph,source,sink,method='dinic')
    flow=result.flow
    assert (flow+flow.T).nnz==0,'flow must be antisymmetric'
    residual=graph-flow
    assert np.all(residual.data>=0),'capacity violation'
    balance=np.asarray(flow.sum(axis=1)).ravel()
    assert balance[source]==result.flow_value and balance[sink]==-result.flow_value
    balance[source]=balance[sink]=0
    assert not np.any(balance),'flow conservation violation'
    residual.eliminate_zeros()
    reachable=breadth_first_order(residual,source,directed=True,return_predecessors=False)
    mask=np.zeros(vertices,dtype=bool);mask[reachable]=True
    assert not mask[sink]
    assert int(graph[mask][:,~mask].sum())==result.flow_value,'cut/flow mismatch'
    state_counts=mask[:N].reshape(L,n).sum(axis=1)
    return dict(lower_scaled=int(result.flow_value)-mu*(N//2),scale=scale,
                flow=int(result.flow_value),mu=mu,states=N,
                minimizing_activity=int(state_counts.sum()),
                minimum_layer_activity=int(state_counts.min()),
                maximum_layer_activity=int(state_counts.max()),
                symmetric_deaths=symmetric,
                verified_integer_flow=True)
