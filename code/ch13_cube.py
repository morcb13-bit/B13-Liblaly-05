"""第13章：1次元から5次元へ ── 回転と時間で描く
n 次元の立方体の頂点を ±1 の列で持つ。すべて整数の数え上げ。
基準（走らせる前に決めた）
  C1 n＝1〜5 で、k 次元の面を全部数えると、一つ下の次元から「写して、つなぐ」足し算
     f_k(n) ＝ 2·f_k(n−1) ＋ f_(k−1)(n−1) とすべて一致する。面の数を符号を交互に足すと、どの n でも 1。
  C2 頂点を長い対角線の向きの段（座標の和）で分けると、段ごとの頂点の数は 1, n, …（二項係数）。
     4次元のまん中の段の頂点6個は正八面体（距離² 8 の組12・距離² 16 の組3、どの頂点も4本）、
     3次元の段の頂点3個は正三角形、5次元の段の頂点5個は4次元の正単体（距離² がすべて 8）。
  C3 立方体を自分に重ねる回転（±1 の置換行列で行列式 ＋1）のうち、動かない向きが一つもないものは、
     n＝2, 4 にあり、n＝1, 3, 5 には一つもない。
     ［知識］4次元で動かない向きのない回転は、直交する二つの面を同時に回す。そのうち −1 倍（向かいへの写し）以外のものがある。
  C4 頂点を結んでできる正方形・立方体のうち、面に乗らない（辺の長さ²が4より大きい）ものは、
     正方形は n＝4 で初めて現れ（12個）、n＝5 で 120 個。立方体は n≤5 で 0 個。
  負の対照 面に乗る正方形・立方体の数は、C(n,k)·2^(n−k) と一致する。
"""
import itertools, collections
from math import comb
def verts(n): return [tuple(v) for v in itertools.product((1,-1),repeat=n)]
def d2(a,b): return sum((x-y)*(x-y) for x,y in zip(a,b))
# C1
def faces(n,k):
    S=set()
    for free in itertools.combinations(range(n),k):
        fixed=[i for i in range(n) if i not in free]
        for sg in itertools.product((1,-1),repeat=len(fixed)):
            S.add((free,tuple(zip(fixed,sg))))
    return len(S)
F={n:[faces(n,k) for k in range(n+1)] for n in range(0,6)}
rule=True
for n in range(1,6):
    for k in range(n+1):
        prev=F[n-1][k] if k<=n-1 else 0
        prevk=F[n-1][k-1] if k>=1 else 0
        if F[n][k]!=2*prev+prevk: rule=False
euler={n:sum((-1)**k*F[n][k] for k in range(n+1)) for n in range(1,6)}
c1= rule and all(v==1 for v in euler.values())
for n in range(1,6): print(f"   {n}次元の立方体 面の数（頂点,辺,…）",F[n])
print("C1 写してつなぐ足し算と一致",rule,"／交互の和",euler,"→","OK" if c1 else "NG")
# C2
def level(n,j): return [v for v in verts(n) if v.count(-1)==j]
lv={n:[len(level(n,j)) for j in range(n+1)] for n in range(1,6)}
binom= all(lv[n]==[comb(n,j) for j in range(n+1)] for n in range(1,6))
def dist_kinds(P): return dict(collections.Counter(d2(a,b) for a,b in itertools.combinations(P,2)))
mid4=level(4,2); k4=dist_kinds(mid4)
deg4=collections.Counter(sum(1 for b in mid4 if d2(a,b)==8) for a in mid4)
tri3=dist_kinds(level(3,1)); sim5=dist_kinds(level(5,1))
c2= binom and k4=={8:12,16:3} and dict(deg4)=={4:6} and tri3=={8:3} and sim5=={8:10}
print("C2 段ごとの頂点",lv,"／4次元のまん中",k4,"次数",dict(deg4),"／3次元の段",tri3,"／5次元の段",sim5,"→","OK" if c2 else "NG")
# C3
def det(M):
    M=[r[:] for r in M]; n=len(M); prev=1; sgn=1
    for k in range(n-1):
        if M[k][k]==0:
            sw=next((i for i in range(k+1,n) if M[i][k]!=0),None)
            if sw is None: return 0
            M[k],M[sw]=M[sw],M[k]; sgn=-sgn
        for i in range(k+1,n):
            for j in range(k+1,n):
                M[i][j]=(M[i][j]*M[k][k]-M[i][k]*M[k][j])//prev
        prev=M[k][k]
    return sgn*M[n-1][n-1]
def rotations(n):
    out=[]
    for perm in itertools.permutations(range(n)):
        for sg in itertools.product((1,-1),repeat=n):
            M=[[0]*n for _ in range(n)]
            for i in range(n): M[i][perm[i]]=sg[i]
            if det(M)==1: out.append(M)
    return out
def nofix(M):
    n=len(M); A=[[M[i][j]-(1 if i==j else 0) for j in range(n)] for i in range(n)]
    return det(A)!=0
cnt={}; free={}
for n in range(1,6):
    R=rotations(n); cnt[n]=len(R); free[n]=sum(1 for M in R if nofix(M))
# 4次元で動かない向きがなく、−1 倍でも M²＝1 でもないもの
R4=rotations(4)
def mul(A,B): n=len(A); return [[sum(A[i][k]*B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
I4=[[1 if i==j else 0 for j in range(4)] for i in range(4)]
neg=[[-x for x in r] for r in I4]
double=sum(1 for M in R4 if nofix(M) and M!=neg and mul(M,M)!=I4)
c3= all((free[n]>0)==(n in (2,4)) for n in range(1,6)) and double>0
print("C3 回転の数",cnt,"／動かない向きがない回転",free,"／4次元で二つの面を同時に回すもの",double,"→","OK" if c3 else "NG")
# C4（k1_inner_cube.py と同じ数え方）
def cubes(n,k):
    V=verts(n); S=set(V); v0=(1,)*n
    D=[tuple(w[i]-v0[i] for i in range(n)) for w in V if w!=v0]
    found={}
    def grow(ch,st):
        if len(ch)==k:
            pts=[]
            for m in range(1<<k):
                p=list(v0)
                for j in range(k):
                    if m>>j&1: p=[p[i]+ch[j][i] for i in range(n)]
                pts.append(tuple(p))
            if all(p in S for p in pts):
                L=sum(x*x for x in ch[0]); found.setdefault(L,set()).add(frozenset(pts))
            return
        for i in range(st,len(D)):
            a=D[i]
            if ch and sum(x*x for x in a)!=sum(x*x for x in ch[0]): continue
            if all(sum(x*y for x,y in zip(a,c))==0 for c in ch): grow(ch+[a],i+1)
    grow([],0)
    return {L:len(s)*(1<<n)//(1<<k) for L,s in found.items()}
inner={}; facec={}; ctrl=True
for k in (2,3):
    for n in range(k,6):
        c=cubes(n,k); facec[(n,k)]=c.get(4,0); inner[(n,k)]=sum(v for L,v in c.items() if L>4)
        if facec[(n,k)]!=comb(n,k)*2**(n-k): ctrl=False
sq={n:inner[(n,2)] for n in range(2,6)}; cu={n:inner[(n,3)] for n in range(3,6)}
c4= sq=={2:0,3:0,4:12,5:120} and all(v==0 for v in cu.values())
print("C4 面に乗らない正方形",sq,"／面に乗らない立方体",cu,"→","OK" if c4 else "NG")
print("負の対照 面に乗る正方形・立方体",{f"{n}次元{k}":facec[(n,k)] for (n,k) in facec},"＝C(n,k)·2^(n−k)",ctrl,"→","OK" if ctrl else "NG")
print("第13章の検定 総合","OK" if (c1 and c2 and c3 and c4 and ctrl) else "NG")
