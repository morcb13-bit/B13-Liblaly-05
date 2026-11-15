"""第14章：6次元と7次元 ── 立方体の中の立方体と、立方体に嵌る正単体
n 次元の立方体の頂点を ±1 の列で持つ。すべて整数の数え上げ。
基準（走らせる前に決めた）
  F1 頂点を結んでできる立方体のうち、面に乗らない（辺の長さ²が4より大きい）ものは、n≤5 で 0、n＝6 で初めて現れる。
     n＝6, 7 の数は記録として出す。
  F2 正単体（n＋1 個の頂点で、どの二つの距離²も等しい）の頂点を、すべて立方体の頂点に置けるのは、
     n≤8 の中では n＝1, 3, 7 だけ。n＝11 でも置ける（一つ見つかればよい）。
     ［知識］置けるのは n＋1 が 4 の倍数のとき（と n＝1）で、アダマール行列の有無と同じ問い。n＝9, 10 は調べない。
  F3 置けるとき、正単体の辺の長さ²は、n＝3 で 8（正方形の対角線）、n＝7 で 16、n＝11 で 24
     （列の半分の (n＋1)/2 か所が違う）。立方体の中に入る正単体の数は、n＝3 で 2、n＝7 は記録。
  負の対照 面に乗る立方体の数は C(n,3)·2^(n−3) と一致する。
"""
import itertools, sys
from math import comb
sys.setrecursionlimit(10000)
def verts(n): return [tuple(v) for v in itertools.product((1,-1),repeat=n)]
def ham(a,b): return sum(1 for x,y in zip(a,b) if x!=y)
# F1（k1_inner_cube.py と同じ数え方）
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
inner={}; face={}; ctrl=True; bylen={}
for n in range(3,8):
    c=cubes(n,3); face[n]=c.get(4,0); inner[n]=sum(v for L,v in c.items() if L>4)
    bylen[n]={L:v for L,v in sorted(c.items()) if L>4}
    if face[n]!=comb(n,3)*2**(n-3): ctrl=False
f1= all(inner[n]==0 for n in (3,4,5)) and inner[6]>0
print("F1 面に乗らない立方体",inner,"／辺の長さ²ごと（n=6,7）",bylen[6],bylen[7],"→","OK" if f1 else "NG")
# F2・F3：v0＝(1,…,1) を含む正単体を数える（立方体はどの頂点も同じ見え方なので、v0 を固定してよい）
def simplices(n,d,count=True):
    V=verts(n); v0=(1,)*n
    C=[v for v in V if ham(v,v0)==d]
    idx={v:i for i,v in enumerate(C)}
    nb=[set(j for j,w in enumerate(C) if j!=i and ham(C[i],w)==d) for i in range(len(C))]
    need=n   # v0 以外に n 個
    total=0; ex=None
    def rec(chosen,cand):
        nonlocal total,ex
        if len(chosen)==need:
            total+=1
            if ex is None: ex=[v0]+[C[i] for i in chosen]
            return not count
        if len(chosen)+len(cand)<need: return False
        for i in sorted(cand):
            if rec(chosen+[i],{j for j in cand if j>i and j in nb[i]}): return True
        return False
    rec([],set(range(len(C))))
    return total,ex
fit={}; edge={}; num={}
for n in range(1,9):
    found=[]
    for d in range(1,n+1):
        t,ex=simplices(n,d,count=(n<=7))
        if t: found.append((d,t,ex))
    fit[n]=bool(found)
    if found:
        d,t,ex=found[0]; edge[n]=4*d
        if n<=7: num[n]=t*2**n//(n+1)      # v0 を含むもの × 頂点の数 ÷ 一つの正単体の頂点の数
t11,ex11=simplices(11,6,count=False); fit[11]=t11>0; edge[11]=24 if t11 else None
f2= [n for n in range(1,9) if fit[n]]==[1,3,7] and fit[11]
f3= edge.get(3)==8 and edge.get(7)==16 and edge.get(11)==24 and num.get(3)==2
print("F2 正単体を立方体の頂点に置ける次元（n≤8）",[n for n in range(1,9) if fit[n]],"／n=11",fit[11],"→","OK" if f2 else "NG")
print("F3 辺の長さ²",edge,"／立方体の中に入る正単体の数",num,"→","OK" if f3 else "NG")
print("   n=7 の一つ",["".join("+" if x>0 else "-" for x in v) for v in simplices(7,4,count=False)[1]])
print("   n=11 の一つ",["".join("+" if x>0 else "-" for x in v) for v in ex11])
print("負の対照 面に乗る立方体",face,"＝C(n,3)·2^(n−3)",ctrl,"→","OK" if ctrl else "NG")
print("第14章の検定 総合","OK" if (f1 and f2 and f3 and ctrl) else "NG")
