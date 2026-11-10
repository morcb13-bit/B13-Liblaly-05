"""第3章：8次元から4次元へ ── E₈からH₄、二つの600胞体
イコシアン（成分が a+bφ の四元数）を、成分を2倍した整数の組で持つ（zphi.py の Z）。
  ノルム N(q)=成分の二乗和。4N = X+Yφ（2倍の座標なので4倍）。
  共役 q*：φ を 1−φ に替える（残渣の側）。
  8次元の長さ² E(q) = X/2（N の「φ のつかない部分」の2倍）。
基準（走らせる前に決めた）
  H0 単位イコシアン120個とその φ 倍120個が張る格子は、整数の基底8本を持ち、長さ² E で測ると
     グラム行列は整数で、正定値（首座小行列式がすべて正）、偶（対角がすべて偶数）、行列式1
     ［知識：8次元でこの性質を持つ格子は E₈ だけ］
  H1 長さ² E が2の240本は、4次元の側で長さ² 1 の120本と φ² の120本に分かれ、
     どちらも600胞体（各頂点の最近接12・辺720）。二つは大きさがちょうど φ 倍
"""
import sys, itertools
from zphi import Z
R2=int(sys.argv[1]) if len(sys.argv)>1 else 6
P=Z(0,1); ONE=Z(1); ZERO=Z(0)
cj=lambda z: Z(z.a+z.b,-z.b)
def q4(q): # 4N
    s=Z()
    for c in q: s=s+c*c
    return s
def E(q):
    X=q4(q).a; assert X%2==0; return X//2
def add(p,q): return tuple(a+b for a,b in zip(p,q))
def sub(p,q): return tuple(a-b for a,b in zip(p,q))
def scl(p,s): return tuple(a*s for a in p)
def star(q): return tuple(cj(c) for c in q)
def dot(p,q):
    s=Z()
    for a,b in zip(p,q): s=s+a*b
    return s
# 単位イコシアン120個（2倍の座標）
U=set()
for i in range(4):
    for s in (2,-2):
        v=[ZERO]*4; v[i]=Z(s); U.add(tuple(v))
for sg in itertools.product((1,-1),repeat=4): U.add(tuple(Z(x) for x in sg))
base=[ZERO,ONE,P,P-ONE]
def even(p):
    inv=sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j]); return inv%2==0
for perm in itertools.permutations(range(4)):
    if not even(perm): continue
    for sg in itertools.product((1,-1),repeat=3):
        v=[ZERO]*4
        vals=[ZERO,ONE*sg[0],P*sg[1],(P-ONE)*sg[2]]
        for k in range(4): v[perm[k]]=vals[k]
        U.add(tuple(v))
U=sorted(U,key=repr)
print("単位イコシアン",len(U),"／4N が 4 でないもの",sum(1 for u in U if q4(u)!=Z(4)))
R=list(U)+[scl(u,P) for u in U]
# ---- H0 格子の基底（整数の掃き出し）とグラム行列 ----
def ints(q): return [x for c in q for x in (c.a,c.b)]
def unints(v): return tuple(Z(v[2*k],v[2*k+1]) for k in range(4))
rows=[ints(r) for r in R]
def hnf(rows):
    rows=[r[:] for r in rows]; basis=[]; col=0; n=len(rows[0])
    for col in range(n):
        piv=[r for r in rows if r[col]!=0]
        if not piv: continue
        while True:
            piv=[r for r in rows if r[col]!=0]
            if len(piv)<=1: break
            piv.sort(key=lambda r:abs(r[col])); p=piv[0]
            for r in rows:
                if r is not p and r[col]!=0:
                    m=r[col]//p[col]
                    for k in range(n): r[k]-=m*p[k]
        piv=[r for r in rows if r[col]!=0]
        if piv:
            p=piv[0]; basis.append(p[:]); rows=[r for r in rows if r is not p]
    return basis
B=[unints(b) for b in hnf(rows)]
def Bf(p,q): return E(add(p,q))-E(p)-E(q)
G=[[Bf(p,q) for q in B] for p in B]
assert all(x%2==0 for r in G for x in r)
G=[[x//2 for x in r] for r in G]      # B(p,p)=E(p) になる対称な形
def minors(M):
    """首座小行列式を整数だけで（Bareiss の割り切れる消去）。途中で 0 が出たら打ち切り"""
    M=[r[:] for r in M]; n=len(M); prev=1; out=[]
    for k in range(n):
        if M[k][k]==0: out.append(0); break
        out.append(M[k][k])
        for i in range(k+1,n):
            for j in range(k+1,n):
                num=M[i][j]*M[k][k]-M[i][k]*M[k][j]
                assert num%prev==0; M[i][j]=num//prev
        prev=M[k][k]
    return out
sym=all(G[i][j]==G[j][i] for i in range(8) for j in range(8))
mn=minors(G)
print("基底",len(B),"本／グラム行列は整数・対称",sym)
print("首座小行列式",mn,"（すべて正なら正定値、最後が行列式）")
print("対角",[G[i][i] for i in range(len(B))])
h0= len(B)==8 and sym and len(mn)==8 and all(x>0 for x in mn) and mn[-1]==1 and all(G[i][i]%2==0 for i in range(8))
print("H0 整数・正定値・偶・行列式1（E₈）","OK" if h0 else "NG")
# ---- H1 根240本と二つの600胞体 ----
print("長さ² E が 2 の本数（構成した240本のうち）",sum(1 for r in R if E(r)==2))
def cell(S):
    d=[q4(sub(a,b)) for a,b in itertools.combinations(S,2)]
    m=min(d); e=sum(1 for x in d if x==m)
    deg={}
    for a,b in itertools.combinations(S,2):
        if q4(sub(a,b))==m: deg[a]=deg.get(a,0)+1; deg[b]=deg.get(b,0)+1
    return m,e,set(deg.values())
S1=list(U); S2=[scl(u,P) for u in U]
m1,e1,g1=cell(S1); m2,e2,g2=cell(S2)
print("4次元の側 半径²1 の組：点",len(S1),"辺²×4",m1,"辺",e1,"次数",g1)
print("4次元の側 半径²φ² の組：点",len(S2),"辺²×4",m2,"辺",e2,"次数",g2)
h1= e1==720 and e2==720 and g1=={12} and g2=={12} and m2==m1*P*P
print("H1 二つの600胞体で大きさは φ 倍","OK" if h1 else "NG")
print("第3章の検定 総合","OK" if (h0 and h1) else "NG")
