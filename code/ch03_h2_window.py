"""第3章 別検定（保留）：H2 窓で切り出す。本には載せない。φ 倍の窓での再検定はここから
イコシアン（成分が a+bφ の四元数）を、成分を2倍した整数の組で持つ（zphi.py の Z）。
  ノルム N(q)=成分の二乗和。4N = X+Yφ（2倍の座標なので4倍）。
  共役 q*：φ を 1−φ に替える（残渣の側）。
  8次元の長さ² E(q) = X/2（N の「φ のつかない部分」の2倍）。
基準（走らせる前に決めた）
  H0 単位イコシアン120個とその φ 倍120個が張る格子は、整数の基底8本を持ち、長さ² E で測ると
     グラム行列の行列式が1、対角がすべて偶数（偶ユニモジュラー）［知識：8次元でこの性質を持つのは E₈ だけ］
  H1 長さ² E が2の240本は、立体の側（4次元）で長さ² 1 の120本と φ² の120本に分かれ、
     どちらも600胞体（各頂点の最近接12・辺720）。二つは大きさがちょうど φ 倍
  H2 窓（残渣の側の600胞体、頂点は単位イコシアンの共役）に入る点だけを採ると、
     採った点どうしの最短距離² は d0 で、窓を 13/10 倍にすると d0 より近い組が出る
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
def det(M):
    from fractions import Fraction as Fr   # 行列式の確認だけ（結果は整数で比べる）
    M=[[Fr(x) for x in r] for r in M]; n=len(M); d=Fr(1)
    for i in range(n):
        k=next((r for r in range(i,n) if M[r][i]!=0),None)
        if k is None: return 0
        if k!=i: M[i],M[k]=M[k],M[i]; d=-d
        d*=M[i][i]
        for r in range(i+1,n):
            f=M[r][i]/M[i][i]
            for c in range(i,n): M[r][c]-=f*M[i][c]
    return d
dg=det(G)
print("基底",len(B),"本／グラム行列の行列式",dg,"／対角",[G[i][i] for i in range(len(B))])
h0= len(B)==8 and dg==1 and all(G[i][i]%2==0 for i in range(8))
print("H0 偶ユニモジュラー（E₈）","OK" if h0 else "NG")
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
print("立体の側 半径²1 の組：点",len(S1),"辺²×4",m1,"辺",e1,"次数",g1)
print("立体の側 半径²φ² の組：点",len(S2),"辺²×4",m2,"辺",e2,"次数",g2)
h1= e1==720 and e2==720 and g1=={12} and g2=={12} and m2==m1*P*P
print("H1 二つの600胞体で大きさは φ 倍","OK" if h1 else "NG")
# ---- H2 窓で切り出す ----
W=[star(u) for u in U]       # 窓の頂点（残渣の側、2倍の座標）
ed=min(q4(sub(a,b)) for a,b in itertools.combinations(W,2))
adj={w:set() for w in W}
for a,b in itertools.combinations(W,2):
    if q4(sub(a,b))==ed: adj[a].add(b); adj[b].add(a)
F=set()
for a in W:
    for b,c,d in itertools.combinations(sorted(adj[a],key=repr),3):
        if c in adj[b] and d in adj[b] and d in adj[c]: F.add(frozenset((a,b,c,d)))
F=[tuple(f) for f in F]
NRM=[(lambda n:(n,dot(n,f[0])))(tuple(sum((v[k] for v in f),Z()) for k in range(4))) for f in F]
print("窓（600胞体）の面",len(F))
D=1000; GAM=(Z(3,1),Z(-2,5),Z(7,-3),Z(1,-4))
def inside(x,num=1,den=1):     # (x−γ) が窓の num/den 倍の内側か（2倍の座標）
    y=tuple(c*(D*den) - g*den for c,g in zip(x,GAM))
    for n,h in NRM:
        if not ((dot(y,n) - h*(D*num)).sign()<0): return False
    return True
# 候補：0 から根の一歩で広げる（立体の側 |q|² ≤ R2、残渣の側は窓の 13/10 倍の外接球まで）
lim_in=Z(4*169)        # 4N*×100 ≤ 4×169
seen={tuple([ZERO]*4)}; todo=list(seen)
while todo:
    p=todo.pop()
    for r in R:
        q=add(p,r)
        if q in seen: continue
        if (q4(q)-Z(4*R2)).sign()>0: continue
        if (q4(star(q))*100-lim_in).sign()>0: continue
        seen.add(q); todo.append(q)
A1=[q for q in seen if inside(star(q))]
A2=[q for q in seen if inside(star(q),13,10)]
print("候補",len(seen),"窓に入る",len(A1),"窓13/10倍に入る",len(A2))
def mind(S):
    best=None; cnt=0
    for a,b in itertools.combinations(S,2):
        d=q4(sub(a,b))
        if best is None or (d-best).sign()<0: best=d; cnt=1
        elif d==best: cnt+=1
    return best,cnt
d0,c0=mind(A1); d1,c1=mind(A2)
print("窓：最短距離²×4",d0,"その組",c0)
print("窓13/10倍：最短距離²×4",d1,"その組",c1)
h2= (d1-d0).sign()<0
print("H2 窓を広げると近すぎる組が出る","OK" if h2 else "NG")
print("第3章の検定 総合","OK" if (h0 and h1 and h2) else "NG")
