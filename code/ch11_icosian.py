"""第11章：イコシアンと600胞体 ── 掛け算で閉じる120個
単位イコシアン120個を、成分を2倍した a+bφ の整数の組で持つ（zphi.py の Z）。
  積は (2p)(2q) = 4pq。表の引き当ては、相手を2倍した 2·(2r) と比べるだけにする（割り算をしない）。
  内積も2倍の座標どうしで取るので、本当の内積の4倍が出る（最近接 φ/2 → 2φ、直交 0）。
基準（走らせる前に決めた）
  G1 120個のどれも 4N=4。14400通りの積がすべて120個の中にあり、1 を含み、どれも共役が逆元として中にある。
  G2 位数の分布は {1:1, 2:1, 3:20, 4:30, 5:24, 6:20, 10:24}。ab≠ba の組は 13320。
  G3 位数10の元 g の累乗10個は大円の十角形（隣どうしの内積が最近接の 2φ）。
     その十角形に120個を左から掛けると、ちょうど12本の十角形になり、互いに共有せず120個を埋める。
  G4 12本の十角形を、二本の間の最大の内積で分けると、どの一本にも、すべての頂点が直交する相手が
     ちょうど1本、近い相手5本、遠い相手5本。近い相手を辺にすると、頂点12・次数5・辺30で、
     直交する相手は辺をたどって3歩目にだけ現れる（正二十面体の頂点と向かいの頂点）。
  G5 120個のどれも、q v q̄ で正二十面体の頂点12個を12個に写す。q と −q は同じ写し方で、写し方は60通り。
  負の対照 偶置換の96個に奇置換の96個を足した216個では、積が集合からはみ出す。
修正の記録：最初の走り（ch11_icosian_v1.py／_v1_out.txt）は正二十面体の頂点を (1,φ,0) の巡回で置き、
  G5 が NG（写し方13通り）。これはイコシアンの向きと鏡映の側の正二十面体だったため、
  (φ,1,0) の巡回に置き直した。基準 G5 そのものは変えていない。鏡映の側の結果は記録として出す。
"""
import itertools, collections
from zphi import Z
P=Z(0,1); ONE=Z(1); ZERO=Z(0)
def q4(q):
    s=Z()
    for c in q: s=s+c*c
    return s
def dot(p,q):
    s=Z()
    for a,b in zip(p,q): s=s+a*b
    return s
def mul(p,q):
    a1,b1,c1,d1=p; a2,b2,c2,d2=q
    return (a1*a2-b1*b2-c1*c2-d1*d2, a1*b2+b1*a2+c1*d2-d1*c2,
            a1*c2-b1*d2+c1*a2+d1*b2, a1*d2+b1*c2-c1*b2+d1*a2)
def conj(q): return (q[0],-q[1],-q[2],-q[3])
def twice(q): return tuple(c*2 for c in q)
def even(p): return sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j])%2==0
def build(parity):
    U=set()
    for i in range(4):
        for s in (2,-2):
            v=[ZERO]*4; v[i]=Z(s); U.add(tuple(v))
    for sg in itertools.product((1,-1),repeat=4): U.add(tuple(Z(x) for x in sg))
    for perm in itertools.permutations(range(4)):
        if even(perm)!=parity and parity is not None: continue
        for sg in itertools.product((1,-1),repeat=3):
            v=[ZERO]*4; vals=[ZERO,ONE*sg[0],P*sg[1],(P-ONE)*sg[2]]
            for k in range(4): v[perm[k]]=vals[k]
            U.add(tuple(v))
    return sorted(U,key=repr)
U=build(True); n=len(U)
idx={twice(u):k for k,u in enumerate(U)}          # 2·(2r) → 番号
E1=(Z(2),ZERO,ZERO,ZERO); one=U.index(E1)
T=[[idx.get(mul(U[i],U[j])) for j in range(n)] for i in range(n)]
closed=all(T[i][j] is not None for i in range(n) for j in range(n))
inv_ok=all(conj(u) in U for u in U) and all(T[k][U.index(conj(U[k]))]==one for k in range(n))
g1= n==120 and all(q4(u)==Z(4) for u in U) and closed and inv_ok
print("G1 単位イコシアン",n,"／4N≠4",sum(1 for u in U if q4(u)!=Z(4)),"／積14400通りがすべて中に",closed,"／共役が逆元",inv_ok,"→","OK" if g1 else "NG")
def order(k):
    m,c=1,k
    while c!=one: c=T[c][k]; m+=1
    return m
od=collections.Counter(order(k) for k in range(n))
nc=sum(1 for i in range(n) for j in range(n) if T[i][j]!=T[j][i])
g2= dict(od)=={1:1,2:1,3:20,4:30,5:24,6:20,10:24} and nc==13320
print("G2 位数の分布",dict(sorted(od.items())),"／ab≠ba の組",nc,"→","OK" if g2 else "NG")
# G3
g=min((k for k in range(n) if order(k)==10),key=lambda k:repr(U[k]))
pw=[one];
for _ in range(9): pw.append(T[pw[-1]][g])
NEAR=P*2
ring_ok=all(dot(U[pw[i]],U[pw[(i+1)%10]])==NEAR for i in range(10))
decs=set(frozenset(T[x][p] for p in pw) for x in range(n))
cover=sorted(set().union(*decs))==list(range(n)) and all(len(d)==10 for d in decs)
g3= ring_ok and len(decs)==12 and cover
print("G3 g の累乗10個：隣の内積がすべて 2φ",ring_ok,"／左から掛けてできる十角形",len(decs),"本／共有なしで120個を埋める",cover,"→","OK" if g3 else "NG")
# G4
D=sorted(decs,key=lambda d:min(d))
def mx(A,B): return max((dot(U[a],U[b]) for a in A for b in B))
def allortho(A,B): return all(dot(U[a],U[b])==ZERO for a in A for b in B)
vals=collections.Counter(); rows=[]
M=[[None]*12 for _ in range(12)]
for i in range(12):
    for j in range(12):
        if i!=j: M[i][j]=mx(D[i],D[j]); vals[repr(M[i][j])]+=1
kinds=sorted({M[i][j] for i in range(12) for j in range(12) if i!=j})
ortho=[[j for j in range(12) if j!=i and allortho(D[i],D[j])] for i in range(12)]
hi=kinds[-1] if kinds else None
nb=[[j for j in range(12) if j!=i and M[i][j]==hi] for i in range(12)]
deg=collections.Counter(len(x) for x in nb); edges=sum(len(x) for x in nb)//2
def dist(i):
    d={i:0}; q=[i]
    for v in q:
        for w in nb[v]:
            if w not in d: d[w]=d[v]+1; q.append(w)
    return d
far3=[[j for j,k in dist(i).items() if k==3] for i in range(12)]
split=collections.Counter((len(ortho[i]),sum(1 for j in range(12) if j!=i and M[i][j]==hi),
      sum(1 for j in range(12) if j!=i and M[i][j]!=hi and j not in ortho[i])) for i in range(12))
g4= len(kinds)==3 and all(len(o)==1 for o in ortho) and dict(split)=={(1,5,5):12} and dict(deg)=={5:12} and edges==30 \
    and all(far3[i]==ortho[i] for i in range(12))
print("G4 二本の間の最大の内積の種類",[repr(k) for k in kinds],"／（直交,近い,遠い）",dict(split),"／近い相手の次数",dict(deg),"辺",edges,
      "／3歩目＝直交の相手",all(far3[i]==ortho[i] for i in range(12)),"→","OK" if g4 else "NG")
# G5
def icosa(a,b):
    S=set()
    for s1,s2 in itertools.product((1,-1),repeat=2):
        for c in range(3):
            v=[ZERO]*3; v[c]=a*s1; v[(c+1)%3]=b*s2; S.add((ZERO,)+tuple(v))
    return sorted(S,key=repr)
def maps(ICO):
    ico4={tuple(c*4 for c in v):k for k,v in enumerate(ICO)}
    perms=[]; okall=True
    for u in U:
        im=[ico4.get(mul(mul(u,v),conj(u))) for v in ICO]
        if None in im or len(set(im))!=12: okall=False
        perms.append(tuple(im))
    return perms,okall
mir,mok=maps(icosa(ONE,P))
print("記録 鏡映の側の正二十面体 (1,φ,0) の巡回：12個を12個に写すもの",sum(1 for p in mir if None not in p and len(set(p))==12),"/120")
ICO=icosa(P,ONE); perms,okall=maps(ICO)
pm={};
for k,u in enumerate(U): pm[u]=perms[k]
pair=all(pm[u]==pm[tuple(-c for c in u)] for u in U)
g5= len(ICO)==12 and okall and pair and len(set(perms))==60
print("G5 正二十面体の頂点",len(ICO),"／120個すべてが12個を12個に写す",okall,"／q と −q は同じ",pair,"／写し方",len(set(perms)),"通り →","OK" if g5 else "NG")
# 対照
W=build(None); widx={twice(u) for u in W}
out=sum(1 for a in W for b in W if mul(a,b) not in widx)
neg= len(W)==216 and out>0
print("負の対照 偶置換と奇置換を合わせた",len(W),"個：積が集合からはみ出す組",out,"→","OK" if neg else "NG")
print("第11章の検定 総合","OK" if (g1 and g2 and g3 and g4 and g5 and neg) else "NG")
