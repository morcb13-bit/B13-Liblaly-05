"""第16章：φ の入る入れ子 ── 整数の格子を半分の次元に写すと、いちばん短い棒が φ² の比で二つに分かれる
整数の格子（2n 次元）を n 次元に写したとき、いちばん短い棒（最短のベクトル）の長さ² を a＋bφ の整数で比べる。
  5→2：Z^5 の単位 e_k を、5乗で1に戻る数 ζ（ζ^5＝1）の ζ^k に写す。ζ の多項式（整数5つ）で計算し、
        実数になったものを a＋bφ に直す（ζ＋ζ⁴ を x と置き、x²＋x−1＝0 を確かめてから x＝φ−1 とする）。
        いちばん短い棒は、和が0の4次元の格子（A4）の e_i − e_j、20本。
  6→3：第2章と同じ生成子 GEN（(±1,φ,0) の巡回）。いちばん短い棒は e_i ± e_j（D6）、60本。単位 ±e_k は12本。
  8→4：第12章の E₈。いちばん短い棒は E＝2 の240本（4次元の長さで層に分ける）。
基準（走らせる前に決めた）
  R0 5次元・6次元の単位 ±e_k（10本・12本）は、写すと長さ² が1種類（正十角形・正二十面体）。
  R1 A4 の20本は、写すと長さ² が2種類で 10本・10本、大きいほうは小さいほうのちょうど φ² 倍。
  R2 D6 の60本は、写すと長さ² が2種類で 30本・30本、比はちょうど φ²。
  R3 E₈ の240本は、写すと長さ² が2種類で 120本・120本、比はちょうど φ²。
  負の対照 Z^4 を4乗で1に戻る i に写すと（e_k → i^k）、e_i − e_j の12本の長さ² の比は 2 で、φ² ではない。
"""
import itertools, io, contextlib, collections
from zphi import Z, GEN, dot
P=Z(0,1)
# ── ζ^5＝1 の多項式（係数5つ、全部に同じ数を足しても同じ値）
def zmul(a,b):
    c=[0]*5
    for i in range(5):
        for j in range(5): c[(i+j)%5]+=a[i]*b[j]
    return c
def zconj(a): return [a[0]]+[a[(5-k)%5] for k in range(1,5)]
def zsub(a,b): return [x-y for x,y in zip(a,b)]
def zeta(k): v=[0]*5; v[k%5]=1; return v
def real_to_phi(c):
    # 実数なら c1=c4, c2=c3（全体に同じ数を足しても変わらない差で見る）
    assert c[1]-c[4]==0 and c[2]-c[3]==0, c
    # c0 + c1·x + c2·(ζ²+ζ³)，ζ²+ζ³ ＝ −1 − x，x ＝ φ−1
    x=P-Z(1)
    return Z(c[0])+x*c[1]+(Z(-1)-x)*c[2]
# x ＝ ζ＋ζ⁴ が x²＋x−1＝0 を満たすこと（ζ の多項式で）
X=[0,1,0,0,1]; q=zmul(X,X); q=[q[k]+X[k]-(1 if k==0 else 0) for k in range(5)]
xok= len(set(q))==1         # 全部同じ数＝0
def norm(z): return real_to_phi(zmul(z,zconj(z)))
def classes(vals):
    c=collections.Counter(vals); return dict(sorted(c.items()))
def two_phi2(cl):
    ks=sorted(cl);
    return len(ks)==2 and cl[ks[0]]==cl[ks[1]] and ks[0]*(P*P)==ks[1]
# R0
u5=classes([norm([s*x for x in zeta(k)]) for k in range(5) for s in (1,-1)])
u6=classes([dot(tuple(c*s for c in GEN[k]),tuple(c*s for c in GEN[k])) for k in range(6) for s in (1,-1)])
r0= xok and len(u5)==1 and sum(u5.values())==10 and len(u6)==1 and sum(u6.values())==12
print("R0 x＝ζ＋ζ⁴ は x²＋x−1＝0",xok,"／5次元の単位10本の長さ²",{repr(k):v for k,v in u5.items()},"／6次元の単位12本",{repr(k):v for k,v in u6.items()},"→","OK" if r0 else "NG")
# R1
a4=classes([norm(zsub(zeta(i),zeta(j))) for i in range(5) for j in range(5) if i!=j])
r1=two_phi2(a4)
print("R1 A4 の20本の長さ²",{repr(k):v for k,v in a4.items()},"／比 φ²",r1,"→","OK" if r1 else "NG")
# R2
def add(a,b,s): return tuple(x+y*s for x,y in zip(a,b))
d6=classes([dot(add(GEN[i],GEN[j],s),add(GEN[i],GEN[j],s)) for i,j in itertools.combinations(range(6),2) for s in (1,-1) for t in (1,-1)])
r2=two_phi2(d6)
print("R2 D6 の60本の長さ²",{repr(k):v for k,v in d6.items()},"／比 φ²",r2,"→","OK" if r2 else "NG")
# R3
with contextlib.redirect_stdout(io.StringIO()):
    import ch12_shells as S
e8=classes([S.q4(q) for q in S.shell[2]])
r3=two_phi2(e8)
print("R3 E₈ の240本の4次元の長さ²（×4）",{repr(k):v for k,v in e8.items()},"／比 φ²",r3,"→","OK" if r3 else "NG")
# 対照：i^4＝1、ガウスの整数 (a,b) で
def g(k): return [(1,0),(0,1),(-1,0),(0,-1)][k%4]
c4=collections.Counter((g(i)[0]-g(j)[0])**2+(g(i)[1]-g(j)[1])**2 for i in range(4) for j in range(4) if i!=j)
ks=sorted(c4); neg= len(ks)==2 and ks[1]==2*ks[0] and not (Z(ks[0])*(P*P)==Z(ks[1]))
print("負の対照 Z^4→i^k の e_i−e_j 12本の長さ²",dict(sorted(c4.items())),"／比",f"{ks[1]}/{ks[0]}","φ² ではない","→","OK" if neg else "NG")
print("第16章の検定 総合","OK" if (r0 and r1 and r2 and r3 and neg) else "NG")
