"""第12章：殻は φ 倍で重なる ── E₈ の殻に入れ子で並ぶ600胞体
第3章と同じく、イコシアンの環（2倍の座標で a+bφ の整数4つ＝整数8個）を使う。
  8次元の長さ² E は、4N = X+Yφ の X の半分。c=a+bφ なら c² の φ のつかない部分は a²+b² なので、
  X は整数8個の二乗和そのもの。E ≤ 10 の点を、整数8個の二乗和 ≤ 20 の中から、格子に入るものとして数える。
  格子に入るかは、第3章の整数の基底（掃き出し）で引き算だけで判定する。
基準（走らせる前に決めた）
  S1 E の殻の点の数は、E=2,4,6,8,10 で 240, 2160, 6720, 17520, 30240
     ［知識：E₈ の殻の数は 240×(E/2 の約数の3乗の和)］
  S2 どの殻も、4次元の長さ² N ごとの層に分けると、どの層の点の数も 120 の倍数
  S3 N=φ^(2k) の層は、ちょうど φ^k 倍の600胞体（単位イコシアン120個を φ^k 倍した集合と一致）。
     その層が現れるのは E=2（k=0,1）、E=4（k=−1,2）、E=10（k=−2,3）だけで、E=6・8 には現れない
     （φ^(2k) の φ のつかない部分は 1,1,2,2,5,5 … と奇数番目のフィボナッチ数になり、E はその2倍）
  負の対照 鏡に映した側の単位イコシアン（奇置換の96個を含む120個）は、N=1 の層に入らない（入るのは24個だけ）
"""
import itertools, collections
from zphi import Z
P=Z(0,1); ONE=Z(1); ZERO=Z(0)
def q4(q):
    s=Z()
    for c in q: s=s+c*c
    return s
def scl(q,s): return tuple(c*s for c in q)
def even(p): return sum(1 for i in range(4) for j in range(i+1,4) if p[i]>p[j])%2==0
def build(parity):
    U=set()
    for i in range(4):
        for s in (2,-2):
            v=[ZERO]*4; v[i]=Z(s); U.add(tuple(v))
    for sg in itertools.product((1,-1),repeat=4): U.add(tuple(Z(x) for x in sg))
    for perm in itertools.permutations(range(4)):
        if even(perm)!=parity: continue
        for sg in itertools.product((1,-1),repeat=3):
            v=[ZERO]*4; vals=[ZERO,ONE*sg[0],P*sg[1],(P-ONE)*sg[2]]
            for k in range(4): v[perm[k]]=vals[k]
            U.add(tuple(v))
    return U
U=build(True); Um=build(False)
def ints(q): return [x for c in q for x in (c.a,c.b)]
def unints(v): return tuple(Z(v[2*k],v[2*k+1]) for k in range(4))
def hnf(rows):
    rows=[r[:] for r in rows]; basis=[]; n=len(rows[0])
    for col in range(n):
        while True:
            piv=[r for r in rows if r[col]!=0]
            if len(piv)<=1: break
            piv.sort(key=lambda r:abs(r[col])); p=piv[0]
            for r in rows:
                if r is not p and r[col]!=0:
                    m=r[col]//p[col]
                    for k in range(n): r[k]-=m*p[k]
        piv=[r for r in rows if r[col]!=0]
        if piv: p=piv[0]; basis.append(p[:]); rows=[r for r in rows if r is not p]
    return basis
BASIS=hnf([ints(u) for u in U]+[ints(scl(u,P)) for u in U])
LEAD=[next(k for k in range(8) if b[k]!=0) for b in BASIS]
def inlat(v):
    v=list(v)
    for b,c in zip(BASIS,LEAD):
        if v[c]%b[c]!=0: return False
        m=v[c]//b[c]
        for k in range(8): v[k]-=m*b[k]
    return not any(v)
# 整数8個で二乗和 ≤ 20 を全部たどる
LIM=20; shell=collections.defaultdict(list)
def rec(pre,s):
    if len(pre)==8:
        if s and inlat(pre): shell[s//2].append(unints(pre))
        return
    k=0
    while k*k+s<=LIM:
        for x in ((k,-k) if k else (0,)): rec(pre+[x],s+k*k)
        k+=1
rec([],0)
cnt={e:len(shell[e]) for e in sorted(shell)}
s1= cnt=={2:240,4:2160,6:6720,8:17520,10:30240}
print("S1 殻の点の数（E=2,4,6,8,10）",[cnt.get(e) for e in (2,4,6,8,10)],"→","OK" if s1 else "NG")
layers={}
for n in sorted(shell):
    c=collections.Counter(q4(q) for q in shell[n]); layers[n]=c
    print("   E=%d の層（4N: 点の数）"%n,", ".join(f"{repr(k)}: {v}" for k,v in sorted(c.items())))
s2= all(v%120==0 for c in layers.values() for v in c.values())
print("S2 どの層の点の数も120の倍数",s2,"→","OK" if s2 else "NG")
# φ^k
def phik(k):
    z=ONE; m=P if k>=0 else P-ONE
    for _ in range(abs(k)): z=z*m
    return z
where={}; same={}
for k in range(-2,4):
    f=phik(k); target=q4((f*2,ZERO,ZERO,ZERO))       # 4φ^(2k)
    hits=[n for n in layers if target in layers[n]]
    where[k]=hits
    S=set(scl(u,f) for u in U)
    n=hits[0] if hits else None
    same[k]= n is not None and set(q for q in shell[n] if q4(q)==target)==S
ladder=[(k,where[k],same[k]) for k in range(-2,4)]
print("   φ^k 倍の600胞体：k, 現れる殻 E, 層＝φ^k×(単位イコシアン120個)",ladder)
# φ^(2k) の形の層が E=6・8 にないこと
pw=set(q4((phik(k)*2,ZERO,ZERO,ZERO)) for k in range(-6,7))
other=[n for n in (6,8) if any(k in pw for k in layers[n])]
s3= where=={-2:[10],-1:[4],0:[2],1:[2],2:[4],3:[10]} and all(same.values()) and other==[]
print("S3 φ^k 倍の600胞体は E=2（k=0,1）・4（k=−1,2）・10（k=−2,3）にだけ現れ、層と一致",s3,"／E=6・8 にある φ^(2k) の層",other,"→","OK" if s3 else "NG")
L1=set(q for q in shell[2] if q4(q)==Z(4))
mi=len(Um & L1)
neg= len(L1)==120 and mi==24
print("負の対照 鏡に映した側の120個のうち N=1 の層に入るもの",mi,"／N=1 の層の点",len(L1),"→","OK" if neg else "NG")
print("第12章の検定 総合","OK" if (s1 and s2 and s3 and neg) else "NG")
