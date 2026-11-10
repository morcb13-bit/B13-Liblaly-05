# 序章：6次元立方体の影は入れ子になっている
# 6次元立方体の頂点（±1 の6つ組）64個を、正二十面体の5回軸6本で3次元に写す。
# 数は a+bφ を (a,b) の整数の組で持つ（φ²=φ+1）。比べるのは整数の等号だけ。
# 基準（走らせる前に決めた）
#   S1 影は殻に分かれ、20点の殻は正十二面体（辺30本・各頂点に3本・五角形の対角線²＝φ²×辺²）、
#      12点の殻は正二十面体（辺30本・各頂点に5本）
#   S2 6次元立方体の中の立方体（面にない立方体）のうち、正十二面体の殻に正しい立方体として写るものが5個、
#      その5個はどれも正四面体と、中心について反転した正四面体の組に分かれる
#   負の対照 φ を 3/2 に替えた写し方では、正十二面体の殻が一つも出ない
from itertools import product, combinations
import json
def mul(x,y): a,b=x; c,d=y; return (a*c+b*d, a*d+b*c+b*d)
def add(x,y): return (x[0]+y[0],x[1]+y[1])
def neg(x): return (-x[0],-x[1])
Z=(0,0)
def less(x,y):            # a+bφ < c+dφ を整数だけで判定
    a=x[0]-y[0]; b=x[1]-y[1]          # a+bφ<0 ?
    # 2(a+bφ) = (2a+b) + b√5
    p=2*a+b; q=b
    if p<=0 and q<=0: return p<0 or q<0
    if p>=0 and q>=0: return False
    if p<0: return p*p>5*q*q          # q>0：p+q√5<0 ⇔ p²>5q²
    return 5*q*q>p*p                  # p>0,q<0
def keysort(xs):
    xs=list(xs)
    for i in range(len(xs)):
        for j in range(len(xs)-1-i):
            if less(xs[j+1],xs[j]): xs[j],xs[j+1]=xs[j+1],xs[j]
    return xs
def minimum(xs):
    m=None
    for x in xs:
        if m is None or less(x,m): m=x
    return m
def setup(o,f):
    m=neg(o)
    return [(o,f,Z),(m,f,Z),(Z,o,f),(Z,m,f),(f,Z,o),(f,Z,m)]
V=list(product((1,-1),repeat=6))
def shadow(GEN):
    P={}
    for s in V:
        p=[Z,Z,Z]
        for si,g in zip(s,GEN):
            for k in range(3): p[k]=add(p[k],(si*g[k][0],si*g[k][1]))
        P[s]=tuple(p)
    return P
def n2(p,q):
    r=Z
    for k in range(3):
        d=(p[k][0]-q[k][0],p[k][1]-q[k][1]); r=add(r,mul(d,d))
    return r
O=(Z,Z,Z)
def shells_of(P):
    sh={}
    for v in V: sh.setdefault(n2(P[v],O),[]).append(v)
    return sh
def shape(pts):
    D={}
    for a,b in combinations(pts,2): D.setdefault(n2(a,b),[]).append((a,b))
    e2=minimum(D); E=D[e2]; deg={p:0 for p in pts}
    for a,b in E: deg[a]+=1; deg[b]+=1
    return e2,E,sorted(set(deg.values())),D
def sub(a,b): return tuple(x-y for x,y in zip(a,b))
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def inner_cubes():
    Vs=set(V); cubes=set()
    for v0 in V:
        Dv=[sub(w,v0) for w in V if w!=v0]
        for a,b,e in combinations(Dv,3):
            L=dot(a,a)
            if L<=4 or dot(b,b)!=L or dot(e,e)!=L or dot(a,b) or dot(a,e) or dot(b,e): continue
            pts=set()
            for m in range(8):
                p=list(v0)
                for j,x in enumerate((a,b,e)):
                    if m>>j&1: p=[u+y for u,y in zip(p,x)]
                pts.add(tuple(p))
            if pts<=Vs and not any(len({p[i] for p in pts})==1 for i in range(6)): cubes.add(frozenset(pts))
    return cubes
def show(x): return f"{x[0]}{'+' if x[1]>=0 else '−'}{abs(x[1])}φ" if x[1] else f"{x[0]}"

print("== 黄金比の写し方（5回軸6本）==")
P=shadow(setup((1,0),(0,1))); sh=shells_of(P); cubes=inner_cubes()
print("6次元立方体の頂点",len(V),"／中の立方体",len(cubes))
ok1=True; ok2=True; out={"shells":[],"edges6":[]}
idx={v:i for i,v in enumerate(V)}
for k in keysort(sh):
    vs=sh[k]; pts=[P[v] for v in vs]
    e2,E,deg,D=shape(pts)
    kind="正十二面体" if (len(pts)==20 and deg==[3]) else "正二十面体" if (len(pts)==12 and deg==[5]) else "?"
    line=f"長さ² {show(k):>8}  点 {len(pts):>2}  辺 {len(E)}  各頂点の辺 {deg[0]}  → {kind}"
    good=[]
    if kind=="正十二面体":
        diag=mul((1,1),e2); line+=f"  対角線²＝φ²×辺² の組 {len(D.get(diag,[]))}"
        ok1&= len(E)==30 and len(D.get(diag,[]))==60
        ptset=set(pts)
        for c in cubes:
            im={P[v] for v in c}
            if len(im)==8 and im<=ptset:
                ds=[n2(a,b) for a,b in combinations(im,2)]
                if ds.count(minimum(ds))==12: good.append(c)
        pair=0
        for c in good:
            im=[P[v] for v in c]; a=im[0]; e=minimum(n2(a,b) for b in im[1:])
            T1={a}|{b for b in im if n2(a,b)==(2*e[0],2*e[1])}; T2=set(im)-T1
            if len(T1)==4 and {tuple(neg(x) for x in p) for p in T1}==T2: pair+=1
        line+=f"\n      写る立方体 {len(good)}（正四面体と反転した正四面体の組 {pair}）"
        ok2&= len(good)==5 and pair==5
    elif kind=="正二十面体": ok1&= len(E)==30
    else: ok1=False
    print(line)
    out["shells"].append({"n2":k,"kind":kind,"pts":[idx[v] for v in vs],
        "edges":[[idx[[v for v in vs if P[v]==a][0]],idx[[v for v in vs if P[v]==b][0]]] for a,b in E],
        "cubes":[[idx[v] for v in c] for c in good]})
print("S1 殻の形","OK" if ok1 else "NG")
print("S2 正十二面体の5つの立方体","OK" if ok2 else "NG")
print("\n== 負の対照：φ を 3/2 に替える（座標を2倍して整数で持つ）==")
Pc=shadow(setup((2,0),(3,0))); shc=shells_of(Pc); nd=0
for k in keysort(shc):
    pts=[Pc[v] for v in shc[k]]; e2,E,deg,D=shape(pts)
    isd=len(pts)==20 and deg==[3]; nd+=isd
    print(f"長さ² {show(k):>6}  点 {len(pts):>2}  辺 {len(E)}  各頂点の辺 {deg}")
print("正十二面体の殻",nd,"→","OK（出ない）" if nd==0 else "NG")
print("\n序章の検定 総合","OK" if ok1 and ok2 and nd==0 else "NG")
# 図のデータ（点は a+bφ の組のまま渡す）
out["pts"]=[[list(c) for c in P[v]] for v in V]
out["edges6"]=[[i,idx[tuple(-x if j==k else x for j,x in enumerate(v))]] for i,v in enumerate(V) for k in range(6) if v[k]==1]
json.dump(out,open("ch00_shadow.json","w"))
