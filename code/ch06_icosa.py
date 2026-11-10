"""第6章：輪の集まりは正二十面体 ── 30面体12個に輪が6つ重なる
タイリングは cp2_tiling.py 900 600 の出力（big/cp2_tiling.json）。数はすべて a+bφ の整数。
基準（走らせる前に決めた）
  I1 タイリングの中の10個の輪どうしは、30面体を 0 個か 8 個だけ共有する。8個共有でつないだ集まりには、
     輪が6つ・30面体が12個あり、どの30面体も集まりの輪のうち5つに入る
  I2 集まりの12個の中心は正二十面体：距離²は3種類で組の数は 30/30/6、比は 1 : φ² : (2+φ)。
     最短の30組はどれも面で接する組（6次元の差が ±1 を4つ）
  I3 集まりの中心（12個の中心の平均）は、窓で採った点（タイリングの頂点）で、そこに30面体はない
  I4 6つの輪はそれぞれ、集まりのうち自分に入らない2個（軸の両端）を結ぶ直線を軸にもつ（両端の中点＝集まりの中心）
  I5 集まりどうしは30面体を共有しない。最も近い集まりの中心どうしの距離²は、中心から30面体までの距離²の φ² 倍
  負の対照 輪1つと、その集まりに入らない30面体2個を合わせた12個は、I2 の距離の組にならない
"""
import json, itertools, io, contextlib
from collections import Counter
from zphi import Z, GEN, dot
with contextlib.redirect_stdout(io.StringIO()):
    import ring5x10 as R
P=Z(0,1); key=lambda z:z.a+1.6180339887*z.b
d=json.load(open("big/cp2_tiling.json")); B=set(tuple(n) for n in d["rts"]); PTS=set(tuple(p) for p in d["pts"])
G=[[dot(GEN[i],GEN[j]) for j in range(6)] for i in range(6)]
SYM=[(p,sg) for p in itertools.permutations(range(6)) for sg in itertools.product((1,-1),repeat=6)
     if all(G[p[i]][p[j]]*(sg[i]*sg[j])==G[i][j] for i in range(6) for j in range(6))]
def app(S,v):
    p,sg=S; o=[0]*6
    for k in range(6): o[p[k]]=sg[k]*v[k]
    return tuple(o)
# 輪の模様（kabai_ring の中心 n+1/2 → 塊の基点 n に直す：中心 c=n+½ なので 2c−1 ＝ 2n）
pat=sorted(R.U); c0=pat[0]; diffs=[tuple(c[t]-c0[t] for t in range(6)) for c in pat]
rings=set()
for S in SYM:
    sd=[app(S,v) for v in diffs]
    for n0 in B:
        s=[tuple(n0[t]+v[t] for t in range(6)) for v in sd]
        if all(x in B for x in s): rings.add(frozenset(s))
rings=list(rings)
sh=Counter(len(a&b) for a,b in itertools.combinations(rings,2))
print("30面体の塊",len(B),"／10個の輪",len(rings),"／輪どうしの共有数",dict(sorted(sh.items())))
# 8個共有でつなぐ
par=list(range(len(rings)))
def f(i):
    while par[i]!=i: par[i]=par[par[i]]; i=par[i]
    return i
for i,j in itertools.combinations(range(len(rings)),2):
    if len(rings[i]&rings[j])==8: par[f(i)]=f(j)
comp={}
for i in range(len(rings)): comp.setdefault(f(i),[]).append(i)
clusters=[]
for c in comp.values():
    if len(c)<2: continue
    rts=set().union(*[rings[i] for i in c])
    cnt=Counter(x for i in c for x in rings[i])
    clusters.append((c,rts,cnt))
print("8個共有でつないだ集まり（輪2つ以上）",len(clusters),"／輪の数ごと",dict(Counter(len(c) for c,_,_ in clusters)))
full=[(c,r,cnt) for c,r,cnt in clusters if len(c)==6]
i1= sh.keys()<={0,8} and len(full)>0 and all(len(r)==12 and set(cnt.values())=={5} for c,r,cnt in full)
print("I1 共有は0か8、輪6つの集まりは30面体12個で各5つの輪に入る","OK" if i1 else "NG", "（輪6つの集まり",len(full),"）")
def phys2(n): return tuple(sum((GEN[k][t]*(2*n[k]+1) for k in range(6)),Z()) for t in range(3))   # 30面体の中心×2
def d2(a,b): return dot(tuple(a[t]-b[t] for t in range(3)),tuple(a[t]-b[t] for t in range(3)))
step=lambda a,b: sorted(abs(x-y) for x,y in zip(a,b))==[0,0,1,1,1,1]
def icosa(ns):
    C=[phys2(n) for n in ns]
    dd=Counter(d2(C[i],C[j]) for i,j in itertools.combinations(range(12),2))
    ks=sorted(dd,key=key)
    if [dd[k] for k in ks]!=[30,30,6]: return False,dd
    a,b,c=ks
    return (b==a*P*P and c==a*(P+2)),dd
i2=True; i3=True; i4=True; centers=[]
for c,r,cnt in full:
    ns=sorted(r); ok,dd=icosa(ns)
    C=[phys2(n) for n in ns]; a=min(dd,key=key)
    edges=[(ns[i],ns[j]) for i,j in itertools.combinations(range(12),2) if d2(C[i],C[j])==a]
    ok=ok and all(step(x,y) for x,y in edges)
    i2&=ok
    s=[sum(n[k] for n in ns)+6 for k in range(6)]      # 中心×12 ＝ Σ(n+½)×12/12 … Σ(2n+1)/2 を12で割る
    cen=tuple(x//12 for x in s) if all(x%12==0 for x in s) else None
    i3&= cen is not None and cen in PTS and cen not in B and tuple(x-1 for x in cen) not in B
    # 中心の物理座標×2 は phys(2·cen)
    cen2=tuple(sum((GEN[k][t]*(2*cen[k]) for k in range(6)),Z()) for t in range(3)) if cen else None
    centers.append((cen,cen2,ns))
    for i in c:
        poles=[n for n in ns if n not in rings[i]]
        mid=tuple(phys2(poles[0])[t]+phys2(poles[1])[t] for t in range(3))
        i4&= len(poles)==2 and mid==tuple(2*x for x in cen2)
print("I2 正二十面体（30/30/6、比 1:φ²:(2+φ)、最短は面で接する組）","OK" if i2 else "NG")
print("I3 集まりの中心はタイリングの頂点で、30面体はない","OK" if i3 else "NG")
print("I4 6つの輪の軸は、入らない2個を結ぶ直線（中点＝中心）","OK" if i4 else "NG")
# I5
shared=sum(1 for x,y in itertools.combinations(full,2) if x[1]&y[1])
rad=d2(centers[0][1],phys2(centers[0][2][0]))          # 中心→30面体（×4）
best=None
for (a,ac,_),(b,bc,_) in itertools.combinations(centers,2):
    dd=d2(ac,bc)
    if best is None or key(dd)<key(best): best=dd
r_all={d2(cc,phys2(n)) for _,cc,ns in centers for n in ns}
i5= shared==0 and len(r_all)==1 and best==rad*P*P if len(centers)>1 else False
print("I5 集まりどうしの共有",shared,"／中心→30面体の距離²×4",rad,"（種類",len(r_all),"）／最も近い中心どうし距離²×4",best,"／φ²倍",best==rad*P*P if best else None)
print("I5","OK" if i5 else "NG")
# 負の対照
c,r,cnt=full[0]; ring0=sorted(rings[c[0]])
others=[n for n in B if n not in r]
ctl_ok=0; tried=0
for x,y in itertools.combinations(sorted(others,key=lambda n:key(d2(phys2(n),centers[0][1])))[:12],2):
    tried+=1; ok,_=icosa(ring0+[x,y]); ctl_ok+=ok
print("負の対照 輪1つ＋集まりの外の30面体2個（近い順に12個から2個、",tried,"通り）で正二十面体になったもの",ctl_ok)
ctl= ctl_ok==0
print("第6章の検定 総合","OK" if (i1 and i2 and i3 and i4 and i5 and ctl) else "NG")
