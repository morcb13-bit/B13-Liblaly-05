"""第7章：切り口に開く第一層 ── 円環5つは φ 倍に開く
タイリングは big/cp2_tiling.json（cp2_tiling.py 900 600）。数はすべて a+bφ の整数。
前の本の「第一層」：円環5つが五芒星を囲み、隣どうしの円環が五角形を2枚ずつ共有する並び（五角形40枚）。
円環の中心どうしの距離² ＝ (φ+2) × 円環の半径²。
基準（走らせる前に決めた）
  K1 軸が5回軸 A の向きの輪について、真ん中の面で切ると、どの輪も五角形10枚の円環になる（同じ面に複数の円環が乗る）
  K2 同じ面に乗る円環5つの中心が正五角形をなす組を、大きさ（辺² ＝ 第一層の辺² × φ^k）で分けると、
     k=0（第一層そのもの）は 0 組、k=2（φ 倍に開いた第一層）は 1 組以上ある
  K3 φ 倍に開いた組の各円環を、組の中心に向かって距離の (1−1/φ) だけ平行に寄せる（中心までの距離を 1/φ 倍にする）と、
     隣どうしの円環がちょうど2枚ずつ五角形を共有し（頂点まで完全に一致）、隣でない円環は共有せず、
     五角形は 50 枚から 40 枚になる ── 第一層と同じ並び
  負の対照 寄せないとき・1/φ² 倍まで寄せるときは、K3 の共有の形にならない
"""
import json, itertools, io, contextlib
from collections import Counter, defaultdict
from zphi import Z, GEN, dot
with contextlib.redirect_stdout(io.StringIO()):
    import ring5x10 as R
    import cut_plane as CP
P=Z(0,1); key=lambda z:z.a+1.6180339887*z.b
A=GEN[2]; AA=dot(A,A); L2=(Z(1)+P*2)*2
B=set(tuple(n) for n in json.load(open("big/cp2_tiling.json"))["rts"])
G=[[dot(GEN[i],GEN[j]) for j in range(6)] for i in range(6)]
SYM=[(p,sg) for p in itertools.permutations(range(6)) for sg in itertools.product((1,-1),repeat=6)
     if p[2]==2 and all(G[p[i]][p[j]]*(sg[i]*sg[j])==G[i][j] for i in range(6) for j in range(6))]
def app(S,v):
    p,sg=S; o=[0]*6
    for k in range(6): o[p[k]]=sg[k]*v[k]
    return tuple(o)
pat=sorted(R.U); c0=pat[0]; diffs=[tuple(c[t]-c0[t] for t in range(6)) for c in pat]
rings=set()
for S in SYM:
    sd=[app(S,v) for v in diffs]
    for n0 in B:
        s=[tuple(n0[t]+v[t] for t in range(6)) for v in sd]
        if all(x in B for x in s): rings.add(frozenset(s))
rings=list(rings)
c2=lambda n: CP.phys(tuple(2*x+1 for x in n))
Sp=lambda p: tuple(AA*p[t]-dot(p,A)*A[t] for t in range(3))      # 軸に垂直な成分 ×|A|²
def d2(a,b): v=tuple(a[t]-b[t] for t in range(3)); return dot(v,v)
def add(a,b): return tuple(x+y for x,y in zip(a,b))
def sc(a,s): return tuple(x*s for x in a)
# K1
plane=defaultdict(list); sect={}; good=True
for r in rings:
    hs=sorted({dot(c2(n),A) for n in r},key=key); h=hs[0]+L2
    secs=[CP.section(c2(n),h) for n in r]
    ok=all(cr==0 and CP.regular_pentagon(on) for on,cr in secs)
    good&=ok and len(hs)==2
    plane[h].append(r); sect[r]=[frozenset(on) for on,cr in secs]
print("軸が A の輪",len(rings),"／乗る面の数",len(plane),"／一枚の面に乗る円環の数",dict(sorted(Counter(len(v) for v in plane.values()).items())))
print("K1 どの輪も真ん中の面で五角形10枚の円環","OK" if good else "NG")
# 円環の中心（×10、軸に垂直な成分）と半径²
cen={r:tuple(sum((Sp(c2(n))[t] for n in r),Z()) for t in range(3)) for r in rings}
r0=rings[0]; R2=d2(sc(Sp(c2(next(iter(r0)))),10),cen[r0])
base=R2*(P+2)
# K2
groups=defaultdict(list)
for h,rs in plane.items():
    for five in itertools.combinations(rs,5):
        pts=[cen[r] for r in five]
        dd=Counter(d2(a,b) for a,b in itertools.combinations(pts,2))
        if len(dd)!=2 or sorted(dd.values())!=[5,5]: continue
        a,b=sorted(dd,key=key)
        if b!=a*P*P: continue
        k=None; pw=Z(1)
        for e in range(12):
            if a==base*pw: k=e
            pw=pw*P
        groups[k].append((h,five))
print("同じ面の円環5つが正五角形：大きさ φ^k ごとの組の数",{k:len(v) for k,v in sorted(groups.items(),key=lambda t:-1 if t[0] is None else t[0])})
k2=groups.get(0,[])==[] and len(groups.get(2,[]))>0
print("K2 第一層そのもの（k=0）は 0 組、φ 倍に開いた組（k=2）がある","OK" if k2 else "NG")
# K3
def close(h,five,f):        # 円環 i を、(c_i−g)×f だけ組の中心へ寄せる（50倍の座標で）
    C=[sc(cen[r],5) for r in five]                       # 円環の中心×50
    g=tuple(sum((cen[r][t] for r in five),Z()) for t in range(3))   # 組の中心×50
    out=[]
    for r,c in zip(five,C):
        sh=sc(tuple(c[t]-g[t] for t in range(3)),f)
        out.append([frozenset(tuple(x-y for x,y in zip(sc(Sp(v),50),sh)) for v in pen) for pen in sect[r]])
    # 円環の隣：中心の正五角形の辺
    dd={(i,j):d2(C[i],C[j]) for i,j in itertools.combinations(range(5),2)}
    m=min(dd.values(),key=key); adj={p for p,v in dd.items() if v==m}
    share={(i,j):len(set(out[i])&set(out[j])) for i,j in itertools.combinations(range(5),2)}
    distinct=len(set(x for o in out for x in o))
    return share,adj,distinct
shrink=Z(2)-P                       # 1−1/φ ＝ 2−φ
k3=True; ctl=True
for h,five in groups.get(2,[]):
    share,adj,dist=close(h,five,shrink)
    ok=all(share[p]==(2 if p in adj else 0) for p in share) and dist==40
    k3&=ok
    s0,_,d0=close(h,five,Z(0)); s1,_,d1=close(h,five,P-1)      # 寄せない／1/φ² 倍まで寄せる（1−1/φ² ＝ φ−1）
    c_ok= not(all(s0[p]==(2 if p in adj else 0) for p in s0) and d0==40) and not(all(s1[p]==(2 if p in adj else 0) for p in s1) and d1==40)
    ctl&=c_ok
    print("  φ倍の組：寄せる→ 隣の共有",sorted(share[p] for p in adj),"隣でない共有",sorted(share[p] for p in share if p not in adj),"五角形",dist,
          "／寄せない→ 五角形",d0,"／1/φ²倍→ 隣の共有",sorted(s1[p] for p in adj),"五角形",d1)
print("K3 1/φ 倍に寄せると第一層（隣は2枚共有・40枚）","OK" if k3 else "NG")
print("負の対照（寄せない・1/φ² 倍では第一層にならない）","OK" if ctl else "NG")
print("第7章の検定 総合","OK" if (good and k2 and k3 and ctl) else "NG")
