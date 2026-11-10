"""第5章：円環は30面体のジグザグの輪 ── 上の段5個・下の段5個
数はすべて a+bφ の整数（zphi.py）。輪は kabai_ring.json の10個（ring5x10.py の U）。
基準（走らせる前に決めた）
  「隣」＝真上（5回軸）から見た影で、角度の順に隣り合う組（10組）
  （初版では「隣」を面で接する組すべてとしたため Z1〜Z3 が NG。面で接する組は20あり、
   上下の両隣と同じ段の一つおきの両方を含んでいた。定義を直して走らせ直した。経緯は引継書）
  Z1 輪：10個の30面体の中心（6次元の整数）。隣どうし10組は面で接し（差が ±1 を4つ持つ）、
     高さ（5回軸方向）は2種類で、隣どうしは必ず上下が入れ替わる（上の段5・下の段5）
     構造として：面で接する組は全部で20（各30面体は上下の両隣と、同じ段の一つおきの4つと接する）
  Z2 影：5回軸の真上から見ると、10個の中心は軸から同じ距離に並び、
     影の半径² ＝ φ² × 隣どうしの影の間隔²（円環の半径と五角形の間隔の比 φ と同じ）
  Z3 切り口：上の段と下の段の真ん中の面で切ると、10個とも正五角形になり、
     隣どうしの五角形はちょうど1本の辺（2頂点）を共有し、それ以外の組は頂点を共有しない ── 五角形10枚の円環
     負の対照：面を少しずらすと正五角形は出ない
  Z4 第2章のタイリング（6次元の窓で切り出した立体のペンローズ）の中に、この輪がそのまま現れる
"""
import json, itertools
from collections import Counter
from zphi import Z, GEN, dot
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):   # 読み込むときの表示を出さない
    import ring5x10 as R
    import cut_plane as CP
P=Z(0,1); U=R.U; A=GEN[2]; AA=dot(A,A)
key=lambda z:z.a+1.6180339887*z.b      # 並べる順番だけ
ok=True
# ---- 隣：影の角度の順（角度は並べる順番を決めるだけに使う） ----
import math
perp=lambda p: tuple(AA*p[t]-dot(p,A)*A[t] for t in range(3))
S={u:perp(CP.phys(u)) for u in U}
cen=tuple(sum((S[u][t] for u in U),Z()) for t in range(3))     # 10倍の中心
ang={u:math.atan2(key(S[u][1])*10-key(cen[1]),key(S[u][0])*10-key(cen[0])) for u in U}
order=sorted(U,key=lambda u:ang[u])
nb=[(order[i],order[(i+1)%10]) for i in range(10)]
# ---- Z1 ----
step=lambda a,b: sorted(abs(x-y) for x,y in zip(a,b))==[0,0,1,1,1,1]
touch=[(a,b) for a,b in itertools.combinations(U,2) if step(a,b)]
deg=Counter(x for p in touch for x in p)
nb_touch=all(step(a,b) for a,b in nb)
h=lambda u: dot(CP.phys(u),A)
hs=sorted({h(u) for u in U},key=key)
up={u for u in U if h(u)==hs[1]}
alt=all((a in up)!=(b in up) for a,b in nb)
same=sum(1 for a,b in touch if (a in up)==(b in up))
print("Z1 隣10組がすべて面で接する",nb_touch,"／隣は必ず上下が入れ替わる",alt,"／高さの種類",len(hs),"差",hs[1]-hs[0],"／上の段",len(up))
print("   構造：面で接する組",len(touch),"（各30面体",set(deg.values()),"、うち同じ段どうし",same,"）")
z1= nb_touch and alt and len(hs)==2 and len(up)==5 and len(touch)==20 and set(deg.values())=={4} and same==10; ok&=z1
# ---- Z2 ----
r2={dot(tuple(10*S[u][t]-cen[t] for t in range(3)),tuple(10*S[u][t]-cen[t] for t in range(3))) for u in U}
sp={dot(tuple(10*(S[a][t]-S[b][t]) for t in range(3)),tuple(10*(S[a][t]-S[b][t]) for t in range(3))) for a,b in nb}
print("Z2 影の半径²の種類",len(r2),"／隣どうしの影の間隔²の種類",len(sp))
z2= len(r2)==1 and len(sp)==1 and list(r2)[0]==list(sp)[0]*P*P
print("   影の半径² ＝ φ² × 間隔²",z2); ok&=z2
# ---- Z3 ----
L2=(Z(1)+P*2)*2
mid=hs[0]+L2
res={u:CP.section(CP.phys(tuple(2*x for x in u)),mid) for u in U}
allreg=all(cr==0 and CP.regular_pentagon(on) for on,cr in res.values())
share=Counter(len(set(res[a][0])&set(res[b][0])) for a,b in itertools.combinations(U,2))
adj_share=all(len(set(res[a][0])&set(res[b][0]))==2 for a,b in nb)
print("Z3 真ん中の面で切る：10個とも正五角形",allreg,"／五角形どうしが共有する頂点の数（全45組）",dict(share),"／隣10組は辺を共有",adj_share)
ctl=[CP.section(CP.phys(tuple(2*x for x in u)),mid+Z(1)) for u in U]
ctl_reg=sum(1 for on,cr in ctl if CP.regular_pentagon(on))
print("   負の対照 面を少しずらす：正五角形",ctl_reg)
z3= allreg and adj_share and share.get(2,0)==10 and ctl_reg==0; ok&=z3
# ---- Z4 ----
d=json.load(open("cp2_tiling.json")); Bset=set(tuple(n) for n in d["rts"])
G=[[dot(GEN[i],GEN[j]) for j in range(6)] for i in range(6)]
SYM=[(p,sg) for p in itertools.permutations(range(6)) for sg in itertools.product((1,-1),repeat=6)
     if all(G[p[i]][p[j]]*(sg[i]*sg[j])==G[i][j] for i in range(6) for j in range(6))]
def app(Sm,v):
    p,sg=Sm; o=[0]*6
    for k in range(6): o[p[k]]=sg[k]*v[k]
    return tuple(o)
pat=sorted(U); c0=pat[0]; diffs=[tuple(c[t]-c0[t] for t in range(6)) for c in pat]
found=set()
for Sm in SYM:
    sd=[app(Sm,v) for v in diffs]
    for n0 in Bset:
        s=[tuple(n0[t]+v[t] for t in range(6)) for v in sd]
        if all(x in Bset for x in s): found.add(frozenset(s))
print("Z4 第2章のタイリングの30面体の塊",len(Bset),"／その中の10個の輪",len(found))
z4= len(found)>0; ok&=z4
print("Z1",("OK" if z1 else "NG"),"Z2",("OK" if z2 else "NG"),"Z3",("OK" if z3 else "NG"),"Z4",("OK" if z4 else "NG"))
print("第5章の検定 総合","OK" if ok else "NG")
from zphi import cross
FACES=[]
for i,j in itertools.combinations(range(6),2):
    n=cross(GEN[i],GEN[j])
    sk=[dot(n,GEN[k]).sign() for k in range(6)]
    for sg in (1,-1):
        quad=[]
        for si,sj in ((1,1),(1,-1),(-1,-1),(-1,1)):
            s_=[sk[k]*sg for k in range(6)]; s_[i]=si; s_[j]=sj
            quad.append([[x.a,x.b] for x in tuple(sum((GEN[k][t]*s_[k] for k in range(6)),Z()) for t in range(3))])
        FACES.append(quad)
json.dump({"faces":FACES,"cen":[[[x.a,x.b] for x in CP.phys(tuple(2*c for c in u))] for u in order],"up":[int(u in up) for u in order],"axis":[[x.a,x.b] for x in A],
           "pent":[[[[c.a,c.b] for c in p] for p in res[u][0]] for u in order]},open("fig05_data.json","w"),separators=(",",":"))
