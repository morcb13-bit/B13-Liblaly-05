"""第9章：正四面体の紐 ── 3次元では閉じず、4次元で30個で閉じる
3次元の座標は「整数の組 ÷ 3^k」で持つ（分子と k の組）。比べるときは 3 の累乗を掛けてそろえるだけ。
4次元（600胞体）の座標は a+bφ の整数の組（cell600_string.py の作り）。
紐の規則：直前の3頂点の三角形を共有するもう一つの正四面体へ進む（3次元では、最古の頂点をその三角形で折り返す）。
  正四面体では、頂点を向かいの面で折り返した点 ＝ 2×(面の重心) − 頂点 ＝ (2(a+b+c) − 3p)/3。
基準（走らせる前に決めた）
  T1 3次元の紐は、1〜299段のどこでも、最初の正四面体を平行移動しただけの形に戻らない（辺の長さ²は1種類のまま）
     対照：立方体を面でつないだ紐は1段目で戻る
  T2 3次元で1本の辺のまわりに正四面体を5個並べると隙間が残る：243×cos(隙間) ＝ 241、243×(半角の sin²) ＝ 1
     対照：立方体4個、正四面体と正八面体を交互に4個は、隙間なく閉じる
  T3 600胞体：頂点120・辺720・正四面体600。どの辺のまわりにも正四面体がちょうど5個
  T4 600胞体の中で同じ規則の紐を作ると、正四面体30個でちょうど最初の正四面体に戻り、紐に現れる頂点も30個
"""
import io, contextlib, json
from itertools import combinations
ok={}
# ---------- 3次元の紐（分子と 3 の指数） ----------
def norm(v,k):
    while k>0 and all(x%3==0 for x in v): v=tuple(x//3 for x in v); k-=1
    return (v,k)
def lift(p,K): v,k=p; return tuple(x*3**(K-k) for x in v)
def refl(p,a,b,c):
    K=max(q[1] for q in (p,a,b,c))
    P_,A,B,C=(lift(q,K) for q in (p,a,b,c))
    return norm(tuple(2*(x+y+z)-3*w for x,y,z,w in zip(A,B,C,P_)),K+1)
def sub(p,q):
    K=max(p[1],q[1]); return (tuple(x-y for x,y in zip(lift(p,K),lift(q,K))),K)
def n2(d): v,k=d; return (sum(x*x for x in v),2*k)      # 長さ² ＝ 整数 ÷ 3^(2k)
T0=[((0,0,0),0),((1,1,0),0),((1,0,1),0),((0,1,1),0)]
v=list(T0)
for _ in range(300): v.append(refl(v[-4],v[-3],v[-2],v[-1]))
def shifted(P,Q):
    for q0 in Q:
        d=sub(q0,P[0])
        K=max([d[1]]+[p[1] for p in P]+[q[1] for q in Q])
        S={tuple(x+y for x,y in zip(lift(p,K),lift(d,K))) for p in P}
        if S=={lift(q,K) for q in Q}: return True
    return False
hit=[n for n in range(1,300) if shifted(v[0:4],v[n:n+4])]
e2=set()
for n in range(300):
    for i,j in combinations(range(n,n+4),2):
        s_,k2=n2(sub(v[i],v[j]))                     # 長さ² ＝ s_ ÷ 3^k2
        while k2>0 and s_%3==0: s_//=3; k2-=1     # 約分してから比べる（整数の割り算は割り切れるときだけ）
        e2.add((s_,k2))
cube=[((i,j,k),0) for i in (0,1) for j in (0,1) for k in (0,1)]
nxt=[((p[0][0]+1,p[0][1],p[0][2]),0) for p in cube]
ok["T1"]= hit==[] and len(e2)==1 and shifted(cube,nxt)
print("T1 3次元の紐：1〜299段で戻る段",hit or "なし","／辺の長さ²の種類",len(e2),"／対照 立方体の紐は1段目で戻る",shifted(cube,nxt),"→","OK" if ok["T1"] else "NG")
print("   構造：n段目に加わる頂点の座標の分母は 3^n で、約分されない（0〜300段のすべて）", all(v[n][1]==max(0,n-3) for n in range(len(v))))
# ---------- T2 辺のまわりの隙間（cos と sin/√8 を、3 の累乗を分母にした整数で持つ） ----------
def mul(x,y): return (x[0]*y[0]-8*x[1]*y[1], x[0]*y[1]+x[1]*y[0])
tet=(1,1); octa=(-1,1)          # cos 70.53° ＝ 1/3、cos 109.47° ＝ −1/3（どちらも sin ＝ √8/3）、分母 3
r=(1,0)
for _ in range(5): r=mul(r,tet)        # 分母 3^5 ＝ 243
cosgap=r[0]                              # 隙間 ＝ 360° − 5θ、cos は同じ
r2=(1,0)
for x in (tet,octa,tet,octa): r2=mul(r2,x)   # 分母 81
ok["T2"]= cosgap==241 and (243-cosgap)==2 and r2==(81,0)
print("T2 正四面体5個：243×cos(隙間) ＝",cosgap,"／243×(半角の sin²) ＝",(243-cosgap)//2,"／対照 正四面体・正八面体の交互4個：81×(cos, sin/√8) ＝",r2,"（(81,0) なら閉じる）／立方体4個は 90°×4 ＝ 360°","→","OK" if ok["T2"] else "NG")
# ---------- T3・T4 600胞体 ----------
src=open("cell600_string.py").read().split("# C2 紐")[0]
g={}
with contextlib.redirect_stdout(io.StringIO()): exec(src,g)
V,adj,tets=g["V"],g["adj"],g["tets"]
nE=sum(len(a) for a in adj)//2
around={}
for t in tets:
    for e in combinations(t,2): around[e]=around.get(e,0)+1
ok["T3"]=(len(V),nE,len(tets))==(120,720,600) and set(around.values())=={5}
print("T3 600胞体：頂点",len(V),"辺",nE,"正四面体",len(tets),"／辺のまわりの正四面体",set(around.values()),"→","OK" if ok["T3"] else "NG")
start=sorted(tets)[0]; w=list(start); close=None
for step in range(1,200):
    a,b,c=w[-3],w[-2],w[-1]
    cand=(adj[a]&adj[b]&adj[c])-{w[-4]}
    assert len(cand)==1
    w.append(cand.pop())
    if set(w[-4:])==set(start): close=step; break
seen=set(w)
ok["T4"]= close==30 and len(seen)==30
print("T4 600胞体の紐：最初の正四面体に戻る段",close,"／紐に現れる頂点",len(seen),"→","OK" if ok["T4"] else "NG")
print("第9章の検定 総合","OK" if all(ok.values()) else "NG")
# 図のデータ：3次元の紐（120段、表示のための小数は図の側で）と、600胞体の紐の30頂点の並び（a+bφ）
json.dump({"s3":[[x/3**p[1] for x in p[0]] for p in v[:64]],   # 表示用の小数（観察者側）
           
           "v4":[[[x.a,x.b] for x in V[i]] for i in w[:33]],
           "all4":[[[x.a,x.b] for x in q] for q in V]},open("ch09_tetra.json","w"),separators=(",",":"))
