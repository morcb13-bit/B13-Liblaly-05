"""第2章：6次元から立体へ ── ひし形30面体の窓と、2種類の黄金のひし面体
数は a+bφ の整数の組（zphi.py）。道具は cp2_tiling.py（11月7日）をそのまま使う。
基準（走らせる前に決めた）
  S0 窓：6次元の単位立方体の頂点64個を残渣の側に写すと、外形は面30・頂点32（5本が集まる頂点12、3本が集まる頂点20）の
     ひし形30面体。どの面にもちょうど4点が乗る
  S1 窓に入る点だけを採り、ひし面体を選ぶと、内側の面はどれもちょうど2つのひし面体に属する
  S2 ひし面体の体積は2種類で、比は φ
  負の対照 窓を 13/10 倍にすると、S1 が崩れる（3つ以上に属する面が出る）
"""
import sys, itertools, json
from zphi import Z
import cp2_tiling as C
R2=int(sys.argv[1]) if len(sys.argv)>1 else 200
R2in=int(sys.argv[2]) if len(sys.argv)>2 else 120
# S0 窓の形
pts=[tuple(C.lin(C.Gp,s)) for s in itertools.product((1,-1),repeat=6)]   # 2倍の座標
on={}
faces=0; ok_face=True
for m,w in zip(C.M,C.W2):
    for sg in (1,-1):
        F=[p for p in pts if C.dot(p,m)*sg==w]
        faces+=1; ok_face&= len(set(F))==4
        for p in set(F): on[p]=on.get(p,0)+1
deg={}
for p,c in on.items(): deg[c]=deg.get(c,0)+1
print("窓：面",faces,"／頂点",len(on),"／頂点に集まる面の数ごとの頂点数",dict(sorted(deg.items())),"／どの面にも4点",ok_face)
s0= faces==30 and len(on)==32 and deg=={3:20,5:12} and ok_face
# 参考：同じ立方体を立体の側に写したときの外形の頂点は、序章のどの殻か
from collections import Counter
pp=[tuple(C.lin(C.GEN,s)) for s in itertools.product((1,-1),repeat=6)]
Mp=[C.cross(C.GEN[i],C.GEN[j]) for i,j in itertools.combinations(range(6),2)]
Wp=[sum((abs(C.dot(G,m)) for G in C.GEN),Z()) for m in Mp]
hull={p for m,w in zip(Mp,Wp) for sg in (1,-1) for p in pp if C.dot(p,m)*sg==w}
print("参考 立体の側の外形の頂点",len(hull),"殻（長さ²）ごと",dict(Counter(str(C.dot(p,p)) for p in hull)))
print("S0 窓はひし形30面体","OK" if s0 else "NG")
pts_t,tiles=C.build(Z(R2))
inner,vols=C.t0(pts_t,tiles,Z(R2in),"窓で選ぶ")
ks=list(vols); s1= set(inner)=={2}
s2= len(ks)==2 and (ks[0]==ks[1]*C.P or ks[1]==ks[0]*C.P)
print("S1 内側の面はすべて2つのひし面体に属する","OK" if s1 else "NG")
print("S2 体積は2種類で比 φ","OK" if s2 else "NG")
pc,tc=C.build(Z(R2),rule="all8",scale=(13,10))
ic,_=C.t0(pc,tc,Z(R2in),"負の対照：窓を 13/10 倍")
ctl= set(ic)!={2}
print("負の対照（崩れる）","OK" if ctl else "NG")
print("第2章の検定 総合","OK" if (s0 and s1 and s2 and ctl) else "NG")
# 図のデータ：ひし面体（中心が近い順）と、採った点の残渣の側の座標
def zz(v): return [[c.a,c.b] for c in v]
T=[]
for n,Tr in tiles:
    c2=tuple(2*n[k]+(1 if k in Tr else 0) for k in range(6)); x=C.phys(c2)
    T.append((C.dot(x,x),n,Tr))
T.sort(key=lambda t:(t[0].a+1.6180339887*t[0].b))
T=[t for t in T if (t[0]-Z(4*60)).sign()<0]
P=sorted(pts_t,key=lambda n:(lambda x:(C.dot(x,x).a+1.6180339887*C.dot(x,x).b))(C.phys(n)))
thick=lambda Tr: abs(C.det3(*[C.GEN[k] for k in Tr]))==Z(2,2)
S=set(pts_t); rej=set()
for n in P:
    if (C.dot(C.phys(n),C.phys(n))-Z(60)).sign()>=0: continue
    for k in range(6):
        for sg in (1,-1):
            m=list(n); m[k]+=sg; m=tuple(m)
            if m not in S: rej.add(m)
json.dump({"tiles":[[list(n),list(Tr),int(thick(Tr))] for _,n,Tr in T],"rej":[[list(m),zz(C.phys(m)),zz(C.perp(m))] for m in rej],
           "pts":[[list(n),zz(C.phys(n)),zz(C.perp(n))] for n in P],
           "gen":[zz(g) for g in C.GEN],"genp":[zz(g) for g in C.Gp],
           "win":[zz(p) for p in on]},open("ch02_solid.json","w"),separators=(",",":"))
print("図のデータ：ひし面体",len(T),"捨てた隣の点",len(rej))
