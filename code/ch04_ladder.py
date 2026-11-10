"""第4章：窓の梯子 ── 5・6・8、その先
数はすべて a+bφ の整数（zphi.py）。第1〜3章の窓をそのまま並べる。
基準（走らせる前に決めた）
  L1 窓の頂点の数：5次元→平面の窓（指数1）は 5、6次元→立体の窓は 32、8次元→4次元の窓は 120
  L2 窓を保つ置き換えの数：
     5次元：5本の向きの並べ替えで、向きどうしの角（内積）を保つもの 10
     6次元：6本の軸の符号つき並べ替えで、軸どうしの内積を保つもの 120
     8次元：イコシアンの左右からの掛け算 q→lqr と q→l q̄ r（l,r は単位イコシアン）で、異なる置き換え 14400
  L3 窓の入れ子：
     (a) 600胞体を頂点の向きに段で切ると、点の数は 1,12,20,12,30,12,20,12,1。
         20点の段と、その一つ内側の12点の段を、切る向きに沿って3次元に落とすと、
         32点の距離の並びが、ひし形30面体の頂点32個と同じ比になる（6次元の窓が現れる）
     (b) ひし形30面体を5回軸の向きに段で切ると、5点が正五角形に並ぶ段があり、
         その中に半径²の比が φ² の組がある（第1章の窓 V1・V2 の比）
  負の対照：L3(a) で、20点の段と一番外側の12点の段を組むと、同じ比にならない
"""
import itertools
from zphi import Z
P=Z(0,1); ONE=Z(1); ZERO=Z(0)
def key(z): return z.a+1.6180339887*z.b        # 並べる順番を決めるだけ（判定には使わない）
ok=True
# ---------- L1 ----------
# 5次元：指数1の切り口 = ζ^(2j) の5点（第1章 W0 と同じ）。凸包の頂点は5
n5=5
# 6次元：第2章 S0（窓の頂点32）
import cp2_tiling as C
pts=[tuple(C.lin(C.Gp,s)) for s in itertools.product((1,-1),repeat=6)]
on=set()
for m,w in zip(C.M,C.W2):
    for sg in (1,-1):
        for p in pts:
            if C.dot(p,m)*sg==w: on.add(p)
n6=len(on)
# 8次元：単位イコシアン（第3章）
src=open("ch03_e8.py").read().split("# ---- H0")[0]
g={}; exec(src.split('"""',2)[2],g)
U=g["U"]; n8=len(U)
print("L1 窓の頂点の数",n5,n6,n8)
l1=(n5,n6,n8)==(5,32,120); ok&=l1
# ---------- L2 ----------
cls=lambda d:min(d%5,(-d)%5)          # 向きの差の種類（0・1・2）で内積が決まる
c5=sum(1 for s in itertools.permutations(range(5)) if all(cls(s[i]-s[j])==cls(i-j) for i in range(5) for j in range(5)))
GEN=C.GEN
Gm=[[C.dot(GEN[i],GEN[j]) for j in range(6)] for i in range(6)]
c6=0
for s in itertools.permutations(range(6)):
    for sg in itertools.product((1,-1),repeat=6):
        if all(Gm[s[i]][s[j]]*(sg[i]*sg[j])==Gm[i][j] for i in range(6) for j in range(i,6)): c6+=1
def qm(p,q):   # 四元数の積（2倍の座標どうしなので、結果は4倍の座標）
    a1,b1,c1,d1=p; a2,b2,c2,d2=q
    return (a1*a2-b1*b2-c1*c2-d1*d2, a1*b2+b1*a2+c1*d2-d1*c2, a1*c2-b1*d2+c1*a2+d1*b2, a1*d2+b1*c2-c1*b2+d1*a2)
def qc(p): return (p[0],-p[1],-p[2],-p[3])
E4=[(Z(2),ZERO,ZERO,ZERO),(ZERO,Z(2),ZERO,ZERO),(ZERO,ZERO,Z(2),ZERO),(ZERO,ZERO,ZERO,Z(2))]
maps=set()
for l in U:
    L=[qm(l,e) for e in E4]; Lc=[qm(l,qc(e)) for e in E4]
    for r in U:
        maps.add(("q",)+tuple(qm(x,r) for x in L))
        maps.add(("c",)+tuple(qm(x,r) for x in Lc))
c8=len({m[1:] for m in maps})
print("L2 窓を保つ置き換えの数",c5,c6,c8)
l2=(c5,c6,c8)==(10,120,14400); ok&=l2
# ---------- L3(a) 600胞体の段 ----------
top=(Z(2),ZERO,ZERO,ZERO)
def dot4(p,q):
    s=Z()
    for a,b in zip(p,q): s=s+a*b
    return s
lev={}
for u in U: lev.setdefault(dot4(u,top),[]).append(u)   # 4×(高さ)
ks=sorted(lev,key=key,reverse=True)
print("L3(a) 600胞体の段（高さ×4 と点の数）",[(str(k),len(lev[k])) for k in ks])
cnt=[len(lev[k]) for k in ks]
def proj(u): return u[1:]                      # 切る向き（第0成分）に沿って落とす
def d2(p,q):
    s=Z()
    for a,b in zip(p,q): s=s+(a-b)*(a-b)
    return s
def dist_list(S):
    return sorted((d2(a,b) for a,b in itertools.combinations(S,2)),key=key)
def same_ratio(A,B):          # 並びが同じ比か：A[i]·B[0] == B[i]·A[0]（掛け算だけで比べる）
    if len(A)!=len(B): return False
    return all(a*B[0]==b*A[0] for a,b in zip(A,B))
RT=dist_list(list(on))
d20=[k for k in ks if len(lev[k])==20 and k.sign()>0][0]
i12=[k for k in ks if len(lev[k])==12 and k.sign()>0]
inner12=min(i12,key=key); outer12=max(i12,key=key)
S_in=[proj(u) for u in lev[d20]+lev[inner12]]
S_out=[proj(u) for u in lev[d20]+lev[outer12]]
a_ok=same_ratio(dist_list(S_in),RT); a_ctl=same_ratio(dist_list(S_out),RT)
print("  20点の段＋内側の12点の段 → ひし形30面体と同じ比",a_ok)
print("  負の対照 20点の段＋外側の12点の段 → 同じ比",a_ctl)
l3a= cnt==[1,12,20,12,30,12,20,12,1] and a_ok and not a_ctl; ok&=l3a
# ---------- L3(b) ひし形30面体の段 ----------
ax=None
for p in on:
    if sum(1 for m,w in zip(C.M,C.W2) for sg in (1,-1) if C.dot(p,m)*sg==w)==5: ax=p; break
lv={}
for p in on: lv.setdefault(C.dot(p,ax),[]).append(p)
ks2=sorted(lv,key=key,reverse=True)
print("L3(b) ひし形30面体の段（5回軸、点の数）",[len(lv[k]) for k in ks2])
A2=C.dot(ax,ax)
pent=[]
for k in ks2:
    S=lv[k]
    if len(S)!=5: continue
    # 軸からの距離²×|ax|²：|p|²|ax|² − (p·ax)²
    r=[C.dot(p,p)*A2 - k*k for p in S]
    side=sorted((d2(a,b) for a,b in itertools.combinations(S,2)),key=key)
    regular= len(set(r))==1 and side[:5]==[side[0]]*5 and side[5:]==[side[5]]*5 and side[5]==side[0]*P*P
    pent.append((k,r[0],regular))
print("  5点の段",len(pent),"／正五角形",sum(1 for x in pent if x[2]))
rs=[x[1] for x in pent if x[2]]
phi2=any(a==b*P*P for a in rs for b in rs)
print("  半径²の比が φ² の組",phi2)
l3b= len(pent)>0 and all(x[2] for x in pent) and phi2; ok&=l3b
print("L1",("OK" if l1 else "NG"),"L2",("OK" if l2 else "NG"),"L3(a)",("OK" if l3a else "NG"),"L3(b)",("OK" if l3b else "NG"))
print("第4章の検定 総合","OK" if ok else "NG")
