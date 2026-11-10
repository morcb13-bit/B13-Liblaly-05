"""第8章：円環で埋めるオートマトン ── 重ならず、5回対称を保ち、φ² で戻る
道具は ring_automaton.py（11月7日）。重なりの判定と、候補を近い順に並べるところを、整数だけの判定に置き換えて走らせる。
  Z[ζ10]（ζ=e^{iπ/5}、基底 1,ζ,ζ²,ζ³）の数 z について
    虚部の符号：Im z / sin36° ＝ a1 + (a2+a3)φ
    実部の2倍：2Re z ＝ (2a0 − a2 + a3) + (a1 + a2 − a3)φ
規則（11月7日の監督の指定どおり）
  既にある円環と五角形2枚を共有し、ほかの円環とは0枚か2枚、五角形どうしは重ならない。
  5回回転で写る5個を一組で置く。中心に近い組から。
基準（走らせる前に決めた）
  A1 60歩（組を60回）置くまで止まらない。どの歩でも、円環の並びは5回回転で自分に写り、
     五角形どうしの重なりは0（整数の判定で全部の組を調べ直す）
  A2 負の対照：重なりの判定を外して同じ規則で置くと、重なる五角形の組が出る
  構造として記録：一枚の五角形が属する円環の数ごとの枚数
  A3・A4（相似）は ch08_retest.py（改訂基準：鏡映と回転を許す）。
  初版（ch08_automaton_v1.py）は A1 の判定に基準にない条件「属する円環は1つか2つ」を入れていた。
  初版の A3（鏡映・回転なし）・A4（回転のみ）の NG は ch08_automaton_v1_out.txt に残す。
"""
import sys, time, io, contextlib, functools, itertools
from collections import Counter
from zphi import Z
t0=time.time()
with contextlib.redirect_stdout(io.StringIO()):
    import ring_automaton as RA
Q=RA.Q
def imq(z): return Z(z[1], z[2]+z[3])                      # Im/sin36
def re2(z): return Z(2*z[0]-z[2]+z[3], z[1]+z[2]-z[3])      # 2Re
def cross(u,v): return imq(Q.zmul(Q.conj(u),v))
def n2(z): return re2(Q.zmul(z,Q.conj(z)))                 # 2|z|²
def edges(P):     # 正五角形の辺：長さ²が最小の5組
    d=[(n2(Q.zsub(a,b)),i,j) for (i,a),(j,b) in itertools.combinations(list(enumerate(P)),2)]
    m=None
    for x,_,_ in d:
        if m is None or (x-m).sign()<0: m=x
    return [(P[i],P[j]) for x,i,j in d if x==m]
def separated_by(P,Qp):
    for a,b in edges(P):
        e=Q.zsub(b,a)
        inside=[cross(e,Q.zsub(x,a)).sign() for x in P if x!=a and x!=b]
        s=inside[0]
        if all(cross(e,Q.zsub(y,a)).sign()*s<=0 for y in Qp): return True
    return False
def overlap(P,Qp): return not(separated_by(P,Qp) or separated_by(Qp,P))
def verts(c,vs): return [Q.zadd(c,v) for v in vs]
def can_add(self,c50,check=True):
    new=RA.ring_pents(c50); share=Counter()
    for p,vs in new.items():
        if p in self.pent:
            for ri in self.pent[p][1]: share[ri]+=1
            continue
        if not check: continue
        A=verts(p,vs)
        for q in self.near(p):
            if q==p or q in new: continue
            if abs(Q.zf(q)-Q.zf(p))>=RA.SPAN: continue          # 遠すぎる組を省くだけ（2倍の外接円より遠い）
            if overlap(A,verts(q,self.pent[q][0])): return None
    if any(v not in (2,) for v in share.values()): return None if check else share
    return share
RA.World.can_add=can_add
_sorted=sorted
def exact_sorted(it,key=None):
    L=list(it)
    if key is not None and len(L) and isinstance(L[0],tuple) and len(L[0])==4:
        return _sorted(L,key=functools.cmp_to_key(lambda a,b:(n2(a)-n2(b)).sign() or ((a>b)-(a<b))))
    return _sorted(L,key=key)
RA.sorted=exact_sorted
STEPS=int(sys.argv[1]) if len(sys.argv)>1 else 60
def run(check,steps,log_sym=False):
    W=RA.start((0,2)); start=list(W.rings); syms=[]
    snaps=[list(W.rings)]
    for st in range(steps):
        n,_=RA.grow(W,10**9,1,check=check,log=False)
        if n<1: return W,start,snaps,syms,st
        rs=set(W.rings); syms.append(all(RA.rot(c,1) in rs for c in W.rings)); snaps.append(list(W.rings))
    return W,start,snaps,syms,steps
W,start,snaps,syms,done=run(True,STEPS)
print("出発：円環",len(start),"／%d歩：進んだ歩 %d・円環 %d・五角形 %d" % (STEPS,done,len(W.rings),len(W.pent)))
print("どの歩でも5回回転で自分に写る",all(syms),"（",len(syms),"歩）")
# 全部の五角形の組を整数で調べ直す
P=list(W.pent.items()); ov=0
for (p,(vs,_)),(q,(ws,_)) in itertools.combinations(P,2):
    if abs(Q.zf(p)-Q.zf(q))>=RA.SPAN: continue
    if overlap(verts(p,vs),verts(q,ws)): ov+=1
print("重なる五角形の組",ov)
print("構造：一枚の五角形が属する円環の数ごとの枚数",dict(sorted(Counter(len(rs) for _,(vs,rs) in W.pent.items()).items())))
a1= done==STEPS and all(syms) and ov==0
print("A1 止まらず・5回対称・重なり0","OK" if a1 else "NG")
# A2
Wc,_,_,_,dc=run(False,8)
Pc=list(Wc.pent.items()); ovc=0
for (p,(vs,_)),(q,(ws,_)) in itertools.combinations(Pc,2):
    if abs(Q.zf(p)-Q.zf(q))>=RA.SPAN: continue
    if overlap(verts(p,vs),verts(q,ws)): ovc+=1
print("負の対照 重なりの判定を外して8歩：円環",len(Wc.rings),"／重なる五角形の組",ovc)
a2= ovc>0
print("A2","OK" if a2 else "NG")
R=max(abs(Q.zf(c)) for c in W.rings)/50
print("いちばん外の円環の半径（表示用の小数）%.1f ／ 時間 %.0f 秒" % (R,time.time()-t0))
print("第8章の検定 A1・A2","OK" if (a1 and a2) else "NG")
import json
json.dump({"rings":[list(c) for c in W.rings],"start":len(start),
           "pents":[[list(p),1 if vs is RA.VLO else 0,rs] for p,(vs,rs) in W.pent.items()],
           "vlo":[list(v) for v in RA.VLO],"vhi":[list(v) for v in RA.VHI],
           "steps":[len(s) for s in snaps]},open("ch08_ra.json","w"))
