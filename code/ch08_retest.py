"""第8章 A3・A4 の再検定（改訂基準）
改訂理由：11月7日の I1 との比較可能性を確保するため（I1 は鏡映と回転を許していた）。
旧基準（回転のみ・A3 は回転も許さない）の結果は ch08_automaton_out.txt に残す：A3 NG、A4 NG。
改訂基準（走らせる前に決めた）
  共通：原点中心のまま、鏡映（複素共役）と回転 ζ^k（k=0,…,9）を許す。拡大率 φ・φ²・φ³ を別々に調べる。
  A3 出発の形（円環10個）：φ² 倍で入り、φ 倍・φ³ 倍では入らない
  A4 1歩目まで進めた形（円環15個）：φ² 倍で入る
データは ch08_automaton.py 60 の出力 ch08_ra.json（円環は置いた順に並んでいる）。"""
import json, io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    import q1_match as Q
D=json.load(open("ch08_ra.json")); C=[tuple(c) for c in D["rings"]]; Cs=set(C)
n0=D["start"]; shapes={"A3 出発の形":C[:n0],"A4 1歩目の形":C[:D["steps"][1]]}
PHI=(1,0,1,-1); P2=Q.zmul(PHI,PHI); P3=Q.zmul(P2,PHI)
SC={"φ":PHI,"φ²":P2,"φ³":P3}
def hits(shape,s):
    out=[]
    for refl in (False,True):
        for k in range(10):
            T=[Q.zmul(Q.zpow(k),Q.conj(Q.zmul(s,t)) if refl else Q.zmul(s,t)) for t in shape]
            if all(x in Cs for x in T): out.append(("鏡映＋" if refl else "")+"ζ^%d"%k)
    return out
res={}
for name,sh in shapes.items():
    for sn,s in SC.items():
        h=hits(sh,s); res[(name,sn)]=h
        print("%s（円環 %d）を %s 倍：入る向き %d 通り %s" % (name,len(sh),sn,len(h),h))
a3= bool(res[("A3 出発の形","φ²")]) and not res[("A3 出発の形","φ")] and not res[("A3 出発の形","φ³")]
a4= bool(res[("A4 1歩目の形","φ²")])
print("A3（改訂基準）","OK" if a3 else "NG")
print("A4（改訂基準）","OK" if a4 else "NG")
print("旧基準の結果（残す）：A3 NG（鏡映・回転なしで φ・φ²・φ³ いずれも入らない）、A4 NG（回転のみで入らない）")
