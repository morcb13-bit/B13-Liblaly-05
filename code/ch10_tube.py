"""第10章：管は4次元で閉じる ── ペンローズの階段
3次元：30面体10個の輪を、5回軸のまわりに回しながら一つずつ積む（helix5x10.py）。
4次元：600胞体の大円の十角形のまわりの頂点を、段に分ける（tube4d.py）。判定はすべて a+bφ の整数。
基準（走らせる前に決めた）
  U1 3次元：隣の輪と30面体を2個共有し、3周（輪15個）積んでも30面体どうしが食い込まない積み方がある。
     そのような積み方はどれも、輪5個で軸方向の平行移動（ずれ ≠ 0）だけに戻る ── 上り続けて閉じない
  U2 4次元：600胞体の大円の十角形 D（10頂点）と、D のどの頂点とも直交する頂点10個 D′。残り100個は D 側50・D′ 側50。
     D 側の50個は十角形の辺の真ん中の向きに5個ずつ10段で、どの段も正五角形。
     どの頂点も下の段と2本・同じ段と2本・上の段と2本の辺でつながり、10段目の上は1段目 ── 管は閉じる
  負の対照：大円ではない10個（正二十面体の両端を除いたジグザグ）には、どれとも直交する頂点がない（管に分けられない）
"""
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    import helix5x10 as H
    import tube4d as T
from zphi import Z
ok={}
good=[o for o in H.out if o["share"]==2 and o["ov"]==0]
pitches=sorted({str(o["pitch"]) for o in good})
ok["U1"]= len(good)>0 and all(o["pitch"]!=Z(0) for o in good)
print("U1 3次元：隣と2個共有・3周で食い込み0の積み方",len(good),"通り／輪5個で上がる高さ（×|A|²）",pitches,"／高さ0（閉じる）もの",sum(1 for o in good if o["pitch"]==Z(0)),"→","OK" if ok["U1"] else "NG")
print("   その積み方の30面体の数（3周）",sorted({o["n"] for o in good}))
lay=[len(T.layers.get(k,[])) for k in range(10)]
ok["U2"]= len(T.D)==10 and len(T.Dp)==10 and (len(T.near),len(T.far),T.tie)==(50,50,0) and lay==[5]*10 and T.okp and set(T.deg)=={(2,2,2)}
print("U2 4次元：大円の十角形",len(T.D),"／直交する頂点",len(T.Dp),"／残りの分かれ方",len(T.near),len(T.far),"／段ごとの頂点",lay,"／どの段も正五角形",T.okp,"／（下,同じ,上）の辺",dict(T.deg),"→","OK" if ok["U2"] else "NG")
ok["対照"]= len(T.ortho)==0
print("負の対照 正二十面体のジグザグ10個：どれとも直交する頂点",len(T.ortho),"→","OK" if ok["対照"] else "NG")
print("第10章の検定 総合","OK" if all(ok.values()) else "NG")
