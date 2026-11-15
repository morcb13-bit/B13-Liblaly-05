"""第17章：13次元までの図形一覧 ── 正多胞体を全部数える
改訂の記録：第4部の範囲を11次元から13次元に改めたため、同じ基準で n＝13 まで延ばした（旧版は ch17_list_v1.py／_v1_out.txt）。
  記号のたどり方は、途中の D が正でなくなった記号の先を打ち切る形に変えた（数える対象は同じ）。
基準（走らせる前に決めた）
  L1 n 次元の正多胞体を、シュレーフリ記号 {p1,…,p(n−1)}（p は 3,4,5,6）で全部たどり、
     「行列 2I − (隣の c) の先頭からの小行列式がすべて正」になるものを数える。
     c² ＝ 4cos²(π/p) は p＝3,4,5,6 で 1, 2, φ²（＝1＋φ）, 3。符号は a＋bφ の整数で判定する。
     小行列式は D_k ＝ 2·D_(k−1) − c²·D_(k−2) の足し算と掛け算だけで出す。
     数は n＝3 で 5、n＝4 で 6、n＝5〜13 で 3。
     ［知識］p≥7 では c² が 3 より大きくなり、n≥3 では正にならない。ここでは p≤6 だけをたどる。
  L2 φ（p＝5）を含むものは、n＝3 で 2、n＝4 で 2、n＝5〜13 で 0。{3,4,3} は n＝4 だけ。
  L3 三つの族を整数の頂点で作り、面を全部数える（n＝1〜13）。
     正単体：頂点 e_0…e_n、面は頂点のどの組でもよい。立方体：±1 の列。正軸体：±e_i、面は向かい合う二つを含まない組。
     立方体と正軸体は、面の数を逆に並べると一致する。外側の面の交互の和は、どの族も 1−(−1)^n。
  記録 行列式が 0 になる（それまでの小行列式は正）もの＝平らに敷き詰める形 を数として出す。
  負の対照 c² の表で φ² を 3 に置き換える（5 を 6 と同じに扱う）と、n＝3 の数は 5 にならない。
"""
import itertools
from math import comb
from zphi import Z
P=Z(0,1)
C2={3:Z(1),4:Z(2),5:P*P,6:Z(3)}
def minors(sym,c2):
    D=[Z(1),Z(2)]
    for k in range(len(sym)):
        D.append(D[-1]*2-c2[sym[k]]*D[-2])
    return D[1:]
def kind(sym,c2):
    D=minors(sym,c2)
    if all(Z(0)<d for d in D): return "球"
    if all(Z(0)<d for d in D[:-1]) and D[-1]==Z(0): return "平"
    return None
def count(n,c2):
    # 記号を左から一つずつ伸ばす。途中の D が正でなくなったら、その先はどれも「球」にも「平」にもならない
    sph=[];flat=[]
    def rec(sym,D2,D1):
        if len(sym)==n-1:
            if Z(0)<D1: sph.append(tuple(sym))
            return
        for p in (3,4,5,6):
            D=D1*2-c2[p]*D2
            if Z(0)<D: rec(sym+[p],D1,D)
            elif D==Z(0) and len(sym)+1==n-1: flat.append(tuple(sym+[p]))
    rec([],Z(1),Z(2))
    return sorted(sph),sorted(flat)
S={};F={}
for n in range(3,14): S[n],F[n]=count(n,C2)
cnt={n:len(S[n]) for n in S}
l1= cnt=={3:5,4:6,**{n:3 for n in range(5,14)}}
print("L1 正多胞体の数",cnt,"→","OK" if l1 else "NG")
for n in (3,4,5,13): print(f"   {n}次元",["{"+",".join(map(str,s))+"}" for s in S[n]])
phi={n:sum(1 for s in S[n] if 5 in s) for n in S}
has343={n:(3,4,3) in S[n] for n in S}
l2= phi=={3:2,4:2,**{n:0 for n in range(5,14)}} and [n for n in has343 if has343[n]]==[4]
print("L2 φ（5）を含むもの",phi,"／{3,4,3} がある次元",[n for n in has343 if has343[n]],"→","OK" if l2 else "NG")
print("記録 平らに敷き詰める形",{n:["{"+",".join(map(str,s))+"}" for s in F[n]] for n in range(3,7)},"／n=7〜13 の数",{n:len(F[n]) for n in range(7,14)})
# L3
def f_simplex(n):
    f=[0]*(n+1)
    for m in range(1,1<<(n+1)):
        f[bin(m).count("1")-1]+=1
    return f
def f_cube(n):
    f=[0]*(n+1)
    for pat in itertools.product((1,-1,0),repeat=n):   # 0 の座標が自由に動く面
        f[pat.count(0)]+=1
    return f
def f_cross(n):
    f=[0]*(n+1)
    for pat in itertools.product((1,-1,0),repeat=n):   # 選んだ頂点：+e_i か −e_i か選ばないか
        k=n-pat.count(0)
        if 1<=k<=n: f[k-1]+=1
    f[n]=1
    return f
ok3=True; rows={}
for n in range(1,14):
    a,b,c=f_simplex(n),f_cube(n),f_cross(n)
    rows[n]=(a,b,c)
    if b[:n]!=c[:n][::-1]: ok3=False
    for f in (a,b,c):
        if sum((-1)**k*f[k] for k in range(n))!=1-(-1)**n: ok3=False
    if a!=[comb(n+1,k+1) for k in range(n+1)] or b!=[comb(n,k)*2**(n-k) for k in range(n+1)]: ok3=False
print("L3 三つの族の面の数（交互の和・立方体と正軸体の逆並び）",ok3,"→","OK" if ok3 else "NG")
for n in (4,7,13):
    a,b,c=rows[n]; print(f"   {n}次元 正単体",a); print(f"   {n}次元 立方体",b); print(f"   {n}次元 正軸体",c)
bad=dict(C2); bad[5]=Z(3)
neg= len(count(3,bad)[0])!=5
print("負の対照 φ² を 3 に置き換えたときの3次元の数",len(count(3,bad)[0]),"→","OK" if neg else "NG")
print("第17章の検定 総合","OK" if (l1 and l2 and ok3 and neg) else "NG")
