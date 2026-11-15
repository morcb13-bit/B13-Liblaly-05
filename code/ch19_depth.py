"""第19章：精度は深さで ── φ 倍を重ねた深さと、近似のずれ
第18章の 26×26 の整数行列 M（番地ごとに (a,b)→(b,a+b)）を k 回掛けることを「深さ k」と呼ぶ。
φ の近似は分数 F_{k+1}/F_k の形をとらず、ずれ e_k ＝ F_{k+1} − F_k·φ を a＋bφ の整数のまま持つ（割り算なし）。
基準（走らせる前に固定。結果を見てから変えない）
  19-1 深さとフィボナッチ：k＝1〜60 で、M^k の番地ごとの区画が [[F_{k−1}, F_k], [F_k, F_{k+1}]]、区画の外は 0。
  19-2 誤差の正体：k＝0〜60 で、φ^k·e_k ＝ (−1)^k が a＋bφ の整数として成り立つ。
  19-3 単調性：k＝0〜59 で |e_{k+1}| ＜ |e_k| が整数の比較で成り立ち、e_k の符号が (−1)^k で交互に変わる。
  19-4 必要な深さ：13^m·|e_k| ＜ 1 となる最小の k を m＝1〜4 で求める。どの m でも k≤60 で見つかり、m とともに減らない。
  負の対照（記録） (a,b)→(b, a＋2b) で作った整数の組 g では、g_{k+1} − g_k·φ の絶対値が深さとともに小さくならない。
"""
from zphi import Z
P=Z(0,1); N=13; D=26; K=60
def matmul(A,B): return [[sum(A[i][t]*B[t][j] for t in range(D) if A[i][t]) for j in range(D)] for i in range(D)]
M=[[0]*D for _ in range(D)]
for i in range(N): M[2*i][2*i+1]=1; M[2*i+1][2*i]=1; M[2*i+1][2*i+1]=1
F=[0,1]
while len(F)<K+3: F.append(F[-1]+F[-2])
# 19-1
Mk=[r[:] for r in M]; ok1=True
for k in range(1,K+1):
    if k>1: Mk=matmul(Mk,M)
    for i in range(D):
        for j in range(D):
            bi,bj=i//2,j//2
            if bi!=bj: want=0
            else: want=[[F[k-1],F[k]],[F[k],F[k+1]]][i%2][j%2]
            if Mk[i][j]!=want: ok1=False
print("19-1 M^k の番地ごとの区画がフィボナッチ（k=1〜60）",ok1,"／k=60 の区画",[[F[59],F[60]],[F[60],F[61]]],"→","OK" if ok1 else "NG")
# 19-2
def phipow(k):
    z=Z(1)
    for _ in range(k): z=z*P
    return z
e=[Z(F[k+1],-F[k]) for k in range(K+1)]
ok2=all(phipow(k)*e[k]==Z((-1)**k) for k in range(K+1))
print("19-2 φ^k·e_k ＝ (−1)^k（k=0〜60）",ok2,"／e_10 ＝",repr(e[10]),"→","OK" if ok2 else "NG")
# 19-3
def ab(z): return z if z.sign()>=0 else Z(0)-z
mono=all(ab(e[k+1])<ab(e[k]) for k in range(K))
alt=all(e[k].sign()==(1 if k%2==0 else -1) for k in range(K+1))
ok3=mono and alt
print("19-3 |e_{k+1}| ＜ |e_k|（k=0〜59）",mono,"／符号が交互",alt,"→","OK" if ok3 else "NG")
# 19-4
need={}
for m in range(1,5):
    Nm=13**m
    need[m]=next((k for k in range(K+1) if ab(e[k])*Nm<Z(1)),None)
ok4=all(need[m] is not None for m in need) and all(need[m]<=need[m+1] for m in range(1,4))
print("19-4 13^m·|e_k| ＜ 1 となる最小の深さ k",need,"→","OK" if ok4 else "NG")
print("   そのときの整数の組 (F_k, F_{k+1})",{m:(F[need[m]],F[need[m]+1]) for m in need})
# 負の対照（記録）
g=[0,1]
while len(g)<K+3: g.append(g[-2]+2*g[-1])
eg=[Z(g[k+1],-g[k]) for k in range(K+1)]
shrink=sum(1 for k in range(1,K) if ab(eg[k+1])<ab(eg[k]))
print("負の対照（記録） (a,b)→(b,a+2b) の組で |g_{k+1}−g_k·φ| が小さくなった段の数（k=1〜59）",shrink,"／例 k=10:",repr(eg[10]))
print("第19章の検定 総合（19-1〜19-4）","OK" if (ok1 and ok2 and ok3 and ok4) else "NG")
