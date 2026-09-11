#!/usr/bin/env python3
"""Preserve the first research round independently of the earlier 27-model audit."""
import csv, hashlib, json, math
from pathlib import Path
from collections import defaultdict
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
OUT=HERE/'search_results'
files=[OUT/'round1_9112026.csv',OUT/'round1_9112027.csv']
rows=[]
for p in files:
    log=p.with_suffix('.log').read_text()
    assert 'PASS rows=1599' in log and 'FAIL:' not in log
    with p.open(newline='') as f: data=list(csv.DictReader(f))
    assert len(data)==1599
    rows.extend(data)
positive=[r for r in rows if r['kind'] in ('search_fields','search_sweep')]
assert len(positive)==3190
assert all(abs(complex(float(r['sign_real']),float(r['sign_imag']))-1)<2e-7 for r in positive)
maxerr=max(float(r['relative_error']) for r in rows)
maxphase=max(abs(complex(float(r['sign_real']),float(r['sign_imag']))-1) for r in positive)
groups=defaultdict(list)
for r in positive:groups[r['model']].append(r)
ctest=(OUT/'ctest_round1.log').read_text()
assert '100% tests passed out of 58' in ctest
registry={
 'scope':'New Hamiltonian constructions found in this goal, excluding obvious duplicates; old catalog and size/parameter variants do not count.',
 'stop_when_verified_count_exceeds':20,
 'verified_count':3,
 'goal_complete':False,
 'results':[
  {'id':'R001','name':'Lorentz local orbital family','status':'verified','counts_as':1,'model_section':'12/R001','test_models':['lorentz4','lorentz5','lorentz8','lorentz4_lattice'],'known_overlap':'m=2 repeats original spin1; tensor d=2 repeats m=3','novelty':'new in this catalog; spectral class IV; no literature-priority claim'},
  {'id':'R002','name':'Conjugate-tensor local orbital family d>=3','status':'verified','counts_as':1,'model_section':'12/R002','test_models':['conjugate_tensor2','conjugate_tensor3','conjugate_tensor3_lattice'],'known_overlap':'d=2 is R001 and counts only as an equivalence check','novelty':'new in this catalog; spectral class IV; no literature-priority claim'},
  {'id':'R003','name':'Symmetric-power TN orbital family d>=3,q>=2','status':'verified','counts_as':1,'model_section':'12/R003','test_models':['symmetric_power3_2','symmetric_power3_3','symmetric_power4_2','symmetric_power3_2_lattice'],'known_overlap':'d=2 spin representations; q=1 open hopping; further Young labels require merging','novelty':'new representation instance of known total-nonnegative semigroup mechanism'}],
 'pending':['quaternion Hermitian (n=2 overlaps Lorentz)','quaternion anti-Hermitian','positive superselection sectors','local Gaussian exponential vertices','dynamic Ising-fermion model'],
 'evidence':['search_results/round1_9112026.csv','search_results/round1_9112027.csv','search_results/ctest_round1.log'],
}
(HERE/'search_registry.json').write_text(json.dumps(registry,ensure_ascii=False,indent=2)+'\n')
lines=['# 持续搜索数值记录','',
 '本页对应 MODELS.md §12 的首轮新增构造，独立于前面的旧目录报告。',
 '**正式计数：3 个新增 Hamiltonian 构造。停止阈值：超过 20 个；当前目标仍进行中。**','',
 '## 本轮验算','',
 '- WSL Ubuntu、Intel oneAPI/MKL、单线程；两个种子 9112026 / 9112027。',
 '- 11 个单胞／阵列预设，每种子每预设 128 个独立实场配置和 16 次真实 PFQMC 左右交替扫描（含初态17个状态）。',
 '- 共 3190 条应为正的配置全部通过 phase=1 检查；另有 4 条预期负号、2 条一般复相位、2 条精确零点控制。',
 f'- 全记录最大权重相对误差 `{maxerr:.6g}`；正配置最大 `|sign-1|={maxphase:.6g}`。',
 '- 完整 CTest 58/58 通过，包含之前 46 项及本轮 12 项。',
 '- 数值样本不替代全参数证明；每个模型的证明、局域性及允许场限制直接写在 [MODELS.md](MODELS.md) §12。','',
 '| 验证预设 | 复模数 | 独立场配置 | 扫描状态 | 最大相对误差 |',
 '| --- | ---: | ---: | ---: | ---: |']
for name,rs in sorted(groups.items()):
    lines.append(f"| `{name}` | {rs[0]['modes']} | {sum(r['kind']=='search_fields' for r in rs)} | {sum(r['kind']=='search_sweep' for r in rs)} | {max(float(r['relative_error']) for r in rs):.3g} |")
lines+=['','## 实际检查覆盖','',
 '- 权重用现有 PFAPACK 的 signed Pfaffian；逐配置与中心化行列式和辅助小矩阵谱公式交叉检查。至多6复模的独立场记录还比较完整 Fock 指数。',
 '- Walker 使用已有 PfQMC 稳定化和真实 Metropolis 更新；逐扫描核对 tracked sign、getSignRaw 和独立 Green 矩阵，并检查有实际提议与接受。',
 '- Tensor d=2 与 Lorentz4 用 Pauli 基逐生成元验证酉等价，不计第二个结果。',
 '- 共同 Wei 排除使用实际允许方向的简单实 Majorana 谱证书；推导针对指定 HS 表示，不是对全部可能 HS 的排除。',
 '- 独立叶势破坏 Lorentz 正性：精确 -2800/81。',
 '- 放开 tensor 共轭相关性：四模权重约 20.473942-1.254887i。',
 '- 放开 symmetric-power 的相关Cartan限制：六模权重约 -7.14635870944e10；高精度参考、完整Fock与signed Pfaffian一致。',
 '- Lorentz/tensor 真正零点单独记nan，不被算成sign1。','',
 '## 复现','',
 '按 [CONTEXT.md](../CONTEXT.md) 激活WSL oneAPI后：','',
 '```bash',
 'build/python-tools/cmake/data/bin/cmake -S . -B build/context-check',
 'build/python-tools/cmake/data/bin/cmake --build build/context-check --target sign_free_search -j 2',
 'build/context-check/sign_free_search --model all --samples 128 --sweeps 16 --seed 9112026 --csv build/search_9112026.csv',
 'build/context-check/sign_free_search --model all --samples 128 --sweeps 16 --seed 9112027 --csv build/search_9112027.csv',
 '```','',
 '也可将 `--model all` 改成上表某一预设或 `controls`。默认 beta=.8、L_tau=3；独立场标准差轮换1.5、.5、1。',
 '源代码：[search_models.h](../inc/catalog/search_models.h)、[catalog_search.cpp](../src/catalog_search.cpp)、[search.cpp](search.cpp)。',
 '局域阵列的胞粒子数仍守恒，没有声称实现跨胞单费米子巡游。',
 '这是验证器与通用dense walker，不据短扫描报告物理热化、自相关或性能提升。','',
 '- [种子9112026 CSV](search_results/round1_9112026.csv) / [日志](search_results/round1_9112026.log)',
 '- [种子9112027 CSV](search_results/round1_9112027.csv) / [日志](search_results/round1_9112027.log)',
 '- [CTest](search_results/ctest_round1.log) / [数据SHA-256](search_results/round1_data.sha256) / [源代码快照SHA-256](search_results/round1_sources.sha256)',
 '- [计数与合并记录](search_registry.json) / [后续研究队列](SEARCH_QUEUE.md)','',
 '本轮是进展，不满足“超过20后停止并最终合并”的完成条件。后续构造必须保留同样的局域性、正号数值和重复性检查。']
(HERE/'SEARCH_RESULTS.md').write_text('\n'.join(lines)+'\n')
paths=[*files,*(p.with_suffix('.log') for p in files),OUT/'ctest_round1.log']
(OUT/'round1_data.sha256').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.name+'\n' for p in paths))
source=[ROOT/'CMakeLists.txt',ROOT/'src/catalog_search.cpp',ROOT/'src/catalog_core.cpp',ROOT/'src/pfqmc.cpp',ROOT/'src/skewMatUtils.cpp',ROOT/'inc/catalog/catalog.h',ROOT/'inc/catalog/search_models.h',ROOT/'inc/qr_udt.h',HERE/'search.cpp',HERE/'checks.cpp',HERE/'summarize_search_round1.py']
(OUT/'round1_sources.sha256').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(ROOT))+'\n' for p in source))
print('Research round1: 3 verified constructions; 3190 positive signed records; 58 tests passed.')
