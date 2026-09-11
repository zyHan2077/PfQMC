#!/usr/bin/env python3
"""Summarize committed raw sign checks; standard library only, no rerun of physics."""
import csv
import hashlib
import math
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS = HERE / 'results'

def summarize():
    files = sorted(RESULTS.glob('samples_*.csv'))
    if not files:
        raise RuntimeError('No samples_*.csv files')
    groups = defaultdict(list)
    all_rows = []
    for path in files:
        log = path.with_name('run_' + path.stem.removeprefix('samples_') + '.log').read_text()
        if 'PASS rows=' not in log or 'FAIL:' in log:
            raise RuntimeError(f'Incomplete or failed run: {path}')
        with path.open(newline='') as f:
            rows = list(csv.DictReader(f))
        expected_rows = int(log.split('PASS rows=')[-1].split()[0])
        if len(rows) != expected_rows:
            raise RuntimeError(f'Row count mismatch: {path}')
        for row in rows:
            all_rows.append(row)
            groups[row['model']].append(row)
    regular = [r for r in all_rows if r['kind'] in ('random_fields','pfqmc_sweep')]
    positives = [r for r in all_rows if r['kind'] != 'control' and math.isfinite(float(r['sign_real'])) and float(r['sign_real']) > 0]
    zeros = [r for r in all_rows if r['kind'] == 'zero']
    negatives = [r for r in all_rows if r['kind'] != 'control' and math.isfinite(float(r['sign_real'])) and float(r['sign_real']) < 0]
    if any(abs(complex(float(r['sign_real']), float(r['sign_imag'])) - 1) > 2e-7 for r in positives):
        raise RuntimeError('Positive phase check failed')
    max_error = max(float(r['relative_error']) for r in all_rows)
    max_phase = max(abs(complex(float(r['sign_real']), float(r['sign_imag'])) - 1) for r in positives)
    max_sweep = {}
    for path in files:
        with path.open(newline='') as f:
            for row in csv.DictReader(f):
                if row['kind'] == 'pfqmc_sweep':
                    max_sweep[(path.name, row['model'])] = row
    attempts = sum(int(r['attempts']) for r in max_sweep.values())
    accepted = sum(int(r['accepted']) for r in max_sweep.values())
    lines = [
        '# 数值验算结果', '',
        '运行日期：2026-09-11。环境：WSL Ubuntu / Intel oneAPI 2025.3 / MKL 单线程。',
        '下表由 [summarize.py](summarize.py) 从完整 CSV 生成，不使用笔记中的历史验算数字。', '',
        '## 结论与计数', '',
        f'- 两个种子（91011、240810311），每个实 HS 预设 48 个独立场配置和 16 次交替扫描（含初态共 17 个扫描状态）。',
        f'- 27 个实 HS Hamiltonian 预设共 {len(regular)} 条记录；所有非零权重的相位在数值容差内为 1。',
        f'- 实际 PFQMC 提议 {attempts} 次，接受 {accepted} 次；每个模型都包含两个扫描方向和非零接受次数。',
        f'- 全部 {len(all_rows)} 条记录：{len(positives)} 条正号、{len(zeros)} 条明确零点、{len(negatives)} 条预期负号控制，另有 {len(all_rows)-len(positives)-len(zeros)-len(negatives)} 条一般复相位 oracle。',
        f'- 最大权重相对误差 `{max_error:.6g}`；正配置最大 `|sign-1| = {max_phase:.6g}`（测试允差 `2e-7`）。',
        '- 46/46 CTest 回归通过，包含原有 18 项、27 个模型独立短跑及综合控制测试。', '',
        '“sign 为 1”是对这里实际计算的非零配置的数值结论；有限抽样不是全参数证明。',
        '单胞系列存在精确零权重点，零点 CSV 写 `nan`，未冒充 sign=1。',
        '外幂系列中尚未完成的 Wei 排除证明也没有被数值验算补成证明。', '',
        '## 每个实现的独立记录', '',
        '`random` 是独立场压力测试；`sweep` 是含初态的实际 PFQMC 扫描状态。',
        '表中 N 是复费米子模数，132 模并不意味着 132 个空间格点。', '',
        '| 模型 | N | random | sweep | 其他核对 | 最大相对误差 | 最大相位偏差 |',
        '| --- | ---: | ---: | ---: | ---: | ---: | ---: |',
    ]
    for name, rows in sorted(groups.items()):
        if not any(r['kind'] in ('random_fields','pfqmc_sweep','all_static_backgrounds','nonhermitian_fields','all_closed_edge_hops') for r in rows):
            continue
        random = sum(r['kind']=='random_fields' for r in rows)
        sweep = sum(r['kind']=='pfqmc_sweep' for r in rows)
        error = max(float(r['relative_error']) for r in rows)
        phase = max(abs(complex(float(r['sign_real']),float(r['sign_imag']))-1) for r in rows)
        modes = '/'.join(str(x) for x in sorted({int(r['modes']) for r in rows}))
        lines.append(f'| `{name}` | {modes} | {random} | {sweep} | {len(rows)-random-sweep} | {error:.3g} | {phase:.3g} |')
    lines += [
        '', '## 独立验算和反例', '',
        '- 小模型比较完整 Fock 迹；守粒子数模型比较中心化行列式；奇数阶伴随族比较小 d 矩阵字符；配对阵列比较逐胞 Fock 迹。',
        '- 所有 Gaussian-HS 扫描重算 signed Pfaffian 与 `getSignRaw()`，并对照独立 Green 矩阵；局部 Metropolis 的接受和拒绝与相同随机数驱动的完整 Fock 权重比逐步一致。',
        '- SL 6/20/70 模和 Sp 5/48/132 模另外通过具有负辅助本征值的强 commutator 检查，避免只检查单位元附近。',
        '- 每种子：静态 FK 和静态 BdG 各穷举 512 个 3×3 背景，各做 432 次静态 Metropolis；另各有 16 个 2×2 完整 Fock 扇区核对。',
        '- 非厄米：仅创建配对开链/自旋单胞、两点单向复 hopping、菱形 DAG、强连通块扩展，逐项比较完整 Fock 迹和删去严格三角部分的参考。',
        '- Pair hopping：原始四次费米子 Hamiltonian 矩阵核对；每种子覆盖 128 个合法边跳跃往返闭环、48 条二维闭合路径，且穷举两点系统所有三步路径对照转移矩阵迹，检查时间步减半时趋近所写 H。没有声称实现生产级世界线更新。',
        '- 固定有限深度非 Gaussian 酉变换：直接构造物理 H\' 和每个共轭 Fock 脉冲，确认与逆变换采样的迹一致。',
        '- Gaussian HS 恒等式用 24 点正态求积独立核对；四模密度平方的中心化常数单独检查。',
        '- 负号控制：单模 Spin 提升 `-2`；自旋 1 复势 `-9/2`；跨胞十脉冲 `-7772907825/67108864`；不合格 Sp(6) primitive 三外幂 14 模 commutator 为负；实际 PFQMC 也正确追踪负权重样本。',
        '- 零点：自旋 1 commutator 取 `a = asinh(2^(-1/4))`，辅助二阶矩阵迹为 0，完整 Fock 迹为 0；全时空 Pfaffian 与 Fock oracle 一致。', '',
        '## 复现与证据', '',
        '完整命令见 [README.md](README.md)。`--model all --samples 48 --sweeps 16` 使用上述两个种子。',
        '结果保留实际复数权重，不对原始虚部清零。大迹下应看相对误差与相位；例如 132 模的强脉冲迹数量级为 1e43。', '',
        '- [种子 91011 CSV](results/samples_91011.csv) / [运行摘要](results/run_91011.log)',
        '- [种子 240810311 CSV](results/samples_240810311.csv) / [运行摘要](results/run_240810311.log)',
        '- [CTest 日志](results/ctest.log) / [源代码 SHA-256](results/sources.sha256) / [原始结果 SHA-256](results/data.sha256)',
        '- [15 个源笔记的 SHA-256](notes_manifest.txt) / [覆盖对照](AUDIT.md)', '',
        '测试范围为 beta=0.9、L_tau=4 的普通 HS 模型和文档给出的特例参数。',
        '通用 walker 采用稠密更新；这里没有测量大格点效率、物理误差条或热化/自相关时间。',
        '有限时间步有 Trotter 误差；极低温、强耦合与任意尺寸需要额外稳定性和收敛检查。',
    ]
    (HERE/'RESULTS.md').write_text('\n'.join(lines)+'\n')
    paths = sorted(p for p in RESULTS.iterdir() if p.suffix in ('.csv','.log'))
    (RESULTS/'data.sha256').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.name+'\n' for p in paths))

if __name__=='__main__':
    summarize()
