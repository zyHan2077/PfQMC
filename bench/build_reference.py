"""Build a baseline with the current harness/toolchain and unoptimized kernels.
Run inside the documented oneAPI shell; never changes the working source files.
"""
import argparse
from pathlib import Path
import shutil
import subprocess

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--ref',default='3832a1f')
p.add_argument('--cmake',default='/home/han/PfQMC/build/python-tools/cmake/data/bin/cmake')
p.add_argument('--eigen',default='/home/han/PfQMC/build/eigen')
a=p.parse_args()
root=Path('build/reference-source').resolve()
root.mkdir(parents=True,exist_ok=True)
for directory in ['inc','src','test','sign_free_verify']:
    shutil.copytree(directory,root/directory,dirs_exist_ok=True)
(root/'bench').mkdir(exist_ok=True)
shutil.copy2('bench/pfqmc_bench.cpp',root/'bench/pfqmc_bench.cpp')
for file in ['CMakeLists.txt','main.cpp']:
    shutil.copy2(file,root/file)
for file in ['inc/qr_udt.h','inc/spinless_tV.h']:
    (root/file).write_bytes(subprocess.check_output(['git','show',f'{a.ref}:{file}']))
subprocess.run([a.cmake,'-S',str(root),'-B','build/reference',
                '-DCMAKE_CXX_COMPILER=mpiicpx',f'-DEIGEN3_INCLUDE_DIR={a.eigen}',
                '-DCMAKE_PREFIX_PATH=/home/han/miniconda3'],check=True)
subprocess.run([a.cmake,'--build','build/reference','--target','pfqmc_bench','-j','2'],check=True)
Path('build/perf').mkdir(exist_ok=True)
shutil.copy2('build/reference/pfqmc_bench','build/perf/pfqmc_bench_baseline')
print('Baseline kernels:',subprocess.check_output(['git','rev-parse',a.ref],text=True).strip())
