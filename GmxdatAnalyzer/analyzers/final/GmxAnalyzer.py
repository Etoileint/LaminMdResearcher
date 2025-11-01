#!/usr/bin/env python3
"""
GROMACS综合批量分析脚本
支持RMSD、RMSF、二级结构、坐标提取、SASA等多种分析类型
具有断点续传、并行处理、灵活的对象选择等功能

使用说明:
1. 确保已安装GROMACS并配置好环境
2. 准备必要的输入文件
3. 运行脚本
输入文件要求:
- 轨迹文件 (.xtc, .trr等)
- 拓扑文件 (.tpr, .gro, .pdb等) 
- 索引文件 (.ndx) - 包含要分析的所有原子组
- 索引组列表文件 (.txt) - 包含组号和组名的对应关系
示例用法:
# RMSD分析 - 所有对象类型
python GmxAnalyzer.py \
    -f trajectory.xtc \
    -s topology.tpr \
    -n index.ndx \
    -nl datanalyse_ndxlist.txt \
    -o results \
    -v pro-h bb \
    -t rmsd \
    -g system chains fragments dimers dimer_fragments \
    -j 4
# RMSF分析 - 计算均方根波动
python GmxAnalyzer.py \
    -f trajectory.xtc \
    -s topology.tpr \
    -n index.ndx \
    -nl datanalyse_ndxlist.txt \
    -o results \
    -v pro-h ca \
    -t rmsf \
    -g system chains \
    --res \
    -j 4
# 二级结构分析 - 只分析单链和片段
python GmxAnalyzer.py \
    -f trajectory.xtc \
    -s topology.tpr \
    -n index.ndx \
    -nl datanalyse_ndxlist.txt \
    -o results \
    -v pro-h \
    -t ss \
    -g chains fragments \
    --hmode dssp
# SASA分析 - 计算溶剂可及表面积
python GmxAnalyzer.py \
    -f trajectory.xtc \
    -s topology.tpr \
    -n index.ndx \
    -nl datanalyse_ndxlist.txt \
    -o results \
    -v pro-h \
    -t sasa \
    -g system chains \
    -j 4
# 坐标提取 - 指定时间范围
python GmxAnalyzer.py \
    -f trajectory.xtc \
    -s topology.tpr \
    -n index.ndx \
    -nl datanalyse_ndxlist.txt \
    -o results \
    -v pro-h bb \
    -t coords \
    -g chains dimers \
    -b 100 -e 1000 -tu ps
# 高级RMSF分析 - 各向异性温度因子
python GmxAnalyzer.py \
    -f trajectory.xtc \
    -s topology.tpr \
    -n index.ndx \
    -nl datanalyse_ndxlist.txt \
    -o results \
    -v pro-h \
    -t rmsf \
    -g system chains \
    --aniso --res
参数说明:
必需参数:
  -f, --trajectory     轨迹文件路径
  -s, --topology       拓扑文件路径  
  -n, --ndx            索引文件路径
  -nl, --ndx-list      索引组列表文件路径
  -o, --output         输出目录
  -v, --versions       分析版本 [pro, pro-h, ca, bb, mc, sc]
  -t, --analysis-type  分析类型 [rmsd, rmsf, ss, sasa, coords]
可选参数:
  -g, --group-types    分析组类型 [system, chains, fragments, dimers, dimer_fragments]
  -b, --begin-time     开始时间 (默认: 0)
  -e, --end-time       结束时间 (默认: -1, 即到最后)
  -tu, --time-unit     时间单位 [ps, ns, us, ms, s] (默认: ps)
  -j, --jobs           并行任务数 (默认: 1)
分析类型特定参数:
RMSF分析:
  --no-fit    不进行最小二乘拟合 (默认进行拟合)
  --res       计算每个残基的平均值
  --aniso     计算各向异性温度因子
二级结构分析:
  --hmode     氢原子处理模式 [gromacs, dssp] (默认: dssp)
  --hbond     氢键定义方式 [energy, geometry] (默认: energy) 
  --no-clear  不清除有缺陷的残基 (默认清除)
SASA分析:
  --probe     溶剂探针半径 (nm) (默认: 0.14)
  --ndots     每个球体的点数 (默认: 24)
  --no-pbc    不使用周期性边界条件 (默认使用)
版本说明:
  pro     蛋白质原子 (包括氢原子)
  pro-h   蛋白质重原子 (不包括氢原子)
  ca      C-alpha原子
  bb      骨架原子 (N, CA, C)
  mc      主链原子 (N, CA, C, O)
  sc      侧链原子
组类型说明:
  system          系统默认组 (整个蛋白质)
  chains          单链组 (如AAA_pro, AAB_pro等)
  fragments       片段组 (如AAE-Head-20505-21022_pro等)  
  dimers          二聚体组 (如AAA_AAB_pro等)
  dimer_fragments 二聚体片段组 (如AAA_AAB-Linker_1_pro等)
输出目录结构:
results/
├── rmsd/
│   ├── pro-h/
│   │   ├── system/
│   │   │   └── whole_protein_pro-h.xvg
│   │   ├── chains/
│   │   │   ├── AAA_pro-h.xvg
│   │   │   └── ...
│   │   ├── fragments/
│   │   │   ├── AAE_Head_pro-h.xvg
│   │   │   └── ...
│   │   ├── dimers/
│   │   │   ├── AAA_AAB_pro-h.xvg
│   │   │   └── ...
│   │   └── dimer_fragments/
│   │       ├── AAA_AAB-Linker_1_pro-h.xvg
│   │       └── ...
│   └── bb/
│       └── ... (类似结构)
├── rmsf/
│   ├── pro-h/
│   │   ├── system/
│   │   │   └── whole_protein_pro-h_rmsf.xvg
│   │   ├── chains/
│   │   │   ├── AAA_pro-h_rmsf.xvg
│   │   │   └── ...
│   │   └── ...
│   └── ... (其他版本)
├── ss/
│   ├── pro-h/
│   │   ├── chains/
│   │   │   ├── AAA_pro-h_dssp.dat
│   │   │   ├── AAA_pro-h_dssp_num.xvg
│   │   │   └── ...
│   │   └── ...
│   └── ... (类似结构)
├── sasa/
│   ├── pro-h/
│   │   ├── system/
│   │   │   ├── whole_protein_pro-h_sasa.xvg
│   │   │   └── whole_protein_pro-h_sasa_resavg.xvg
│   │   ├── chains/
│   │   │   ├── AAA_pro-h_sasa.xvg
│   │   │   ├── AAA_pro-h_sasa_resavg.xvg
│   │   │   └── ...
│   │   └── ...
│   └── ... (类似结构)
└── coords/
    └── ... (类似结构)
状态管理:
- 脚本支持断点续传，会自动记录已完成的任务
- 状态文件保存在输出目录中: {analysis_type}_status.json
- 如果中途中断，重新运行会跳过已完成的任务
- 支持从现有文件恢复状态记录
二级结构符号说明:
H — alpha-helix
B — residue in isolated beta-bridge  
E — extended strand that participates in beta-ladder
G — 3_10-helix
I — pi-helix
P — kappa-helix (poly-proline II helix)
S — bend
T — hydrogen-bonded turn
= — break
~ — loop (no special secondary structure designation)
RMSF分析说明:
RMSF (Root Mean Square Fluctuation) 表示原子位置在轨迹中的均方根波动
- 输出文件包含每个原子的RMSF值 (单位: nm)
- 使用--res参数可计算每个残基的平均RMSF值
- 使用--aniso参数可计算各向异性温度因子
- 使用--no-fit参数可禁用结构拟合 (需确保参考结构和轨迹匹配)
SASA分析说明:
SASA (Solvent Accessible Surface Area) 计算溶剂可及表面积
- 输出总SASA随时间变化文件 (*_sasa.xvg)
- 输出残基平均SASA文件 (*_sasa_resavg.xvg)
- 使用--probe参数可调整溶剂探针半径 (默认: 0.14 nm)
- 使用--ndots参数可调整精度 (默认: 24)
- 使用--no-pbc参数可禁用周期性边界条件
"""

import os
import sys
import argparse
import subprocess
import re
import json
import time
import threading
from pathlib import Path
import logging
from typing import Dict, List, Tuple, Set, Optional, Callable, Any
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime

# 设置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

class MultiFileStatusManager:
    """多文件状态管理器，支持不同分析类型的多文件输出"""
    
    def __init__(self, output_dir: str, analysis_type: str):
        """
        初始化状态管理器
        
        Args:
            output_dir: 输出目录
            analysis_type: 分析类型
        """
        self.status_file = os.path.join(output_dir, f"{analysis_type}_status.json")
        self.lock = threading.Lock()
        self.analysis_type = analysis_type
        
        # 定义不同分析类型的文件模式
        self.file_patterns = {
            'rmsd': {
                'pattern': 'single',  # 单文件输出
                'extensions': ['.xvg']
            },
            'rmsf': {
                'pattern': 'single',  # 单文件输出
                'extensions': ['.xvg']
            },
            'ss': {
                'pattern': 'double',  # 双文件输出
                'extensions': ['.dat', '_num.xvg']
            },
            'sasa': {
                'pattern': 'double',  # 双文件输出：总SASA和残基平均SASA
                'extensions': ['_sasa.xvg', '_sasa_resavg.xvg']
            },
            'coords': {
                'pattern': 'single',  # 单文件输出
                'extensions': ['.xvg']
            }
        }
        
        # 默认状态数据结构
        self.status_data = {
            'analysis_type': analysis_type,
            'version': '2.0',  # 状态文件版本
            'start_time': datetime.now().isoformat(),
            'last_update': datetime.now().isoformat(),
            'completed_tasks': {},  # task_id -> task_info
            'failed_tasks': {},     # task_id -> error_info
            'parameters': {}
        }
        self.loaded = False
    
    def load_status(self):
        """加载状态文件"""
        if os.path.exists(self.status_file):
            try:
                with self.lock:
                    with open(self.status_file, 'r', encoding='utf-8') as f:
                        loaded_data = json.load(f)
                        
                        # 验证分析类型匹配
                        if loaded_data.get('analysis_type') == self.analysis_type:
                            self.status_data = loaded_data
                            self.loaded = True
                            logger.info(f"成功加载状态文件: {self.status_file}")
                        else:
                            logger.warning(f"状态文件分析类型不匹配: {loaded_data.get('analysis_type')} != {self.analysis_type}")
            except Exception as e:
                logger.warning(f"加载状态文件失败，将创建新状态文件: {e}")
        else:
            logger.info("未找到状态文件，将创建新状态文件")
    
    def save_status(self):
        """保存状态文件（覆盖模式）"""
        try:
            with self.lock:
                self.status_data['last_update'] = datetime.now().isoformat()
                # 使用临时文件避免写入过程中损坏
                temp_file = self.status_file + '.tmp'
                with open(temp_file, 'w', encoding='utf-8') as f:
                    json.dump(self.status_data, f, ensure_ascii=False, indent=2)
                
                # 原子性替换
                os.replace(temp_file, self.status_file)
        except Exception as e:
            logger.error(f"保存状态文件失败: {e}")
    
    def get_task_id(self, output_file: str) -> str:
        """根据输出文件路径生成任务ID"""
        # 对于二级结构分析，使用基础路径作为任务ID
        if self.analysis_type == 'ss':
            return output_file  # output_file已经是基础路径
        elif self.analysis_type == 'sasa':
            # 对于SASA分析，使用基础路径作为任务ID
            return output_file
        else:
            return output_file
    
    def get_output_files(self, task_id: str) -> List[str]:
        """获取任务对应的所有输出文件路径"""
        pattern_info = self.file_patterns.get(self.analysis_type)
        if not pattern_info:
            return [task_id]
        
        if pattern_info['pattern'] == 'single':
            # 单文件模式
            return [task_id]
        elif pattern_info['pattern'] == 'double':
            # 双文件模式
            extensions = pattern_info['extensions']
            return [task_id + ext for ext in extensions]
        else:
            # 默认单文件模式
            return [task_id]
    
    def add_completed_task(self, task_id: str, task_info: dict):
        """添加已完成的任务记录"""
        with self.lock:
            task_info['completed_time'] = datetime.now().isoformat()
            task_info['output_files'] = self.get_output_files(task_id)
            
            # 记录文件大小和修改时间
            file_sizes = {}
            file_mtimes = {}
            for file_path in task_info['output_files']:
                if os.path.exists(file_path):
                    file_sizes[file_path] = os.path.getsize(file_path)
                    file_mtimes[file_path] = os.path.getmtime(file_path)
            
            task_info['file_sizes'] = file_sizes
            task_info['file_mtimes'] = file_mtimes
            
            self.status_data['completed_tasks'][task_id] = task_info
            # 如果之前失败过，移除失败记录
            if task_id in self.status_data['failed_tasks']:
                del self.status_data['failed_tasks'][task_id]
        
        self.save_status()
    
    def add_failed_task(self, task_id: str, error_info: dict):
        """添加失败的任务记录"""
        with self.lock:
            error_info['failed_time'] = datetime.now().isoformat()
            self.status_data['failed_tasks'][task_id] = error_info
        self.save_status()
    
    def is_task_completed(self, task_id: str) -> bool:
        """检查任务是否已完成"""
        return task_id in self.status_data['completed_tasks']
    
    def is_task_failed(self, task_id: str) -> bool:
        """检查任务是否失败"""
        return task_id in self.status_data['failed_tasks']
    
    def get_task_info(self, task_id: str) -> dict:
        """获取任务信息"""
        return self.status_data['completed_tasks'].get(task_id, {})
    
    def set_parameters(self, parameters: dict):
        """设置运行参数"""
        with self.lock:
            self.status_data['parameters'] = parameters
        self.save_status()
    
    def verify_task_integrity(self, task_id: str) -> bool:
        """
        验证任务输出文件的完整性
        
        Args:
            task_id: 任务ID
            
        Returns:
            任务输出是否完整
        """
        task_info = self.get_task_info(task_id)
        if not task_info:
            return False
        
        output_files = self.get_output_files(task_id)
        
        # 检查所有文件是否存在
        for file_path in output_files:
            if not os.path.exists(file_path):
                logger.warning(f"输出文件不存在: {file_path}")
                return False
        
        # 验证文件大小（如果记录过）
        recorded_sizes = task_info.get('file_sizes', {})
        for file_path in output_files:
            if file_path in recorded_sizes:
                current_size = os.path.getsize(file_path)
                recorded_size = recorded_sizes[file_path]
                
                # 如果文件大小显著小于记录大小，认为文件损坏
                if current_size < recorded_size * 0.9:
                    logger.warning(f"文件可能损坏: {file_path} (当前大小: {current_size}, 记录大小: {recorded_size})")
                    return False
        
        return True
    
    def get_completed_tasks(self) -> List[str]:
        """获取所有已完成的任务ID"""
        return list(self.status_data['completed_tasks'].keys())
    
    def get_failed_tasks(self) -> List[str]:
        """获取所有失败的任务ID"""
        return list(self.status_data['failed_tasks'].keys())

class GromacsGroupAnalyzer:
    """GROMACS组分析器基类"""
    
    def __init__(self, ndx_file: str, ndx_list_file: str):
        """
        初始化分析器
        
        Args:
            ndx_file: GROMACS索引文件路径
            ndx_list_file: 索引组列表文件路径
        """
        self.ndx_file = ndx_file
        self.ndx_list_file = ndx_list_file
        self.group_mapping = {}  # 组名到编号的映射
        self.group_info = {}     # 编号到组名的映射
        self.system_groups = {}  # 系统默认组
        self.chain_groups = {}   # 单链组
        self.fragment_groups = {} # 片段组
        self.dimer_groups = {}   # 二聚体组
        self.dimer_fragment_groups = {} # 二聚体片段组
        
        self._load_ndx_list()
        self._categorize_groups()
    
    def _load_ndx_list(self):
        """加载索引组列表文件"""
        try:
            with open(self.ndx_list_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split()
                    if len(parts) >= 2:
                        group_num = int(parts[0])
                        group_name = parts[1]
                        self.group_mapping[group_name] = group_num
                        self.group_info[group_num] = group_name
            logger.info(f"成功加载 {len(self.group_mapping)} 个索引组")
        except Exception as e:
            logger.error(f"加载索引组列表失败: {e}")
            raise
    
    def _categorize_groups(self):
        """对索引组进行分类"""
        for group_name, group_num in self.group_mapping.items():
            # 系统默认组（数字编号，没有前缀）
            if group_num <= 9:
                self.system_groups[group_name] = group_num
            
            # 单链组（格式：AAA_pro, AAA_pro-h等）
            elif re.match(r'^[A-Z]{3}_(pro|pro-h|ca|bb|mc|sc)$', group_name):
                self.chain_groups[group_name] = group_num
            
            # 片段组（格式：AAE-Head-20505-21022_pro等）
            elif re.match(r'^[A-Z]{3}-[A-Za-z0-9_]+-\d+-\d+_(pro|pro-h|ca|bb|mc|sc)$', group_name):
                self.fragment_groups[group_name] = group_num
            
            # 二聚体组（格式：AAA_AAB_pro, AAA_AAB_pro-h等）
            elif re.match(r'^[A-Z]{3}_[A-Z]{3}_(pro|pro-h|ca|bb|mc|sc)$', group_name):
                self.dimer_groups[group_name] = group_num
            
            # 二聚体片段组（格式：AAA_AAB-Linker_1_pro, AAA_AAB-Coil_1B_pro等）
            elif re.match(r'^[A-Z]{3}_[A-Z]{3}-[A-Za-z0-9_]+_(pro|pro-h|ca|bb|mc|sc)$', group_name):
                self.dimer_fragment_groups[group_name] = group_num
        
        logger.info(f"分类完成: 系统组 {len(self.system_groups)} 个, 单链组 {len(self.chain_groups)} 个, "
                   f"片段组 {len(self.fragment_groups)} 个, 二聚体组 {len(self.dimer_groups)} 个, "
                   f"二聚体片段组 {len(self.dimer_fragment_groups)} 个")
    
    def get_system_group_for_version(self, version: str) -> int:
        """
        获取系统默认组中对应版本的组号
        
        Args:
            version: 版本后缀 (pro, pro-h, ca, bb, mc, sc)
            
        Returns:
            组号
        """
        version_mapping = {
            'pro': 'Protein',
            'pro-h': 'Protein-H', 
            'ca': 'C-alpha',
            'bb': 'Backbone',
            'mc': 'MainChain',
            'sc': 'SideChain'
        }
        
        system_group_name = version_mapping.get(version)
        if system_group_name and system_group_name in self.system_groups:
            return self.system_groups[system_group_name]
        else:
            raise ValueError(f"未找到系统默认组对应的版本: {version}")
    
    def get_chain_groups_for_version(self, version: str) -> Dict[str, int]:
        """
        获取指定版本的所有单链组
        
        Args:
            version: 版本后缀
            
        Returns:
            字典: {链名: 组号}
        """
        chain_groups = {}
        pattern = re.compile(rf'^([A-Z]{{3}})_{version}$')
        
        for group_name, group_num in self.chain_groups.items():
            match = pattern.match(group_name)
            if match:
                chain_name = match.group(1)
                chain_groups[chain_name] = group_num
        
        return chain_groups
    
    def get_fragment_groups_for_version(self, version: str) -> Dict[str, Dict[str, int]]:
        """
        获取指定版本的所有片段组
        
        Args:
            version: 版本后缀
            
        Returns:
            字典: {链名: {片段名: 组号}}
        """
        fragment_groups = {}
        pattern = re.compile(rf'^([A-Z]{{3}})-([A-Za-z0-9_]+)-\d+-\d+_{version}$')
        
        for group_name, group_num in self.fragment_groups.items():
            match = pattern.match(group_name)
            if match:
                chain_name = match.group(1)
                fragment_name = match.group(2)
                
                if chain_name not in fragment_groups:
                    fragment_groups[chain_name] = {}
                
                fragment_groups[chain_name][fragment_name] = group_num
        
        return fragment_groups
    
    def get_dimer_groups_for_version(self, version: str) -> Dict[str, int]:
        """
        获取指定版本的所有二聚体组
        
        Args:
            version: 版本后缀
            
        Returns:
            字典: {二聚体名: 组号}
        """
        dimer_groups = {}
        pattern = re.compile(rf'^([A-Z]{{3}}_[A-Z]{{3}})_{version}$')
        
        for group_name, group_num in self.dimer_groups.items():
            match = pattern.match(group_name)
            if match:
                dimer_name = match.group(1)
                dimer_groups[dimer_name] = group_num
        
        return dimer_groups
    
    def get_dimer_fragment_groups_for_version(self, version: str) -> Dict[str, Dict[str, int]]:
        """
        获取指定版本的所有二聚体片段组
        
        Args:
            version: 版本后缀
            
        Returns:
            字典: {二聚体名: {片段名: 组号}}
        """
        dimer_fragment_groups = {}
        pattern = re.compile(rf'^([A-Z]{{3}}_[A-Z]{{3}})-([A-Za-z0-9_]+)_{version}$')
        
        for group_name, group_num in self.dimer_fragment_groups.items():
            match = pattern.match(group_name)
            if match:
                dimer_name = match.group(1)
                fragment_name = match.group(2)
                
                if dimer_name not in dimer_fragment_groups:
                    dimer_fragment_groups[dimer_name] = {}
                
                dimer_fragment_groups[dimer_name][fragment_name] = group_num
        
        return dimer_fragment_groups

class RMSDAnalyzer(GromacsGroupAnalyzer):
    """RMSD分析器"""
    
    def analyze(self, trajectory: str, topology: str, output_file: str, 
                group_num: int, begin_time: float = 0, end_time: float = -1,
                time_unit: str = 'ps', **kwargs) -> bool:
        """
        运行GROMACS RMSD计算
        
        Args:
            trajectory: 轨迹文件路径
            topology: 拓扑文件路径
            output_file: 输出文件路径
            group_num: 索引组号
            begin_time: 开始时间
            end_time: 结束时间
            time_unit: 时间单位
            
        Returns:
            是否成功
        """
        try:
            # 构建命令
            cmd = [
                'gmx', 'rms', 
                '-f', trajectory,
                '-s', topology, 
                '-n', self.ndx_file,
                '-o', output_file,
                '-tu', time_unit
            ]
            
            # 如果指定了时间范围，添加时间参数
            if begin_time > 0:
                cmd.extend(['-b', str(begin_time)])
            if end_time > 0:
                cmd.extend(['-e', str(end_time)])
            
            # 准备输入（选择两次相同的组进行RMSD计算）
            input_str = f"{group_num}\n{group_num}\n"
            
            logger.info(f"计算RMSD: 组 {group_num} -> {output_file}")
            
            # 执行命令
            result = subprocess.run(
                cmd, 
                input=input_str, 
                text=True, 
                capture_output=True,
                check=True
            )
            
            return True
            
        except subprocess.CalledProcessError as e:
            logger.error(f"GROMACS RMSD计算失败: {e}")
            logger.error(f"错误输出: {e.stderr}")
            return False
        except Exception as e:
            logger.error(f"执行RMSD计算时发生错误: {e}")
            return False

class RMSFAnalyzer(GromacsGroupAnalyzer):
    """RMSF分析器"""
    
    def analyze(self, trajectory: str, topology: str, output_file: str, 
                group_num: int, begin_time: float = 0, end_time: float = -1,
                time_unit: str = 'ps', fit: bool = True, res: bool = False,
                aniso: bool = False, **kwargs) -> bool:
        """
        运行GROMACS RMSF计算
        
        Args:
            trajectory: 轨迹文件路径
            topology: 拓扑文件路径
            output_file: 输出文件路径
            group_num: 索引组号
            begin_time: 开始时间
            end_time: 结束时间
            time_unit: 时间单位（不支持）
            fit: 是否进行最小二乘拟合
            res: 是否计算每个残基的平均值
            aniso: 是否计算各向异性温度因子
            
        Returns:
            是否成功
        """
        try:
            # 构建命令
            cmd = [
                'gmx', 'rmsf', 
                '-f', trajectory,
                '-s', topology, 
                '-n', self.ndx_file,
                '-o', output_file
            ]
            
            # '-tu', time_unit
            
            # 如果指定了时间范围，添加时间参数
            if begin_time > 0:
                cmd.extend(['-b', str(begin_time)])
            if end_time > 0:
                cmd.extend(['-e', str(end_time)])
            
            # 添加可选参数
            if not fit:
                cmd.append('-nofit')
            if res:
                cmd.append('-res')
            if aniso:
                cmd.append('-aniso')
            
            # 准备输入（选择要分析的组）
            input_str = f"{group_num}\n"
            
            logger.info(f"计算RMSF: 组 {group_num} -> {output_file}")
            
            # 执行命令
            result = subprocess.run(
                cmd, 
                input=input_str, 
                text=True, 
                capture_output=True,
                check=True
            )
            
            return True
            
        except subprocess.CalledProcessError as e:
            logger.error(f"GROMACS RMSF计算失败: {e}")
            logger.error(f"错误输出: {e.stderr}")
            return False
        except Exception as e:
            logger.error(f"执行RMSF计算时发生错误: {e}")
            return False

class SecondaryStructureAnalyzer(GromacsGroupAnalyzer):
    """二级结构分析器"""
    
    def analyze(self, trajectory: str, topology: str, output_file: str, 
                group_num: int, begin_time: float = 0, end_time: float = -1,
                time_unit: str = 'ps', hmode: str = 'dssp', hbond: str = 'energy',
                clear_defective: bool = True, **kwargs) -> bool:
        """
        运行GROMACS DSSP二级结构分析
        
        Args:
            trajectory: 轨迹文件路径
            topology: 拓扑文件路径
            output_file: 输出文件基础名（不含扩展名）
            group_num: 索引组号
            begin_time: 开始时间
            end_time: 结束时间
            time_unit: 时间单位
            hmode: 氢原子处理模式 (gromacs, dssp)
            hbond: 氢键定义方式 (energy, geometry)
            clear_defective: 是否清除有缺陷的残基
            
        Returns:
            是否成功
        """
        try:
            # 获取组名
            group_name = self.group_info.get(group_num)
            if not group_name:
                logger.error(f"未找到组号 {group_num} 对应的组名")
                return False
            
            # 构建命令 - 使用 -sel 参数指定索引组名称
            cmd = [
                'gmx', 'dssp', 
                '-f', trajectory,
                '-s', topology, 
                '-n', self.ndx_file,
                '-o', f"{output_file}.dat",
                '-num', f"{output_file}_num.xvg",
                '-hmode', hmode,
                '-hbond', hbond,
                '-tu', time_unit,
                '-sel', f'"{group_name}"'  # 使用 -sel 参数指定组名而不是交互式选择
            ]
            
            # 如果指定了时间范围，添加时间参数
            if begin_time > 0:
                cmd.extend(['-b', str(begin_time)])
            if end_time > 0:
                cmd.extend(['-e', str(end_time)])
            
            # 添加清除有缺陷残基选项
            if clear_defective:
                cmd.append('-clear')
            
            # 注意：不再需要标准输入，因为已经通过 -sel 参数指定了组
            
            logger.info(f"分析二级结构: 组 {group_num} ({group_name}) -> {output_file}.*")
            logger.debug(f"执行命令: {' '.join(cmd)}")
            
            # 执行命令 - 不再需要输入重定向
            result = subprocess.run(
                cmd, 
                text=True, 
                capture_output=True,
                check=True
            )
            
            # 检查输出文件是否生成
            dat_file = f"{output_file}.dat"
            num_file = f"{output_file}_num.xvg"
            
            if os.path.exists(dat_file) and os.path.exists(num_file):
                logger.info(f"二级结构分析成功: {dat_file}, {num_file}")
                return True
            else:
                logger.error(f"二级结构分析输出文件未生成: {dat_file}, {num_file}")
                return False
            
        except subprocess.CalledProcessError as e:
            logger.error(f"GROMACS DSSP分析失败: {e}")
            logger.error(f"错误输出: {e.stderr}")
            logger.error(f"标准输出: {e.stdout}")
            return False
        except Exception as e:
            logger.error(f"执行DSSP分析时发生错误: {e}")
            return False

class SASAAnalyzer(GromacsGroupAnalyzer):
    """SASA分析器"""
    
    def analyze(self, trajectory: str, topology: str, output_file: str, 
                group_num: int, begin_time: float = 0, end_time: float = -1,
                time_unit: str = 'ps', probe: float = 0.14, ndots: int = 24,
                pbc: bool = True, **kwargs) -> bool:
        """
        运行GROMACS SASA计算
        
        Args:
            trajectory: 轨迹文件路径
            topology: 拓扑文件路径
            output_file: 输出文件基础名（不含扩展名）
            group_num: 索引组号
            begin_time: 开始时间
            end_time: 结束时间
            time_unit: 时间单位
            probe: 溶剂探针半径 (nm)
            ndots: 每个球体的点数
            pbc: 是否使用周期性边界条件
            
        Returns:
            是否成功
        """
        try:
            # 构建命令 - 输出总SASA和残基平均SASA
            cmd = [
                'gmx', 'sasa', 
                '-f', trajectory,
                '-s', topology, 
                '-n', self.ndx_file,
                '-o', f"{output_file}_sasa.xvg",      # 总SASA
                '-or', f"{output_file}_sasa_resavg.xvg", # 残基平均SASA
                '-probe', str(probe),
                '-ndots', str(ndots),
                '-tu', time_unit
            ]
            
            # 如果指定了时间范围，添加时间参数
            if begin_time > 0:
                cmd.extend(['-b', str(begin_time)])
            if end_time > 0:
                cmd.extend(['-e', str(end_time)])
            
            # 添加周期性边界条件选项
            if not pbc:
                cmd.append('-nopbc')
            
            # 准备输入（选择表面计算组）
            input_str = f"{group_num}\n"
            
            logger.info(f"计算SASA: 组 {group_num} -> {output_file}_sasa.xvg 和 {output_file}_sasa_resavg.xvg")
            
            # 执行命令
            result = subprocess.run(
                cmd, 
                input=input_str, 
                text=True, 
                capture_output=True,
                check=True
            )
            
            return True
            
        except subprocess.CalledProcessError as e:
            logger.error(f"GROMACS SASA计算失败: {e}")
            logger.error(f"错误输出: {e.stderr}")
            return False
        except Exception as e:
            logger.error(f"执行SASA计算时发生错误: {e}")
            return False

class CoordsAnalyzer(GromacsGroupAnalyzer):
    """坐标提取分析器"""
    
    def analyze(self, trajectory: str, topology: str, output_file: str, 
                group_num: int, begin_time: float = 0, end_time: float = -1,
                time_unit: str = 'ps', **kwargs) -> bool:
        """
        运行GROMACS坐标提取
        
        Args:
            trajectory: 轨迹文件路径
            topology: 拓扑文件路径
            output_file: 输出文件路径
            group_num: 索引组号
            begin_time: 开始时间
            end_time: 结束时间
            time_unit: 时间单位
            
        Returns:
            是否成功
        """
        try:
            # 构建命令
            cmd = [
                'gmx', 'traj', 
                '-f', trajectory,
                '-s', topology, 
                '-n', self.ndx_file,
                '-ox', output_file,
                '-tu', time_unit
            ]
            
            # 如果指定了时间范围，添加时间参数
            if begin_time > 0:
                cmd.extend(['-b', str(begin_time)])
            if end_time > 0:
                cmd.extend(['-e', str(end_time)])
            
            # 准备输入（选择要提取坐标的组）
            input_str = f"{group_num}\n"
            
            logger.info(f"提取坐标: 组 {group_num} -> {output_file}")
            
            # 执行命令
            result = subprocess.run(
                cmd, 
                input=input_str, 
                text=True, 
                capture_output=True,
                check=True
            )
            
            return True
            
        except subprocess.CalledProcessError as e:
            logger.error(f"GROMACS坐标提取失败: {e}")
            logger.error(f"错误输出: {e.stderr}")
            return False
        except Exception as e:
            logger.error(f"执行坐标提取时发生错误: {e}")
            return False

class GromacsBatchAnalyzer:
    """GROMACS批量分析器主类"""
    
    def __init__(self, analysis_type: str, ndx_file: str, ndx_list_file: str):
        """
        初始化批量分析器
        
        Args:
            analysis_type: 分析类型
            ndx_file: 索引文件路径
            ndx_list_file: 索引列表文件路径
        """
        self.analysis_type = analysis_type
        self.analyzer = self._create_analyzer(analysis_type, ndx_file, ndx_list_file)
    
    def _create_analyzer(self, analysis_type: str, ndx_file: str, ndx_list_file: str):
        """创建对应的分析器"""
        analyzers = {
            'rmsd': RMSDAnalyzer,
            'rmsf': RMSFAnalyzer,
            'ss': SecondaryStructureAnalyzer,
            'sasa': SASAAnalyzer,
            'coords': CoordsAnalyzer
        }
        
        if analysis_type not in analyzers:
            raise ValueError(f"不支持的分析类型: {analysis_type}")
        
        return analyzers[analysis_type](ndx_file, ndx_list_file)
    
    def create_output_directories(self, base_dir: str, versions: List[str], 
                                 group_types: List[str]) -> Dict[str, Dict[str, str]]:
        """
        创建输出目录结构
        
        Args:
            base_dir: 基础输出目录
            versions: 版本列表
            group_types: 组类型列表
            
        Returns:
            目录路径字典
        """
        base_path = Path(base_dir) / self.analysis_type
        base_path.mkdir(parents=True, exist_ok=True)
        
        dirs = {}
        for version in versions:
            version_dir = base_path / version
            version_dir.mkdir(exist_ok=True)
            
            version_dirs = {'base': str(version_dir)}
            
            # 根据选择的组类型创建子目录
            for group_type in group_types:
                group_dir = version_dir / group_type
                group_dir.mkdir(exist_ok=True)
                version_dirs[group_type] = str(group_dir)
            
            dirs[version] = version_dirs
        
        return dirs
    
    def get_output_filename(self, group_type: str, name: str, version: str, 
                           output_dirs: Dict[str, Dict[str, str]]) -> str:
        """
        获取输出文件名
        
        Args:
            group_type: 组类型
            name: 组名称
            version: 版本
            output_dirs: 输出目录字典
            
        Returns:
            输出文件路径
        """
        version_dirs = output_dirs[version]
        
        if self.analysis_type == 'rmsd':
            return os.path.join(version_dirs[group_type], f"{name}_{version}.xvg")
        elif self.analysis_type == 'rmsf':
            return os.path.join(version_dirs[group_type], f"{name}_{version}_rmsf.xvg")
        elif self.analysis_type == 'ss':
            # 对于二级结构分析，返回基础文件名（不含扩展名）
            return os.path.join(version_dirs[group_type], f"{name}_{version}_dssp")
        elif self.analysis_type == 'sasa':
            # 对于SASA分析，返回基础文件名（不含扩展名）
            return os.path.join(version_dirs[group_type], f"{name}_{version}")
        elif self.analysis_type == 'coords':
            return os.path.join(version_dirs[group_type], f"{name}_{version}_coords.xvg")
        else:
            raise ValueError(f"未知的分析类型: {self.analysis_type}")
    
    def generate_analysis_tasks(self, versions: List[str], group_types: List[str],
                               output_dirs: Dict[str, Dict[str, str]], 
                               begin_time: float, end_time: float, time_unit: str,
                               **kwargs) -> List[Dict]:
        """
        生成分析任务列表
        
        Args:
            versions: 版本列表
            group_types: 组类型列表
            output_dirs: 输出目录字典
            begin_time: 开始时间
            end_time: 结束时间
            time_unit: 时间单位
            **kwargs: 其他参数
            
        Returns:
            任务列表
        """
        tasks = []
        
        for version in versions:
            version_dirs = output_dirs[version]
            
            # 系统组分析
            if 'system' in group_types:
                try:
                    group_num = self.analyzer.get_system_group_for_version(version)
                    output_file = self.get_output_filename('system', 'whole_protein', version, output_dirs)
                    tasks.append({
                        'type': 'system',
                        'version': version,
                        'name': 'whole_protein',
                        'group_num': group_num,
                        'output': output_file,
                        'begin_time': begin_time,
                        'end_time': end_time,
                        'time_unit': time_unit,
                        'kwargs': kwargs
                    })
                except ValueError as e:
                    logger.warning(f"跳过系统组分析 {version}: {e}")
            
            # 单链组分析
            if 'chains' in group_types:
                chain_groups = self.analyzer.get_chain_groups_for_version(version)
                for chain_name, group_num in chain_groups.items():
                    output_file = self.get_output_filename('chains', chain_name, version, output_dirs)
                    tasks.append({
                        'type': 'chains',
                        'version': version,
                        'name': chain_name,
                        'group_num': group_num,
                        'output': output_file,
                        'begin_time': begin_time,
                        'end_time': end_time,
                        'time_unit': time_unit,
                        'kwargs': kwargs
                    })
            
            # 片段组分析
            if 'fragments' in group_types:
                fragment_groups = self.analyzer.get_fragment_groups_for_version(version)
                for chain_name, fragments in fragment_groups.items():
                    for fragment_name, group_num in fragments.items():
                        full_name = f"{chain_name}_{fragment_name}"
                        output_file = self.get_output_filename('fragments', full_name, version, output_dirs)
                        tasks.append({
                            'type': 'fragments',
                            'version': version,
                            'name': full_name,
                            'group_num': group_num,
                            'output': output_file,
                            'begin_time': begin_time,
                            'end_time': end_time,
                            'time_unit': time_unit,
                            'kwargs': kwargs
                        })
            
            # 二聚体组分析
            if 'dimers' in group_types:
                dimer_groups = self.analyzer.get_dimer_groups_for_version(version)
                for dimer_name, group_num in dimer_groups.items():
                    output_file = self.get_output_filename('dimers', dimer_name, version, output_dirs)
                    tasks.append({
                        'type': 'dimers',
                        'version': version,
                        'name': dimer_name,
                        'group_num': group_num,
                        'output': output_file,
                        'begin_time': begin_time,
                        'end_time': end_time,
                        'time_unit': time_unit,
                        'kwargs': kwargs
                    })
            
            # 二聚体片段组分析
            if 'dimer_fragments' in group_types:
                dimer_fragment_groups = self.analyzer.get_dimer_fragment_groups_for_version(version)
                for dimer_name, fragments in dimer_fragment_groups.items():
                    for fragment_name, group_num in fragments.items():
                        full_name = f"{dimer_name}-{fragment_name}"
                        output_file = self.get_output_filename('dimer_fragments', full_name, version, output_dirs)
                        tasks.append({
                            'type': 'dimer_fragments',
                            'version': version,
                            'name': full_name,
                            'group_num': group_num,
                            'output': output_file,
                            'begin_time': begin_time,
                            'end_time': end_time,
                            'time_unit': time_unit,
                            'kwargs': kwargs
                        })
        
        return tasks
    
    def execute_analysis(self, trajectory: str, topology: str, tasks: List[Dict],
                        status_manager: MultiFileStatusManager, max_workers: int = 1) -> Dict[str, int]:
        """
        执行分析任务
        
        Args:
            trajectory: 轨迹文件路径
            topology: 拓扑文件路径
            tasks: 任务列表
            status_manager: 状态管理器
            max_workers: 最大并行任务数
            
        Returns:
            统计信息字典
        """
        stats = {
            'total': len(tasks),
            'completed': 0,
            'skipped': 0,
            'failed': 0
        }
        
        def should_process_task(task: Dict) -> bool:
            """判断任务是否需要处理"""
            output_file = task['output']
            task_id = status_manager.get_task_id(output_file)
            
            # 检查任务是否已完成且文件完整
            if status_manager.is_task_completed(task_id):
                if status_manager.verify_task_integrity(task_id):
                    return False  # 任务已完成且文件完整，跳过
                else:
                    logger.warning(f"任务文件不完整，重新处理: {task_id}")
                    return True
            
            # 检查任务是否失败
            if status_manager.is_task_failed(task_id):
                logger.info(f"重试之前失败的任务: {task_id}")
                return True
            
            # 新任务或需要重新处理的任务
            return True
        
        def process_task(task: Dict) -> Tuple[bool, str, Dict]:
            """处理单个任务"""
            output_file = task['output']
            task_id = status_manager.get_task_id(output_file)
            
            try:
                # 执行分析
                success = self.analyzer.analyze(
                    trajectory=trajectory,
                    topology=topology,
                    output_file=output_file,
                    group_num=task['group_num'],
                    begin_time=task['begin_time'],
                    end_time=task['end_time'],
                    time_unit=task['time_unit'],
                    **task['kwargs']
                )
                
                if success:
                    # 记录成功状态
                    task_info = {
                        'group_num': task['group_num'],
                        'name': task['name'],
                        'version': task['version'],
                        'group_type': task['type'],
                        'analysis_type': self.analysis_type,
                        'output_file': output_file
                    }
                    
                    status_manager.add_completed_task(task_id, task_info)
                    return True, task_id, task_info
                else:
                    # 记录失败状态
                    error_info = {
                        'error': 'Analysis failed',
                        'group_num': task['group_num'],
                        'name': task['name'],
                        'version': task['version'],
                        'group_type': task['type']
                    }
                    status_manager.add_failed_task(task_id, error_info)
                    return False, task_id, error_info
                    
            except Exception as e:
                # 记录异常状态
                error_info = {
                    'error': str(e),
                    'group_num': task['group_num'],
                    'name': task['name'],
                    'version': task['version'],
                    'group_type': task['type']
                }
                status_manager.add_failed_task(task_id, error_info)
                return False, task_id, error_info
        
        # 过滤需要处理的任务
        tasks_to_process = [task for task in tasks if should_process_task(task)]
        stats['skipped'] = stats['total'] - len(tasks_to_process)
        
        if not tasks_to_process:
            logger.info("所有任务已完成，无需处理")
            return stats
        
        logger.info(f"需要处理 {len(tasks_to_process)} 个任务，使用 {max_workers} 个并行工作线程")
        
        # 使用线程池执行任务
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # 提交所有任务
            future_to_task = {
                executor.submit(process_task, task): task for task in tasks_to_process
            }
            
            # 处理完成的任务
            for future in as_completed(future_to_task):
                task = future_to_task[future]
                try:
                    success, task_id, info = future.result()
                    if success:
                        stats['completed'] += 1
                        logger.info(f"任务完成: {task_id}")
                    else:
                        stats['failed'] += 1
                        logger.error(f"任务失败: {task_id} - {info.get('error', 'Unknown error')}")
                except Exception as e:
                    stats['failed'] += 1
                    logger.error(f"任务异常: {task['output']} - {e}")
        
        return stats

def scan_existing_files(output_dirs: Dict[str, Dict[str, str]], analysis_type: str) -> Set[str]:
    """
    扫描现有的输出文件
    
    Args:
        output_dirs: 输出目录字典
        analysis_type: 分析类型
        
    Returns:
        现有文件路径集合
    """
    existing_files = set()
    
    for version, dirs in output_dirs.items():
        for dir_type, dir_path in dirs.items():
            if dir_type == 'base':
                continue
                
            dir_obj = Path(dir_path)
            if dir_obj.exists():
                # 根据分析类型确定文件扩展名
                if analysis_type == 'rmsd':
                    for file_path in dir_obj.glob("*.xvg"):
                        existing_files.add(str(file_path))
                elif analysis_type == 'rmsf':
                    for file_path in dir_obj.glob("*_rmsf.xvg"):
                        existing_files.add(str(file_path))
                elif analysis_type == 'ss':
                    # 二级结构分析有两个文件
                    for dat_file in dir_obj.glob("*.dat"):
                        existing_files.add(str(dat_file))
                    for num_file in dir_obj.glob("*_num.xvg"):
                        existing_files.add(str(num_file))
                elif analysis_type == 'sasa':
                    # SASA分析有两个文件
                    for sasa_file in dir_obj.glob("*_sasa.xvg"):
                        existing_files.add(str(sasa_file))
                    for resavg_file in dir_obj.glob("*_sasa_resavg.xvg"):
                        existing_files.add(str(resavg_file))
                elif analysis_type == 'coords':
                    for file_path in dir_obj.glob("*_coords.xvg"):
                        existing_files.add(str(file_path))
    
    return existing_files

def recover_status_from_files(status_manager: MultiFileStatusManager, existing_files: Set[str]):
    """
    从现有文件恢复状态
    
    Args:
        status_manager: 状态管理器
        existing_files: 现有文件集合
    """
    logger.info("检测到状态文件丢失，从现有输出文件重新生成状态记录...")
    logger.warning("注意：请检查现有文件是否完整，状态文件丢失可能导致文件完整性无法验证")
    
    # 按任务ID分组文件
    task_files = {}
    
    for file_path in existing_files:
        if status_manager.analysis_type == 'ss':
            # 对于二级结构分析，提取基础路径作为任务ID
            if file_path.endswith('.dat'):
                task_id = file_path[:-4]  # 移除 .dat 扩展名
            elif file_path.endswith('_num.xvg'):
                task_id = file_path[:-8]  # 移除 _num.xvg 扩展名
            else:
                continue
            
            if task_id not in task_files:
                task_files[task_id] = []
            task_files[task_id].append(file_path)
        elif status_manager.analysis_type == 'sasa':
            # 对于SASA分析，提取基础路径作为任务ID
            if file_path.endswith('_sasa.xvg'):
                task_id = file_path[:-9]  # 移除 _sasa.xvg 扩展名
            elif file_path.endswith('_sasa_resavg.xvg'):
                task_id = file_path[:-16]  # 移除 _sasa_resavg.xvg 扩展名
            else:
                continue
            
            if task_id not in task_files:
                task_files[task_id] = []
            task_files[task_id].append(file_path)
        else:
            # 对于其他分析类型，文件路径就是任务ID
            task_id = file_path
            task_files[task_id] = [file_path]
    
    # 为每个完整的任务组创建状态记录
    for task_id, files in task_files.items():
        # 检查任务是否完整
        expected_files = status_manager.get_output_files(task_id)
        if set(files) == set(expected_files):
            # 所有预期文件都存在，创建状态记录
            task_info = {
                'recovered': True,
                'recovery_time': datetime.now().isoformat(),
                'output_files': files
            }
            
            # 记录文件大小和修改时间
            file_sizes = {}
            file_mtimes = {}
            for file_path in files:
                if os.path.exists(file_path):
                    file_sizes[file_path] = os.path.getsize(file_path)
                    file_mtimes[file_path] = os.path.getmtime(file_path)
            
            task_info['file_sizes'] = file_sizes
            task_info['file_mtimes'] = file_mtimes
            
            status_manager.add_completed_task(task_id, task_info)
    
    logger.info(f"已从 {len(task_files)} 个任务恢复状态记录")

def reconcile_files_and_status(status_manager: MultiFileStatusManager, output_dirs: Dict[str, Dict[str, str]]):
    """
    协调状态文件和输出文件的一致性
    
    Args:
        status_manager: 状态管理器
        output_dirs: 输出目录字典
    """
    # 扫描现有文件
    existing_files = scan_existing_files(output_dirs, status_manager.analysis_type)
    
    # 获取状态文件中记录的任务
    recorded_tasks = set(status_manager.get_completed_tasks())
    
    # 按任务ID分组现有文件
    existing_tasks = set()
    for file_path in existing_files:
        if status_manager.analysis_type == 'ss':
            if file_path.endswith('.dat'):
                task_id = file_path[:-4]
            elif file_path.endswith('_num.xvg'):
                task_id = file_path[:-8]
            else:
                continue
        elif status_manager.analysis_type == 'sasa':
            if file_path.endswith('_sasa.xvg'):
                task_id = file_path[:-9]
            elif file_path.endswith('_sasa_resavg.xvg'):
                task_id = file_path[:-16]
            else:
                continue
        else:
            task_id = file_path
        existing_tasks.add(task_id)
    
    # 情况1: 状态文件中有记录但实际文件不存在
    missing_tasks = recorded_tasks - existing_tasks
    if missing_tasks:
        logger.warning(f"发现 {len(missing_tasks)} 个在状态文件中记录但实际不存在的任务，将从状态记录中移除")
        for task_id in missing_tasks:
            if task_id in status_manager.status_data['completed_tasks']:
                del status_manager.status_data['completed_tasks'][task_id]
        status_manager.save_status()
    
    # 情况2: 实际文件存在但状态文件中无记录
    unrecorded_tasks = existing_tasks - recorded_tasks
    if unrecorded_tasks:
        logger.info(f"发现 {len(unrecorded_tasks)} 个状态文件中未记录的任务")
        # 这些任务将在后续处理中被重新生成
    
    # 情况3: 验证已记录任务的完整性
    valid_tasks = recorded_tasks & existing_tasks
    invalid_tasks = []
    
    for task_id in valid_tasks:
        if not status_manager.verify_task_integrity(task_id):
            invalid_tasks.append(task_id)
    
    if invalid_tasks:
        logger.warning(f"发现 {len(invalid_tasks)} 个可能损坏的任务，将重新生成")
        # 这些任务将在后续处理中被重新生成

def main():
    parser = argparse.ArgumentParser(description='GROMACS综合批量分析工具')
    
    # 必需参数
    parser.add_argument('-f', '--trajectory', required=True, help='轨迹文件路径 (.xtc, .trr)')
    parser.add_argument('-s', '--topology', required=True, help='拓扑文件路径 (.tpr, .gro, .pdb)')
    parser.add_argument('-n', '--ndx', required=True, help='GROMACS索引文件路径')
    parser.add_argument('-nl', '--ndx-list', required=True, help='索引组列表文件路径')
    parser.add_argument('-o', '--output', required=True, help='输出目录')
    parser.add_argument('-v', '--versions', nargs='+', required=True, 
                       choices=['pro', 'pro-h', 'ca', 'bb', 'mc', 'sc'],
                       help='要分析的版本列表')
    parser.add_argument('-t', '--analysis-type', required=True,
                       choices=['rmsd', 'rmsf', 'ss', 'sasa', 'coords'],
                       help='分析类型: rmsd(RMSD分析), rmsf(RMSF分析), ss(二级结构分析), sasa(SASA分析), coords(坐标提取)')
    
    # 可选参数
    parser.add_argument('-g', '--group-types', nargs='+', 
                       choices=['system', 'chains', 'fragments', 'dimers', 'dimer_fragments'],
                       default=['system', 'chains', 'fragments', 'dimers', 'dimer_fragments'],
                       help='要分析的组类别')
    parser.add_argument('-b', '--begin-time', type=float, default=0, help='开始时间')
    parser.add_argument('-e', '--end-time', type=float, default=-1, help='结束时间')
    parser.add_argument('-tu', '--time-unit', default='ps', choices=['ps', 'ns', 'us', 'ms', 's'],
                       help='时间单位')
    parser.add_argument('-j', '--jobs', type=int, default=1, 
                       help='并行任务数 (默认: 1，即串行执行)')
    
    # 分析类型特定参数
    parser.add_argument('--hmode', default='dssp', choices=['gromacs', 'dssp'],
                       help='二级结构分析氢原子处理模式（默认：dssp，适用于缺少氢原子的结构）')
    parser.add_argument('--hbond', default='energy', choices=['energy', 'geometry'],
                       help='二级结构分析氢键定义方式（默认：energy）')
    parser.add_argument('--no-clear', action='store_true', 
                       help='二级结构分析不清除有缺陷的残基（默认清除）')
    parser.add_argument('--no-fit', action='store_true', 
                       help='RMSF分析不进行最小二乘拟合（默认进行拟合）')
    parser.add_argument('--res', action='store_true', 
                       help='RMSF分析计算每个残基的平均值')
    parser.add_argument('--aniso', action='store_true', 
                       help='RMSF分析计算各向异性温度因子')
    parser.add_argument('--probe', type=float, default=0.14,
                       help='SASA分析溶剂探针半径 (nm)（默认：0.14）')
    parser.add_argument('--ndots', type=int, default=24,
                       help='SASA分析每个球体的点数（默认：24）')
    parser.add_argument('--no-pbc', action='store_true',
                       help='SASA分析不使用周期性边界条件（默认使用）')
    
    args = parser.parse_args()
    
    # 检查输入文件是否存在
    for file_path in [args.trajectory, args.topology, args.ndx, args.ndx_list]:
        if not os.path.exists(file_path):
            logger.error(f"文件不存在: {file_path}")
            sys.exit(1)
    
    # 初始化批量分析器
    try:
        batch_analyzer = GromacsBatchAnalyzer(args.analysis_type, args.ndx, args.ndx_list)
    except ValueError as e:
        logger.error(f"初始化分析器失败: {e}")
        sys.exit(1)
    
    # 创建输出目录
    output_dirs = batch_analyzer.create_output_directories(args.output, args.versions, args.group_types)
    
    # 初始化状态管理器
    status_manager = MultiFileStatusManager(args.output, args.analysis_type)
    status_manager.load_status()
    
    # 扫描现有输出文件
    existing_files = scan_existing_files(output_dirs, args.analysis_type)
    
    # 断点续跑逻辑处理
    if not status_manager.loaded:
        if not existing_files:
            logger.info("检测到新项目，从头开始运行")
        else:
            logger.warning("检测到状态文件丢失，但存在输出文件")
            recover_status_from_files(status_manager, existing_files)
    else:
        logger.info("检测到状态文件，进行文件和状态协调")
        reconcile_files_and_status(status_manager, output_dirs)
    
    # 准备分析参数
    analysis_kwargs = {}
    if args.analysis_type == 'ss':
        analysis_kwargs.update({
            'hmode': args.hmode,
            'hbond': args.hbond,
            'clear_defective': not args.no_clear
        })
    elif args.analysis_type == 'rmsf':
        analysis_kwargs.update({
            'fit': not args.no_fit,
            'res': args.res,
            'aniso': args.aniso
        })
    elif args.analysis_type == 'sasa':
        analysis_kwargs.update({
            'probe': args.probe,
            'ndots': args.ndots,
            'pbc': not args.no_pbc
        })
    
    # 保存运行参数
    parameters = {
        'trajectory': args.trajectory,
        'topology': args.topology,
        'ndx': args.ndx,
        'ndx_list': args.ndx_list,
        'versions': args.versions,
        'analysis_type': args.analysis_type,
        'group_types': args.group_types,
        'begin_time': args.begin_time,
        'end_time': args.end_time,
        'time_unit': args.time_unit,
        'jobs': args.jobs,
        'analysis_kwargs': analysis_kwargs
    }
    status_manager.set_parameters(parameters)
    
    # 生成分析任务
    tasks = batch_analyzer.generate_analysis_tasks(
        versions=args.versions,
        group_types=args.group_types,
        output_dirs=output_dirs,
        begin_time=args.begin_time,
        end_time=args.end_time,
        time_unit=args.time_unit,
        **analysis_kwargs
    )
    
    logger.info(f"生成了 {len(tasks)} 个分析任务")
    
    # 执行分析
    stats = batch_analyzer.execute_analysis(
        trajectory=args.trajectory,
        topology=args.topology,
        tasks=tasks,
        status_manager=status_manager,
        max_workers=args.jobs
    )
    
    # 输出统计信息
    logger.info("\n=== 处理完成汇总 ===")
    logger.info(f"总任务数: {stats['total']}")
    logger.info(f"跳过任务: {stats['skipped']}")
    logger.info(f"成功任务: {stats['completed']}")
    logger.info(f"失败任务: {stats['failed']}")
    logger.info(f"结果保存在: {os.path.join(args.output, args.analysis_type)}")
    logger.info(f"状态日志: {status_manager.status_file}")
    
    # 打印汇总信息
    print("\n=== 分析汇总 ===")
    print(f"分析类型: {args.analysis_type}")
    print(f"输出目录: {os.path.join(args.output, args.analysis_type)}")
    
    for version in args.versions:
        version_dirs = output_dirs[version]
        print(f"\n版本: {version}")
        for group_type in args.group_types:
            if group_type in version_dirs:
                print(f"  {group_type}: {version_dirs[group_type]}/")
    
    # 如果是二级结构分析，打印符号说明
    if args.analysis_type == 'ss':
        print("\n=== 二级结构符号说明 ===")
        print("H — alpha-helix")
        print("B — residue in isolated beta-bridge")
        print("E — extended strand that participates in beta-ladder")
        print("G — 3_10-helix")
        print("I — pi-helix")
        print("P — kappa-helix (poly-proline II helix)")
        print("S — bend")
        print("T — hydrogen-bonded turn")
        print("= — break")
        print("~ — loop (no special secondary structure designation)")
    
    # 如果是RMSF分析，打印说明
    elif args.analysis_type == 'rmsf':
        print("\n=== RMSF分析说明 ===")
        print("RMSF (Root Mean Square Fluctuation) 表示原子位置在轨迹中的均方根波动")
        print("输出文件包含每个原子的RMSF值 (nm)")
        if args.res:
            print("已启用残基平均模式：输出每个残基的平均RMSF值")
        if args.aniso:
            print("已启用各向异性温度因子计算")
        if args.no_fit:
            print("注意：未进行结构拟合，请确保参考结构和轨迹匹配")
    
    # 如果是SASA分析，打印说明
    elif args.analysis_type == 'sasa':
        print("\n=== SASA分析说明 ===")
        print("SASA (Solvent Accessible Surface Area) 计算溶剂可及表面积")
        print("输出文件包含:")
        print("  - 总SASA随时间变化 (*_sasa.xvg)")
        print("  - 残基平均SASA (*_sasa_resavg.xvg)")
        print(f"分析参数: 探针半径={args.probe} nm, 点数={args.ndots}")
        if args.no_pbc:
            print("注意：未使用周期性边界条件")
        else:
            print("已使用周期性边界条件")

if __name__ == "__main__":
    main()

