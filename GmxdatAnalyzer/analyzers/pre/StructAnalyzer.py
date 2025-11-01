#!/usr/bin/env python3
"""
从PDB/GRO文件中根据用户定义序列找出指定片段
支持PDB和GRO格式，可生成片段标识和PyMol染色命令
新增二聚体支持功能

基本用法：
python StructAnalyzer.py -f input.pdb -s sequence_definitions.txt -o results.txt -pm my_protein -ndx
注意：使用索引组功能时，优先手动make_ndx看一下具体的默认组是什么样，再在脚本中调整对应关系和起始组号
"""

import argparse
import os
import sys
import subprocess
from typing import List, Dict, Tuple, Optional

# 氨基酸三字母到单字母转换字典
THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'HID': 'H', 'HIE': 'H', 'HIP': 'H',  # 组氨酸变体
    'ASH': 'D', 'GLH': 'E',  # 质子化形式
}

def get_three_letter_chain_id(index: int) -> str:
    """
    将索引转换为三字母链标识符 (AAA, AAB, ..., ZZZ)
    支持最多 26*26*26 = 17576 条链
    """
    if index < 0 or index >= 26*26*26:
        raise ValueError(f"链索引 {index} 超出范围 (0-{26*26*26-1})")
    
    first = index // (26 * 26)
    second = (index % (26 * 26)) // 26
    third = index % 26
    
    return chr(65 + first) + chr(65 + second) + chr(65 + third)

def validate_file_extension(filepath: str, expected_extensions: List[str]) -> str:
    """验证文件扩展名并返回完整路径"""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"文件不存在: {filepath}")
    
    _, ext = os.path.splitext(filepath)
    if ext.lower() not in expected_extensions:
        raise ValueError(f"文件扩展名错误: {filepath}，期望: {expected_extensions}")
    
    return filepath

def get_output_path(input_file: str, user_specified: Optional[str], 
                   default_name: str, expected_ext: str = None) -> str:
    """
    统一处理输出路径
    
    Args:
        input_file: 输入文件路径
        user_specified: 用户指定的输出路径
        default_name: 默认文件名
        expected_ext: 期望的文件扩展名（可选）
    """
    if user_specified is None:
        # 使用默认名，输出到输入文件所在文件夹
        input_dir = os.path.dirname(input_file)
        if input_dir == '':
            input_dir = '.'
        output_path = os.path.join(input_dir, default_name)
    else:
        if os.path.isdir(user_specified):
            # 用户指定的是文件夹，使用默认名
            output_path = os.path.join(user_specified, default_name)
        else:
            # 用户指定的是完整路径
            output_path = user_specified
    
    # 检查文件扩展名
    if expected_ext:
        _, ext = os.path.splitext(output_path)
        if ext.lower() != expected_ext.lower():
            raise ValueError(f"输出文件扩展名错误: {output_path}，期望: {expected_ext}")
    
    # 确保输出目录存在
    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    return output_path

class StructureParser:
    """结构文件解析器基类"""
    
    def __init__(self, filename: str):
        self.filename = filename
        self.chains = {}
        self.dimers = []
        
    def parse(self):
        raise NotImplementedError("子类必须实现parse方法")
    
    def get_residue_sequence(self, chain_id: str) -> str:
        """获取指定链的氨基酸序列"""
        raise NotImplementedError("子类必须实现get_residue_sequence方法")
    
    def get_chain_residues(self, chain_id: str):
        """获取指定链的残基信息"""
        return self.chains.get(chain_id, [])
    
    def generate_dimers(self):
        """生成二聚体列表"""
        chain_ids = list(self.chains.keys())
        num_chains = len(chain_ids)
        
        if num_chains == 1:
            print("只有一条链，不生成二聚体")
            self.dimers = []
        elif num_chains % 2 == 0:
            # 偶数条链，按顺序两两配对
            self.dimers = []
            for i in range(0, num_chains, 2):
                dimer_name = f"{chain_ids[i]}_{chain_ids[i+1]}"
                self.dimers.append({
                    'name': dimer_name,
                    'chain1': chain_ids[i],
                    'chain2': chain_ids[i+1]
                })
            print(f"生成 {len(self.dimers)} 个二聚体: {[d['name'] for d in self.dimers]}")
        else:
            print(f"警告: 链数为奇数 ({num_chains})，无法完全配对生成二聚体")
            # 仍然生成前n-1条链的配对
            self.dimers = []
            for i in range(0, num_chains-1, 2):
                dimer_name = f"{chain_ids[i]}_{chain_ids[i+1]}"
                self.dimers.append({
                    'name': dimer_name,
                    'chain1': chain_ids[i],
                    'chain2': chain_ids[i+1]
                })
            print(f"生成 {len(self.dimers)} 个二聚体 (最后一条链未配对): {[d['name'] for d in self.dimers]}")

class PDBParser(StructureParser):
    """PDB文件解析器"""
    
    def parse(self):
        with open(self.filename, 'r') as f:
            current_chain = None
            current_residue = None
            residues = []
            
            for line in f:
                if line.startswith('ATOM'):
                    chain_id = line[21].strip()
                    if chain_id == '':
                        chain_id = 'A'  # 默认链ID
                    
                    res_name = line[17:20].strip()
                    res_seq = int(line[22:26].strip())
                    atom_name = line[12:16].strip()
                    
                    # 新链开始
                    if current_chain != chain_id:
                        if current_chain is not None and residues:
                            self.chains[current_chain] = residues
                        current_chain = chain_id
                        residues = []
                        current_residue = None
                    
                    # 新残基开始
                    if current_residue != res_seq:
                        current_residue = res_seq
                        residues.append({
                            'res_seq': res_seq,
                            'res_name': res_name,
                            'atom_name': atom_name,
                            'line': line.strip()
                        })
                
                elif line.startswith('TER'):
                    # 链结束
                    if current_chain is not None and residues:
                        self.chains[current_chain] = residues
                        residues = []
                        current_chain = None
                        current_residue = None
                
                elif line.startswith('END'):
                    # 文件结束
                    if current_chain is not None and residues:
                        self.chains[current_chain] = residues
                    break
        
        # 生成二聚体
        self.generate_dimers()
    
    def get_residue_sequence(self, chain_id: str) -> str:
        """获取PDB链的氨基酸序列"""
        if chain_id not in self.chains:
            return ""
        
        sequence = ""
        last_res_seq = None
        
        for residue in self.chains[chain_id]:
            # 避免重复残基
            if residue['res_seq'] != last_res_seq:
                one_letter = THREE_TO_ONE.get(residue['res_name'], 'X')
                sequence += one_letter
                last_res_seq = residue['res_seq']
        
        return sequence

class GROParser(StructureParser):
    """GRO文件解析器"""
    
    def parse(self):
        with open(self.filename, 'r') as f:
            # 跳过标题行和原子数行
            title = f.readline().strip()
            num_atoms = int(f.readline().strip())
            
            current_chain_index = 0
            current_chain = get_three_letter_chain_id(current_chain_index)
            residues = []
            last_res_seq = None
            
            # 生成绝对原子序号
            absolute_atom_seq = 1
            
            for i in range(num_atoms):
                line = f.readline()
                if not line:
                    break
                
                # GRO格式解析
                res_seq = int(line[0:5].strip())
                res_name = line[5:10].strip()
                atom_name = line[10:15].strip()
                # 忽略文件中的原子序号，使用绝对原子序号
                file_atom_seq = int(line[15:20].strip())
                
                # 检测新链：残基序号不连续
                if last_res_seq is not None and res_seq < last_res_seq:
                    # 新链开始
                    if residues:
                        self.chains[current_chain] = residues
                        residues = []
                    # 更新链ID
                    current_chain_index += 1
                    current_chain = get_three_letter_chain_id(current_chain_index)
                
                residues.append({
                    'res_seq': res_seq,
                    'res_name': res_name,
                    'atom_seq': absolute_atom_seq,  # 使用绝对原子序号
                    'atom_name': atom_name,
                    'line': line.strip()
                })
                
                last_res_seq = res_seq
                absolute_atom_seq += 1  # 递增绝对原子序号
            
            # 添加最后一链
            if residues:
                self.chains[current_chain] = residues
        
        # 生成二聚体
        self.generate_dimers()
    
    def get_residue_sequence(self, chain_id: str) -> str:
        """获取GRO链的氨基酸序列"""
        if chain_id not in self.chains:
            return ""
        
        sequence = ""
        last_res_seq = None
        
        for residue in self.chains[chain_id]:
            # 避免重复残基
            if residue['res_seq'] != last_res_seq:
                one_letter = THREE_TO_ONE.get(residue['res_name'], 'X')
                sequence += one_letter
                last_res_seq = residue['res_seq']
        
        return sequence

class SequenceSearcher:
    """序列搜索器"""
    
    def __init__(self, structure_parser: StructureParser):
        self.parser = structure_parser
        self.found_fragments = {}
        self.found_dimer_fragments = {}
    
    def search_sequence(self, front_seq: str, target_seq: str, back_seq: str, 
                       fragment_name: str, color: str = "") -> Dict:
        """
        在所有链中搜索指定序列模式
        
        Args:
            front_seq: 前段序列，'-'表示任意，'N*'表示N端
            target_seq: 目标序列
            back_seq: 后段序列，'-'表示任意，'C*'表示C端
            fragment_name: 片段名称
            color: 颜色（可选）
        """
        results = []
        
        for chain_id in self.parser.chains:
            chain_sequence = self.parser.get_residue_sequence(chain_id)
            
            if not chain_sequence:
                continue
            
            # 搜索序列模式
            positions = self._find_sequence_positions(
                chain_sequence, front_seq, target_seq, back_seq
            )
            
            for start_pos, end_pos in positions:
                # 获取实际的残基/原子序号
                if isinstance(self.parser, PDBParser):
                    start_num, end_num = self._get_pdb_residue_numbers(
                        chain_id, start_pos, end_pos
                    )
                else:  # GROParser
                    start_num, end_num = self._get_gro_atom_numbers(
                        chain_id, start_pos, end_pos
                    )
                
                fragment_id = f"{chain_id}-{fragment_name}-{start_num}-{end_num}"
                results.append({
                    'fragment_id': fragment_id,
                    'chain': chain_id,
                    'fragment_name': fragment_name,
                    'start': start_num,
                    'end': end_num,
                    'color': color,
                    'sequence': target_seq
                })
        
        self.found_fragments[fragment_name] = results
        return results
    
    def search_dimer_fragments(self):
        """在二聚体中搜索片段（基于单链片段结果组合）"""
        dimer_results = {}
        
        for dimer in self.parser.dimers:
            dimer_name = dimer['name']
            chain1 = dimer['chain1']
            chain2 = dimer['chain2']
            
            dimer_fragments = []
            
            # 对于每个片段类型，检查两条链是否都有该片段
            for fragment_name, fragments in self.found_fragments.items():
                # 查找两条链上的该片段
                chain1_fragments = [f for f in fragments if f['chain'] == chain1]
                chain2_fragments = [f for f in fragments if f['chain'] == chain2]
                
                # 如果两条链都有该片段，则组合成二聚体片段
                if chain1_fragments and chain2_fragments:
                    # 假设一种片段在一条链上只出现一次
                    chain1_frag = chain1_fragments[0]
                    chain2_frag = chain2_fragments[0]
                    
                    dimer_fragment_id = f"{dimer_name}-{fragment_name}"
                    
                    dimer_fragments.append({
                        'fragment_id': dimer_fragment_id,
                        'dimer': dimer_name,
                        'fragment_name': fragment_name,
                        'chain1': chain1,
                        'chain2': chain2,
                        'chain1_start': chain1_frag['start'],
                        'chain1_end': chain1_frag['end'],
                        'chain2_start': chain2_frag['start'],
                        'chain2_end': chain2_frag['end'],
                        'color': chain1_frag['color'],  # 使用第一条链的颜色
                        'sequence': chain1_frag['sequence']
                    })
            
            if dimer_fragments:
                dimer_results[dimer_name] = dimer_fragments
        
        self.found_dimer_fragments = dimer_results
        return dimer_results
    
    def _find_sequence_positions(self, sequence: str, front_seq: str, 
                                target_seq: str, back_seq: str) -> List[Tuple[int, int]]:
        """在序列中查找匹配的位置"""
        positions = []
        seq_len = len(sequence)
        target_len = len(target_seq)
        
        # 处理特殊标记
        if front_seq == 'N*':
            # N端序列
            if sequence.startswith(target_seq):
                if back_seq == 'C*':
                    # 整个链都是目标序列
                    if sequence == target_seq:
                        positions.append((0, target_len))
                elif back_seq == '-' or sequence[target_len:].startswith(back_seq):
                    positions.append((0, target_len))
        
        elif back_seq == 'C*':
            # C端序列
            if sequence.endswith(target_seq):
                start_pos = seq_len - target_len
                if front_seq == '-' or sequence[:start_pos].endswith(front_seq):
                    positions.append((start_pos, seq_len))
        
        else:
            # 普通序列搜索
            search_pattern = ""
            if front_seq != '-':
                search_pattern += front_seq
            search_pattern += target_seq
            if back_seq != '-':
                search_pattern += back_seq
            
            pattern_len = len(search_pattern)
            front_len = len(front_seq) if front_seq != '-' else 0
            
            i = 0
            while i <= len(sequence) - pattern_len:
                if sequence[i:i+pattern_len] == search_pattern:
                    start_pos = i + front_len
                    end_pos = start_pos + len(target_seq)
                    positions.append((start_pos, end_pos))
                    i = end_pos  # 跳过已匹配的部分
                else:
                    i += 1
        
        return positions
    
    def _get_pdb_residue_numbers(self, chain_id: str, start_pos: int, end_pos: int) -> Tuple[int, int]:
        """获取PDB文件中实际的残基序号"""
        residues = self.parser.get_chain_residues(chain_id)
        
        # 获取唯一的残基
        unique_residues = []
        last_res_seq = None
        for residue in residues:
            if residue['res_seq'] != last_res_seq:
                unique_residues.append(residue)
                last_res_seq = residue['res_seq']
        
        start_residue = unique_residues[start_pos]
        end_residue = unique_residues[end_pos - 1]
        
        return start_residue['res_seq'], end_residue['res_seq']
    
    def _get_gro_atom_numbers(self, chain_id: str, start_pos: int, end_pos: int) -> Tuple[int, int]:
        """获取GRO文件中实际的原子序号"""
        residues = self.parser.get_chain_residues(chain_id)
        
        # 获取唯一的残基
        unique_residues = []
        last_res_seq = None
        for residue in residues:
            if residue['res_seq'] != last_res_seq:
                unique_residues.append(residue)
                last_res_seq = residue['res_seq']
        
        # 找到起始和结束残基的第一个原子
        start_residue = unique_residues[start_pos]
        end_residue = unique_residues[end_pos - 1]
        
        # 在原始数据中查找这些残基的原子范围
        start_atom = None
        end_atom = None
        
        for residue in residues:
            if residue['res_seq'] == start_residue['res_seq'] and start_atom is None:
                start_atom = residue['atom_seq']
            if residue['res_seq'] == end_residue['res_seq']:
                end_atom = residue['atom_seq']
        
        return start_atom, end_atom

class OutputGenerator:
    """输出生成器"""
    
    def __init__(self, structure_parser: StructureParser, is_gro: bool = False):
        self.parser = structure_parser
        self.is_gro = is_gro
    
    def generate_statistics(self, all_results: Dict, dimer_results: Dict) -> str:
        """生成统计结果"""
        output = "=== 片段搜索统计结果 ===\n\n"
        
        # 单链片段统计
        total_found = 0
        for fragment_name, results in all_results.items():
            count = len(results)
            total_found += count
            chains = set(result['chain'] for result in results)
            
            output += f"单链片段 '{fragment_name}':\n"
            output += f"  找到数量: {count}\n"
            output += f"  所在链: {', '.join(sorted(chains))}\n\n"
        
        # 二聚体片段统计
        dimer_fragment_count = 0
        for dimer_name, fragments in dimer_results.items():
            count = len(fragments)
            dimer_fragment_count += count
            fragment_names = set(frag['fragment_name'] for frag in fragments)
            
            output += f"二聚体 '{dimer_name}':\n"
            output += f"  片段数量: {count}\n"
            output += f"  片段类型: {', '.join(sorted(fragment_names))}\n\n"
        
        output += f"总计找到单链片段: {total_found}\n"
        output += f"总计找到二聚体片段: {dimer_fragment_count}\n"
        output += "=" * 30 + "\n\n"
        
        return output
    
    def generate_fragment_list(self, all_results: Dict, dimer_results: Dict) -> str:
        """生成片段列表"""
        output = "=== 找到的片段列表 ===\n\n"
        
        # 单链片段
        output += "单链片段:\n"
        for fragment_name, results in all_results.items():
            if results:
                output += f"片段 '{fragment_name}':\n"
                for result in results:
                    output += f"  {result['fragment_id']}\n"
                output += "\n"
        
        # 二聚体片段
        if dimer_results:
            output += "二聚体片段:\n"
            for dimer_name, fragments in dimer_results.items():
                if fragments:
                    output += f"二聚体 '{dimer_name}':\n"
                    for fragment in fragments:
                        output += f"  {fragment['fragment_id']}\n"
                    output += "\n"
        
        return output
    
    def generate_pymol_commands(self, all_results: Dict, dimer_results: Dict, object_name: str) -> str:
        """生成PyMOL染色命令"""
        output = f"# PyMOL染色命令 - 对象: {object_name}\n"
        output += f"# 生成自: {self.parser.filename}\n\n"
        
        # 对于GRO文件，先定义链
        if self.is_gro:
            output += "# 定义链（GRO文件无链标识符）\n"
            for chain_id in self.parser.chains:
                residues = self.parser.get_chain_residues(chain_id)
                if residues:
                    start_atom = residues[0]['atom_seq']
                    end_atom = residues[-1]['atom_seq']
                    output += f"select chain_{chain_id}, id {start_atom}-{end_atom}\n"
                    output += f"color gray, chain_{chain_id}\n"
            output += "\n"
        
        # 为单链片段生成选择命令
        output += "# 单链片段染色\n"
        for fragment_name, results in all_results.items():
            for result in results:
                chain_id = result['chain']
                start_num = result['start']
                end_num = result['end']
                color = result['color']
                
                if color:  # 如果有指定颜色
                    if self.is_gro:
                        # GRO文件使用原子序号
                        output += f"select {result['fragment_id']}, id {start_num}-{end_num}\n"
                        output += f"color {color}, {result['fragment_id']}\n"
                    else:
                        # PDB文件使用残基序号
                        output += f"select {result['fragment_id']}, chain {chain_id} and resi {start_num}-{end_num}\n"
                        output += f"color {color}, {result['fragment_id']}\n"
        
        # 为二聚体片段生成选择命令
        if dimer_results:
            output += "\n# 二聚体片段染色\n"
            for dimer_name, fragments in dimer_results.items():
                for fragment in fragments:
                    color = fragment['color']
                    
                    if color:  # 如果有指定颜色
                        if self.is_gro:
                            # GRO文件使用原子序号
                            chain1_atoms = f"id {fragment['chain1_start']}-{fragment['chain1_end']}"
                            chain2_atoms = f"id {fragment['chain2_start']}-{fragment['chain2_end']}"
                            output += f"select {fragment['fragment_id']}, {chain1_atoms} | {chain2_atoms}\n"
                            output += f"color {color}, {fragment['fragment_id']}\n"
                        else:
                            # PDB文件使用残基序号
                            chain1_residues = f"chain {fragment['chain1']} and resi {fragment['chain1_start']}-{fragment['chain1_end']}"
                            chain2_residues = f"chain {fragment['chain2']} and resi {fragment['chain2_start']}-{fragment['chain2_end']}"
                            output += f"select {fragment['fragment_id']}, {chain1_residues} | {chain2_residues}\n"
                            output += f"color {color}, {fragment['fragment_id']}\n"
        
        # 为二聚体生成分组命令
        if self.parser.dimers:
            output += "\n# 二聚体分组\n"
            for dimer in self.parser.dimers:
                dimer_name = dimer['name']
                chain1 = dimer['chain1']
                chain2 = dimer['chain2']
                
                if self.is_gro:
                    chain1_atoms = f"chain_{chain1}"
                    chain2_atoms = f"chain_{chain2}"
                    output += f"select dimer_{dimer_name}, {chain1_atoms} | {chain2_atoms}\n"
                else:
                    chain1_residues = f"chain {chain1}"
                    chain2_residues = f"chain {chain2}"
                    output += f"select dimer_{dimer_name}, {chain1_residues} | {chain2_residues}\n"
        
        output += "\n# 清理选择\nselect none\n"
        return output

class IndexFileGenerator:
    """索引文件生成器"""
    
    def __init__(self, structure_parser: StructureParser, all_results: Dict, dimer_results: Dict):
        self.parser = structure_parser
        self.all_results = all_results
        self.dimer_results = dimer_results
    
    def generate_ndx_commands(self) -> str:
        """生成创建索引文件的GROMACS命令"""
        commands = []
        
        # 获取所有链
        chains = list(self.parser.chains.keys())
        
        # 基础组索引映射
        base_group_map = {
            'pro': 1,      # Protein
            'pro-h': 2,    # Protein-H
            'ca': 3,       # C-alpha
            'bb': 4,       # Backbone
            'mc': 5,       # MainChain
            'sc': 8        # SideChain
        }
        
        # 当前组计数器，从10/17开始（默认组0-9无溶液/0-16全体系）
        group_counter = 10
        
        # 1. 为每条链创建单链分组
        for chain in chains:
            # 获取链的原子范围
            residues = self.parser.get_chain_residues(chain)
            if not residues:
                continue
                
            start_atom = residues[0]['atom_seq']
            end_atom = residues[-1]['atom_seq']
            
            # 为每条链创建各种分组
            for suffix, base_group in base_group_map.items():
                # 选择基础组和原子范围
                commands.append(f"{base_group} & a {start_atom}-{end_atom}")
                # 重命名新组
                commands.append(f"name {group_counter} {chain}_{suffix}")
                group_counter += 1
        
        # 2. 为每个单链片段创建分组
        for fragment_name, results in self.all_results.items():
            for result in results:
                chain = result['chain']
                start_atom = result['start']
                end_atom = result['end']
                fragment_id = result['fragment_id']
                
                # 为每个片段创建各种分组
                for suffix, base_group in base_group_map.items():
                    # 选择基础组和原子范围
                    commands.append(f"{base_group} & a {start_atom}-{end_atom}")
                    # 重命名新组
                    commands.append(f"name {group_counter} {fragment_id}_{suffix}")
                    group_counter += 1
        
        # 3. 为每个二聚体创建分组
        for dimer in self.parser.dimers:
            dimer_name = dimer['name']
            chain1 = dimer['chain1']
            chain2 = dimer['chain2']
            
            # 获取两条链的原子范围
            residues1 = self.parser.get_chain_residues(chain1)
            residues2 = self.parser.get_chain_residues(chain2)
            
            if not residues1 or not residues2:
                continue
                
            start_atom1 = residues1[0]['atom_seq']
            end_atom1 = residues1[-1]['atom_seq']
            start_atom2 = residues2[0]['atom_seq']
            end_atom2 = residues2[-1]['atom_seq']
            
            # 为每个二聚体创建各种分组
            for suffix, base_group in base_group_map.items():
                # 选择基础组和原子范围
                commands.append(f"a {start_atom1}-{end_atom1} | a {start_atom2}-{end_atom2} & {base_group}")
                # 重命名新组
                commands.append(f"name {group_counter} {dimer_name}_{suffix}")
                group_counter += 1
        
        # 4. 为每个二聚体片段创建分组
        for dimer_name, fragments in self.dimer_results.items():
            for fragment in fragments:
                fragment_id = fragment['fragment_id']
                start_atom1 = fragment['chain1_start']
                end_atom1 = fragment['chain1_end']
                start_atom2 = fragment['chain2_start']
                end_atom2 = fragment['chain2_end']
                
                # 为每个二聚体片段创建各种分组
                for suffix, base_group in base_group_map.items():
                    # 选择基础组和原子范围
                    commands.append(f"a {start_atom1}-{end_atom1} | a {start_atom2}-{end_atom2} & {base_group}")
                    # 重命名新组
                    commands.append(f"name {group_counter} {fragment_id}_{suffix}")
                    group_counter += 1
        
        # 保存命令
        commands.append('q\n')
        
        return '\n'.join(commands)
    
    def create_ndx_file(self, output_ndx: str = "datanalyse.ndx", execute: bool = False) -> str:
        """创建索引文件"""
        
        # 生成GROMACS命令
        commands = self.generate_ndx_commands()
        
        if execute:
            # 直接执行GROMACS命令
            try:
                # 创建临时命令文件
                with open("temp_ndx_commands.txt", "w") as f:
                    f.write(commands)
                
                # 执行gmx make_ndx命令，不指定输入索引文件
                cmd = [
                    "gmx", "make_ndx",
                    "-f", self.parser.filename,
                    "-o", output_ndx
                ]
                
                result = subprocess.run(
                    cmd,
                    input=commands,
                    text=True,
                    capture_output=True
                )
                
                # 清理临时文件
                if os.path.exists("temp_ndx_commands.txt"):
                    os.remove("temp_ndx_commands.txt")
                
                if result.returncode == 0:
                    return f"成功创建索引文件: {output_ndx}"
                else:
                    return f"GROMACS命令执行失败: {result.stderr}"
                    
            except Exception as e:
                return f"执行GROMACS命令时出错: {str(e)}"
        else:
            # 将命令写入文件
            script_filename = "create_datanalyse_ndx.sh"
            with open(script_filename, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"# 自动生成的索引文件创建脚本\n")
                f.write(f"# 使用: bash {script_filename}\n\n")
                f.write(f"gmx make_ndx -f {self.parser.filename} -o {output_ndx} << EOF\n")
                f.write(commands)
                f.write("\nEOF\n")
            
            # 设置执行权限
            os.chmod(script_filename, 0o755)
            return f"已生成执行脚本: {script_filename}\n运行: bash {script_filename}"

def read_sequence_definitions(filename: str) -> List[Tuple]:
    """读取用户定义的序列文件"""
    definitions = []
    
    with open(filename, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split('/')
            if len(parts) < 4:
                print(f"警告: 第{line_num}行格式错误，跳过: {line}")
                continue
            
            front_seq = parts[0].strip().upper()
            target_seq = parts[1].strip().upper()
            back_seq = parts[2].strip().upper()
            fragment_name = parts[3].strip()
            color = parts[4].strip() if len(parts) > 4 else ""
            
            definitions.append((front_seq, target_seq, back_seq, fragment_name, color))
    
    return definitions

def main():
    parser = argparse.ArgumentParser(description='从PDB/GRO文件中搜索指定序列片段')
    parser.add_argument('-f', '--file', required=True, help='输入PDB/GRO文件')
    parser.add_argument('-s', '--standard', required=True, help='用户定义序列文件')
    parser.add_argument('-o', '--output', help='输出文件路径 (默认: 输入文件目录/fragment_results.txt)')
    parser.add_argument('-pm', '--pymol_commands', nargs='?', const='', 
                       help='生成PyMOL命令。如果指定对象名则使用，否则使用输入文件名')
    parser.add_argument('-opm', '--output_pymol_commands', 
                       help='PyMOL命令输出文件路径 (仅在-pm存在时使用)')
    parser.add_argument('-ndx', '--make_ndx', nargs='?', const='script', 
                       help='创建datanalyse.ndx文件。如果指定"run"则直接执行，否则生成脚本文件')
    parser.add_argument('-ondx', '--output_ndx', 
                       help='输出索引文件路径 (默认: 输入文件目录/datanalyse.ndx)')
    
    args = parser.parse_args()
    
    # 检查输入文件
    if not os.path.exists(args.file):
        print(f"错误: 输入文件不存在: {args.file}")
        sys.exit(1)
    
    if not os.path.exists(args.standard):
        print(f"错误: 序列定义文件不存在: {args.standard}")
        sys.exit(1)
    
    # 根据文件扩展名选择解析器
    file_ext = os.path.splitext(args.file)[1].lower()
    is_gro = False
    
    if file_ext == '.pdb':
        structure_parser = PDBParser(args.file)
        # 检查PDB文件是否使用了-ndx选项
        if args.make_ndx:
            print("错误: PDB文件不支持-ndx选项")
            sys.exit(1)
    elif file_ext == '.gro':
        structure_parser = GROParser(args.file)
        is_gro = True
    else:
        print(f"错误: 不支持的文件格式: {file_ext}")
        sys.exit(1)
       
    # 统一处理输出路径
    try:
        # 获取输入文件的基础名称（不含扩展名）
        input_basename = os.path.splitext(os.path.basename(args.file))[0]
        
        # 处理主输出文件
        output_path = get_output_path(
            args.file, 
            args.output, 
            "fragment_results.txt", 
            ".txt"
        )
        
        # 处理PyMOL命令输出文件
        pymol_output_path = None
        if args.pymol_commands is not None:
            pymol_output_path = get_output_path(
                args.file,
                args.output_pymol_commands,
                f"{input_basename}_pymol.txt",
                ".txt"
            )
        
        # 处理索引文件输出路径
        ndx_output_path = None
        if args.make_ndx:
            if args.output_ndx:
                # 用户指定了输出路径
                ndx_output_path = get_output_path(
                    args.file,
                    args.output_ndx,
                    "datanalyse.ndx",
                    ".ndx"
                )
            else:
                # 使用默认路径（输入文件目录）
                input_dir = os.path.dirname(args.file)
                if input_dir == '':
                    input_dir = '.'
                ndx_output_path = os.path.join(input_dir, "datanalyse.ndx")
                
    except ValueError as e:
        print(f"错误: {e}")
        sys.exit(1)
    
    # 解析结构文件
    print("正在解析结构文件...")
    structure_parser.parse()
    print(f"找到 {len(structure_parser.chains)} 条链: {list(structure_parser.chains.keys())}")
    
    # 读取序列定义
    print("正在读取序列定义...")
    sequence_definitions = read_sequence_definitions(args.standard)
    print(f"读取到 {len(sequence_definitions)} 个序列定义")
    
    # 搜索序列
    print("正在搜索序列...")
    searcher = SequenceSearcher(structure_parser)
    all_results = {}
    
    for i, (front_seq, target_seq, back_seq, fragment_name, color) in enumerate(sequence_definitions, 1):
        print(f"搜索序列 {i}/{len(sequence_definitions)}: {fragment_name}")
        results = searcher.search_sequence(front_seq, target_seq, back_seq, fragment_name, color)
        all_results[fragment_name] = results
    
    # 搜索二聚体片段
    print("正在搜索二聚体片段...")
    dimer_results = searcher.search_dimer_fragments()
    print(f"找到 {len(dimer_results)} 个二聚体的片段")
    
    # 生成输出
    output_generator = OutputGenerator(structure_parser, is_gro)
    
    output_content = output_generator.generate_statistics(all_results, dimer_results)
    output_content += output_generator.generate_fragment_list(all_results, dimer_results)
    
    # 如果指定了PyMOL选项，生成PyMOL命令
    if args.pymol_commands is not None:
        # 确定对象名：如果用户指定了对象名则使用，否则使用输入文件名
        object_name = args.pymol_commands if args.pymol_commands != '' else input_basename
        pymol_commands = output_generator.generate_pymol_commands(all_results, dimer_results, object_name)
        output_content += pymol_commands
        
        # 写入PyMOL命令文件
        with open(pymol_output_path, 'w') as f:
            f.write(pymol_commands)
        print(f"PyMOL命令已写入: {pymol_output_path}")
    
    # 写入主输出文件
    with open(output_path, 'w') as f:
        f.write(output_content)
    
    print(f"结果已写入: {output_path}")
    
    # 在控制台也显示统计结果
    print("\n" + output_generator.generate_statistics(all_results, dimer_results))
    
    # 如果指定了创建索引文件
    if args.make_ndx:
        print("\n正在生成索引文件...")
        ndx_generator = IndexFileGenerator(structure_parser, all_results, dimer_results)
        
        execute = (args.make_ndx == 'run')
        result = ndx_generator.create_ndx_file(
            output_ndx=ndx_output_path,
            execute=execute
        )
        
        print(result)

if __name__ == "__main__":
    main()

