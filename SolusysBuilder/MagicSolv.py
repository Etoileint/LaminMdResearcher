#!/usr/bin/env python3
"""
GROMACS体系准备自动化脚本
作者：Etoileint
功能：自动处理多链蛋白体系，生成完整的GROMACS模拟文件
"""

import os
import sys
import shutil
import argparse
import subprocess
import glob
import re
from collections import defaultdict, OrderedDict
from itertools import groupby

# ========================= 用户配置区域 =========================
# 力场和水模型选择 (运行gmx pdb2gmx时使用)
FORCEFIELD = 8      # 8: Charmm36
WATER_MODEL = 1     # 1: TIP3P-Charmm36

# 端基处理模式
AUTO_TERMINI = False  # True: 全自动模式; False: 半自动模式(需要用户交互)

# editconf参数
EDITCONF_GROUP = 1   # 默认选择Protein组 (1)

# genion参数
GENION_GROUP = 13    # 默认选择SOL组 (13)

# make_ndx索引文件生成命令
MAKE_NDX_COMMANDS = """1
11
name 17 SOLU
name 18 SOLV
q

"""

# ======================== 全局变量和常量 ========================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PDB2GMX = "gmx pdb2gmx"
MAKE_NDX = "gmx make_ndx"
EDITCONF = "gmx editconf"
SOLVATE = "gmx solvate"
GENION = "gmx genion"
GROMPP = "gmx grompp"

# ========================= 核心功能函数 =========================
def copy_standard_files():
    """复制standard_files文件夹内容到当前目录"""
    std_dir = "standard_files"
    if not os.path.exists(std_dir):
        print(f"错误: 未找到{std_dir}文件夹!")
        sys.exit(1)
    
    print(f"复制{std_dir}中的文件到当前目录...")
    for item in os.listdir(std_dir):
        src = os.path.join(std_dir, item)
        if os.path.isfile(src):
            shutil.copy2(src, ".")

def setup_directories():
    """创建必要的目录结构"""
    dirs = ["single_chains", "chain_types", "toppar"]
    for d in dirs:
        if os.path.exists(d):
            shutil.rmtree(d)
        os.makedirs(d)
        print(f"创建目录: {d}")

def merge_pdb_files(pdb_files, output_file="combined.pdb"):
    """合并多个PDB文件并进行预处理"""
    with open(output_file, 'w') as out:
        for pdb_file in pdb_files:
            with open(pdb_file, 'r') as f:
                # 读取并过滤空行
                lines = [line for line in f if line.strip()]
            
            # 处理文件末尾的END行
            if lines and lines[-1].startswith(('END', 'ENDMDL')):
                end_line = lines.pop()
                # 检查是否需要添加TER行
                if lines and not lines[-1].startswith('TER'):
                    lines.append('TER\n')
            # 写入处理后的内容
            out.writelines(lines)
        # 在合并文件末尾添加END行
        out.write("END\n")
    return output_file

def split_pdb_chains(pdb_file):
    """拆分多链PDB文件为单链文件（正确处理TER和END）"""
    print(f"拆分PDB文件: {os.path.basename(pdb_file)}")
    
    with open(pdb_file, 'r') as f:
        pdb_lines = f.readlines()
    
    chains = []
    current_chain = []
    chain_count = 0
    found_end = False
    
    for line in pdb_lines:
        # 遇到END停止处理
        if line.startswith("END") or line.startswith("ENDMDL"):
            found_end = True
            # 如果当前链非空，添加到链列表
            if current_chain:
                chains.append(current_chain)
                chain_count += 1
            break
        
        # 处理ATOM/HETATM行
        if line.startswith("ATOM") or line.startswith("HETATM"):
            current_chain.append(line)
        
        # 处理TER行
        elif line.startswith("TER"):
            current_chain.append(line)
            chains.append(current_chain)
            chain_count += 1
            current_chain = []
    
    # 处理文件末尾没有TER的情况
    if not found_end and current_chain:
        chains.append(current_chain)
        chain_count += 1
    
    # 确定文件名格式
    num_width = len(str(chain_count))
    chain_format = f"chain_{{:0{num_width}d}}.pdb"
    
    # 写入单链文件
    for i, chain in enumerate(chains):
        chain_file = os.path.join("single_chains", chain_format.format(i+1))
        with open(chain_file, 'w') as f:
            f.writelines(chain)
            # 确保每条链都有END
            if not chain[-1].startswith("END"):
                f.write("END\n")
            else:
                # 如果最后一行已经是END，不需要额外添加
                pass
        print(f"  生成链文件: {os.path.basename(chain_file)}")
    
    return chain_count, num_width

def process_single_chain(chain_file):
    """处理单链PDB文件: pdb2gmx转换"""
    print(f"处理链: {os.path.basename(chain_file)}")
    
    base_name = os.path.splitext(os.path.basename(chain_file))[0]
    gro_file = os.path.join("single_chains", f"{base_name}.gro")
    top_file = os.path.join("single_chains", f"{base_name}.top")
    
    # 构建pdb2gmx命令
    command = f"{PDB2GMX} -f {chain_file} -o {gro_file} -p {top_file} -ignh"
    
    if not AUTO_TERMINI:
        command += " -ter"
    
    # 准备自动输入
    input_str = f"{FORCEFIELD}\n{WATER_MODEL}\n"
    
    if not AUTO_TERMINI:
        print("="*50)
        print(" 进入半自动模式 - 请手动选择端基状态")
        print("="*50)
        # 用户需要手动交互
        subprocess.run(command, shell=True, check=True)
    else:
        # 全自动处理
        try:
            subprocess.run(command, input=input_str.encode(), shell=True, check=True)
        except subprocess.CalledProcessError as e:
            if "termini" in str(e).lower():
                print("="*50)
                print(" 错误: 检测到端基选择问题")
                print(" 请在配置中设置 AUTO_TERMINI = False 并重新运行脚本")
                print("="*50)
                sys.exit(1)
            else:
                raise
    
    # 清理生成的itp文件
    for itp in glob.glob(os.path.join("single_chains", "*.itp")):
        os.remove(itp)
    
    return gro_file, top_file

def merge_gro_files(pdb_name, chain_count):
    """合并所有链的GRO文件"""
    print("合并GRO文件...")
    
    # 获取所有链的gro文件
    gro_files = sorted(glob.glob(os.path.join("single_chains", "chain_*.gro")))
    
    # 计算总原子数
    total_atoms = 0
    chain_data = []
    
    for gro_file in gro_files:
        with open(gro_file, 'r') as f:
            lines = f.readlines()
        
        # 跳过前2行和最后1行
        atom_lines = lines[2:-1]
        chain_data.append(atom_lines)
        total_atoms += len(atom_lines)
    
    # 创建合并后的GRO文件
    merged_gro = f"{pdb_name}.gro"
    with open(merged_gro, 'w') as f:
        # 标题行
        f.write(f"{pdb_name} by Etoileint\n")
        # 原子数行
        f.write(f"{total_atoms}\n")
        # 所有原子坐标
        for data in chain_data:
            f.writelines(data)
        # 盒子尺寸 (临时值，后续会修改)
        f.write("   0.0000   0.0000   0.0000\n")
    
    print(f"生成合并GRO文件: {merged_gro}")
    return merged_gro

def create_main_topology(pdb_name, first_top):
    """创建主拓扑文件topol.top（修复块划分问题）"""
    print("创建主拓扑文件: topol.top")
    
    # 定义要保留的部分及其标识
    keep_sections = {
        "forcefield": "; Include forcefield parameters",
        "water": "; Include water topology",
        "ions": "; Include topology for ions",
        "system": "[ system ]",
        "molecules": "[ molecules ]"
    }
    
    # 读取文件内容并按空行分割成块
    blocks = []
    current_block = []
    with open(first_top, 'r') as f:
        for line in f:
            stripped = line.strip()
            
            # 遇到空行表示块结束
            if stripped == "":
                if current_block:
                    blocks.append(current_block)
                    current_block = []
                continue
            
            current_block.append(line)
        
        # 添加最后一个块
        if current_block:
            blocks.append(current_block)
    
    output_blocks = []
    found_molecules = False
    
    # 处理每个块
    for block in blocks:
        block_header = block[0].strip()
        matched_section = None
        
        # 检查块是否匹配任何需要保留的节
        for section_name, section_header in keep_sections.items():
            if section_header == block_header:
                matched_section = section_name
                break
        
        if matched_section:
            # 处理system节：修改系统名称 - 修复多余标题行问题
            if matched_section == "system":
                new_block = [block[0]]  # 保留节标题 [ system ]
                name_line_found = False
                
                # 只保留注释行，移除非注释内容行
                for line in block[1:]:
                    if line.strip().startswith(';'):  # 保留注释行
                        new_block.append(line)
                    elif line.strip().startswith('name'):  # 处理name行
                        new_block.append(f'name  = "{pdb_name}"\n')
                        name_line_found = True
                
                # 如果没有找到name行，添加新的name行
                if not name_line_found:
                    new_block.append(f'{pdb_name}\n')
                    
                output_blocks.append(new_block)
            
            # 处理molecules节：只保留标题
            elif matched_section == "molecules":
                output_blocks.append([block[0]])
                found_molecules = True
            
            # 处理其他保留节
            else:
                output_blocks.append(block)
    
    # 确保molecules节存在
    if not found_molecules:
        output_blocks.append([keep_sections["molecules"] + "\n"])
    
    # 写入主拓扑文件（块之间用空行分隔）
    with open("topol.top", 'w') as f:
        for i, block in enumerate(output_blocks):
            for line in block:
                f.write(line)
            # 在块之间添加空行（最后一个块除外）
            if i < len(output_blocks) - 1:
                f.write("\n")
    
    return "topol.top"

def classify_chains():
    """根据氨基酸序列分类链类型"""
    print("根据氨基酸序列分类链类型...")
    
    # 获取所有链的GRO文件
    gro_files = sorted(glob.glob(os.path.join("single_chains", "chain_*.gro")))
    chain_types = {}
    type_counter = defaultdict(int)
    type_representatives = {}
    
    # 读取并分类链
    for gro_file in gro_files:
        with open(gro_file, 'r') as f:
            lines = f.readlines()
        
        # 提取残基序列 - 修正后的解析方式
        residues = OrderedDict()
        for line in lines[2:-1]:  # 跳过头两行和最后一行
            # 正确解析GRO文件格式：
            # 第1-5位: 残基编号 (5字符)
            # 第6-10位: 残基名称 (5字符)
            res_id = line[0:5].strip()
            res_name = line[5:10].strip()
            
            # 只保留每个残基的第一次出现
            if res_id not in residues:
                residues[res_id] = res_name
        
        # 创建序列签名
        seq = tuple(residues.values())
        seq_len = len(seq)
        
        # 分类链类型
        if seq in chain_types:
            chain_type = chain_types[seq]
        else:
            # 生成新类型名称
            base_name = f"C{seq_len}"
            suffix = ""
            
            if base_name in type_counter:
                count = type_counter[base_name]
                # 生成字母后缀 (aa, ab, ..., zz)
                first_char = chr(97 + count // 26)
                second_char = chr(97 + count % 26)
                suffix = f"_{first_char}{second_char}"
            
            chain_type = f"{base_name}{suffix}"
            chain_types[seq] = chain_type
            type_counter[base_name] += 1
            type_representatives[chain_type] = gro_file
    
    # 记录链类型到日志文件
    with open("chain_types.log", 'w') as log:
        for i, gro_file in enumerate(gro_files):
            chain_name = os.path.basename(gro_file).replace(".gro", "")
            with open(gro_file, 'r') as f:
                lines = f.readlines()
            
            residues = OrderedDict()
            for line in lines[2:-1]:
                res_id = line[0:5].strip()
                res_name = line[5:10].strip()
                if res_id not in residues:
                    residues[res_id] = res_name
            
            seq = tuple(residues.values())
            chain_type = chain_types[seq]
            log.write(f"{chain_name}: {chain_type}\n")
    
    # 复制代表链到chain_types目录
    for chain_type, gro_file in type_representatives.items():
        # 复制GRO文件
        dest_gro = os.path.join("chain_types", f"{chain_type}.gro")
        shutil.copy2(gro_file, dest_gro)
        
        # 复制对应的TOP文件
        src_top = gro_file.replace(".gro", ".top")
        dest_top = os.path.join("chain_types", f"{chain_type}.top")
        if os.path.exists(src_top):
            shutil.copy2(src_top, dest_top)
        
        print(f"  链类型: {chain_type} -> 代表文件: {os.path.basename(dest_gro)}")
    
    return type_representatives

def process_chain_types():
    """处理每种链类型：生成NDX文件并修改ITP文件"""
    print("处理链类型...")
    
    for top_file in glob.glob(os.path.join("chain_types", "*.top")):
        base_name = os.path.splitext(os.path.basename(top_file))[0]
        gro_file = os.path.join("chain_types", f"{base_name}.gro")
        ndx_file = os.path.join("chain_types", f"{base_name}_bb_sc.ndx")
        itp_file = os.path.join("chain_types", f"{base_name}.itp")
        final_itp = os.path.join("toppar", f"{base_name}.itp")
        
        # 生成NDX文件
        print(f"  为 {base_name} 生成索引文件")
        command = f"{MAKE_NDX} -f {gro_file} -o {ndx_file}"
        subprocess.run(command, input="q\n".encode(), shell=True, check=True)
        
        # 转换TOP为ITP
        convert_top_to_itp(top_file, itp_file, base_name)
        
        # 添加位置限制
        add_posres_to_itp(itp_file, ndx_file, final_itp)
        
        # 清理临时文件
        os.remove(itp_file)
    
    print("链类型处理完成")

def convert_top_to_itp(top_file, itp_file, chain_type):
    """将TOP文件转换为ITP文件格式"""
    # 定义要删除的部分
    exclude_sections = [
        "; Include forcefield parameters",
        "; Include chain topologies",
        "; Include water topology",
        "; Include topology for ions",
        "[ system ]",
        "[ molecules ]"
    ]
    
    # 定义要保留的部分
    keep_sections = [
        "[ moleculetype ]",
        "[ atoms ]",
        "[ bonds ]",
        "[ pairs ]",
        "[ angles ]",
        "[ dihedrals ]",
        "[ cmap ]"
    ]
    
    current_section = None
    output_lines = []
    skip_section = False
    
    with open(top_file, 'r') as f:
        for line in f:
            stripped = line.strip()
            
            # 检查排除部分
            if any(stripped.startswith(sec) for sec in exclude_sections):
                current_section = stripped
                skip_section = True
                continue
            
            # 检查保留部分
            if any(stripped.startswith(sec) for sec in keep_sections):
                current_section = stripped
                skip_section = False
                output_lines.append(line)
                
                # 修改moleculetype名称
                if stripped == "[ moleculetype ]":
                    # 下一行是注释
                    next_line = next(f)
                    output_lines.append(next_line)
                    # 名称行
                    name_line = next(f)
                    parts = name_line.split()
                    parts[0] = chain_type
                    output_lines.append(" ".join(parts) + "\n")
                continue
            
            # 处理内容行
            if not skip_section:
                output_lines.append(line)
    
    # 写入ITP文件
    with open(itp_file, 'w') as f:
        f.writelines(output_lines)

def add_posres_to_itp(itp_file, ndx_file, output_file):
    """为ITP文件添加位置限制（集成add_posres.py功能）"""
    # 读取NDX文件
    backbone_atoms = []
    sidechain_h_atoms = []
    current_section = None
    
    with open(ndx_file, 'r') as f:
        for line in f:
            stripped = line.strip()
            
            # 检查节标题
            if stripped == '[ Backbone ]':
                current_section = 'Backbone'
                continue
            elif stripped == '[ SideChain-H ]':
                current_section = 'SideChain-H'
                continue
            elif stripped.startswith('[') and stripped.endswith(']'):
                current_section = None
                continue
            
            # 收集原子序号
            if current_section == 'Backbone' and stripped:
                backbone_atoms.extend(map(int, stripped.split()))
            elif current_section == 'SideChain-H' and stripped:
                sidechain_h_atoms.extend(map(int, stripped.split()))
    
    # 检查重复原子
    common = set(backbone_atoms) & set(sidechain_h_atoms)
    if common:
        raise ValueError(f"错误: 发现重复原子索引: {sorted(common)}")
    
    # 生成位置限制行
    all_atoms = sorted(set(backbone_atoms + sidechain_h_atoms))
    posres_lines = ["[ position_restraints ]"]
    for atom in all_atoms:
        fc_type = "POSRES_FC_BB" if atom in backbone_atoms else "POSRES_FC_SC"
        posres_lines.append(f"{atom:6d}     1    {fc_type:15s} {fc_type:15s} {fc_type:15s}")
    
    # 处理ITP文件
    posres_block = ["#ifdef POSRES"] + posres_lines + ["#endif"]
    
    with open(itp_file, 'r') as f:
        itp_lines = f.readlines()
    
    # 查找现有POSRES块
    posres_start = -1
    posres_end = -1
    for i, line in enumerate(itp_lines):
        stripped = line.strip()
        if stripped == "#ifdef POSRES":
            posres_start = i
        elif stripped == "#endif" and posres_start != -1:
            posres_end = i
            break
    
    # 插入或替换POSRES块
    if posres_start != -1 and posres_end != -1:
        new_content = itp_lines[:posres_start] + posres_block + itp_lines[posres_end+1:]
    else:
        new_content = itp_lines + ["\n"] + posres_block
    
    # 写入输出文件
    with open(output_file, 'w') as f:
        f.writelines([line if line.endswith('\n') else line + '\n' for line in new_content])

def update_main_topology():
    """更新主拓扑文件：包含链类型ITP文件并构建[molecules]部分"""
    print("更新主拓扑文件...")
    
    # 读取链类型日志
    chain_mapping = {}
    with open("chain_types.log", 'r') as f:
        for line in f:
            if ":" in line:
                chain, ctype = line.split(":", 1)
                chain_mapping[chain.strip()] = ctype.strip()
    
    # 按原始顺序收集链类型
    chain_types = []
    for chain in sorted(chain_mapping.keys()):
        chain_types.append(chain_mapping[chain])
    
    # 合并连续相同的链类型
    merged_types = []
    for key, group in groupby(chain_types):
        count = len(list(group))
        merged_types.append((key, count))
    
    # 在拓扑文件中包含ITP文件
    itp_files = sorted(glob.glob(os.path.join("toppar", "*.itp")))
    
    with open("topol.top", 'r') as f:
        top_lines = f.readlines()
    
    new_top_lines = []
    chain_topology_found = False
    molecules_found = False
    forcefield_block_end = -1  # 记录forcefield块结束位置
    
    # 第一遍：检查是否存在chain topologies块，并记录forcefield块结束位置
    for i, line in enumerate(top_lines):
        stripped = line.strip()
        if stripped == "; Include chain topologies":
            chain_topology_found = True
        if stripped == "; Include forcefield parameters":
            # 找到forcefield块结束位置（下一个空行）
            j = i
            while j < len(top_lines) and top_lines[j].strip() != "":
                j += 1
            forcefield_block_end = j  # 记录空行位置
    
    # 如果未找到chain topologies块，且找到forcefield块结束位置
    if not chain_topology_found and forcefield_block_end != -1:
        # 在forcefield块结束后添加chain topologies块
        top_lines.insert(forcefield_block_end + 1, "; Include chain topologies\n")
    
    # 第二遍：处理更新后的内容
    for line in top_lines:
        stripped = line.strip()
        
        # 在链拓扑部分后添加ITP包含
        if stripped == "; Include chain topologies":
            new_top_lines.append(line)
            chain_topology_found = True
            for itp in itp_files:
                itp_base = os.path.basename(itp)
                new_top_lines.append(f'#include "toppar/{itp_base}"\n')
            # 确保块结束后有一个空行
            new_top_lines.append("\n")
            continue
        
        # 替换[molecules]部分
        if stripped == "[ molecules ]":
            new_top_lines.append(line)
            new_top_lines.append("; Compound\t#mols\n")
            molecules_found = True
            for ctype, count in merged_types:
                new_top_lines.append(f"{ctype}\t\t{count}\n")
            continue
        
        # 跳过旧的molecules内容
        if molecules_found and not stripped.startswith('[') and stripped != "":
            continue
        
        new_top_lines.append(line)
    
    # 写入更新的拓扑文件
    with open("topol.top", 'w') as f:
        f.writelines(new_top_lines)
    
    print("主拓扑文件更新完成")

def run_editconf(gro_file, box_distance=None, box_size=None):
    """运行editconf创建模拟盒子"""
    print("创建模拟盒子...")
    box_gro = gro_file.replace(".gro", "_box.gro")
    
    # 构建基础命令
    command = f"{EDITCONF} -f {gro_file} -o {box_gro} -bt triclinic -princ -c"
    
    # 根据参数选择盒子生成方式
    if box_distance is not None:
        # 使用边界距离模式
        command += f" -d {box_distance}"
        print(f"  使用边界距离模式: {box_distance} nm")
    elif box_size is not None:
        # 使用直接指定盒子尺寸模式
        x, y, z = box_size
        command += f" -box {x} {y} {z}"
        print(f"  使用盒子尺寸模式: {x} x {y} x {z} nm")
    else:
        raise ValueError("必须指定盒子生成方式：-bd/--box_distance 或 -bs/--box_size")
    
    # 使用预定义的组选择
    subprocess.run(command, input=f"{EDITCONF_GROUP}\n".encode(), shell=True, check=True)
    return box_gro

def run_solvate(box_gro):
    """运行solvate添加水分子"""
    print("添加水分子...")
    sol_gro = box_gro.replace("_box.gro", "_sol.gro")
    
    command = f"{SOLVATE} -cp {box_gro} -o {sol_gro} -p topol.top"
    subprocess.run(command, shell=True, check=True)
    return sol_gro

def add_ions(sol_gro):
    """添加离子平衡电荷"""
    print("添加离子...")
    
    # 生成离子tpr文件
    grompp_cmd = f"{GROMPP} -f ions.mdp -c {sol_gro} -p topol.top -o ions.tpr -maxwarn 1"
    subprocess.run(grompp_cmd, shell=True, check=True)
    
    # 添加离子
    genion_cmd = f"{GENION} -s ions.tpr -p topol.top -o step3_input.gro -neutral -conc 0.15"
    subprocess.run(genion_cmd, input=f"{GENION_GROUP}\n".encode(), shell=True, check=True)
    
    return "step3_input.gro"

def create_final_index(final_gro):
    """创建最终的索引文件"""
    print("创建最终索引文件...")
    command = f"{MAKE_NDX} -f {final_gro} -o index.ndx"
    subprocess.run(command, input=MAKE_NDX_COMMANDS.encode(), shell=True, check=True)

# ========================= 主程序 =========================
def main():
    # 解析命令行参数
    parser = argparse.ArgumentParser(description="GROMACS多链蛋白体系准备脚本")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-p", "--pdb", help="输入PDB文件")
    group.add_argument("-d", "--directory", help="包含单链PDB文件的目录")
    
    # 盒子参数组（互斥选择）
    box_group = parser.add_mutually_exclusive_group(required=True)
    box_group.add_argument("-bd", "--box_distance", type=float, 
                         help="使用边界距离生成盒子（单位：nm）")
    box_group.add_argument("-bs", "--box_size", nargs=3, type=float, 
                         help="直接指定盒子尺寸（单位：nm），格式：x y z")    
        
    args = parser.parse_args()
    
    try:
        # 步骤1: 准备目录和标准文件
        copy_standard_files()
        setup_directories()
        
        # 步骤2: 处理输入文件
        num_width = None
        if args.pdb:
            pdb_name = os.path.splitext(os.path.basename(args.pdb))[0]
            chain_count, num_width = split_pdb_chains(args.pdb)
        else:
            # 使用传入的目录名作为系统名
            dir_name = os.path.basename(os.path.normpath(args.directory))
            if not dir_name or dir_name == '.':
                # 如果目录名为空或是当前目录，使用默认名称
                pdb_name = "multi_chain_system"
            else:
                pdb_name = dir_name

            chain_files = sorted(glob.glob(os.path.join(args.directory, "*.pdb")))
            chain_count = len(chain_files)
            
            # 合并预处理所有PDB文件
            combined_pdb = os.path.join(args.directory, "combined.pdb")
            with open(combined_pdb, 'w') as out:
                for pdb_file in chain_files:
                    with open(pdb_file, 'r') as f:
                        # 读取并过滤空行
                        lines = [line for line in f if line.strip()]
                    
                    # 处理文件末尾的END行
                    if lines:
                        # 检查最后一行是否是END/ENDMDL
                        if lines[-1].startswith(('END', 'ENDMDL')):
                            end_line = lines.pop()
                            # 检查是否需要添加TER行
                            if lines and not lines[-1].startswith('TER'):
                                lines.append('TER\n')
                    
                    # 写入处理后的内容
                    out.writelines(lines)
                
                # 在合并文件末尾添加END行
                out.write("END\n")
            
            # 拆分合并后的文件
            chain_count, num_width = split_pdb_chains(combined_pdb)
            
            # 删除临时合并文件
            os.remove(combined_pdb)
        
        # 步骤3: 处理每条链 - 使用正确的文件名格式
        first_top = None
        for i in range(1, chain_count + 1):
            # 使用num_width确保文件名格式一致
            chain_file = os.path.join("single_chains", f"chain_{i:0{num_width}d}.pdb")
            gro, top = process_single_chain(chain_file)
            if i == 1:
                first_top = top
        
        # 步骤4: 创建主GRO和TOP文件
        merged_gro = merge_gro_files(pdb_name, chain_count)
        main_top = create_main_topology(pdb_name, first_top)
        
        # 步骤5: 分类链类型
        type_reps = classify_chains()
        
        # 步骤6: 处理链类型
        process_chain_types()
        
        # 步骤7: 更新主拓扑
        update_main_topology()
        
        # 步骤8: 体系准备流程
        box_gro = run_editconf(merged_gro, args.box_distance, args.box_size)
        sol_gro = run_solvate(box_gro)
        final_gro = add_ions(sol_gro)
        create_final_index(final_gro)
        
        print("\n" + "="*50)
        print("体系准备完成！")
        print(f"最终结构文件: {final_gro}")
        print(f"拓扑文件: topol.top")
        print(f"索引文件: index.ndx")
        print("="*50)
    
    except Exception as e:
        print(f"\n错误: {str(e)}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()

