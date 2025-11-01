#!/usr/bin/env python3
"""
提取.ndx文件中的组名并生成序号
"""

import argparse
import os
import re
import sys

def extract_group_names(ndx_file):
    """
    从.ndx文件中提取所有组名
    
    Args:
        ndx_file: .ndx文件路径
        
    Returns:
        list: 组名列表
    """
    group_names = []
    pattern = r'^\s*\[\s*(.+?)\s*\]\s*$'  # 匹配 [ Group Name ] 格式
    
    try:
        with open(ndx_file, 'r', encoding='utf-8') as f:
            for line in f:
                match = re.match(pattern, line)
                if match:
                    group_name = match.group(1).strip()
                    group_names.append(group_name)
    except Exception as e:
        print(f"读取文件时出错: {e}")
        sys.exit(1)
        
    return group_names

def write_group_names_to_file(group_names, output_file):
    """
    将组名和序号写入输出文件
    
    Args:
        group_names: 组名列表
        output_file: 输出文件路径
    """
    try:
        with open(output_file, 'w', encoding='utf-8') as f:
            for i, group_name in enumerate(group_names):
                f.write(f"{i} {group_name}\n")
        print(f"成功生成输出文件: {output_file}")
        print(f"共找到 {len(group_names)} 个组")
    except Exception as e:
        print(f"写入文件时出错: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description='提取.ndx文件中的组名并生成序号',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例:
  python ndxlister.py -i index.ndx
  python ndxlister.py -i index.ndx -o groups.txt
        '''
    )
    
    parser.add_argument('-i', '--input', required=True,
                       help='输入的.ndx文件路径')
    parser.add_argument('-o', '--output',
                       help='输出的txt文件路径（可选）')
    
    args = parser.parse_args()
    
    # 检查输入文件是否存在
    if not os.path.isfile(args.input):
        print(f"错误: 输入文件 '{args.input}' 不存在")
        sys.exit(1)
    
    # 生成默认输出文件名
    if args.output is None:
        input_dir = os.path.dirname(args.input)
        input_basename = os.path.basename(args.input)
        name_without_ext = os.path.splitext(input_basename)[0]
        default_output = f"{name_without_ext}_ndxlist.txt"
        args.output = os.path.join(input_dir, default_output) if input_dir else default_output
    
    # 提取组名
    group_names = extract_group_names(args.input)
    
    if not group_names:
        print("警告: 未在文件中找到任何组名")
        print("请确认文件格式正确，组名格式为: [ Group Name ]")
    
    # 写入输出文件
    write_group_names_to_file(group_names, args.output)

if __name__ == "__main__":
    main()

