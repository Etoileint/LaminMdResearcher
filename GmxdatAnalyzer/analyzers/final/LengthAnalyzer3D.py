'''
参数说明：
输入输出：
-i md_trajectory.xvg：指定GROMACS生成的轨迹文件
-o 后接三个输出文件：长度数据、时程图、动画文件（必须按顺序）

核心参数：
-w/--window_size：滑动窗口包含的原子数
-cp/--curve_points：控制曲线精度的离散点数（特征点间插值数，值越大曲线越平滑，默认为20）
-cc/--curve_coords：保存曲线坐标

动画优化：
-a/--animate：启用3D动画生成（会显著增加计算时间）
-fs/--frame_scale：帧采样率（10表示只处理1/10的帧）

示例：
python LengthAnalyzer3Dv4.1.py -i "C:/Users/Magichoco/Desktop/backbone_coords.xvg" -o "C:/Users/Magichoco/Desktop/len.xvg" "C:/Users/Magichoco/Desktop/plot.svg" "C:/Users/Magichoco/Desktop/anim.html" -w 12 -a -fs 10 -cp 10

python LengthAnalyzer3Dv4.1.py -i backbone_coords_cA.xvg -o len_cA.xvg plot_cA.svg anim_cA.html -w 12 -a -fs 10 -cp 10 -cc curve_cord_cA.npz
python LengthAnalyzer3Dv4.1.py -i backbone_coords_cB.xvg -o len_cB.xvg plot_cB.svg anim_cB.html -w 12 -a -fs 10 -cp 10 -cc curve_cord_cB.npz

'''


import argparse
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import plotly.graph_objects as go
from scipy.interpolate import CubicSpline
import os

def parse_xvg(file_path):
    # ...保持原有parse_xvg函数不变...  
    times = []
    all_coords = []
    num_atoms = None
    
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(('@', '#')):
                continue
            
            parts = list(map(float, line.split()))
            if len(parts) < 4:
                continue
            
            time = parts[0]
            coords = np.array(parts[1:])
            
            if (len(coords) % 3) != 0:
                raise ValueError("Invalid coordinate count in line")
            
            current_atoms = len(coords) // 3
            if num_atoms is None:
                num_atoms = current_atoms
            elif current_atoms != num_atoms:
                raise ValueError("Inconsistent number of atoms between frames")
            
            times.append(time)
            all_coords.append(coords.reshape(-1, 3))
    
    return times, all_coords, num_atoms

def generate_windows(num_atoms, window_size):
    if window_size > num_atoms:
        raise ValueError("Window size cannot be larger than total number of atoms")
    return [(i, i + window_size - 1) for i in range(num_atoms - window_size + 1)]

def main():
    parser = argparse.ArgumentParser(description='Analyze protein length using sliding window approach')
    parser.add_argument('-i', '--input', required=True, help='Input .xvg file')
    parser.add_argument('-o', '--output', nargs='*', default=[], help='Output files (length data, plot, animation)')
    parser.add_argument('-a', '--animate', action='store_true', help='Generate 3D animation')
    parser.add_argument('-w', '--window_size', type=int, required=True,
                       help='Number of atoms in each sliding window')
    parser.add_argument('-cp', '--curve_points', type=int, default=20,
                       help='Number of points between adjacent feature points (default: 20)')
    parser.add_argument('-fs', '--frame_scale', type=int, default=1,
                       help='Frame scaling factor (e.g. 10 reduces frames to 1/10)')
    # 新增参数：保存曲线坐标
    parser.add_argument('-cc', '--curve_coords', help='Save curve coordinates to .npz file')  # <--- 新增参数
    
    args = parser.parse_args()

    default_outputs = ['protein_length.xvg', 'protein_length.svg', 'protein_animation.html']
    outputs = args.output + default_outputs[len(args.output):]
    length_path = outputs[0]
    plot_path = outputs[1]
    animation_path = outputs[2] if args.animate else None

    times, all_coords, num_atoms = parse_xvg(args.input)
    print(f"Total atoms: {num_atoms}, Total frames: {len(times)}")
    print(f"Window size: {args.window_size}")
    
    try:
        windows = generate_windows(num_atoms, args.window_size)
    except ValueError as e:
        print(f"Error: {e}")
        return

    print(f"Number of windows: {len(windows)}")

    lengths = []
    animation_data = [] if args.animate else None
    # 新增曲线坐标存储容器
    curve_coords_list = [] if args.curve_coords else None  # <--- 新增变量

    for frame_idx in tqdm(range(len(times)), desc='Processing frames'):
        coords = all_coords[frame_idx]
        start_points = []
        end_points = []
        mid_points = []
        frame_curve = None

        # Process each window
        for start, end in windows:
            seg_coords = coords[start:end+1]
            centroid = np.mean(seg_coords, axis=0)
            centered = seg_coords - centroid
            
            cov = np.cov(centered, rowvar=False)
            eigenvalues, eigenvectors = np.linalg.eigh(cov)
            principal_axis = eigenvectors[:, np.argmax(eigenvalues)]
            
            projections = np.dot(centered, principal_axis)
            min_p = np.min(projections)
            max_p = np.max(projections)
            start_point = centroid + principal_axis * min_p
            end_point = centroid + principal_axis * max_p
            mid_point = (start_point + end_point) / 2
            
            start_points.append(start_point)
            end_points.append(end_point)
            mid_points.append(mid_point)
        
        # Create curve from feature points
        if len(mid_points) == 0:
            curve_length = 0.0
            curve = np.empty((0,3))
        else:
            # Construct feature points (first start, mids, last end)
            first_start = start_points[0]
            first_end = end_points[0]
            last_start = start_points[-1]
            last_end = end_points[-1]
            feature_points = [first_start] + mid_points + [last_end]
            feature_points = np.array(feature_points)
            
            if len(feature_points) < 2:
                curve_length = 0.0
                curve = np.empty((0,3))
            else:
                # Cubic spline interpolation
                k = len(feature_points)
                t = np.linspace(0, 1, k)
                try:
                    cs_x = CubicSpline(t, feature_points[:, 0])
                    cs_y = CubicSpline(t, feature_points[:, 1])
                    cs_z = CubicSpline(t, feature_points[:, 2])
                except:
                    curve_length = 0.0
                    curve = np.empty((0,3))
                else:
                    # Generate interpolation points
                    n = args.curve_points
                    total_points = (k-1)*n + 1
                    t_curve = np.linspace(0, 1, total_points)
                    curve = np.column_stack([cs_x(t_curve), cs_y(t_curve), cs_z(t_curve)])
                    
                    # Calculate total curve length
                    if len(curve) < 2:
                        curve_length = 0.0
                    else:
                        curve_length = np.sum(np.linalg.norm(np.diff(curve, axis=0), axis=1))
        
        lengths.append(curve_length)
        
        # 保存曲线坐标数据
        if curve_coords_list is not None:  # <--- 新增保存逻辑
            curve_coords_list.append(curve)
        
        if args.animate:
            animation_data.append({
                'coords': coords,
                'curve': curve
            })

    # Save length data
    with open(length_path, 'w') as f:
        f.write("# Time (ps)\tLength (nm)\n")
        for t, l in zip(times, lengths):
            f.write(f"{t}\t{l}\n")

    # Generate plot
    plt.figure(figsize=(10, 6))
    plt.plot(times, lengths)
    plt.xlabel('Time (ps)')
    plt.ylabel('Length (nm)')
    plt.title('Protein Length Over Time')
    plt.savefig(plot_path, format='svg')
    plt.close()

    # 保存曲线坐标到文件
    if args.curve_coords:  # <--- 新增保存功能
        print(f"Saving curve coordinates to {args.curve_coords}")
        np.savez(args.curve_coords, curves=curve_coords_list)

   # Generate animation
    if args.animate:
        print("Generating animation...")
        
        if args.frame_scale > 1:
            animation_data = animation_data[::args.frame_scale]
            times_animation = times[::args.frame_scale]
        else:
            times_animation = times

        # Calculate global coordinate range
        all_coords_list = [data['coords'] for data in animation_data]
        all_curves_list = [data['curve'] for data in animation_data]
        all_points = np.concatenate(all_coords_list + all_curves_list, axis=0)
        
        x_min, x_max = np.min(all_points[:, 0]), np.max(all_points[:, 0])
        y_min, y_max = np.min(all_points[:, 1]), np.max(all_points[:, 1])
        z_min, z_max = np.min(all_points[:, 2]), np.max(all_points[:, 2])
        
        padding = 1.0
        x_range = [x_min - padding, x_max + padding]
        y_range = [y_min - padding, y_max + padding]
        z_range = [z_min - padding, z_max + padding]
        
        # Calculate aspect ratio
        x_len = x_range[1] - x_range[0]
        y_len = y_range[1] - y_range[0]
        z_len = z_range[1] - z_range[0]
        max_len = max(x_len, y_len, z_len)
        aspect_ratio = {
            'x': x_len/max_len,
            'y': y_len/max_len,
            'z': z_len/max_len
        }

        fig = go.Figure()

        # Initial frame
        initial = animation_data[0]
        # Atoms
        fig.add_trace(go.Scatter3d(
            x=initial['coords'][:, 0],
            y=initial['coords'][:, 1],
            z=initial['coords'][:, 2],
            mode='markers',
            marker=dict(size=2, color='blue', symbol='circle'),
            name='Atoms'
        ))
        # Curve
        fig.add_trace(go.Scatter3d(
            x=initial['curve'][:, 0],
            y=initial['curve'][:, 1],
            z=initial['curve'][:, 2],
            mode='lines',
            line=dict(color='red', width=6),
            name='Curve'
        ))

        # Create animation frames
        frames = []
        for idx, data in enumerate(tqdm(animation_data, desc='Preparing animation')):
            traces = [
                # Atoms
                go.Scatter3d(
                    x=data['coords'][:, 0],
                    y=data['coords'][:, 1],
                    z=data['coords'][:, 2],
                    mode='markers',
                    marker=dict(size=1, color='blue', symbol='circle'),
                    showlegend=False
                ),
                # Curve
                go.Scatter3d(
                    x=data['curve'][:, 0],
                    y=data['curve'][:, 1],
                    z=data['curve'][:, 2],
                    mode='lines',
                    line=dict(color='red', width=6),
                    showlegend=False
                )
            ]
            frames.append(go.Frame(
                data=traces,
                name=str(idx),
                traces=[0, 1]  # 对应两个trace的顺序
            ))

        # Set layout
        fig.update_layout(
            title='Protein Backbone with Interpolated Curve',
            scene=dict(
                xaxis=dict(
                    showgrid=True,
                    gridcolor='rgba(150,150,150,0.5)',
                    range=x_range,
                    title='X (nm)'
                ),
                yaxis=dict(
                    showgrid=True,
                    gridcolor='rgba(150,150,150,0.5)',
                    range=y_range,
                    title='Y (nm)'
                ),
                zaxis=dict(
                    showgrid=True,
                    gridcolor='rgba(150,150,150,0.5)',
                    range=z_range,
                    title='Z (nm)'
                ),
                aspectmode='manual',
                aspectratio=aspect_ratio,
                bgcolor='white'
            ),
            updatemenus=[{
                'type': 'buttons',
                'buttons': [
                    {
                        'label': '▶ Play',
                        'method': 'animate',
                        'args': [None, {
                            'frame': {'duration': 50, 'redraw': True},
                            'fromcurrent': True,
                            'mode': 'immediate',
                            'transition': {'duration': 0}
                        }]
                    },
                    {
                        'label': '⏸ Pause',
                        'method': 'animate',
                        'args': [[None], {
                            'frame': {'duration': 0, 'redraw': False},
                            'mode': 'immediate',
                            'transition': {'duration': 0}
                        }]
                    }
                ],
                'x': 0.1,
                'xanchor': 'right',
                'y': 0,
                'yanchor': 'top'
            }],
            sliders=[{
                'steps': [{
                    'args': [[f.name], {
                        'frame': {'duration': 50, 'redraw': True},
                        'mode': 'immediate',
                        'transition': {'duration': 0}
                    }],
                    'label': f'{t:.1f} ps',
                    'method': 'animate'
                } for t, f in zip(times_animation, frames)],
                'transition': {'duration': 0},
                'x': 0.1,
                'len': 0.9,
                'currentvalue': {'prefix': 'Time: ', 'suffix': ' ps'}
            }]
        )

        fig.frames = frames
        fig.write_html(animation_path)
        print(f"Animation saved to {animation_path}")

if __name__ == '__main__':
    main()
