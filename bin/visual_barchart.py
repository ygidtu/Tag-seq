#!/usr/bin/env python3
"""
Visualize CRISPR off-target sites with coverage depth bars.
Author: Converted from original Perl/Python2 workflow
Description: Generate SVG visualization for off-target sites showing coverage 
             from different strand and primer combinations.
"""

import svgwrite
import sys
import os
import math
from pathlib import Path
from typing import List, Tuple, Dict, Any

# Constants
BOX_WIDTH = 10
BOX_SIZE = 15
V_SPACING = 3

# Color mappings
COLORS = {
    'G': '#F5F500',  # Yellow
    'A': '#FF5454',  # Red
    'T': '#00D118',  # Green
    'C': '#26A8FF',  # Blue
    'N': '#B3B3B3'   # Gray
}

CHAIN_COLORS = {
    'pp': '#70CDE2',  # Light blue - (+) Strand, forward primer
    'pm': '#BE70E2',  # Purple - (+) Strand, reverse primer
    'mp': '#E28570',  # Orange - (-) Strand, forward primer
    'mm': '#94E270'   # Light green - (-) Strand, reverse primer
}

def parse_depth_file(infile: str, strand: str) -> Tuple[List[List[str]], List[int]]:
    """
    Parse depth file containing coverage data for each position.
    
    Args:
        infile: Path to depth file (tab-separated)
        strand: Strand direction (+ or -)
    
    Returns:
        Tuple of (depth_array, sorted_depth_sums)
    """
    depth = []
    depth2 = []
    
    try:
        with open(infile, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:  # Skip empty lines
                    continue
                    
                items = line.split('\t')
                if len(items) < 7:
                    print(f"Warning: Invalid line format in {infile}: {line}")
                    continue
                    
                # Extract coverage values for each chain
                d_arr = [items[3], items[4], items[5], items[6]]
                depth.append(d_arr)
                
                # Calculate sum of all coverages
                pos_sum = sum(int(x) for x in d_arr)
                depth2.append(pos_sum)
    
    except FileNotFoundError:
        print(f"Error: Cannot find depth file {infile}")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading depth file {infile}: {e}")
        sys.exit(1)
    
    # Reverse for negative strand
    if strand == "-":
        depth = depth[::-1]
    
    depth2.sort()
    return depth, depth2

def calculate_y_axis_scale(depth_values: List[int]) -> Tuple[int, int, int]:
    """
    Calculate appropriate Y-axis scale based on depth values.
    
    Args:
        depth_values: Sorted list of depth sums
    
    Returns:
        Tuple of (max_value, unit, num_units)
    """
    if not depth_values:
        return 0, 1, 0
    
    max_val = depth_values[-1]
    
    # Calculate unit size (round to nice numbers)
    if max_val == 0:
        return 0, 1, 1
    
    # Find appropriate step size
    devide6 = max_val / 6
    devide6_len = len(str(int(devide6))) - 1
    base = 10 ** max(0, devide6_len)
    unit = (int(devide6 / base) + 1) * base
    
    # Calculate number of units
    num_units = int(max_val / unit) + 2  # +2 to include zero and one extra
    
    return max_val, unit, num_units

def draw_axes(dwg: svgwrite.Drawing, x_offset: int, y_offset: int, 
              seq_length: int, max_val: int, unit: int, num_units: int) -> None:
    """
    Draw X and Y axes with ticks and labels.
    
    Args:
        dwg: SVG drawing object
        x_offset: X-axis starting offset
        y_offset: Y-axis starting offset
        seq_length: Length of sequence
        max_val: Maximum coverage value
        unit: Unit size for Y-axis
        num_units: Number of units to display
    """
    # Calculate axis positions
    bar_height = 180
    x_axis_start = (x_offset, y_offset + bar_height)
    x_axis_end = (x_offset + seq_length * BOX_SIZE, y_offset + bar_height)
    
    y_axis_start = (x_offset, y_offset)
    y_axis_end = (x_offset, y_offset + bar_height)
    
    # Draw axes
    dwg.add(dwg.line(start=x_axis_start, end=x_axis_end, stroke='black', stroke_width=1))
    dwg.add(dwg.line(start=y_axis_start, end=y_axis_end, stroke='black', stroke_width=1))
    
    # Draw Y-axis ticks and labels
    for i in range(num_units):
        y_val = i * unit
        y_pos = y_offset + bar_height - i * (bar_height / (num_units - 1))
        
        # Add tick mark
        dwg.add(dwg.line(start=(x_offset - 3, y_pos), 
                        end=(x_offset, y_pos), 
                        stroke='black', 
                        stroke_width=1))
        
        # Add label
        dwg.add(dwg.text(str(y_val), 
                        insert=(x_offset - 32, y_pos + 3), 
                        style="font-size:10px; font-family:Courier"))

def draw_coverage_bars(dwg: svgwrite.Drawing, depth_arr: List[List[str]], 
                      x_offset: int, y_offset: int, unit: int) -> None:
    """
    Draw stacked coverage bars for each position.
    
    Args:
        dwg: SVG drawing object
        depth_arr: Array of coverage values for each chain
        x_offset: X-axis starting offset
        y_offset: Y-axis starting offset
        unit: Unit size for scaling bars
    """
    bar_height = 180
    y_base = y_offset + bar_height
    
    for k, coverage in enumerate(depth_arr):
        x = x_offset + k * BOX_SIZE
        
        # Convert to integers and scale to pixels
        pp = float(coverage[0])
        pm = float(coverage[1])
        mp = float(coverage[2])
        mm = float(coverage[3])
        
        # Scale heights (30 pixels per unit)
        pp_height = (pp / unit) * 30 if unit > 0 else 0
        pm_height = (pm / unit) * 30 if unit > 0 else 0
        mp_height = (mp / unit) * 30 if unit > 0 else 0
        mm_height = (mm / unit) * 30 if unit > 0 else 0
        
        # Draw stacked bars (bottom to top: pp, pm, mp, mm)
        current_y = y_base
        
        # mm (top)
        if mm_height > 0:
            dwg.add(dwg.rect(insert=(x + 2, current_y - mm_height - pp_height - pm_height - mp_height),
                            size=(BOX_SIZE - 4, mm_height),
                            fill=CHAIN_COLORS["mm"],
                            stroke="black",
                            stroke_width=0.5))
        
        # mp
        if mp_height > 0:
            dwg.add(dwg.rect(insert=(x + 2, current_y - mp_height - pp_height - pm_height),
                            size=(BOX_SIZE - 4, mp_height),
                            fill=CHAIN_COLORS["mp"],
                            stroke="black",
                            stroke_width=0.5))
        
        # pm
        if pm_height > 0:
            dwg.add(dwg.rect(insert=(x + 2, current_y - pm_height - pp_height),
                            size=(BOX_SIZE - 4, pm_height),
                            fill=CHAIN_COLORS["pm"],
                            stroke="black",
                            stroke_width=0.5))
        
        # pp (bottom)
        if pp_height > 0:
            dwg.add(dwg.rect(insert=(x + 2, current_y - pp_height),
                            size=(BOX_SIZE - 4, pp_height),
                            fill=CHAIN_COLORS["pp"],
                            stroke="black",
                            stroke_width=0.5))

def draw_legend(dwg: svgwrite.Drawing, legend_x: int, legend_y: int) -> None:
    """
    Draw legend explaining bar colors.
    
    Args:
        dwg: SVG drawing object
        legend_x: X-coordinate for legend
        legend_y: Y-coordinate for legend
    """
    legend_items = [
        ('pp', '(+) Strand, forward primer'),
        ('pm', '(+) Strand, reverse primer'),
        ('mp', '(-) Strand, forward primer'),
        ('mm', '(-) Strand, reverse primer')
    ]
    
    for i, (key, label) in enumerate(legend_items):
        y_pos = legend_y + i * (BOX_SIZE + 5)
        dwg.add(dwg.rect(insert=(legend_x, y_pos),
                        size=(BOX_SIZE, BOX_SIZE * 0.8),
                        fill=CHAIN_COLORS[key],
                        stroke="black",
                        stroke_width=0.5))
        dwg.add(dwg.text(label,
                        insert=(legend_x + BOX_SIZE + 5, y_pos + BOX_SIZE * 0.6),
                        style="font-size:10px; font-family:Courier"))

def draw_sequences(dwg: svgwrite.Drawing, seq: str, y_offset: int, 
                  x_offset: int, label: str, legend_x: int) -> None:
    """
    Draw sequence row with colored bases.
    
    Args:
        dwg: SVG drawing object
        seq: Sequence string
        y_offset: Y-coordinate offset
        x_offset: X-coordinate offset
        label: Label for the sequence
        legend_x: X-coordinate for legend
    """
    y = y_offset
    
    for i, base in enumerate(seq.upper()):
        x = x_offset + i * BOX_SIZE
        color = COLORS.get(base, COLORS['N'])
        
        # Draw background
        dwg.add(dwg.rect(insert=(x, y),
                        size=(BOX_SIZE, BOX_SIZE),
                        fill=color,
                        stroke="black",
                        stroke_width=0.5))
        
        # Draw base character
        dwg.add(dwg.text(base,
                        insert=(x + 3, y + BOX_SIZE - 3),
                        fill='black',
                        style="font-size:12px; font-family:Courier"))
    
    # Add label
    dwg.add(dwg.text(label,
                    insert=(legend_x, y + BOX_SIZE * 0.8),
                    style="font-size:12px; font-family:Courier"))

def visualize_sites(offtarget_seq: str, ref_seq: str, strand: str, 
                   depth_file: str, prefix: str) -> None:
    """
    Main visualization function.
    
    Args:
        offtarget_seq: Off-target sequence
        ref_seq: Reference sequence
        strand: Strand direction (+ or -)
        depth_file: Path to depth matrix file
        prefix: Output file prefix
    """
    # Validate inputs
    if not offtarget_seq or not ref_seq:
        print("Error: Empty sequence provided")
        sys.exit(1)
    
    if len(offtarget_seq) != len(ref_seq):
        print(f"Warning: Sequence lengths differ: offtarget={len(offtarget_seq)}, reference={len(ref_seq)}")
    
    # Parse depth file
    depth_arr, depth_sums = parse_depth_file(depth_file, strand)
    
    # Calculate Y-axis scale
    max_val, unit, num_units = calculate_y_axis_scale(depth_sums)
    
    if max_val == 0:
        print("Warning: All coverage values are zero")
    
    # Initialize SVG drawing
    output_file = Path(f'./{prefix}.svg')
    dwg = svgwrite.Drawing(str(output_file), profile='full', size=('100%', '100%'))
    
    # Set margins
    if prefix:
        x_offset = 50
        y_offset = 50
        dwg.add(dwg.text(prefix, 
                        insert=(x_offset, 30), 
                        style="font-size:16px; font-family:Courier; font-weight:bold"))
    else:
        x_offset = 20
        y_offset = 20
    
    # Draw axes and coverage bars
    seq_length = len(offtarget_seq)
    draw_axes(dwg, x_offset, y_offset, seq_length, max_val, unit, num_units)
    draw_coverage_bars(dwg, depth_arr, x_offset, y_offset, unit)
    
    # Draw legend
    legend_x = x_offset + seq_length * BOX_SIZE + 20
    legend_y = y_offset
    draw_legend(dwg, legend_x, legend_y)
    
    # Draw X-axis ticks
    tick_y = y_offset + 180 + 11
    tick_locations = [1, seq_length] + list(range(seq_length + 1)[::10][1:])
    
    for pos in tick_locations:
        if 1 <= pos <= seq_length:
            dwg.add(dwg.text(str(pos),
                            insert=(x_offset + (pos - 1) * BOX_SIZE + 2, tick_y - 2),
                            style="font-size:9px; font-family:Courier"))
    
    # Draw sequences
    seq_y_offset = y_offset + 180 + 11
    draw_sequences(dwg, offtarget_seq, seq_y_offset, x_offset, 
                  'Off-target', legend_x)
    draw_sequences(dwg, ref_seq, seq_y_offset + BOX_SIZE + 5, x_offset,
                  'Reference', legend_x)
    
    # Save SVG file
    try:
        dwg.save()
        print(f"Successfully created: {output_file}")
    except Exception as e:
        print(f"Error saving SVG file: {e}")
        sys.exit(1)

def main():
    """Main entry point."""
    if len(sys.argv) < 5:
        print("Usage: python visualization.py offtarget_seq ref_seq strand count_matrix.txt [id]")
        print("\nArguments:")
        print("  offtarget_seq  - Off-target sequence")
        print("  ref_seq        - Reference sequence")
        print("  strand         - Strand direction (+ or -)")
        print("  count_matrix.txt - Depth matrix file (tab-separated)")
        print("  id             - Optional output file identifier")
        sys.exit(1)
    
    # Parse arguments
    offtarget_seq = sys.argv[1]
    ref_seq = sys.argv[2]
    strand = sys.argv[3]
    depth_file = sys.argv[4]
    prefix = sys.argv[5] if len(sys.argv) > 5 else "output"
    
    # Validate strand
    if strand not in ['+', '-']:
        print(f"Error: Invalid strand '{strand}'. Must be '+' or '-'")
        sys.exit(1)
    
    # Check if depth file exists
    if not Path(depth_file).exists():
        print(f"Error: Depth file '{depth_file}' not found")
        sys.exit(1)
    
    # Generate visualization
    try:
        visualize_sites(offtarget_seq, ref_seq, strand, depth_file, prefix)
    except Exception as e:
        print(f"Error during visualization: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()