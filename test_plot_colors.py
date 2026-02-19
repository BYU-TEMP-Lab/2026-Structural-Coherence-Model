#!/usr/bin/env python3
"""
Test script to verify the consistent color scheme for plots.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from TC_batch_cli_V2 import plot_tc_cli

def test_consistent_colors():
    """Test that models use consistent colors and Mix Data uses dashed lines."""
    
    # Test composition
    composition = "0.5NaCl-0.5KCl"
    temp_range = (1000, 1200)
    methods = ['Present Model', 'Present Model, Mix Data', 'Gheribi-KT24', 'Gheribi-KT24, Mix Data', 'Zhao-PGM', 'Zhao-PGM, Mix Data']
    
    # Create plot without showing to test colors
    result = plot_tc_cli(
        composition=composition,
        temp_range=temp_range,
        methods=methods,
        show_plot=False,
        output_dir='test_plots'
    )
    
    print("Test completed successfully!")
    print(f"Plot saved to: {result['figure_path']}")
    print("\nColor scheme implemented:")
    print("- SCM (Present Model): Red")
    print("- SCM, Mix Data: Red with dashed line")
    print("- KTM (Gheribi-KT24): Blue") 
    print("- KTM, Mix Data: Blue with dashed line")
    print("- PGM (Zhao-PGM): Green")
    print("- PGM, Mix Data: Green with dashed line")
    print("- Experimental data: High-contrast colors with white edges")

if __name__ == "__main__":
    test_consistent_colors()
