# Test composition parsing and extraction logic
from tc_batch_cli import _parse_composition
import re

# Test composition parsing
comp1 = '0.465LiF-0.115NaF-0.42KF'
comp2 = '0.42KF-0.465LiF-0.115NaF'

_, _, label1 = _parse_composition(comp1)
_, _, label2 = _parse_composition(comp2)

print(f'Composition 1: {comp1}')
print(f'Parsed label 1: {label1}')
print(f'Composition 2: {comp2}')
print(f'Parsed label 2: {label2}')
print(f'Labels match: {label1 == label2}')

# Test composition extraction from source
src = '0.465LiF-0.115NaF-0.42KF (Merritt, 2022)'

if '(' in src and ')' in src:
    exp_composition_str = src.split('(')[0].strip()
    print(f'\nSource: {src}')
    print(f'Extracted composition: {exp_composition_str}')

    _, _, exp_comp_label = _parse_composition(exp_composition_str)
    print(f'Parsed label: {exp_comp_label}')
    print(f'Matches comp1 label: {exp_comp_label == label1}')
