import os
import numpy as np

# Load DNA structure
path = '/home/ava/Dropbox (ASU)/temp/Metabackbone/structure_files/six_helix_oxdna_file/unmodified/1512_bp'
dat_path = os.path.join(path,'1512_bp.dat')
top_path = os.path.join(path,'1512_bp.top')
dna = load_dna_structure(top_path, dat_path) 

# Indices for left and right sides
left_indices = [...]
right_indices = [...]

# Function to find point on the axis
def find_axis_point(dna, left_indices, right_indices, t):
    left_pos = [base.pos for strand in dna.strands for base in strand if base.uid in left_indices]
    right_pos = [base.pos for strand in dna.strands for base in strand if base.uid in right_indices]

    if left_pos:
        cms_left_side = np.mean(left_pos, axis=0)
    if right_pos:
        cms_right_side = np.mean(right_pos, axis=0)
    
    P = cms_left_side + t * (cms_right_side - cms_left_side)   
    return P

# Function to find bases in a sphere
def find_bases_in_sphere(dna, P, sphere_radius):
    bases_in_sphere = []
    for strand in dna.strands:
        for base in strand:
            base_position = np.array(base.pos)
            distance = np.linalg.norm(base_position - P)
            if distance < sphere_radius:
                bases_in_sphere.append(base.uid)
    return bases_in_sphere

# Function to create new DNA structures by removing staples one by one
def create_new_structures(dna, staples_in_sphere):
    new_structures = []
    for staple_uid in staples_in_sphere:
        new_strands = []
        for strand in dna.strands:
            new_bases = [base for base in strand if base.uid != staple_uid]
            if new_bases:
                new_strands.append(strand_from_info(new_bases))
        
        new_dna_structure = DNAStructure(new_strands, dna.time, dna.box, dna.energy)
        new_structures.append(new_dna_structure)
    return new_structures

# Main process
P = find_axis_point(dna, left_indices, right_indices, 0.5)
sphere_radius = 2.00
staples_in_sphere = find_bases_in_sphere(dna, P, sphere_radius)
new_dna_structures = create_new_structures(dna, staples_in_sphere)

print("Number of new DNA structures:", len(new_dna_structures))



# #example
# new_dna_structures = []
# t_values = np.linspace(0, 1, 5)
# sphere_radius = 2.7
# safe_points = find_safe_points(dna, left_indices, right_indices, t_values, sphere_radius)
# # print(safe_points)
# for P in safe_points:
#     staples_in_sphere = find_bases_in_sphere(dna, P, sphere_radius)
#     new_dna_structures.extend(dna_structures(dna, staples_in_sphere))

# print("Number of new DNA structures:", len(new_dna_structures))


# for i, structure in enumerate(new_dna_structures):
#     print(f"Structure {i + 1}: Number of bases = {structure.get_num_bases()}")

# print("Original DNA number of bases:", dna.get_num_bases())

# path_struct = Path('/home/ava/Dropbox (ASU)/temp/Metabackbone/structure_files/six_helix_oxdna_file/unmodified')
# for i, new_dna_structure in enumerate(new_dna_structures):
#     dat_path = path_struct / f"1512_bp_mod_{i + 1}.dat"
#     top_path = path_struct / f"1512_bp_mod_{i + 1}.top"
#     new_dna_structure.export_top_conf(top_path, dat_path)