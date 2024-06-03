from itertools import combinations
import numpy as np
import sys
from oxDNA_analysis_tools.UTILS.oxview import oxdna_conf, from_path
from oxDNA_analysis_tools.UTILS.RyeReader import describe, get_confs, inbox
from oxDNA_analysis_tools.UTILS.data_structures import TopInfo, TrajInfo
from pathlib import Path
import os
from ipy_oxdna.dna_structure import DNAStructure, DNAStructureStrand, load_dna_structure, DNABase, strand_from_info
from copy import deepcopy
from ipy_oxdna.oxdna_simulation import Simulation 
import copy
from tqdm.auto import tqdm

path = '/home/ava/Dropbox (ASU)/temp/Metabackbone/structure_files/six_helix_oxdna_file/unmodified/1512_bp'
dat_path = os.path.join(path,'1512_bp.dat')
top_path = os.path.join(path,'1512_bp.top')

dna = load_dna_structure(top_path, dat_path) 

 
left_indices = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,391,392,393,394,395,396,397,398,399,400,401,402,403,404,405,406,407,408,409,410,411,412,413,414,415,416,417,418,419,420,421,422,423,424,425,426,427,428,429,430,431,432,433,434,435,436,437,438,1390,1391,1392,1393,1394,1395,1396,1397,1398,1399,1400,1401,1402,1403,1404,1405,1406,1407,1408,1409,1410,1411,1412,1413,1414,1415,1416,1417,1418,1419,1420,1421,1422,1423,1424,1425,1426,1427,1428,1429,1430,1431,1432,1433,1434,1435,1436,1437,1438,1439,1440,1441,1442,1443,1444,1445,1446,1447,1448,1449,1450,1451,1452,1453,1454,2293,2294,2295,2301,2302,2303,2304,2305,2306,2307,2308,2309,2310,2311,2312,2313,2314,2315,2316,2317,2318,2319,2320,2321,2322,2323,2324,2325,2326,2327,2328,2329,2330,2331,2332,2333,2334,2335,2336,2337,2338,2339,2340,2341,2342,2343,2344,2345,2346,2347,2348,2349,2350,2351,2352,2353,2354,2355,2356,2357,2358,2359,2360,2361,2362,2363,2364,2365,2366,2367,2368,2369,2370,2371,2372,2373,2374,2375,2376,2377,2378,2379,2380,2381,2783,2784,2785,2786,2787,2788,2789,2802,2803,2804,2805,2806,3011,3012,3013,3014,3015,3016,3017,3018,3019,3020,3021,3022,3023,3024,3025,3026,3027,3028,3029,3030,3031,3032,3033,3034,3035,3036,3037,3038,3039,3040,3041,3042,3043,3044,3045,3046,3047,3048,3049,3050,3051,3052,3053,3054,3055,3056,3057,3058,3059,3060,3061,3062,3063,3064,3065,3066,3067]
right_indices = [20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,631,632,633,634,635,636,637,638,639,640,641,642,643,644,645,646,647,648,649,650,651,652,653,654,655,656,657,658,659,660,661,662,663,664,665,666,667,668,669,670,671,672,673,674,675,676,677,678,679,680,681,682,683,684,685,686,687,688,689,690,691,692,693,694,695,696,697,698,699,700,701,875,876,877,878,879,880,881,882,883,884,885,886,887,888,889,890,891,892,893,894,895,896,897,898,899,900,901,902,903,904,905,906,907,908,909,910,911,912,913,914,915,916,917,918,919,920,921,922,923,924,925,926,927,928,929,930,931,932,933,934,935,936,937,938,939,940,941,942,943,944,945,946,947,948,949,950,951,952,953,954,955,956,957,958,959,960,961,962,1139,1140,1141,1142,1143,1144,1145,1146,1147,1148,1149,1150,1151,1152,1153,1154,1155,1156,1157,1158,1159,1160,1161,1162,1163,1164,1165,1166,1167,1168,1169,1170,1171,1172,1173,1174,1175,1176,1177,1178,1179,1180,1181,1182,1183,1184,1185,1186,1187,1188,1189,1190,1191,1192,1193,1194,1195,1196,1197,1198,1199,1555,1556,1557,1558,1559,1560,1561,1562,1563,1564,1565,1566,1567,1568,1569,1570,1571,1572,1573,1574,1575,1576,1577,1578,1579,1580,1581,1582,1583,1584,1585,1586,1587,1588,1589,1590,1591,1592,1593,1594,1595,1596,1597,1598,1599,1600,1601,1602,1603,1604,1605,1606,1607,1608,1609,1610,1611,1612,1613,1614,1615,1616,1617,1618,1619,1620,1621,1622,1623,1624,1625,1626,1627,1628,1629,1630,1631,1632,1633,1634,1635,1636,1637,1638,1639,1640,1641,1642,1643,1644,1645,1646,1647,1648,1649,1650,1651,1652,1653,1654,1655,1656,1657,1658,1667,1668,1669,1670,2382,2383,2384,2385,2386,2387,2388,2389,2390,2391,2392,2393,2394,2395,2396,2397,2398,2399,2400,2401,2402,2403,2404,2405,2406,2407,2408,2409,2410,2411,2412,2413,2414,2415,2416,2417,2418,2419,2420,2421,2422,2423,2424,2425,2426,2427,2428,2429,2430,2431,2432,2433,2434,2435,2436,2437,2438,2439,2440,2441,2442,2443,2444,2445,2446,2447,2448,2449,2450,2451,2452,2453,2454,2455,2456,2457,2458,2459,2460,2461,2462,2463,2464,2465,2466,2467,2479,2480,2481,2482,2483,2484,2491,2492,2493,2494,2495,2496,2497,2498,2499,2500,2501,2502,2818,2819,2820,2821,2822,2823]

# Find the longest strand in the DNA structure
def find_longest_strand(dna):
    longest_strand = None
    max_length = 0
    
    for strand in dna.strands:
        if len(strand.bases) > max_length:
            max_length = len(strand.bases)
            longest_strand = strand
    
    return longest_strand
    
longest_strand = list(find_longest_strand(dna)) 
  
# find the cross over in the scaffold   
def find_cross_over_in_longest_strand(dna):
    min_distance = float('inf')        # Initialize the minimum distance to infinity
    max_index_difference = 0           # Initialize the maximum index difference to 0
    cross_over_bases = (None, None)    # Initialize the crossover bases to (None, None)
    num_bases = len(longest_strand)
    
    for i in range(num_bases):
        for j in range(i + 1, num_bases):
            base_i = longest_strand[i] 
            base_j = longest_strand[j]
            index_difference = abs(base_i.uid - base_j.uid)
            distance = np.linalg.norm(np.array(base_i.pos) - np.array(base_j.pos))
            
            if index_difference > max_index_difference or (index_difference == max_index_difference and distance < min_distance):
                max_index_difference = index_difference
                min_distance = distance
                cross_over_bases = (base_i, base_j)
                
    return cross_over_bases, max_index_difference, min_distance
    
    
    
    
# Function to find safe points P without scaffold crossover bases
def find_safe_points(dna, left_indices, right_indices, t_values, sphere_radius):
    
    longest_strand = find_longest_strand(dna)
    cross_over_bases, max_index_difference, min_distance = find_cross_over_in_longest_strand(longest_strand)
    
    left_pos = []
    right_pos = []  
    
    for strand in dna.strands:
        for base in strand: 
            if base.uid in left_indices:
                left_pos.append(base.pos)
            elif base.uid in right_indices:  
                right_pos.append(base.pos)
                
    if left_pos:
        cms_left_side = np.mean(left_pos, axis=0)
    if right_pos:
        cms_right_side = np.mean(right_pos, axis=0)
    
    safe_points = []
    
    for t in t_values:
        P = np.array(cms_left_side + t * (cms_right_side - cms_left_side))
         
        safe = True
        for crossover_base in cross_over_bases:
            distance = np.linalg.norm(np.array(crossover_base.pos) - P)
            if distance < sphere_radius:
                safe = False
                break

        if safe:
            safe_points.append(P)
    
    return safe_points 


def find_bases_in_sphere(dna, P, sphere_radius):
    staples_in_sphere = []
    for strand in dna.strands:
        for base in strand:
            base_position = np.array(base.pos)
            distance = np.linalg.norm(base_position - P)
            if distance < sphere_radius:
                staples_in_sphere.append(base.uid)
    return staples_in_sphere


def dna_structures(dna, staples_in_sphere):
    strand_list = []
    sphere_list = []
    new_structures = []

    for strand in dna.strands:
        bases_in_strand = []
        bases_in_sphere = []

        for base in strand:
            if base.uid in staples_in_sphere:
                base_info = (base.base, base.pos, base.a1, base.a3)
                bases_in_sphere.append(base_info)
            else:
                base_info = (base.base, base.pos, base.a1, base.a3)
                bases_in_strand.append(base_info)

        if bases_in_strand:
            new_strands = strand_from_info(bases_in_strand)
            strand_list.append(new_strands)
            
        if bases_in_sphere:
            new_strands_sphere = strand_from_info(bases_in_sphere)
            sphere_list.append(new_strands_sphere)
            
    while sphere_list:
        rmv_staple = sphere_list.pop()
        new_dna_structure = DNAStructure(strand_list + [rmv_staple], dna.time, dna.box, dna.energy)
        new_structures.append(new_dna_structure)

    return new_structures


#example
new_dna_structures = []
t_values = np.linspace(0, 1, 5)
sphere_radius = 2.7
safe_points = find_safe_points(dna, left_indices, right_indices, t_values, sphere_radius)
# print(safe_points)
for P in safe_points:
    staples_in_sphere = find_bases_in_sphere(dna, P, sphere_radius)
    new_dna_structures.extend(dna_structures(dna, staples_in_sphere))

print("Number of new DNA structures:", len(new_dna_structures))


for i, structure in enumerate(new_dna_structures):
    print(f"Structure {i + 1}: Number of bases = {structure.get_num_bases()}")

print("Original DNA number of bases:", dna.get_num_bases())

path_struct = Path('/home/ava/Dropbox (ASU)/temp/Metabackbone/structure_files/six_helix_oxdna_file/unmodified')
for i, new_dna_structure in enumerate(new_dna_structures):
    dat_path = path_struct / f"1512_bp_mod_{i + 1}.dat"
    top_path = path_struct / f"1512_bp_mod_{i + 1}.top"
    new_dna_structure.export_top_conf(top_path, dat_path)