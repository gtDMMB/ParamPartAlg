import os
import RNA
from scipy.spatial import ConvexHull
from geometry import intersect_half_planes
import numpy as np

def find_nucleotide_pairs(dot_brac, sequence=None, save_nc=False):
    con_pairs = ["AU", "CG", "GU", "UA", "GC", "UG"]
    bracket_stack = []
    nucleotide_pairs = set()
    for i, n in enumerate(dot_brac):
        if n == "(":
            bracket_stack.append(i) 
        if n == ")":
            start = bracket_stack.pop()
            if save_nc or f"{sequence[start]}{sequence[i]}" in con_pairs: #Save connonical paris
                nucleotide_pairs.add((start+1, i+1)) #Reindex
    return nucleotide_pairs

def write_sig_struct_file(save_path, name, interface):
   with open(os.path.join(save_path, f"{name}_sig_structs.txt"), "w") as f:
      sig_structs = set(interface.computed.values())
      f.write(f"> NAME: {name}\n")
      f.write(f"> NUM STRUCTS: {len(sig_structs)}\n")
      f.write(f"> SEQ: {interface.seq}\n")
      f.write("\n")
      for sig, struct in sig_structs:
         f.write(f"{sig}: {struct}\n")



def write_geometry_file(save_path, name, sig_map):
   with open(os.path.join(save_path, f"{name}_geometry.txt"), "w") as f: 
      for r in sig_map.values():
         full_sig = r.full_sig
         points = r.vertices
         
         f.write(f"{full_sig}:{';'.join(str(p) for p in points)}\n")

def create_prediction_file(save_path, name, sig_map, interface, noLP, dangles):
   import pandas as pd
   sequence = interface.seq
   sig_structs = {sig: struct for sig, struct in set(interface.computed.values())}

   area_map = {}
   #Ax + b <= 0
   bounding_box = np.array([[1, 0, -800],[-1, 0, 300],[0, 1, -100],[0, -1, -350]])
   bounding_points = np.array([[800,100],[300,100],[300,-350],[800,-350]])


   for sig, reg in sig_map.items():
      array_points = np.array([np.array(p) for p in reg.vertices])
      if reg.zero_area:
         area_map[sig] = 0
      else:
         n = array_points.shape[0]
         m = 4
         n1 = reg.hps.halfspaces.shape[0]
         contains = np.any(np.all(bounding_box[:,:-1].dot(array_points.T) + np.vstack([bounding_box[:,-1]]*n).T <= np.zeros((m,n)), axis=0))
         is_contained = np.any(np.all(reg.hps.halfspaces[:,:-1].dot(bounding_points.T) + np.vstack([reg.hps.halfspaces[:,-1]]*m).T <= np.zeros((n1,m)), axis=0))

         if contains or is_contained:
            hpi, _, zero_area = intersect_half_planes(np.array([0,0]), np.vstack((reg.hps.halfspaces, bounding_box)))
            if zero_area:
               area_map[sig] = 0
            else:
               area_map[sig] = ConvexHull(hpi.intersections).volume
         else:
            area_map[sig] = 0


   ssa_list = [] 
   for sig, struct in sig_structs.items():
      ssa_list.append([sig, struct, area_map[sig[::2]]])
   
   RNA.params_load_RNA_Turner2004()
   md = RNA.md()
   md.noLP = noLP
   md.dangles = dangles
   fc = RNA.fold_compound(sequence, md)
   MFE_struct, _ = fc.mfe()
   fc.pf()
   CMP_struct, _ = fc.centroid()

   ssa = pd.DataFrame(ssa_list, columns=["sig","struct","area"])
   ssa["normed_area"] = (ssa["area"] / ssa["area"].sum())

   MAP_struct = ssa[ssa["area"] == ssa["area"].max()].iloc[0]["struct"]
   non_zero = ssa[ssa["area"] > 0]

   pair_probabilities = {}
   centroid_pairs = set()
   for row in non_zero.itertuples():
      if row.area > 0:
         pairs = find_nucleotide_pairs(row.struct, sequence, save_nc=True)
         for p in pairs:
            try:
               pair_probabilities[p] += row.normed_area
            except KeyError:
               pair_probabilities[p] = row.normed_area

            if pair_probabilities[p] > 0.5:
               centroid_pairs.add(p)

   centroid_area_struct = ["."]*len(sequence)
   for pi, pj in centroid_pairs:
      centroid_area_struct[pi - 1] = "("
      centroid_area_struct[pj - 1] = ")"

   CAP_struct = "".join(centroid_area_struct)

   print(non_zero.head())

   with open(os.path.join(save_path, f"{name}_predictions.txt"), "w") as f:
      f.write(f"> NAME: {name} \n")
      f.write(f"> SEQ: {sequence} \n")
      f.write(f"> MFE: {MFE_struct} \n")
      f.write(f"> MAP: {MAP_struct} \n")
      f.write(f"> CMP: {CMP_struct} \n")
      f.write(f"> CAP: {CAP_struct} \n")
      f.write(f"> MFE = default vienna MFE, NAT = Native Structure, MAP = Maximum Area Prediction, CMP = Centroid MFE Prediction, CAP = Centroid Area Prediction")



def write_timing_file(save_path, name, INITIALIZATION_TIME, add_new_region_times, oracle_times, iteration_times, region_added, oracle_call_counts, oracle_points):
      with open(os.path.join(save_path, f"{name}_timing.txt"), "w") as f:
         f.write(f"INIT_TIME: {INITIALIZATION_TIME[0]}\n")
         f.write(f"ADD_NEW_REGION: {','.join(str(t) for t in add_new_region_times)} \n")
         f.write(f"oracle_times: {','.join(str(t) for t in oracle_times)}\n")
         f.write(f"iteration_times: {','.join(str(t) for t in iteration_times)}\n")         
         f.write(f"region_added: {','.join(str(t) for t in region_added)}\n")
         f.write(f"oracle_call_counts: {','.join(str(t) for t in oracle_call_counts)}\n")
         f.write(f"oracle_call_gridpoints: {';'.join(str(t[1]) for t in oracle_points)}\n")
         f.write(f"oracle_call_points: {';'.join(str(t[0]) for t in oracle_points)}")
