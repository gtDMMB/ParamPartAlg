from math import floor, ceil
import time

import os
import numpy as np

from geometry import intersect_half_planes, compute_TLs
from ViennaInterface import ViennaInterface
from debug_funcs import nice_print_HP, nice_print_points
from saving import *

ROUND_PRECISION = 5
COMPARE_PRECISION = 2**(-20)

INITIALIZATION_TIME = []

add_new_region_times = []
oracle_times = []
iteration_times = []
region_added = []
oracle_call_counts = []
oracle_points = []

class Point():
   def __init__(self, ac, regs):
      self.ac = ac
      self.regs = regs
      self.full_sigs = np.array([r.full_sig for r in self.regs])
      self.sigs = [r.xz for r in self.regs]

   def rounded(self):
      ac = self.ac
      return set([
            (ceil(ac[0]), ceil(ac[1])), 
            (floor(ac[0]), ceil(ac[1])), 
            (floor(ac[0]), floor(ac[1])), 
            (ceil(ac[0]), floor(ac[1]))
            ])

   def update_regs(self, new_reg):
      if new_reg not in self.regs:
         self.regs.append(new_reg)
         self.full_sigs = np.array([r.full_sig for r in self.regs])
         self.sigs = [r.xz for r in self.regs]
         return True
      
      else:
         return False
      
class Region():
   def __init__(self, full_sig, hp_intersection, base_point, zero_area = False):
      #STATIC
      self.full_sig = full_sig
      self.xz = tuple(int(i) for i in full_sig[::2])

      #DYNAMIC
      self.hps = hp_intersection
      self.bp = base_point
      self.zero_area = zero_area
      if not self.zero_area:
         self.vertices = set(tuple(float(i) for i in np.round(x, ROUND_PRECISION)) for x in self.hps.intersections)
      else:
         if hp_intersection != None:
            self.hps = None
            self.vertices = set([tuple(float(j) for j in i) for i in np.round(hp_intersection, ROUND_PRECISION)])
         else:
            self.vertices = set([tuple(float(i) for i in np.round(self.bp, ROUND_PRECISION))])

   def update_region(self, new_sig, b_val):
      new_TLs = compute_TLs(np.array([self.full_sig]), np.array([new_sig]), b_val=b_val)

      if not self.zero_area:
         if self.full_sig[::2] == new_sig:
            raise Exception(f"New signature {self.new_sig} found for existing region, {self.full_sig}:{self.vertices}")
         old_hps = self.hps.halfspaces[self.hps.dual_vertices]
         all_TLs = np.vstack((new_TLs, old_hps))

         self.hps, self.bp, self.zero_area = intersect_half_planes(self.bp, all_TLs)
         if self.zero_area:
            rounded_bp = tuple(float(x) for x in np.round(self.bp, ROUND_PRECISION))
            
            if self.hps != None: #If region is 1D then self.hps contains the 2 endpoints of the segment
               new_vertices = set([tuple(float(p) for p in np.round(x)) for x in self.hps])
               self.hps = None
               # print("SHRUNK TO LINE", self.xz)
            else:
               new_vertices = set([rounded_bp])
               # print("SHRUNK TO POINT", self.xz)

            lost_vertices = self.vertices - new_vertices
            if len(new_vertices) == 0:
               raise ValueError(f"Intersection of {self.vertices} with {all_TLs[0]} is empty.")
            
            self.vertices = new_vertices
         else:
            new_vertices = set(tuple(float(i) for i in np.round(x, ROUND_PRECISION)) for x in self.hps.intersections)
            lost_vertices = self.vertices - new_vertices

            if len(new_vertices) == 0:
               raise ValueError(f"Intersection of {self.vertices} with {all_TLs[0]} is empty.")

            self.vertices = new_vertices

         return lost_vertices
      else:
         ls = self.lost_vertices_zero_area(new_TLs[0])
         return ls
      
   def lost_vertices_zero_area(self, new_TL):
      if len(self.vertices) == 1:
         if (new_TL[:-1].dot(self.bp) + new_TL[-1]) > COMPARE_PRECISION:
            raise Exception(f"New halfplane {new_TL}, doesn't contain 1D region {self.bp, self.full_sig}")
         return set()
      else:
         lost_vertices = set()
         for v in self.vertices:
            if (new_TL[:-1].dot(v) + new_TL[-1]) > COMPARE_PRECISION:
               lost_vertices.add(v)
         
         if len(lost_vertices) == 0:
            return lost_vertices
         elif len(lost_vertices) == 2:
            raise ValueError(f"Intersection of {self.vertices} with {new_TL} is empty.")
         else:
            v1, v2 = self.vertices
            dx0, dz0, b0 = new_TL
            da = v1[0] - v2[0]
            dc = v1[1] - v2[1]
            dx1, dz1, b1 = dc, -da, v1[0]*dc - v1[1]*da

            A = np.array([[dx0, dz0],[dx1, dz1]])
            b = np.array([-b0,b1])
            try:
               intersection_point = np.linalg.solve(A,b)
            except np.linalg.LinAlgError as e:
               print(self.full_sig, self.vertices, new_TL, A, b)
               raise(e)

            new_vertex = set([tuple(float(x) for x in np.round(intersection_point, ROUND_PRECISION))])
            self.vertices -= lost_vertices
            self.vertices.update(new_vertex)
            self.bp = np.array(list(new_vertex)[0])

            return lost_vertices

def new_0_area_region(hp_intersection, sig_reg, sig_map, TLS, intersection_sigs):
   if hp_intersection == None:
      new_point = tuple(float(i) for i in np.round(sig_reg.bp, ROUND_PRECISION))

      #Find TLs that intersect new_point. 
      reg_idx = (np.abs(TLS[:,:-1].dot(new_point) + TLS[:,-1]) < COMPARE_PRECISION)[:-4]
      reg_sigs = [sig_map[tuple(int(x) for x in s[::2])] for s in intersection_sigs[reg_idx]]

      return [Point(new_point, regs=[sig_reg] + reg_sigs)]
   else:
      new_points = []
      for v in hp_intersection:
         #Find TLs that intersect v. 
         reg_idx = (np.abs(TLS[:,:-1].dot(v) + TLS[:,-1]) < COMPARE_PRECISION)[:-4]
         reg_sigs = [sig_map[tuple(int(x) for x in s[::2])] for s in intersection_sigs[reg_idx]]

         point = tuple(float(i) for i in np.round(v, ROUND_PRECISION))
         new_points.append(Point(point, regs=[sig_reg] + reg_sigs))
      
      return new_points


#Given a point, signature, and regions update all the regions around the point and add the the region R(sig) to REGION
def add_new_region(REGIONS, sig_map, all_sigs, frame_HPs, point, full_sig, b_val):
   if full_sig[::2] in sig_map.keys():
      print("PRINTING INTERSECTIONS")
      for r in sig_map.values():
         if r.zero_area:
            print(*r.vertices, sep=",")
         else:
            nice_print_points(r.hps.intersections) 
      raise Exception(f"Signature found twice, {(full_sig, point.ac, point.sigs, point.full_sigs, (sig_map[full_sig[::2]].vertices))}")


   intersection_sigs = np.array(all_sigs)
   #Update all regions
   points_to_remove = set()
   for r in REGIONS:
      region = sig_map[r]
      points_to_remove.update(region.update_region(full_sig, b_val))

   all_sigs.append(full_sig)

   TLS = np.vstack((compute_TLs(np.array(full_sig), intersection_sigs, b_val=b_val), frame_HPs))
   hp_intersection, bp, zero_area = intersect_half_planes(np.array(point.ac), TLS)

   sig_reg = Region(full_sig, hp_intersection, bp, zero_area)
   REGIONS.append(sig_reg.xz)
   sig_map[sig_reg.xz] = sig_reg

   if zero_area:
      return new_0_area_region(hp_intersection, sig_reg, sig_map, TLS, intersection_sigs), points_to_remove
      
   new_points = []
   if len(sig_reg.hps.dual_facets) != sig_reg.hps.intersections.shape[0]:
      raise Exception("Number of intersection points and number of dual facets does not match.")
   for i, sigidxs in enumerate(sig_reg.hps.dual_facets):
      point = tuple(float(i) for i in np.round(sig_reg.hps.intersections[i], ROUND_PRECISION))

      true_idxs = [i for i in sigidxs if i < len(intersection_sigs)]
      adjacent_sigs = list((int(s[0]), int(s[1]), int(s[2]), float(np.round(s[3], ROUND_PRECISION))) for s in intersection_sigs[true_idxs])
      adjacent_sigs += [full_sig]

      adjacent_regs = [sig_map[s[::2]] for s in adjacent_sigs]
      new_points.append(Point(point, adjacent_regs))

   return new_points, points_to_remove

"""
params: Points p1, p2 with scores s1, s2
return: All tipping points beteween p1 and p2
restrictions: p1 and p2 must be collinear (integer points) on a horizontal or vertical line
"""
def find_collinear_TPs(p1, p2, s1, s2, oracle, b_val):
   intersection = find_tipping_line_intersection(p1, p2, s1, s2, b_val)
   
   if (intersection == None):
      return [],[]
   
   #DOUBLE CHECK THIS CONDITION
   if (intersection == p1):
      return [],[]
   
   elif (intersection == p2):
      return [],[]
   
   #Find the ceiling and score point (since one of the cordinates is equal)
   pc, pf = p2, p1
   sc, sf = s2, s1
   if (p1[0] > p2[0] or  p1[1] > p2[1]):
      pc, pf = p1, p2
      sc, sf = s1, s2
   
   # print(intersection)
   intersection_ceil = tuple(ceil(x) for x in intersection)
   intersection_floor = tuple(floor(x) for x in intersection)
   ceil_score = oracle(intersection_ceil[0], b_val, intersection_ceil[1])
   floor_score = oracle(intersection_floor[0], b_val, intersection_floor[1])

   colinear_TPs = []
   colinear_scores = []
   if ((sc[::2] == ceil_score[::2]) or (sf[::2] == floor_score[::2])): #ONLY COMPARE xz sig not entire sig :)
      colinear_TPs.append(intersection)
      colinear_scores.append(sc)

   A = find_collinear_TPs(intersection_floor, pf, floor_score, sf, oracle, b_val)
   B = find_collinear_TPs(intersection_ceil, pc, ceil_score, sc, oracle, b_val)

   colinear_TPs.extend(A[0])
   colinear_TPs.extend(B[0])
   colinear_scores.extend(A[1])
   colinear_scores.extend(B[1])

   return colinear_TPs, colinear_scores

"""
Requirements p1, p2 are horizonatl or vertical
"""
def find_tipping_line_intersection(p1, p2, s1, s2, b_val):
   #No intersection between points with the same signature
   if s1[::2] == s2[::2]:
      return None
      
   TLs = compute_TLs(s1, np.array([s2]), b_val=b_val)
   delta_a, delta_c, const = TLs[0]

   if p1[0] == p2[0]:
      a0 = p1[0]
      return ((a0, (delta_a * a0 + const) / (-delta_c)) if (delta_c != 0) else None)
   elif p1[1] == p2[1]:
      c0 = p1[1]
      return (((delta_c * c0 + const) / (-delta_a), c0)  if (delta_a != 0) else None)
   else:
      raise Exception("Tipping line intersection for non-horizontal/vertical p1, p2 is not implemented.")

"""
A complex initization procedure which is not used.
"""
def run_WBP_HP(oracle, a=-10000, A=10000, c=-10000, C=10000, b_val=0):
   INIT_START = time.process_time_ns()
   #We need the internal frame points, boundry points, and regions induced by frame points
   corners = [(A, C), (a, C), (a, c), (A, c)]
   corner_scores = [oracle(corner[0], b_val, corner[1]) for corner in corners]
   
   boundary_tps = []
   all_boundary_scores = []
   for i in range(4):
      # print(f"FRAME{i}")
      p1 = corners[i]
      p2 = corners[(i+1)%4]
      colinear_tps, colinear_scores = find_collinear_TPs(p1, p2, corner_scores[i], corner_scores[(i+1)%4],oracle, b_val)
      boundary_tps.extend(colinear_tps)
      all_boundary_scores.extend(colinear_scores)
   
   #Create list of scores without xz duplicates
   boundary_sigs = set()
   boundary_scores = []
   for s in all_boundary_scores:
      if s[::2] not in boundary_sigs:
         boundary_sigs.add(s[::2])
         boundary_scores.append(s)
   
   regions = []
   internal_points = {} #Map of (a, c): (point at ac)
   boundary_points = {}
   boundary_half_planes = np.array([[0,1,-C],[-1,0,a],[0,-1,c],[1,0,-A]])
   for i in range(len(boundary_scores)):
      other_scores = boundary_scores[:i] + (boundary_scores[i+1:] if i < len(boundary_tps) - 1 else [])
      other_scores = np.array(other_scores)
      TLs = compute_TLs(np.array(boundary_scores[i]), other_scores, b_val)
      half_planes = np.vstack((TLs, boundary_half_planes))
      
      try:
         hp_intersection, base_point = intersect_half_planes(np.array([0,0]), half_planes)
      except Exception as e:
         """""DEBUG"""
         print(boundary_scores, boundary_scores[:i], boundary_scores[i], other_scores)
         nice_print_HP(half_planes)
         """"END DEBUG"""
         raise (e) 

      region = Region(boundary_scores[i], hp_intersection, base_point)
      regions.append(region)

      # vertices = set(hp_intersection.intersections)
      reg_points = region.hps.intersections
      internal_condition = (reg_points[:,0] != A) & (reg_points[:,0] != a) & (reg_points[:,1] != C) & (reg_points[:,1] != c)
      internal_reg_points = reg_points[internal_condition]
      boundary_reg_points = reg_points[~internal_condition]
      
      for p in internal_reg_points:
         p = tuple(float(x) for x in np.round(p, ROUND_PRECISION))
         try:
            internal_points[p].update_regs(region)
         except KeyError:
            internal_points[p] = Point(p, [region])

      for p in boundary_reg_points:
         p = tuple(float(x) for x in np.round(p, ROUND_PRECISION))
         try:
            boundary_points[p].update_regs(region)
         except KeyError:
            boundary_points[p] = Point(p, [region])

   INIT_END = time.process_time_ns()
   global INITIALIZATION_TIME
   INITIALIZATION_TIME = INIT_END - INIT_START

   return WBP_HP_main(oracle, list(internal_points.values()), list(boundary_points.values()), regions, boundary_half_planes, b_val)


def simple_WBP_init(oracle, a=-10000, A=10000, c=-10000, C=10000, b_val=0):
   INIT_START = time.process_time_ns()
   all_corners = [(A, C), (a, c)] + [(a, C), (A, c)]
   corners = [(A, C), (a, c)]

   boundary_half_planes = np.array([[0,1,-C],[-1,0,a],[0,-1,c],[1,0,-A]])

   corner_scores = [oracle(corner[0], b_val, corner[1]) for corner in corners]
   if corner_scores[0] == corner_scores[-1]:
      corners = [(a, C), (A, c)]
      corner_scores = [oracle(corner[0], b_val, corner[1]) for corner in corners]
      if corner_scores[0] == corner_scores[-1]:
         #Bounding box is only contains 1 score
         hp_intersection, base_point, zero_area = intersect_half_planes(np.array(corners).mean(axis=1), boundary_half_planes)

         if zero_area:
            raise Exception(f"Frame halfplanes returned zero area")
         
         reg = Region(corner_scores[0], hp_intersection, base_point, zero_area)
         INIT_END = time.process_time_ns()
         INITIALIZATION_TIME.append(INIT_END - INIT_START)
         return WBP_HP_main(oracle, [], [Point(p, [reg]) for p in corners], [reg], boundary_half_planes, b_val)
   

   regions = []
   for i in range(2):
      TLs = compute_TLs(np.array(corner_scores[i]), np.array(corner_scores[(i+1)%2::2]), b_val)
      half_planes = np.vstack((TLs, boundary_half_planes))
      hp_intersection, base_point, zero_area = intersect_half_planes(np.array(corners[i])/2, half_planes)
      
      regions.append(Region(corner_scores[i], hp_intersection, base_point, zero_area = zero_area))

   reg = regions[0]
   points = [Point(p, regions) for p in reg.vertices if (p not in  all_corners)]
   confirmed = [Point(p, [r]) for p, r in zip(corners, regions)]

   INIT_END = time.process_time_ns()
   INITIALIZATION_TIME.append(INIT_END - INIT_START)

   return WBP_HP_main(oracle, points, confirmed, regions, boundary_half_planes, b_val)
   # reg

# Main algorithm
def WBP_HP_main(oracle, points, completed_points, regions, frame_HPS, b_val):
   point_oracle = lambda p: oracle(p[0], b_val, p[1])

   point_map = {}
   for p in points + completed_points:
      point_map[p.ac] = p
   
   sig_map = {}
   all_sigs = []
   for r in regions:
      sig_map[r.xz] = r
      all_sigs.append(r.full_sig)

   POINTS = set(p.ac for p in points)
   COMPLETED = set(p.ac for p in completed_points)
   REGIONS = [s for s in sig_map.keys()]

   # print(POINTS, COMPLETED, REGIONS)
   while len(POINTS) != 0:
      iteration_start = time.process_time_ns()

      test_point = point_map[POINTS.pop()]
      
      #Find signatures at all surrounding lattice points
      new_points = []
      oracle_calls = 0
      for p_lattice in test_point.rounded():
         added_region = 0
         oracle_start = time.process_time_ns()
         full_sig = point_oracle(p_lattice)
         oracle_calls += 1
         oracle_points.append((test_point.ac,p_lattice))
         oracle_end = time.process_time_ns()
         oracle_times.append((oracle_end - oracle_start))

         sig = full_sig[::2]

         if sig not in sig_map.keys():
            anr_start = time.process_time_ns()
            added_region = 1
            new_points, points_to_remove = add_new_region(REGIONS, sig_map, all_sigs, frame_HPS, test_point, full_sig,  b_val)
            anr_end = time.process_time_ns()
            
            #Remove points that are contined in the new region. 
            POINTS = POINTS - points_to_remove
            add_new_region_times.append((anr_end - anr_start))
            break
      else:
         COMPLETED.add(test_point.ac)

      for p in new_points:
         if (p.ac == test_point.ac): #If the point is the test point, add it back to points
            POINTS.add(p.ac)
         if (p.ac in POINTS) or (p.ac in COMPLETED): #If we have already found the point, update the regions at the point
            for r in p.sigs:
               point_map[p.ac].update_regs(sig_map[r])
         else: #Otherwise add the new point
            POINTS.add(p.ac)
            point_map[p.ac] = p

      iteration_end = time.process_time_ns()
      iteration_times.append((iteration_end - iteration_start))
      region_added.append(added_region)
      oracle_call_counts.append(oracle_calls)

      
   return COMPLETED, REGIONS, sig_map


if __name__ == "__main__":
   import argparse
   parser = argparse.ArgumentParser(
            prog='TL_HPI',
            description='Computes a,c slice for an RNA polytope')

   parser.add_argument("fasta_path", metavar="FASTA_PATH", help="Path to fasta file.")
   parser.add_argument("save_path", metavar="SAVE_PATH", help="Directory to save output files.")
   parser.add_argument("-b", type = int, default = 0, help="b value for computation. Default 0.")
   parser.add_argument("--geometry", action="store_true", help="Save geometry of regions as SAVE_PATH/seqname_geometry.txt.")
   parser.add_argument("--timing", action="store_true", help="Save timing data as SAVE_PATH/seqname_timing.txt")
   parser.add_argument("--unskewed", action="store_true", help="Compute in unskewed space (i.e. excess branching)")
   parser.add_argument("--bounds", metavar="n", nargs=4, default=[0, 10000, -10000, 10000], type=int, help="Integer bounds in dckals/mol for slice in format a_min a_max c_min c_max. Default 0 10000 -10000 10000.")
   parser.add_argument("--LP", action="store_true", help="Compute without --noLP option in vienna")
   parser.add_argument("-d", metavar="n", type=int, default=2, help="Set dandle mode (1 or 2).")

   args = parser.parse_args()

   path = args.fasta_path
   name = os.path.splitext(os.path.split(args.fasta_path)[-1])[-2]
   b_val = args.b
   interface = ViennaInterface(path, transform=(not args.unskewed), dangles=args.d, noLP = (not args.LP))
   vienna_oracle = lambda a, b, c : interface.vertex_oracle(a, c, b, 1)

   a, A, c, C = args.bounds

   COMPLETED, REGIONS, sig_map = simple_WBP_init(vienna_oracle, a=a, A=A, c=c, C=C, b_val=b_val)
   write_sig_struct_file(args.save_path, name, interface)

   if args.geometry:
      write_geometry_file(args.save_path, name, sig_map)

   if args.timing:
      write_timing_file(args.save_path, name, INITIALIZATION_TIME, add_new_region_times, oracle_times, iteration_times, region_added, oracle_call_counts, oracle_points)