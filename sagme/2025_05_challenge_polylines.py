# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 08:04:59 2024

@author: jpn3t
"""
# https://oeis.org/A066372 diff shapes always bending
# https://oeis.org/A001444 diff shapes but diff reflections count
# https://oeis.org/A001997 diff shapes bending n-1 times
#   https://oeis.org/A001997/a001997.pdf

# but none are exactly what we want (!)

from itertools import product

reversal   = lambda p: p[::-1]
reflection = lambda p: p.replace('R','*').replace('L','R').replace('*','L')

def canonize(polyline):
  """ paths are invariant across reversals, reflections, or both """
  return min(polyline, 
             reversal(polyline),
             reflection(polyline),
             reflection(reversal(polyline)))


def path(polyline):
  """ return the points occupied by a given polyline """
  d = 1j # all polylines start at 0 going north
  res = [0, cur := 0+d]
  for p in polyline:
    if p == 'R': d *= -1j
    if p == 'L': d *=  1j
    res.append(cur := cur+d)
  return res


# def canonize_path(polyline):
#   """ shift bottom-left point to origin """
#   p_path = path(polyline)
#   corner = min(p_path, key=lambda z: (z.imag, z.real))
#   where = p_path.index(corner)
#   shift = [x-corner for x in p_path]
#   return tuple(shift[where:] + shift[:where])

# def canonize_path(polyline):
#   """ shift bottom-left most point to origin """
#   p_path = path(polyline)
#   corner = min(p_path, key=lambda z: (z.imag, z.real))
#   shift = [x-corner for x in p_path]
#   return frozenset({(x,y) for a,b in zip(shift, shift[1:])
#                           for x,y in ((a,b),(b,a))})


def canonize_path_aux(p_path):
  """ shift bottom-left most point to origin """
  corner = min(p_path, key=lambda z: (z.imag, z.real))
  shift = [x-corner for x in p_path]
  return tuple(sorted({(x.real, x.imag, y.real, y.imag) 
                        for a,b in zip(shift, shift[1:])
                        for x,y in ((a,b),(b,a))}))


def canonize_path(polyline): # need to consider rotations
  p_path = path(polyline)
  path0 = canonize_path_aux(p_path)
  path1 = canonize_path_aux([ p.imag + p.real*1j for p in p_path])
  path2 = canonize_path_aux([-p.imag + p.real*1j for p in p_path])
  path3 = canonize_path_aux([ p.imag - p.real*1j for p in p_path])
  # return min(path0, path1, path2, path3)
  path4 = canonize_path_aux([-p.imag - p.real*1j for p in p_path])
  path5 = canonize_path_aux([-p.real + p.imag*1j for p in p_path])
  path6 = canonize_path_aux([-p.real - p.imag*1j for p in p_path])
  path7 = canonize_path_aux([+p.real - p.imag*1j for p in p_path])
  return min(path0, path1, path2, path3, path4, path5, path6, path7)


def is_non_coincidental(polyline):
  """ check if any pair of lines don't overlap """
  d = 1j # all polylines start at 0 going north
  seen = { (0, cur:=0+d) } 
  
  for p in polyline:
    if p == 'R': d *= -1j
    if p == 'L': d *=  1j
    if (cur, cur+d) in seen or (cur+d, cur) in seen:
      return False
    seen.add( (cur, cur+d) )
    cur += d
  return True


def polylines(m):
  """ set of distinct polylines """
  res, seen_paths = [], set()
  for polyline in product('0RL', repeat=m-1):
    p = canonize(''.join(polyline))
    if is_non_coincidental(p):
      path = canonize_path(p)
      if path not in seen_paths:
        seen_paths.add(path)
        res.append(p)
  return res

############

import matplotlib.pyplot as plt

def plot_poly(p, i, sz=20, scale=3, per_line=6):
  p_path = path(p)
  # path = shift_path(path)
  for a,b in zip(p_path, p_path[1:]):
    dx, dy = sz*(i%per_line), -sz*(i//per_line)
    x = (scale*a.real + dx, scale*b.real + dx)
    y = (scale*a.imag + dy, scale*b.imag + dy)
    plt.plot(x, y, 'k-', linewidth=.5)


def draw_polylines(ps):
  plt.figure(figsize=(10,5))
  plt.axes().set_aspect('equal')
  plt.axis('off')
  for i,p in enumerate(ps):
    plot_poly(p, i)
  plt.savefig("polylines.svg")
  plt.show()

# ps = polylines(m=5)
# draw_polylines(ps)

####################################################
# def canonize_path_aux(polyline):
#   """ shift bottom-left point to origin """
#   p_path = path(polyline)
#   corner = min(p_path, key=lambda z: (z.imag, z.real))
#   where = p_path.index(corner)
#   shift = [x-corner for x in p_path]
#   return tuple(shift[where:] + shift[:where])

def all_paths(polyline):
  p_path = path(polyline)
  path0 = tuple( p_path )
  path1 = tuple( p.imag + p.real*1j for p in p_path)
  path2 = tuple(-p.imag + p.real*1j for p in p_path)
  path3 = tuple( p.imag - p.real*1j for p in p_path)
  path4 = tuple(-p.imag - p.real*1j for p in p_path)
  path5 = tuple(-p.real + p.imag*1j for p in p_path)
  path6 = tuple(-p.real - p.imag*1j for p in p_path)
  path7 = tuple(+p.real - p.imag*1j for p in p_path)
  return {path0, path1, path2, path3, path4, path5, path6, path7}

def shift_path(path, dz):
  return tuple(z+dz for z in path)

def overlap(path1, path2):
  set1 = {(x,y) for a,b in zip(path1, path1[1:]) for x,y in ((a,b), (b,a))}
  set2 = {(x,y) for a,b in zip(path2, path2[1:]) for x,y in ((a,b), (b,a))}
  return set1 & set2

def distance(z, center=22):
  z -= center
  return (z.real**2 + z.imag**2)**.5


def plot_piece(p):
  for a,b in zip(p, p[1:]):
    x = (a.real, b.real)
    y = (a.imag, b.imag)
    plt.plot(x, y, 'k-', linewidth=.5)
  plt.plot(x[1], y[1], 'o')
    
def draw_sol(sol):
  plt.figure(figsize=(10,6))
  plt.axes().set_aspect('equal')
  plt.axis('off')
  for p in sol:
    plot_piece(p)
  plt.savefig("polylines.svg")
  plt.show()
  
  
  

def backtrack(ps, cur=0, sol=None, seen=None, best=None):
  if best is None:
    best = []
    sol = []
    seen = set()
  if not ps: # placed all polylines!
    if any(z.imag == 0 for z in sol[-1]):
      draw_sol(sol) # (need to compute score)
    return
  
  for p in ps:
    for path in all_paths(p):
      new_path = shift_path(path, cur)
      if any(z.real < -1 or z.imag < 0 for z in new_path): 
        continue
      lines = {(x,y) for a,b in zip(new_path, new_path[1:]) 
                     for x,y in ((a,b), (b,a))}
      if lines & seen: # no overlap with prev moves needed
        # TODO: however we must deal with polylines that are squares...
        continue
      if any(z in seen for z in new_path[1:]): # polylines should not touch
        continue
      if len(sol)>1 and new_path[-1].real + 1 < sol[-2][-1].real:
        continue
      # if 2 <= distance(new_path[-1], 3) <= 5.5: # close enough to circle
      if 6 <= distance(new_path[-1], 8) <= 12: # close enough to circle
        sol.append(new_path)
        backtrack(ps-{p}, sol[-1][-1], sol, seen|lines|{*new_path}, best)
        sol.pop()

ps = polylines(m=4)
draw_polylines(ps)
# backtrack(set(ps))



  