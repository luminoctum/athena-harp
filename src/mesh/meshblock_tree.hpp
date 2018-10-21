#ifndef MESHBLOCK_TREE_HPP
#define MESHBLOCK_TREE_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file meshblock_tree.hpp
//  \brief defines the LogicalLocation structure and MeshBlockTree class
//======================================================================================

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../defs.hpp"
#include "../bvals/bvals.hpp"

//--------------------------------------------------------------------------------------
//! \struct LogicalLocation
//  \brief stores logical location and level of meshblock

typedef struct LogicalLocation {
  long int lx1, lx2, lx3;
  int level;

  LogicalLocation() : lx1(-1), lx2(-1), lx3(-1), level(-1) {};

  // operators useful for sorting
  bool operator==(LogicalLocation &ll)
    { return ((ll.level==level) && (ll.lx1==lx1) && (ll.lx2==lx2) && (ll.lx3==lx3)); }
  static bool Lesser(const LogicalLocation &left, const LogicalLocation &right)
    { return left.level < right.level; };
  static bool Greater(const LogicalLocation & left, const LogicalLocation &right)
    { return left.level > right.level; };

} LogicalLocation;

//--------------------------------------------------------------------------------------
//! \class MeshBlockTree
//  \brief Objects are nodes in an AMR MeshBlock tree structure

class MeshBlockTree {
  friend class Mesh;
  friend class MeshBlock;
public:
  MeshBlockTree();
  MeshBlockTree(MeshBlockTree *parent, int ox, int oy, int oz);
  ~MeshBlockTree();

  // accessor
  MeshBlockTree* GetLeaf(int ox, int oy, int oz) {return pleaf[oz][oy][ox];}

  // functions
  void CreateRootGrid(long int nx, long int ny, long int nz, int nl);
  void AddMeshBlock(MeshBlockTree& root, LogicalLocation rloc, int dim,
       enum BoundaryFlag* mesh_bcs, long int rbx, long int rby, long int rbz,
       int rl, int &nnew);
  void AddMeshBlockWithoutRefine(LogicalLocation rloc, 
                                 long int rbx, long int rby, long int rbz, int rl);
  void Refine(MeshBlockTree& root, int dim, enum BoundaryFlag* mesh_bcs,
              long int rbx, long int rby, long int rbz, int rl, int &nnew);
  void Derefine(MeshBlockTree& root, int dim, enum BoundaryFlag* mesh_bcs,
              long int rbx, long int rby, long int rbz, int rl, int &ndel);
  MeshBlockTree* FindMeshBlock(LogicalLocation tloc);
  void CountMeshBlock(int& count);
  void GetMeshBlockList(LogicalLocation *list, int *pglist, int& count);
  MeshBlockTree* FindNeighbor(LogicalLocation myloc, int ox1, int ox2, int ox3,
                 enum BoundaryFlag* bcs, long int rbx, long int rby, long int rbz,
                 int rl, bool amrflag=false);

private:
  // data
  bool flag; // false: virtual node, has leaves; true: real node, is a leaf
  MeshBlockTree* pparent;
  MeshBlockTree* pleaf[2][2][2];
  LogicalLocation loc;
  int gid;
};

#endif // MESHBLOCK_TREE_HPP
