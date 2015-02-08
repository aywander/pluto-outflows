/************************************/
/************************************/
/************************************/
class EBISBox
{
public:
  //public interface can remain the same except for construction

private:
  
  //region over which this EBISBox exists
  Box m_region;
  
  //domain of the calculation at this level of refinement
  Box m_domain;

  //graph connectivity
  EBGraph m_graph;

  //volume fraction
  BaseIVFAB<Real> m_volFrac;

  //normal to the boundary (spacedim vars)
  BaseIVFAB<Real> m_normal;

  //centroid of the volume (spacedim vars)
  BaseIVFAB<Real> m_volCentroid;

  //area of the boundry 
  BaseIVFAB<Real> m_bndryArea;

  //centroid of the boundary (spacedim vars)
  BaseIVFAB<Real> m_bndryCentroid;

  //area fraction in each direction
  BaseIFFAB<Real> m_areaFrac[SpaceDim];

  //centroid of the faces  in each direction
  BaseIFFAB<Real> m_faceCentroid[SpaceDim];
}
/************************************/
/************************************/
/************************************/

//this is the graph data held for each irregular node
class GraphNodeImplem
{
  vector<int>      m_arcs[2*SpaceDim];
  int              m_coarserNode;
  Vector<VolIndex> m_finerNodes;
}

/************************************/
//i made this a doubly-linked list 
class GraphNode
{
private:
  //the connectivity data at this point
  //if m_thisNode == 0, node is covered 
  //if m_thisNode == 1, node is regular
  GraphNodeImplem* m_thisNode;

  //the next node in this cell.  
  GraphNode*       m_nextNode;

  //the previous node to this in the cell
  GraphNode*       m_prevNode;
}

/************************************/
/************************************/
/************************************/
class EBGraph
{
public: 
  //only graph questions can be answered by this class
  /// I list the ones from EBISBox that I think are relevant
  Vector<VolIndex> getVoFs(const IntVect& a_iv) const;
  bool isRegular(const IntVect& a_iv) const;
  bool isAllRegular() const;
  bool isAllCovered() const;
  bool isIrregular(const IntVect& a_iv) const;
  bool isCovered(const Box& a_box) const;
  bool isRegular(const Box& a_box) const;
  Vector<FaceIndex> getFaces(const VolIndex& a_vof, 
                             const int& a_idir, 
                             const Side::LoHiSide& a_sd) const;
  Vector<VolIndex> refine(const VolIndex& a_coarVoF) const;
  VolIndex coarsen(const VolIndex& a_fineVoF) const;

  enum TAG {AllRegular = 0, AllCovered, HasIrregular};
private:
  //   if this is allregular or allcovered, the rest of stuff below is undefined.
  TAG m_tag;

  //the graph description of the box.
  //if the pointer GraphNode.m_thisNode == 0, the cell is covered.
  //if the pointer GraphNode.m_thisNode == 1, the cell is regular
  //otherwise, the pointer goes to the node of the graph
  BaseFab<GraphNode>  m_graphData;
}
/************************************/
/************************************/
/************************************/
