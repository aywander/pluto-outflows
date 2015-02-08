Here is the new API for construction.  Note that I did not include any
file I/O stuff.   We did discuss using the Vector<Node> stuff for 
file I/O but I think this is largely orthogonal to our May goals.


dtg


/************************************/
//ebindexspace's construction interface can remain unchanged
//though I would love to remove the dx argument.
/************************************/
class  EBIndexSpace
{
public:  
  void
  define(const Box&             a_domain,
         const RealVect&        a_origin,
         const Real&            a_dx,
         const GeometryService& a_geoserver);
}

/************************************/
//GeometryService  uses the linear node structure to build an EBISBox
//we can also rip out all the access functions to priviate EBISBox data
//As far as I can tell, this is the whole class.
//we could put the translation from Vector<node> to EBISBox either
//here or in EBIndexSpace.
/************************************/

class GeometryService
{
public:
  virtual bool isRegular(const Box&      a_region,
                         const Box&      a_domain,
                         const RealVect& a_origin,
                         const Real&     a_dx) const = 0;

  virtual bool isCovered(const Box&      a_region,
                         const Box&      a_domain,
                         const RealVect& a_origin,
                         const Real&     a_dx) const = 0;

  // In valid region---all nodes are complete.
  // In complement of valid region inside ghost region,
  // the nodes only have arcs that point  into the 
  // valid region.
  virtual void fillGraph(Vector<Node>&   a_nodes,
                         const Box&      a_validRegion,
                         const Box&      a_ghostRegion,
                         const Box&      a_domain,
                         const RealVect& a_origin,
                         const Real&     a_dx) const = 0;
}
/*********************************************/
/*********************************************/
class Node
{
  enum TAG {Regular = 0, Covered, Irregular};

private:
  //if the flavor is regular or covered, none of the
  //other data is valid.
  TAG          m_flavor;
  //for an irregular node, this is a 1x1 box.
  //for covered or regular, the region can be larger
  Box          m_region;

  //data for irregular nodes
  //volume centered data
  Real         m_volFrac;
  Real         m_boundaryArea;

  RealVect     m_normal;
  RealVect     m_volCentroid
  RealVect     m_boundaryCentroid;

  //face centered irregular data
  //indicies into a_nodes to show connectivity
  Vector<int>      m_arcs;
  Vector<Real>     m_areaFracs;
  Vector<RealVect> m_faceCentroid;
  
}
/*********************************************/
/*********************************************/
