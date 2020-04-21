// MAUS WARNING: THIS IS LEGACY CODE.
#ifndef Mesh_hh
#define Mesh_hh

#include <vector>
#include <algorithm>
#include <math.h>
#include <ostream>

#include "src/legacy/Interface/Meshing/Mesh.hh"
#include "src/legacy/Interface/Meshing/ThreeDGrid.hh"

//Holds grid info for 2dTo1d interpolation algorithm; 
//also controls memory usage (essentially a boost::shared_pointer but not as elegent);
class TwoDGrid : public Mesh
{
public:
    class Iterator;
    //get value
    Mesh * Clone() {return new TwoDGrid(*this);}
    Mesh * Dual () {return NULL;}
    //ERROR SHOULD NOT ALLOW MESH WITH 1 GRID POINT (causes error in LowerBound)
    TwoDGrid();
    TwoDGrid(double dX, double dY, double minX, double minY, int numberOfXCoords, int numberOfYCoords);
    TwoDGrid(int xSize, const double *x, int ySize, const double *y);
    TwoDGrid(std::vector<double> x, std::vector<double> y);
    ~TwoDGrid();
    //get coordinate; round bracket indexing starts at 1 goes to nCoords
    inline double& x    (const int& i) {return _x[i-1];}
    inline double& y    (const int& j) {return _y[j-1];}
    inline const double& x    (const int& i) const {return _x[i-1];}
    inline const double& y    (const int& j) const {return _y[j-1];}
    inline int     xSize() const       {return int(_x.size());}
    inline int     ySize() const       {return int(_y.size());}

    std::vector<double> xVector()      {return std::vector<double>(_x);}
    std::vector<double> yVector()      {return std::vector<double>(_y);}

    double* newXArray()                {double *x = new double[_x.size()]; for(unsigned int i=0; i<_x.size(); i++) x[i]=_x[i]; return x;}
    double* newYArray()                {double *y = new double[_y.size()]; for(unsigned int i=0; i<_y.size(); i++) y[i]=_y[i]; return y;}

    //if you are sure the grid has constant spacing ConstSpacing is quicker; VarSpacing in either case
    //return lower bound on x for insertion i.e. 0 if x < x[1], 1 if  x[1] < x < x[2], ... , x.size()-1 if x > x.back()
    inline void xLowerBound            (const double& x, int& xIndex) const 
    {if(_constantSpacing) xIndex = static_cast<int>(std::floor( (x - _x[0])/(_x[1]-_x[0]) )); else xIndex = std::lower_bound(_x.begin(), _x.end(), x)-_x.begin()-1;}
    inline void yLowerBound            (const double& y, int& yIndex) const
    {if(_constantSpacing) yIndex = static_cast<int>(std::floor( (y - _y[0])/(_y[1]-_y[0]) )); else yIndex = std::lower_bound(_y.begin(), _y.end(), y)-_y.begin()-1;}
    inline void LowerBound             (const double& x, int& xIndex, const double& y, int& yIndex) const
    {xLowerBound(x, xIndex); yLowerBound(y, yIndex);}

    inline double    MinX()            const {return _x[0];}
    inline double    MaxX()            const {return _x[_xSize-1];}
    inline double    MinY()            const {return _y[0];}
    inline double    MaxY()            const {return _y[_ySize-1];}

    void    Add   (VectorMap* map); //add *map if it has not already been added
    void    Remove(VectorMap* map); //remove *map if it exists; delete this if there are no more VectorMaps

    void    SetX(int nXCoords, double * x) {_x = std::vector<double>(x,x+nXCoords);}
    void    SetY(int nYCoords, double * y) {_y = std::vector<double>(y,y+nYCoords);}

    Mesh::Iterator Begin() const;
    Mesh::Iterator End()   const;
    //position of mesh point
    virtual void   Position(const Mesh::Iterator& it, double * position) const;
    //position of centre of mesh point at it and opposite diagonal in +ve index direction
    virtual void   CentrePosition(const Mesh::Iterator& it, double * position) const;
    int            PositionDimension() const {return 2;}
    int            ToInteger(const Mesh::Iterator& lhs) const {return (lhs.State()[0]-1)*(_ySize)+lhs.State()[1]-1;}

    void           SetConstantSpacing(bool spacing) {_constantSpacing = spacing;}
    void           SetConstantSpacing();
    bool           GetConstantSpacing() const       {return _constantSpacing; }

    Mesh::Iterator Nearest(const double* position) const;

protected:

    //Change position
    Mesh::Iterator& AddEquals(Mesh::Iterator& lhs, int difference) const;
    Mesh::Iterator& SubEquals(Mesh::Iterator& lhs, int difference) const;
    Mesh::Iterator& AddEquals(Mesh::Iterator& lhs, const Mesh::Iterator& rhs) const;
    Mesh::Iterator& SubEquals(Mesh::Iterator& lhs, const Mesh::Iterator& rhs) const;
    Mesh::Iterator& AddOne   (Mesh::Iterator& lhs) const;
    Mesh::Iterator& SubOne   (Mesh::Iterator& lhs) const;
    //Check relative position
    bool IsGreater(const Mesh::Iterator& lhs, const Mesh::Iterator& rhs) const;

private:
    std::vector<double>     _x;
    std::vector<double>     _y;
    int                     _xSize;
    int                     _ySize;
    std::vector<VectorMap*> _maps;
    bool                    _constantSpacing;

    friend Mesh::Iterator  operator++(Mesh::Iterator& lhs, int);
    friend Mesh::Iterator  operator--(Mesh::Iterator& lhs, int);
    friend Mesh::Iterator& operator++(Mesh::Iterator& lhs);
    friend Mesh::Iterator& operator--(Mesh::Iterator& lhs);
    friend Mesh::Iterator  operator- (const Mesh::Iterator& lhs, const Mesh::Iterator& rhs);
    friend Mesh::Iterator  operator+ (const Mesh::Iterator& lhs, const Mesh::Iterator& rhs);
    friend Mesh::Iterator& operator-=(Mesh::Iterator& lhs,       const Mesh::Iterator& rhs);
    friend Mesh::Iterator& operator+=(Mesh::Iterator& lhs,       const Mesh::Iterator& rhs);
    friend Mesh::Iterator  operator- (const Mesh::Iterator&, const int&);
    friend Mesh::Iterator  operator+ (const Mesh::Iterator&, const int&);
    friend Mesh::Iterator& operator-=(Mesh::Iterator&,       const int&);
    friend Mesh::Iterator& operator+=(Mesh::Iterator&,       const int&);

    friend bool operator==(const Mesh::Iterator& lhs, const Mesh::Iterator& rhs);
    friend bool operator!=(const Mesh::Iterator& lhs, const Mesh::Iterator& rhs);
    friend bool operator>=(const Mesh::Iterator& lhs, const Mesh::Iterator& rhs);
    friend bool operator<=(const Mesh::Iterator& lhs, const Mesh::Iterator& rhs);
    friend bool operator< (const Mesh::Iterator& lhs, const Mesh::Iterator& rhs);
    friend bool operator> (const Mesh::Iterator& lhs, const Mesh::Iterator& rhs);

};

//Holds grid info for some vector mapping; 
//also controls memory usage (essentially a boost::shared_pointer but not as elegent);
class NDGrid : public Mesh
{
public:
    class Iterator;

    Mesh * Clone() {return new NDGrid(*this);}
    Mesh * Dual () {return NULL;}
    //get value
    NDGrid();
    NDGrid(int nDimensions, int* size, double* spacing, double* min);
    NDGrid(std::vector<int> size, std::vector<const double *> gridCoordinates);
    NDGrid(std::vector< std::vector<double> > gridCoordinates);
    ~NDGrid() {;}

    //get coordinate; round bracket indexing starts at 1 goes to nCoords
    inline double&       coord    ( const int& index, const int& dimension)       {return _coord[dimension][index-1];}
    inline const double& coord    ( const int& index, const int& dimension) const {return _coord[dimension][index-1];}
    inline int           size     ( const int& dimension ) const                  {return _coord[dimension].size();}
    std::vector<double>  coordVector( const int& dimension ) const {return _coord[dimension];}
    double*              newCoordArray  ( const int& dimension ) const;

    //if you are sure the grid has constant spacing ConstSpacing is quicker; VarSpacing in either case
    //indexing starts at 0 and goes to nCoords+1
    inline void    coordLowerBound (const double& x, const int& dimension, int& xIndex)       const
    {
        if(_constantSpacing) xIndex = static_cast<int>(std::floor( (x - _coord[dimension][0])/(_coord[dimension][1]-_coord[dimension][0]) )); 
        else xIndex = std::lower_bound(_coord[dimension].begin(), _coord[dimension].end(), x)-_coord[dimension].begin()-1;
    }

    inline void    LowerBound      (const std::vector<double>& pos, std::vector<int>& xIndex) const
    {xIndex = std::vector<int>(pos.size()); for(unsigned int i=0; i<pos.size(); i++) coordLowerBound(pos[i], i, xIndex[i]);}
    inline double  Min             (const int& dimension) const {return _coord[dimension][0];}
    inline double  Max             (const int& dimension) const {return _coord[dimension][_coord[dimension].size()-1];}

    void           SetCoord        ( int dimension, int nCoords, double * x) {_coord[dimension] = std::vector<double>(x, x+nCoords);}

    Mesh::Iterator Begin() const  {return Mesh::Iterator(std::vector<int>(_coord.size(), 1), this);}
    Mesh::Iterator End()   const  {std::vector<int> end(_coord.size(), 1); end[0] = _coord[0].size()+1; return Mesh::Iterator(end, this);}

    void           Position(const Mesh::Iterator& it, double * position) const {for(unsigned int i=0; i<it.State().size(); i++) position[i] = _coord[i][it[i]-1];}
    int            PositionDimension() const {return _coord.size();}

    bool           GetConstantSpacing() const       {return _constantSpacing;}
    void           SetConstantSpacing(bool spacing) {_constantSpacing = spacing;}
    void           SetConstantSpacing();
    int            ToInteger         (const Mesh::Iterator& lhs) const;
    Mesh::Iterator Nearest           (const double* position)          const;

protected:

    //Change position
    virtual Mesh::Iterator& AddEquals(Mesh::Iterator& lhs, int difference) const;
    virtual Mesh::Iterator& SubEquals(Mesh::Iterator& lhs, int difference) const;
    virtual Mesh::Iterator& AddEquals(Mesh::Iterator& lhs, const Mesh::Iterator& rhs) const;
    virtual Mesh::Iterator& SubEquals(Mesh::Iterator& lhs, const Mesh::Iterator& rhs) const;
    virtual Mesh::Iterator& AddOne   (Mesh::Iterator& lhs) const;
    virtual Mesh::Iterator& SubOne   (Mesh::Iterator& lhs) const;
    //Check relative position
    virtual bool IsGreater(const Mesh::Iterator& lhs, const Mesh::Iterator& rhs) const;

private:
    std::vector< std::vector<double> > _coord;
    std::vector<VectorMap*>            _maps;
    bool                               _constantSpacing;

    friend Mesh::Iterator  operator++(Mesh::Iterator& lhs, int);
    friend Mesh::Iterator  operator--(Mesh::Iterator& lhs, int);
    friend Mesh::Iterator& operator++(Mesh::Iterator& lhs);
    friend Mesh::Iterator& operator--(Mesh::Iterator& lhs);
    friend Mesh::Iterator  operator- (const Mesh::Iterator& lhs, const Mesh::Iterator& rhs);
    friend Mesh::Iterator  operator+ (const Mesh::Iterator& lhs, const Mesh::Iterator& rhs);
    friend Mesh::Iterator& operator-=(Mesh::Iterator& lhs,       const Mesh::Iterator& rhs);
    friend Mesh::Iterator& operator+=(Mesh::Iterator& lhs,       const Mesh::Iterator& rhs);

    friend bool operator==(const Mesh::Iterator& lhs, const Mesh::Iterator& rhs);
    friend bool operator!=(const Mesh::Iterator& lhs, const Mesh::Iterator& rhs);
    friend bool operator>=(const Mesh::Iterator& lhs, const Mesh::Iterator& rhs);
    friend bool operator<=(const Mesh::Iterator& lhs, const Mesh::Iterator& rhs);
    friend bool operator< (const Mesh::Iterator& lhs, const Mesh::Iterator& rhs);
    friend bool operator> (const Mesh::Iterator& lhs, const Mesh::Iterator& rhs);
};

#endif
