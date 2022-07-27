#include "vtkIntArray.h"
#include "vtkDoubleArray.h"
#include "vtkPolyData.h"
#include "vtkCellData.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkTetra.h"
#include "vtkGenericCell.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLGenericDataObjectReader.h"
#include <stdlib.h>
void add_int_PointData( vtkPointSet * const &grid_w,
    const std::vector<int> &ptdata, const std::string &dataname );

void add_int_CellData( vtkPointSet * const &grid_w,
    const std::vector<int> &cldata, const std::string &dataname );

int main(int argc, char * argv[])
{
  int n = 2, id;
  double x,y,z;

  const bool isXML = false;
  const std::string filename("pro4");

  vtkUnstructuredGrid * grid_w = vtkUnstructuredGrid::New();

  vtkPoints * ppt = vtkPoints::New();
  ppt->SetDataTypeToDouble();
 
  for (int i = 0; i < n+1; i++)
    for (int j = 0; j < n+1; j++)
      for (int k = 0; k < n+1; k++){
	id = (n+1)*(n+1)*k + (n+1)*j + i;
	x = 1.0/n*i;
	y = 1.0/n*j;
	z = 1.0/n*k;

  	ppt -> InsertPoint(id, x, y, z);
      }
  grid_w -> SetPoints(ppt);
  ppt -> Delete();

  std::vector<int> IEN;

  vtkCell * cl = vtkTetra::New();
  for(int i=0; i<n; ++i)
    for(int j=0; j<n; ++j)
      for(int k=0; k<n; ++k){
	id = (n+1)*(n+1)*k + (n+1)*j + i;
	IEN = {id, id+1, id+n+1, id+(n+1)*(n+1)+1,
	       id, id+n+1, id+(n+1)*(n+1), id+(n+1)*(n+1)+1,
               id+n+1,id+(n+1)*(n+1), id+(n+1)*(n+1)+1, id+(n+1)*(n+2),
	       id+1, id+n+1, id+n+2, id+(n+1)*(n+1)+1,
               id+n+1, id+n+2, id+(n+1)*(n+1)+1, id+(n+1)*(n+2)+1,
               id+n+1, id+(n+1)*(n+2), id+(n+1)*(n+1)+1, id+(n+1)*(n+2)+1};

	for(int ii=0; ii<6; ++ii){
          for(int jj=0; jj<4; ++jj)
      	    cl->GetPointIds()->SetId( jj, IEN[4*ii + jj] );
    	  grid_w->InsertNextCell( cl->GetCellType(), cl->GetPointIds() );
  	}
      }

  cl -> Delete();


  /*std::vector<int> node_id = {1,2,3,4,5,6,7,8};
  add_int_PointData(grid_w, node_id, "GlobalNodeID");

  std::vector<int> cell_id = {1,2,3,4,5,6};
  add_int_CellData(grid_w, cell_id, "GlobalCellID");*/

  if( isXML )
  {
    vtkXMLUnstructuredGridWriter * writer = vtkXMLUnstructuredGridWriter::New();
    std::string name_to_write(filename);
    name_to_write.append(".vtu");
    writer -> SetFileName( name_to_write.c_str() );

    writer->SetInputData(grid_w);
    writer->Write();
    writer->Delete();
  }
  else
  {
    vtkUnstructuredGridWriter * writer = vtkUnstructuredGridWriter::New();
    std::string name_to_write(filename);
    name_to_write.append(".vtk");
    writer -> SetFileName( name_to_write.c_str() );

    writer->SetInputData(grid_w);
    writer->Write();
    writer->Delete();
  }


  grid_w->Delete();
  return 0;
}

void add_int_PointData( vtkPointSet * const &grid_w,
    const std::vector<int> &ptdata, const std::string &dataname )
{
  vtkIntArray * data = vtkIntArray::New();
  data -> SetNumberOfComponents(1);
  data -> SetName(dataname.c_str());

  for(unsigned int ii=0; ii<ptdata.size(); ++ii)
    data -> InsertComponent(ii, 0, ptdata[ii]);

  grid_w -> GetPointData() -> AddArray( data );
  data -> Delete();
}

void add_int_CellData( vtkPointSet * const &grid_w,
    const std::vector<int> &cldata, const std::string &dataname )
{
  vtkIntArray * data = vtkIntArray::New();
  data -> SetNumberOfComponents(1);
  data -> SetName(dataname.c_str());

  for(unsigned int ii=0; ii<cldata.size(); ++ii)
    data -> InsertComponent(ii, 0, cldata[ii]);

  grid_w -> GetCellData() -> AddArray( data );
  data -> Delete();
}

// EOF
