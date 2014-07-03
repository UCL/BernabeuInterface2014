#include "mex.h"
#include <math.h>

#include <vtkVersion.h>
#include <vtkPoints.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkTubeFilter.h>
#include <vtkDataObject.h>
#include <vtkCleanPolyData.h>

#define acces_2d(matrix,rowIndex,columnIndex,numRows) (matrix[(columnIndex)*(numRows) + (rowIndex)])

double DistanceTwoVertices(double* vertices, unsigned numVertices, vtkIdType* edgeEnds)
{
    double x_diff = acces_2d(vertices, edgeEnds[0], 0, numVertices) - acces_2d(vertices, edgeEnds[1], 0, numVertices);
    double y_diff = acces_2d(vertices, edgeEnds[0], 1, numVertices) - acces_2d(vertices, edgeEnds[1], 1, numVertices);
    
    double dist = sqrt(x_diff*x_diff + y_diff*y_diff);
    return dist;
}

void WritePolyDataToDisk(const vtkSmartPointer<vtkPolyData> polyData, const std::string& fileName)
{
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(fileName.c_str());
    writer->SetInput(polyData);

    // Optional - set the mode. The default is binary.
    //writer->SetDataModeToBinary();
    //writer->SetDataModeToAscii();

    writer->Write();    
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
    // Check for proper number of arguments
    if (nrhs!=6) 
    { 
	    mexErrMsgTxt("Six input argument required."); 
    } 
    else if (nlhs>1) 
    {
	    mexErrMsgTxt("Too many output arguments provided."); 
    } 

    // Get access to MATLAB data through C data types
    double *vertices = mxGetPr(prhs[0]);
    unsigned num_vertices = mxGetM(prhs[0]);
    
    double *edges = mxGetPr(prhs[1]);
    unsigned num_edges = mxGetM(prhs[1]);

    double *radii = mxGetPr(prhs[2]);
    unsigned num_radii_values = mxGetM(prhs[2]);
    
    std::string dataset_name(mxArrayToString(prhs[3]));
    
    double pixel_to_um = mxGetScalar(prhs[4]);
    
    double radii_fudge_factor = mxGetScalar(prhs[5]);
    
    if (num_vertices != num_radii_values)
    {
        mexErrMsgTxt("Number of radii values provided differs from the number of vertices specified.");
    }
    
    // Create the output variable that will store the edge lengths and radii required for one of the histograms
    plhs[0] = mxCreateDoubleMatrix(num_edges, 2, mxREAL);
    double *edge_length_radii = mxGetPr(plhs[0]);
    
    // Create and fill a point container
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (unsigned vertex=0; vertex<num_vertices; vertex++)
    {
        points->InsertNextPoint(acces_2d(vertices, vertex, 0, num_vertices) * pixel_to_um, 
                                acces_2d(vertices, vertex, 1, num_vertices) * pixel_to_um, 
                                0);
    }
 
    // Create a polydata object and add the points to it.
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);

    // Create the array containing edge lengths 
    vtkSmartPointer<vtkDoubleArray> lengths_array = vtkSmartPointer<vtkDoubleArray>::New();
    lengths_array->SetNumberOfComponents(1);
    lengths_array->SetName("EdgeLengths");
    
    // Allocate memory and add the definition of edges to the polydata object
    polydata->Allocate();
    vtkIdType edge_ends[2];    
    for (unsigned edge_id=0; edge_id<num_edges; edge_id++)
    {
        // The definition of edges comes from MATLAB and and the vertices are one-based indexed. Offset them.
        edge_ends[0] = acces_2d(edges, edge_id, 0, num_edges) - 1; 
        edge_ends[1] = acces_2d(edges, edge_id, 1, num_edges) - 1;
        
        polydata->InsertNextCell(VTK_LINE, 2, edge_ends);
        double edge_length = DistanceTwoVertices(vertices, num_vertices, edge_ends) * pixel_to_um;
        lengths_array->InsertNextValue(edge_length);
        
        acces_2d(edge_length_radii, edge_id, 0, num_edges) = edge_length;
        acces_2d(edge_length_radii, edge_id, 1, num_edges) = radii_fudge_factor * pixel_to_um * (radii[edge_ends[0]] + radii[edge_ends[1]])/2;
    }
    //polydata->GetCellData()->AddArray(lengths_array);
    
    // Create and fill the array containing the plexus radii at each vertex.
    vtkSmartPointer<vtkDoubleArray> radii_array = vtkSmartPointer<vtkDoubleArray>::New();
    radii_array->SetNumberOfComponents(1);
    radii_array->SetName("Radius");    
    for (unsigned vertex_id=0; vertex_id<num_vertices; vertex_id++)
    {
        radii_array->InsertNextValue(radii_fudge_factor * radii[vertex_id] * pixel_to_um);
    }
    polydata->GetPointData()->AddArray(radii_array);
    polydata->GetPointData()->SetActiveScalars("Radius");
     
    WritePolyDataToDisk(polydata, dataset_name + ".vtp");
}
