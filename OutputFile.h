#pragma once
#include "Node.h"
#include <fstream>
#include <string>
#include <mpi.h>
#include "Utility.h"
#include <filesystem>

#ifndef _WIN32
#include <sys/stat.h>
#else
#include <windows.h>
#endif

#include <memory> 

#include "Struct.h"
#include "HaloLayer.h"
#include "Subdomain.h"

#include <vtk-9.2/vtkSmartPointer.h>
#include <vtk-9.2/vtkUnstructuredGrid.h>
#include <vtk-9.2/vtkPoints.h>
#include <vtk-9.2/vtkDoubleArray.h>
#include <vtk-9.2/vtkPointData.h>
#include <vtk-9.2/vtkDelaunay3D.h>
#include <vtk-9.2/vtkXMLUnstructuredGridWriter.h>

#include <vtk-9.2/vtkPolyData.h>
#include <vtk-9.2/vtkCellArray.h>
#include <vtk-9.2/vtkXMLPolyDataWriter.h>
#include <vtk-9.2/vtkFloatArray.h>
#include <vtk-9.2/vtkTriangleFilter.h>
#include <vtk-9.2/vtkPolyDataWriter.h>



class OutputFile
{
private:

	std::mutex m_fileMutex;
	int m_outputIdx;
	std::mutex m_dataMutex;
	Utility Utility_;

public:

	OutputFile(int outputIdx);
	~OutputFile();


	void WriteVTK_LocalData(SubDomain &subdomain);
	void WriteWallShearStressVTK(SubDomain& subdomain, int iter);

};
