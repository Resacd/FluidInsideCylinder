#include "OutputFile.h"


OutputFile::OutputFile(int outputIdx)
    : m_outputIdx(outputIdx)
{
}

OutputFile::~OutputFile()
{
}



void OutputFile::WriteVTK_LocalData(SubDomain& subdomain)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    std::string folder_name = "Fluid_Variables_LVTK";

#ifdef _WIN32
    std::filesystem::create_directory(folder_name);
#else 
    mkdir(folder_name.c_str(), 0777);
#endif

    std::string file_name = folder_name + "/" + std::to_string(m_outputIdx) +
        "_rank_" + std::to_string(rank) + ".vtu";  // Use .vtu extension for unstructured grid

    vtkSmartPointer<vtkUnstructuredGrid> uGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();

    vtkSmartPointer<vtkDoubleArray> velocities = vtkSmartPointer<vtkDoubleArray>::New();
    velocities->SetNumberOfComponents(3);
    velocities->SetName("Velocity");
    vtkSmartPointer<vtkDoubleArray> pressures = vtkSmartPointer<vtkDoubleArray>::New();
    pressures->SetName("Pressure");
    vtkSmartPointer<vtkDoubleArray> velMags = vtkSmartPointer<vtkDoubleArray>::New();
    velMags->SetName("VelocityMagnitude");
    vtkSmartPointer<vtkDoubleArray> Density = vtkSmartPointer<vtkDoubleArray>::New();
    Density->SetName("Density");

    int pointId = 0;

    for (int z = 0; z < subdomain.my_Nz_; ++z) {
        for (int y = 0; y < subdomain.my_Ny_; ++y) {
            for (int x = 0; x < subdomain.my_Nx_; ++x) {
                auto& node = subdomain.m_lattice[x][y][z];
                if (node != nullptr && node->NodeType_ != NodeType::Halo /*node->NodeType_ == NodeType::Fluid*/) {
                    double velMag = sqrt(node->velocity[0] * node->velocity[0] +
                        node->velocity[1] * node->velocity[1] +
                        node->velocity[2] * node->velocity[2]);

                    std::vector<int> globalIndices = subdomain.Local2Global(node->x_position, node->y_position, node->z_position);
                    points->InsertNextPoint(globalIndices[0], globalIndices[1], globalIndices[2]);
                    velocities->InsertNextTuple(node->velocity.data());
                    pressures->InsertNextValue(node->Pressure(subdomain.m_velSet));
                    velMags->InsertNextValue(velMag);
                    Density->InsertNextValue(node->density);

                    cells->InsertNextCell(1);
                    cells->InsertCellPoint(pointId);
                    ++pointId;
                }
            }
        }
    }

    uGrid->SetPoints(points);
    uGrid->SetCells(VTK_VERTEX, cells);

    // Add arrays to the point data
    uGrid->GetPointData()->AddArray(velocities);
    uGrid->GetPointData()->AddArray(pressures);
    uGrid->GetPointData()->AddArray(velMags);
    uGrid->GetPointData()->AddArray(Density);

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(file_name.c_str());
    writer->SetInputData(uGrid);
    writer->Write();
}




//void OutputFile::WriteVTK_LocalData(SubDomain& subdomain)
//{
//    MPI_Comm comm = MPI_COMM_WORLD;
//    int rank, size;
//    MPI_Comm_rank(comm, &rank);
//    MPI_Comm_size(comm, &size);
//
//    std::string folder_name = "Fluid_Variables_LVTK";
//
//#ifdef _WIN32
//    std::filesystem::create_directory(folder_name);
//#else 
//    mkdir(folder_name.c_str(), 0777);
//#endif
//
//    std::string file_name = folder_name + "/" + std::to_string(m_outputIdx) +
//        "_rank_" + std::to_string(rank) + ".vtu";  // Use .vtu extension for unstructured grid
//
//    vtkSmartPointer<vtkUnstructuredGrid> uGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
//    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
//    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
//
//    vtkSmartPointer<vtkDoubleArray> velocities = vtkSmartPointer<vtkDoubleArray>::New();
//    velocities->SetNumberOfComponents(3);
//    velocities->SetName("Velocity");
//    vtkSmartPointer<vtkDoubleArray> pressures = vtkSmartPointer<vtkDoubleArray>::New();
//    pressures->SetName("Pressure");
//    vtkSmartPointer<vtkDoubleArray> velMags = vtkSmartPointer<vtkDoubleArray>::New();
//    velMags->SetName("VelocityMagnitude");
//    vtkSmartPointer<vtkDoubleArray> Density = vtkSmartPointer<vtkDoubleArray>::New();
//    Density->SetName("Density");
//
//    int pointId = 0;
//
//    for (int z = 0; z < subdomain.my_Nz_; ++z) {
//        for (int y = 0; y < subdomain.my_Ny_; ++y) {
//            for (int x = 0; x < subdomain.my_Nx_; ++x) {
//                auto& node = subdomain.m_lattice[x][y][z];
//                if (node != nullptr && node->NodeType_ != NodeType::Halo /*node->NodeType_ == NodeType::Fluid*/) {
//                    double velMag = sqrt(node->Velocity(subdomain.m_velSet)[0] * node->Velocity(subdomain.m_velSet)[0] +
//                        node->Velocity(subdomain.m_velSet)[1] * node->Velocity(subdomain.m_velSet)[1] +
//                        node->Velocity(subdomain.m_velSet)[2] * node->Velocity(subdomain.m_velSet)[2]);
//
//                    std::vector<int> globalIndices = subdomain.Local2Global(node->x_position, node->y_position, node->z_position);
//                    points->InsertNextPoint(globalIndices[0], globalIndices[1], globalIndices[2]);
//                    velocities->InsertNextTuple(node->Velocity(subdomain.m_velSet).data());
//                    pressures->InsertNextValue(node->Pressure(subdomain.m_velSet));
//                    velMags->InsertNextValue(velMag);
//                    Density->InsertNextValue(node->Density());
//
//                    cells->InsertNextCell(1);
//                    cells->InsertCellPoint(pointId);
//                    ++pointId;
//                }
//            }
//        }
//    }
//
//    uGrid->SetPoints(points);
//    uGrid->SetCells(VTK_VERTEX, cells);
//
//    // Add arrays to the point data
//    uGrid->GetPointData()->AddArray(velocities);
//    uGrid->GetPointData()->AddArray(pressures);
//    uGrid->GetPointData()->AddArray(velMags);
//    uGrid->GetPointData()->AddArray(Density);
//
//    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
//    writer->SetFileName(file_name.c_str());
//    writer->SetInputData(uGrid);
//    writer->Write();
//}

