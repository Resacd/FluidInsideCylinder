#include "Subdomain.h"

SubDomain::SubDomain(Geometry& Geometry) {

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank_);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs_);
    Geometry_ = Geometry;
    m_velSet = Geometry_.velSet_;
    gridSize_ = Geometry_.GeomSize;
    Create_cart_comm();
    CreateNodes();
    Connectivity();
    GetInvolevedLat();
    Initialize();
    Need_Comm();
    MPI_Comm_free(&cart_comm_); 

    std::string wallFileName = "wallLattice" + std::to_string(myrank_) + ".vtk";
    std::string inletFileName = "inletLattice" + std::to_string(myrank_) + ".vtk";
    std::string outletFileName = "outletLattice" + std::to_string(myrank_) + ".vtk";


    WriteBoundaryLatticeToVTK(inletFileName.c_str(),"Inlet");
    WriteBoundaryLatticeToVTK(outletFileName.c_str(),"Outlet");
    WriteBoundaryLatticeToVTK(wallFileName.c_str(), "Wall");
    WriteFluidLatticeToVTK();

}


SubDomain::~SubDomain() {

}



void SubDomain::getNeighborInfo(int neighbor_rank, Neighbors& neighbor) {
    neighbor.rank = neighbor_rank;
    if (neighbor_rank != MPI_PROC_NULL && neighbor_rank >= 0) {
        MPI_Cart_coords(cart_comm_, neighbor_rank, 3, neighbor.coords);
    }
    else {
        neighbor.coords[0] = -1; 
        neighbor.coords[1] = -1; 
        neighbor.coords[2] = -1; 
    }
}


//void SubDomain::getNeighborInfo(int neighbor_rank, Neighbors& neighbor) {
//    neighbor.rank = neighbor_rank;
//    if (neighbor_rank != MPI_PROC_NULL) {
//        MPI_Cart_coords(cart_comm_, neighbor_rank, 3, neighbor.coords);
//    }
//}



bool SubDomain::is_valid_coordinate(int coords[]) {
    return coords[0] >= 0 && coords[0] < dims_[0] && coords[1] >= 0 && coords[1] < dims_[1] && coords[2] >= 0 && coords[2] < dims_[2];
}

int SubDomain::get_rank(MPI_Comm comm, int coords[]) {
    int rank;
    MPI_Cart_rank(comm, coords, &rank);
    return rank;
}

void SubDomain::Create_cart_comm() {
    int const ndims = 3;
    dims_ = { 0, 0, 0 };
    int periods[ndims] = { 0, 0, 0 };
    int reorder = 1;
    coords_ = { 0, 0, 0 };
    int remain_dims[3] = { 1, 1, 1 };
    int Nx = gridSize_[0];
    int Ny = gridSize_[1];
    int Nz = gridSize_[2];
    MPI_Dims_create(num_procs_, ndims, dims_.data());
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims_.data(), periods, reorder, &cart_comm_);

    MPI_Cart_shift(cart_comm_, 0, 1, &my_left_, &my_right_);
    MPI_Cart_shift(cart_comm_, 1, 1, &my_front_, &my_back_);
    MPI_Cart_shift(cart_comm_, 2, 1, &my_down_, &my_up_);
    MPI_Comm_rank(cart_comm_, &myrank_);
    MPI_Cart_coords(cart_comm_, myrank_, ndims, coords_.data());
    int cart_coords[ndims];
    MPI_Cart_get(cart_comm_, ndims, dims_.data(), periods, cart_coords);

    neighbors_.resize(18);

    getNeighborInfo(my_left_, neighbors_[0]);   
    getNeighborInfo(my_right_, neighbors_[1]);  
    getNeighborInfo(my_front_, neighbors_[2]);  
    getNeighborInfo(my_back_, neighbors_[3]);   
    getNeighborInfo(my_up_, neighbors_[4]);     
    getNeighborInfo(my_down_, neighbors_[5]);   

    int right_up_coords[3] = { coords_[0] + 1, coords_[1], coords_[2] + 1 };
    int right_down_coords[3] = { coords_[0] + 1, coords_[1], coords_[2] - 1 };
    int left_up_coords[3] = { coords_[0] - 1, coords_[1], coords_[2] + 1 };
    int left_down_coords[3] = { coords_[0] - 1, coords_[1], coords_[2] - 1 };
    int up_front_coords[3] = { coords_[0], coords_[1] - 1, coords_[2] + 1 };
    int up_back_coords[3] = { coords_[0], coords_[1] + 1, coords_[2] + 1 };
    int down_front_coords[3] = { coords_[0], coords_[1] - 1, coords_[2] - 1 };
    int down_back_coords[3] = { coords_[0], coords_[1] + 1, coords_[2] - 1 };
    int right_front_coords[3] = { coords_[0] + 1, coords_[1] - 1, coords_[2] };
    int right_back_coords[3] = { coords_[0] + 1, coords_[1] + 1, coords_[2] };
    int left_front_coords[3] = { coords_[0] - 1, coords_[1] - 1, coords_[2] };
    int left_back_coords[3] = { coords_[0] - 1, coords_[1] + 1, coords_[2] };

    right_up_ = is_valid_coordinate(right_up_coords) ? get_rank(cart_comm_, right_up_coords) : -1;
    right_down_ = is_valid_coordinate(right_down_coords) ? get_rank(cart_comm_, right_down_coords) : -1;
    left_up_ = is_valid_coordinate(left_up_coords) ? get_rank(cart_comm_, left_up_coords) : -1;
    left_down_ = is_valid_coordinate(left_down_coords) ? get_rank(cart_comm_, left_down_coords) : -1;
    up_front_ = is_valid_coordinate(up_front_coords) ? get_rank(cart_comm_, up_front_coords) : -1;
    up_back_ = is_valid_coordinate(up_back_coords) ? get_rank(cart_comm_, up_back_coords) : -1;
    down_front_ = is_valid_coordinate(down_front_coords) ? get_rank(cart_comm_, down_front_coords) : -1;
    down_back_ = is_valid_coordinate(down_back_coords) ? get_rank(cart_comm_, down_back_coords) : -1;
    right_front_ = is_valid_coordinate(right_front_coords) ? get_rank(cart_comm_, right_front_coords) : -1;
    right_back_ = is_valid_coordinate(right_back_coords) ? get_rank(cart_comm_, right_back_coords) : -1;
    left_front_ = is_valid_coordinate(left_front_coords) ? get_rank(cart_comm_, left_front_coords) : -1;
    left_back_ = is_valid_coordinate(left_back_coords) ? get_rank(cart_comm_, left_back_coords) : -1;


    getNeighborInfo(right_up_, neighbors_[6]);     
    getNeighborInfo(right_down_, neighbors_[7]);   
    getNeighborInfo(left_up_, neighbors_[8]);      
    getNeighborInfo(left_down_, neighbors_[9]);    
    getNeighborInfo(up_front_, neighbors_[10]);    
    getNeighborInfo(up_back_, neighbors_[11]);     
    getNeighborInfo(down_front_, neighbors_[12]);  
    getNeighborInfo(down_back_, neighbors_[13]);   
    getNeighborInfo(right_front_, neighbors_[14]); 
    getNeighborInfo(right_back_, neighbors_[15]);  
    getNeighborInfo(left_front_, neighbors_[16]);  
    getNeighborInfo(left_back_, neighbors_[17]);   

    halo_size = 1;
    my_Nx_ = Nx / dims_[0];
    my_Ny_ = Ny / dims_[1];
    my_Nz_ = Nz / dims_[2];
    start_x_ = coords_[0] * my_Nx_;
    start_y_ = coords_[1] * my_Ny_;
    start_z_ = coords_[2] * my_Nz_;
    end_x_ = start_x_ + my_Nx_ - 1;
    end_y_ = start_y_ + my_Ny_ - 1;
    end_z_ = start_z_ + my_Nz_ - 1;

    if (Nx % dims_[0] > 0 && coords_[0] == dims_[0] - 1) {
        my_Nx_ += Nx % dims_[0];
    }
    if (Ny % dims_[1] > 0 && coords_[1] == dims_[1] - 1) {
        my_Ny_ += Ny % dims_[1];
    }
    if (Nz % dims_[2] > 0 && coords_[2] == dims_[2] - 1) {
        my_Nz_ += Nz % dims_[2];
    }


    if (dims_[0] > 1) {
        if (coords_[0] == 0) {
            my_Nx_ += halo_size;
            end_x_ = start_x_ + my_Nx_ - 1;
        }
        else if (coords_[0] == dims_[0] - 1) {
            my_Nx_ += halo_size;
            start_x_ = start_x_ + 2 * halo_size * coords_[0] - halo_size;
            end_x_ = start_x_ + my_Nx_ - 1;
        }
        else {
            my_Nx_ += 2 * halo_size;
            start_x_ = start_x_ + 2 * halo_size * coords_[0] - halo_size;
            end_x_ = start_x_ + my_Nx_ - 1;
        }
    }

    if (dims_[1] > 1) {
        if (coords_[1] == 0) {
            my_Ny_ += halo_size;
            end_y_ = start_y_ + my_Ny_ - 1;
        }
        else if (coords_[1] == dims_[1] - 1) {
            my_Ny_ += halo_size;
            start_y_ = start_y_ + 2 * halo_size * coords_[1] - halo_size;
            end_y_ = start_y_ + my_Ny_ - 1;
        }
        else {
            my_Ny_ += 2 * halo_size;
            start_y_ = start_y_ + 2 * halo_size * coords_[1] - halo_size;
            end_y_ = start_y_ + my_Ny_ - 1;
        }
    }

    if (dims_[2] > 1) {
        if (coords_[2] == 0) {
            my_Nz_ += halo_size;
            end_z_ = start_z_ + my_Nz_ - 1;
        }
        else if (coords_[2] == dims_[2] - 1) {
            my_Nz_ += halo_size;
            start_z_ = start_z_ + 2 * halo_size * coords_[2] - halo_size;
            end_z_ = start_z_ + my_Nz_ - 1;
        }
        else {
            my_Nz_ += 2 * halo_size;
            start_z_ = start_z_ + 2 * halo_size * coords_[2] - halo_size;
            end_z_ = start_z_ + my_Nz_ - 1;
        }
    }

    NewgridSize_.resize(3);

    NewgridSize_[0] = (dims_[0] - 1) * 2 * halo_size + gridSize_[0];
    NewgridSize_[1] = (dims_[1] - 1) * 2 * halo_size + gridSize_[1];
    NewgridSize_[2] = (dims_[2] - 1) * 2 * halo_size + gridSize_[2];

    for (int p = 0; p < num_procs_; p++) {
        if (myrank_ == p) {

            //printf("Processor %d:\n", myrank_);
            //printf("my_Nx_ = %d, my_Ny_ = %d, my_Nz_ = %d\n", my_Nx_, my_Ny_, my_Nz_);
            //printf("start_x_ = %d, end_x_ = %d\n", start_x_, end_x_);
            //printf("start_y_ = %d, end_y_ = %d\n", start_y_, end_y_);
            //printf("start_z_ = %d, end_z_ = %d\n", start_z_, end_z_);
            //fflush(stdout);
        }
        MPI_Barrier(cart_comm_);
    }


}


bool SubDomain::IsInSubDomain(int coord_x, int coord_y, int coord_z) {

    if (coord_x >= 0 && coord_x < my_Nx_
        && coord_y >= 0 && coord_y < my_Ny_
        && coord_z >= 0 && coord_z < my_Nz_)
        return true;
    else
        return false;

}

std::vector<int>SubDomain::Local2Global(int coord_x, int coord_y, int coord_z) {

    std::vector<int> L2G;
    int global_coordX = coord_x;
    int global_coordY = coord_y;
    int global_coordZ = coord_z;

    if (dims_[0] > 1) {
        if (coords_[0] == 0) {
            global_coordX = coord_x;
        }
        else if (coords_[0] == dims_[0] - 1) {

            global_coordX = start_x_ + coord_x - 2 * coords_[0] * halo_size;
        }
        else {

            global_coordX = start_x_ + coord_x - 2 * coords_[0] * halo_size;

        }
    }


    if (dims_[1] > 1) {
        if (coords_[1] == 0) {
            ;
            global_coordY = coord_y;
        }
        else if (coords_[1] == dims_[1] - 1) {

            global_coordY = start_y_ + coord_y - 2 * coords_[1] * halo_size;
        }
        else {

            global_coordY = start_y_ + coord_y - 2 * coords_[1] * halo_size;

        }
    }

    if (dims_[2] > 1) {
        if (coords_[2] == 0) {
            ;
            global_coordZ = coord_z;
        }
        else if (coords_[2] == dims_[2] - 1) {

            global_coordZ = start_z_ + coord_z - 2 * coords_[2] * halo_size;
        }
        else {

            global_coordZ = start_z_ + coord_z - 2 * coords_[2] * halo_size;

        }
    }

    L2G = { global_coordX,global_coordY,global_coordZ };
    return L2G;

}


bool SubDomain::IsHaloLayer(int coordX, int coordY, int coordZ) {

    if (my_left_ != MPI_PROC_NULL && coordX == 0) {
        return true;
    }

    if (my_right_ != MPI_PROC_NULL && coordX == my_Nx_ - 1) {
        return true;
    }

    if (my_front_ != MPI_PROC_NULL && coordY == 0) {
        return true;
    }

    if (my_back_ != MPI_PROC_NULL && coordY == my_Ny_ - 1) {
        return true;
    }

    if (my_up_ != MPI_PROC_NULL && coordZ == my_Nz_ - 1) {
        return true;
    }

    if (my_down_ != MPI_PROC_NULL && coordZ == 0) {
        return true;
    }

    return false;
}

bool SubDomain::IsFluidInSubDomain(int coord_x, int coord_y, int coord_z) {
    std::vector<int> L2G = Local2Global(coord_x, coord_y, coord_z);

    if (Geometry_.pointsBox[L2G[0]][L2G[1]][L2G[2]]!=nullptr) {

        return true;
    }
    else {
        return false;
    }
}



void SubDomain::CreateNodes() {
    std::vector<Point> InletBoundary = Geometry_.InletboundaryPoints;
    std::vector<Point> OutletBoundary = Geometry_.OutletboundaryPoints;
    std::vector<Point> WallBoundary = Geometry_.WallboundaryPoints;

    m_lattice.resize(my_Nx_);

    for (int i = 0; i < my_Nx_; i++) {
        m_lattice[i].resize(my_Ny_);
        for (int j = 0; j < my_Ny_; j++) {
            m_lattice[i][j].resize(my_Nz_);
            for (int k = 0; k < my_Nz_; k++) {
                m_lattice[i][j][k] = nullptr;
            }
        }
    }

    for (int i = 0; i < my_Nx_; i++) {
        for (int j = 0; j < my_Ny_; j++) {
            for (int k = 0; k < my_Nz_; k++) {
                if (IsHaloLayer(i, j, k)) {
                    m_lattice[i][j][k] = std::make_shared<Node>(i, j, k, NodeType::Halo);
                }
                if (IsFluidInSubDomain(i, j, k)) {
                    std::vector<int> L2GIndex = Local2Global(i, j, k);
                    int global_coordX = L2GIndex[0];
                    int global_coordY = L2GIndex[1];
                    int global_coordZ = L2GIndex[2];

                    std::shared_ptr<Point> GlobalcurrentPointPtr = Geometry_.pointsBox[global_coordX][global_coordY][global_coordZ];
                    Point& GlobalcurrentPoint = *GlobalcurrentPointPtr; 

                    if (Utility_.isPointInBoundary(GlobalcurrentPoint, InletBoundary)) {
                        m_lattice[i][j][k] = std::make_shared<ZouHe>(i, j, k, m_velSet, NodeType::Boundary);
                        m_lattice[i][j][k]->BoundaryType_ = BoundaryType::Inlet;
                        InletBoundaryNode_.push_back(m_lattice[i][j][k]);
                    }
                    else if (Utility_.isPointInBoundary(GlobalcurrentPoint, OutletBoundary)) {
                        m_lattice[i][j][k] = std::make_shared<Node>(i, j, k, NodeType::Boundary);
                        m_lattice[i][j][k]->BoundaryType_ = BoundaryType::Outlet;
                        OutletBoundaryNode_.push_back(m_lattice[i][j][k]);
                    }
                    else if (Utility_.isPointInBoundary(GlobalcurrentPoint, WallBoundary)) {
                        m_lattice[i][j][k] = std::make_shared<Node>(i, j, k, NodeType::Boundary);
                        m_lattice[i][j][k]->BoundaryType_ = BoundaryType::Wall;
                        WallBoundaryNode_.push_back(m_lattice[i][j][k]);
                    }
                    else {
                        m_lattice[i][j][k] = std::make_shared<Node>(i, j, k, NodeType::Fluid);
                    }
                }
            }
        }
    }
}



void SubDomain::Connectivity() {

    for (int i = 0; i < my_Nx_; i++) {
        for (int j = 0; j < my_Ny_; j++) {
            for (int k = 0; k < my_Nz_; k++) {
                if (m_lattice[i][j][k] != nullptr) {


                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        int x_neighbor = i + m_velSet->GetDirection(dir)[0];
                        int y_neighbor = j + m_velSet->GetDirection(dir)[1];
                        int z_neighbor = k + m_velSet->GetDirection(dir)[2];

                        if (IsInSubDomain(x_neighbor, y_neighbor, z_neighbor) &&
                            m_lattice[x_neighbor][y_neighbor][z_neighbor] != nullptr) {
                            m_lattice[i][j][k]->m_neighbours[dir] = m_lattice[x_neighbor][y_neighbor][z_neighbor];
                        }
                        else {
                            m_lattice[i][j][k]->m_neighbours[dir] = nullptr;
                        }
                    }
                }
            }
        }
    }
}



void SubDomain::OutputConnectivityToVTK(const std::string& filename) {
    std::ofstream vtkFile(filename);
    vtkFile << "# vtk DataFile Version 3.0" << std::endl;
    vtkFile << "Lattice Connectivity" << std::endl;
    vtkFile << "ASCII" << std::endl;
    vtkFile << "DATASET POLYDATA" << std::endl;

    vtkFile << "POINTS " << my_Nx_ * my_Ny_ * my_Nz_ << " float" << std::endl;
    for (int i = 0; i < my_Nx_; i++) {
        for (int j = 0; j < my_Ny_; j++) {
            for (int k = 0; k < my_Nz_; k++) {
                vtkFile << i << " " << j << " " << k << std::endl;
            }
        }
    }

    int totalLines = 0;
    for (int i = 0; i < my_Nx_; i++) {
        for (int j = 0; j < my_Ny_; j++) {
            for (int k = 0; k < my_Nz_; k++) {
                if (m_lattice[i][j][k] != nullptr) {
                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        if (m_lattice[i][j][k]->m_neighbours[dir] != nullptr) {
                            totalLines++;
                        }
                    }
                }
            }
        }
    }
    vtkFile << "LINES " << totalLines << " " << totalLines * 3 << std::endl;
    for (int i = 0; i < my_Nx_; i++) {
        for (int j = 0; j < my_Ny_; j++) {
            for (int k = 0; k < my_Nz_; k++) {
                if (m_lattice[i][j][k] != nullptr) {
                    int baseIndex = (i * my_Ny_ + j) * my_Nz_ + k;
                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        if (m_lattice[i][j][k]->m_neighbours[dir] != nullptr) {
                            int nx = i + m_velSet->GetDirection(dir)[0];
                            int ny = j + m_velSet->GetDirection(dir)[1];
                            int nz = k + m_velSet->GetDirection(dir)[2];
                            int neighborIndex = (nx * my_Ny_ + ny) * my_Nz_ + nz;
                            vtkFile << "2 " << baseIndex << " " << neighborIndex << std::endl;
                        }
                    }
                }
            }
        }
    }
}


void SubDomain::Need_Comm() {

    need_communication_.resize(18, 0);

    bool fluidFound[18] = { false }; 

    for (int i = 0; i < my_Nx_; i++) {
        for (int j = 0; j < my_Ny_; j++) {
            for (int k = 0; k < my_Nz_; k++) {

                //left
                if (!fluidFound[0]  && m_lattice[1][j][k]!=nullptr) {
                    need_communication_[0] = 1;
                    fluidFound[0] = true;
                }

                //right
                if (!fluidFound[1] && m_lattice[my_Nx_ - 2][j][k]!=nullptr) {
                    need_communication_[1] = 1;
                    fluidFound[1] = true;
                }
                //front
                if (!fluidFound[2] && m_lattice[i][1][k]!=nullptr) {
                    need_communication_[2] = 1;
                    fluidFound[2] = true;
                }
                //back
                if (!fluidFound[3] && m_lattice[i][my_Ny_ - 2][k] !=nullptr) {
                    need_communication_[3] = 1;
                    fluidFound[3] = true;
                }
                //up
                if (!fluidFound[4] && m_lattice[i][j][my_Nz_ - 2]!=nullptr) {
                    need_communication_[4] = 1;
                    fluidFound[4] = true;
                }
                //down
                if (!fluidFound[5] && m_lattice[i][j][1]!=nullptr) {
                    need_communication_[5] = 1;
                    fluidFound[5] = true;
                }
                //right_up
                if (!fluidFound[6] && m_lattice[my_Nx_ - 2][j][my_Nz_ - 2]!=nullptr) {
                    need_communication_[6] = 1;
                    fluidFound[6] = true;
                }
                //right_down
                if (!fluidFound[7] && m_lattice[my_Nx_ - 2][j][1]!=nullptr) {
                    need_communication_[7] = 1;
                    fluidFound[7] = true;
                }
                //left_up
                if (!fluidFound[8] && m_lattice[1][j][my_Nz_ - 2]!=nullptr) {
                    need_communication_[8] = 1;
                    fluidFound[8] = true;
                }
                //left_down
                if (!fluidFound[9] && m_lattice[1][j][1]!=nullptr) {
                    need_communication_[9] = 1;
                    fluidFound[9] = true;
                }
                //up_front
                if (!fluidFound[10] && m_lattice[i][1][my_Nz_ - 2]!=nullptr) {
                    need_communication_[10] = 1;
                    fluidFound[10] = true;
                }
                //up_back
                if (!fluidFound[11] && m_lattice[i][my_Ny_ - 2][my_Nz_ - 2]!=nullptr) {
                    need_communication_[11] = 1;
                    fluidFound[11] = true;
                }
                //down_front
                if (!fluidFound[12] && m_lattice[i][1][1]!=nullptr) {
                    need_communication_[12] = 1;
                    fluidFound[12] = true;
                }
                //down_back
                if (!fluidFound[13] && m_lattice[i][my_Ny_ - 2][1]!=nullptr) {
                    need_communication_[13] = 1;
                    fluidFound[13] = true;
                }
                //right_front
                if (!fluidFound[14] && m_lattice[my_Nx_ - 2][1][k]!=nullptr) {
                    need_communication_[14] = 1;
                    fluidFound[14] = true;
                }
                //right_back
                if (!fluidFound[15] && m_lattice[my_Nx_ - 2][my_Ny_ - 2][k]!=nullptr) {
                    need_communication_[15] = 1;
                    fluidFound[15] = true;
                }
                //left_front
                if (!fluidFound[16] && m_lattice[1][1][k]!=nullptr) {
                    need_communication_[16] = 1;
                    fluidFound[16] = true;
                }
                //left_back
                if (!fluidFound[17] && m_lattice[1][my_Ny_ - 2][k]!=nullptr) {
                    need_communication_[17] = 1;
                    fluidFound[17] = true;
                }
            }
        }
    }
}





void SubDomain::Initialize()
{
    for (int i = 0; i < my_Nx_; i++) {
        for (int j = 0; j < my_Ny_; j++) {
            for (int k = 0; k < my_Nz_; k++) {

                if (m_lattice[i][j][k] != nullptr) {

                    for (int dir = 0; dir < _nLatNodes; dir++) {
                        m_lattice[i][j][k]->Initialize(m_velSet, dir);
                    }

                }
            }
        }
    }
}





const std::shared_ptr<VelocitySet> SubDomain::GetVelSet()
{
    return m_velSet;
}

void SubDomain::GetInvolevedLat() {

    for (int i = 0; i < my_Nx_; i++) {
        for (int j = 0; j < my_Ny_; j++) {
            for (int k = 0; k < my_Nz_; k++) {

                if (m_lattice[i][j][k] != nullptr) {

                        Involved_lat_.push_back(m_lattice[i][j][k]);

                }
            }
        }
    }
}



void SubDomain::WriteBoundaryLatticeToVTK(const std::string& filename, const std::string& boundaryType) {
    std::ofstream vtkFile(filename);

    if (!vtkFile.is_open()) {
        std::cerr << "Error: Unable to open VTK file for writing." << std::endl;
        return;
    }

    vtkFile << "# vtk DataFile Version 3.0" << std::endl;
    vtkFile << "Lattice Data" << std::endl;
    vtkFile << "ASCII" << std::endl;
    vtkFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

    std::vector<std::shared_ptr<Node>> boundaryNodes;

    if (boundaryType == "Inlet") {
        boundaryNodes = InletBoundaryNode_;
    }
    else if (boundaryType == "Outlet") {
        boundaryNodes = OutletBoundaryNode_;
    }
    else if (boundaryType == "Wall") {
        boundaryNodes = WallBoundaryNode_;
    }
    else {
        std::cerr << "Error: Unknown boundary type." << std::endl;
        return;
    }

    vtkFile << "POINTS " << boundaryNodes.size() << " double" << std::endl;
    for (const auto& nodePtr : boundaryNodes) {
        int xpos = nodePtr->x_position;
        int ypos = nodePtr->y_position;
        int zpos = nodePtr->z_position;

        std::vector<int> globalInd = Local2Global(xpos, ypos, zpos);
        vtkFile << globalInd[0] << " " << globalInd[1] << " " << globalInd[2] << std::endl;
    }

    vtkFile << "CELLS " << boundaryNodes.size() << " " << 2 * boundaryNodes.size() << std::endl;
    for (size_t i = 0; i < boundaryNodes.size(); ++i) {
        vtkFile << "1 " << i << std::endl;
    }

    vtkFile << "CELL_TYPES " << boundaryNodes.size() << std::endl;
    for (size_t i = 0; i < boundaryNodes.size(); ++i) {
        vtkFile << "1" << std::endl;
    }

    vtkFile.close();
    std::cout << "VTK file " << filename << " written successfully." << std::endl;
}



void SubDomain::WriteFluidLatticeToVTK() {
    std::string filename = "fluidLattice_" + std::to_string(myrank_) + ".vtk"; 

    std::ofstream vtkFile(filename);

    if (!vtkFile.is_open()) {
        std::cerr << "Error: Unable to open VTK file for writing." << std::endl;
        return;
    }

    int localFluidNodesCount = 0;
    for (int i = 0; i < my_Nx_; i++) {
        for (int j = 0; j < my_Ny_; j++) {
            for (int k = 0; k < my_Nz_; k++) {
                if (m_lattice[i][j][k] != nullptr
                    && m_lattice[i][j][k]->NodeType_ != NodeType::Halo) {
                    localFluidNodesCount++;
                }
            }
        }
    }

    int totalFluidNodesCount = 0;
    MPI_Allreduce(&localFluidNodesCount, &totalFluidNodesCount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    vtkFile << "# vtk DataFile Version 3.0" << std::endl;
    vtkFile << "Fluid Lattice Data for Processor " << myrank_ << std::endl;
    vtkFile << "ASCII" << std::endl;
    vtkFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

    std::vector<std::shared_ptr<Node>> fluidNodes;

    for (int i = 0; i < my_Nx_; i++) {
        for (int j = 0; j < my_Ny_; j++) {
            for (int k = 0; k < my_Nz_; k++) {
                if (m_lattice[i][j][k] != nullptr
                    && m_lattice[i][j][k]->NodeType_ != NodeType::Halo) {
                    fluidNodes.push_back(m_lattice[i][j][k]);
                }
            }
        }
    }

    vtkFile << "POINTS " << fluidNodes.size() << " double" << std::endl;
    for (const auto& nodePtr : fluidNodes) {

        int i = nodePtr->x_position;
        int j = nodePtr->y_position;
        int k = nodePtr->z_position;

        std::vector<int> globInd = Local2Global(i, j, k);


        vtkFile << globInd[0] << " " << globInd[1] << " " << globInd[2] << std::endl;
    }

    vtkFile << "CELLS " << fluidNodes.size() << " " << 2 * fluidNodes.size() << std::endl;
    for (size_t i = 0; i < fluidNodes.size(); ++i) {
        vtkFile << "1 " << i << std::endl;
    }

    vtkFile << "CELL_TYPES " << fluidNodes.size() << std::endl;
    for (size_t i = 0; i < fluidNodes.size(); ++i) {
        vtkFile << "1" << std::endl;
    }

    vtkFile.close();

    int globalRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);
    if (globalRank == 0) {
        std::cout << "Total fluid lattice nodes across all processors: " << totalFluidNodesCount << std::endl;
    }
}

