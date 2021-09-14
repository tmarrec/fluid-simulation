#include "MarchingCube.h"

#define  TINYPLY_IMPLEMENTATION
#include "tinyply.h"

inline double MarchingCube::interp(const Field<double,std::uint16_t>& F, double x, double y, double z) const
{
    const std::uint16_t i0 = static_cast<std::uint16_t>(x);
    const std::uint16_t i1 = i0 + 1;
    const std::uint16_t j0 = static_cast<std::uint16_t>(y);
    const std::uint16_t j1 = j0 + 1;
    const std::uint16_t k0 = static_cast<std::uint16_t>(z);
    const std::uint16_t k1 = k0 + 1;

    const double s1 = x - i0;
    const double s0 = 1.0 - s1;
    const double t1 = y - j0;
    const double t0 = 1.0 - t1;
    const double u1 = z - k0;
    const double u0 = 1.0 - u1;

    return s0 * (     t0 * ( u0 * F(i0,j0,k0) + u1 * F(i0,j0,k1) )
                    + t1 * ( u0 * F(i0,j1,k0) + u1 * F(i0,j1,k1) )
                )
         + s1 * (     t0 * ( u0 * F(i1,j0,k0) + u1 * F(i1,j0,k1) )
                    + t1 * ( u0 * F(i1,j1,k0) + u1 * F(i1,j1,k1) )
                );
}

inline std::array<glm::vec3, 3> MarchingCube::computeNormal(std::array<glm::vec3, 3> p, float eps) const
{
}

inline bool MarchingCube::check(const glm::vec3 &left, const glm::vec3 &right) const
{
    if (left.x < right.x)
        return true;
    else if (left.x > right.x)
        return false;
    if (left.y < right.y)
        return true;
    else if (left.y > right.y)
        return false;
    if (left.z < right.z)
        return true;
    else if (left.z > right.z)
        return false;
    return false;
}

glm::vec3 MarchingCube::VertexInterp(glm::vec3 p1, glm::vec3 p2, float valp1, float valp2) const
{
    if (check(p2, p1))
    {
        glm::vec3 temp = p1;
        p1 = p2;
        p2 = temp;
        float temp2 = valp1;
        valp1 = valp2;
        valp2 = temp2;
    }
    glm::vec3 p;
    if(std::fabs(valp1 - valp2) > 0.00001f)
    {
        p = p1 + (p2 - p1)/(valp2 - valp1)*(0.0f - valp1);
    }
    else
    {
        p = p1;
    }
    return p;
}

void MarchingCube::run(const Field<double, std::uint16_t>& F, std::uint64_t iteration)
{
    std::vector<float> vertices;
    std::vector<float> normals;
    std::vector<std::uint32_t> indices;
    const std::uint64_t nbEchant = 32;
    for (std::uint16_t k = 0; k < nbEchant; ++k)
    {
        for (std::uint16_t j = 0; j < nbEchant; ++j)
        {
            for (std::uint16_t i = 0; i < nbEchant; ++i)
            {
                const double x = (static_cast<double>(i)/static_cast<double>(nbEchant-1))*(F.x()-1);
                const double y = (static_cast<double>(j)/static_cast<double>(nbEchant-1))*(F.y()-1);
                const double z = (static_cast<double>(k)/static_cast<double>(nbEchant-1))*(F.z()-1);

                const double s = 0.5;
                const glm::vec3 p[8] =
                {
                    {x-s,y+s,z-s},
                    {x+s,y+s,z-s},
                    {x+s,y-s,z-s},
                    {x-s,y-s,z-s},
                    {x-s,y+s,z+s},
                    {x+s,y+s,z+s},
                    {x+s,y-s,z+s},
                    {x-s,y-s,z+s},
                };
                const double cell[8] =
                {
                    interp(F, p[0].x, p[0].y, p[0].z),
                    interp(F, p[1].x, p[1].y, p[1].z),
                    interp(F, p[2].x, p[2].y, p[2].z),
                    interp(F, p[3].x, p[3].y, p[3].z),
                    interp(F, p[4].x, p[4].y, p[4].z),
                    interp(F, p[5].x, p[5].y, p[5].z),
                    interp(F, p[6].x, p[6].y, p[6].z),
                    interp(F, p[7].x, p[7].y, p[7].z),
                };

                std::uint16_t cubeindex = 0;
                if (cell[0] < 0.0) cubeindex |= 1;
                if (cell[1] < 0.0) cubeindex |= 2;
                if (cell[2] < 0.0) cubeindex |= 4;
                if (cell[3] < 0.0) cubeindex |= 8;
                if (cell[4] < 0.0) cubeindex |= 16;
                if (cell[5] < 0.0) cubeindex |= 32;
                if (cell[6] < 0.0) cubeindex |= 64;
                if (cell[7] < 0.0) cubeindex |= 128;

                // Cube is entirely in/out of the surface
                if (_edgeTable[cubeindex] != 0)
                {
                    glm::vec3 vertlist[12];
                    // Find the vertices where the surface intersects the cube
                    if (_edgeTable[cubeindex] & 1)
                        vertlist[0] = VertexInterp(p[0],p[1],cell[0],cell[1]);
                    if (_edgeTable[cubeindex] & 2)
                        vertlist[1] = VertexInterp(p[1],p[2],cell[1],cell[2]);
                    if (_edgeTable[cubeindex] & 4)
                        vertlist[2] = VertexInterp(p[2],p[3],cell[2],cell[3]);
                    if (_edgeTable[cubeindex] & 8)
                        vertlist[3] = VertexInterp(p[3],p[0],cell[3],cell[0]);
                    if (_edgeTable[cubeindex] & 16)
                        vertlist[4] = VertexInterp(p[4],p[5],cell[4],cell[5]);
                    if (_edgeTable[cubeindex] & 32)
                        vertlist[5] = VertexInterp(p[5],p[6],cell[5],cell[6]);
                    if (_edgeTable[cubeindex] & 64)
                        vertlist[6] = VertexInterp(p[6],p[7],cell[6],cell[7]);
                    if (_edgeTable[cubeindex] & 128)
                        vertlist[7] = VertexInterp(p[7],p[4],cell[7],cell[4]);
                    if (_edgeTable[cubeindex] & 256)
                        vertlist[8] = VertexInterp(p[0],p[4],cell[0],cell[4]);
                    if (_edgeTable[cubeindex] & 512)
                        vertlist[9] = VertexInterp(p[1],p[5],cell[1],cell[5]);
                    if (_edgeTable[cubeindex] & 1024)
                        vertlist[10] = VertexInterp(p[2],p[6],cell[2],cell[6]);
                    if (_edgeTable[cubeindex] & 2048)
                        vertlist[11] = VertexInterp(p[3],p[7],cell[3],cell[7]);

                    // Create the triangle
                    for (std::uint64_t n = 0; _triTable[cubeindex][n] != -1; n += 3)
                    {
                        std::array<glm::vec3, 3> triangle;
                        triangle[0] = vertlist[static_cast<std::uint8_t>(_triTable[cubeindex][n  ])];
                        triangle[1] = vertlist[static_cast<std::uint8_t>(_triTable[cubeindex][n+1])];
                        triangle[2] = vertlist[static_cast<std::uint8_t>(_triTable[cubeindex][n+2])];

                        std::uint64_t oldSize = vertices.size();
                        vertices.resize(oldSize+9);
                        memcpy(&vertices[oldSize], &triangle, sizeof(float)*9);

                        std::array<glm::vec3, 3> ns;
                        ns[0] = glm::vec3{1,1,1};
                        ns[1] = glm::vec3{1,1,1};
                        ns[2] = glm::vec3{1,1,1};
                        normals.resize(oldSize+9);
                        memcpy(&normals[oldSize], &ns, sizeof(float)*9);
                    }
                }
            }
        }
    }
    indices.resize(vertices.size()/3);
	std::iota(indices.begin(), indices.end(), 0);

    std::filesystem::create_directories("result");
    std::filebuf fb_binary;
    std::string path = "result/";
    path += std::to_string(iteration);
    path += ".ply";
    fb_binary.open(path, std::ios::out | std::ios::trunc | std::ios::binary);
    std::ostream outstream_binary(&fb_binary);
    if (outstream_binary.fail())
    {
        ERROR("Failed to open res.ply file");
    }

    tinyply::PlyFile meshFile;
    meshFile.add_properties_to_element("vertex", {"x","y","z"}, tinyply::Type::FLOAT32, vertices.size()/3, reinterpret_cast<std::uint8_t*>(vertices.data()), tinyply::Type::INVALID, 0);
    meshFile.add_properties_to_element("vertex", {"nx","ny","nz"}, tinyply::Type::FLOAT32, normals.size()/3, reinterpret_cast<std::uint8_t*>(normals.data()), tinyply::Type::INVALID, 0);
    meshFile.add_properties_to_element("face", {"vertex_indices"}, tinyply::Type::UINT32, indices.size()/3, reinterpret_cast<std::uint8_t*>(indices.data()), tinyply::Type::UINT8, 3);
    meshFile.get_comments().push_back("Generated by tmarrec.dev fluid-simulation project.");

    meshFile.write(outstream_binary, true);
    fb_binary.close();
}
