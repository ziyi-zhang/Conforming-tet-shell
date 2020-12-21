#include <shell/TetMeshCheck.h>
#include <shell/common.hpp>
#include <shell/Shell.h>
#include <shell/Utils.h>
#include <tetwild/CGALTypes.h>
#include <tetwild/Logger.h>
#include <igl/boundary_facets.h>

#include <Eigen/Dense>
#include <unordered_set>


namespace tetshell {

using namespace std;
using tetwild::logger;
using tetwild::TetVertex;

namespace {

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////

bool TetMeshCheck::SanityCheck(int eulerNumber) {

    logger().info("======================================");
    logger().info("Tetrahedral mesh sanity check results:");
    bool result = true;

    if (eulerNumber != std::numeric_limits<int>::max()) {
        if (eulerNumber == 1) {
            logger().debug(">>> Euler number == 1 >>>");
            result = result && true;
        } else {
            logger().warn("Euler number = {}", eulerNumber);
            result = result && false;
        }
    }

    if (args.positiveTet) {
        bool res = PositiveTetCheck();
        result = result && res;
    }

    if (args.conform) {
        bool res = ConformCheck();
        result = result && res;
    }

    if (args.vertexAttri) {
        bool res = VertexAttriCheck();
        result = result && res;
    }

    if (args.boundary) {
        bool res = BoundaryCheck();
        result = result && res;
    }

    if (args.euler) {
        bool res = EulerCheck();
        result = result && res;
    }

    logger().info("Tetrahedral mesh sanity check done.");
    logger().info("======================================");
    return result;
}


bool TetMeshCheck::ConformityCheck() {

    logger().info("======================================");
    logger().debug(">>> Post-refinement Conformity Check >>>");
    bool result = true;

    ConvertDouble();

    int idx1, idx2, idx3;
    int count = 0;
    for (int i=0; i<TO.size(); i++)
        for (int j=0; j<4; j++) {

            idx1 = TO[i][(j+1) % 4];
            idx2 = TO[i][(j+2) % 4];
            idx3 = TO[i][(j+3) % 4];
            if (VO[idx1].is_frozen && VO[idx2].is_frozen && VO[idx3].is_frozen) {
                std::unordered_set<int> sharedInputFace;
                UnorderedsetIntersection(VO[idx1].on_face, VO[idx2].on_face, VO[idx3].on_face, sharedInputFace);

                if (sharedInputFace.empty()) {
                    // if (VO[idx1].on_face.size() > 0) countNeg++;  // DEBUG PURPOSE
                    continue;  // this is the bounding box, skip
                }
                for (auto it=sharedInputFace.begin(); it!=sharedInputFace.end(); it++) {

                    int FI_idx = *it;
                    Eigen::Matrix<double, 3, 3> inputTriangle, outputTriangle;
                    inputTriangle << VI.row(FI(FI_idx, 0)), VI.row(FI(FI_idx, 1)), VI.row(FI(FI_idx, 2));
                    outputTriangle << VO[idx1].posf.x(), VO[idx1].posf.y(), VO[idx1].posf.z(), 
                                      VO[idx2].posf.x(), VO[idx2].posf.y(), VO[idx2].posf.z(), 
                                      VO[idx3].posf.x(), VO[idx3].posf.y(), VO[idx3].posf.z();
                    if (SameUnorderedTriangle(inputTriangle, outputTriangle, 0.0)) {
                        count++;
                    } else {
                        result = false;
                        std::cout << "Input triangle:" << std::endl;
                        std::cout << inputTriangle << std::endl;
                        std::cout << "Output triangle:" << std::endl;
                        std::cout << outputTriangle << std::endl;
                        logger().warn("Triangle coordinates not match: FI_idx={}", FI_idx);
                    }
                }
            }
        }

    if (count != 2 * (FI.rows()/4)) {
        result = false;
        logger().warn("#[locked faces] after refinement is wrong: count={} FI/4*2={}", count, FI.rows()/2);
    }

    logger().info("======================================");
    return result;
}


bool TetMeshCheck::PositiveTetCheck() {

    logger().debug(">>> Positive Tet Check >>>");
    bool result = true;

    for (auto it=TO.begin(); it!=TO.end(); it++) {
        if (!IsTetPositive(VO[it->at(0)].pos, VO[it->at(1)].pos, VO[it->at(2)].pos, VO[it->at(3)].pos)) {
            logger().warn("Flipped or degenerated tetrahedron. #TO = {}", it-TO.begin());
            result = false;
        }
    }

    return result;
}


bool TetMeshCheck::ConformCheck() {

    logger().debug(">>> Conform Check >>>");
    bool result = true;

    // Whether input faces are in tet mesh 
    // (not rigorous: the input face might be subdivided into many smaller triangles)
    int N = FI.rows();
    std::vector<int> on_input_face(N, 0);
    std::unordered_set<int> sharedInputFace;
    for (int i=0; i<TO.size(); i++) {
        // for the i-th tet
        for (int j=0; j<4; j++) {
            // for the j-th face
            UnorderedsetIntersection(VO[TO[i][j]].on_face, VO[TO[i][(j+1)%4]].on_face, VO[TO[i][(j+2)%4]].on_face, sharedInputFace);

            if (sharedInputFace.empty()) continue;
            if (sharedInputFace.size() > 1) tetwild::log_and_throw("ConformCheck: shared_input_face size > 1");
            int sharedInputFaceIdx = *(sharedInputFace.begin());
            on_input_face[sharedInputFaceIdx]++;
        }
    }
    for (int i=0; i<N; i++) {
        if (on_input_face[i] == 0) {
            logger().warn("Input triangle face #{} not found in tetMesh: cnt == 0", i);
            result = false;
        }
        // DEBUG PURPOSE
        // if (on_input_face[i] > 4)
        //     std::cerr << on_input_face[i] << std::endl;
    }

    return result;
}


bool TetMeshCheck::VertexAttriCheck() {

    logger().debug(">>> Vertex Attributes Check >>>");
    bool result = true;

    // make sure VO.posf has been computed
    ConvertDouble();

    // check tets are indeed in "conn_tets"
    std::vector<int> count_conn_tets(VO.size(), 0);
    for (int i=0; i<TO.size(); i++) {
        // for the i-th tet
        for (int j=0; j<4; j++) {
            // for the j-th vertex
            const TetVertex &vertex = VO[TO[i][j]];
            if (vertex.conn_tets.count(i) == 0) {
                logger().warn("conn_tets attribute of vertex {} does not contain tet {}", TO[i][j], i);
                result = false;
            } else {
                count_conn_tets[TO[i][j]]++;
            }
        }
    }
    // check "conn_tets" do not have deprecated tets
    for (int i=0; i<VO.size(); i++) {
        if (VO[i].conn_tets.size() != count_conn_tets[i]) {
            logger().warn("conn_tets attribute of vertex {} is wrong. conn_tet size = {}, count = {}", i, VO[i].conn_tets.size(), count_conn_tets[i]);
            result = false;
        }
    }

    // check "on_face": which input trianglular mesh faces are this vertex on
    for (int i=0; i<VO.size(); i++) {
        // for the i-th output vertex
        for (auto it=VO[i].on_face.begin(); it!=VO[i].on_face.end(); it++) {
            // for it-th input face
            int faceIdx = *it;
            // this output vertex must match one of the following three points
            bool matched = false;
            for (int j=0; j<3; j++) {
                if (VI(FI(faceIdx, j), 0) == VO[i].posf[0] && VI(FI(faceIdx, j), 1) == VO[i].posf[1] && VI(FI(faceIdx, j), 2) == VO[i].posf[2]) {
                    matched = true;
                    break;
                }
            }
            if (!matched) {
                logger().warn("on_face attribute of vertex {} contains wrong face index {}", i, faceIdx);
                result = false;
            }
        }
    }

    // what about on_edge and on_face?

    return result;
}


bool TetMeshCheck::BoundaryCheck() {

    logger().debug(">>> Boundary Check >>>");
    bool result = true;

    Eigen::MatrixXi T_temp, F_boundary;
    Eigen::VectorXi T_temp2TO, F_index, tet_local_idx;
    int cnt_surface_inner = 0, cnt_surface_bottom = 0, cnt_surface_top = 0, cnt_surface_outer = 0;
    int cnt_face_on_shell_inner = 0, cnt_face_on_shell_bottom = 0, cnt_face_on_shell_top = 0, cnt_face_on_shell_outer = 0;
    int cnt = 0;

    for (int i=0; i<face_on_shell.size(); i++)
        for (int j=0; j<4; j++) {

            switch (face_on_shell[i][j]) {
            case SURFACE_INNER:
                cnt_face_on_shell_inner++;
                break;
            case SURFACE_BOTTOM:
                cnt_face_on_shell_bottom++;
                break;
            case SURFACE_TOP:
                cnt_face_on_shell_top++;
                break;
            case SURFACE_OUTER:
                cnt_face_on_shell_outer++;
                break;
            default:
                break;
            }
        }
    if ((cnt_face_on_shell_top != FI.rows()/4) || (cnt_face_on_shell_bottom != FI.rows()/4)) {
        logger().warn("face_on_shell error: cnt_face_on_shell_top={}, cnt_face_on_shell_bottom={}, FI.rows()/4={}", cnt_face_on_shell_top, cnt_face_on_shell_bottom, FI.rows()/4);
        result = false;
    }

    //////////////////////////////////////////////////
    // Check region-1 is bounded by surface 1 and 2 //
    //////////////////////////////////////////////////
    // retrieve region-1 tets in F_temp
    for (int i=0; i<labels.rows(); i++) {
        if (labels(i) == SHELL_INNER_BOTTOM) cnt++;
    }
    T_temp.resize(cnt, 4);
    T_temp2TO.resize(cnt, 1);  // index in T_temp -> index in TO
    cnt = 0;
    for (int i=0; i<labels.rows(); i++) {
        if (labels(i) == SHELL_INNER_BOTTOM) {
            T_temp(cnt, 0) = TO[i][0];
            T_temp(cnt, 1) = TO[i][1];
            T_temp(cnt, 2) = TO[i][2];
            T_temp(cnt, 3) = TO[i][3];
            T_temp2TO(cnt) = i;
            cnt++;
        }
    }
    // find boundary
    igl::boundary_facets(T_temp, F_boundary, F_index, tet_local_idx);
    // check 
    for (int i=0; i<F_index.rows(); i++) {

        int TO_boundary_idx = T_temp2TO[F_index[i]];
        if (face_on_shell[TO_boundary_idx][tet_local_idx(i)] == SURFACE_INNER) {
            cnt_surface_inner++;
        } else if (face_on_shell[TO_boundary_idx][tet_local_idx(i)] == SURFACE_BOTTOM) {
            cnt_surface_bottom++;
        } else {
            logger().warn("Boundary of SHELL_INNER_BOTTOM not valid. TO index = {}, face_on_shell = {}", TO_boundary_idx, face_on_shell[TO_boundary_idx][tet_local_idx(i)]);
            result = false;
        }
    }
    // check the number matches
    /// NOTE: Why "cnt_surface_inner * 2" multiplied by 2? Because each face is counted twice in face_on_shell (two sides)
    ///       Why "cnt_surface_bottom" not multiplied? Because the middle region is hallow
    if (cnt_surface_inner * 2 != cnt_face_on_shell_inner) {
        logger().warn("Boundary of SHELL_INNER_BOTTOM: cnts of SURFACE_INNER do not match. cnt_surface_inner={}, cnt_face_on_shell_inner={}", cnt_surface_inner, cnt_face_on_shell_inner);
    }
    if (cnt_surface_bottom != cnt_face_on_shell_bottom) {
        logger().warn("Boundary of SHELL_INNER_BOTTOM: cnts of SURFACE_BOTTOM do not match. cnt_surface_bottom={}, cnt_face_on_shell_bottom={}", cnt_surface_bottom, cnt_face_on_shell_bottom);
    }

    //////////////////////////////////////////////////
    // Check region-3 is bounded by surface 3 and 4 //
    //////////////////////////////////////////////////
    // retrieve region-3 tets in F_temp
    cnt = 0;
    for (int i=0; i<labels.rows(); i++) {
        if (labels(i) == SHELL_TOP_OUTER) cnt++;
    }
    T_temp.resize(cnt, 4);
    T_temp2TO.resize(cnt, 1);  // index in T_temp -> index in TO
    cnt = 0;
    for (int i=0; i<labels.rows(); i++) {
        if (labels(i) == SHELL_TOP_OUTER) {
            T_temp(cnt, 0) = TO[i][0];
            T_temp(cnt, 1) = TO[i][1];
            T_temp(cnt, 2) = TO[i][2];
            T_temp(cnt, 3) = TO[i][3];
            T_temp2TO(cnt) = i;
            cnt++;
        }
    }
    // find boundary
    igl::boundary_facets(T_temp, F_boundary, F_index, tet_local_idx);
    // check 
    for (int i=0; i<F_index.rows(); i++) {

        int TO_boundary_idx = T_temp2TO[F_index[i]];
        if (face_on_shell[TO_boundary_idx][tet_local_idx(i)] == SURFACE_TOP) {
            cnt_surface_top++;
        } else if (face_on_shell[TO_boundary_idx][tet_local_idx(i)] == SURFACE_OUTER) {
            cnt_surface_outer++;
        } else {
            logger().warn("Boundary of SHELL_TOP_OUTER not valid. TO index = {}, face_on_shell = {}", TO_boundary_idx, face_on_shell[TO_boundary_idx][tet_local_idx(i)]);
            result = false;
        }
    }
    // check the number matches
    if (cnt_surface_top != cnt_face_on_shell_top) {
        logger().warn("Boundary of SHELL_TOP_OUTER: cnts of SURFACE_TOP do not match. cnt_surface_top={}, cnt_face_on_shell_top={}", cnt_surface_top, cnt_face_on_shell_top);
    }
    if (cnt_surface_outer * 2 != cnt_face_on_shell_outer) {
        logger().warn("Boundary of SHELL_TOP_OUTER: cnts of SURFACE_OUTER do not match. cnt_surface_outer={}, cnt_face_on_shell_outer={}", cnt_surface_outer, cnt_face_on_shell_outer);
    }

    return result;
}


bool TetMeshCheck::EulerCheck() {

    return true;
}


void TetMeshCheck::ConvertDouble() {

    if (VO.size() == 0 || VO[0].is_rounded) return;

    for (auto it=VO.begin(); it!=VO.end(); it++) {
        it->is_rounded = true;
        it->round();
    }
}

}  // namespace tetshell
