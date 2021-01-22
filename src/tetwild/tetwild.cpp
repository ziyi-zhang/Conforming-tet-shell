// This file is part of TetWild, a software for generating tetrahedral meshes.
//
// Copyright (C) 2018 Yixin Hu <yixin.hu@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
//
// Created by Yixin Hu on 5/31/18.
//

#include <tetwild/tetwild.h>
#include <tetwild/Common.h>
#include <tetwild/Logger.h>
#include <tetwild/Preprocess.h>
#include <tetwild/DelaunayTetrahedralization.h>
#include <tetwild/BSPSubdivision.h>
#include <tetwild/SimpleTetrahedralization.h>
#include <tetwild/MeshRefinement.h>
#include <tetwild/InoutFiltering.h>
#include <shell/common.hpp>
#include <shell/Label.h>
#include <shell/Shell.h>
#include <shell/Utils.h>
#include <shell/TetMeshCheck.h>
#include <igl/boundary_facets.h>
#include <igl/remove_unreferenced.h>
#include <pymesh/MshSaver.h>
#include <geogram/mesh/mesh.h>


namespace tetwild {

////////////////////////////////////////////////////////////////////////////////

void printFinalQuality(double time, const std::vector<TetVertex>& tet_vertices,
                       // const std::vector<std::array<int, 4>>& tets,
                       const std::vector<bool> &t_is_removed,
                       const std::vector<TetQuality>& tet_qualities,
                       const std::vector<int>& v_ids,
                       const Args &args, const State &state)
{
    logger().debug("final quality:");
    double min = 10, max = 0;
    double min_avg = 0, max_avg = 0;
    // double max_asp_ratio = 0, avg_asp_ratio = 0;
    double max_slim_energy = 0, avg_slim_energy = 0;
    std::array<double, 6> cmp_cnt = {{0, 0, 0, 0, 0, 0}};
    std::array<double, 6> cmp_d_angles = {{6 / 180.0 * M_PI, 12 / 180.0 * M_PI, 18 / 180.0 * M_PI,
                                           162 / 180.0 * M_PI, 168 / 180.0 * M_PI, 174 / 180.0 * M_PI}};
    int cnt = 0;
    for (int i = 0; i < tet_qualities.size(); i++) {
        if (t_is_removed[i])
            continue;
        cnt++;
        if (tet_qualities[i].min_d_angle < min)
            min = tet_qualities[i].min_d_angle;
        if (tet_qualities[i].max_d_angle > max)
            max = tet_qualities[i].max_d_angle;
        // if (tet_qualities[i].asp_ratio_2 > max_asp_ratio)
            // max_asp_ratio = tet_qualities[i].asp_ratio_2;
        if (tet_qualities[i].slim_energy > max_slim_energy)
            max_slim_energy = tet_qualities[i].slim_energy;
        min_avg += tet_qualities[i].min_d_angle;
        max_avg += tet_qualities[i].max_d_angle;
        // avg_asp_ratio += tet_qualities[i].asp_ratio_2;
        avg_slim_energy += tet_qualities[i].slim_energy;

        for (int j = 0; j < 3; j++) {
            if (tet_qualities[i].min_d_angle < cmp_d_angles[j])
                cmp_cnt[j]++;
        }
        for (int j = 0; j < 3; j++) {
            if (tet_qualities[i].max_d_angle > cmp_d_angles[j + 3])
                cmp_cnt[j + 3]++;
        }
    }
    logger().debug("min_d_angle = {}, max_d_angle = {}, max_slim_energy = {}", min, max, max_slim_energy);
    logger().debug("avg_min_d_angle = {}, avg_max_d_angle = {}, avg_slim_energy = {}", min_avg / cnt, max_avg / cnt, avg_slim_energy / cnt);
    logger().debug("min_d_angle: <6 {};   <12 {};  <18 {}", cmp_cnt[0] / cnt, cmp_cnt[1] / cnt, cmp_cnt[2] / cnt);
    logger().debug("max_d_angle: >174 {}; >168 {}; >162 {}", cmp_cnt[5] / cnt, cmp_cnt[4] / cnt, cmp_cnt[3] / cnt);

    addRecord(MeshRecord(MeshRecord::OpType::OP_WN, time, v_ids.size(), cnt,
                         min, min_avg / cnt, max, max_avg / cnt, max_slim_energy, avg_slim_energy / cnt), args, state);

    // output unrounded vertices:
    cnt = 0;
    for (int v_id: v_ids) {
        if (!tet_vertices[v_id].is_rounded) {
            cnt++;
        }
    }
    logger().debug("{}/{} vertices are unrounded!!!", cnt, v_ids.size());
    addRecord(MeshRecord(MeshRecord::OpType::OP_UNROUNDED, -1, cnt, -1), args, state);
}


void extractSurfaceMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T,
    Eigen::MatrixXd &VS, Eigen::MatrixXi &FS)
{
    Eigen::VectorXi I;
    igl::boundary_facets(T, FS);
    igl::remove_unreferenced(V, FS, VS, FS, I);
    for(int i=0;i < FS.rows();i++){
        int tmp = FS(i, 0);
        FS(i, 0) = FS(i, 2);
        FS(i, 2) = tmp;
    }
}


void extractFinalTetmesh(MeshRefinement& MR,
    Eigen::MatrixXd &V_out, Eigen::MatrixXi &T_out, Eigen::VectorXd &A_out,
    const Args &args, const State &state) {

    std::vector<TetVertex> &tet_vertices = MR.tet_vertices;
    std::vector<std::array<int, 4>> &tets = MR.tets;
//    std::vector<bool> &v_is_removed = MR.v_is_removed;
    std::vector<bool> &t_is_removed = MR.t_is_removed;
    std::vector<TetQuality> &tet_qualities = MR.tet_qualities;
    int t_cnt = std::count(t_is_removed.begin(), t_is_removed.end(), false);
    double tmp_time = 0;

/*
    if (!args.smooth_open_boundary) {

        InoutFiltering IOF(tet_vertices, tets, MR.is_surface_fs, v_is_removed, t_is_removed, tet_qualities, state);

        igl::Timer igl_timer;
        igl_timer.start();
        IOF.filter();  // do in-out filter
        t_cnt = std::count(t_is_removed.begin(), t_is_removed.end(), false);
        tmp_time = igl_timer.getElapsedTime();
        logger().info("time = {}s", tmp_time);
        logger().debug("{} tets inside!", t_cnt);
    }
*/

    // v_ids is the vector of index of vertices
    std::vector<int> v_ids;
    for (int i = 0; i < tets.size(); i++) {
        if (t_is_removed[i])
            continue;
        for (int j = 0; j < 4; j++)
            v_ids.push_back(tets[i][j]);
    }
    // unique "v_ids"
    std::sort(v_ids.begin(), v_ids.end());
    v_ids.erase(std::unique(v_ids.begin(), v_ids.end()), v_ids.end());
    // "map_ids" is defined as this for all used vertices:
    // Key: index in V
    // Value: new index i
    std::unordered_map<int, int> map_ids;
    for (int i = 0; i < v_ids.size(); i++)
        map_ids[v_ids[i]] = i;

    // Prepare V, T, A
    V_out.resize(v_ids.size(), 3);
    T_out.resize(t_cnt, 4);
    A_out.resize(t_cnt);
    // Fill V
    for (int i = 0; i < v_ids.size(); i++) {
        for (int j = 0; j < 3; j++) {
            V_out(i, j) = tet_vertices[v_ids[i]].posf[j];
        }
    }
    // Fill T & A
    int cnt = 0;
    for (int i = 0; i < tets.size(); i++) {
        if (t_is_removed[i]) {
            continue;
        }
        for (int j = 0; j < 4; j++) {
            T_out(cnt, j) = map_ids[tets[i][j]];  // the index is new as defined in map_ids
        }
        A_out(cnt) = tet_qualities[i].min_d_angle;
        cnt++;
    }
    logger().debug("final output #v = {}", V_out.rows());
    logger().debug("final output #t = {}", T_out.rows());

    if (args.is_quiet) {
        return;
    }
    printFinalQuality(tmp_time, tet_vertices, t_is_removed, tet_qualities, v_ids, args, state);
}


////////////////////////////////////////////////////////////////////////////////


// Simplify the input surface by swapping and removing edges, while staying within the envelope
double tetwild_stage_one_preprocess(
    const Eigen::MatrixXd &VI,
    const Eigen::MatrixXi &FI,
    const Args &args,
    State &state,
    GEO::Mesh &geo_sf_mesh,
    GEO::Mesh &geo_b_mesh,
    std::vector<Point_3> &m_vertices,
    std::vector<std::array<int, 3>> &m_faces)
{
    igl::Timer igl_timer;
    igl_timer.start();
    logger().info("Preprocessing...");
    Preprocess pp(state);
    if (!pp.init(VI, FI, geo_b_mesh, geo_sf_mesh, args)) {
        //todo: output a empty tetmesh
        PyMesh::MshSaver mSaver(state.working_dir + state.postfix + ".msh", true);
        Eigen::VectorXd oV;
        Eigen::VectorXi oT;
        oV.resize(0);
        oT.resize(0);
        mSaver.save_mesh(oV, oT, 3, mSaver.TET);
        log_and_throw("Empty mesh!");
    }
    addRecord(MeshRecord(MeshRecord::OpType::OP_INIT, 0, geo_sf_mesh.vertices.nb(), geo_sf_mesh.facets.nb()), args, state);

    m_vertices.clear();
    m_faces.clear();
    pp.process(geo_sf_mesh, m_vertices, m_faces, args);
    double tmp_time = igl_timer.getElapsedTime();
    addRecord(MeshRecord(MeshRecord::OpType::OP_PREPROCESSING, tmp_time, m_vertices.size(), m_faces.size()), args, state);
    logger().info("time = {}s", tmp_time);
    return tmp_time;
}

// -----------------------------------------------------------------------------

// Compute an initial Delaunay triangulation of the input triangle soup
double tetwild_stage_one_delaunay(
    const Args &args,
    const State &state,
    GEO::Mesh &geo_sf_mesh,
    const std::vector<Point_3> &m_vertices,
    const std::vector<std::array<int, 3>> &m_faces,
    std::vector<Point_3> &bsp_vertices,
    std::vector<BSPEdge> &bsp_edges,
    std::vector<BSPFace> &bsp_faces,
    std::vector<BSPtreeNode> &bsp_nodes,
    std::vector<int> &m_f_tags,
    std::vector<int> &raw_e_tags,
    std::vector<std::vector<int>> &raw_conn_e4v)
{
    igl::Timer igl_timer;
    igl_timer.start();
    logger().info("Delaunay tetrahedralizing...");
    DelaunayTetrahedralization DT;
    m_f_tags.clear();
    raw_e_tags.clear();
    raw_conn_e4v.clear();
    DT.init(m_vertices, m_faces, m_f_tags, raw_e_tags, raw_conn_e4v);
    bsp_vertices.clear();
    bsp_edges.clear();
    bsp_faces.clear();
    bsp_nodes.clear();
    DT.tetra(m_vertices, m_faces, geo_sf_mesh, bsp_vertices, bsp_edges, bsp_faces, bsp_nodes, args, state);
    logger().debug("# bsp_vertices = {}", bsp_vertices.size());
    logger().debug("# bsp_edges = {}", bsp_edges.size());
    logger().debug("# bsp_faces = {}", bsp_faces.size());
    logger().debug("# bsp_nodes = {}", bsp_nodes.size());
    logger().info("Delaunay tetrahedralization done!");
    double tmp_time = igl_timer.getElapsedTime();
    addRecord(MeshRecord(MeshRecord::OpType::OP_DELAUNEY_TETRA, tmp_time, bsp_vertices.size(), bsp_nodes.size()), args, state);
    logger().info("time = {}s", tmp_time);
    return tmp_time;
}

// -----------------------------------------------------------------------------

// Match faces of the Delaunay tetrahedralization with faces from the input mesh
double tetwild_stage_one_mc(
    const Args &args,
    const State &state,
    MeshConformer &MC)
{
    igl::Timer igl_timer;
    igl_timer.start();
    logger().info("Divfaces matching...");
    MC.match();
    logger().info("Divfaces matching done!");
    double tmp_time = igl_timer.getElapsedTime();
    addRecord(MeshRecord(MeshRecord::OpType::OP_DIVFACE_MATCH, tmp_time, MC.bsp_vertices.size(), MC.bsp_nodes.size()), args, state);
    logger().info("time = {}s", tmp_time);
    return tmp_time;
}

// -----------------------------------------------------------------------------

// Compute BSP partition of the domain
double tetwild_stage_one_bsp(
    const Args &args,
    const State &state,
    MeshConformer &MC)
{
    igl::Timer igl_timer;
    igl_timer.start();
    logger().info("BSP subdivision ...");
    BSPSubdivision BS(MC);
    BS.init();
    BS.subdivideBSPNodes();
    logger().debug("Output: ");
    logger().debug("# node = {}", MC.bsp_nodes.size());
    logger().debug("# face = {}", MC.bsp_faces.size());
    logger().debug("# edge = {}", MC.bsp_edges.size());
    logger().debug("# vertex = {}", MC.bsp_vertices.size());
    logger().info("BSP subdivision done!");
    double tmp_time = igl_timer.getElapsedTime();
    addRecord(MeshRecord(MeshRecord::OpType::OP_BSP, tmp_time, MC.bsp_vertices.size(), MC.bsp_nodes.size()), args, state);
    logger().info("time = {}s", tmp_time);
    return tmp_time;
}

// -----------------------------------------------------------------------------

// Compute an initial tetrahedral mesh from the BSP partition
double tetwild_stage_one_tetra(
    const Args &args,
    const State &state,
    MeshConformer &MC,
    const std::vector<int> &m_f_tags,
    const std::vector<int> &raw_e_tags,
    const std::vector<std::vector<int>> &raw_conn_e4v,
    std::vector<TetVertex> &tet_vertices,
    std::vector<std::array<int, 4>> &tet_indices,
    std::vector<std::array<int, 4>> &is_surface_facet, 
    std::vector<std::array<int, 4>>& face_on_shell,
    std::vector<std::array<int, 4>>& faceIdx_on_shell) {

    igl::Timer igl_timer;
    igl_timer.start();
    logger().info("Tetrehedralizing ...");

    SimpleTetrahedralization ST(state, MC);
    tet_vertices.clear();
    tet_indices.clear();
    is_surface_facet.clear();
    face_on_shell.clear();
    faceIdx_on_shell.clear();

    ST.tetra(tet_vertices, tet_indices);
    ST.labelSurface(m_f_tags, raw_e_tags, raw_conn_e4v, tet_vertices, tet_indices, is_surface_facet, face_on_shell, faceIdx_on_shell);
    ST.labelBbox(tet_vertices, tet_indices);
    if (!state.is_mesh_closed)//if input is an open mesh
        ST.labelBoundary(tet_vertices, tet_indices, is_surface_facet);

    logger().debug("# tet_vertices = {}", tet_vertices.size());
    logger().debug("# tets = {}", tet_indices.size());
    logger().info("Tetrahedralization done!");
    double tmp_time = igl_timer.getElapsedTime();
    addRecord(MeshRecord(MeshRecord::OpType::OP_SIMPLE_TETRA, tmp_time, tet_vertices.size(), tet_indices.size()), args, state);
    logger().info("time = {}s", tmp_time);

    return tmp_time;
}

////////////////////////////////////////////////////////////////////////////////




///
/// Generate tet mesh from input
///
/// TODO
///
void tetwild_stage_one(
    const Eigen::MatrixXd &VI,
    const Eigen::MatrixXi &FI,
    const Args &args,
    State &state,
    GEO::Mesh &geo_sf_mesh,
    GEO::Mesh &geo_b_mesh,
    std::vector<TetVertex> &tet_vertices,
    std::vector<std::array<int, 4>> &tet_indices,
    std::vector<std::array<int, 4>> &is_surface_facet, 
    std::vector<std::array<int, 4>>& face_on_shell, 
    std::vector<std::array<int, 4>>& faceIdx_on_shell) {

    igl::Timer igl_timer;
    double sum_time = 0;

    // preprocess
    std::vector<Point_3> m_vertices;
    std::vector<std::array<int, 3>> m_faces;  // m_faces is F_in after deleting all degenerate triangles
    sum_time += tetwild_stage_one_preprocess(VI, FI, args, state, geo_sf_mesh, geo_b_mesh, m_vertices, m_faces);

    // delaunay tetrahedralization
    std::vector<Point_3> bsp_vertices;
    std::vector<BSPEdge> bsp_edges;
    std::vector<BSPFace> bsp_faces;
    std::vector<BSPtreeNode> bsp_nodes;
    std::vector<int> m_f_tags;
    std::vector<int> raw_e_tags;
    std::vector<std::vector<int>> raw_conn_e4v;
    sum_time += tetwild_stage_one_delaunay(args, state, geo_sf_mesh, m_vertices, m_faces,
        bsp_vertices, bsp_edges, bsp_faces, bsp_nodes, m_f_tags, raw_e_tags, raw_conn_e4v);

    // mesh conforming
    MeshConformer MC(m_vertices, m_faces, bsp_vertices, bsp_edges, bsp_faces, bsp_nodes);
    sum_time += tetwild_stage_one_mc(args, state, MC);

    // bsp subdivision
    sum_time += tetwild_stage_one_bsp(args, state, MC);

    // simple tetrahedralization
    sum_time += tetwild_stage_one_tetra(args, state, MC, m_f_tags, raw_e_tags, raw_conn_e4v,
        tet_vertices, tet_indices, is_surface_facet, face_on_shell, faceIdx_on_shell);

    logger().info("Total time for the first stage = {}s", sum_time);
}

// -----------------------------------------------------------------------------

///
/// Shell related operations
///
/// TODO
///
void tetwild_stage_shell(
    const Args &args,
    const Eigen::MatrixXd &VI,
    const Eigen::MatrixXi &FI,
    std::vector<TetVertex> &VO,
    std::vector<std::array<int, 4>> &TO,
    std::vector<std::array<int, 4>> &is_surface_facet, 
    std::vector<std::array<int, 4>> &face_on_shell,
    std::vector<std::array<int, 4>> &faceIdx_on_shell,
    Eigen::VectorXi &labels, 
    int &eulerNumber) {

    igl::Timer igl_timer;
    igl_timer.start();

    // the input mesh size
    logger().debug("VI size = {} x {}", VI.rows(), VI.cols());
    logger().debug("FI size = {} x {}", FI.rows(), FI.cols());
    logger().debug("VO size = {} TetVertex", VO.size());
    logger().debug("TO size = {} x 4", TO.size());

    // Report Euler number (this is immediately after tetwild stage 1)
    tetshell::EulerNumber(TO, "Pre-shell");

    // Convert VI to CGAL rational
    std::vector<Point_3> VI_cgal;
    for (int i=0; i<VI.rows(); i++) {
        VI_cgal.push_back(Point_3(VI(i, 0), VI(i, 1), VI(i, 2)));
    }

    // Construct dualShell & label tets
    tetshell::DualShell_t dualShell;
    tetshell::LabelTet(args, VI_cgal, FI, VO, TO, face_on_shell,  // input
                       dualShell, labels);  // output

    // Replace some tets with prism tets
    if (!args.skip_prism) {
        tetshell::ReplaceWithPrismTet(args, dualShell,  // input
                                      VO, TO, labels, is_surface_facet, face_on_shell, faceIdx_on_shell);  // output
    }

    // Freeze vertices on bottom and top surfaces
    tetshell::FreezeVertices(face_on_shell, TO,  // input
                             VO, is_surface_facet);  // output

    // Report Euler number again (with the hallow region filled with pseudo-tets)
    std::vector<std::array<int, 4>> TO_with_pseudo_tets;
    tetshell::GetMeshWithPseudoTets(dualShell, VO, TO,  // input
                                    TO_with_pseudo_tets);  // output
    eulerNumber = tetshell::EulerNumber(TO_with_pseudo_tets, "Post-shell with pseudo tets");

    double stageTime = igl_timer.getElapsedTime();
    logger().info("Total time for the shell stage = {}s", stageTime);
}

// -----------------------------------------------------------------------------

///
/// Mesh refinement
///
/// TODO
///
void tetwild_stage_two(const Args &args, State &state,
    GEO::Mesh &geo_sf_mesh,
    GEO::Mesh &geo_b_mesh,
    std::vector<TetVertex> &tet_vertices,
    std::vector<std::array<int, 4>> &tet_indices,
    std::vector<std::array<int, 4>> &is_surface_facet, 
    std::vector<bool> &t_is_removed, 
    std::vector<TetQuality> &tetQuality) {

    igl::Timer igl_timer;
    igl_timer.start();

    spdlog::level::level_enum verbose_level = logger().level();
    logger().set_level(static_cast<spdlog::level::level_enum>(spdlog::level::debug));
    // init
    logger().info("Refinement initializing...");
    MeshRefinement MR(geo_sf_mesh, geo_b_mesh, args, state);
    MR.tet_vertices = std::move(tet_vertices);
    MR.tets = std::move(tet_indices);
    MR.is_surface_fs = std::move(is_surface_facet);
    MR.prepareData();
    logger().info("Refinement initialization done!");

    // improvement
    MR.refine(state.ENERGY_AMIPS);

    tet_vertices = std::move(MR.tet_vertices);
    tet_indices = std::move(MR.tets);
    is_surface_facet = std::move(MR.is_surface_fs);
    t_is_removed = std::move(MR.t_is_removed);
    tetQuality = std::move(MR.tet_qualities);

    double stageTime = igl_timer.getElapsedTime();
    logger().info("Total time for the optimization stage = {}s", stageTime);
    logger().set_level(verbose_level);
}

////////////////////////////////////////////////////////////////////////////////

void tetrahedralization(const Eigen::MatrixXd &VI, const Eigen::MatrixXi &FI,
                        Eigen::MatrixXd &VO, Eigen::MatrixXi &TO, Eigen::VectorXd &AO, Eigen::VectorXi &LO,
                        const Args &args) {

    GEO::initialize();

    igl::Timer igl_timer;
    igl_timer.start();

    State state(args, VI);
    GEO::Mesh geo_sf_mesh;  // surface
    GEO::Mesh geo_b_mesh;  // open boundary
    std::vector<TetVertex> tet_vertices;
    std::vector<std::array<int, 4>> tet_indices;
    std::vector<std::array<int, 4>> is_surface_facet;
    std::vector<std::array<int, 4>> face_on_shell;  // 0 for not on shell, records which surface this face is on
    std::vector<std::array<int, 4>> faceIdx_on_shell;  // -1 for not on shell, records the face index in per-layer-F
    int eulerNumber;  // Post-shell with pseudo tets

    /// STAGE 1: Preprocess & Generate tet mesh
    tetwild_stage_one(VI, FI, args, state, geo_sf_mesh, geo_b_mesh,
        tet_vertices, tet_indices, is_surface_facet, face_on_shell, faceIdx_on_shell);

    /// STAGE 1.5: Shell
    Eigen::VectorXi labels;
    tetwild_stage_shell(args, VI, FI, tet_vertices, tet_indices, is_surface_facet, face_on_shell, faceIdx_on_shell, labels, eulerNumber);
    if (args.tet_mesh_sanity_check) {
        tetshell::TetMeshCheckArgs_t tetMeshCheckArgs;
        tetMeshCheckArgs.vertexAttri = false;
        tetMeshCheckArgs.conform = false;
        tetshell::TetMeshCheck tetMeshCheck(VI, FI, tet_vertices, tet_indices, labels, face_on_shell, tetMeshCheckArgs);
        tetMeshCheck.SanityCheck(eulerNumber);
    }

    // DEBUG PURPOSE
    if (args.tet_mesh_sanity_check) {
        tetshell::TetMeshCheckArgs_t tetMeshCheckArgs;
        tetshell::TetMeshCheck tetMeshCheck(VI, FI, tet_vertices, tet_indices, labels, face_on_shell, tetMeshCheckArgs);
    }
    /// STAGE 2: Mesh refinement
    std::vector<bool> t_is_removed(tet_indices.size(), false);
    std::vector<TetQuality> tetQuality;
    if (!args.skip_optim) {
        tetwild_stage_two(args, state, geo_sf_mesh, geo_b_mesh,
                          tet_vertices, tet_indices, is_surface_facet, t_is_removed, tetQuality);
        // Post-refinement check
        if (args.tet_mesh_sanity_check) {
            tetshell::TetMeshCheckArgs_t tetMeshCheckArgs;
            tetshell::TetMeshCheck tetMeshCheck(VI, FI, tet_vertices, tet_indices, t_is_removed, labels, face_on_shell, tetMeshCheckArgs);
            tetMeshCheck.ConformityCheck();
        }
    }

    // Extract to VO TO AO LO
    tetshell::ReorderVertices(VI, tet_vertices, tet_indices);
    tetshell::LabelInOut(tet_vertices, tet_indices, t_is_removed, labels);
    tetshell::ExtractMesh(args, tet_vertices, tet_indices, tetQuality, labels, t_is_removed, VO, TO, AO, LO);

    double total_time = igl_timer.getElapsedTime();
    logger().info("Total time for all stages = {}s", total_time);
}

}  // namespace tetwild
