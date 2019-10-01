#include <ecell4/core/Polygon.hpp>
#include <ecell4/core/exceptions.hpp>
#include <boost/container/flat_map.hpp>
#include <boost/container/static_vector.hpp>
#include <boost/format.hpp>

namespace ecell4
{

const Real Polygon::absolute_tolerance = 1e-12;
const Real Polygon::relative_tolerance = 1e-8;

void Polygon::assign(const std::vector<Triangle>& ts)
{
    constexpr Real pi = boost::math::constants::pi<Real>();
    const Real tol_abs2 = absolute_tolerance * absolute_tolerance;
    const Real tol_rel2 = relative_tolerance * relative_tolerance;

    vertices_.clear();
       faces_.clear();
       edges_.clear();
    this->total_area_ = 0.0;

    // prepair temporal data storage
    typedef std::pair<FaceID, std::size_t>                     fid_vidx_pair;
    typedef std::pair<Real3, std::vector<fid_vidx_pair>>       tmp_vtx_type;
    typedef boost::container::flat_map<VertexID, tmp_vtx_type> tmp_vertex_map;
    tmp_vertex_map tmp_vtxs;

    // first, generate (FaceIDs for all triangles) and (EdgeIDs for all Edges).
    // and collect vertices that are at the same position.
    for(const Triangle& triangle : ts)
    {
        this->total_area_ += triangle.area();

        const FaceID fid = FaceID(faces_.size());
        face_data fd;
        fd.triangle = triangle;

        for(std::size_t i=0; i<3; ++i)
        {
            const Real3& v1 = triangle.vertices()[i];
            boost::optional<VertexID> found_vtx = boost::none;

            // find near vertex
            for(auto& vid_vtx : tmp_vtxs)
            {
                const Real3&  v2 = vid_vtx.second.first;
                const Real dist2 =
                    length_sq(this->periodic_transpose(v1, v2) - v2);

                if(dist2 < tol_abs2 || dist2 < tol_rel2 * length_sq(v1))
                {
                    // vertex that locates near the vertex found
                    found_vtx = vid_vtx.first;
                    auto& vtx = vid_vtx.second;

                    // calculating mean position...
                    vtx.first = (v2 * vtx.second.size() + this->apply_boundary(v1)) /
                                (vtx.second.size() + 1);
                    // assign face-id to the vertex
                    vtx.second.push_back(std::make_pair(fid, i));
                    break;
                }
            }
            if(!found_vtx) // new vertices! add VertexID.
            {
                const VertexID new_vid = VertexID(tmp_vtxs.size());
                tmp_vtxs[new_vid] = std::make_pair(v1,
                        std::vector<fid_vidx_pair>(1, std::make_pair(fid, i)));
                found_vtx = new_vid;
            }
            fd.vertices[i] = *found_vtx; // store vertex id to face data
        }

        // make 3 edges around the face
        for(std::size_t i=0; i<3; ++i)
        {
            // in this point, edge length and direction are not fixed (because
            // vertex positions are corrected after all the faces are assigned).
            const EdgeID eid = EdgeID(edges_.size());
            edge_data ed;
            ed.face   = fid;
            ed.target = fd.vertices[i==2?0:i+1];
            this->edges_.push_back(ed);

            fd.edges[i] = eid;
        }
        // set `next` of these 3 edges
        for(std::size_t i=0; i<3; ++i)
        {
            this->edge_at(fd.edges[i]).next = fd.edges[i==2?0:i+1];
        }
        faces_.push_back(fd);
    }

    // * assign tmp_vtxs to this->vertices_
    // * set outgoing_edges without order
    for(const auto& vid_vtx : tmp_vtxs)
    {
        const VertexID                         vid = vid_vtx.first;
        const Real3                            pos = vid_vtx.second.first;
        const std::vector<fid_vidx_pair>& face_pos = vid_vtx.second.second;

        vertex_data vd;
        vd.position = pos;

        // * set vertex.outgoing_edges, but not sorted.
        for(const auto& fid_vidx : face_pos)
        {
            const FaceID      fid = fid_vidx.first;
            const std::size_t idx = fid_vidx.second;
            face_data& fd = this->face_at(fid);

            assert(vid == fd.vertices[idx]);
            vd.outgoing_edges.push_back(std::make_pair(fd.edges[idx], 0.0));
        }
        this->vertices_.push_back(vd);
    }

    // * refine vertex positions
    for(face_data& fd : this->faces_)
    {
        boost::array<Real3, 3> vs = fd.triangle.vertices();
        vs[0] = this->periodic_transpose(
                this->vertex_at(fd.vertices[0]).position, vs[0]);
        vs[1] = this->periodic_transpose(
                this->vertex_at(fd.vertices[1]).position, vs[1]);
        vs[2] = this->periodic_transpose(
                this->vertex_at(fd.vertices[2]).position, vs[2]);
        fd.triangle = Triangle(vs);
    }

    // set edge.length, edge.direction by using face.traingle
    for(const face_data& fd : this->faces_)
    {
        for(std::size_t i=0; i<3; ++i)
        {
            const EdgeID eid = fd.edges[i];
            this->edge_at(eid).length    = fd.triangle.length_of_edge_at(i);
            this->edge_at(eid).direction = fd.triangle.edge_at(i);
        }
    }

    // search pairs of opposite edges & calculate edge.tilt.
    for(std::size_t i=0; i<edges_.size(); ++i)
    {
        const EdgeID   eid = EdgeID(i);
        const VertexID start  = this->target_of(eid);
        const VertexID target = this->target_of(
                this->next_of(this->next_of(eid)));

        bool opposite_found = false;
        for(const auto& eid_angle : this->vertex_at(start).outgoing_edges)
        {
            const EdgeID outgoing = eid_angle.first;
            if(this->target_of(outgoing) == target)
            {
                // found opposite edge! calculate tilt...
                this->edge_at(eid).opposite_edge = outgoing;

                const FaceID fid1 = face_of(eid);
                const FaceID fid2 = face_of(outgoing);
                const Real3 n1 = this->face_at(fid1).triangle.normal();
                const Real3 n2 = this->face_at(fid2).triangle.normal();
                const Real3 cr = cross_product(this->edge_at(eid).direction, n1);
                const Real  sg = dot_product(cr, n2);
                const Real ang = calc_angle(n1, n2) * (sg > 0 ? 1 : -1);

                this->edges_.at(i     ).tilt = ang;
                this->edge_at(outgoing).tilt = ang;

                opposite_found = true;
                break;
            }
        }
        if(!opposite_found)
        {
            throw ecell4::NotSupported("The given polygon is not closed.");
        }
    }

    // set vertex_data.angle by traversing edges.
    for(vertex_data& vtx : this->vertices_)
    {
        const std::size_t num_edges = vtx.outgoing_edges.size();
        std::vector<EdgeID> outgoing_edges_tmp(vtx.outgoing_edges.size());
        for(std::size_t idx=0; idx<vtx.outgoing_edges.size(); ++idx)
        {
            outgoing_edges_tmp[idx] = vtx.outgoing_edges[idx].first;
        }
        vtx.outgoing_edges.clear();

        Real total_angle = 0.0;
        const EdgeID start = outgoing_edges_tmp.front();
        EdgeID current = start;
        do
        {
            {
                const auto found = std::find(outgoing_edges_tmp.begin(),
                        outgoing_edges_tmp.end(), current);
                assert(found != outgoing_edges_tmp.end());
                outgoing_edges_tmp.erase(found);
            }
            const face_data& f = this->face_at(this->face_of(current));
            Real angle = std::numeric_limits<Real>::max();
            for(std::size_t idx=0; idx<3; ++idx)
            {
                if(f.edges[idx] == current)
                {
                    angle = f.triangle.angle_at(idx);
                    break;
                }
            }
            assert(angle != std::numeric_limits<Real>::max());

            total_angle += angle;
            vtx.outgoing_edges.push_back(std::make_pair(current, angle));
            current = this->opposite_of(this->next_of(this->next_of(current)));
        }
        while(current != start);

        vtx.apex_angle = total_angle;

        if(!outgoing_edges_tmp.empty())
        {
            std::cout << outgoing_edges_tmp.size() << std::endl;
            for(const auto& oet: outgoing_edges_tmp)
            {
                const auto fid = this->face_of(oet);
                std::cout << oet << ": on " << fid << ", {";
                const face_data& f = this->face_at(fid);
                for(std::size_t idx=0; idx<3; ++idx)
                {
                    std::cout << f.triangle.vertex_at(idx);
                    if(idx != 2) {std::cout << ", ";}
                }
                std::cout << "}, ";
                for(std::size_t idx=0; idx<3; ++idx)
                {
                    if(f.edges[idx] == oet)
                    {
                        std::cout << f.triangle.vertex_at(idx) << " -> ";
                        std::cout << f.triangle.vertex_at(idx==2?0:idx+1);
                        std::cout << std::endl;
                        break;
                    }
                }
            }
        }
        if(not outgoing_edges_tmp.empty())
        {
            throw std::runtime_error("Polygon::assign: internal error: "
                    "cannot traverse all the outgoing edges from a vertex");
        }
        if(vtx.outgoing_edges.size() != num_edges)
        {
            throw std::runtime_error("Polygon::assign: internal error: "
                    "inconsistent number of outgoing edges");
        }
    }

    // make neighbor list for faces!
    for(face_data& face : this->faces_)
    {
        for(std::size_t i=0; i<3; ++i)
        {
            const VertexID vid = face.vertices[i];
            const Real3 v_pos  = face.triangle.vertex_at(i);
            const Real3 normal = face.triangle.normal();

            // counter clock wise
            {
                const Real3 ref_edge = face.triangle.edge_at(i==0?2:i-1) /
                    (-length(face.triangle.edge_at(i==0?2:i-1)));

                face.neighbor_ccw[i].clear();
                const auto start_edge  = face.edges[i];
                EdgeID   current_edge  = opposite_of(next_of(next_of(start_edge)));
                Real     current_angle = 0.0;
                while(current_edge != start_edge)
                {
                    const FaceID fid  = face_of(current_edge);

                    const std::size_t vidx0 = this->face_at(fid).index_of(vid);
                    const std::size_t vidx1 = this->face_at(fid).index_of(
                                              target_of(current_edge));
                    const std::size_t vidx2 = this->face_at(fid).index_of(
                                              target_of(next_of(current_edge)));
                    assert(vidx0 < 3);
                    assert(vidx1 < 3);
                    assert(vidx2 < 3);

                    const Real next_angle = current_angle +
                        this->face_at(fid).triangle.angle_at(vidx0);

                    boost::array<Real3, 3> unfolded;
                    unfolded[vidx0] = v_pos;
                    unfolded[vidx1] = v_pos +
                        rotate(current_angle, normal, ref_edge) *
                        length_of(current_edge);
                    unfolded[vidx2] = v_pos +
                        rotate(next_angle,    normal, ref_edge) *
                        length_of(next_of(next_of(current_edge)));

                    face.neighbors.push_back(fid);
                    face.neighbor_ccw[i].emplace_back(fid, Triangle(unfolded));

                    current_angle = next_angle;
                    current_edge  = opposite_of(next_of(next_of(current_edge)));
                }
            }

            // clock wise
            {
                const Real3 ref_edge = face.triangle.edge_at(i) /
                    length(face.triangle.edge_at(i));

                face.neighbor_cw[i].clear();
                const auto start_edge  = face.edges[i];
                EdgeID   current_edge  = next_of(opposite_of(start_edge));
                Real     current_angle = 0.0;
                while(current_edge != start_edge)
                {
                    const FaceID fid  = face_of(current_edge);

                    const std::size_t vidx0 = this->face_at(fid).index_of(vid);
                    const std::size_t vidx1 = this->face_at(fid).index_of(
                                              target_of(current_edge));
                    const std::size_t vidx2 = this->face_at(fid).index_of(
                                              target_of(next_of(current_edge)));
                    assert(vidx0 < 3);
                    assert(vidx1 < 3);
                    assert(vidx2 < 3);

                    const Real next_angle = current_angle +
                        this->face_at(fid).triangle.angle_at(vidx0);

                    boost::array<Real3, 3> unfolded;
                    unfolded[vidx0] = v_pos;
                    unfolded[vidx1] = v_pos +
                        rotate(-next_angle, normal, ref_edge) *
                        length_of(current_edge);
                    unfolded[vidx2] = v_pos +
                        rotate(-current_angle, normal, ref_edge) *
                        length_of(next_of(next_of(current_edge)));

                    face.neighbors.push_back(fid);
                    face.neighbor_cw[i].emplace_back(fid, Triangle(unfolded));

                    current_angle = next_angle;
                    current_edge  = next_of(opposite_of(current_edge));
                }
            }
        }
        std::sort(face.neighbors.begin(), face.neighbors.end());
        const auto last = std::unique(face.neighbors.begin(), face.neighbors.end());
        face.neighbors.erase(last, face.neighbors.end());
    }
    return;
}

Real Polygon::distance_sq(const std::pair<Real3, FaceID>& pos1,
                          const std::pair<Real3, FaceID>& pos2) const
{
    typedef utils::pair_first_element_unary_predicator<FaceID, Triangle>
            face_finder_type;
    constexpr Real pi = boost::math::constants::pi<Real>();

    // if two particles are on the same face, return just a 3D distance.
    if(pos1.second == pos2.second)
    {
        return length_sq(pos2.first - pos1.first);
    }

    // If positions are on different faces, there can be several cases.
    // 1.)  ______
    //     /\    /\   | The positions are on the faces connected by a vertex.
    //    /  \  /p2\  | using p1-vtx-p2 angle and the low-of-cosines, the
    //   /____\/____\ | minimum distance on the surface can be calculated.
    //        /\    / |
    //       /p1\  /  |
    //      /____\/   |
    //
    // 2.)     ______
    //        /\ p2 / | The positions are on the faces that are not connected
    //       /  \  /  | by any vertex. There can be several options to unfold
    //      /____\/   | the polygon to make the particles on the same plane.
    //     /\    /    | In this case, finding the minimum path is too difficult
    //    /p1\  /     | to use in simulation, so just returns inf. In the SGFRD
    //   /____\/      | simulation, movement from p1 to p2 is inhibited.
    //
    // 3.)  ______
    //     /\    /\   | The positions are on the faces connected by a vertex
    //    /p2\  /  \  | and the apex angle exceeds 360 degree. There can be 2
    //   /____\/____\ | pathes, counter-clockwise and clockwise, and both angle
    //        /\    / | exceeds 180 degree. In this case, the minimum distance
    //       /p1\  /  | pathway goes across the vertex.
    //      /____\/   |
    //
    // 4.)  ......... ______
    //     /\connected\    /\   | The particles are on the faces connected by a
    //    /  \ <=====> \  /  \  | vertex and the apex angle exceeds 360 degree.
    //   /____\.........\/____\ | And the triangles overlaps when unfolded.
    //   \    /\        /\    / | In this case, we must use the minimum angle
    //    \  /p2\      /p1\  /  | because just unfolding triangle makes minimum
    //     \/____\    /____\/   | pathway shorter than the 3D distance.
    //
    // 5.) TODO
    //     \`.p1           | If a polygon is severely deformed, the minimum
    //      \o`.           | angle pathway can protrude the face. To avoid this,
    //       ^  `.         | It is required to check the minimum angle pathway
    //       |\___`.vertex | goes across all the edges.
    //       |/   .'       |
    //       v  .'         |
    //      /o.'           |
    //     /.'p2           |
    //

    const Real3& p1 = pos1.first;
    const Real3& p2 = pos2.first;
    const FaceID f1 = pos1.second;
    const FaceID f2 = pos2.second;

    // for comparison
    const auto min_edge_length = std::min(edge_length_[0],
            std::min(edge_length_[1], edge_length_[2]));
    const auto rel_tol = relative_tolerance * min_edge_length;

    const face_data& face = face_at(f1);

    boost::container::static_vector<VertexID, 3> connecting_vtxs;

    Real distance_sq = std::numeric_limits<Real>::infinity();
    for(std::size_t i=0; i<3; ++i)
    {
        // counter clockwise
        //
        //          ^ vertices[i]
        // edge[i] /|\
        //        / |~\
        //       /  o  \ next(next(egde[i]))
        //      v______>\
        //    next(edge[i])

        const VertexID    vid = face.vertices[i];
        const Real3&     vpos = position_at(vid);
        const Real3     p1tov = this->periodic_transpose(vpos, p1) - p1;
        const Real3     vtop2 = this->periodic_transpose(p2, vpos) - vpos;

        // check p1 or p2 are exactly on the vertex.
        // If they are on, it causes NaN because the length of v->p vector is 0.
        const Real p1tov_len = length(p1tov);
        const Real p2tov_len = length(vtop2);
        if(p1tov_len < relative_tolerance * min_edge_length)
        {
            return this->distance_sq(std::make_pair(p1, vid), pos2);
        }
        if(p2tov_len < relative_tolerance * min_edge_length)
        {
            return this->distance_sq(pos1, std::make_pair(p2, vid));
        }

        // calc the initial angle
        Real angle = calc_angle(p1tov, face.triangle.edge_at(i==0?2:i-1));
        const Real apex_angle = apex_angle_at(vid);

        // ------------------------------------------------------------------
        // search `f2`
        bool connected = false;
        for(const auto& neighbor : face.neighbor_ccw[i])
        {
            assert(neighbor.first != f1);

            const auto& nface = face_at(neighbor.first);
            const auto  vidx  = nface.index_of(vid);
            if(neighbor.first == f2)
            {
                angle += calc_angle(vtop2, nface.triangle.edge_at(vidx));
                connected = true;
                break;
            }
            angle += nface.triangle.angle_at(vidx);
        }
        if(!connected)
        {
            continue;
        }
        connecting_vtxs.push_back(vid);

        // ------------------------------------------------------------------
        // calculate the minimum angle

        const Real angle_ccw = angle;
        const Real angle_cw  = apex_angle - angle;
        assert(angle_cw >= 0.0);
        const Real min_angle = std::min(angle_ccw, angle_cw);

        // skip case 3.
        if(min_angle > pi) {continue;}

        // if the minimum angle < 180 degree, its case 1.
        // calculate distance using the low of cosine.

        // TODO
        // But... before calculating the distance, check whether this is the
        // case 5.
        // TODO

        const auto p1tov_lensq = length_sq(p1tov);
        const auto p2tov_lensq = length_sq(vtop2);
        const auto dist_sq = p1tov_lensq + p2tov_lensq -
            2 * std::sqrt(p1tov_lensq * p2tov_lensq) * std::cos(min_angle);

        distance_sq = std::min(distance_sq, dist_sq);
    }
    if(distance_sq != std::numeric_limits<Real>::infinity())
    {
        // distance_sq is updated! It founds the minimum path (case 1).
        return distance_sq;
    }
    if(connecting_vtxs.empty())
    {
        // if `f1` and `f2` are not connected via any vertex, the distance
        // cannot be calculated (case 2).
        return std::numeric_limits<Real>::infinity();
    }

    // Here, the positions are connected via some vertices but the minimum
    // distance is non-trivial. It is case 3. return distance that passes
    // through the vertex that connects `f1` and `f2`.

    for(const VertexID vtx : connecting_vtxs)
    {
        const Real3 vpos = position_at(vtx);
        const Real  l1   = length(this->periodic_transpose(p1, vpos) - vpos);
        const Real  l2   = length(this->periodic_transpose(p2, vpos) - vpos);
        distance_sq = std::min(distance_sq, (l1 + l2) * (l1 + l2));
    }
    return distance_sq;
}

Real Polygon::distance_sq(const std::pair<Real3, VertexID>& pos1,
                          const std::pair<Real3, FaceID>&   pos2) const
{
    const std::vector<std::pair<EdgeID, Real> >& outs =
        this->vertex_at(pos1.second).outgoing_edges;

    for(std::vector<std::pair<EdgeID, Real> >::const_iterator
            i(outs.begin()), e(outs.end()); i!=e; ++i)
    {
        const EdgeID eid(i->first);
        const FaceID fid = this->face_of(eid); // face that connects to the vtx
        if(fid == pos2.second)
        {
            // on the same face.
            return length_sq(
                this->periodic_transpose(pos1.first, pos2.first) - pos2.first);
        }
        const FaceID adj = this->face_of(this->opposite_of(this->next_of(eid)));
        if(adj == pos2.second)
        {
            const face_data& fd = this->face_at(fid);
            const std::pair<FaceID, Triangle>& fp =
                fd.neighbor_cw.at(fd.index_of(target_of(eid))).front();
            assert(fp.first == pos2.second);

            const Barycentric b2 =
                to_barycentric(pos2.first, this->face_at(pos2.second).triangle);
            const Real3 pos2rot = to_absolute(b2, fp.second);

            return length_sq(
                    this->periodic_transpose(pos1.first, pos2rot) - pos2rot);
        }
    }
    return std::numeric_limits<Real>::infinity();
}

Real Polygon::distance_sq(const std::pair<Real3, VertexID>& pos1,
                          const std::pair<Real3, VertexID>& pos2) const
{
    if(pos1.second == pos2.second)
    {
        return 0.0;
    }

    const std::vector<std::pair<EdgeID, Real> >& outs =
        this->vertex_at(pos1.second).outgoing_edges;

    for(std::vector<std::pair<EdgeID, Real> >::const_iterator
            i(outs.begin()), e(outs.end()); i!=e; ++i)
    {
        const EdgeID eid(i->first);
        if(this->target_of(eid) == pos2.second)
        {
            // directly connected.
            const Real l = length_of(eid);
            return l * l;
        }
        const EdgeID inbetween = this->opposite_of(this->next_of(eid));
        if(this->target_of(this->next_of(inbetween)) == pos2.second)
        {
            const FaceID fid = this->face_of(eid);
            const FaceID adj = this->face_of(inbetween);

            const face_data& fd = this->face_at(fid);
            const std::pair<FaceID, Triangle>& fp =
                fd.neighbor_cw.at(fd.index_of(target_of(eid))).front();
            assert(fp.first == adj);

            const Barycentric b2 =
                to_barycentric(pos2.first, this->face_at(adj).triangle);
            const Real3 pos2rot = to_absolute(b2, fp.second);

            return length_sq(
                    this->periodic_transpose(pos1.first, pos2rot) - pos2rot);
        }
    }
    return std::numeric_limits<Real>::infinity();
}

Real Polygon::distance(const std::pair<Real3, VertexID>& pos1,
                       const std::pair<Real3, VertexID>& pos2) const
{
    if(pos1.second == pos2.second)
    {
        return 0.0;
    }

    const std::vector<std::pair<EdgeID, Real> >& outs =
        this->vertex_at(pos1.second).outgoing_edges;

    for(std::vector<std::pair<EdgeID, Real> >::const_iterator
            i(outs.begin()), e(outs.end()); i!=e; ++i)
    {
        const EdgeID eid(i->first);
        if(this->target_of(eid) == pos2.second)
        {
            // directly connected.
            return length_of(eid);
        }
        const EdgeID inbetween = this->opposite_of(this->next_of(eid));
        if(this->target_of(this->next_of(inbetween)) == pos2.second)
        {
            const FaceID fid = this->face_of(eid);
            const FaceID adj = this->face_of(inbetween);

            const face_data& fd = this->face_at(fid);
            const std::pair<FaceID, Triangle>& fp =
                fd.neighbor_cw[fd.index_of(target_of(eid))].front();
            assert(fp.first == adj);

            const Barycentric b2 =
                to_barycentric(pos2.first, this->face_at(adj).triangle);
            const Real3 pos2rot = to_absolute(b2, fp.second);

            return length(
                    this->periodic_transpose(pos1.first, pos2rot) - pos2rot);
        }
    }
    return std::numeric_limits<Real>::infinity();
}


Real3 Polygon::direction(
        const std::pair<Real3, FaceID>& pos1,
        const std::pair<Real3, FaceID>& pos2) const
{
    typedef utils::pair_first_element_unary_predicator<FaceID, Triangle>
            face_finder_type;

    if(pos1.second == pos2.second)
    {
        return pos2.first - pos1.first;
    }

    const Real3& p1 = pos1.first;
    const FaceID f1 = pos1.second;
    const FaceID f2 = pos2.second;
    const Barycentric b2 = to_barycentric(pos2.first, face_at(f2).triangle);

    const face_data& face = face_at(f1);
    const Real3&   normal = face.triangle.normal();

    Real mindist2 = std::numeric_limits<Real>::max();
    Real3 direction(0,0,0);
    for(std::size_t i=0; i<3; ++i)
    {
        const VertexID vid = face.vertices[i];
        const Real3& vpos(position_at(vid));
        const Real3 vtop1(this->periodic_transpose(p1, vpos) - vpos);

        { // counter clock wise
            const std::vector<std::pair<FaceID, Triangle> >::const_iterator fi =
                std::find_if(face.neighbor_ccw[i].begin(),
                             face.neighbor_ccw[i].end(), face_finder_type(f2));
            if(fi != face.neighbor_ccw[i].end())
            {
                // unfolded place of p2
                const Real3 p2 = to_absolute(b2, fi->second);
                const Real3 vtop2(this->periodic_transpose(p2, vpos) - vpos);
                // check the angle between p1-v-p2 does not exceeds PI
                if(dot_product(normal, cross_product(vtop1, vtop2)) >= 0)
                {
                    const Real3 dr = this->periodic_transpose(p2, p1) - p1;
                    const Real  d2 = length_sq(dr);
                    if(d2 < mindist2)
                    {
                        mindist2  = d2;
                        direction = dr;
                    }
                    continue;
                }
            }
        }
        { // clock wise
            const std::vector<std::pair<FaceID, Triangle> >::const_iterator fi =
                std::find_if(face.neighbor_cw[i].begin(),
                             face.neighbor_cw[i].end(), face_finder_type(f2));
            if(fi != face.neighbor_cw[i].end())
            {
                // unfolded place of p2
                const Real3 p2 = to_absolute(b2, fi->second);
                const Real3 vtop2(this->periodic_transpose(p2, vpos) - vpos);
                // check the angle between p1-v-p2 does not exceeds PI
                if(dot_product(normal, cross_product(vtop1, vtop2)) <= 0)
                {
                    const Real3 dr = this->periodic_transpose(p2, p1) - p1;
                    const Real  d2 = length_sq(dr);
                    if(d2 < mindist2)
                    {
                        mindist2  = d2;
                        direction = dr;
                    }
                }
            }
        }
    }
    if(mindist2 == std::numeric_limits<Real>::max())
    {
        throw std::runtime_error((boost::format(
            "polygon::direction: couldn't find the min path between "
            "%1% on %2% <-> %3% on %4%") % pos1.first % pos1.second %
            pos2.first % pos2.second).str());
    }
    return direction;
}

std::pair<Real3, Polygon::FaceID> Polygon::travel(
        const std::pair<Real3, FaceID>& pos, const Real3& disp) const
{
    const Real3& p = pos.first;
    const FaceID f = pos.second;
    const Real3 np = p + disp;

    const face_data& fd = this->face_at(f);
    const Triangle& tri = fd.triangle;
    const Barycentric b2(to_barycentric(np, tri));

    // if pos + disp is inside of the current face, just return the sum.
    if(::ecell4::is_inside(b2))
    {
        // to avoid numerical error that make the particle goes outside of the
        // face that the particle belongs, use `to_absolute`.
        return std::make_pair(to_absolute(b2, tri), f);
    }

    const Barycentric b1(to_barycentric(p, tri));
    const Barycentric db(b2 - b1);
    const std::pair<std::size_t, Real> cs   = first_cross_edge(b1, db);
    const std::pair<FaceID, Triangle>& next = fd.neighbor_cw[cs.first].front();

    // if the position is inside of the adjacent face, return the position
    // reconstructed from unfolded-Barycentric by using folded Triangle.
    const Barycentric unfolded_b(to_barycentric(np, next.second));
    if(::ecell4::is_inside(unfolded_b, 1e-8))
    {
        const Real3 nxt = to_absolute(
                force_put_inside(unfolded_b), this->triangle_at(next.first));

        if(!this->is_inside_of_boundary(nxt))
        {
            std::cerr << "travel: initial pos          = " << p   << " on " << f            << std::endl;
            std::cerr << "travel: initial face         = " << f << " -> " << this->triangle_at(f) << std::endl;
            std::cerr << "travel: initial disp         = " << disp                          << std::endl;
            std::cerr << "travel: next face (unfolded) = " << next.second                   << std::endl;
            std::cerr << "travel: next barycentric crd = " << unfolded_b                    << std::endl;
            std::cerr << "travel: next bary (inside)   = " << force_put_inside(unfolded_b)  << std::endl;
            std::cerr << "travel: next face            = " << this->triangle_at(next.first) << std::endl;
            std::cerr << "travel: next pos             = " << nxt << " on " << next.first   << std::endl;
        }
        return std::make_pair(nxt, next.first);
        // use folded (normal) Triangle, NOT next.second
    }

    // stride over not only the edge but adjacent face.
    // XXX to make it sure that `on_edge` should be on the edge under the PBC,
    //     to_absolute is used with the next triangle.
    const Barycentric on_edge_b(to_barycentric(p + disp * cs.second, next.second));
    const Real3 edge_over = direction_of(fd.edges[cs.first]);
    const Real3 next_pos  = to_absolute(on_edge_b, this->triangle_at(next.first));

    if(!this->is_inside_of_boundary(next_pos))
    {
        std::cerr << "travel: initial  pos               = " << p        << " on " << f          << std::endl;
        std::cerr << "travel: initial face               = " << f << " -> " << this->triangle_at(f)             << std::endl;
        std::cerr << "travel: initial disp               = " << disp                          << std::endl;
        std::cerr << "travel: next face (unfolded)       = " << next.second                      << std::endl;
        std::cerr << "travel: next barycnetric (on edge) = " << on_edge_b                        << std::endl;
        std::cerr << "travel: next face                  = " << this->triangle_at(next.first)    << std::endl;
        std::cerr << "travel: next  pos                  = " << next_pos << " on " << next.first << std::endl;
    }

    return this->travel(std::make_pair(next_pos, next.first),
        rotate(tilt_angle_at(fd.edges[cs.first]),     // rotate disp by tilt_angle
               edge_over * (1.0 / length(edge_over)), // around the edge
               disp * (1 - cs.second)));              // the rest of displacement
}

std::pair<Real3, Polygon::FaceID> Polygon::travel(
        const std::pair<Real3, FaceID>& pos, const Real3& disp,
        const std::size_t restraint) const
{
    if(restraint == 0)
    {
        std::cerr << "movement along surface of a Polygon: "
                     "restraint hits 0. tolerance violated!" << std::endl;
        std::cerr << "the rest of displacement: " << disp    << std::endl;
        return pos;
    }
    const Real3& p = pos.first;
    const FaceID f = pos.second;
    const Real3 np = p + disp;

    const face_data& fd = this->face_at(f);
    const Triangle& tri = fd.triangle;
    const Barycentric b2(to_barycentric(np, tri));

    // if pos + disp is inside of the current face, just return the sum.
    if(::ecell4::is_inside(b2))
    {
        // to avoid numerical error that make the particle goes outside of the
        // face that the particle belongs, use `to_absolute`.
        return std::make_pair(to_absolute(b2, tri), f);
    }

    const Barycentric b1(to_barycentric(p, tri));
    const Barycentric db(b2 - b1);
    const std::pair<std::size_t, Real> cs   = first_cross_edge(b1, db);
    const std::pair<FaceID, Triangle>& next = fd.neighbor_cw[cs.first].front();

    // if the position is inside of the adjacent face, return the position
    // reconstructed from unfolded-Barycentric by using folded Triangle.
    const Barycentric unfolded_b(to_barycentric(np, next.second));
    if(::ecell4::is_inside(unfolded_b, 1e-8))
    {
        const Real3 nxt = to_absolute(
                force_put_inside(unfolded_b), this->triangle_at(next.first));

        if(!this->is_inside_of_boundary(nxt))
        {
            std::cerr << "travel: initial  pos         = " << p   << " on " << f            << std::endl;
            std::cerr << "travel: initial face         = " << f << " -> " << this->triangle_at(f)          << std::endl;
            std::cerr << "travel: initial disp         = " << disp                          << std::endl;
            std::cerr << "travel: next face (unfolded) = " << next.second                   << std::endl;
            std::cerr << "travel: next barycentric crd = " << unfolded_b                    << std::endl;
            std::cerr << "travel: next bary (inside)   = " << force_put_inside(unfolded_b)  << std::endl;
            std::cerr << "travel: next face            = " << this->triangle_at(next.first) << std::endl;
            std::cerr << "travel: next  pos            = " << nxt << " on " << next.first   << std::endl;
        }
        return std::make_pair(nxt, next.first);
        // use folded (normal) Triangle, NOT next.second
    }

    // stride over not only the edge but adjacent face.
    // XXX to make it sure that `on_edge` should be on the edge under the PBC,
    //     to_absolute is used with the next triangle.
    const Barycentric on_edge_b(to_barycentric(p + disp * cs.second, next.second));
    const Real3 edge_over = direction_of(fd.edges[cs.first]);
    const Real3 next_pos  = to_absolute(on_edge_b, this->triangle_at(next.first));

    if(!this->is_inside_of_boundary(next_pos))
    {
        std::cerr << "travel: initial  pos               = " << p        << " on " << f          << std::endl;
        std::cerr << "travel: initial face               = " << f << " -> " << this->triangle_at(f)             << std::endl;
        std::cerr << "travel: initial disp               = " << disp                          << std::endl;
        std::cerr << "travel: next face (unfolded)       = " << next.second                      << std::endl;
        std::cerr << "travel: next barycnetric (on edge) = " << on_edge_b                        << std::endl;
        std::cerr << "travel: next face                  = " << this->triangle_at(next.first)    << std::endl;
        std::cerr << "travel: next  pos                  = " << next_pos << " on " << next.first << std::endl;
    }
    return this->travel(std::make_pair(next_pos, next.first),
        rotate(tilt_angle_at(fd.edges[cs.first]),     // rotate disp by tilt_angle
               edge_over * (1.0 / length(edge_over)), // around the edge
               disp * (1 - cs.second)),               // the rest of displacement
        restraint - 1);
}

} // ecell4
