#include "LatticeSpace.hpp"

namespace ecell4
{

LatticeSpace::LatticeSpace(const Position3& edge_lengths)
{
    this->edge_lengths_ = edge_lengths;
    vacant_ = new VacantType();
    set_lattice_properties();
}

LatticeSpace::~LatticeSpace()
{
    delete vacant_;
}

Integer LatticeSpace::num_species() const
{
    return spmap_.size();
}

bool LatticeSpace::has_species(const Species& sp) const
{
    return spmap_.find(sp) != spmap_.end();
}

Integer LatticeSpace::num_molecules(const Species& sp) const
{
    return num_particles(sp) * 1;
}

const Position3& LatticeSpace::edge_lengths() const
{
    return edge_lengths_;
}

Integer LatticeSpace::num_particles() const
{
    Integer count(0);
    for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
    {
        count += (*itr).second.voxels().size();
    }

    return count;
}

Integer LatticeSpace::num_particles(const Species& sp) const
{
    spmap::const_iterator itr(spmap_.find(sp));
    Integer count(0);
    if (itr != spmap_.end())
    {
        count = (*itr).second.voxels().size();
    }
    return count;
}

bool LatticeSpace::has_particle(const ParticleID& pid) const
{
    bool flg(false);
    for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
    {
        const MolecularTypeBase& mt((*itr).second);
        if (mt.is_vacant())
        {
            return false;
        }
        for (MolecularType::container_type::const_iterator vitr(mt.begin());
                vitr != mt.end(); ++vitr)
        {
            if ((*vitr).second == pid)
            {
                flg = true;
                break;
            }
        }
    }

    return flg;
}

std::vector<std::pair<ParticleID, Particle> >
    LatticeSpace::list_particles() const
{
    std::vector<std::pair<ParticleID, Particle> > retval;
    for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
    {
        const MolecularTypeBase& mt((*itr).second);
        if (mt.is_vacant())
        {
            continue;
        }
        for (MolecularType::container_type::const_iterator vitr(mt.begin());
                vitr != mt.end(); ++vitr)
        {
            retval.push_back(std::pair<ParticleID, Particle>(
                        (*vitr).second, particle_at((*vitr).first)));
        }
    }

    return retval;
}

std::vector<std::pair<ParticleID, Particle> >
    LatticeSpace::list_particles(const Species& sp) const
{
    std::vector<std::pair<ParticleID, Particle> > retval;
    spmap::const_iterator itr(spmap_.find(sp));
    if (itr != spmap_.end())
    {
        const MolecularTypeBase& mt((*itr).second);
        if (mt.is_vacant())
        {
            return retval;
        }
        for (MolecularType::container_type::const_iterator vitr(mt.begin());
                vitr != mt.end(); ++vitr)
        {
            retval.push_back(std::pair<ParticleID, Particle>(
                        (*vitr).second, particle_at((*vitr).first)));
        }
    }

    return retval;
}

bool LatticeSpace::update_particle(const ParticleID& pid, const Particle& p)
{
    Integer coord(position2coord(p.position()));
    MolecularTypeBase* dest_mt = get_molecular_type(p.species());

    if (voxels_.size() - coord <= 0)
    {
        return false;
    }
    if (has_particle(pid))
    {
        Integer coord2(get_coord(pid));
        MolecularTypeBase* src_ptr_mt(voxels_.at(coord));
        src_ptr_mt->removeVoxel(coord2);
        voxel_container::iterator itr(voxels_.begin() + coord2);
        voxels_.erase(itr);
        voxels_.insert(itr, NULL);
    }
    dest_mt->addVoxel(MolecularType::particle_info(coord, pid));

    voxel_container::iterator itr(voxels_.begin() + coord);
    (*itr) = dest_mt;

    return true;
}


/*
 * derived from SpatiocyteStepper::setLatticeProperties()
 */
void LatticeSpace::set_lattice_properties()
{
    lattice_type_ = HCP_LATTICE;

    theNormalizedVoxelRadius = 2.5e-9;

    HCP_L = theNormalizedVoxelRadius/sqrt(3);
    HCP_X = theNormalizedVoxelRadius*sqrt(8.0/3); //Lx
    HCP_Y = theNormalizedVoxelRadius*sqrt(3); //Ly

    Real lengthX = edge_lengths_[0];
    Real lengthY = edge_lengths_[1];
    Real lengthZ = edge_lengths_[2];

    /*
    row_size_ = (Integer)rint((lengthZ/2)/theNormalizedVoxelRadius) + 4;
    layer_size_ = (Integer)rint(lengthY/HCP_Y) + 4;
    col_size_ = (Integer)rint(lengthX/HCP_X) + 4;
    */
    row_size_ = (Integer)rint((lengthZ/2)/theNormalizedVoxelRadius) + 2;
    layer_size_ = (Integer)rint(lengthY/HCP_Y) + 2;
    col_size_ = (Integer)rint(lengthX/HCP_X) + 2;

    MolecularTypeBase* border = get_molecular_type(Species("Border", "0"));
    for (Integer coord(0); coord < row_size_ * layer_size_ * col_size_; ++coord)
    {
        Global global(coord2global(coord));
        if (global.col == 0 || global.col == col_size_-1 ||
                global.row == 0 || global.row == row_size_-1 ||
                global.layer == 0 || global.layer == layer_size_-1)
        {
            //voxels_.push_back(border_);
            voxels_.push_back(border);
        }
        else
        {
            voxels_.push_back(vacant_);
        }
    }
}

/*
 * original methods
 */


std::vector<Species> LatticeSpace::list_species() const
{
    std::vector<Species> keys;
    for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
    {
        keys.push_back((*itr).first);
    }
    return keys;
}

MolecularTypeBase* LatticeSpace::get_molecular_type(const Species& sp)
{
    spmap::iterator itr(spmap_.find(sp));
    if (itr == spmap_.end())
    {
        MolecularType mt(sp);
        std::pair<spmap::iterator, bool> result = spmap_.insert(spmap::value_type(sp, mt));
        if (result.second)
            itr = result.first;
        else
            throw "insert error";
    }
    return &((*itr).second);
}

const Global LatticeSpace::position2global(const Position3& pos) const
{
    Global global;
    switch(lattice_type_)
    {
        case HCP_LATTICE:
            global.col = (Integer)(pos[0] / HCP_X);
            global.layer = (Integer)((pos[1] - (global.col % 2) * HCP_L) / HCP_Y);
            global.row = (Integer)((pos[2] - (((global.layer + global.col) % 2) *
                            theNormalizedVoxelRadius)) /
                    theNormalizedVoxelRadius / 2);
            break;

        case CUBIC_LATTICE:
            global.layer = (Integer)(pos[1] / theNormalizedVoxelRadius / 2);
            global.row = (Integer)(pos[2] / theNormalizedVoxelRadius / 2);
            global.col = (Integer)(pos[0] / theNormalizedVoxelRadius / 2);
            break;
    }
    return global;
}

const Global LatticeSpace::coord2global(Integer aCoord) const
{
    Global retval;
    retval.col = aCoord / (row_size_ * layer_size_);
    retval.layer = (aCoord % (row_size_ * layer_size_)) / row_size_;
    retval.row = (aCoord % (row_size_ * layer_size_)) % row_size_;
    return retval;
}

const Position3 LatticeSpace::global2position(const Global& global) const
{
    //the center point of a voxel
    Position3 position;
    switch(lattice_type_)
    {
        case HCP_LATTICE:
            position[0] = global.col * HCP_X;
            position[1] = (global.col % 2) * HCP_L + HCP_Y * global.layer;
            position[2] = (global.row * 2 + (global.layer + global.col) % 2)
                * theNormalizedVoxelRadius;
            break;
        case CUBIC_LATTICE:
            position[0] = global.col * 2 * theNormalizedVoxelRadius;
            position[1] = global.layer * 2 * theNormalizedVoxelRadius;
            position[2] = global.row * 2 * theNormalizedVoxelRadius;
            break;
    }
    return position;
}

const Position3 LatticeSpace::coord2position(Integer coord) const
{
    return global2position(coord2global(coord));
}

Integer LatticeSpace::position2coord(const Position3& pos) const
{
    return global2coord(position2global(pos));
}

Integer LatticeSpace::get_coord(const ParticleID& pid) const
{
    for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
    {
        const MolecularTypeBase& mt((*itr).second);
        if (mt.is_vacant())
        {
            return false;
        }
        for (MolecularType::container_type::const_iterator vitr(mt.begin());
                vitr != mt.end(); ++vitr)
        {
            if ((*vitr).second == pid)
            {
                return (*vitr).first;
            }
        }
    }
    throw "Exception: Not in lattice";
}

MolecularTypeBase* LatticeSpace::get_molecular_type(Integer coord) const
{
    return voxels_.at(coord);
}

bool LatticeSpace::add(const Species& sp)
{
    if (has_species(sp))
    {
        return false;
    }
    MolecularType mt(sp);
    std::pair<spmap::iterator, bool> retval = spmap_.insert(spmap::value_type(sp, mt));
    return retval.second;
}

bool LatticeSpace::add(const Species& sp, Coord coord, const ParticleID& pid) throw(std::out_of_range)
{
    if (!is_in_range(coord))
    {
        throw std::out_of_range("");
    }

    MolecularTypeBase* mt(get_molecular_type(sp));
    if (mt->is_vacant())
    {
        return false;
    }
    MolecularTypeBase* mt_at(get_molecular_type(coord));
    if (!mt_at->is_vacant())
    {
        return false;
    }
    /*
     * Warning: Not Checking duplication of ParticleID
     */
    mt->addVoxel(std::pair<Coord, ParticleID>(coord, pid));

    voxel_container::iterator itr(voxels_.begin() + coord);
    (*itr) = mt;

    return true;
}

bool LatticeSpace::move(Coord from, Coord to) throw(std::out_of_range)
{
    if (!is_in_range(from) || !is_in_range(to))
    {
        throw std::out_of_range("");
    }

    MolecularTypeBase* to_mt(get_molecular_type(to));

    if (!to_mt->is_vacant()) {
        return false;
    }

    MolecularTypeBase* from_mt(get_molecular_type(from));

    if (!from_mt->is_vacant()) {
        MolecularType::container_type::iterator itr(from_mt->find(from));
        (*itr).first = to;
        voxel_container::iterator fitr(voxels_.begin() + from);
        (*fitr) = vacant_;
        voxel_container::iterator titr(voxels_.begin() + to);
        (*titr) = from_mt;
    }

    return true;
}

bool LatticeSpace::react(Coord coord, const Species& species) throw(std::out_of_range)
{
    if (!is_in_range(coord))
    {
        throw std::out_of_range("");
    }

    MolecularTypeBase* mt(get_molecular_type(coord));
    MolecularTypeBase* new_mt(get_molecular_type(species));
    if (mt->is_vacant())
    {
        return false;
    }
    else
    {
        MolecularType::particle_info info(*mt->find(coord));
        mt->removeVoxel(coord);
        new_mt->addVoxel(info);
        voxel_container::iterator itr(voxels_.begin() + coord);
        (*itr) = new_mt;
    }
    return true;
}

Integer LatticeSpace::global2coord(const Global& global) const
{
    return global.row +
        row_size_ * global.layer +
        row_size_ * layer_size_ * global.col;
}

const Particle LatticeSpace::particle_at(Integer coord) const
{
    const MolecularTypeBase* ptr_mt(get_molecular_type(coord));
    const Species& sp = ptr_mt->species();
    const Position3& pos = coord2position(coord);
    const Real& radius = 0;
    const Real& D = 0;
    Particle particle(sp, pos, radius, D);
    return particle;
}

bool LatticeSpace::is_in_range(Coord coord) const
{
    return coord >= 0 && coord < row_size_ * layer_size_ * col_size_;
}

} // ecell4
