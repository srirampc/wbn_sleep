#ifndef CELL_TYPES_HPP
#define CELL_TYPES_HPP

#include <cassert>
#include <iostream>
#include <string>

#include <utils.hpp>

class CellTypes {
  public:
    static constexpr int MAXTYPES = 33;

    struct cell_t {
        int x, y;
        std::string name;
        //
        cell_t(int ix, int iy, std::string nm) : x(ix), y(iy), name(nm) {}
        cell_t() : x(-1), y(-1) {}
    };

  private:
    int n;                  // number of all types we have
    cell_t cell[MAXTYPES];  // info for each type
    double MRI_ratio;  // mainly for debugging purposes when we want simulate
                       // smaller network from hardcoded MRI input data, we
                       // shrink number of neurons by MRU_ratio factor

  public:
    CellTypes() : n(0), MRI_ratio(1) {}

    inline void set_mri_ratio(double rx) { MRI_ratio = rx; }

    inline int size() const { return n; }

    inline const CellTypes::cell_t& get_cell(std::size_t cell_id) const {
        return cell[cell_id];
    }

    inline const CellTypes::cell_t& operator[](std::size_t cell_id) const {
        return cell[cell_id];
    }

    inline CellTypes::cell_t& operator[](std::size_t cell_id) {
        return cell[cell_id];
    }

    inline void add_type(const CellTypes::cell_t& pcx) {
        assert(n < CellTypes::MAXTYPES);
        cell[n] = pcx;
        n++;
    }

    inline void apply_reduction(double ratio_x, double ratio_y) {
        // ratios after loop because MRI scenario
        for (int i = 0; i < n; i++) {
            cell[n].x = MathUtil::round(cell[n].x * ratio_x);
            cell[n].y = MathUtil::round(cell[n].y * ratio_y);
        }
    }

    inline bool right_hemisphere(int cellt) const {
        // is the cell in left or right hemisphere (left is
        // default for networks without hemispheric difference)
        return (cell[cellt].name.find("r") == 0);
    }  // right_hemisphere

    inline int get_type(const std::string& t) const {  // map name into index
        for (int i = 0; i < n; i++)
            if (cell[i].name == t)
                return i;
        std::cerr << "CAN NOT FIND NEURON : " << t << std::endl;
        assert(0);  // no such neuron type exists
        return -1;
    }  // get_type

    inline void dump(std::ostream& ConnSummaryFile) const {
        for (int i = 0; i < n; i++) {
            dump(i, ConnSummaryFile);
            ConnSummaryFile << "\n";
        }
    }
    inline void dump(int i, std::ostream& ConnSummaryFile) const {
        assert(i < n);
        ConnSummaryFile << cell[i].name << " " << cell[i].x << " " << cell[i].y;
    }
    inline void output(std::ostream& out_stream)
        const {  // fixed output format used by simulation
        out_stream << n << "\n";
        for (int i = 0; i < n; i++)
            out_stream << cell[i].x << " " << cell[i].y << "\n";
    }

    inline void MRI_reduce() {
        for (int i = 0; i < n; i++)
            cell[i].x *= MRI_ratio;
    }

    inline bool has_both_hemispheres() const {
        for (int i = 0; i < n; i++)
            if (cell[i].name.find("r") == 0)
                return true;
        return false;
    }
};  // CellTypes

#endif
