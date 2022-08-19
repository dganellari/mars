#ifndef MARS_TRACKER_HPP
#define MARS_TRACKER_HPP

#include "mars_mesh.hpp"

namespace mars {

    class Tracker {
    public:
        Tracker() : current_iterate_(1), is_tracking_(false) {}

        void begin_iterate() { is_tracking_ = true; }

        void element_refined(const Integer element_id) {
            if (!is_tracking_) return;

            if (iterates_.size() <= element_id) {
                iterates_.resize(element_id + 1, 0);
            }

            iterates_[element_id] = current_iterate_;
        }

        inline Integer get_iterate(const Integer element_id) const {
            assert(element_id >= 0);
            assert(element_id < iterates_.size());

            return iterates_[element_id];
        }

        template <Integer Dim, Integer ManifoldDim>
        void undo_last_iterate(SimplicialMesh<Dim, ManifoldDim> &mesh) {
            const Integer last_iterate = current_iterate_ - 1;
            const Integer iter_size = iterates_.size();

            std::vector<Integer> elements_to_remove;
            for (Integer i = 0; i < iter_size; ++i) {
                if (get_iterate(i) == last_iterate) {
                    assert(!mesh.is_active(i));
                    mesh.set_active(i, true);
                    iterates_[i] = INVALID_INDEX;

                    for (auto c : mesh.elem(i).children) {
                        mesh.set_active(c, false);
                        elements_to_remove.push_back(c);
                    }
                }
            }

            mesh.remove_elements(elements_to_remove);
            --current_iterate_;
        }

        void end_iterate() {
            ++current_iterate_;
            is_tracking_ = false;
        }

    private:
        bool is_tracking_;
        Integer current_iterate_;
        std::vector<Integer> iterates_;
    };
}  // namespace mars

#endif  // MARS_TRACKER_HPP
