	#ifndef MARS_TRACKER_HPP
	#define MARS_TRACKER_HPP

	#include "mars_mesh.hpp"
    

namespace mars {

	class Tracker {
	public:
		Tracker()
		: current_iterate_(1), is_tracking_(false)
		{}

		void begin_iterate()
		{
			is_tracking_ = true;
		}

		void element_refined(const Integer element_id)
		{
			if(!is_tracking_) return;

			if(iterates_.size() <= element_id) {
				iterates_.resize(element_id + 1, 0);
			}

			iterates_[element_id] = current_iterate_;
		}

		void element_deactivated(const Integer element_id)
		{	
			if(iterates_.size() <= element_id) {
				return;
			}

			if(iterates_[element_id] == current_iterate_) {
				iterates_[element_id] = -current_iterate_;
			}
		}

		inline Integer get_iterate(const Integer element_id) const
		{
			assert(element_id >= 0);
			assert(element_id < iterates_.size());

			return iterates_[element_id];
		}

			template<Integer Dim, Integer ManifoldDim>
		void undo_last_iterate(Mesh<Dim, ManifoldDim> &mesh)
		{
			const Integer last_iterate = current_iterate_ - 1;
			const Integer iter_size = iterates_.size();

			std::vector<Integer> elements_to_remove;
			for(Integer i = 0; i < iter_size; ++i) {
				if(get_iterate(i) == last_iterate) {
					assert(!mesh.is_active(i));
					mesh.set_active(i, true);
					iterates_[i] = INVALID_INDEX;

					for(auto c : mesh.elem(i).children) {
						mesh.set_active(c, false);
						elements_to_remove.push_back(c);
					}
				}
			}

			mesh.remove_elements(elements_to_remove);
			--current_iterate_;
		}

		void end_iterate()
		{
			++current_iterate_;
			is_tracking_ = false;
		}

		auto current_iterate(){return current_iterate_;}

		auto has_tracked()const {if(iterates_.size()==0)
			                       return false;
			                     else return true;}

		inline Integer get_level(const Integer element_id) const
		{
			if(iterates_.size()==0)
				{

					return 0;
				}
			else
				{
					return get_iterate(element_id);
				}
		}


	private:
		bool is_tracking_;
		Integer current_iterate_;
		std::vector<Integer> iterates_;
	};


	template<typename MeshT>
	bool elem_belongs_to_level(const MeshT& mesh, const Integer i, const Integer level, const Tracker& track)
	{
		return  (!track.has_tracked()) || (track.get_iterate(i)==level || (mesh.elem(i).children.size()==0 && track.get_iterate(i)<=level));
		// return   (track.get_iterate(i)==level || (mesh.elem(i).children.size()==0 && track.get_iterate(i)<=level));
	}

	template<typename MeshT>
	bool elem_belongs_to_level(const std::shared_ptr<MeshT> mesh_ptr, const Integer i, const Integer level, const Tracker& track)
	{

		// std::cout<<"!track.has_tracked()->"<<!track.has_tracked()<<" track.get_iterate(i)==level-> "<<(track.get_iterate(i)==level)<<"children.size()==0 && track.get_iterate(i)<=level)->"<<  (mesh_ptr->elem(i).children.size()==0 && track.get_iterate(i)<=level)<<std::endl;
		return  (!track.has_tracked())||(track.get_iterate(i)==level || (mesh_ptr->elem(i).children.size()==0 && track.get_iterate(i)<=level));
		// return (track.get_iterate(i)==level || (mesh_ptr->elem(i).children.size()==0 && track.get_iterate(i)<=level));

	}


	template<class MeshT>
	class Bisection;

	// template<typename MeshT>
	// bool elem_belongs_to_level(const std::shared_ptr<MeshT> mesh_ptr, const Integer i, const Integer level, const std::shared_ptr<Bisection<MeshT>>& bisection_ptr)
	// {
	// 	if(level==-1 )
	// 	     return mesh_ptr->is_active(i);
	// 	// otherwise we check if it belongs to the given level l!=-1
	// 	else return elem_belongs_to_level(mesh_ptr,i,level,bisection_ptr->tracker());
	// }

	// template<typename MeshT>
	// bool elem_belongs_to_level(const MeshT& mesh, const Integer i, const Integer level, const std::shared_ptr<Bisection<MeshT>>& bisection_ptr)
	// {

	// 	if(level==-1 )
	// 	     return mesh.is_active(i);
	// 	// otherwise we check if it belongs to the given level l!=-1
	// 	else return elem_belongs_to_level(mesh,i,level,bisection_ptr->tracker());
	// }


	template<typename MeshT>
	void get_elems_of_level(std::vector<Integer>& list_of_elems,const MeshT& mesh,const Integer level, const Tracker& track)
	{
		list_of_elems.reserve(mesh.n_elements());

		for(Integer i=0;i<mesh.n_elements();i++)
		{
	            if(elem_belongs_to_level(mesh,i,level,track))//track.get_iterate(i)==given_level || (mesh.elem(i).children.size()==0 && track.get_iterate(i)<=given_level))
	            {
	            	list_of_elems.push_back(i);
	            }
	        }
	    }


	

	template<typename MeshT>
	void get_elems_of_level(std::vector<Integer>& list_of_elems,const std::shared_ptr<MeshT> mesh_ptr,const Integer level, const Tracker& track)
	{
		list_of_elems.reserve(mesh_ptr->n_elements());

		for(Integer i=0;i<mesh_ptr->n_elements();i++)
		{
	            if(elem_belongs_to_level(mesh_ptr,i,level,track))//track.get_iterate(i)==given_level || (mesh.elem(i).children.size()==0 && track.get_iterate(i)<=given_level))
	            {
	            	list_of_elems.push_back(i);
	            }
	        }
	}

     template<typename MeshT>
	 inline Real get_used_parent_elem_id(const Integer child_id,const MeshT & mesh,const Tracker& tracker)
	 {
	 	Integer parent_id=mesh.elem(child_id).parent_id;
        if(tracker.get_level(parent_id)<0)
        	return get_used_parent_elem_id(parent_id,mesh,tracker);
        else
        	return parent_id;
	 }

     template<typename MeshT>
	 inline void  get_used_children_elem_id_aux(std::vector<Integer>& children_full, const Integer parent_id,const MeshT & mesh,const Tracker& tracker)
	 {
	 	auto& children=mesh.elem(parent_id).children;

	 	for(Integer i=0;i<children.size();i++)
	 	{
	 		if(tracker.get_level(children[i])>=0)
	 		{
             children_full.push_back(children[i]);
	 		}
	 		else
	 		{
	 			get_used_children_elem_id_aux(children_full, children[i],mesh,tracker);
	 		}

	 	}
	 }

	 template<typename MeshT>
	 inline std::vector<Integer>  get_used_children_elem_id(const Integer parent_id,const MeshT & mesh,const Tracker& tracker)
	 {
	 	std::vector<Integer> children_full;
	 	get_used_children_elem_id_aux(children_full,parent_id,mesh,tracker);
	 	return children_full;
	 }


}

	#endif //MARS_TRACKER_HPP