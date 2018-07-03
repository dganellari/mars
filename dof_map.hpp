#ifndef MARS_DOF_MAP_HPP
#define MARS_DOF_MAP_HPP

namespace mars {
	class Map {
	public:
		Map(const Integer partition_id)
		: partition_id_(partition_id)
		{}

		void add(
			const Integer local_id,
			const Integer global_id,
			const Integer owner = INVALID_INDEX)
		{
			if(local_id >= global_id_.size()) {
				global_id_.resize(local_id, INVALID_INDEX);
			}

			set_global_id(local_id, global_id);
		}

		void set_owner(
			const Integer local_id,
			const Integer owner)
		{
			if(local_id >= owner_id_.size()) {
				owner_id_.resize(local_id, INVALID_INDEX);
			}

			owner_id_[local_id] = INVALID_INDEX;
		}

		inline Integer node_owner(const Integer local_id) const
		{
			assert(local_id < owner_id_.size());
			return owner_id_[local_id];
		}

		void set_global_id(const Integer local_id, const Integer global_id)
		{
			assert(local_id < global_id_.size());
			global_id_[local_id] = global_id;

			if(global_id != INVALID_INDEX) {
				global_to_local_[global_id] = local_id;
			}
		}

		inline Integer local(const Integer global_id) const
		{
			auto it = global_to_local_.find(global_id);
			if(it == global_to_local_.end()) return INVALID_INDEX;

			return it->second;
		}

		inline Integer global(const Integer local_id) const
		{
			assert(local_id < global_id_.size());
			return global_id_[local_id];
		}

		inline Integer n_owned() const
		{
			Integer ret = 0;
			for(auto n : owner_id_) {
				ret += (n == partition_id());
			}

			return ret;
		}

		inline Integer n_unassigned() const
		{
			Integer ret = 0;
			for(auto n : owner_id_) {
				ret += (n == INVALID_INDEX);
			}

			return ret;
		}

		inline Integer partition_id() const
		{
			return partition_id_;
		}

	private:
		Integer partition_id_;
		std::vector<Integer> global_id_;
		std::vector<Integer> owner_id_;

		//global to local
		std::map<Integer, Integer> global_to_local_;
	};
}

#endif //MARS_DOF_MAP_HPP
